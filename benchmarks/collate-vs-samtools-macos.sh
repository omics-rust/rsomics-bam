#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

if [[ $# -lt 4 || $# -gt 6 ]]; then
    echo "usage: $0 RSOMICS_BAM SAMTOOLS INPUT OUTPUT_DIR [REPEATS] [MODE]" >&2
    exit 2
fi

ours=$1
samtools=$2
input=$3
output_dir=$4
repeats=${5:-20}
mode=${6:-default}
ours_output=$output_dir/ours.bam
samtools_output=$output_dir/samtools.bam
timing=$output_dir/.timing
times=$output_dir/times.tsv

[[ $(uname -s) == Darwin ]]
[[ $repeats =~ ^[1-9][0-9]*$ ]]
((repeats >= 2))
[[ $mode == default || $mode == equal-workers ]]
mkdir -p "$output_dir"

if [[ $mode == default ]]; then
    ours_threads=automatic
    samtools_threads=0
else
    ours_threads=4
    samtools_threads=4
fi

{
    date -u '+%Y-%m-%dT%H:%M:%SZ'
    printf 'rsomics_commit=%s\n' "${RSOMICS_COMMIT:-unknown}"
    uname -a
    sw_vers
    sysctl -n hw.model machdep.cpu.brand_string hw.memsize
    "$ours" --version
    "$samtools" --version
    python3 --version
    shasum -a 256 "$ours" "$samtools" "$input"
    stat -L -f 'input_bytes=%z' "$input"
    printf 'input_records=%s\n' "$("$samtools" view -c "$input")"
    printf 'mode=%s\nrepeats=%s\n' "$mode" "$repeats"
    printf 'rsomics_memory=128M\n'
    printf 'rsomics_additional_threads=%s\n' "$ours_threads"
    printf 'samtools_additional_threads=%s\n' "$samtools_threads"
} > "$output_dir/environment.txt"
printf 'round\torder\ttool\treal_s\tuser_s\tsystem_s\tmax_rss_bytes\n' > "$times"

record_timing() {
    local round=$1
    local order=$2
    local tool=$3
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$round" "$order" "$tool" \
        "$(awk '$1 == "real" { print $2 }' "$timing")" \
        "$(awk '$1 == "user" { print $2 }' "$timing")" \
        "$(awk '$1 == "sys" { print $2 }' "$timing")" \
        "$(awk '/maximum resident set size/ { print $1 }' "$timing")" >> "$times"
}

run_ours() {
    local round=$1
    local run_order=$2
    local command=("$ours" collate --no-PG -m 128M)
    if [[ $mode == equal-workers ]]; then
        command+=(-@ 4)
    fi
    command+=("$input" -o "$ours_output")
    /usr/bin/time -p -l -o "$timing" "${command[@]}"
    record_timing "$round" "$run_order" rsomics
}

run_samtools() {
    local round=$1
    local run_order=$2
    local command=("$samtools" collate --no-PG -@ "$samtools_threads" -o "$samtools_output" "$input")
    /usr/bin/time -p -l -o "$timing" "${command[@]}"
    record_timing "$round" "$run_order" samtools
}

validate_pair() {
    "$samtools" quickcheck "$ours_output" "$samtools_output"
    "$samtools" view -H --no-PG "$ours_output" > "$output_dir/ours.header"
    "$samtools" view -H --no-PG "$samtools_output" > "$output_dir/samtools.header"
    cmp "$output_dir/ours.header" "$output_dir/samtools.header"
    record_fingerprint "$ours_output" > "$output_dir/ours.checksum"
    record_fingerprint "$samtools_output" > "$output_dir/samtools.checksum"
    cmp "$output_dir/ours.checksum" "$output_dir/samtools.checksum"
}

record_fingerprint() {
    "$samtools" view "$1" | python3 -c '
import hashlib
import sys

modulus = 1 << 256
count = 0
total = 0
squares = 0
xor = 0
for line in sys.stdin.buffer:
    fields = line.rstrip(b"\n").split(b"\t")
    fields[11:] = sorted(fields[11:])
    digest = int.from_bytes(hashlib.sha256(b"\t".join(fields)).digest(), "big")
    count += 1
    total = (total + digest) % modulus
    squares = (squares + digest * digest) % modulus
    xor ^= digest
print(f"{count}\t{total:064x}\t{squares:064x}\t{xor:064x}")
'
}

if [[ $mode == default ]]; then
    "$ours" --json collate --no-PG -m 128M "$input" -o "$ours_output" \
        > "$output_dir/rsomics-summary.json"
else
    "$ours" --json collate --no-PG -m 128M -@ 4 "$input" -o "$ours_output" \
        > "$output_dir/rsomics-summary.json"
fi
run_samtools 0 warmup
validate_pair
printf 'round\torder\ttool\treal_s\tuser_s\tsystem_s\tmax_rss_bytes\n' > "$times"

for ((round = 1; round <= repeats; round++)); do
    if ((round % 2 == 1)); then
        run_order=rsomics-samtools
        run_ours "$round" "$run_order"
        run_samtools "$round" "$run_order"
    else
        run_order=samtools-rsomics
        run_samtools "$round" "$run_order"
        run_ours "$round" "$run_order"
    fi
    validate_pair
done

awk -F '\t' -v repeats="$repeats" '
    NR > 1 {
        count[$3]++
        real[$3] += $4
        user[$3] += $5
        system_time[$3] += $6
        rss[$3] += $7
        paired[$1, $3] = $4
    }
    END {
        for (round = 1; round <= repeats; round++) {
            difference = paired[round, "samtools"] - paired[round, "rsomics"]
            sum += difference
            sum_squares += difference * difference
        }
        mean = sum / repeats
        deviation = sqrt((sum_squares - repeats * mean * mean) / (repeats - 1))
        statistic = deviation > 0 ? mean / (deviation / sqrt(repeats)) : 0
        for (tool in count) {
            printf "%s\t%.6f\t%.6f\t%.6f\t%.1f\n", tool,
                real[tool] / count[tool], user[tool] / count[tool],
                system_time[tool] / count[tool], rss[tool] / count[tool]
        }
        printf "paired_samtools_minus_rsomics\t%.6f\t%.6f\t%.6f\n",
            mean, deviation, statistic
    }
' "$times" > "$output_dir/summary.tsv"

shasum -a 256 \
    "$output_dir/environment.txt" \
    "$output_dir/rsomics-summary.json" \
    "$times" \
    "$output_dir/summary.tsv" \
    "$output_dir/ours.header" \
    "$output_dir/ours.checksum" > "$output_dir/artifacts.sha256"
rm -f "$timing"
