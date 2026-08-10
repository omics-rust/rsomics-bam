#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

if [[ $# -lt 7 ]]; then
    echo "usage: $0 RSOMICS_BAM SAMTOOLS OPERATION OUTPUT_DIR REPEATS INPUT..." >&2
    exit 2
fi

ours=$1
samtools=$2
operation=$3
output_dir=$4
repeats=$5
shift 5
inputs=("$@")
ours_output=$output_dir/ours.bam
samtools_output=$output_dir/samtools.bam
timing=$output_dir/.timing
times=$output_dir/times.tsv

[[ $(uname -s) == Darwin ]]
[[ $repeats =~ ^[1-9][0-9]*$ ]]
((repeats >= 2))
[[ $operation == cat || $operation == reheader ]]
if [[ $operation == cat ]]; then
    ((${#inputs[@]} >= 2))
else
    ((${#inputs[@]} == 2))
fi
mkdir -p "$output_dir"

{
    date -u '+%Y-%m-%dT%H:%M:%SZ'
    printf 'rsomics_commit=%s\n' "${RSOMICS_COMMIT:-unknown}"
    uname -a
    sw_vers
    sysctl -n hw.model machdep.cpu.brand_string hw.memsize
    "$ours" --version
    "$samtools" --version
    shasum -a 256 "$ours" "$samtools" "${inputs[@]}"
    for input in "${inputs[@]}"; do
        stat -L -f 'input=%N bytes=%z' "$input"
    done
    printf 'operation=%s\nrepeats=%s\n' "$operation" "$repeats"
    if [[ $operation == cat ]]; then
        for input in "${inputs[@]}"; do
            printf 'records=%s input=%s\n' "$("$samtools" view -c "$input")" "$input"
        done
    else
        printf 'records=%s input=%s\n' \
            "$("$samtools" view -c "${inputs[1]}")" "${inputs[1]}"
    fi
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
    local order=$2
    /usr/bin/time -p -l -o "$timing" \
        "$ours" "$operation" --no-pg "${inputs[@]}" -o "$ours_output"
    record_timing "$round" "$order" rsomics
}

run_samtools() {
    local round=$1
    local order=$2
    if [[ $operation == cat ]]; then
        /usr/bin/time -p -l -o "$timing" \
            "$samtools" cat --no-PG -o "$samtools_output" "${inputs[@]}"
    else
        /usr/bin/time -p -l -o "$timing" \
            "$samtools" reheader --no-PG "${inputs[@]}" > "$samtools_output"
    fi
    record_timing "$round" "$order" samtools
}

validate_pair() {
    "$samtools" quickcheck "$ours_output" "$samtools_output"
    "$samtools" view -h --no-PG "$ours_output" | shasum -a 256 > "$output_dir/ours.checksum"
    "$samtools" view -h --no-PG "$samtools_output" | shasum -a 256 > "$output_dir/samtools.checksum"
    cmp "$output_dir/ours.checksum" "$output_dir/samtools.checksum"
}

"$ours" --json "$operation" --no-pg "${inputs[@]}" -o "$ours_output" \
    > "$output_dir/rsomics-summary.json"
if [[ $operation == cat ]]; then
    "$samtools" cat --no-PG -o "$samtools_output" "${inputs[@]}"
else
    "$samtools" reheader --no-PG "${inputs[@]}" > "$samtools_output"
fi
validate_pair

for ((round = 1; round <= repeats; round++)); do
    if ((round % 2 == 1)); then
        order=rsomics-samtools
        run_ours "$round" "$order"
        run_samtools "$round" "$order"
    else
        order=samtools-rsomics
        run_samtools "$round" "$order"
        run_ours "$round" "$order"
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
            difference = paired[round, "rsomics"] - paired[round, "samtools"]
            sum += difference
            sum_squares += difference * difference
            if (difference < 0) wins++
        }
        mean = sum / repeats
        deviation = sqrt((sum_squares - repeats * mean * mean) / (repeats - 1))
        statistic = deviation > 0 ? mean / (deviation / sqrt(repeats)) : 0
        for (tool in count) {
            printf "%s\t%.6f\t%.6f\t%.6f\t%.1f\n", tool,
                real[tool] / count[tool], user[tool] / count[tool],
                system_time[tool] / count[tool], rss[tool] / count[tool]
        }
        printf "paired\t%.6f\t%.6f\t%.6f\t%d/%d\n", mean, deviation,
            statistic, wins, repeats
    }
' "$times" > "$output_dir/summary.tsv"

shasum -a 256 \
    "$output_dir/environment.txt" \
    "$times" \
    "$output_dir/summary.tsv" \
    "$output_dir/rsomics-summary.json" \
    "$output_dir/ours.checksum" > "$output_dir/artifacts.sha256"
rm -f "$timing"
