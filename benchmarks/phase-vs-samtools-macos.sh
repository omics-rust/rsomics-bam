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
repeats=${5:-12}
mode=${6:-phase}
ours_output=$output_dir/rsomics.phase
samtools_output=$output_dir/samtools.phase
timing=$output_dir/.timing
times=$output_dir/times.tsv

die() {
    echo "$1" >&2
    exit 2
}

[[ $(uname -s) == Darwin ]] || die "this benchmark requires macOS"
[[ -x $ours ]] || die "rsomics executable is not executable: $ours"
[[ -x $samtools ]] || die "samtools executable is not executable: $samtools"
[[ -f $input ]] || die "input is not a regular file: $input"
[[ $repeats =~ ^[1-9][0-9]*$ ]] || die "REPEATS must be a positive integer"
((repeats >= 2)) || die "REPEATS must be at least 2"
[[ $mode == phase || $mode == scan ]] || die "MODE must be phase or scan"
for artifact in .timing artifacts.sha256 environment.txt rsomics.checksum rsomics.phase \
    samtools.checksum samtools.phase summary.tsv times.tsv; do
    [[ ! -e $output_dir/$artifact ]] || die "output directory contains prior benchmark data: $output_dir"
done
mkdir -p "$output_dir"

if [[ $mode == phase ]]; then
    minimum_lod=37
else
    minimum_lod=1000
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
    printf 'mode=%s\nminimum_lod=%s\nrepeats=%s\n' "$mode" "$minimum_lod" "$repeats"
    "$samtools" view -H "$input"
} > "$output_dir/environment.txt"

fingerprint() {
    python3 -c '
import hashlib
import sys

digest = hashlib.sha256()
evidence = []
phase_sets = 0
markers = 0

def flush():
    for line in sorted(evidence):
        digest.update(line)
    evidence.clear()

for line in sys.stdin.buffer:
    if line.startswith(b"EV\t"):
        evidence.append(line)
        continue
    flush()
    digest.update(line)
    phase_sets += line.startswith(b"PS\t")
    markers += line.startswith((b"M0\t", b"M1\t", b"M2\t"))
flush()
print(f"{phase_sets}\t{markers}\t{digest.hexdigest()}")
' < "$1"
}

validate() {
    fingerprint "$ours_output" > "$output_dir/rsomics.checksum"
    fingerprint "$samtools_output" > "$output_dir/samtools.checksum"
    cmp "$output_dir/rsomics.checksum" "$output_dir/samtools.checksum"
}

"$ours" phase -q "$minimum_lod" "$input" > "$ours_output"
"$samtools" phase -q "$minimum_lod" "$input" > "$samtools_output"
validate
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

measure() {
    local round=$1
    local order=$2
    local tool=$3
    if [[ $tool == rsomics ]]; then
        /usr/bin/time -p -l -o "$timing" "$ours" phase \
            -q "$minimum_lod" "$input" > /dev/null
    else
        /usr/bin/time -p -l -o "$timing" "$samtools" phase \
            -q "$minimum_lod" "$input" > /dev/null
    fi
    record_timing "$round" "$order" "$tool"
}

for ((round = 1; round <= repeats; round++)); do
    if ((round % 2 == 1)); then
        order=rsomics-samtools
        measure "$round" "$order" rsomics
        measure "$round" "$order" samtools
    else
        order=samtools-rsomics
        measure "$round" "$order" samtools
        measure "$round" "$order" rsomics
    fi
done
validate

awk -F '\t' -v repeats="$repeats" '
    function median(tool, n,    i, j, key, value) {
        for (i = 1; i <= n; i++) ordered[i] = real_value[tool, i]
        for (i = 2; i <= n; i++) {
            key = ordered[i]
            j = i - 1
            while (j >= 1 && ordered[j] > key) {
                ordered[j + 1] = ordered[j]
                j--
            }
            ordered[j + 1] = key
        }
        value = n % 2 ? ordered[(n + 1) / 2] : (ordered[n / 2] + ordered[n / 2 + 1]) / 2
        for (i = 1; i <= n; i++) delete ordered[i]
        return value
    }
    NR > 1 {
        count[$3]++
        real_value[$3, count[$3]] = $4
        real[$3] += $4
        real_squared[$3] += $4 * $4
        user[$3] += $5
        system_time[$3] += $6
        rss[$3] += $7
        paired[$1, $3] = $4
    }
    END {
        print "tool\tn\tmean_real_s\tmedian_real_s\tstddev_real_s\tmean_user_s\tmean_system_s\tmean_rss_bytes"
        for (tool in count) {
            mean = real[tool] / count[tool]
            variance = (real_squared[tool] - count[tool] * mean * mean) / (count[tool] - 1)
            if (variance < 0) variance = 0
            printf "%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.1f\n", tool,
                count[tool], mean, median(tool, count[tool]), sqrt(variance),
                user[tool] / count[tool], system_time[tool] / count[tool],
                rss[tool] / count[tool]
        }
        sum = squared = wins = 0
        for (round = 1; round <= repeats; round++) {
            difference = paired[round, "samtools"] - paired[round, "rsomics"]
            sum += difference
            squared += difference * difference
            if (difference > 0) wins++
        }
        mean = sum / repeats
        variance = (squared - repeats * mean * mean) / (repeats - 1)
        if (variance < 0) variance = 0
        deviation = sqrt(variance)
        statistic = deviation > 0 ? mean / (deviation / sqrt(repeats)) : 0
        printf "paired\t%d\t%.6f\t%.6f\t%.6f\t%d\n", repeats, mean,
            deviation, statistic, wins
    }
' "$times" > "$output_dir/summary.tsv"

shasum -a 256 \
    "$output_dir/environment.txt" \
    "$times" \
    "$output_dir/summary.tsv" \
    "$output_dir/rsomics.checksum" \
    "$output_dir/samtools.checksum" > "$output_dir/artifacts.sha256"
rm -f "$timing"
