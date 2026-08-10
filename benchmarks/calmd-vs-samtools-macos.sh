#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

if [[ $# -lt 5 || $# -gt 7 ]]; then
    echo "usage: $0 RSOMICS_BAM SAMTOOLS INPUT REFERENCE OUTPUT_DIR [REPEATS] [THREADS]" >&2
    exit 2
fi

ours=$1
samtools=$2
input=$3
reference=$4
output_dir=$5
repeats=${6:-20}
threads=${7:-4}
ours_output=$output_dir/rsomics.bam
samtools_output=$output_dir/samtools.bam
timing=$output_dir/.timing
times=$output_dir/times.tsv

[[ $(uname -s) == Darwin ]]
[[ $repeats =~ ^[1-9][0-9]*$ ]]
[[ $threads =~ ^[1-9][0-9]*$ ]]
((repeats >= 2))
mkdir -p "$output_dir"

{
    date -u '+%Y-%m-%dT%H:%M:%SZ'
    printf 'rsomics_commit=%s\n' "${RSOMICS_COMMIT:-unknown}"
    uname -a
    sw_vers
    sysctl -n hw.model machdep.cpu.brand_string hw.memsize
    "$ours" --version
    "$samtools" --version
    shasum -a 256 "$ours" "$samtools" "$input" "$reference"
    stat -L -f 'input=%N bytes=%z' "$input"
    stat -L -f 'reference=%N bytes=%z' "$reference"
    printf 'input_records=%s\n' "$("$samtools" view -c "$input")"
    printf 'repeats=%s\nthreads=%s\n' "$repeats" "$threads"
} > "$output_dir/environment.txt"

run_ours() {
    "$ours" calmd --no-PG -b -@ "$threads" "$input" "$reference" > "$ours_output"
}

run_samtools() {
    "$samtools" calmd --no-PG -b -@ "$threads" "$input" "$reference" \
        > "$samtools_output"
}

validate() {
    "$samtools" quickcheck -v "$ours_output" "$samtools_output"
    "$samtools" view --no-PG -h "$ours_output" | shasum -a 256 \
        > "$output_dir/rsomics-view.sha256"
    "$samtools" view --no-PG -h "$samtools_output" | shasum -a 256 \
        > "$output_dir/samtools-view.sha256"
    cmp "$output_dir/rsomics-view.sha256" "$output_dir/samtools-view.sha256"
}

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
        /usr/bin/time -p -l -o "$timing" "$ours" calmd --no-PG -b \
            -@ "$threads" "$input" "$reference" > "$ours_output"
    else
        /usr/bin/time -p -l -o "$timing" "$samtools" calmd --no-PG -b \
            -@ "$threads" "$input" "$reference" > "$samtools_output"
    fi
    record_timing "$round" "$order" "$tool"
}

run_ours
run_samtools
validate
printf 'round\torder\ttool\treal_s\tuser_s\tsystem_s\tmax_rss_bytes\n' > "$times"

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
    NR > 1 {
        count[$3]++
        real[$3] += $4
        real_squared[$3] += $4 * $4
        user[$3] += $5
        system_time[$3] += $6
        rss[$3] += $7
        paired[$1, $3] = $4
    }
    END {
        print "tool\tn\tmean_real_s\tstddev_real_s\tmean_user_s\tmean_system_s\tmean_rss_bytes"
        for (tool in count) {
            mean = real[tool] / count[tool]
            variance = (real_squared[tool] - count[tool] * mean * mean) / (count[tool] - 1)
            if (variance < 0) variance = 0
            printf "%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.1f\n", tool,
                count[tool], mean, sqrt(variance), user[tool] / count[tool],
                system_time[tool] / count[tool], rss[tool] / count[tool]
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
    "$output_dir/rsomics-view.sha256" > "$output_dir/artifacts.sha256"
rm -f "$timing"
