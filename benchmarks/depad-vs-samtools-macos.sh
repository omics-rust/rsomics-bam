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
    "$ours" depad --no-pg -@ "$threads" -T "$reference" -o "$ours_output" "$input"
}

run_samtools() {
    "$samtools" depad --no-PG --threads "$threads" -T "$reference" \
        -o "$samtools_output" "$input"
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
        /usr/bin/time -p -l -o "$timing" "$ours" depad --no-pg \
            -@ "$threads" -T "$reference" -o "$ours_output" "$input"
    else
        /usr/bin/time -p -l -o "$timing" "$samtools" depad --no-PG \
            --threads "$threads" -T "$reference" -o "$samtools_output" "$input"
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
    function median(tool, n,    i, j, key) {
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
    "$output_dir/rsomics-view.sha256" > "$output_dir/artifacts.sha256"
rm -f "$timing"
