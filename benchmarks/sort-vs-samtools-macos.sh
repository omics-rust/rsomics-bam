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
ours_prefix=$output_dir/ours-tmp
samtools_prefix=$output_dir/samtools-tmp
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
    samtools_memory=768M
else
    ours_threads=4
    samtools_threads=4
    samtools_memory=192M
fi

{
    date -u '+%Y-%m-%dT%H:%M:%SZ'
    printf 'rsomics_commit=%s\n' "${RSOMICS_COMMIT:-unknown}"
    uname -a
    sw_vers
    sysctl -n hw.model machdep.cpu.brand_string hw.memsize
    "$ours" --version
    "$samtools" --version
    shasum -a 256 "$ours" "$samtools" "$input"
    stat -L -f 'input_bytes=%z' "$input"
    printf 'input_records=%s\n' "$("$samtools" view -c "$input")"
    printf 'mode=%s\nrepeats=%s\n' "$mode" "$repeats"
    printf 'rsomics_memory=768M\nrsomics_additional_threads=%s\n' "$ours_threads"
    printf 'samtools_memory_per_thread=%s\nsamtools_additional_threads=%s\n' \
        "$samtools_memory" "$samtools_threads"
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
    if [[ $mode == default ]]; then
        /usr/bin/time -p -l -o "$timing" "$ours" sort --no-PG -m 768M \
            -T "$ours_prefix" -o "$ours_output" "$input"
    else
        /usr/bin/time -p -l -o "$timing" "$ours" sort --no-PG -@ 4 -m 768M \
            -T "$ours_prefix" -o "$ours_output" "$input"
    fi
    record_timing "$round" "$order" rsomics
}

run_samtools() {
    local round=$1
    local order=$2
    /usr/bin/time -p -l -o "$timing" "$samtools" sort --no-PG \
        -@ "$samtools_threads" -m "$samtools_memory" -T "$samtools_prefix" \
        -o "$samtools_output" "$input"
    record_timing "$round" "$order" samtools
}

validate_pair() {
    "$samtools" quickcheck "$ours_output" "$samtools_output"
    "$samtools" view -H --no-PG "$ours_output" > "$output_dir/ours.header"
    "$samtools" view -H --no-PG "$samtools_output" > "$output_dir/samtools.header"
    cmp "$output_dir/ours.header" "$output_dir/samtools.header"
    "$samtools" checksum -a -O -T "$ours_output" \
        | awk '!/^#/' > "$output_dir/ours.checksum"
    "$samtools" checksum -a -O -T "$samtools_output" \
        | awk '!/^#/' > "$output_dir/samtools.checksum"
    cmp "$output_dir/ours.checksum" "$output_dir/samtools.checksum"
}

if [[ $mode == default ]]; then
    "$ours" --json sort --no-PG -m 768M -T "$ours_prefix" \
        -o "$ours_output" "$input" > "$output_dir/rsomics-summary.json"
else
    "$ours" --json sort --no-PG -@ 4 -m 768M -T "$ours_prefix" \
        -o "$ours_output" "$input" > "$output_dir/rsomics-summary.json"
fi
"$samtools" sort --no-PG -@ "$samtools_threads" -m "$samtools_memory" \
    -T "$samtools_prefix" -o "$samtools_output" "$input"
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
    "$output_dir/ours.header" \
    "$output_dir/ours.checksum" > "$output_dir/artifacts.sha256"
rm -f "$timing"
