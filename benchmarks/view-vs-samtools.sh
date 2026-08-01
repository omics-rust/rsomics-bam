#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

if [[ $# -lt 4 || $# -gt 7 ]]; then
    echo "usage: $0 RSOMICS_BAM SAMTOOLS INPUT OUTPUT_DIR [REPEATS] [THREADS] [CPUSET]" >&2
    exit 2
fi

ours=$1
samtools=$2
input=$3
output_dir=$4
repeats=${5:-5}
threads=${6:-4}
cpuset=${7:-0-4}
ours_output=$output_dir/ours.bam
samtools_output=$output_dir/samtools.bam
times=$output_dir/times.tsv
checksums=$output_dir/checksums.tsv

[[ $(uname -s) == Linux ]]
command -v taskset >/dev/null
/usr/bin/time --version 2>&1 | grep -q GNU
[[ $repeats =~ ^[1-9][0-9]*$ ]]
[[ $threads =~ ^[1-9][0-9]*$ ]]
[[ $cpuset =~ ^[0-9,-]+$ ]]
mkdir -p "$output_dir"

{
    date --utc --iso-8601=seconds
    printf 'rsomics_commit=%s\n' "${RSOMICS_COMMIT:-unknown}"
    uname -a
    uptime
    lscpu
    "$ours" --version
    "$samtools" --version
    sha256sum "$ours" "$samtools" "$input"
    stat --format='input_bytes=%s' "$input"
    printf 'input_records=%s\n' "$("$samtools" view -c "$input")"
    printf 'threads=%s\ncpuset=%s\nrepeats=%s\n' "$threads" "$cpuset" "$repeats"
} > "$output_dir/environment.txt"
printf 'round\ttool\treal_s\tuser_s\tsystem_s\tmax_rss_kib\n' > "$times"
printf 'round\theader_sha256\trecord_sha256\trecords\n' > "$checksums"

run_ours() {
    local round=$1
    /usr/bin/time -a -o "$times" -f "$round\tours\t%e\t%U\t%S\t%M" \
        taskset -c "$cpuset" "$ours" view -b -@ "$threads" --no-pg \
        -o "$ours_output" "$input"
}

run_samtools() {
    local round=$1
    /usr/bin/time -a -o "$times" -f "$round\tsamtools\t%e\t%U\t%S\t%M" \
        taskset -c "$cpuset" "$samtools" view -b -@ "$threads" --no-PG \
        -o "$samtools_output" "$input"
}

validate_pair() {
    local round=$1
    "$samtools" quickcheck "$ours_output" "$samtools_output"

    local ours_header samtools_header ours_records samtools_records ours_count samtools_count
    ours_header=$("$samtools" view -H --no-PG "$ours_output" | sha256sum | cut -d' ' -f1)
    samtools_header=$("$samtools" view -H --no-PG "$samtools_output" | sha256sum | cut -d' ' -f1)
    ours_records=$("$samtools" view "$ours_output" | sha256sum | cut -d' ' -f1)
    samtools_records=$("$samtools" view "$samtools_output" | sha256sum | cut -d' ' -f1)
    ours_count=$("$samtools" view -c "$ours_output")
    samtools_count=$("$samtools" view -c "$samtools_output")

    [[ $ours_header == "$samtools_header" ]]
    [[ $ours_records == "$samtools_records" ]]
    [[ $ours_count == "$samtools_count" ]]
    printf '%s\t%s\t%s\t%s\n' "$round" "$ours_header" "$ours_records" "$ours_count" \
        >> "$checksums"
}

taskset -c "$cpuset" "$ours" view -b -@ "$threads" --no-pg \
    -o "$ours_output" "$input"
taskset -c "$cpuset" "$samtools" view -b -@ "$threads" --no-PG \
    -o "$samtools_output" "$input"
validate_pair warmup

for ((round = 1; round <= repeats; round++)); do
    if ((round % 2 == 1)); then
        run_ours "$round"
        run_samtools "$round"
    else
        run_samtools "$round"
        run_ours "$round"
    fi
    validate_pair "$round"
done
