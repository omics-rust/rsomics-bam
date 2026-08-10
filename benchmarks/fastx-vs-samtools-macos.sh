#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

if [[ $# -ne 5 ]]; then
    echo "usage: $0 RSOMICS_BAM SAMTOOLS INPUT_BAM OUTPUT_DIR REPEATS" >&2
    exit 2
fi

ours=$1
samtools=$2
input=$3
output_dir=$4
repeats=$5
timing=$output_dir/.timing
times=$output_dir/times.tsv

[[ $(uname -s) == Darwin ]]
[[ $repeats =~ ^[1-9][0-9]*$ ]]
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
    shasum -a 256 "$ours" "$samtools" "$input"
    stat -L -f 'input=%N bytes=%z' "$input"
    printf 'records=%s\nrepeats=%s\n' "$($samtools view -c "$input")" "$repeats"
} > "$output_dir/environment.txt"

for format in fasta fastq; do
    "$ours" "$format" "$input" | shasum -a 256 > "$output_dir/rsomics-$format.sha256"
    "$samtools" "$format" "$input" 2>/dev/null | shasum -a 256 \
        > "$output_dir/samtools-$format.sha256"
    cmp "$output_dir/rsomics-$format.sha256" "$output_dir/samtools-$format.sha256"
done

printf 'format\tround\torder\ttool\treal_s\tuser_s\tsystem_s\tmax_rss_bytes\n' > "$times"

record_timing() {
    local format=$1
    local round=$2
    local order=$3
    local tool=$4
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$format" "$round" "$order" "$tool" \
        "$(awk '$1 == "real" { print $2 }' "$timing")" \
        "$(awk '$1 == "user" { print $2 }' "$timing")" \
        "$(awk '$1 == "sys" { print $2 }' "$timing")" \
        "$(awk '/maximum resident set size/ { print $1 }' "$timing")" >> "$times"
}

run_ours() {
    local format=$1
    local round=$2
    local order=$3
    /usr/bin/time -p -l -o "$timing" "$ours" "$format" "$input" >/dev/null
    record_timing "$format" "$round" "$order" rsomics
}

run_samtools() {
    local format=$1
    local round=$2
    local order=$3
    /usr/bin/time -p -l -o "$timing" "$samtools" "$format" "$input" \
        >/dev/null 2>/dev/null
    record_timing "$format" "$round" "$order" samtools
}

for format in fasta fastq; do
    run_ours "$format" 0 warmup rsomics
    run_samtools "$format" 0 warmup samtools
    for ((round = 1; round <= repeats; round++)); do
        if ((round % 2 == 1)); then
            order=rsomics-samtools
            run_ours "$format" "$round" "$order"
            run_samtools "$format" "$round" "$order"
        else
            order=samtools-rsomics
            run_samtools "$format" "$round" "$order"
            run_ours "$format" "$round" "$order"
        fi
    done
done

awk -F '\t' '
    NR > 1 && $2 != 0 {
        key = $1 SUBSEP $4
        count[key]++
        real[key] += $5
        real_squared[key] += $5 * $5
        user[key] += $6
        system_time[key] += $7
        rss[key] += $8
    }
    END {
        print "format\ttool\tmean_real_s\tstddev_real_s\tmean_user_s\tmean_system_s\tmean_rss_bytes"
        for (key in count) {
            split(key, fields, SUBSEP)
            mean = real[key] / count[key]
            variance = count[key] > 1 \
                ? (real_squared[key] - count[key] * mean * mean) / (count[key] - 1) \
                : 0
            printf "%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.1f\n", \
                fields[1], fields[2], mean, sqrt(variance), \
                user[key] / count[key], system_time[key] / count[key], \
                rss[key] / count[key]
        }
    }
' "$times" > "$output_dir/summary.tsv"

shasum -a 256 \
    "$output_dir/environment.txt" \
    "$times" \
    "$output_dir/summary.tsv" \
    "$output_dir/rsomics-fasta.sha256" \
    "$output_dir/rsomics-fastq.sha256" > "$output_dir/artifacts.sha256"
rm -f "$timing"
