#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

if [[ $# -ne 7 ]]; then
    echo "usage: $0 RSOMICS_BAM SAMTOOLS INPUT OUTPUT_DIR REPEATS THREADS RG_LINE" >&2
    exit 2
fi

ours=$1
samtools=$2
input=$3
output_dir=$4
repeats=$5
threads=$6
rg_line=$7
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
    shasum -a 256 "$ours" "$samtools" "$input"
    stat -L -f 'input=%N bytes=%z' "$input"
    printf 'input_records=%s\n' "$("$samtools" view -c "$input")"
    printf 'input_records_with_rg=%s\n' \
        "$("$samtools" view -d RG "$input" | wc -l | tr -d ' ')"
    printf 'repeats=%s\nthreads=%s\nrg_line=%s\n' "$repeats" "$threads" "$rg_line"
} > "$output_dir/environment.txt"

run_edit() {
    local tool=$1
    local context=$2
    local destination=$3
    if [[ $tool == rsomics ]]; then
        "$ours" addreplacerg -r "$rg_line" -m "$context" --no-PG \
            -@ "$threads" -o "$destination" "$input"
    else
        "$samtools" addreplacerg -r "$rg_line" -m "$context" --no-PG \
            -@ "$threads" -o "$destination" "$input"
    fi
}

fingerprint() {
    local context=$1
    local tool=$2
    local bam=$3
    local raw_header=$output_dir/.$context-$tool.header
    "$samtools" view -c "$bam" >/dev/null
    "$samtools" view -H "$bam" > "$raw_header"
    {
        grep -E '^@HD[[:space:]]' "$raw_header"
        grep -E '^@SQ[[:space:]]' "$raw_header"
        grep -E '^@RG[[:space:]]' "$raw_header"
        grep -Ev '^@(HD|SQ|RG|PG)[[:space:]]' "$raw_header" | sort
    } | shasum -a 256 > "$output_dir/$context-$tool-header.sha256"
    "$samtools" view "$bam" | shasum -a 256 \
        > "$output_dir/$context-$tool-records.sha256"
    rm -f "$raw_header"
}

validate() {
    local context=$1
    fingerprint "$context" rsomics "$output_dir/$context-rsomics.bam"
    fingerprint "$context" samtools "$output_dir/$context-samtools.bam"
    cmp "$output_dir/$context-rsomics-header.sha256" \
        "$output_dir/$context-samtools-header.sha256"
    cmp "$output_dir/$context-rsomics-records.sha256" \
        "$output_dir/$context-samtools-records.sha256"
}

for context in overwrite_all orphan_only; do
    run_edit rsomics "$context" "$output_dir/$context-rsomics.bam"
    run_edit samtools "$context" "$output_dir/$context-samtools.bam"
    validate "$context"
done

printf 'context\tround\torder\ttool\treal_s\tuser_s\tsystem_s\tmax_rss_bytes\n' > "$times"

record_timing() {
    local context=$1
    local round=$2
    local order=$3
    local tool=$4
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$context" "$round" "$order" "$tool" \
        "$(awk '$1 == "real" { print $2 }' "$timing")" \
        "$(awk '$1 == "user" { print $2 }' "$timing")" \
        "$(awk '$1 == "sys" { print $2 }' "$timing")" \
        "$(awk '/maximum resident set size/ { print $1 }' "$timing")" >> "$times"
}

measure() {
    local context=$1
    local round=$2
    local order=$3
    local tool=$4
    local destination=$output_dir/$context-$tool.bam
    if [[ $tool == rsomics ]]; then
        /usr/bin/time -p -l -o "$timing" "$ours" addreplacerg \
            -r "$rg_line" -m "$context" --no-PG -@ "$threads" \
            -o "$destination" "$input"
    else
        /usr/bin/time -p -l -o "$timing" "$samtools" addreplacerg \
            -r "$rg_line" -m "$context" --no-PG -@ "$threads" \
            -o "$destination" "$input"
    fi
    "$samtools" view -c "$destination" >/dev/null
    record_timing "$context" "$round" "$order" "$tool"
}

for context in overwrite_all orphan_only; do
    measure "$context" 0 warmup rsomics
    measure "$context" 0 warmup samtools
    for ((round = 1; round <= repeats; round++)); do
        if ((round % 2 == 1)); then
            order=rsomics-samtools
            measure "$context" "$round" "$order" rsomics
            measure "$context" "$round" "$order" samtools
        else
            order=samtools-rsomics
            measure "$context" "$round" "$order" samtools
            measure "$context" "$round" "$order" rsomics
        fi
    done
    validate "$context"
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
        paired[$1 SUBSEP $2 SUBSEP $4] = $5
        contexts[$1] = 1
    }
    END {
        print "context\ttool\tmean_real_s\tstddev_real_s\tmean_user_s\tmean_system_s\tmean_rss_bytes"
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

awk -F '\t' '
    NR > 1 && $2 != 0 {
        key = $1 SUBSEP $2
        if ($4 == "rsomics") ours[key] = $5
        else oracle[key] = $5
        contexts[$1] = 1
    }
    END {
        print "context\tpairs\tmean_rsomics_minus_samtools_s\tstddev_s\tt_statistic\trsomics_wins"
        for (context in contexts) {
            n = sum = squared = wins = 0
            for (key in ours) {
                split(key, fields, SUBSEP)
                if (fields[1] != context) continue
                difference = ours[key] - oracle[key]
                sum += difference
                squared += difference * difference
                if (difference < 0) wins++
                n++
            }
            mean = sum / n
            variance = n > 1 ? (squared - n * mean * mean) / (n - 1) : 0
            stddev = sqrt(variance)
            statistic = stddev > 0 ? mean / (stddev / sqrt(n)) : 0
            printf "%s\t%d\t%.6f\t%.6f\t%.6f\t%d\n", \
                context, n, mean, stddev, statistic, wins
        }
    }
' "$times" > "$output_dir/paired.tsv"

shasum -a 256 \
    "$output_dir/environment.txt" \
    "$times" \
    "$output_dir/summary.tsv" \
    "$output_dir/paired.tsv" \
    "$output_dir/overwrite_all-rsomics-header.sha256" \
    "$output_dir/overwrite_all-rsomics-records.sha256" \
    "$output_dir/orphan_only-rsomics-header.sha256" \
    "$output_dir/orphan_only-rsomics-records.sha256" > "$output_dir/artifacts.sha256"
rm -f "$timing"
