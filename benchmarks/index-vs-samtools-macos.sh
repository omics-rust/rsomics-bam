#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

if [[ $# -lt 4 || $# -gt 5 ]]; then
    echo "usage: $0 RSOMICS_BAM SAMTOOLS INPUT OUTPUT_DIR [REPEATS]" >&2
    exit 2
fi

ours=$1
samtools=$2
input=$3
output_dir=$4
repeats=${5:-20}
ours_index=$output_dir/ours.bai
samtools_index=$output_dir/samtools.bai
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
    stat -f 'input_bytes=%z' "$input"
    printf 'input_records=%s\n' "$("$samtools" view -c "$input")"
    printf 'repeats=%s\n' "$repeats"
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
    /usr/bin/time -p -l -o "$timing" "$ours" index -o "$ours_index" "$input"
    record_timing "$round" "$order" rsomics
}

run_samtools() {
    local round=$1
    local order=$2
    /usr/bin/time -p -l -o "$timing" "$samtools" index -o "$samtools_index" "$input"
    record_timing "$round" "$order" samtools
}

validate_pair() {
    cmp "$ours_index" "$samtools_index"
    cmp \
        <("$samtools" idxstats -X "$input" "$ours_index") \
        <("$samtools" idxstats -X "$input" "$samtools_index")
}

"$ours" index -o "$ours_index" "$input"
"$samtools" index -o "$samtools_index" "$input"
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

shasum -a 256 "$ours_index" "$samtools_index" > "$output_dir/indexes.sha256"
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
        for (tool in count) {
            printf "%s\t%.6f\t%.6f\t%.6f\t%.1f\n", tool,
                real[tool] / count[tool], user[tool] / count[tool],
                system_time[tool] / count[tool], rss[tool] / count[tool]
        }
        printf "paired\t%.6f\t%.6f\t%.6f\t%d/%d\n", mean, deviation,
            mean / (deviation / sqrt(repeats)), wins, repeats
    }
' "$times" > "$output_dir/summary.tsv"

rm -f "$timing"
