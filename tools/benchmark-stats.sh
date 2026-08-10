#!/bin/sh
set -eu

if [ "$#" -ne 6 ]; then
    printf 'usage: %s RSOMICS_BAM SAMTOOLS INPUT OUTPUT_DIR ROUNDS THREADS\n' "$0" >&2
    exit 2
fi

ours=$1
oracle=$2
input=$3
output=$4
rounds=$5
threads=$6

mkdir -p "$output"
ledger="$output/timings.tsv"
summary="$output/summary.tsv"
environment="$output/environment.txt"

"$ours" stats -@ "$threads" "$input" > "$output/rsomics.report"
"$oracle" stats -@ "$threads" "$input" > "$output/samtools.report"
tail -n +4 "$output/rsomics.report" > "$output/rsomics.stable"
tail -n +4 "$output/samtools.report" > "$output/samtools.stable"
cmp "$output/rsomics.stable" "$output/samtools.stable"

{
    uname -a
    sw_vers
    "$ours" --version
    "$oracle" --version
    printf 'threads %s\n' "$threads"
    printf 'rounds %s\n' "$rounds"
    printf 'input_bytes %s\n' "$(stat -f %z "$input")"
    printf 'input_records %s\n' "$("$oracle" view -c "$input")"
    shasum -a 256 "$ours" "$oracle" "$input" "$output/rsomics.stable" "$output/samtools.stable"
} > "$environment"

printf 'tool\tround\treal_seconds\tuser_seconds\tsystem_seconds\tmax_rss_bytes\n' > "$ledger"

measure() {
    tool=$1
    round=$2
    shift 2
    timing=$(mktemp "$output/time.XXXXXX")
    /usr/bin/time -lp "$@" > /dev/null 2> "$timing"
    real=$(awk '$1 == "real" { print $2 }' "$timing")
    user=$(awk '$1 == "user" { print $2 }' "$timing")
    system=$(awk '$1 == "sys" { print $2 }' "$timing")
    rss=$(awk '$2 == "maximum" && $3 == "resident" { print $1 }' "$timing")
    rm "$timing"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$tool" "$round" "$real" "$user" "$system" "$rss" >> "$ledger"
}

"$ours" stats -@ "$threads" "$input" > /dev/null
"$oracle" stats -@ "$threads" "$input" > /dev/null

round=1
while [ "$round" -le "$rounds" ]; do
    if [ $((round % 2)) -eq 1 ]; then
        measure rsomics "$round" "$ours" stats -@ "$threads" "$input"
        measure samtools "$round" "$oracle" stats -@ "$threads" "$input"
    else
        measure samtools "$round" "$oracle" stats -@ "$threads" "$input"
        measure rsomics "$round" "$ours" stats -@ "$threads" "$input"
    fi
    round=$((round + 1))
done

printf 'tool\tmean_wall_seconds\tmean_cpu_seconds\tmean_max_rss_bytes\tmedian_wall_seconds\tmedian_max_rss_bytes\n' > "$summary"
for tool in rsomics samtools; do
    means=$(awk -F '\t' -v tool="$tool" '$1 == tool { wall += $3; cpu += $4 + $5; rss += $6; n += 1 } END { printf "%.6f\t%.6f\t%.0f", wall / n, cpu / n, rss / n }' "$ledger")
    median_wall=$(awk -F '\t' -v tool="$tool" '$1 == tool { print $3 }' "$ledger" | sort -n | awk '{ value[NR] = $1 } END { if (NR % 2) print value[(NR + 1) / 2]; else printf "%.6f", (value[NR / 2] + value[NR / 2 + 1]) / 2 }')
    median_rss=$(awk -F '\t' -v tool="$tool" '$1 == tool { print $6 }' "$ledger" | sort -n | awk '{ value[NR] = $1 } END { if (NR % 2) print value[(NR + 1) / 2]; else printf "%.0f", (value[NR / 2] + value[NR / 2 + 1]) / 2 }')
    printf '%s\t%s\t%s\t%s\n' "$tool" "$means" "$median_wall" "$median_rss" >> "$summary"
done
