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

if [ -e "$output" ]; then
    printf 'output directory already exists: %s\n' "$output" >&2
    exit 2
fi
mkdir -p "$output"

version=$($oracle --version | sed -n '1p')
if [ "$version" != 'samtools 1.24' ]; then
    printf 'expected samtools 1.24, found %s\n' "$version" >&2
    exit 2
fi

ledger="$output/timings.tsv"
summary="$output/summary.tsv"
environment="$output/environment.txt"

"$ours" checksum -@ "$threads" "$input" > "$output/rsomics.report"
"$oracle" checksum -@ "$threads" "$input" > "$output/samtools.report"
cmp "$output/rsomics.report" "$output/samtools.report"

{
    uname -a
    if command -v sw_vers >/dev/null 2>&1; then
        sw_vers
    fi
    "$ours" --version
    "$oracle" --version
    printf 'additional_bam_workers %s\n' "$threads"
    printf 'rounds %s\n' "$rounds"
    printf 'input_bytes %s\n' "$(stat -f %z "$input")"
    printf 'input_records %s\n' "$("$oracle" view -c "$input")"
    shasum -a 256 "$ours" "$oracle" "$input" "$0" "$output/rsomics.report" "$output/samtools.report"
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

"$ours" checksum -@ "$threads" "$input" > /dev/null
"$oracle" checksum -@ "$threads" "$input" > /dev/null

round=1
while [ "$round" -le "$rounds" ]; do
    if [ $((round % 2)) -eq 1 ]; then
        measure rsomics "$round" "$ours" checksum -@ "$threads" "$input"
        measure samtools "$round" "$oracle" checksum -@ "$threads" "$input"
    else
        measure samtools "$round" "$oracle" checksum -@ "$threads" "$input"
        measure rsomics "$round" "$ours" checksum -@ "$threads" "$input"
    fi
    round=$((round + 1))
done

printf 'tool\tmean_wall_seconds\tmedian_wall_seconds\tmean_user_seconds\tmean_system_seconds\tmean_max_rss_bytes\tmedian_max_rss_bytes\n' > "$summary"
for tool in rsomics samtools; do
    means=$(awk -F '\t' -v tool="$tool" '$1 == tool { wall += $3; user += $4; sys += $5; rss += $6; n += 1 } END { printf "%.6f\t%.6f\t%.6f\t%.0f", wall / n, user / n, sys / n, rss / n }' "$ledger")
    median_wall=$(awk -F '\t' -v tool="$tool" '$1 == tool { print $3 }' "$ledger" | sort -n | awk '{ value[NR] = $1 } END { if (NR % 2) print value[(NR + 1) / 2]; else printf "%.6f", (value[NR / 2] + value[NR / 2 + 1]) / 2 }')
    median_rss=$(awk -F '\t' -v tool="$tool" '$1 == tool { print $6 }' "$ledger" | sort -n | awk '{ value[NR] = $1 } END { if (NR % 2) print value[(NR + 1) / 2]; else printf "%.0f", (value[NR / 2] + value[NR / 2 + 1]) / 2 }')
    mean_wall=$(printf '%s\n' "$means" | cut -f1)
    rest=$(printf '%s\n' "$means" | cut -f2-)
    printf '%s\t%s\t%s\t%s\t%s\n' "$tool" "$mean_wall" "$median_wall" "$rest" "$median_rss" >> "$summary"
done

awk -F '\t' '
    NR == 1 { next }
    $1 == "rsomics" { ours[$2] = $3 }
    $1 == "samtools" { oracle[$2] = $3 }
    END {
        for (round in ours) {
            difference = ours[round] - oracle[round]
            sum += difference
            squared += difference * difference
            if (difference < 0) wins += 1
            n += 1
        }
        mean = sum / n
        variance = n > 1 ? (squared - n * mean * mean) / (n - 1) : 0
        statistic = variance > 0 ? mean / sqrt(variance / n) : 0
        printf "paired_rounds\t%d\nrsomics_wins\t%d\npaired_mean_difference_seconds\t%.6f\npaired_sample_sd_seconds\t%.6f\npaired_t_statistic\t%.6f\n", n, wins, mean, sqrt(variance), statistic
    }
' "$ledger" > "$output/paired.tsv"
