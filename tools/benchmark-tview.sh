#!/bin/sh
set -eu

if [ "$#" -ne 9 ]; then
    printf 'usage: %s RSOMICS_BAM SAMTOOLS INPUT REFERENCE OUTPUT_DIR ROUNDS THREADS POSITION WIDTH\n' "$0" >&2
    exit 2
fi

ours=$1
oracle=$2
input=$3
reference=$4
output=$5
rounds=$6
threads=$7
position=$8
width=$9

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

"$ours" tview -@ "$threads" -d T -w "$width" -p "$position" -T "$reference" "$input" > "$output/rsomics.txt"
"$oracle" tview -d T -w "$width" -p "$position" "$input" "$reference" > "$output/samtools.txt"
cmp "$output/rsomics.txt" "$output/samtools.txt"

"$ours" tview -@ "$threads" -d H -w "$width" -p "$position" -T "$reference" "$input" > "$output/rsomics.html"
"$oracle" tview -d H -w "$width" -p "$position" "$input" "$reference" > "$output/samtools.html"

ledger="$output/timings.tsv"
summary="$output/summary.tsv"
paired="$output/paired.tsv"
interactive="$output/interactive.tsv"
environment="$output/environment.txt"

{
    uname -a
    sw_vers 2>/dev/null || true
    "$ours" --version
    "$oracle" --version
    printf 'additional_bam_workers %s\n' "$threads"
    printf 'position %s\n' "$position"
    printf 'width %s\n' "$width"
    printf 'rounds %s\n' "$rounds"
    printf 'input_bytes %s\n' "$(stat -f %z "$input")"
    printf 'input_records %s\n' "$("$oracle" view -c "$input")"
    shasum -a 256 "$ours" "$oracle" "$input" "$reference" "$0" \
        "$(dirname "$0")/tview-pty.py" "$output/rsomics.txt" "$output/samtools.txt" \
        "$output/rsomics.html" "$output/samtools.html"
} > "$environment"

printf 'mode\ttool\tround\treal_seconds\tuser_seconds\tsystem_seconds\tmax_rss_bytes\n' > "$ledger"

measure() {
    mode=$1
    tool=$2
    round=$3
    shift 3
    timing=$(mktemp "$output/time.XXXXXX")
    /usr/bin/time -lp "$@" > /dev/null 2> "$timing"
    real=$(awk '$1 == "real" { print $2 }' "$timing")
    user=$(awk '$1 == "user" { print $2 }' "$timing")
    system=$(awk '$1 == "sys" { print $2 }' "$timing")
    rss=$(awk '$2 == "maximum" && $3 == "resident" { print $1 }' "$timing")
    rm "$timing"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$mode" "$tool" "$round" "$real" "$user" "$system" "$rss" >> "$ledger"
}

round=1
while [ "$round" -le "$rounds" ]; do
    if [ $((round % 2)) -eq 1 ]; then
        measure text rsomics "$round" "$ours" tview -@ "$threads" -d T -w "$width" -p "$position" -T "$reference" "$input"
        measure text samtools "$round" "$oracle" tview -d T -w "$width" -p "$position" "$input" "$reference"
        measure html rsomics "$round" "$ours" tview -@ "$threads" -d H -w "$width" -p "$position" -T "$reference" "$input"
        measure html samtools "$round" "$oracle" tview -d H -w "$width" -p "$position" "$input" "$reference"
    else
        measure text samtools "$round" "$oracle" tview -d T -w "$width" -p "$position" "$input" "$reference"
        measure text rsomics "$round" "$ours" tview -@ "$threads" -d T -w "$width" -p "$position" -T "$reference" "$input"
        measure html samtools "$round" "$oracle" tview -d H -w "$width" -p "$position" "$input" "$reference"
        measure html rsomics "$round" "$ours" tview -@ "$threads" -d H -w "$width" -p "$position" -T "$reference" "$input"
    fi
    round=$((round + 1))
done

printf 'mode\ttool\tmean_wall_seconds\tmedian_wall_seconds\tmean_user_seconds\tmean_system_seconds\tmean_max_rss_bytes\tmedian_max_rss_bytes\n' > "$summary"
for mode in text html; do
    for tool in rsomics samtools; do
        awk -F '\t' -v mode="$mode" -v tool="$tool" '$1 == mode && $2 == tool { print $4 }' "$ledger" | sort -n > "$output/wall.values"
        awk -F '\t' -v mode="$mode" -v tool="$tool" '$1 == mode && $2 == tool { print $7 }' "$ledger" | sort -n > "$output/rss.values"
        median_wall=$(awk '{ value[NR]=$1 } END { if (NR%2) print value[(NR+1)/2]; else printf "%.6f", (value[NR/2]+value[NR/2+1])/2 }' "$output/wall.values")
        median_rss=$(awk '{ value[NR]=$1 } END { if (NR%2) print value[(NR+1)/2]; else printf "%.0f", (value[NR/2]+value[NR/2+1])/2 }' "$output/rss.values")
        means=$(awk -F '\t' -v mode="$mode" -v tool="$tool" '$1 == mode && $2 == tool { wall+=$4; user+=$5; sys+=$6; rss+=$7; n++ } END { printf "%.6f\t%.6f\t%.6f\t%.0f", wall/n, user/n, sys/n, rss/n }' "$ledger")
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$mode" "$tool" "$(printf '%s\n' "$means" | cut -f1)" "$median_wall" "$(printf '%s\n' "$means" | cut -f2-)" "$median_rss" >> "$summary"
    done
done
rm "$output/wall.values" "$output/rss.values"

printf 'mode\trounds\trsomics_wins\tties\tmean_delta_seconds\tsample_sd_delta_seconds\tmean_ratio\n' > "$paired"
for mode in text html; do
    awk -F '\t' -v mode="$mode" '
        $1 == mode && $2 == "rsomics" { ours[$3]=$4 }
        $1 == mode && $2 == "samtools" { oracle[$3]=$4 }
        END {
            for (round in ours) {
                delta=ours[round]-oracle[round]
                sum+=delta
                sumsq+=delta*delta
                ratio+=ours[round]/oracle[round]
                if (delta < 0) wins++
                else if (delta == 0) ties++
                n++
            }
            mean=sum/n
            sd=n > 1 ? sqrt((sumsq-sum*sum/n)/(n-1)) : 0
            printf "%s\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\n", mode, n, wins, ties, mean, sd, ratio/n
        }
    ' "$ledger" >> "$paired"
done

printf 'round\tseconds\n' > "$interactive"
round=1
while [ "$round" -le "$rounds" ]; do
    seconds=$("$(dirname "$0")/tview-pty.py" "$ours" tview -@ "$threads" -p "$position" -T "$reference" "$input")
    printf '%s\t%s\n' "$round" "$seconds" >> "$interactive"
    round=$((round + 1))
done
