#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

if [[ $# -lt 6 || $# -gt 8 ]]; then
    echo "usage: $0 RSOMICS_BAM SAMTOOLS MODE INPUT OUTPUT_DIR REPEATS [BED] [INDEX]" >&2
    exit 2
fi

ours=$1
samtools=$2
mode=$3
input=$4
output_dir=$5
repeats=$6
bed=${7:-}
index=${8:-}
batch=${BENCH_BATCH:-1}
ours_output=$output_dir/ours.tsv
samtools_output=$output_dir/samtools.tsv
timing=$output_dir/.timing
times=$output_dir/times.tsv

[[ $(uname -s) == Darwin ]]
[[ $repeats =~ ^[1-9][0-9]*$ ]]
[[ $batch =~ ^[1-9][0-9]*$ ]]
((repeats >= 2))
[[ $mode == coverage || $mode == idxstats || $mode == bedcov ]]
if [[ $mode == idxstats ]]; then
    [[ -n $index && -z $bed ]]
elif [[ $mode == bedcov ]]; then
    [[ -n $bed && -n $index ]]
else
    [[ -z $bed && -z $index ]]
fi
mkdir -p "$output_dir"

ours_command=("$ours" "$mode")
samtools_command=("$samtools" "$mode")
if [[ $mode == coverage ]]; then
    ours_command+=("$input")
    samtools_command+=("$input")
elif [[ $mode == idxstats ]]; then
    ours_command+=(-X "$input" "$index")
    samtools_command+=(-X "$input" "$index")
else
    ours_command+=(-X "$bed" "$input" "$index")
    samtools_command+=(-X "$bed" "$input" "$index")
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
    printf 'mode=%s\nrepeats=%s\nbatch=%s\n' "$mode" "$repeats" "$batch"
    if [[ -n $bed ]]; then
        shasum -a 256 "$bed"
        printf 'bed_rows=%s\n' "$(wc -l < "$bed" | tr -d ' ')"
    fi
    if [[ -n $index ]]; then
        shasum -a 256 "$index"
    fi
} > "$output_dir/environment.txt"
printf 'round\torder\ttool\treal_s\tuser_s\tsystem_s\tmax_rss_bytes\n' > "$times"

record_timing() {
    local round=$1
    local order=$2
    local tool=$3
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$round" "$order" "$tool" \
        "$(awk -v batch="$batch" '$1 == "real" { print $2 / batch }' "$timing")" \
        "$(awk -v batch="$batch" '$1 == "user" { print $2 / batch }' "$timing")" \
        "$(awk -v batch="$batch" '$1 == "sys" { print $2 / batch }' "$timing")" \
        "$(awk '/maximum resident set size/ { print $1 }' "$timing")" >> "$times"
}

run_ours() {
    local round=$1
    local order=$2
    /usr/bin/time -p -l -o "$timing" /bin/zsh -c \
        'count=$1; shift; repeat $count { "$@" > /dev/null }' \
        zsh "$batch" "${ours_command[@]}"
    record_timing "$round" "$order" rsomics
}

run_samtools() {
    local round=$1
    local order=$2
    /usr/bin/time -p -l -o "$timing" /bin/zsh -c \
        'count=$1; shift; repeat $count { "$@" > /dev/null }' \
        zsh "$batch" "${samtools_command[@]}"
    record_timing "$round" "$order" samtools
}

"${ours_command[@]}" > "$ours_output"
"${samtools_command[@]}" > "$samtools_output"
cmp "$ours_output" "$samtools_output"
shasum -a 256 "$ours_output" > "$output_dir/output.sha256"

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
    "$output_dir/output.sha256" > "$output_dir/artifacts.sha256"
rm -f "$timing"
