#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

if [[ $# -lt 5 || $# -gt 7 ]]; then
    echo "usage: $0 RSOMICS_BAM SAMTOOLS INPUT OUTPUT_DIR REPEATS [BATCH] [MODE]" >&2
    exit 2
fi

ours=$1
samtools=$2
input=$3
output_dir=$4
repeats=$5
batch=${6:-20}
mode=${7:-default}
timing=$output_dir/.timing
times=$output_dir/times.tsv

[[ $(uname -s) == Darwin ]]
[[ $repeats =~ ^[1-9][0-9]*$ ]]
[[ $batch =~ ^[1-9][0-9]*$ ]]
((repeats >= 2))
[[ $mode == default || $mode == encodings ]]
mkdir -p "$output_dir"

ours_command=("$ours" cram-size)
samtools_command=("$samtools" cram-size)
if [[ $mode == encodings ]]; then
    ours_command+=(-e)
    samtools_command+=(-e)
fi
ours_command+=("$input")
samtools_command+=("$input")

"${ours_command[@]}" > "$output_dir/rsomics.txt"
"${samtools_command[@]}" > "$output_dir/samtools.txt"
cmp "$output_dir/rsomics.txt" "$output_dir/samtools.txt"

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
    sed -n '/^Number of containers/,$p' "$output_dir/samtools.txt"
    printf 'repeats=%s\nbatch=%s\nmode=%s\n' "$repeats" "$batch" "$mode"
} > "$output_dir/environment.txt"
shasum -a 256 "$output_dir/rsomics.txt" > "$output_dir/output.sha256"
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

measure() {
    local round=$1
    local order=$2
    local tool=$3
    shift 3
    /usr/bin/time -p -l -o "$timing" /bin/zsh -c \
        'count=$1; shift; repeat $count { "$@" > /dev/null }' \
        zsh "$batch" "$@"
    record_timing "$round" "$order" "$tool"
}

measure 0 warmup rsomics "${ours_command[@]}"
measure 0 warmup samtools "${samtools_command[@]}"
printf 'round\torder\ttool\treal_s\tuser_s\tsystem_s\tmax_rss_bytes\n' > "$times"
for ((round = 1; round <= repeats; round++)); do
    if ((round % 2 == 1)); then
        order=rsomics-samtools
        measure "$round" "$order" rsomics "${ours_command[@]}"
        measure "$round" "$order" samtools "${samtools_command[@]}"
    else
        order=samtools-rsomics
        measure "$round" "$order" samtools "${samtools_command[@]}"
        measure "$round" "$order" rsomics "${ours_command[@]}"
    fi
done

awk -F '\t' -v repeats="$repeats" '
    function median(tool, n,    i, j, key, value) {
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
                user[tool] / count[tool], system_time[tool] / count[tool], rss[tool] / count[tool]
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
    "$output_dir/output.sha256" > "$output_dir/artifacts.sha256"
rm -f "$timing"
