#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

if [[ $# -ne 6 ]]; then
    echo "usage: $0 RSOMICS_BAM SAMTOOLS INPUT OUTPUT_DIR REPEATS THREADS" >&2
    exit 2
fi

ours=$1
samtools=$2
input=$3
output_dir=$4
repeats=$5
threads=$6
ours_output=$output_dir/ours.bam
samtools_output=$output_dir/samtools.bam
timing=$output_dir/.timing
times=$output_dir/times.tsv

[[ $(uname -s) == Darwin ]]
[[ $repeats =~ ^[1-9][0-9]*$ ]]
((repeats >= 5))
[[ $threads =~ ^[0-9]+$ ]]
[[ $($samtools --version | sed -n '1p') == "samtools 1.24" ]]
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
    printf 'records=%s\n' "$("$samtools" view -c "$input")"
    printf 'threads=%s\nrepeats=%s\n' "$threads" "$repeats"
} > "$output_dir/environment.txt"
printf 'round\torder\ttool\treal_s\tuser_s\tsystem_s\tmax_rss_bytes\n' > "$times"

record_timing() {
    local round=$1
    local order=$2
    local tool=$3
    awk -v round="$round" -v order="$order" -v tool="$tool" '
        $1 == "real" { real = $2 }
        $1 == "user" { user = $2 }
        $1 == "sys" { system_time = $2 }
        /maximum resident set size/ { rss = $1 }
        END {
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", round, order, tool,
                real, user, system_time, rss
        }
    ' "$timing" >> "$times"
}

run_ours() {
    local round=$1
    local order=$2
    /usr/bin/time -lp -o "$timing" \
        "$ours" reset --no-PG -O bam -@ "$threads" -o "$ours_output" "$input"
    record_timing "$round" "$order" rsomics
}

run_samtools() {
    local round=$1
    local order=$2
    /usr/bin/time -lp -o "$timing" \
        "$samtools" reset --no-PG -O bam -@ "$threads" -o "$samtools_output" "$input"
    record_timing "$round" "$order" samtools
}

validate_pair() {
    "$samtools" view -h --no-PG "$ours_output" | shasum -a 256 > "$output_dir/ours.sha256"
    "$samtools" view -h --no-PG "$samtools_output" | shasum -a 256 \
        > "$output_dir/samtools.sha256"
    cmp "$output_dir/ours.sha256" "$output_dir/samtools.sha256"
}

"$ours" --json reset --no-PG -O bam -@ "$threads" -o "$ours_output" "$input" \
    > "$output_dir/rsomics-summary.json"
"$samtools" reset --no-PG -O bam -@ "$threads" -o "$samtools_output" "$input"
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
    "$output_dir/ours.sha256" > "$output_dir/artifacts.sha256"
rm -f "$timing"
