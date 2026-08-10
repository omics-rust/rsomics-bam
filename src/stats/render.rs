use std::fmt::Write as _;
use std::io::Write;

use rsomics_common::{Result, RsomicsError};

use super::record::{Accumulator, BaseCounts, GcDepth, QualityCycles};
use super::{CoverageBins, Options};

pub(crate) fn write(
    mut output: impl Write,
    stats: &Accumulator,
    options: Options<'_>,
    scope: Option<([u8; 2], &[u8])>,
) -> Result<()> {
    let mut text = String::new();
    writeln!(
        text,
        "# This file was produced by rsomics-bam stats ({}) and can be plotted using plot-bamstats",
        env!("CARGO_PKG_VERSION")
    )
    .unwrap();
    if let Some((tag, value)) = scope {
        writeln!(
            text,
            "# This file contains statistics only for reads with tag: {}={}",
            String::from_utf8_lossy(&tag),
            String::from_utf8_lossy(value)
        )
        .unwrap();
    } else {
        writeln!(text, "# This file contains statistics for all reads.").unwrap();
    }
    writeln!(text, "# The command line was:  rsomics-bam stats").unwrap();
    write_body(&mut text, stats, options);
    output
        .write_all(text.as_bytes())
        .map_err(RsomicsError::Io)?;
    output.flush().map_err(RsomicsError::Io)
}

fn write_body(text: &mut String, stats: &Accumulator, options: Options<'_>) {
    let summary = &stats.summary;
    let sequences = summary.first + summary.last + summary.other;
    let (insert_sizes, insert_end, insert_average, insert_deviation, orientations) =
        insert_statistics(stats, options.insert_bulk);
    writeln!(
        text,
        "# CHK, Checksum\t[2]Read Names\t[3]Sequences\t[4]Qualities"
    )
    .unwrap();
    writeln!(
        text,
        "# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)"
    )
    .unwrap();
    writeln!(
        text,
        "CHK\t{:08x}\t{:08x}\t{:08x}",
        stats.checksum.names, stats.checksum.sequences, stats.checksum.qualities
    )
    .unwrap();
    writeln!(
        text,
        "# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part."
    )
    .unwrap();
    writeln!(
        text,
        "SN\traw total sequences:\t{}\t# excluding supplementary and secondary reads",
        summary.filtered + sequences
    )
    .unwrap();
    writeln!(text, "SN\tfiltered sequences:\t{}", summary.filtered).unwrap();
    writeln!(text, "SN\tsequences:\t{sequences}").unwrap();
    writeln!(
        text,
        "SN\tis sorted:\t{}\t# {} by coordinate",
        usize::from(stats.sorted),
        if stats.sorted { "sorted" } else { "not sorted" }
    )
    .unwrap();
    writeln!(text, "SN\t1st fragments:\t{}", summary.first).unwrap();
    writeln!(text, "SN\tlast fragments:\t{}", summary.last).unwrap();
    writeln!(
        text,
        "SN\treads mapped:\t{}",
        summary.mapped_paired + summary.mapped_single
    )
    .unwrap();
    writeln!(
        text,
        "SN\treads mapped and paired:\t{}\t# paired-end technology bit set + both mates mapped",
        summary.mapped_paired
    )
    .unwrap();
    writeln!(text, "SN\treads unmapped:\t{}", summary.unmapped).unwrap();
    writeln!(
        text,
        "SN\treads properly paired:\t{}\t# proper-pair bit set",
        summary.properly_paired
    )
    .unwrap();
    writeln!(
        text,
        "SN\treads paired:\t{}\t# paired-end technology bit set",
        summary.paired
    )
    .unwrap();
    writeln!(
        text,
        "SN\treads duplicated:\t{}\t# PCR or optical duplicate bit set",
        summary.duplicated
    )
    .unwrap();
    writeln!(text, "SN\treads MQ0:\t{}\t# mapped and MQ=0", summary.mq0).unwrap();
    writeln!(text, "SN\treads QC failed:\t{}", summary.qc_failed).unwrap();
    writeln!(text, "SN\tnon-primary alignments:\t{}", summary.secondary).unwrap();
    writeln!(
        text,
        "SN\tsupplementary alignments:\t{}",
        summary.supplementary
    )
    .unwrap();
    writeln!(
        text,
        "SN\ttotal length:\t{}\t# ignores clipping",
        summary.total_length
    )
    .unwrap();
    writeln!(
        text,
        "SN\ttotal first fragment length:\t{}\t# ignores clipping",
        summary.first_length
    )
    .unwrap();
    writeln!(
        text,
        "SN\ttotal last fragment length:\t{}\t# ignores clipping",
        summary.last_length
    )
    .unwrap();
    writeln!(
        text,
        "SN\tbases mapped:\t{}\t# ignores clipping",
        summary.mapped_bases
    )
    .unwrap();
    writeln!(
        text,
        "SN\tbases mapped (cigar):\t{}\t# more accurate",
        summary.cigar_bases
    )
    .unwrap();
    writeln!(text, "SN\tbases trimmed:\t{}", summary.trimmed_bases).unwrap();
    writeln!(text, "SN\tbases duplicated:\t{}", summary.duplicated_bases).unwrap();
    writeln!(
        text,
        "SN\tmismatches:\t{}\t# from NM fields",
        summary.mismatches
    )
    .unwrap();
    let error_rate = if summary.cigar_bases == 0 {
        0.0f32
    } else {
        summary.mismatches as f32 / summary.cigar_bases as f32
    };
    writeln!(
        text,
        "SN\terror rate:\t{}\t# mismatches / bases mapped (cigar)",
        scientific(error_rate)
    )
    .unwrap();
    let average_length = ratio(summary.total_length, sequences);
    writeln!(text, "SN\taverage length:\t{average_length:.0}").unwrap();
    writeln!(
        text,
        "SN\taverage first fragment length:\t{:.0}",
        ratio(summary.first_length, summary.first)
    )
    .unwrap();
    writeln!(
        text,
        "SN\taverage last fragment length:\t{:.0}",
        ratio(summary.last_length, summary.last)
    )
    .unwrap();
    writeln!(text, "SN\tmaximum length:\t{}", stats.max_length).unwrap();
    writeln!(
        text,
        "SN\tmaximum first fragment length:\t{}",
        stats.max_first_length
    )
    .unwrap();
    writeln!(
        text,
        "SN\tmaximum last fragment length:\t{}",
        stats.max_last_length
    )
    .unwrap();
    writeln!(
        text,
        "SN\taverage quality:\t{:.1}",
        if summary.total_length == 0 {
            0.0
        } else {
            summary.quality_sum / summary.total_length as f64
        }
    )
    .unwrap();
    writeln!(text, "SN\tinsert size average:\t{insert_average:.1}").unwrap();
    writeln!(
        text,
        "SN\tinsert size standard deviation:\t{insert_deviation:.1}"
    )
    .unwrap();
    writeln!(text, "SN\tinward oriented pairs:\t{}", orientations[0]).unwrap();
    writeln!(text, "SN\toutward oriented pairs:\t{}", orientations[1]).unwrap();
    writeln!(
        text,
        "SN\tpairs with other orientation:\t{}",
        orientations[2]
    )
    .unwrap();
    writeln!(
        text,
        "SN\tpairs on different chromosomes:\t{}",
        summary.anomalous / 2
    )
    .unwrap();
    writeln!(
        text,
        "SN\tpercentage of properly paired reads (%):\t{:.1}",
        if sequences == 0 {
            0.0
        } else {
            100.0 * summary.properly_paired as f64 / sequences as f64
        }
    )
    .unwrap();
    if let Some(target_bases) = stats.target_bases {
        writeln!(text, "SN\tbases inside the target:\t{target_bases}").unwrap();
        writeln!(
            text,
            "SN\tpercentage of target genome with coverage > {} (%):\t{:.2}",
            options.coverage_threshold,
            100.0 * stats.coverage.bases_above(options.coverage_threshold) as f64
                / target_bases as f64
        )
        .unwrap();
    }

    let maximum_quality = if stats.max_quality < 255 {
        stats.max_quality + 1
    } else {
        stats.max_quality
    };
    write_quality(
        text,
        "First Fragment",
        "FFQ",
        &stats.first_qualities,
        stats.max_first_length,
        maximum_quality,
    );
    write_quality(
        text,
        "Last Fragment",
        "LFQ",
        &stats.last_qualities,
        stats.max_last_length,
        maximum_quality,
    );
    if let Some(cycles) = &stats.mismatch_cycles {
        writeln!(
            text,
            "# Mismatches per cycle and quality. Use `grep ^MPC | cut -f 2-` to extract this part."
        )
        .unwrap();
        writeln!(text, "# Columns correspond to qualities, rows to cycles. First column is the cycle number, second").unwrap();
        writeln!(
            text,
            "# is the number of N's and the rest is the number of mismatches"
        )
        .unwrap();
        for cycle in 0..=stats.max_length {
            write!(text, "MPC\t{}", cycle + 1).unwrap();
            for quality in 0..=maximum_quality {
                let value = cycles.get(cycle, quality);
                write!(text, "\t{value}").unwrap();
            }
            text.push('\n');
        }
    }
    write_gc(text, "first", "GCF", &stats.first_gc);
    write_gc(text, "last", "GCL", &stats.last_gc);
    writeln!(text, "# ACGT content per cycle. Use `grep ^GCC | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N and O counts as a percentage of all A/C/G/T bases [%]").unwrap();
    for cycle in 0..stats.max_length {
        let counts = value(&stats.first_bases, cycle).add(value(&stats.last_bases, cycle));
        write_base_percentages(text, "GCC", cycle, counts, true);
    }
    writeln!(text, "# ACGT content per cycle, read oriented. Use `grep ^GCT | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]").unwrap();
    for cycle in 0..stats.max_length {
        write_base_percentages(
            text,
            "GCT",
            cycle,
            value(&stats.oriented_bases, cycle),
            false,
        );
    }
    write_fragment_bases(
        text,
        "first",
        "FBC",
        "FTC",
        &stats.first_bases,
        stats.max_length,
    );
    write_fragment_bases(
        text,
        "last",
        "LBC",
        "LTC",
        &stats.last_bases,
        stats.max_length,
    );
    write_barcodes(text, stats);

    writeln!(text, "# Insert sizes. Use `grep ^IS | cut -f 2-` to extract this part. The columns are: insert size, pairs total, inward oriented pairs, outward oriented pairs, other pairs").unwrap();
    for size in 0..insert_end {
        let values = insert_sizes.get(&size).copied().unwrap_or_default();
        let total = values.iter().sum::<u64>();
        if !options.sparse || total != 0 {
            writeln!(
                text,
                "IS\t{size}\t{total}\t{}\t{}\t{}",
                values[0], values[1], values[2]
            )
            .unwrap();
        }
    }
    write_lengths(text, "Read lengths", "RL", &stats.read_lengths);
    write_lengths(
        text,
        "Read lengths - first fragments",
        "FRL",
        &stats.first_lengths,
    );
    write_lengths(
        text,
        "Read lengths - last fragments",
        "LRL",
        &stats.last_lengths,
    );
    writeln!(text, "# Mapping qualities for reads !(UNMAP|SECOND|SUPPL|QCFAIL|DUP). Use `grep ^MAPQ | cut -f 2-` to extract this part. The columns are: mapq, count").unwrap();
    for (quality, &count) in stats.mapping_qualities.iter().enumerate() {
        if count != 0 {
            writeln!(text, "MAPQ\t{quality}\t{count}").unwrap();
        }
    }
    writeln!(text, "# Indel distribution. Use `grep ^ID | cut -f 2-` to extract this part. The columns are: length, number of insertions, number of deletions").unwrap();
    let maximum_indel = stats
        .insertions
        .keys()
        .chain(stats.deletions.keys())
        .copied()
        .max()
        .unwrap_or(0);
    for length in 1..=maximum_indel {
        let insertions = stats.insertions.get(&length).copied().unwrap_or(0);
        let deletions = stats.deletions.get(&length).copied().unwrap_or(0);
        if insertions != 0 || deletions != 0 {
            writeln!(text, "ID\t{length}\t{insertions}\t{deletions}").unwrap();
        }
    }
    writeln!(text, "# Indels per cycle. Use `grep ^IC | cut -f 2-` to extract this part. The columns are: cycle, number of insertions (fwd), .. (rev) , number of deletions (fwd), .. (rev)").unwrap();
    for cycle in 0..stats
        .insertion_cycles
        .len()
        .max(stats.deletion_cycles.len())
    {
        let insertions = stats.insertion_cycles.get(cycle).copied().unwrap_or([0; 2]);
        let deletions = stats.deletion_cycles.get(cycle).copied().unwrap_or([0; 2]);
        if insertions != [0; 2] || deletions != [0; 2] {
            writeln!(
                text,
                "IC\t{}\t{}\t{}\t{}\t{}",
                cycle + 1,
                insertions[0],
                insertions[1],
                deletions[0],
                deletions[1]
            )
            .unwrap();
        }
    }
    if stats.sorted {
        write_coverage(text, stats, options.coverage);
        writeln!(text, "# GC-depth. Use `grep ^GCD | cut -f 2-` to extract this part. The columns are: GC%, unique sequence percentiles, 10th, 25th, 50th, 75th and 90th depth percentile").unwrap();
        write_gc_depth(text, stats, options.gc_depth, average_length);
    }
    if let Some(reference) = &stats.reference_stats {
        writeln!(
            text,
            "# Reference statistics. Use `grep ^RFS | cut -f 2-` to extract this part."
        )
        .unwrap();
        writeln!(text, "# Total count, Output count, Average GC, Min length, Max length, Average length, Total length in first row.").unwrap();
        writeln!(
            text,
            "# Sequence name, Length, GC content, Unknown count in following rows."
        )
        .unwrap();
        writeln!(
            text,
            "RFS\t{}\t{}\t{:.2}\t{}\t{}\t{:.2}\t{}",
            reference.total_count,
            reference.output_count,
            reference.average_gc,
            reference.minimum_length,
            reference.maximum_length,
            reference.average_length,
            reference.total_length
        )
        .unwrap();
        for sequence in &reference.sequences {
            writeln!(
                text,
                "RFS\t{}\t{}\t{:.2}\t{}",
                sequence.name, sequence.length, sequence.gc, sequence.unknown
            )
            .unwrap();
        }
    }
}

fn write_gc_depth(text: &mut String, stats: &Accumulator, raw_bin_size: f64, average: f64) {
    let count = stats.gc_depth.len().saturating_sub(1);
    if count == 0 {
        return;
    }
    let reference = stats.mismatch_cycles.is_some();
    let mut bins = stats.gc_depth.clone();
    for bin in bins.iter_mut().take(count) {
        bin.gc = if reference {
            (100.0 * bin.gc).round_ties_even()
        } else if bin.depth == 0 {
            bin.gc
        } else {
            (100.0 * bin.gc / bin.depth as f32).round_ties_even()
        };
    }
    bins.sort_by(|left, right| {
        left.gc
            .total_cmp(&right.gc)
            .then(left.depth.cmp(&right.depth))
    });
    let mut index = 0usize;
    while index < count {
        let gc = bins[index].gc;
        let mut end = index;
        while end < count && (bins[end].gc - gc).abs() < 0.1 {
            end += 1;
        }
        let group = &bins[index..end];
        let scale = average / raw_bin_size;
        writeln!(
            text,
            "GCD\t{gc:.1}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
            (end + 1) as f64 * 100.0 / (count + 1) as f64,
            gc_percentile(group, 10) * scale,
            gc_percentile(group, 25) * scale,
            gc_percentile(group, 50) * scale,
            gc_percentile(group, 75) * scale,
            gc_percentile(group, 90) * scale,
        )
        .unwrap();
        index = end;
    }
}

fn gc_percentile(values: &[GcDepth], percentile: usize) -> f64 {
    let n = percentile as f64 * (values.len() + 1) as f64 / 100.0;
    let k = n as usize;
    if k == 0 {
        return values[0].depth as f64;
    }
    if k >= values.len() {
        return values.last().unwrap().depth as f64;
    }
    let fraction = n - k as f64;
    values[k - 1].depth as f64 + fraction * (values[k].depth as f64 - values[k - 1].depth as f64)
}

fn write_barcodes(text: &mut String, stats: &Accumulator) {
    for barcode in &stats.barcodes {
        if barcode.bases.is_empty() {
            continue;
        }
        let sequence_tag = String::from_utf8_lossy(&barcode.sequence_tag);
        let quality_tag = String::from_utf8_lossy(&barcode.quality_tag);
        writeln!(text, "# ACGT content per cycle for barcodes. Use `grep ^{sequence_tag}C | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N counts as a percentage of all A/C/G/T bases [%]").unwrap();
        for (index, &counts) in barcode.bases.iter().enumerate() {
            if barcode.separator == Some(index) || counts.acgt() == 0 {
                continue;
            }
            let (part, cycle) = barcode_cycle(index, barcode.separator);
            writeln!(
                text,
                "{sequence_tag}C{part}\t{cycle}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}",
                percent(counts.a, counts.acgt()),
                percent(counts.c, counts.acgt()),
                percent(counts.g, counts.acgt()),
                percent(counts.t, counts.acgt()),
                percent(counts.n, counts.acgt())
            )
            .unwrap();
        }
        writeln!(
            text,
            "# Barcode Qualities. Use `grep ^{quality_tag}Q | cut -f 2-` to extract this part."
        )
        .unwrap();
        writeln!(text, "# Columns correspond to qualities and rows to barcode cycles. First column is the cycle number.").unwrap();
        for index in 0..barcode.qualities.len() {
            if barcode.separator == Some(index) {
                continue;
            }
            let (part, cycle) = barcode_cycle(index, barcode.separator);
            write!(text, "{quality_tag}Q{part}\t{cycle}").unwrap();
            if let Some(maximum) = barcode.maximum_quality {
                for quality in 0..=maximum {
                    write!(text, "\t{}", barcode.qualities.get(index, quality)).unwrap();
                }
            }
            text.push('\n');
        }
    }
}

fn barcode_cycle(index: usize, separator: Option<usize>) -> (usize, usize) {
    match separator {
        Some(separator) if index > separator => (2, index - separator),
        _ => (1, index + 1),
    }
}

fn write_quality(
    text: &mut String,
    name: &str,
    prefix: &str,
    values: &QualityCycles,
    length: usize,
    maximum_quality: usize,
) {
    writeln!(
        text,
        "# {name} Qualities. Use `grep ^{prefix} | cut -f 2-` to extract this part."
    )
    .unwrap();
    writeln!(
        text,
        "# Columns correspond to qualities and rows to cycles. First column is the cycle number."
    )
    .unwrap();
    for cycle in 0..length {
        write!(text, "{prefix}\t{}", cycle + 1).unwrap();
        for quality in 0..=maximum_quality {
            let count = values.get(cycle, quality);
            write!(text, "\t{count}").unwrap();
        }
        text.push('\n');
    }
}

fn write_gc(text: &mut String, name: &str, prefix: &str, values: &[u64]) {
    writeln!(
        text,
        "# GC Content of {name} fragments. Use `grep ^{prefix} | cut -f 2-` to extract this part."
    )
    .unwrap();
    let mut previous = 0;
    for index in 0..values.len() {
        if values[index] == values[previous] {
            continue;
        }
        writeln!(
            text,
            "{prefix}\t{:.2}\t{}",
            (index + previous) as f64 * 50.0 / 199.0,
            values[previous]
        )
        .unwrap();
        previous = index;
    }
}

fn write_base_percentages(
    text: &mut String,
    prefix: &str,
    cycle: usize,
    counts: BaseCounts,
    include_unknown: bool,
) {
    let total = counts.acgt();
    if total == 0 {
        return;
    }
    write!(
        text,
        "{prefix}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}",
        cycle + 1,
        percent(counts.a, total),
        percent(counts.c, total),
        percent(counts.g, total),
        percent(counts.t, total)
    )
    .unwrap();
    if include_unknown {
        write!(
            text,
            "\t{:.2}\t{:.2}",
            percent(counts.n, total),
            percent(counts.other, total)
        )
        .unwrap();
    }
    text.push('\n');
}

fn write_fragment_bases(
    text: &mut String,
    name: &str,
    prefix: &str,
    total_prefix: &str,
    values: &[BaseCounts],
    maximum_length: usize,
) {
    writeln!(text, "# ACGT content per cycle for {name} fragments. Use `grep ^{prefix} | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N and O counts as a percentage of all A/C/G/T bases [%]").unwrap();
    let mut totals = BaseCounts::default();
    for cycle in 0..maximum_length {
        let counts = value(values, cycle);
        totals = totals.add(counts);
        write_base_percentages(text, prefix, cycle, counts, true);
    }
    writeln!(text, "# ACGT raw counters for {name} fragments. Use `grep ^{total_prefix} | cut -f 2-` to extract this part. The columns are: A,C,G,T,N base counters").unwrap();
    writeln!(
        text,
        "{total_prefix}\t{}\t{}\t{}\t{}\t{}",
        totals.a, totals.c, totals.g, totals.t, totals.n
    )
    .unwrap();
}

fn write_lengths(
    text: &mut String,
    title: &str,
    prefix: &str,
    values: &std::collections::BTreeMap<usize, u64>,
) {
    writeln!(text, "# {title}. Use `grep ^{prefix} | cut -f 2-` to extract this part. The columns are: read length, count").unwrap();
    for (length, count) in values {
        writeln!(text, "{prefix}\t{length}\t{count}").unwrap();
    }
}

fn write_coverage(text: &mut String, stats: &Accumulator, bins: CoverageBins) {
    writeln!(
        text,
        "# Coverage distribution. Use `grep ^COV | cut -f 2-` to extract this part."
    )
    .unwrap();
    let histogram = stats.coverage_histogram(bins);
    if histogram[0] != 0 {
        writeln!(
            text,
            "COV\t[<{}]\t{}\t{}",
            bins.minimum,
            bins.minimum - 1,
            histogram[0]
        )
        .unwrap();
    }
    for (index, count) in histogram
        .iter()
        .enumerate()
        .take(histogram.len() - 1)
        .skip(1)
    {
        if *count != 0 {
            let start = bins.minimum + (index - 1) * bins.step;
            let end = bins.minimum + index * bins.step - 1;
            writeln!(text, "COV\t[{start}-{end}]\t{end}\t{count}").unwrap();
        }
    }
    let last = histogram.len() - 1;
    if histogram[last] != 0 {
        let boundary = bins.minimum + (last - 1) * bins.step - 1;
        writeln!(text, "COV\t[{boundary}<]\t{boundary}\t{}", histogram[last]).unwrap();
    }
}

fn insert_statistics(
    stats: &Accumulator,
    bulk_fraction: f64,
) -> (
    std::collections::BTreeMap<usize, [u64; 3]>,
    usize,
    f64,
    f64,
    [u64; 3],
) {
    let values = stats
        .insert_sizes
        .iter()
        .map(|(&size, value)| (size, [value[0] / 2, value[1] / 2, value[2] / 2]))
        .collect::<std::collections::BTreeMap<_, _>>();
    let orientations = values.values().fold([0u64; 3], |mut total, value| {
        for index in 0..3 {
            total[index] += value[index];
        }
        total
    });
    let total = orientations.iter().sum::<u64>();
    let mut cumulative = 0u64;
    let mut weighted = 0f64;
    let mut end = 0;
    let mut selected = total;
    for (&size, value) in &values {
        let count = value.iter().sum::<u64>();
        if count != 0 {
            end = size + 1;
        }
        cumulative += count;
        weighted += size as f64 * count as f64;
        if total != 0 && cumulative as f64 / total as f64 > bulk_fraction {
            end = size + 1;
            selected = cumulative;
            break;
        }
    }
    let average = weighted / selected.max(1) as f64;
    let variance = values
        .iter()
        .filter(|(size, _)| **size > 0 && **size < end)
        .map(|(&size, value)| {
            let count = value.iter().sum::<u64>() as f64;
            count * (size as f64 - average).powi(2) / selected.max(1) as f64
        })
        .sum::<f64>();
    let deviation = if variance == 0.0 {
        0.0
    } else {
        variance.sqrt()
    };
    (values, end, average, deviation, orientations)
}

fn ratio(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

fn percent(value: u64, total: u64) -> f64 {
    100.0 * value as f64 / total as f64
}

fn value(values: &[BaseCounts], index: usize) -> BaseCounts {
    values.get(index).copied().unwrap_or_default()
}

fn scientific(value: f32) -> String {
    let raw = format!("{value:.6e}");
    let (mantissa, exponent) = raw.split_once('e').unwrap();
    let exponent = exponent.parse::<i32>().unwrap();
    format!("{mantissa}e{exponent:+03}")
}
