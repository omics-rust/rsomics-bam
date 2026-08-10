use std::io::{self, Write};

use super::{BlockSummary, MethodSummary, Options, Report};

pub(super) fn write(report: &Report, options: Options, mut writer: impl Write) -> io::Result<()> {
    if options.encodings {
        for container in &report.encodings {
            writeln!(writer, "Container encodings")?;
            for entry in &container.entries {
                writeln!(writer, "\t{}\t{}", entry.data_series, entry.encoding)?;
            }
            writeln!(writer)?;
        }
    }

    if options.verbose {
        writeln!(
            writer,
            "#   Content_ID  Uncomp.size    Comp.size   Ratio Method      Data_series"
        )?;
        for block in &report.blocks {
            write_verbose_block(&mut writer, block)?;
        }
    } else {
        writeln!(
            writer,
            "#   Content_ID  Uncomp.size    Comp.size   Ratio Method  Data_series"
        )?;
        for block in &report.blocks {
            write_compact_block(&mut writer, block)?;
        }
    }

    writeln!(writer)?;
    writeln!(writer, "Number of containers  {:18}", report.containers)?;
    writeln!(writer, "Number of slices      {:18}", report.slices)?;
    writeln!(writer, "Number of sequences   {:18}", report.sequences)?;
    writeln!(writer, "Number of bases       {:18}", report.bases)?;
    writeln!(writer, "Total file size       {:18}", report.file_size)?;
    writeln!(
        writer,
        "Format overhead size  {:18}",
        report.format_overhead_size
    )?;
    writer.flush()
}

fn write_compact_block(writer: &mut impl Write, block: &BlockSummary) -> io::Result<()> {
    let uncompressed = block
        .methods
        .iter()
        .map(|method| method.uncompressed_size)
        .sum::<u64>();
    let compressed = block
        .methods
        .iter()
        .map(|method| method.compressed_size)
        .sum::<u64>();
    let method = block
        .methods
        .iter()
        .filter(|method| method.compressed_size != 0)
        .map(|method| method.short.as_str())
        .collect::<String>();
    let method = if method.is_empty() { "." } else { &method };
    write_prefix(writer, block, uncompressed, compressed)?;
    write_ratio(writer, compressed, uncompressed)?;
    write!(writer, " {method:<7}")?;
    write_series(writer, block, true)?;
    writeln!(writer)
}

fn write_verbose_block(writer: &mut impl Write, block: &BlockSummary) -> io::Result<()> {
    let visible_count = block
        .methods
        .iter()
        .enumerate()
        .take_while(|(index, method)| *index == 0 || method.compressed_size != 0)
        .count();
    for (index, method) in visible_methods(&block.methods).enumerate() {
        write_prefix(
            writer,
            block,
            method.uncompressed_size,
            method.compressed_size,
        )?;
        write_ratio(writer, method.compressed_size, method.uncompressed_size)?;
        write!(writer, " {:<11}", method.method)?;
        write_series(writer, block, index + 1 == visible_count)?;
        writeln!(writer)?;
    }
    Ok(())
}

fn visible_methods(methods: &[MethodSummary]) -> impl Iterator<Item = &MethodSummary> {
    methods
        .iter()
        .enumerate()
        .take_while(|(index, method)| *index == 0 || method.compressed_size != 0)
        .map(|(_, method)| method)
}

fn write_prefix(
    writer: &mut impl Write,
    block: &BlockSummary,
    uncompressed: u64,
    compressed: u64,
) -> io::Result<()> {
    match block.content_id {
        Some(content_id) => write!(writer, "BLOCK {content_id:>8}"),
        None => write!(writer, "BLOCK {:>8}", "CORE"),
    }?;
    write!(writer, " {uncompressed:>12} {compressed:>12}")
}

fn write_ratio(writer: &mut impl Write, compressed: u64, uncompressed: u64) -> io::Result<()> {
    let ratio = 100.0 * (compressed as f64 + 0.0001) / (uncompressed as f64 + 0.0001);
    if ratio > 999.0 {
        write!(writer, "   >999%")
    } else {
        write!(writer, " {ratio:6.2}%")
    }
}

fn write_series(
    writer: &mut impl Write,
    block: &BlockSummary,
    include_embedded_reference: bool,
) -> io::Result<()> {
    for data_series in &block.data_series {
        write!(writer, " {data_series}")?;
    }
    if include_embedded_reference && block.embedded_reference {
        write!(writer, " embedded_ref")?;
    }
    Ok(())
}
