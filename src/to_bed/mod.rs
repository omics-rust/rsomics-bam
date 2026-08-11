mod pair;
mod record;
mod render;

use std::io::{BufWriter, Write};
use std::path::Path;

use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::input;

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub reference: Option<&'a Path>,
    pub additional_threads: usize,
    pub layout: Layout,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Score {
    MappingQuality,
    EditDistance,
    Tag([u8; 2]),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PairScore {
    MappingQuality,
    EditDistance,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Layout {
    Records(RecordLayout),
    Bedpe { score: PairScore, mate1_first: bool },
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RecordLayout {
    Bed6 {
        score: Score,
        cigar: bool,
    },
    SplitBed6 {
        score: Score,
        split_deletions: bool,
    },
    Bed12 {
        score: Score,
        split_deletions: bool,
        color: [u8; 3],
    },
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub format: Format,
    pub records_read: u64,
    pub records_mapped: u64,
    pub records_skipped: u64,
    pub pairs_written: u64,
    pub rows_written: u64,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Format {
    #[default]
    Bed6,
    Bed12,
    Bedpe,
}

pub fn write<W: Write>(input_path: &Path, options: Options<'_>, output: W) -> Result<Summary> {
    if options.additional_threads > 256 {
        return Err(RsomicsError::ConfigError(
            "to-bed additional thread count cannot exceed 256".to_owned(),
        ));
    }
    let mut reader = input::open(input_path, options.reference, options.additional_threads)?;
    let header = reader.read_header(input_path)?;
    let references: Vec<_> = header
        .reference_sequences()
        .keys()
        .map(ToString::to_string)
        .collect();
    let mut output = BufWriter::with_capacity(256 * 1024, output);
    let mut summary = Summary {
        format: match options.layout {
            Layout::Records(RecordLayout::Bed6 { .. } | RecordLayout::SplitBed6 { .. }) => {
                Format::Bed6
            }
            Layout::Records(RecordLayout::Bed12 { .. }) => Format::Bed12,
            Layout::Bedpe { .. } => Format::Bedpe,
        },
        ..Summary::default()
    };
    let mut context = WriteContext {
        input_path,
        references: &references,
        output: &mut output,
        summary: &mut summary,
    };

    match options.layout {
        Layout::Bedpe { score, mate1_first } => {
            write_pairs(&mut reader, &header, &mut context, score, mate1_first)?
        }
        Layout::Records(layout) => write_records(&mut reader, &header, &mut context, layout)?,
    }
    output.flush().map_err(RsomicsError::Io)?;
    Ok(summary)
}

struct WriteContext<'a, W> {
    input_path: &'a Path,
    references: &'a [String],
    output: &'a mut W,
    summary: &'a mut Summary,
}

fn write_pairs<W: Write>(
    reader: &mut input::Reader,
    header: &noodles::sam::Header,
    context: &mut WriteContext<'_, W>,
    score: PairScore,
    mate1_first: bool,
) -> Result<()> {
    let mut pairs = pair::State::new();
    if reader.has_reusable_raw_bam_path() {
        reader.visit_mut_raw_bam_records(context.input_path, |raw| {
            visit_pair(
                raw,
                context.references,
                score,
                mate1_first,
                context.output,
                &mut pairs,
                context.summary,
            )
        })?;
    } else {
        reader.visit_owned_raw_records(header, context.input_path, |raw| {
            visit_pair(
                &raw,
                context.references,
                score,
                mate1_first,
                context.output,
                &mut pairs,
                context.summary,
            )
        })?;
    }
    pairs.finish()
}

fn write_records<W: Write>(
    reader: &mut input::Reader,
    header: &noodles::sam::Header,
    context: &mut WriteContext<'_, W>,
    layout: RecordLayout,
) -> Result<()> {
    let mut cigar_ops = Vec::new();
    if reader.has_reusable_raw_bam_path() {
        reader.visit_mut_raw_bam_records(context.input_path, |raw| {
            visit_record(
                raw,
                context.references,
                layout,
                context.output,
                &mut cigar_ops,
                context.summary,
            )
        })?;
    } else {
        reader.visit_owned_raw_records(header, context.input_path, |raw| {
            visit_record(
                &raw,
                context.references,
                layout,
                context.output,
                &mut cigar_ops,
                context.summary,
            )
        })?;
    }
    Ok(())
}

fn visit_pair(
    raw: &RawRecord,
    references: &[String],
    score: PairScore,
    mate1_first: bool,
    output: &mut impl Write,
    pairs: &mut pair::State,
    summary: &mut Summary,
) -> Result<bool> {
    summary.records_read = increment(summary.records_read)?;
    if raw.flags() & record::UNMAPPED == 0 {
        summary.records_mapped = increment(summary.records_mapped)?;
    }
    if pairs.push(output, references, raw, score, mate1_first)? {
        summary.pairs_written = increment(summary.pairs_written)?;
        summary.rows_written = increment(summary.rows_written)?;
    }
    Ok(true)
}

fn visit_record(
    raw: &RawRecord,
    references: &[String],
    layout: RecordLayout,
    output: &mut impl Write,
    cigar_ops: &mut Vec<(u8, u32)>,
    summary: &mut Summary,
) -> Result<bool> {
    summary.records_read = increment(summary.records_read)?;
    let Some(mapped) = record::project(raw, references)? else {
        summary.records_skipped = increment(summary.records_skipped)?;
        return Ok(true);
    };
    summary.records_mapped = increment(summary.records_mapped)?;
    raw.decode_cigar_into(cigar_ops)?;
    match layout {
        RecordLayout::Bed6 { score, cigar } => {
            let score = record::score(raw, score)?;
            let cigar = cigar
                .then(|| record::cigar_text(cigar_ops, mapped.name))
                .transpose()?;
            let end = record::reference_end(cigar_ops, mapped.start, mapped.name)?;
            render::bed6(
                output,
                render::Bed {
                    reference: mapped.reference,
                    start: mapped.start,
                    end,
                    name: mapped.name,
                    flags: mapped.flags,
                    score,
                },
                cigar.as_deref(),
            )?;
            summary.rows_written = increment(summary.rows_written)?;
        }
        RecordLayout::SplitBed6 {
            score,
            split_deletions,
        } => {
            let score = record::score(raw, score)?;
            for (block_start, block_end) in
                record::blocks(cigar_ops, mapped.start, split_deletions)?
            {
                render::bed6(
                    output,
                    render::Bed {
                        reference: mapped.reference,
                        start: block_start,
                        end: block_end,
                        name: mapped.name,
                        flags: mapped.flags,
                        score,
                    },
                    None,
                )?;
                summary.rows_written = increment(summary.rows_written)?;
            }
        }
        RecordLayout::Bed12 {
            score,
            split_deletions,
            color,
        } => {
            let score = record::score(raw, score)?;
            let blocks = record::blocks(cigar_ops, mapped.start, split_deletions)?;
            render::bed12(
                output,
                render::Bed {
                    reference: mapped.reference,
                    start: mapped.start,
                    end: record::reference_end(cigar_ops, mapped.start, mapped.name)?,
                    name: mapped.name,
                    flags: mapped.flags,
                    score,
                },
                color,
                &blocks,
            )?;
            summary.rows_written = increment(summary.rows_written)?;
        }
    }
    Ok(true)
}

fn increment(value: u64) -> Result<u64> {
    value
        .checked_add(1)
        .ok_or_else(|| RsomicsError::InvalidInput("record count overflows".to_owned()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writer_rejects_excess_workers_before_opening_input() {
        let error = write(
            Path::new("missing.sam"),
            Options {
                reference: None,
                additional_threads: 257,
                layout: Layout::Records(RecordLayout::Bed6 {
                    score: Score::MappingQuality,
                    cigar: false,
                }),
            },
            Vec::new(),
        )
        .unwrap_err();

        assert!(error.to_string().contains("cannot exceed 256"), "{error}");
    }
}
