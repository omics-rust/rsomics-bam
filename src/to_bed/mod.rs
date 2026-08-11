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
    pub bed12: bool,
    pub split: bool,
    pub split_deletions: bool,
    pub color: &'a str,
    pub score_tag: Option<[u8; 2]>,
    pub cigar: bool,
    pub bedpe: bool,
    pub mate1_first: bool,
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
    let mut reader = input::open(input_path, options.reference, options.additional_threads)?;
    let header = reader.read_header(input_path)?;
    let references: Vec<_> = header
        .reference_sequences()
        .keys()
        .map(ToString::to_string)
        .collect();
    let mut output = BufWriter::with_capacity(256 * 1024, output);
    let mut summary = Summary {
        format: if options.bedpe {
            Format::Bedpe
        } else if options.bed12 {
            Format::Bed12
        } else {
            Format::Bed6
        },
        ..Summary::default()
    };

    if options.bedpe {
        let mut pairs = pair::State::new();
        if reader.has_reusable_raw_bam_path() {
            reader.visit_mut_raw_bam_records(input_path, |raw| {
                visit_pair(
                    raw,
                    &references,
                    options,
                    &mut output,
                    &mut pairs,
                    &mut summary,
                )
            })?;
        } else {
            reader.visit_owned_raw_records(&header, input_path, |raw| {
                visit_pair(
                    &raw,
                    &references,
                    options,
                    &mut output,
                    &mut pairs,
                    &mut summary,
                )
            })?;
        }
        pairs.finish()?;
        output.flush().map_err(RsomicsError::Io)?;
        return Ok(summary);
    }

    let mut cigar_ops = Vec::new();
    if reader.has_reusable_raw_bam_path() {
        reader.visit_mut_raw_bam_records(input_path, |raw| {
            visit_record(
                raw,
                &references,
                options,
                &mut output,
                &mut cigar_ops,
                &mut summary,
            )
        })?;
    } else {
        reader.visit_owned_raw_records(&header, input_path, |raw| {
            visit_record(
                &raw,
                &references,
                options,
                &mut output,
                &mut cigar_ops,
                &mut summary,
            )
        })?;
    }
    output.flush().map_err(RsomicsError::Io)?;
    Ok(summary)
}

fn visit_pair(
    raw: &RawRecord,
    references: &[String],
    options: Options<'_>,
    output: &mut impl Write,
    pairs: &mut pair::State,
    summary: &mut Summary,
) -> Result<bool> {
    summary.records_read = increment(summary.records_read)?;
    if raw.flags() & record::UNMAPPED == 0 {
        summary.records_mapped = increment(summary.records_mapped)?;
    }
    if pairs.push(output, references, raw, options)? {
        summary.pairs_written = increment(summary.pairs_written)?;
        summary.rows_written = increment(summary.rows_written)?;
    }
    Ok(true)
}

fn visit_record(
    raw: &RawRecord,
    references: &[String],
    options: Options<'_>,
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
    let score = record::score(raw, options.score_tag)?;
    let cigar = options
        .cigar
        .then(|| record::cigar_text(cigar_ops, mapped.name))
        .transpose()?;
    if options.bed12 {
        let blocks = record::blocks(cigar_ops, mapped.start, options.split_deletions)?;
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
            options.color,
            &blocks,
        )?;
        summary.rows_written = increment(summary.rows_written)?;
    } else if options.split {
        for (block_start, block_end) in
            record::blocks(cigar_ops, mapped.start, options.split_deletions)?
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
    } else {
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
    Ok(true)
}

fn increment(value: u64) -> Result<u64> {
    value
        .checked_add(1)
        .ok_or_else(|| RsomicsError::InvalidInput("record count overflows".to_owned()))
}
