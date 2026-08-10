use std::io::Write;
use std::path::Path;

use noodles::sam;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::{Program, input, md, output};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Format {
    Sam,
    Bam,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Compression {
    #[default]
    Default,
    Uncompressed,
}

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub format: Format,
    pub compression: Compression,
    pub use_equal: bool,
    pub additional_threads: usize,
    pub program: Option<Program<'a>>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub records_read: u64,
    pub records_recalculated: u64,
    pub records_with_corrected_tags: u64,
    pub records_without_sequence: u64,
}

pub fn write<W>(
    input_path: &Path,
    reference_path: &Path,
    options: Options<'_>,
    output: W,
) -> Result<Summary>
where
    W: Write + Send + 'static,
{
    if options.additional_threads > 256 {
        return Err(RsomicsError::ConfigError(
            "calmd additional thread count cannot exceed 256".to_owned(),
        ));
    }

    let mut reader = input::open(input_path, Some(reference_path), options.additional_threads)?;
    let input_format = reader.format();
    let input_header = reader.read_header(input_path)?;
    let mut output_header = input_header.clone();
    if let Some(program) = options.program {
        program.add_to(&mut output_header)?;
    }

    let output_format = match options.format {
        Format::Sam => output::Format::Sam,
        Format::Bam => output::Format::Bam,
    };
    let compression = match options.compression {
        Compression::Default => output::Compression::Default,
        Compression::Uncompressed => output::Compression::Uncompressed,
    };
    let mut writer = output::Writer::new(
        output_format,
        compression,
        options.additional_threads,
        output,
    );
    writer.write_header(&output_header)?;

    let mut reference = md::ReferenceCache::open(reference_path)?;
    let mut summary = Summary::default();
    if input_format == input::Format::Bam && options.format == Format::Bam {
        let mut cigar = Vec::new();
        let mut mismatch_positions = Vec::new();
        reader.visit_mut_raw_bam_records(input_path, |record| {
            let result = md::recalculate_raw(
                &input_header,
                record,
                &mut reference,
                options.use_equal,
                &mut cigar,
                &mut mismatch_positions,
            )?;
            update_summary(&mut summary, result)?;
            writer.write_owned_raw_record(record)?;
            Ok(true)
        })?;
    } else {
        reader.visit_records(&input_header, input_path, |record| {
            let mut record = if input_format == input::Format::Cram {
                md::complete(&input_header, record, Some(&mut reference))?
            } else {
                sam::alignment::RecordBuf::try_from_alignment_record(&input_header, record)
                    .map_err(RsomicsError::Io)?
            };
            let result = md::recalculate_record(
                &input_header,
                &mut record,
                &mut reference,
                options.use_equal,
            )?;
            update_summary(&mut summary, result)?;
            writer.write_record(&output_header, &record)?;
            Ok(true)
        })?;
    }
    writer.finish(&output_header)?;
    Ok(summary)
}

fn update_summary(summary: &mut Summary, result: md::Recalculation) -> Result<()> {
    summary.records_read = summary
        .records_read
        .checked_add(1)
        .ok_or_else(summary_overflow)?;
    match result {
        md::Recalculation::Unchanged => {}
        md::Recalculation::MissingSequence => {
            summary.records_without_sequence = summary
                .records_without_sequence
                .checked_add(1)
                .ok_or_else(summary_overflow)?;
        }
        md::Recalculation::Updated { corrected_existing } => {
            summary.records_recalculated = summary
                .records_recalculated
                .checked_add(1)
                .ok_or_else(summary_overflow)?;
            if corrected_existing {
                summary.records_with_corrected_tags = summary
                    .records_with_corrected_tags
                    .checked_add(1)
                    .ok_or_else(summary_overflow)?;
            }
        }
    }
    Ok(())
}

fn summary_overflow() -> RsomicsError {
    RsomicsError::InvalidInput("calmd record count overflow".to_owned())
}
