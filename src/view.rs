use std::io::Write;
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::{filter::Filter, input, md, output};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Format {
    #[default]
    Sam,
    Bam,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Options<'a> {
    pub with_header: bool,
    pub header_only: bool,
    pub count_only: bool,
    pub require_flags: u16,
    pub exclude_flags: u16,
    pub include_flags: u16,
    pub exclude_all_flags: u16,
    pub minimum_mapping_quality: u8,
    pub output_format: Format,
    pub reference: Option<&'a Path>,
    pub additional_threads: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub selected: u64,
    pub rejected: u64,
}

pub fn write(input_path: &Path, options: Options<'_>, mut output: impl Write) -> Result<Summary> {
    let mut reader = input::open(input_path, options.reference, options.additional_threads)?;
    let header = reader.read_header(input_path)?;
    let format = reader.format();

    let filter = Filter {
        require_all: options.require_flags,
        exclude_any: options.exclude_flags,
        include_any: options.include_flags,
        exclude_all: options.exclude_all_flags,
        minimum_mapping_quality: options.minimum_mapping_quality,
    };
    let mut selected = 0u64;
    let mut rejected = 0u64;

    if options.count_only {
        reader.visit_records(&header, input_path, |record| {
            if filter.accepts(record)? {
                selected = selected.checked_add(1).ok_or_else(count_overflow)?;
            } else {
                rejected = rejected.checked_add(1).ok_or_else(count_overflow)?;
            }
            Ok(true)
        })?;
    } else {
        let output_format = match options.output_format {
            Format::Sam => output::Format::Sam,
            Format::Bam => output::Format::Bam,
        };
        let mut writer = output::Writer::new(output_format, &mut output);
        if options.with_header || options.header_only || options.output_format != Format::Sam {
            writer.write_header(&header)?;
        }

        let mut reference = if format == input::Format::Cram {
            options
                .reference
                .map(md::ReferenceCache::open)
                .transpose()?
        } else {
            None
        };

        if !options.header_only {
            reader.visit_records(&header, input_path, |record| {
                if filter.accepts(record)? {
                    selected = selected.checked_add(1).ok_or_else(count_overflow)?;
                    if format == input::Format::Cram {
                        let record = md::complete(&header, record, reference.as_mut())?;
                        writer.write_record(&header, &record)?;
                    } else {
                        writer.write_record(&header, record)?;
                    }
                } else {
                    rejected = rejected.checked_add(1).ok_or_else(count_overflow)?;
                }
                Ok(true)
            })?;
        }
        writer.finish(&header)?;
    }

    Ok(Summary { selected, rejected })
}

fn count_overflow() -> RsomicsError {
    RsomicsError::InvalidInput("alignment record count exceeds u64".to_owned())
}
