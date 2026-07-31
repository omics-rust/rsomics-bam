use std::io::Write;
use std::path::Path;

use noodles::sam::{self, alignment::io::Write as _};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::{filter::Filter, input, md};

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

    if options.with_header || options.header_only {
        let mut writer = sam::io::Writer::new(&mut output);
        writer.write_header(&header).map_err(RsomicsError::Io)?;
    }

    let filter = Filter {
        require_all: options.require_flags,
        exclude_any: options.exclude_flags,
        include_any: options.include_flags,
        exclude_all: options.exclude_all_flags,
        minimum_mapping_quality: options.minimum_mapping_quality,
    };
    let mut selected = 0u64;
    let mut rejected = 0u64;

    if !options.header_only {
        let mut reference = if format == input::Format::Cram && !options.count_only {
            options
                .reference
                .map(md::ReferenceCache::open)
                .transpose()?
        } else {
            None
        };
        let mut writer = sam::io::Writer::new(&mut output);

        reader.visit_records(&header, input_path, |record| {
            if filter.accepts(record)? {
                selected = selected.checked_add(1).ok_or_else(count_overflow)?;
                if !options.count_only {
                    if format == input::Format::Cram {
                        let record = md::complete(&header, record, reference.as_mut())?;
                        writer
                            .write_alignment_record(&header, &record)
                            .map_err(RsomicsError::Io)?;
                    } else {
                        writer
                            .write_alignment_record(&header, record)
                            .map_err(RsomicsError::Io)?;
                    }
                }
            } else {
                rejected = rejected.checked_add(1).ok_or_else(count_overflow)?;
            }
            Ok(true)
        })?;
    }

    output.flush().map_err(RsomicsError::Io)?;
    Ok(Summary { selected, rejected })
}

fn count_overflow() -> RsomicsError {
    RsomicsError::InvalidInput("alignment record count exceeds u64".to_owned())
}
