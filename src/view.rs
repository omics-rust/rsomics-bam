use std::io::Write;
use std::path::Path;
use std::str::FromStr;

use noodles::core;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::{filter::Filter, input, md, output};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Format {
    #[default]
    Sam,
    Bam,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Compression {
    #[default]
    Default,
    Fast,
    Uncompressed,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Region {
    Mapped(core::Region),
    Unmapped,
}

impl FromStr for Region {
    type Err = String;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        if value == "*" {
            Ok(Self::Unmapped)
        } else {
            value
                .parse()
                .map(Self::Mapped)
                .map_err(|error| format!("invalid region {value}: {error}"))
        }
    }
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
    pub compression: Compression,
    pub reference: Option<&'a Path>,
    pub additional_threads: usize,
    pub regions: &'a [Region],
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub selected: u64,
    pub rejected: u64,
}

pub fn write(input_path: &Path, options: Options<'_>, mut output: impl Write) -> Result<Summary> {
    if !options.regions.is_empty() && options.additional_threads > 0 && !options.header_only {
        return Err(RsomicsError::ConfigError(
            "additional decoding threads are not available for indexed region queries yet"
                .to_owned(),
        ));
    }

    let mut reader = if options.regions.is_empty() || options.header_only {
        input::open(input_path, options.reference, options.additional_threads)?
    } else {
        input::open_indexed(input_path, options.reference)?
    };
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
        visit_records(
            &mut reader,
            &header,
            input_path,
            options.regions,
            |record| {
                if filter.accepts(record)? {
                    selected = selected.checked_add(1).ok_or_else(count_overflow)?;
                } else {
                    rejected = rejected.checked_add(1).ok_or_else(count_overflow)?;
                }
                Ok(true)
            },
        )?;
    } else {
        let output_format = match options.output_format {
            Format::Sam => output::Format::Sam,
            Format::Bam => output::Format::Bam,
        };
        let compression = match options.compression {
            Compression::Default => output::Compression::Default,
            Compression::Fast => output::Compression::Fast,
            Compression::Uncompressed => output::Compression::Uncompressed,
        };
        let mut writer = output::Writer::new(output_format, compression, &mut output);
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
            visit_records(
                &mut reader,
                &header,
                input_path,
                options.regions,
                |record| {
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
                },
            )?;
        }
        writer.finish(&header)?;
    }

    Ok(Summary { selected, rejected })
}

fn visit_records(
    reader: &mut input::Reader,
    header: &noodles::sam::Header,
    input_path: &Path,
    regions: &[Region],
    mut visit: impl FnMut(&dyn noodles::sam::alignment::Record) -> Result<bool>,
) -> Result<()> {
    if regions.is_empty() {
        return reader.visit_records(header, input_path, visit);
    }

    for region in regions {
        match region {
            Region::Mapped(region) => {
                reader.visit_region(header, input_path, Some(region), &mut visit)?;
            }
            Region::Unmapped => {
                reader.visit_region(header, input_path, None, &mut visit)?;
            }
        }
    }
    Ok(())
}

fn count_overflow() -> RsomicsError {
    RsomicsError::InvalidInput("alignment record count exceeds u64".to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_mapped_and_unmapped_regions() {
        assert_eq!("*".parse(), Ok(Region::Unmapped));
        assert!(matches!(
            "chr1:8-13".parse(),
            Ok(Region::Mapped(region)) if region.to_string() == "chr1:8-13"
        ));
        assert!("chr1:0".parse::<Region>().is_err());
    }
}
