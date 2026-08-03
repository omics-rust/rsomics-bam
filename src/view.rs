use std::io::Write;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use noodles::core;
use noodles::sam;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::{
    Program,
    filter::{Filter, LibraryFilter, QnameFilter, ReadGroupFilter},
    input, md, output,
};

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
    pub read_groups: &'a [String],
    pub qname_files: &'a [PathBuf],
    pub library: Option<&'a str>,
    pub minimum_mapping_quality: u8,
    pub minimum_query_length: u64,
    pub add_flags: u16,
    pub remove_flags: u16,
    pub remove_tags: Option<&'a [[u8; 2]]>,
    pub keep_tags: Option<&'a [[u8; 2]]>,
    pub output_format: Format,
    pub compression: Compression,
    pub reference: Option<&'a Path>,
    pub additional_threads: usize,
    pub regions: &'a [Region],
    pub program: Option<Program<'a>>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub selected: u64,
    pub rejected: u64,
}

pub fn write<W>(input_path: &Path, options: Options<'_>, output: W) -> Result<Summary>
where
    W: Write + Send + 'static,
{
    if options.remove_tags.is_some() && options.keep_tags.is_some() {
        return Err(RsomicsError::ConfigError(
            "--remove-tag and --keep-tag are mutually exclusive".to_owned(),
        ));
    }
    let parallel_output = !options.count_only
        && options.output_format == Format::Bam
        && options.additional_threads > 0;
    if !options.regions.is_empty()
        && options.additional_threads > 0
        && !options.header_only
        && !parallel_output
    {
        return Err(RsomicsError::ConfigError(
            "additional decoding threads are not available for indexed region queries yet"
                .to_owned(),
        ));
    }
    let qnames = QnameFilter::from_files(options.qname_files)?;

    let input_threads = if parallel_output {
        0
    } else {
        options.additional_threads
    };
    let mut reader = if options.regions.is_empty() || options.header_only {
        input::open(input_path, options.reference, input_threads)?
    } else {
        input::open_indexed(input_path, options.reference)?
    };
    let header = reader.read_header(input_path)?;
    let mut output_header = header.clone();
    let format = reader.format();

    let read_groups = ReadGroupFilter::new(options.read_groups);
    let library = LibraryFilter::new(&header, options.library);
    if let Some(read_groups) = &read_groups {
        read_groups.retain_header(&mut output_header);
    }
    let filter = Filter {
        require_all: options.require_flags,
        exclude_any: options.exclude_flags,
        include_any: options.include_flags,
        exclude_all: options.exclude_all_flags,
        read_groups: read_groups.as_ref(),
        qnames: qnames.as_ref(),
        library: library.as_ref(),
        minimum_mapping_quality: options.minimum_mapping_quality,
        minimum_query_length: options.minimum_query_length,
    };
    let mut selected = 0u64;
    let mut rejected = 0u64;

    if options.count_only {
        if format == input::Format::Bam && options.regions.is_empty() {
            reader.visit_raw_bam_records(input_path, |record| {
                if filter.accepts_raw(&record)? {
                    selected = selected.checked_add(1).ok_or_else(count_overflow)?;
                } else {
                    rejected = rejected.checked_add(1).ok_or_else(count_overflow)?;
                }
                Ok(true)
            })?;
        } else {
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
        }
    } else {
        if let Some(program) = options.program {
            program.add_to(&mut output_header)?;
        }
        let output_format = match options.output_format {
            Format::Sam => output::Format::Sam,
            Format::Bam => output::Format::Bam,
        };
        let compression = match options.compression {
            Compression::Default => output::Compression::Default,
            Compression::Fast => output::Compression::Fast,
            Compression::Uncompressed => output::Compression::Uncompressed,
        };
        let mut writer = output::Writer::new(
            output_format,
            compression,
            options.additional_threads,
            output,
        );
        if options.with_header || options.header_only || options.output_format != Format::Sam {
            writer.write_header(&output_header)?;
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
            if format == input::Format::Bam
                && options.output_format == Format::Bam
                && options.regions.is_empty()
                && options.add_flags == 0
                && options.remove_flags == 0
                && options.remove_tags.is_none()
                && options.keep_tags.is_none()
            {
                reader.visit_raw_bam_records(input_path, |record| {
                    if filter.accepts_raw(&record)? {
                        selected = selected.checked_add(1).ok_or_else(count_overflow)?;
                        writer.write_raw_record(&record)?;
                    } else {
                        rejected = rejected.checked_add(1).ok_or_else(count_overflow)?;
                    }
                    Ok(true)
                })?;
            } else {
                visit_records(
                    &mut reader,
                    &header,
                    input_path,
                    options.regions,
                    |record| {
                        if filter.accepts(record)? {
                            selected = selected.checked_add(1).ok_or_else(count_overflow)?;
                            if format == input::Format::Cram {
                                let mut record = md::complete(&header, record, reference.as_mut())?;
                                transform_record(&mut record, options);
                                writer.write_record(&output_header, &record)?;
                            } else if has_record_transform(options) {
                                let mut record =
                                    sam::alignment::RecordBuf::try_from_alignment_record(
                                        &header, record,
                                    )
                                    .map_err(RsomicsError::Io)?;
                                transform_record(&mut record, options);
                                writer.write_record(&output_header, &record)?;
                            } else {
                                writer.write_record(&output_header, record)?;
                            }
                        } else {
                            rejected = rejected.checked_add(1).ok_or_else(count_overflow)?;
                        }
                        Ok(true)
                    },
                )?;
            }
        }
        writer.finish(&output_header)?;
    }

    Ok(Summary { selected, rejected })
}

fn has_record_transform(options: Options<'_>) -> bool {
    options.add_flags != 0
        || options.remove_flags != 0
        || options.remove_tags.is_some()
        || options.keep_tags.is_some()
}

fn transform_record(record: &mut sam::alignment::RecordBuf, options: Options<'_>) {
    let flags = (u16::from(record.flags()) | options.add_flags) & !options.remove_flags;
    *record.flags_mut() = sam::alignment::record::Flags::from_bits_retain(flags);

    let tags = match (options.keep_tags, options.remove_tags) {
        (Some(tags), None) => Some((tags, true)),
        (None, Some(tags)) => Some((tags, false)),
        _ => None,
    };
    let Some((tags, keep)) = tags else {
        return;
    };

    let data = record.data_mut();
    let removed = data
        .keys()
        .filter(|tag| {
            let present = tags.contains(tag.as_ref());
            present != keep
        })
        .collect::<Vec<_>>();
    for tag in removed {
        data.remove(&tag);
    }
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
