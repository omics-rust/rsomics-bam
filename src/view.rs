use std::io::Write;
use std::path::Path;
use std::str::FromStr;

use noodles::core;
use noodles::sam;
use noodles::sam::header::record::value::{
    Map,
    map::{Program as SamProgram, program::tag},
};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::{
    filter::{Filter, ReadGroupFilter},
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

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Program<'a> {
    name: &'a str,
    version: &'a str,
    command_line: &'a str,
}

impl<'a> Program<'a> {
    pub fn new(name: &'a str, version: &'a str, command_line: &'a str) -> Result<Self> {
        if name.is_empty() || version.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "program name and version cannot be empty".to_owned(),
            ));
        }
        if [name, version, command_line]
            .iter()
            .any(|value| value.contains(['\t', '\r', '\n']))
        {
            return Err(RsomicsError::InvalidInput(
                "program header fields cannot contain tabs or line breaks".to_owned(),
            ));
        }
        Ok(Self {
            name,
            version,
            command_line,
        })
    }
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
    pub minimum_mapping_quality: u8,
    pub minimum_query_length: u64,
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
    if let Some(read_groups) = &read_groups {
        read_groups.retain_header(&mut output_header);
    }
    let filter = Filter {
        require_all: options.require_flags,
        exclude_any: options.exclude_flags,
        include_any: options.include_flags,
        exclude_all: options.exclude_all_flags,
        read_groups: read_groups.as_ref(),
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
            add_program(&mut output_header, program)?;
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
                                let record = md::complete(&header, record, reference.as_mut())?;
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

fn add_program(header: &mut sam::Header, program: Program<'_>) -> Result<()> {
    let programs = header.programs_mut().as_mut();
    let previous = programs.last().map(|(id, _)| id.clone());
    let mut id = program.name.to_owned();
    let mut suffix = 0u64;

    while programs.contains_key(id.as_bytes()) {
        suffix = suffix.checked_add(1).ok_or_else(|| {
            RsomicsError::InvalidInput("program ID suffix exceeds u64".to_owned())
        })?;
        id = format!("{}.{}", program.name, suffix);
    }

    let mut builder = Map::<SamProgram>::builder()
        .insert(tag::NAME, program.name)
        .insert(tag::VERSION, program.version)
        .insert(tag::COMMAND_LINE, program.command_line);
    if let Some(previous) = previous {
        builder = builder.insert(tag::PREVIOUS_PROGRAM_ID, previous);
    }
    let map = builder
        .build()
        .map_err(|error| RsomicsError::InvalidInput(format!("building program record: {error}")))?;
    programs.insert(id.into(), map);
    Ok(())
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
    use noodles::sam::header::record::value::map::program::tag;

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

    #[test]
    fn adds_program_after_the_last_existing_record_with_a_unique_id() {
        let mut header = sam::Header::builder()
            .add_program("rsomics-bam", Map::default())
            .add_program("aligner", Map::default())
            .build();
        add_program(
            &mut header,
            Program::new("rsomics-bam", "1.2.3", "rsomics-bam view input.bam").unwrap(),
        )
        .unwrap();

        let program = &header.programs().as_ref()[b"rsomics-bam.1".as_slice()];
        assert_eq!(
            program
                .other_fields()
                .get(&tag::PREVIOUS_PROGRAM_ID)
                .map(|value| value.as_ref()),
            Some(b"aligner".as_slice())
        );
        assert_eq!(
            program
                .other_fields()
                .get(&tag::NAME)
                .map(|value| value.as_ref()),
            Some(b"rsomics-bam".as_slice())
        );
        assert_eq!(
            program
                .other_fields()
                .get(&tag::VERSION)
                .map(|value| value.as_ref()),
            Some(b"1.2.3".as_slice())
        );
    }

    #[test]
    fn program_fields_cannot_create_header_fields_or_lines() {
        assert!(Program::new("rsomics-bam", "1.2.3", "view\tinput.bam").is_err());
        assert!(Program::new("", "1.2.3", "view input.bam").is_err());
    }
}
