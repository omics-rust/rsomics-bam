use std::path::{Path, PathBuf};

use noodles::sam;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::{Program, input};

mod label;
mod output;
mod tag;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    Sam,
    Bam,
    Cram,
}

impl Format {
    fn extension(self) -> &'static str {
        match self {
            Self::Sam => "sam",
            Self::Bam => "bam",
            Self::Cram => "cram",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Mode {
    ReadGroup,
    Tag([u8; 2]),
}

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub mode: Mode,
    pub output_prefix: &'a Path,
    pub unaccounted: Option<&'a Path>,
    pub unaccounted_header: Option<&'a Path>,
    pub format: Format,
    pub maximum_outputs: usize,
    pub zero_pad: usize,
    pub reference: Option<&'a Path>,
    pub additional_threads: usize,
    pub program: Option<Program<'a>>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct OutputSummary {
    pub label: String,
    pub path: PathBuf,
    pub records: u64,
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub records: u64,
    pub outputs: Vec<OutputSummary>,
    pub skipped: u64,
}

pub fn run(input_path: &Path, options: Options<'_>) -> Result<Summary> {
    validate(input_path, options)?;
    let mut reader = input::open(input_path, options.reference, options.additional_threads)?;
    let mut header = reader.read_header(input_path)?;
    let mut unaccounted_header = options
        .unaccounted_header
        .map(|path| read_header(path, options.reference))
        .transpose()?;
    if unaccounted_header
        .as_ref()
        .is_some_and(|candidate| !same_dictionary(&header, candidate))
    {
        return Err(RsomicsError::InvalidInput(
            "unaccounted header reference dictionary differs from the input".to_owned(),
        ));
    }
    if let Some(program) = options.program {
        program.add_to(&mut header)?;
        if let Some(candidate) = &mut unaccounted_header {
            program.add_to(candidate)?;
        }
    }

    let mut router = output::Router::new(input_path, &header, unaccounted_header, options)?;
    let tag = match options.mode {
        Mode::ReadGroup => *b"RG",
        Mode::Tag(tag) => tag,
    };
    let require_string = options.mode == Mode::ReadGroup || tag == *b"RG";
    let mut records = 0u64;
    if options.format == Format::Bam {
        reader.visit_owned_raw_records(&header, input_path, |record| {
            let outcome = if require_string {
                tag::read_string(&record, tag)
            } else {
                tag::read(&record, tag, options.zero_pad)
            };
            router.write_raw(&record, outcome)?;
            records = increment(records)?;
            Ok(true)
        })?;
    } else {
        reader.visit_records(&header, input_path, |record| {
            let record = sam::alignment::RecordBuf::try_from_alignment_record(&header, record)
                .map_err(RsomicsError::Io)?;
            let outcome = if require_string {
                tag::read_string_record(&record, tag)?
            } else {
                tag::read_record(&record, tag, options.zero_pad)?
            };
            router.write_record(&record, outcome)?;
            records = increment(records)?;
            Ok(true)
        })?;
    }
    let mut summary = router.finish()?;
    summary.records = records;
    Ok(summary)
}

fn validate(input_path: &Path, options: Options<'_>) -> Result<()> {
    if options.additional_threads > 256 {
        return Err(RsomicsError::ConfigError(
            "split additional thread count cannot exceed 256".to_owned(),
        ));
    }
    if options.maximum_outputs == 0 {
        return Err(RsomicsError::ConfigError(
            "split maximum output count must be positive".to_owned(),
        ));
    }
    if options.zero_pad > 128 {
        return Err(RsomicsError::ConfigError(
            "split integer padding cannot exceed 128".to_owned(),
        ));
    }
    if input_path == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "split requires a file-backed input".to_owned(),
        ));
    }
    if options.format == Format::Cram && options.reference.is_none() {
        return Err(RsomicsError::ConfigError(
            "CRAM split output requires an indexed reference".to_owned(),
        ));
    }
    if options.unaccounted_header.is_some() && options.unaccounted.is_none() {
        return Err(RsomicsError::ConfigError(
            "an unaccounted header requires an unaccounted output".to_owned(),
        ));
    }
    Ok(())
}

fn increment(value: u64) -> Result<u64> {
    value
        .checked_add(1)
        .ok_or_else(|| RsomicsError::InvalidInput("split record count exceeds u64".to_owned()))
}

fn read_header(path: &Path, reference: Option<&Path>) -> Result<sam::Header> {
    let mut reader = input::open(path, reference, 0)?;
    reader.read_header(path)
}

fn same_dictionary(left: &sam::Header, right: &sam::Header) -> bool {
    left.reference_sequences().len() == right.reference_sequences().len()
        && left
            .reference_sequences()
            .iter()
            .zip(right.reference_sequences())
            .all(|((left_name, left_map), (right_name, right_map))| {
                left_name == right_name && left_map.length() == right_map.length()
            })
}
