use std::io::{BufWriter, Write};
use std::path::Path;
use std::str::FromStr;

use noodles::core::Region;
use noodles::sam;
use rsomics_common::{Result, RsomicsError};
use rust_htslib::bam::{self, Read as _};
use rust_htslib::htslib;
use serde::Serialize;

use crate::input;

mod builder;
mod embedded;
mod md;

use builder::{Builder, Reference, Selection};

#[derive(Clone, Debug, Default)]
/// Reference reconstruction configuration.
pub struct Options {
    /// Read embedded reference blocks from CRAM instead of alignment MD fields.
    pub embedded: bool,
    /// Optional indexed region using samtools region syntax.
    pub region: Option<String>,
    /// Additional alignment I/O workers.
    pub additional_threads: usize,
}

#[derive(Clone, Debug, Serialize)]
/// Recovery statistics for one FASTA record.
pub struct ReferenceSummary {
    /// FASTA record name.
    pub name: String,
    /// Number of emitted bases.
    pub bases: u64,
    /// Number of emitted bases other than `N`.
    pub known_bases: u64,
    /// Percentage of emitted bases other than `N`.
    pub coverage: f64,
}

#[derive(Clone, Debug, Default, Serialize)]
/// Aggregate reference recovery statistics.
pub struct Summary {
    /// Number of emitted FASTA records.
    pub references: u64,
    /// Total number of emitted bases.
    pub bases: u64,
    /// Total number of emitted bases other than `N`.
    pub known_bases: u64,
    /// Per-record statistics in output order.
    pub items: Vec<ReferenceSummary>,
}

/// Recover reference FASTA into `output` from a SAM, BAM, or CRAM input.
pub fn write(input_path: &Path, options: Options, output: impl Write) -> Result<Summary> {
    if options.embedded {
        return embedded::write(
            input_path,
            options.region.as_deref(),
            options.additional_threads,
            output,
        );
    }
    let region = options
        .region
        .as_deref()
        .map(Region::from_str)
        .transpose()
        .map_err(|error| RsomicsError::ConfigError(format!("invalid region: {error}")))?;
    let use_hts =
        input_path == Path::new("-") || input::detect_format(input_path)? == input::Format::Cram;
    if use_hts {
        return write_hts(
            input_path,
            region.as_ref(),
            options.additional_threads,
            output,
        );
    }
    let mut reader = if region.is_some() {
        input::open_indexed(input_path, None)?
    } else {
        input::open(input_path, None, options.additional_threads)?
    };
    let header = reader.read_header(input_path)?;
    let references = sam_references(&header);
    let selection = region
        .as_ref()
        .map(|region| Selection::new(&references, region))
        .transpose()?;
    let mut builder = Builder::new(references, selection, BufWriter::new(output));

    if let Some(region) = region.as_ref() {
        reader
            .visit_owned_raw_region(&header, input_path, region, |record| builder.add(&record))?;
    } else if reader.has_reusable_raw_bam_path() {
        reader.visit_raw_bam_records(input_path, |record| builder.add(&record))?;
    } else {
        reader.visit_owned_raw_records(&header, input_path, |record| builder.add(&record))?;
    }
    builder.finish()
}

fn write_hts(
    input_path: &Path,
    region: Option<&Region>,
    additional_threads: usize,
    output: impl Write,
) -> Result<Summary> {
    if let Some(region) = region {
        let mut reader = bam::IndexedReader::from_path(input_path)
            .map_err(|error| hts_error("opening indexed alignment", input_path, error))?;
        configure_hts_reader(&mut reader, additional_threads, input_path)?;
        let references = hts_references(reader.header())?;
        let selection = Selection::new(&references, region)?;
        reader
            .fetch((
                i32::try_from(selection.reference_id).unwrap(),
                i64::try_from(selection.start).unwrap(),
                i64::try_from(selection.end).unwrap(),
            ))
            .map_err(|error| hts_error("querying alignment region", input_path, error))?;
        let mut builder = Builder::new(references, Some(selection), BufWriter::new(output));
        for result in reader.records() {
            let record =
                result.map_err(|error| hts_error("reading alignment record", input_path, error))?;
            builder.add(&record)?;
        }
        builder.finish()
    } else {
        let mut reader = if input_path == Path::new("-") {
            bam::Reader::from_stdin()
        } else {
            bam::Reader::from_path(input_path)
        }
        .map_err(|error| hts_error("opening alignment", input_path, error))?;
        configure_hts_reader(&mut reader, additional_threads, input_path)?;
        let references = hts_references(reader.header())?;
        let mut builder = Builder::new(references, None, BufWriter::new(output));
        for result in reader.records() {
            let record =
                result.map_err(|error| hts_error("reading alignment record", input_path, error))?;
            builder.add(&record)?;
        }
        builder.finish()
    }
}

fn configure_hts_reader(
    reader: &mut impl bam::Read,
    additional_threads: usize,
    input_path: &Path,
) -> Result<()> {
    if additional_threads > 0 {
        reader
            .set_threads(additional_threads)
            .map_err(|error| hts_error("configuring alignment threads", input_path, error))?;
    }
    let fields = htslib::sam_fields_SAM_FLAG
        | htslib::sam_fields_SAM_RNAME
        | htslib::sam_fields_SAM_POS
        | htslib::sam_fields_SAM_CIGAR
        | htslib::sam_fields_SAM_SEQ
        | htslib::sam_fields_SAM_AUX;
    reader
        .set_cram_options(htslib::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS, fields)
        .map_err(|error| hts_error("configuring CRAM fields", input_path, error))?;
    reader
        .set_cram_options(htslib::hts_fmt_option_CRAM_OPT_DECODE_MD, 1)
        .map_err(|error| hts_error("enabling CRAM MD decoding", input_path, error))?;
    Ok(())
}

fn sam_references(header: &sam::Header) -> Vec<Reference> {
    header
        .reference_sequences()
        .iter()
        .map(|(name, reference)| Reference {
            name: name.to_vec().into_boxed_slice(),
            length: usize::from(reference.length()),
        })
        .collect()
}

fn hts_references(header: &bam::HeaderView) -> Result<Vec<Reference>> {
    header
        .target_names()
        .iter()
        .enumerate()
        .map(|(reference_id, name)| {
            let length = header
                .target_len(u32::try_from(reference_id).unwrap())
                .ok_or_else(|| {
                    RsomicsError::InvalidInput(format!(
                        "reference {reference_id} has no declared length"
                    ))
                })?;
            Ok(Reference {
                name: name.to_vec().into_boxed_slice(),
                length: usize::try_from(length).map_err(|_| {
                    RsomicsError::InvalidInput(format!(
                        "reference {reference_id} length exceeds this platform"
                    ))
                })?,
            })
        })
        .collect()
}

fn hts_error(action: &str, input_path: &Path, error: rust_htslib::errors::Error) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{action} {}: {error}", input_path.display()))
}
