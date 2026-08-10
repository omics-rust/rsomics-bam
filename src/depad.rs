use std::fs::File;
use std::io::Write;
use std::num::NonZero;
use std::path::{Path, PathBuf};

use noodles::core::Region;
use noodles::fasta;
use noodles::sam::{
    self,
    alignment::{RecordBuf, record::MappingQuality},
};
use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::{Program, input, output};

mod bam_record;
mod cigar;

use self::bam_record::replace_raw_cigar;
use self::cigar::{
    CigarOp, MATCH, SKIP, build_position_map, decode_padded, decode_padded_raw, position,
    project_cigar, record_name, reference_code, reference_name, typed_cigar, typed_projected_cigar,
    unknown_reference,
};

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
    Fast,
    Uncompressed,
}

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub format: Format,
    pub compression: Compression,
    pub reference: Option<&'a Path>,
    pub additional_threads: usize,
    pub program: Option<Program<'a>>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub records_read: u64,
    pub records_projected: u64,
    pub records_preserved_unmapped: u64,
    pub embedded_references: u64,
    pub records_with_reference_skips: u64,
}

pub fn write<W>(input_path: &Path, options: Options<'_>, output: W) -> Result<Summary>
where
    W: Write + Send + 'static,
{
    if options.additional_threads > 256 {
        return Err(RsomicsError::ConfigError(
            "depad additional thread count cannot exceed 256".to_owned(),
        ));
    }

    let mut reader = input::open(input_path, None, 0)?;
    let input_format = reader.format();
    let input_header = reader.read_header(input_path)?;
    let (mut output_header, reference) = match options.reference {
        Some(path) => PaddedFasta::open(path, &input_header)?,
        None => (input_header.clone(), None),
    };
    if let Some(program) = options.program {
        program.add_to(&mut output_header)?;
    }

    let format = match options.format {
        Format::Sam => output::Format::Sam,
        Format::Bam => output::Format::Bam,
    };
    let compression = match options.compression {
        Compression::Default => output::Compression::Default,
        Compression::Fast => output::Compression::Fast,
        Compression::Uncompressed => output::Compression::Uncompressed,
    };
    let mut writer = output::Writer::new(format, compression, options.additional_threads, output);
    writer.write_header(&output_header)?;

    let mut projector = Projector::new(reference);
    let mut summary = Summary::default();
    if input_format == input::Format::Bam && options.format == Format::Bam {
        reader.visit_mut_raw_bam_records(input_path, |record| {
            let result = projector.project_raw(&input_header, record)?;
            update_summary(&mut summary, result)?;
            writer.write_owned_raw_record(record)?;
            Ok(true)
        })?;
    } else {
        reader.visit_records(&input_header, input_path, |record| {
            let mut record = RecordBuf::try_from_alignment_record(&input_header, record)
                .map_err(RsomicsError::Io)?;
            if input_format == input::Format::Cram && record.flags().is_unmapped() {
                normalize_unmapped_cram(&mut record);
            }
            let result = projector.project(&input_header, &mut record)?;
            update_summary(&mut summary, result)?;
            writer.write_record(&output_header, &record)?;
            Ok(true)
        })?;
    }
    writer.finish(&output_header)?;
    Ok(summary)
}

fn update_summary(summary: &mut Summary, result: ProjectionOutcome) -> Result<()> {
    summary.records_read = increment(summary.records_read)?;
    if result.reference_skip {
        summary.records_with_reference_skips = increment(summary.records_with_reference_skips)?;
    }
    match result.kind {
        Projection::EmbeddedReference => {
            summary.records_projected = increment(summary.records_projected)?;
            summary.embedded_references = increment(summary.embedded_references)?;
        }
        Projection::Projected => summary.records_projected = increment(summary.records_projected)?,
        Projection::Unmapped => {
            summary.records_preserved_unmapped = increment(summary.records_preserved_unmapped)?;
        }
    }
    Ok(())
}

fn normalize_unmapped_cram(record: &mut RecordBuf) {
    if record.mapping_quality().is_none() {
        *record.mapping_quality_mut() = Some(MappingQuality::MIN);
    }
    *record.mate_reference_sequence_id_mut() = None;
}

fn increment(value: u64) -> Result<u64> {
    value
        .checked_add(1)
        .ok_or_else(|| RsomicsError::InvalidInput("depad record count overflow".to_owned()))
}

struct PaddedFasta {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    path: PathBuf,
}

impl PaddedFasta {
    fn open(path: &Path, input_header: &sam::Header) -> Result<(sam::Header, Option<Self>)> {
        let index = fasta::fs::index(path).map_err(|error| {
            RsomicsError::ConfigError(format!(
                "indexing padded reference {}: {error}",
                path.display()
            ))
        })?;
        let mut reference = Self {
            reader: fasta::io::indexed_reader::Builder::default()
                .set_index(index)
                .build_from_path(path)
                .map_err(|error| {
                    RsomicsError::ConfigError(format!(
                        "opening padded reference {}: {error}",
                        path.display()
                    ))
                })?,
            path: path.to_path_buf(),
        };
        let mut output_header = input_header.clone();
        let mut sequence = Vec::new();

        for (name, map) in input_header.reference_sequences() {
            reference.read(name.as_ref(), usize::from(map.length()), &mut sequence)?;
            let unpadded = sequence.iter().filter(|&&base| base != 0).count();
            let length = NonZero::new(unpadded).ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "padded reference {} contains no bases",
                    String::from_utf8_lossy(name)
                ))
            })?;
            let output_map = output_header
                .reference_sequences_mut()
                .get_mut(name)
                .expect("output header is cloned from input header");
            *output_map.length_mut() = length;
        }

        Ok((output_header, Some(reference)))
    }

    fn read(&mut self, name: &[u8], expected: usize, sequence: &mut Vec<u8>) -> Result<()> {
        let record = self.reader.query(&Region::new(name, ..)).map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "reading padded reference {} from {}: {error}",
                String::from_utf8_lossy(name),
                self.path.display()
            ))
        })?;
        if record.sequence().len() != expected {
            return Err(RsomicsError::InvalidInput(format!(
                "padded reference {} has length {}, expected {expected}",
                String::from_utf8_lossy(name),
                record.sequence().len()
            )));
        }
        sequence.clear();
        sequence.reserve(expected);
        for &base in record.sequence().as_ref() {
            sequence.push(reference_code(base).ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "invalid base '{}' in padded reference {}",
                    char::from(base),
                    String::from_utf8_lossy(name)
                ))
            })?);
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Projection {
    EmbeddedReference,
    Projected,
    Unmapped,
}

struct ProjectionOutcome {
    kind: Projection,
    reference_skip: bool,
}

struct Projector {
    reference: Option<PaddedFasta>,
    reference_id: Option<usize>,
    sequence: Vec<u8>,
    position_map: Vec<usize>,
    padded_query: Vec<u8>,
    raw_cigar: Vec<(u8, u32)>,
    original_cigar: Vec<CigarOp>,
    projected_cigar: Vec<CigarOp>,
}

impl Projector {
    fn new(reference: Option<PaddedFasta>) -> Self {
        Self {
            reference,
            reference_id: None,
            sequence: Vec::new(),
            position_map: Vec::new(),
            padded_query: Vec::new(),
            raw_cigar: Vec::new(),
            original_cigar: Vec::new(),
            projected_cigar: Vec::new(),
        }
    }

    fn project(
        &mut self,
        header: &sam::Header,
        record: &mut RecordBuf,
    ) -> Result<ProjectionOutcome> {
        if record.flags().is_unmapped() {
            return Ok(ProjectionOutcome {
                kind: Projection::Unmapped,
                reference_skip: false,
            });
        }

        let name = record_name(record);
        let reference_id = record.reference_sequence_id().ok_or_else(|| {
            RsomicsError::InvalidInput(format!("mapped record {name} has no reference sequence"))
        })?;
        let start = record.alignment_start().ok_or_else(|| {
            RsomicsError::InvalidInput(format!("mapped record {name} has no alignment start"))
        })?;
        let padded_start = usize::from(start) - 1;
        let reference_name = reference_name(header, reference_id)?;
        let cigar = typed_cigar(record.cigar())?;
        if cigar.is_empty() {
            return Err(RsomicsError::InvalidInput(format!(
                "mapped record {name} has no CIGAR"
            )));
        }

        decode_padded(
            &cigar,
            record.sequence().as_ref(),
            &name,
            &mut self.padded_query,
        )?;
        let embedded =
            padded_start == 0 && record.name().is_some_and(|value| value == reference_name);
        if embedded {
            self.set_embedded(header, reference_id, &name)?;
            let length = record.sequence().len();
            let length = u32::try_from(length).map_err(|_| {
                RsomicsError::InvalidInput(format!("record {name} sequence is too long"))
            })?;
            self.projected_cigar.clear();
            self.projected_cigar.push(CigarOp::new(MATCH, length)?);
        } else {
            self.ensure_reference(header, reference_id)?;
            project_cigar(
                &cigar,
                &self.padded_query,
                &self.sequence,
                padded_start,
                &name,
                &mut self.projected_cigar,
            )?;
        }

        let new_start = self.mapped_position(header, reference_id, padded_start, &name)?;
        *record.alignment_start_mut() = Some(position(new_start, &name)?);
        *record.cigar_mut() = typed_projected_cigar(&self.projected_cigar);

        match (
            record.mate_reference_sequence_id(),
            record.mate_alignment_start(),
        ) {
            (Some(mate_reference_id), Some(mate_start)) => {
                let mate_start = usize::from(mate_start) - 1;
                let mapped = self.mapped_position(header, mate_reference_id, mate_start, &name)?;
                *record.mate_alignment_start_mut() = Some(position(mapped, &name)?);
                if mate_reference_id != reference_id {
                    self.ensure_reference(header, reference_id)?;
                }
            }
            _ => {
                *record.mate_reference_sequence_id_mut() = None;
                *record.mate_alignment_start_mut() = None;
            }
        }

        Ok(ProjectionOutcome {
            kind: if embedded {
                Projection::EmbeddedReference
            } else {
                Projection::Projected
            },
            reference_skip: cigar.iter().any(|op| op.kind == SKIP),
        })
    }

    fn project_raw(
        &mut self,
        header: &sam::Header,
        record: &mut RawRecord,
    ) -> Result<ProjectionOutcome> {
        if record.flags() & 0x04 != 0 {
            return Ok(ProjectionOutcome {
                kind: Projection::Unmapped,
                reference_skip: false,
            });
        }

        let name = String::from_utf8_lossy(record.name());
        let reference_id = usize::try_from(record.reference_sequence_id()).map_err(|_| {
            RsomicsError::InvalidInput(format!("mapped record {name} has no reference sequence"))
        })?;
        let padded_start = usize::try_from(record.alignment_start()).map_err(|_| {
            RsomicsError::InvalidInput(format!("mapped record {name} has no alignment start"))
        })?;
        let reference_name = reference_name(header, reference_id)?;
        record.decode_cigar_into(&mut self.raw_cigar)?;
        self.original_cigar.clear();
        self.original_cigar.reserve(self.raw_cigar.len());
        for &(kind, len) in &self.raw_cigar {
            self.original_cigar.push(CigarOp::new(kind, len)?);
        }
        if self.original_cigar.is_empty() {
            return Err(RsomicsError::InvalidInput(format!(
                "mapped record {name} has no CIGAR"
            )));
        }
        decode_padded_raw(&self.original_cigar, record, &name, &mut self.padded_query)?;

        let embedded = padded_start == 0 && record.name() == reference_name;
        if embedded {
            self.set_embedded(header, reference_id, &name)?;
            let length = u32::try_from(record.sequence_len()).map_err(|_| {
                RsomicsError::InvalidInput(format!("record {name} sequence is too long"))
            })?;
            self.projected_cigar.clear();
            self.projected_cigar.push(CigarOp::new(MATCH, length)?);
        } else {
            self.ensure_reference(header, reference_id)?;
            project_cigar(
                &self.original_cigar,
                &self.padded_query,
                &self.sequence,
                padded_start,
                &name,
                &mut self.projected_cigar,
            )?;
        }

        let new_start = self.mapped_position(header, reference_id, padded_start, &name)?;
        let new_start = i32::try_from(new_start).map_err(|_| {
            RsomicsError::InvalidInput(format!("record {name} position exceeds BAM limits"))
        })?;
        let mate_reference_id = record.mate_reference_sequence_id();
        let mate_start = record.mate_alignment_start();
        let new_mate_start = if mate_reference_id >= 0 && mate_start >= 0 {
            let mate_reference_id = usize::try_from(mate_reference_id).map_err(|_| {
                RsomicsError::InvalidInput(format!("record {name} has an invalid mate reference"))
            })?;
            let mate_start = usize::try_from(mate_start).map_err(|_| {
                RsomicsError::InvalidInput(format!("record {name} has an invalid mate position"))
            })?;
            let mapped = self.mapped_position(header, mate_reference_id, mate_start, &name)?;
            let mapped = i32::try_from(mapped).map_err(|_| {
                RsomicsError::InvalidInput(format!(
                    "record {name} mate position exceeds BAM limits"
                ))
            })?;
            if mate_reference_id != reference_id {
                self.ensure_reference(header, reference_id)?;
            }
            Some(mapped)
        } else {
            None
        };
        let source_long_cigar = self.raw_cigar.len() != record.cigar_ops().count();
        let reference_skip = self.original_cigar.iter().any(|op| op.kind == SKIP);
        drop(name);

        record.set_alignment_start(new_start);
        if let Some(mate_start) = new_mate_start {
            record.set_mate_alignment_start(mate_start);
        } else {
            record.set_mate_reference_sequence_id(-1);
            record.set_mate_alignment_start(-1);
        }
        replace_raw_cigar(record, &self.projected_cigar, source_long_cigar)?;
        Ok(ProjectionOutcome {
            kind: if embedded {
                Projection::EmbeddedReference
            } else {
                Projection::Projected
            },
            reference_skip,
        })
    }

    fn set_embedded(
        &mut self,
        header: &sam::Header,
        reference_id: usize,
        name: &str,
    ) -> Result<()> {
        let expected = usize::from(
            header
                .reference_sequences()
                .get_index(reference_id)
                .ok_or_else(|| unknown_reference(reference_id))?
                .1
                .length(),
        );
        if self.padded_query.len() != expected {
            return Err(RsomicsError::InvalidInput(format!(
                "embedded reference {name} has padded length {}, expected {expected}",
                self.padded_query.len()
            )));
        }
        if self.reference.is_some() {
            self.load_reference(header, reference_id)?;
            if self.sequence != self.padded_query {
                return Err(RsomicsError::InvalidInput(format!(
                    "embedded reference {name} does not match the padded FASTA"
                )));
            }
        } else {
            self.sequence.clear();
            self.sequence.extend_from_slice(&self.padded_query);
            self.reference_id = Some(reference_id);
            build_position_map(&self.sequence, &mut self.position_map);
        }
        Ok(())
    }

    fn ensure_reference(&mut self, header: &sam::Header, reference_id: usize) -> Result<()> {
        if self.reference_id == Some(reference_id) {
            return Ok(());
        }
        if self.reference.is_none() {
            let name = reference_name(header, reference_id)?;
            return Err(RsomicsError::InvalidInput(format!(
                "missing {} embedded reference sequence and no padded FASTA was supplied",
                String::from_utf8_lossy(name)
            )));
        }
        self.load_reference(header, reference_id)
    }

    fn load_reference(&mut self, header: &sam::Header, reference_id: usize) -> Result<()> {
        let (name, map) = header
            .reference_sequences()
            .get_index(reference_id)
            .ok_or_else(|| unknown_reference(reference_id))?;
        self.reference
            .as_mut()
            .expect("reference source was checked")
            .read(name.as_ref(), usize::from(map.length()), &mut self.sequence)?;
        self.reference_id = Some(reference_id);
        build_position_map(&self.sequence, &mut self.position_map);
        Ok(())
    }

    fn mapped_position(
        &mut self,
        header: &sam::Header,
        reference_id: usize,
        padded: usize,
        name: &str,
    ) -> Result<usize> {
        if self.reference_id != Some(reference_id) && self.reference.is_none() {
            return Err(RsomicsError::InvalidInput(format!(
                "record {name} has a mate on another reference; --reference is required"
            )));
        }
        self.ensure_reference(header, reference_id)?;
        self.position_map.get(padded).copied().ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "record {name} position {} exceeds the padded reference",
                padded + 1
            ))
        })
    }
}
