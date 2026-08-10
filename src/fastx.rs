use std::io::Write;
use std::path::Path;

use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::{Compression, Format as SequenceFormat, OutputWriter, Record};
use serde::Serialize;

use crate::input;

const REVERSE: u16 = 0x10;
const READ1: u16 = 0x40;
const READ2: u16 = 0x80;
const BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
const COMPLEMENT: [u8; 16] = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    Fasta,
    Fastq,
}

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub format: Format,
    pub compression: Compression,
    pub append_mate_suffix: bool,
    pub use_original_quality: bool,
    pub default_quality: u8,
    pub require_flags: u16,
    pub exclude_flags: u16,
    pub include_flags: u16,
    pub exclude_all_flags: u16,
    pub reference: Option<&'a Path>,
    pub additional_threads: usize,
}

impl Default for Options<'_> {
    fn default() -> Self {
        Self {
            format: Format::Fastq,
            compression: Compression::Plain,
            append_mate_suffix: true,
            use_original_quality: false,
            default_quality: 1,
            require_flags: 0,
            exclude_flags: 0x900,
            include_flags: 0,
            exclude_all_flags: 0,
            reference: None,
            additional_threads: 0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub records_read: u64,
    pub records_filtered: u64,
    pub records_written: u64,
}

pub fn write<W: Write>(input_path: &Path, options: Options<'_>, output: W) -> Result<Summary> {
    if options.default_quality > 93 {
        return Err(RsomicsError::ConfigError(format!(
            "default quality must be between 0 and 93, got {}",
            options.default_quality
        )));
    }

    let format = match options.format {
        Format::Fasta => SequenceFormat::Fasta,
        Format::Fastq => SequenceFormat::Fastq,
    };
    let writer = OutputWriter::new(output, format, options.compression)?;
    let mut emitter = Emitter::new(writer, options);
    let mut reader = input::open(input_path, options.reference, options.additional_threads)?;
    let header = reader.read_header(input_path)?;
    reader.visit_owned_raw_records(&header, input_path, |record| {
        emitter.push(record)?;
        Ok(true)
    })?;
    emitter.finish()
}

struct Candidate {
    record: RawRecord,
    has_quality: bool,
}

struct Emitter<'a, W: Write> {
    writer: OutputWriter<W>,
    options: Options<'a>,
    name: Vec<u8>,
    candidates: [Option<Candidate>; 3],
    output_name: Vec<u8>,
    sequence: Vec<u8>,
    quality: Vec<u8>,
    summary: Summary,
}

impl<'a, W: Write> Emitter<'a, W> {
    fn new(writer: OutputWriter<W>, options: Options<'a>) -> Self {
        Self {
            writer,
            options,
            name: Vec::new(),
            candidates: std::array::from_fn(|_| None),
            output_name: Vec::new(),
            sequence: Vec::new(),
            quality: Vec::new(),
            summary: Summary::default(),
        }
    }

    fn push(&mut self, record: RawRecord) -> Result<()> {
        self.summary.records_read = checked_increment(self.summary.records_read)?;
        if !accepts(record.flags(), self.options) {
            self.summary.records_filtered = checked_increment(self.summary.records_filtered)?;
            return Ok(());
        }

        if !self.name.is_empty() && self.name != record.name() {
            self.flush()?;
        }
        if self.name.is_empty() {
            self.name.extend_from_slice(record.name());
        }

        let category = category(record.flags());
        let has_quality = quality_source(&record, self.options.use_original_quality)?.is_some();
        let replace = self.candidates[category]
            .as_ref()
            .is_none_or(|candidate| !candidate.has_quality && has_quality);
        if replace {
            self.candidates[category] = Some(Candidate {
                record,
                has_quality,
            });
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        for category in [1, 2, 0] {
            if let Some(candidate) = self.candidates[category].take() {
                self.write_candidate(&candidate.record, category)?;
            }
        }
        self.name.clear();
        Ok(())
    }

    fn write_candidate(&mut self, record: &RawRecord, category: usize) -> Result<()> {
        self.output_name.clear();
        self.output_name.extend_from_slice(record.name());
        if self.options.append_mate_suffix {
            match category {
                1 => self.output_name.extend_from_slice(b"/1"),
                2 => self.output_name.extend_from_slice(b"/2"),
                _ => {}
            }
        }

        decode_sequence(record, &mut self.sequence);
        let quality = if self.options.format == Format::Fastq {
            encode_quality(record, self.options, &mut self.quality)?;
            Some(self.quality.as_slice())
        } else {
            None
        };
        self.writer.write_record(Record {
            id: &self.output_name,
            seq: &self.sequence,
            qual: quality,
        })?;
        self.summary.records_written = checked_increment(self.summary.records_written)?;
        Ok(())
    }

    fn finish(mut self) -> Result<Summary> {
        self.flush()?;
        self.writer.finish()?;
        Ok(self.summary)
    }
}

fn accepts(flags: u16, options: Options<'_>) -> bool {
    flags & options.require_flags == options.require_flags
        && flags & options.exclude_flags == 0
        && (options.include_flags == 0 || flags & options.include_flags != 0)
        && (options.exclude_all_flags == 0
            || flags & options.exclude_all_flags != options.exclude_all_flags)
}

fn category(flags: u16) -> usize {
    match (flags & READ1 != 0, flags & READ2 != 0) {
        (true, false) => 1,
        (false, true) => 2,
        _ => 0,
    }
}

enum Quality<'a> {
    Phred(&'a [u8]),
    Ascii(&'a [u8]),
}

fn quality_source(record: &RawRecord, original: bool) -> Result<Option<Quality<'_>>> {
    if original && record.aux_type(*b"OQ") == Some(b'Z') {
        let quality = record
            .aux_value(*b"OQ")
            .and_then(|value| value.strip_suffix(&[0]))
            .ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "record {} has a malformed OQ tag",
                    String::from_utf8_lossy(record.name())
                ))
            })?;
        if quality.len() != record.sequence_len() {
            return Err(RsomicsError::InvalidInput(format!(
                "record {} has OQ length {}, expected {}",
                String::from_utf8_lossy(record.name()),
                quality.len(),
                record.sequence_len()
            )));
        }
        return Ok(Some(Quality::Ascii(quality)));
    }

    let quality = record.quality_scores();
    if quality.is_empty() && record.sequence_len() > 0 {
        Ok(None)
    } else {
        Ok(Some(Quality::Phred(quality)))
    }
}

fn decode_sequence(record: &RawRecord, output: &mut Vec<u8>) {
    output.clear();
    output.reserve(record.sequence_len());
    if record.flags() & REVERSE == 0 {
        output.extend((0..record.sequence_len()).map(|i| BASES[usize::from(record.seq_nibble(i))]));
    } else {
        output.extend((0..record.sequence_len()).rev().map(|i| {
            let code = COMPLEMENT[usize::from(record.seq_nibble(i))];
            BASES[usize::from(code)]
        }));
    }
}

fn encode_quality(record: &RawRecord, options: Options<'_>, output: &mut Vec<u8>) -> Result<()> {
    output.clear();
    output.reserve(record.sequence_len());
    match quality_source(record, options.use_original_quality)? {
        Some(Quality::Phred(scores)) => {
            for &score in scores {
                if score > 93 {
                    return Err(RsomicsError::InvalidInput(format!(
                        "record {} has quality score {score} above 93",
                        String::from_utf8_lossy(record.name())
                    )));
                }
                output.push(score + b'!');
            }
        }
        Some(Quality::Ascii(scores)) => {
            if scores.iter().any(|&score| !(b'!'..=b'~').contains(&score)) {
                return Err(RsomicsError::InvalidInput(format!(
                    "record {} has an OQ byte outside ASCII 33..=126",
                    String::from_utf8_lossy(record.name())
                )));
            }
            output.extend_from_slice(scores);
        }
        None => output.resize(record.sequence_len(), options.default_quality + b'!'),
    }
    if record.flags() & REVERSE != 0 {
        output.reverse();
    }
    Ok(())
}

fn checked_increment(value: u64) -> Result<u64> {
    value
        .checked_add(1)
        .ok_or_else(|| RsomicsError::InvalidInput("alignment record count exceeds u64".to_owned()))
}
