use std::fs::File;
use std::io::{BufWriter, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};

use noodles::core::{Position, Region};
use noodles::fasta;
use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use rsomics_pileup::{
    BaqOptions, Column, FlagFilter, PileupEngine, PileupError, PileupOptions, PileupRead,
    RecordFilter,
};
use serde::Serialize;

use crate::input;

const REVERSE: u16 = 0x10;
const DEFAULT_EXCLUDED_FLAGS: u16 = 0x704;
const UPPER_BASES: &[u8; 16] = b".ACMGRSVTWYHKDBN";
const LOWER_BASES: &[u8; 16] = b",acmgrsvtwyhkdbn";
const INSERTION_BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum BaqMode {
    Disabled,
    #[default]
    Calculate,
    Recalculate,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum PositionMode {
    #[default]
    Covered,
    UsedReferences,
    AllReferences,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Options<'a> {
    pub reference: Option<&'a Path>,
    pub additional_threads: usize,
    pub minimum_base_quality: u8,
    pub minimum_mapping_quality: u8,
    pub maximum_depth: usize,
    pub adjust_overlaps: bool,
    pub include_anomalous_pairs: bool,
    pub excluded_flags: u16,
    pub required_flags: u16,
    pub baq: BaqMode,
    pub positions: PositionMode,
}

impl Default for Options<'_> {
    fn default() -> Self {
        Self {
            reference: None,
            additional_threads: 0,
            minimum_base_quality: 13,
            minimum_mapping_quality: 0,
            maximum_depth: 8000,
            adjust_overlaps: true,
            include_anomalous_pairs: false,
            excluded_flags: DEFAULT_EXCLUDED_FLAGS,
            required_flags: 0,
            baq: BaqMode::Calculate,
            positions: PositionMode::Covered,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub positions: u64,
}

pub fn write(input_path: &Path, options: Options<'_>, output: impl Write) -> Result<Summary> {
    let mut reader = input::open(input_path, options.reference, options.additional_threads)?;
    let header = reader.read_header(input_path)?;
    let references = header
        .reference_sequences()
        .iter()
        .map(|(name, reference)| {
            (
                name.to_vec().into_boxed_slice(),
                usize::from(reference.length()) as u64,
            )
        })
        .collect::<Vec<_>>();
    let pileup_options = PileupOptions {
        filter: RecordFilter {
            flags: FlagFilter {
                skip_any_set: options.excluded_flags,
                skip_any_unset: options.required_flags,
                ..FlagFilter::default()
            },
            minimum_mapping_quality: options.minimum_mapping_quality,
            include_anomalous_pairs: options.include_anomalous_pairs,
        },
        adjust_overlaps: options.adjust_overlaps,
        maximum_depth_per_source: (options.maximum_depth != 0).then_some(options.maximum_depth),
    };
    let mut pileup =
        PileupEngine::new(references.iter().map(|(_, length)| *length), pileup_options);
    let mut reference = options
        .reference
        .map(|path| ReferenceCache::open(path, &references))
        .transpose()?;
    let names = references
        .iter()
        .map(|(name, _)| name.as_ref())
        .collect::<Vec<_>>();
    let lengths = references
        .iter()
        .map(|(_, length)| {
            i64::try_from(*length).map_err(|_| {
                RsomicsError::InvalidInput(
                    "reference length exceeds the supported coordinate range".to_owned(),
                )
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let mut output = BufWriter::with_capacity(256 * 1024, output);
    let mut buffers = TextBuffers::default();
    let mut summary = Summary::default();
    let mut last_reference_id = -1;
    let mut last_position = -1;

    reader.visit_owned_raw_records(&header, input_path, |record| {
        pileup
            .push(record)
            .map_err(|error| pileup_error(input_path, error))?;
        drain(
            &mut pileup,
            &names,
            &lengths,
            &mut reference,
            options,
            &mut output,
            &mut buffers,
            &mut summary,
            &mut last_reference_id,
            &mut last_position,
        )?;
        Ok(true)
    })?;
    pileup
        .finish()
        .map_err(|error| pileup_error(input_path, error))?;
    drain(
        &mut pileup,
        &names,
        &lengths,
        &mut reference,
        options,
        &mut output,
        &mut buffers,
        &mut summary,
        &mut last_reference_id,
        &mut last_position,
    )?;
    emit_trailing(
        &mut output,
        &mut buffers.line,
        &names,
        &lengths,
        &mut reference,
        options.positions,
        &mut last_reference_id,
        &mut last_position,
        &mut summary,
    )?;
    output.flush().map_err(RsomicsError::Io)?;
    Ok(summary)
}

#[allow(clippy::too_many_arguments)]
fn drain(
    pileup: &mut PileupEngine,
    names: &[&[u8]],
    lengths: &[i64],
    reference: &mut Option<ReferenceCache>,
    options: Options<'_>,
    output: &mut impl Write,
    buffers: &mut TextBuffers,
    summary: &mut Summary,
    last_reference_id: &mut i32,
    last_position: &mut i64,
) -> Result<()> {
    pileup
        .drain_with(|context| {
            if let Some(reference) = reference.as_mut() {
                let baq = BaqOptions {
                    adjust_qualities: true,
                    extended: true,
                    redo: options.baq == BaqMode::Recalculate,
                };
                if options.baq != BaqMode::Disabled {
                    context.apply_full_baq(
                        usize::MAX,
                        baq,
                        |reference_id, range, buffer: &mut Vec<u8>| {
                            buffer.extend_from_slice(reference.sequence(reference_id, range)?);
                            Ok::<_, StreamError>(())
                        },
                    )?;
                }
            }
            let column = context.column();
            emit_gaps(
                output,
                &mut buffers.line,
                names,
                lengths,
                reference,
                options.positions,
                column.reference_id(),
                column.position(),
                last_reference_id,
                last_position,
                summary,
            )?;
            encode_column(
                buffers,
                names,
                reference,
                options.minimum_base_quality,
                &column,
            )?;
            output.write_all(&buffers.line).map_err(RsomicsError::Io)?;
            summary.positions = summary
                .positions
                .checked_add(1)
                .ok_or_else(position_overflow)?;
            *last_reference_id = column.reference_id();
            *last_position = column.position();
            Ok::<_, StreamError>(())
        })
        .map_err(StreamError::into_rsomics)
}

fn encode_column(
    buffers: &mut TextBuffers,
    names: &[&[u8]],
    reference: &mut Option<ReferenceCache>,
    minimum_base_quality: u8,
    column: &Column<'_>,
) -> Result<()> {
    buffers.line.clear();
    buffers.bases.clear();
    buffers.qualities.clear();
    let reference_id = usize::try_from(column.reference_id()).unwrap();
    let name = names.get(reference_id).copied().unwrap_or(b"*");
    let reference_base = reference_base(reference, column.reference_id(), column.position())?;
    buffers.line.extend_from_slice(name);
    buffers.line.push(b'\t');
    push_i64(&mut buffers.line, column.position() + 1);
    buffers.line.push(b'\t');
    buffers.line.push(reference_base);

    buffers.bases.reserve(column.len().saturating_mul(2));
    buffers.qualities.reserve(column.len());
    for entry in column.entries() {
        let quality = base_quality(entry.record(), entry.projection());
        if quality < minimum_base_quality {
            continue;
        }
        encode_base(
            &mut buffers.bases,
            entry.record(),
            entry.projection(),
            reference_base,
            reference,
            column.reference_id(),
            column.position(),
        )?;
        buffers.qualities.push(quality.saturating_add(33).min(126));
    }

    buffers.line.push(b'\t');
    push_usize(&mut buffers.line, buffers.qualities.len());
    buffers.line.push(b'\t');
    if buffers.bases.is_empty() {
        buffers.line.push(b'*');
    } else {
        buffers.line.extend_from_slice(&buffers.bases);
    }
    buffers.line.push(b'\t');
    if buffers.qualities.is_empty() {
        buffers.line.push(b'*');
    } else {
        buffers.line.extend_from_slice(&buffers.qualities);
    }
    buffers.line.push(b'\n');
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn encode_base(
    output: &mut Vec<u8>,
    record: &RawRecord,
    read: &PileupRead,
    reference_nucleotide: u8,
    reference: &mut Option<ReferenceCache>,
    reference_id: i32,
    position: i64,
) -> Result<()> {
    let reverse = record.flags() & REVERSE != 0;
    if read.is_head {
        output.push(b'^');
        output.push(record.mapping_quality().saturating_add(33).min(126));
    }

    if read.is_deletion {
        output.push(if read.is_reference_skip {
            if reverse { b'<' } else { b'>' }
        } else {
            b'*'
        });
    } else {
        let mut base = if read.qpos < record.sequence_len() {
            record.seq_nibble(read.qpos)
        } else {
            15
        };
        if reference.is_some() && base == reference_code(reference_nucleotide) {
            base = 0;
        }
        output.push(if reverse {
            LOWER_BASES[usize::from(base)]
        } else {
            UPPER_BASES[usize::from(base)]
        });
    }

    if read.indel > 0 {
        output.push(b'+');
        push_i64(output, read.indel);
        let deletion_offset = usize::from(read.is_deletion);
        for offset in 1..=usize::try_from(read.indel).unwrap() {
            let query_position = read.qpos + offset - deletion_offset;
            let base = if query_position < record.sequence_len() {
                INSERTION_BASES[usize::from(record.seq_nibble(query_position))]
            } else {
                b'N'
            };
            output.push(if reverse {
                base.to_ascii_lowercase()
            } else {
                base
            });
        }
    } else if read.indel < 0 {
        output.push(b'-');
        let length = read.indel.unsigned_abs();
        push_u64(output, length);
        for offset in 1..=length {
            let offset = i64::try_from(offset).unwrap();
            let base = reference_base(reference, reference_id, position + offset)?;
            output.push(if reverse {
                base.to_ascii_lowercase()
            } else {
                base.to_ascii_uppercase()
            });
        }
    }

    if read.is_tail {
        output.push(b'$');
    }
    Ok(())
}

fn base_quality(record: &RawRecord, read: &PileupRead) -> u8 {
    record.quality_scores().get(read.qpos).copied().unwrap_or(0)
}

#[allow(clippy::too_many_arguments)]
fn emit_gaps(
    output: &mut impl Write,
    line: &mut Vec<u8>,
    names: &[&[u8]],
    lengths: &[i64],
    reference: &mut Option<ReferenceCache>,
    mode: PositionMode,
    reference_id: i32,
    position: i64,
    last_reference_id: &mut i32,
    last_position: &mut i64,
    summary: &mut Summary,
) -> Result<()> {
    if mode == PositionMode::Covered {
        return Ok(());
    }
    while reference_id > *last_reference_id {
        if *last_reference_id >= 0 {
            let length = lengths
                .get(usize::try_from(*last_reference_id).unwrap())
                .copied()
                .unwrap_or(0);
            while *last_position + 1 < length {
                *last_position += 1;
                emit_empty(
                    output,
                    line,
                    names,
                    reference,
                    *last_reference_id,
                    *last_position,
                    summary,
                )?;
            }
        }
        *last_reference_id += 1;
        *last_position = -1;
        if mode == PositionMode::UsedReferences {
            break;
        }
    }
    if *last_reference_id != reference_id {
        *last_reference_id = reference_id;
        *last_position = -1;
    }
    while *last_position + 1 < position {
        *last_position += 1;
        emit_empty(
            output,
            line,
            names,
            reference,
            reference_id,
            *last_position,
            summary,
        )?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn emit_trailing(
    output: &mut impl Write,
    line: &mut Vec<u8>,
    names: &[&[u8]],
    lengths: &[i64],
    reference: &mut Option<ReferenceCache>,
    mode: PositionMode,
    last_reference_id: &mut i32,
    last_position: &mut i64,
    summary: &mut Summary,
) -> Result<()> {
    if mode == PositionMode::Covered {
        return Ok(());
    }
    if *last_reference_id < 0 && mode == PositionMode::AllReferences {
        *last_reference_id = 0;
    }
    let reference_count = i32::try_from(names.len()).unwrap_or(i32::MAX);
    while *last_reference_id >= 0 && *last_reference_id < reference_count {
        let length = lengths
            .get(usize::try_from(*last_reference_id).unwrap())
            .copied()
            .unwrap_or(0);
        while *last_position + 1 < length {
            *last_position += 1;
            emit_empty(
                output,
                line,
                names,
                reference,
                *last_reference_id,
                *last_position,
                summary,
            )?;
        }
        *last_reference_id += 1;
        *last_position = -1;
        if mode == PositionMode::UsedReferences {
            break;
        }
    }
    Ok(())
}

fn emit_empty(
    output: &mut impl Write,
    line: &mut Vec<u8>,
    names: &[&[u8]],
    reference: &mut Option<ReferenceCache>,
    reference_id: i32,
    position: i64,
    summary: &mut Summary,
) -> Result<()> {
    line.clear();
    let name = names
        .get(usize::try_from(reference_id).unwrap())
        .copied()
        .unwrap_or(b"*");
    line.extend_from_slice(name);
    line.push(b'\t');
    push_i64(line, position + 1);
    line.push(b'\t');
    line.push(reference_base(reference, reference_id, position)?);
    line.extend_from_slice(b"\t0\t*\t*\n");
    output.write_all(line).map_err(RsomicsError::Io)?;
    summary.positions = summary
        .positions
        .checked_add(1)
        .ok_or_else(position_overflow)?;
    Ok(())
}

fn reference_base(
    reference: &mut Option<ReferenceCache>,
    reference_id: i32,
    position: i64,
) -> Result<u8> {
    let Some(reference) = reference else {
        return Ok(b'N');
    };
    let position = usize::try_from(position).map_err(|error| reference.error(error))?;
    let end = position
        .checked_add(1)
        .ok_or_else(|| reference.error("reference position overflows"))?;
    Ok(reference.sequence(reference_id, position..end)?[0])
}

fn reference_code(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'=' => 0,
        b'A' => 1,
        b'C' => 2,
        b'M' => 3,
        b'G' => 4,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'T' => 8,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        _ => 15,
    }
}

fn push_i64(output: &mut Vec<u8>, value: i64) {
    let mut buffer = itoa::Buffer::new();
    output.extend_from_slice(buffer.format(value).as_bytes());
}

fn push_u64(output: &mut Vec<u8>, value: u64) {
    let mut buffer = itoa::Buffer::new();
    output.extend_from_slice(buffer.format(value).as_bytes());
}

fn push_usize(output: &mut Vec<u8>, value: usize) {
    let mut buffer = itoa::Buffer::new();
    output.extend_from_slice(buffer.format(value).as_bytes());
}

struct TextBuffers {
    line: Vec<u8>,
    bases: Vec<u8>,
    qualities: Vec<u8>,
}

impl Default for TextBuffers {
    fn default() -> Self {
        Self {
            line: Vec::with_capacity(256),
            bases: Vec::with_capacity(256),
            qualities: Vec::with_capacity(128),
        }
    }
}

struct ReferenceCache {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    path: PathBuf,
    references: Vec<(Box<[u8]>, u64)>,
    reference_id: Option<usize>,
    sequence_start: usize,
    sequence: Vec<u8>,
}

impl ReferenceCache {
    fn open(path: &Path, references: &[(Box<[u8]>, u64)]) -> Result<Self> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .map_err(|error| reference_error(path, error))?;
        Ok(Self {
            reader,
            path: path.to_path_buf(),
            references: references.to_vec(),
            reference_id: None,
            sequence_start: 0,
            sequence: Vec::new(),
        })
    }

    fn sequence(&mut self, reference_id: i32, range: Range<usize>) -> Result<&[u8]> {
        const CHUNK_SIZE: usize = 1024 * 1024;

        let reference_id = usize::try_from(reference_id).map_err(|error| self.error(error))?;
        let (name, length) = self
            .references
            .get(reference_id)
            .map(|(name, length)| (name.to_vec(), *length))
            .ok_or_else(|| self.error("reference ID is absent"))?;
        let length = usize::try_from(length).map_err(|error| self.error(error))?;
        if range.start >= range.end || range.end > length {
            return Err(self.error("requested reference range is invalid"));
        }
        if self.reference_id != Some(reference_id)
            || range.start < self.sequence_start
            || range.end > self.sequence_start + self.sequence.len()
        {
            let start = range.start / CHUNK_SIZE * CHUNK_SIZE;
            let chunk_end = start
                .checked_add(CHUNK_SIZE)
                .map_or(length, |end| end.min(length));
            let end = chunk_end.max(range.end);
            let interval_start =
                Position::try_from(start + 1).map_err(|error| self.error(error))?;
            let interval_end = Position::try_from(end).map_err(|error| self.error(error))?;
            let record = self
                .reader
                .query(&Region::new(name, interval_start..=interval_end))
                .map_err(|error| self.error(error))?;
            self.sequence.clear();
            self.sequence.extend_from_slice(record.sequence().as_ref());
            self.reference_id = Some(reference_id);
            self.sequence_start = start;
        }
        let start = range
            .start
            .checked_sub(self.sequence_start)
            .ok_or_else(|| self.error("invalid reference cache position"))?;
        let end = range
            .end
            .checked_sub(self.sequence_start)
            .ok_or_else(|| self.error("invalid reference cache position"))?;
        self.sequence
            .get(start..end)
            .ok_or_else(|| self.error("reference range is outside the cache"))
    }

    fn error(&self, error: impl std::fmt::Display) -> RsomicsError {
        reference_error(&self.path, error)
    }
}

enum StreamError {
    Pileup(PileupError),
    Product(RsomicsError),
}

impl From<PileupError> for StreamError {
    fn from(error: PileupError) -> Self {
        Self::Pileup(error)
    }
}

impl From<RsomicsError> for StreamError {
    fn from(error: RsomicsError) -> Self {
        Self::Product(error)
    }
}

impl StreamError {
    fn into_rsomics(self) -> RsomicsError {
        match self {
            Self::Pileup(error) => RsomicsError::InvalidInput(format!("building pileup: {error}")),
            Self::Product(error) => error,
        }
    }
}

fn pileup_error(path: &Path, error: PileupError) -> RsomicsError {
    RsomicsError::InvalidInput(format!("building pileup from {}: {error}", path.display()))
}

fn reference_error(path: &Path, error: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!("reading reference {}: {error}", path.display()))
}

fn position_overflow() -> RsomicsError {
    RsomicsError::InvalidInput("pileup position count exceeds u64".to_owned())
}
