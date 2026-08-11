use std::{
    io::{BufWriter, Write},
    sync::Arc,
};

use noodles::fasta::record::Sequence;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use super::{call::BayesianObservation, walker::CalledColumn};

const BASES: &[u8; 17] = b"NACMGRSVTWYHKDBN*";

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub(crate) struct Summary {
    pub(crate) sequence_records: u64,
    pub(crate) sequence_symbols: u64,
    pub(crate) pileup_rows: u64,
}

impl Summary {
    pub(super) fn add(&mut self, other: Self) {
        self.sequence_records += other.sequence_records;
        self.sequence_symbols += other.sequence_symbols;
        self.pileup_rows += other.pileup_rows;
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Format {
    Pileup,
    Fasta,
    Fastq,
}

#[derive(Clone)]
pub(super) struct Reference {
    pub(super) name: Box<[u8]>,
    pub(super) label: Box<[u8]>,
    pub(super) length: u64,
    pub(super) start: i64,
    pub(super) end: i64,
    pub(super) enabled: bool,
}

pub(super) struct Configuration {
    pub(super) format: Format,
    pub(super) show_deletions: bool,
    pub(super) show_insertions: bool,
    pub(super) mark_insertions: bool,
    pub(super) all_positions: u8,
    pub(super) reference_sequences: Option<Vec<Arc<Sequence>>>,
    pub(super) reference_quality: u8,
    pub(super) line_width: usize,
}

pub(super) struct Output<W: Write> {
    writer: BufWriter<W>,
    format: Format,
    show_deletions: bool,
    show_insertions: bool,
    mark_insertions: bool,
    all_positions: u8,
    reference_sequences: Option<Vec<Arc<Sequence>>>,
    reference_quality: u8,
    line_width: usize,
    next_reference_id: i32,
    reference_id: Option<i32>,
    position: Option<i64>,
    sequence: Vec<u8>,
    qualities: Vec<u8>,
    summary: Summary,
}

impl<W: Write> Output<W> {
    pub(super) fn new(writer: W, configuration: Configuration) -> Self {
        Self {
            writer: BufWriter::with_capacity(256 * 1024, writer),
            format: configuration.format,
            show_deletions: configuration.show_deletions,
            show_insertions: configuration.show_insertions,
            mark_insertions: configuration.mark_insertions,
            all_positions: configuration.all_positions,
            reference_sequences: configuration.reference_sequences,
            reference_quality: configuration.reference_quality,
            line_width: configuration.line_width,
            next_reference_id: 0,
            reference_id: None,
            position: None,
            sequence: Vec::new(),
            qualities: Vec::new(),
            summary: Summary::default(),
        }
    }

    pub(super) fn write(
        &mut self,
        references: &[Reference],
        call: CalledColumn,
        observations: &[BayesianObservation],
    ) -> Result<()> {
        let reference = &references[call.reference_id as usize];
        if !reference.enabled || call.position < reference.start || call.position >= reference.end {
            return Ok(());
        }
        self.prepare(references, call)?;
        match self.format {
            Format::Pileup => self.write_pileup(references, call, observations),
            Format::Fasta | Format::Fastq => self.write_sequence(call),
        }
    }

    pub(super) fn begin_selected_reference(&mut self, references: &[Reference]) {
        if self.all_positions == 0 {
            return;
        }
        let reference_id = references
            .iter()
            .position(|reference| reference.enabled)
            .unwrap() as i32;
        self.reference_id = Some(reference_id);
        self.next_reference_id = reference_id + 1;
    }

    pub(super) fn finish(mut self, references: &[Reference]) -> Result<Summary> {
        self.finish_reference(references)?;
        if self.all_positions == 2 {
            for reference_id in self.next_reference_id..references.len() as i32 {
                if !references[reference_id as usize].enabled {
                    continue;
                }
                self.reference_id = Some(reference_id);
                self.finish_reference(references)?;
            }
        }
        self.writer.flush().map_err(RsomicsError::Io)?;
        Ok(self.summary)
    }

    fn prepare(&mut self, references: &[Reference], call: CalledColumn) -> Result<()> {
        if self.reference_id != Some(call.reference_id) {
            self.finish_reference(references)?;
            if self.all_positions == 2 {
                for reference_id in self.next_reference_id..call.reference_id {
                    if !references[reference_id as usize].enabled {
                        continue;
                    }
                    self.reference_id = Some(reference_id);
                    self.finish_reference(references)?;
                }
            }
            self.reference_id = Some(call.reference_id);
            self.next_reference_id = call.reference_id + 1;
        }
        if call.offset == 0 {
            let start = match self.position {
                Some(position) => position + 1,
                None if self.all_positions > 0 => references[call.reference_id as usize].start,
                None => call.position,
            };
            if self.all_positions > 0 || !matches!(self.format, Format::Pileup) {
                for position in start..call.position {
                    self.write_uncovered(references, position)?;
                }
            }
            self.position = Some(call.position);
        }
        Ok(())
    }

    fn finish_reference(&mut self, references: &[Reference]) -> Result<()> {
        if self.all_positions > 0
            && let Some(reference_id) = self.reference_id
        {
            let reference = &references[reference_id as usize];
            if reference.enabled {
                let start = self
                    .position
                    .map_or(reference.start, |position| position + 1);
                for position in start..reference.end {
                    self.write_uncovered(references, position)?;
                }
            }
        }
        if matches!(self.format, Format::Pileup) {
            self.reference_id = None;
            self.position = None;
            Ok(())
        } else {
            self.flush_record(references)
        }
    }

    fn write_uncovered(&mut self, references: &[Reference], position: i64) -> Result<()> {
        let reference_id = self.reference_id.unwrap() as usize;
        let reference_base = self.reference_sequences.as_ref().map(|sequences| {
            let sequence: &[u8] = sequences[reference_id].as_ref().as_ref();
            sequence[position as usize]
        });
        let base = reference_base.unwrap_or(b'N');
        if matches!(self.format, Format::Pileup) {
            self.summary.pileup_rows += 1;
            self.writer
                .write_all(&references[reference_id].name)
                .map_err(RsomicsError::Io)?;
            writeln!(
                self.writer,
                "\t{}\t0\t0\t{}\t0\t*\t*",
                position + 1,
                char::from(base)
            )
            .map_err(RsomicsError::Io)
        } else {
            let quality = reference_base.map_or(0, |_| self.reference_quality);
            self.push(base, i32::from(quality));
            Ok(())
        }
    }

    fn write_pileup(
        &mut self,
        references: &[Reference],
        call: CalledColumn,
        observations: &[BayesianObservation],
    ) -> Result<()> {
        if call.offset > 0 && !self.show_insertions {
            return Ok(());
        }
        if call.base == b'*' && !self.show_deletions {
            return Ok(());
        }
        self.summary.pileup_rows += 1;
        self.writer
            .write_all(&references[call.reference_id as usize].name)
            .map_err(RsomicsError::Io)?;
        write!(
            self.writer,
            "\t{}\t{}\t{}\t{}\t{}\t",
            call.position + 1,
            call.offset,
            call.depth,
            char::from(call.base),
            call.quality
        )
        .map_err(RsomicsError::Io)?;
        for observation in observations {
            let base = if observation.reference_skip {
                b'.'
            } else {
                BASES[usize::from(observation.base.min(16))]
            };
            self.writer.write_all(&[base]).map_err(RsomicsError::Io)?;
        }
        self.writer.write_all(b"\t").map_err(RsomicsError::Io)?;
        for observation in observations {
            let quality = observation.quality.saturating_add(b'!').min(b'~');
            self.writer
                .write_all(&[quality])
                .map_err(RsomicsError::Io)?;
        }
        self.writer.write_all(b"\n").map_err(RsomicsError::Io)
    }

    fn write_sequence(&mut self, call: CalledColumn) -> Result<()> {
        if call.offset > 0 && (!self.show_insertions || call.base == b'*') {
            return Ok(());
        }
        if call.offset == 0 && call.base == b'*' && !self.show_deletions {
            return Ok(());
        }
        if call.offset > 0 && self.mark_insertions {
            self.push(b'_', 62);
        }
        self.push(call.base, call.quality);
        Ok(())
    }

    fn push(&mut self, base: u8, quality: i32) {
        self.sequence.push(base);
        self.summary.sequence_symbols += 1;
        if matches!(self.format, Format::Fastq) {
            self.qualities.push((quality.clamp(0, 93) as u8) + b'!');
        }
    }

    fn flush_record(&mut self, references: &[Reference]) -> Result<()> {
        let Some(reference_id) = self.reference_id.take() else {
            return Ok(());
        };
        if self.sequence.is_empty() {
            self.position = None;
            return Ok(());
        }
        let marker = if matches!(self.format, Format::Fastq) {
            b'@'
        } else {
            b'>'
        };
        self.writer.write_all(&[marker]).map_err(RsomicsError::Io)?;
        self.writer
            .write_all(&references[reference_id as usize].label)
            .map_err(RsomicsError::Io)?;
        self.writer.write_all(b"\n").map_err(RsomicsError::Io)?;
        self.summary.sequence_records += 1;
        write_wrapped(&mut self.writer, &self.sequence, self.line_width)?;
        if matches!(self.format, Format::Fastq) {
            self.writer.write_all(b"+\n").map_err(RsomicsError::Io)?;
            write_wrapped(&mut self.writer, &self.qualities, self.line_width)?;
        }
        self.sequence.clear();
        self.qualities.clear();
        self.position = None;
        Ok(())
    }
}

fn write_wrapped(writer: &mut impl Write, values: &[u8], width: usize) -> Result<()> {
    for line in values.chunks(width) {
        writer.write_all(line).map_err(RsomicsError::Io)?;
        writer.write_all(b"\n").map_err(RsomicsError::Io)?;
    }
    Ok(())
}
