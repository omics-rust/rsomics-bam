use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};

use noodles::core::Region;
use noodles::sam;
use rsomics_bamio::raw::{RawRecord, RawRecordEncoder, RecordRef};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::input;

const FLAG_PAIRED: u16 = 0x1;
const FLAG_MATE_UNMAPPED: u16 = 0x8;
const MATCH: u8 = 0;
const INSERTION: u8 = 1;
const DELETION: u8 = 2;
const REFERENCE_SKIP: u8 = 3;
const SOFT_CLIP: u8 = 4;
const HARD_CLIP: u8 = 5;
const PADDING: u8 = 6;
const SEQUENCE_MATCH: u8 = 7;
const SEQUENCE_MISMATCH: u8 = 8;

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
    pub minimum_read_length: usize,
    pub excluded_flags: u16,
    pub included_flags: u16,
    pub required_flags: u16,
    pub include_deletions: bool,
    pub remove_overlaps: bool,
    pub positions: PositionMode,
    pub region: Option<&'a str>,
    pub bed: Option<&'a Path>,
    pub header: bool,
}

impl Default for Options<'_> {
    fn default() -> Self {
        Self {
            reference: None,
            additional_threads: 0,
            minimum_base_quality: 0,
            minimum_mapping_quality: 0,
            minimum_read_length: 0,
            excluded_flags: 0x704,
            included_flags: 0,
            required_flags: 0,
            include_deletions: false,
            remove_overlaps: false,
            positions: PositionMode::Covered,
            region: None,
            bed: None,
            header: false,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub inputs: usize,
    pub positions: u64,
}

enum Message {
    Header(Box<sam::Header>),
    Record(RawRecord),
    Finished,
    Error(RsomicsError),
}

struct Stream {
    input: PathBuf,
    receiver: Receiver<Message>,
    next: Option<RawRecord>,
    finished: bool,
    last_coordinate: Option<(i32, i32)>,
}

pub fn write(inputs: &[PathBuf], options: Options<'_>, output: impl Write) -> Result<Summary> {
    if inputs.is_empty() {
        return Err(RsomicsError::ConfigError(
            "depth requires at least one alignment input".to_owned(),
        ));
    }
    if inputs.len() > 1 && inputs.iter().any(|input| input == Path::new("-")) {
        return Err(RsomicsError::ConfigError(
            "standard input cannot be combined with other depth inputs".to_owned(),
        ));
    }
    let region = options
        .region
        .map(str::parse::<Region>)
        .transpose()
        .map_err(|error| RsomicsError::ConfigError(format!("invalid region: {error}")))?;

    if let [input] = inputs {
        return write_single(input, options, region.as_ref(), output);
    }
    write_multiple(inputs, options, region.as_ref(), output)
}

fn write_single(
    path: &Path,
    options: Options<'_>,
    region: Option<&Region>,
    output: impl Write,
) -> Result<Summary> {
    let mut reader = if region.is_some() {
        input::open_indexed(path, options.reference)?
    } else {
        input::open(path, options.reference, options.additional_threads)?
    };
    let header = reader.read_header(path)?;
    let references = reference_dictionary(&header);
    let selection = Selection::resolve(region, &header, &references)?;
    let bed = options
        .bed
        .map(|bed| Bed::read(bed, &references))
        .transpose()?;
    let mut accumulator = Accumulator::new(references, selection, bed, 1, options, output);
    accumulator.write_header(&[path.to_string_lossy().into_owned()])?;

    if let Some(region) = region {
        let mut encoder = RawRecordEncoder::new();
        reader.visit_region(&header, path, Some(region), |record| {
            let record = encoder.encode(&header, record)?;
            accumulator.add(0, &record)?;
            Ok(true)
        })?;
    } else if reader.format() == input::Format::Bam {
        reader.visit_raw_bam_records(path, |record| {
            accumulator.add(0, &record)?;
            Ok(true)
        })?;
    } else {
        reader.visit_owned_raw_records(&header, path, |record| {
            accumulator.add(0, &record)?;
            Ok(true)
        })?;
    }
    accumulator.finish()
}

fn write_multiple(
    inputs: &[PathBuf],
    options: Options<'_>,
    region: Option<&Region>,
    output: impl Write,
) -> Result<Summary> {
    std::thread::scope(|scope| {
        let mut streams = Vec::with_capacity(inputs.len());
        for input in inputs {
            let (sender, receiver) = sync_channel(1);
            let input = input.clone();
            let worker_input = input.clone();
            let worker_region = region.cloned();
            scope.spawn(move || {
                if let Err(error) = read_input(
                    &worker_input,
                    options.reference,
                    options.additional_threads,
                    worker_region.as_ref(),
                    &sender,
                ) {
                    let _ = sender.send(Message::Error(error));
                }
            });
            streams.push(Stream {
                input,
                receiver,
                next: None,
                finished: false,
                last_coordinate: None,
            });
        }

        let mut headers = Vec::with_capacity(streams.len());
        for stream in &mut streams {
            headers.push(receive_header(stream)?);
        }
        let references = reference_dictionary(&headers[0]);
        for (index, header) in headers.iter().enumerate().skip(1) {
            if reference_dictionary(header) != references {
                return Err(RsomicsError::InvalidInput(format!(
                    "reference dictionary in {} differs from {}",
                    streams[index].input.display(),
                    streams[0].input.display()
                )));
            }
        }

        let selection = Selection::resolve(region, &headers[0], &references)?;
        let bed = options
            .bed
            .map(|path| Bed::read(path, &references))
            .transpose()?;
        let names = inputs
            .iter()
            .map(|input| input.to_string_lossy().into_owned())
            .collect::<Vec<_>>();
        let mut accumulator =
            Accumulator::new(references, selection, bed, inputs.len(), options, output);
        accumulator.write_header(&names)?;

        for stream in &mut streams {
            receive_next(stream)?;
        }
        while let Some(index) = streams
            .iter()
            .enumerate()
            .filter_map(|(index, stream)| {
                stream.next.as_ref().map(|record| {
                    (
                        index,
                        record.reference_sequence_id(),
                        record.alignment_start(),
                    )
                })
            })
            .min_by_key(|&(index, reference_id, position)| (reference_id, position, index))
            .map(|(index, _, _)| index)
        {
            let record = streams[index].next.take().unwrap();
            let coordinate = (record.reference_sequence_id(), record.alignment_start());
            if coordinate.0 >= 0 {
                if streams[index]
                    .last_coordinate
                    .is_some_and(|previous| coordinate < previous)
                {
                    return Err(RsomicsError::InvalidInput(format!(
                        "alignment input is not coordinate sorted: {}",
                        streams[index].input.display()
                    )));
                }
                streams[index].last_coordinate = Some(coordinate);
            }
            accumulator.add(index, &record)?;
            receive_next(&mut streams[index])?;
        }
        accumulator.finish()
    })
}

fn read_input(
    path: &Path,
    reference: Option<&Path>,
    additional_threads: usize,
    region: Option<&Region>,
    sender: &SyncSender<Message>,
) -> Result<()> {
    let mut reader = if region.is_some() {
        input::open_indexed(path, reference)?
    } else {
        input::open(path, reference, additional_threads)?
    };
    let header = reader.read_header(path)?;
    if sender
        .send(Message::Header(Box::new(header.clone())))
        .is_err()
    {
        return Ok(());
    }

    if let Some(region) = region {
        let mut encoder = RawRecordEncoder::new();
        reader.visit_region(&header, path, Some(region), |record| {
            let record = encoder.encode(&header, record)?;
            Ok(sender.send(Message::Record(record)).is_ok())
        })?;
    } else {
        reader.visit_owned_raw_records(&header, path, |record| {
            Ok(sender.send(Message::Record(record)).is_ok())
        })?;
    }
    let _ = sender.send(Message::Finished);
    Ok(())
}

fn receive_header(stream: &mut Stream) -> Result<sam::Header> {
    match stream.receiver.recv() {
        Ok(Message::Header(header)) => Ok(*header),
        Ok(Message::Error(error)) => Err(error),
        Ok(Message::Record(_) | Message::Finished) => Err(RsomicsError::InvalidInput(format!(
            "alignment stream ended before its header: {}",
            stream.input.display()
        ))),
        Err(_) => Err(RsomicsError::InvalidInput(format!(
            "alignment reader stopped unexpectedly: {}",
            stream.input.display()
        ))),
    }
}

fn receive_next(stream: &mut Stream) -> Result<()> {
    if stream.finished {
        stream.next = None;
        return Ok(());
    }
    match stream.receiver.recv() {
        Ok(Message::Record(record)) => stream.next = Some(record),
        Ok(Message::Finished) => {
            stream.next = None;
            stream.finished = true;
        }
        Ok(Message::Error(error)) => return Err(error),
        Ok(Message::Header(_)) => {
            return Err(RsomicsError::InvalidInput(format!(
                "alignment stream emitted more than one header: {}",
                stream.input.display()
            )));
        }
        Err(_) => {
            return Err(RsomicsError::InvalidInput(format!(
                "alignment reader stopped unexpectedly: {}",
                stream.input.display()
            )));
        }
    }
    Ok(())
}

fn reference_dictionary(header: &sam::Header) -> Vec<(Vec<u8>, u64)> {
    header
        .reference_sequences()
        .iter()
        .map(|(name, reference)| {
            (
                name.to_vec(),
                u64::try_from(usize::from(reference.length())).unwrap(),
            )
        })
        .collect()
}

#[derive(Clone, Copy)]
struct Selection {
    reference_id: Option<usize>,
    start: u64,
    end: u64,
}

impl Selection {
    fn resolve(
        region: Option<&Region>,
        header: &sam::Header,
        references: &[(Vec<u8>, u64)],
    ) -> Result<Self> {
        let Some(region) = region else {
            return Ok(Self {
                reference_id: None,
                start: 0,
                end: 0,
            });
        };
        let reference_id = header
            .reference_sequences()
            .get_index_of(region.name())
            .ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "region reference is absent from the alignment header: {}",
                    String::from_utf8_lossy(region.name())
                ))
            })?;
        let reference_length = references[reference_id].1;
        let start = region
            .interval()
            .start()
            .map(|position| u64::try_from(usize::from(position) - 1).unwrap())
            .unwrap_or(0);
        let end = region
            .interval()
            .end()
            .map(|position| u64::try_from(usize::from(position)).unwrap())
            .unwrap_or(reference_length)
            .min(reference_length);
        if start >= end {
            return Err(RsomicsError::InvalidInput(format!(
                "region starts outside reference length {reference_length}"
            )));
        }
        Ok(Self {
            reference_id: Some(reference_id),
            start,
            end,
        })
    }

    fn bounds(self, reference_id: usize, reference_length: u64) -> Option<Range<u64>> {
        match self.reference_id {
            Some(selected) if selected == reference_id => Some(self.start..self.end),
            Some(_) => None,
            None => Some(0..reference_length),
        }
    }
}

struct Bed {
    intervals: HashMap<usize, Vec<Range<u64>>>,
}

impl Bed {
    fn read(path: &Path, references: &[(Vec<u8>, u64)]) -> Result<Self> {
        let names = references
            .iter()
            .enumerate()
            .map(|(index, (name, _))| (name.as_slice(), index))
            .collect::<HashMap<_, _>>();
        let file = std::fs::File::open(path).map_err(RsomicsError::Io)?;
        let mut intervals: HashMap<usize, Vec<Range<u64>>> = HashMap::new();
        for (index, result) in BufReader::new(file).lines().enumerate() {
            let line = result.map_err(RsomicsError::Io)?;
            if line.is_empty() || line.starts_with('#') || line.starts_with("track ") {
                continue;
            }
            let mut fields = line.split('\t');
            let name = fields.next().unwrap();
            let start = fields
                .next()
                .ok_or_else(|| invalid_bed(path, index + 1, "missing start"))?
                .parse::<u64>()
                .map_err(|_| invalid_bed(path, index + 1, "invalid start"))?;
            let end = fields
                .next()
                .ok_or_else(|| invalid_bed(path, index + 1, "missing end"))?
                .parse::<u64>()
                .map_err(|_| invalid_bed(path, index + 1, "invalid end"))?;
            if start >= end {
                return Err(invalid_bed(
                    path,
                    index + 1,
                    "start must be smaller than end",
                ));
            }
            let Some(&reference_id) = names.get(name.as_bytes()) else {
                continue;
            };
            let reference_length = references[reference_id].1;
            if start >= reference_length {
                continue;
            }
            intervals
                .entry(reference_id)
                .or_default()
                .push(start..end.min(reference_length));
        }
        for ranges in intervals.values_mut() {
            ranges.sort_unstable_by_key(|range| (range.start, range.end));
            let mut merged: Vec<Range<u64>> = Vec::with_capacity(ranges.len());
            for range in ranges.drain(..) {
                if let Some(previous) = merged.last_mut()
                    && range.start <= previous.end
                {
                    previous.end = previous.end.max(range.end);
                } else {
                    merged.push(range);
                }
            }
            *ranges = merged;
        }
        Ok(Self { intervals })
    }

    fn contains(&self, reference_id: usize, position: u64) -> bool {
        let Some(ranges) = self.intervals.get(&reference_id) else {
            return false;
        };
        let index = ranges.partition_point(|range| range.end <= position);
        ranges
            .get(index)
            .is_some_and(|range| range.start <= position)
    }
}

fn invalid_bed(path: &Path, line: usize, reason: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{}:{line}: {reason}", path.display()))
}

trait DepthRecord {
    fn flags(&self) -> u16;
    fn reference_sequence_id(&self) -> i32;
    fn alignment_start(&self) -> i32;
    fn mapping_quality(&self) -> u8;
    fn mate_reference_sequence_id(&self) -> i32;
    fn mate_alignment_start(&self) -> i32;
    fn name(&self) -> &[u8];
    fn copy_cigar(&self, buffer: &mut Vec<(u8, u32)>) -> Result<()>;
    fn sequence_len(&self) -> usize;
    fn quality_scores(&self) -> &[u8];
}

macro_rules! impl_depth_record {
    ($ty:ty) => {
        impl DepthRecord for $ty {
            fn flags(&self) -> u16 {
                self.flags()
            }

            fn reference_sequence_id(&self) -> i32 {
                self.reference_sequence_id()
            }

            fn alignment_start(&self) -> i32 {
                self.alignment_start()
            }

            fn mapping_quality(&self) -> u8 {
                self.mapping_quality()
            }

            fn mate_reference_sequence_id(&self) -> i32 {
                self.mate_reference_sequence_id()
            }

            fn mate_alignment_start(&self) -> i32 {
                self.mate_alignment_start()
            }

            fn name(&self) -> &[u8] {
                self.name()
            }

            fn copy_cigar(&self, buffer: &mut Vec<(u8, u32)>) -> Result<()> {
                buffer.clear();
                buffer.extend(self.cigar_ops());
                if buffer.len() == 2
                    && buffer[0] == (SOFT_CLIP, u32::try_from(self.sequence_len()).unwrap())
                    && buffer[1].0 == REFERENCE_SKIP
                    && self.aux_type(*b"CG") == Some(b'B')
                {
                    *buffer = self.decoded_cigar()?;
                }
                if buffer.iter().any(|&(kind, length)| kind > 8 || length == 0) {
                    return Err(RsomicsError::InvalidInput(format!(
                        "read {}: invalid CIGAR operation",
                        String::from_utf8_lossy(self.name())
                    )));
                }
                Ok(())
            }

            fn sequence_len(&self) -> usize {
                self.sequence_len()
            }

            fn quality_scores(&self) -> &[u8] {
                self.quality_scores()
            }
        }
    };
}

impl_depth_record!(RawRecord);
impl_depth_record!(RecordRef<'_>);

struct Accumulator<'a, W: Write> {
    references: Vec<(Vec<u8>, u64)>,
    selection: Selection,
    bed: Option<Bed>,
    options: Options<'a>,
    output: BufWriter<W>,
    depths: Vec<Vec<u64>>,
    ends: Vec<u64>,
    overlaps: Vec<HashMap<Vec<u8>, u64>>,
    source_coordinates: Vec<Option<(usize, u64)>>,
    cigar: Vec<(u8, u32)>,
    current_reference: Option<usize>,
    output_position: u64,
    line: Vec<u8>,
    summary: Summary,
}

impl<'a, W: Write> Accumulator<'a, W> {
    fn new(
        references: Vec<(Vec<u8>, u64)>,
        selection: Selection,
        bed: Option<Bed>,
        input_count: usize,
        options: Options<'a>,
        output: W,
    ) -> Self {
        Self {
            references,
            selection,
            bed,
            options,
            output: BufWriter::with_capacity(256 * 1024, output),
            depths: vec![vec![0; 2048]; input_count],
            ends: vec![0; input_count],
            overlaps: (0..input_count).map(|_| HashMap::new()).collect(),
            source_coordinates: vec![None; input_count],
            cigar: Vec::with_capacity(16),
            current_reference: None,
            output_position: 0,
            line: Vec::with_capacity(128),
            summary: Summary {
                inputs: input_count,
                positions: 0,
            },
        }
    }

    fn write_header(&mut self, inputs: &[String]) -> Result<()> {
        if self.options.header {
            self.output
                .write_all(b"#CHROM\tPOS")
                .map_err(RsomicsError::Io)?;
            for input in inputs {
                self.output.write_all(b"\t").map_err(RsomicsError::Io)?;
                self.output
                    .write_all(input.as_bytes())
                    .map_err(RsomicsError::Io)?;
            }
            self.output.write_all(b"\n").map_err(RsomicsError::Io)?;
        }
        Ok(())
    }

    fn add(&mut self, source: usize, record: &impl DepthRecord) -> Result<()> {
        let flags = record.flags();
        if flags & self.options.excluded_flags != 0
            || (self.options.included_flags != 0 && flags & self.options.included_flags == 0)
            || flags & self.options.required_flags != self.options.required_flags
            || record.mapping_quality() < self.options.minimum_mapping_quality
            || record.reference_sequence_id() < 0
            || record.alignment_start() < 0
        {
            return Ok(());
        }

        let reference_id = usize::try_from(record.reference_sequence_id()).unwrap();
        let Some((_, reference_length)) = self.references.get(reference_id) else {
            return Err(RsomicsError::InvalidInput(format!(
                "read {} has reference ID {} absent from the header",
                String::from_utf8_lossy(record.name()),
                record.reference_sequence_id()
            )));
        };
        if self
            .selection
            .bounds(reference_id, *reference_length)
            .is_none()
        {
            return Ok(());
        }
        record.copy_cigar(&mut self.cigar)?;
        let used_length = self
            .cigar
            .iter()
            .try_fold(0usize, |length, &(kind, count)| {
                if matches!(kind, MATCH | INSERTION | SEQUENCE_MATCH | SEQUENCE_MISMATCH) {
                    length.checked_add(usize::try_from(count).ok()?)
                } else {
                    Some(length)
                }
            });
        let Some(used_length) = used_length else {
            return Err(invalid_read(record, "used read length overflows"));
        };
        if used_length < self.options.minimum_read_length {
            return Ok(());
        }
        let query_span = self
            .cigar
            .iter()
            .try_fold(0usize, |length, &(kind, count)| {
                if matches!(
                    kind,
                    MATCH | INSERTION | SOFT_CLIP | SEQUENCE_MATCH | SEQUENCE_MISMATCH
                ) {
                    length.checked_add(usize::try_from(count).ok()?)
                } else {
                    Some(length)
                }
            });
        if !self.cigar.is_empty() && query_span != Some(record.sequence_len()) {
            return Err(invalid_read(
                record,
                "CIGAR query span differs from sequence length",
            ));
        }
        let start = u64::try_from(record.alignment_start()).unwrap();
        let coordinate = (reference_id, start);
        if self.source_coordinates[source].is_some_and(|previous| coordinate < previous) {
            return Err(RsomicsError::InvalidInput(
                "alignment input is not coordinate sorted".to_owned(),
            ));
        }
        self.source_coordinates[source] = Some(coordinate);
        let reference_span = self.cigar.iter().try_fold(0u64, |length, &(kind, count)| {
            if matches!(
                kind,
                MATCH | DELETION | REFERENCE_SKIP | SEQUENCE_MATCH | SEQUENCE_MISMATCH
            ) {
                length.checked_add(u64::from(count))
            } else {
                Some(length)
            }
        });
        let Some(reference_span) = reference_span else {
            return Err(invalid_read(record, "CIGAR reference span overflows"));
        };
        let end = start
            .checked_add(reference_span.max(1))
            .ok_or_else(|| invalid_read(record, "alignment end overflows"))?;
        if end > *reference_length {
            return Err(invalid_read(
                record,
                "alignment extends beyond its reference",
            ));
        }

        self.move_to(reference_id, start)?;
        self.flush_until(start)?;
        self.ensure_capacity(end)?;
        self.clear_extension(source, end);
        let overlap_clip = self.overlap_clip(source, record, end);
        self.ends[source] = self.ends[source].max(end);

        let qualities = record.quality_scores();
        let mut reference_position = start;
        let mut query_position = 0usize;
        for cigar_index in 0..self.cigar.len() {
            let (kind, count) = self.cigar[cigar_index];
            let count = usize::try_from(count).unwrap();
            match kind {
                MATCH | SEQUENCE_MATCH | SEQUENCE_MISMATCH => {
                    let run_end = reference_position + u64::try_from(count).unwrap();
                    let run_start = reference_position.max(overlap_clip);
                    if run_start < run_end {
                        let query_start = query_position
                            + usize::try_from(run_start - reference_position).unwrap();
                        if self.options.minimum_base_quality == 0 || qualities.is_empty() {
                            self.increment_range(source, run_start, run_end)?;
                        } else {
                            for (offset, quality) in qualities[query_start..query_position + count]
                                .iter()
                                .copied()
                                .enumerate()
                            {
                                if quality >= self.options.minimum_base_quality {
                                    self.increment(
                                        source,
                                        run_start + u64::try_from(offset).unwrap(),
                                    )?;
                                }
                            }
                        }
                    }
                    reference_position = run_end;
                    query_position += count;
                }
                DELETION => {
                    if self.options.include_deletions {
                        let quality = qualities.get(query_position).copied().unwrap_or(u8::MAX);
                        if quality >= self.options.minimum_base_quality {
                            let run_end = reference_position + u64::try_from(count).unwrap();
                            let run_start = reference_position.max(overlap_clip);
                            if run_start < run_end {
                                self.increment_range(source, run_start, run_end)?;
                            }
                        }
                    }
                    reference_position += u64::try_from(count).unwrap();
                }
                REFERENCE_SKIP => reference_position += u64::try_from(count).unwrap(),
                INSERTION | SOFT_CLIP => query_position += count,
                HARD_CLIP | PADDING => {}
                _ => return Err(invalid_read(record, "unsupported CIGAR operation")),
            }
        }
        Ok(())
    }

    fn move_to(&mut self, reference_id: usize, start: u64) -> Result<()> {
        if let Some(current) = self.current_reference {
            if reference_id < current {
                return Err(RsomicsError::InvalidInput(
                    "merged alignment inputs are not coordinate sorted".to_owned(),
                ));
            }
            if reference_id == current {
                return Ok(());
            }
            self.finish_reference()?;
        }
        if self.options.positions == PositionMode::AllReferences {
            let first = self.current_reference.map_or(0, |current| current + 1);
            for unused in first..reference_id {
                self.emit_empty_reference(unused)?;
            }
        }
        self.start_reference(reference_id, start);
        Ok(())
    }

    fn start_reference(&mut self, reference_id: usize, first_start: u64) {
        self.current_reference = Some(reference_id);
        let reference_length = self.references[reference_id].1;
        let bounds = self
            .selection
            .bounds(reference_id, reference_length)
            .unwrap();
        self.output_position = if self.options.positions == PositionMode::Covered {
            first_start.max(bounds.start)
        } else {
            bounds.start
        };
        self.ends.fill(self.output_position);
        for depth in &mut self.depths {
            depth.fill(0);
        }
        for overlaps in &mut self.overlaps {
            overlaps.clear();
        }
    }

    fn ensure_capacity(&mut self, end: u64) -> Result<()> {
        let required = end.saturating_sub(self.output_position);
        let required = usize::try_from(required).map_err(|_| {
            RsomicsError::InvalidInput("depth buffer exceeds addressable memory".to_owned())
        })?;
        let old_capacity = self.depths[0].len();
        if required < old_capacity {
            return Ok(());
        }
        let mut capacity = old_capacity;
        while required >= capacity {
            capacity = capacity.checked_mul(2).ok_or_else(|| {
                RsomicsError::InvalidInput("depth buffer capacity overflows".to_owned())
            })?;
        }
        let active_end = self
            .ends
            .iter()
            .copied()
            .max()
            .unwrap_or(self.output_position);
        for depth in &mut self.depths {
            let mut expanded = vec![0; capacity];
            for position in self.output_position..active_end {
                expanded[ring_index(position, capacity)] =
                    depth[ring_index(position, old_capacity)];
            }
            *depth = expanded;
        }
        Ok(())
    }

    fn clear_extension(&mut self, source: usize, end: u64) {
        let start = self.ends[source].max(self.output_position);
        let capacity = self.depths[source].len();
        for position in start..end {
            self.depths[source][ring_index(position, capacity)] = 0;
        }
    }

    fn increment(&mut self, source: usize, position: u64) -> Result<()> {
        if position < self.output_position {
            return Ok(());
        }
        let index = ring_index(position, self.depths[source].len());
        self.depths[source][index] += 1;
        Ok(())
    }

    fn increment_range(&mut self, source: usize, start: u64, end: u64) -> Result<()> {
        let length = usize::try_from(end - start).unwrap();
        let capacity = self.depths[source].len();
        let index = ring_index(start, capacity);
        let first_length = length.min(capacity - index);
        for value in &mut self.depths[source][index..index + first_length] {
            *value += 1;
        }
        for value in &mut self.depths[source][..length - first_length] {
            *value += 1;
        }
        Ok(())
    }

    fn overlap_clip(&mut self, source: usize, record: &impl DepthRecord, end: u64) -> u64 {
        if !self.options.remove_overlaps
            || record.flags() & FLAG_PAIRED == 0
            || record.flags() & FLAG_MATE_UNMAPPED != 0
        {
            return 0;
        }
        if let Some(first_end) = self.overlaps[source].remove(record.name()) {
            return first_end;
        }
        let mate_position = record.mate_alignment_start();
        if mate_position < 0
            || (record.mate_reference_sequence_id() == record.reference_sequence_id()
                && u64::try_from(mate_position).unwrap() <= end)
        {
            self.overlaps[source].insert(record.name().to_vec(), end);
        }
        0
    }

    fn flush_until(&mut self, limit: u64) -> Result<()> {
        let reference_id = self.current_reference.unwrap();
        let reference_length = self.references[reference_id].1;
        let bounds = self
            .selection
            .bounds(reference_id, reference_length)
            .unwrap();
        let limit = limit.min(bounds.end);
        while self.output_position < limit {
            let active = self.ends.iter().any(|&end| self.output_position < end);
            if self.options.positions != PositionMode::Covered || active {
                self.emit_position(reference_id, self.output_position)?;
            }
            self.output_position += 1;
        }
        Ok(())
    }

    fn emit_position(&mut self, reference_id: usize, position: u64) -> Result<()> {
        if self
            .bed
            .as_ref()
            .is_some_and(|bed| !bed.contains(reference_id, position))
        {
            return Ok(());
        }
        self.line.clear();
        self.line
            .extend_from_slice(&self.references[reference_id].0);
        self.line.push(b'\t');
        let mut number = itoa::Buffer::new();
        self.line
            .extend_from_slice(number.format(position + 1).as_bytes());
        for (source, depth) in self.depths.iter().enumerate() {
            self.line.push(b'\t');
            let value = if position < self.ends[source] {
                depth[ring_index(position, depth.len())]
            } else {
                0
            };
            self.line.extend_from_slice(number.format(value).as_bytes());
        }
        self.line.push(b'\n');
        self.output
            .write_all(&self.line)
            .map_err(RsomicsError::Io)?;
        self.summary.positions =
            self.summary.positions.checked_add(1).ok_or_else(|| {
                RsomicsError::InvalidInput("position count exceeds u64".to_owned())
            })?;
        Ok(())
    }

    fn finish_reference(&mut self) -> Result<()> {
        let reference_id = self.current_reference.unwrap();
        let reference_length = self.references[reference_id].1;
        let bounds = self
            .selection
            .bounds(reference_id, reference_length)
            .unwrap();
        let end = if self.options.positions == PositionMode::Covered {
            self.ends
                .iter()
                .copied()
                .max()
                .unwrap_or(self.output_position)
        } else {
            bounds.end
        };
        self.flush_until(end)?;
        Ok(())
    }

    fn emit_empty_reference(&mut self, reference_id: usize) -> Result<()> {
        let Some(bounds) = self
            .selection
            .bounds(reference_id, self.references[reference_id].1)
        else {
            return Ok(());
        };
        self.current_reference = Some(reference_id);
        self.output_position = bounds.start;
        self.ends.fill(bounds.start);
        for depth in &mut self.depths {
            depth.fill(0);
        }
        self.flush_until(bounds.end)
    }

    fn finish(mut self) -> Result<Summary> {
        if self.current_reference.is_some() {
            self.finish_reference()?;
        } else if self.options.positions != PositionMode::Covered
            && let Some(reference_id) = self.selection.reference_id
        {
            self.emit_empty_reference(reference_id)?;
        }
        if self.options.positions == PositionMode::AllReferences
            && self.selection.reference_id.is_none()
        {
            let first = self.current_reference.map_or(0, |current| current + 1);
            for reference_id in first..self.references.len() {
                self.emit_empty_reference(reference_id)?;
            }
        }
        self.output.flush().map_err(RsomicsError::Io)?;
        Ok(self.summary)
    }
}

fn invalid_read(record: &impl DepthRecord, reason: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "read {}: {reason}",
        String::from_utf8_lossy(record.name())
    ))
}

fn ring_index(position: u64, capacity: usize) -> usize {
    usize::try_from(position & u64::try_from(capacity - 1).unwrap()).unwrap()
}
