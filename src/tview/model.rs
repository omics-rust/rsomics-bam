use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashSet};
use std::ops::RangeInclusive;
use std::path::Path;

use noodles::core::{Position, Region};
use noodles::fasta;
use noodles::sam;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::cigar::{Op, op::Kind};
use noodles::sam::alignment::record::{Flags, MappingQuality};
use noodles::sam::alignment::record_buf::{
    Cigar, Data, QualityScores, Sequence, data::field::Value,
};
use noodles::sam::header::record::value::map::read_group::tag;
use rsomics_bamio::raw::{RawRecord, RawRecordEncoder};
use rsomics_common::{Result, RsomicsError};
use rsomics_pileup::{PileupEngine, PileupOptions};
use rust_htslib::bam as hts_bam;
use rust_htslib::bam::Read as _;

use super::grid::{GridBuilder, own_column};
use super::rows::{ReadState, RowPacker};
use super::{Options, Settings, Viewport};
use crate::input;

const MAXIMUM_DEPTH: usize = 8000;

#[derive(Default)]
struct DepthLimiter {
    ends: BinaryHeap<Reverse<i64>>,
    last_start: Option<i64>,
}

impl DepthLimiter {
    fn accept(&mut self, start: i64, end: i64) -> bool {
        while self.ends.peek().is_some_and(|end| end.0 < start) {
            self.ends.pop();
        }
        let repeated_start = self.last_start == Some(start);
        self.last_start = Some(start);
        if repeated_start && self.ends.len() >= MAXIMUM_DEPTH {
            return false;
        }
        self.ends.push(Reverse(end));
        true
    }
}

pub(super) fn load(
    input_path: &Path,
    options: Options<'_>,
    settings: Settings,
) -> Result<Viewport> {
    let mut hts_reader = None;
    let mut reader = None;
    let header = if options.additional_threads > 0 {
        let (source, header) = open_threaded(
            input_path,
            options.index,
            options.reference,
            options.additional_threads,
        )?;
        hts_reader = Some(source);
        header
    } else {
        let mut source = if let Some(index) = options.index {
            input::open_indexed_with_index(input_path, index, options.reference)?
        } else {
            input::open_indexed(input_path, options.reference)?
        };
        let header = source.read_header(input_path)?;
        reader = Some(source);
        header
    };
    let references = header
        .reference_sequences()
        .iter()
        .map(|(name, sequence)| {
            (
                name.to_vec(),
                u64::try_from(usize::from(sequence.length())).unwrap(),
            )
        })
        .collect::<Vec<_>>();
    let (reference_id, reference_name, start) = resolve_position(options.position, &references)?;
    let reference_length = references[reference_id].1;
    let start_zero = start - 1;
    if start_zero >= reference_length {
        return Err(RsomicsError::InvalidInput(format!(
            "tview position {reference_name}:{start} is outside the reference"
        )));
    }
    let query_end = start_zero
        .saturating_add(u64::try_from(options.width).unwrap())
        .min(reference_length - 1);
    let region = Region::new(
        reference_name.as_bytes(),
        position_range(start_zero, query_end)?,
    );
    let selected_groups = options
        .sample
        .map(|sample| sample_groups(&header, sample))
        .transpose()?;
    let mut pileup = PileupEngine::with_record_state(
        references.iter().map(|(_, length)| *length),
        PileupOptions::default(),
    );
    let mut records = if let Some(source) = hts_reader.as_mut() {
        threaded_records(source, input_path, reference_id, start_zero, query_end)?
    } else {
        let mut records = Vec::new();
        reader
            .as_mut()
            .expect("one indexed reader is configured")
            .visit_region(&header, input_path, Some(&region), |record| {
                records.push(
                    RecordBuf::try_from_alignment_record(&header, record).map_err(|error| {
                        RsomicsError::InvalidInput(format!(
                            "reading tview record from {}: {error}",
                            input_path.display()
                        ))
                    })?,
                );
                Ok(true)
            })?;
        records
    };
    let coordinate_origin = records
        .iter()
        .filter_map(|record| record.alignment_start())
        .map(|position| u64::try_from(usize::from(position) - 1).unwrap())
        .min()
        .unwrap_or(start_zero);
    let mut encoder = RawRecordEncoder::new();
    let mut depth = DepthLimiter::default();
    for mut source in records.drain(..) {
        translate_record(&mut source, coordinate_origin)?;
        let record = encoder.encode(&header, &source)?;
        if selected_groups
            .as_ref()
            .is_some_and(|groups| !record_in_sample(&record, groups))
        {
            continue;
        }
        let end = alignment_end(&record)?;
        let start = i64::from(record.alignment_start());
        if !depth.accept(start, end) {
            continue;
        }
        pileup
            .push_with_state(record, ReadState { start, end })
            .map_err(|error| {
                RsomicsError::InvalidInput(format!(
                    "building tview pileup from {}: {error}",
                    input_path.display()
                ))
            })?;
    }
    pileup.finish().map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "finishing tview pileup from {}: {error}",
            input_path.display()
        ))
    })?;

    let reference = load_reference(options.reference, &reference_name, start_zero, query_end)?;
    let mut grid = GridBuilder::new(
        reference_name,
        reference_length,
        references
            .iter()
            .map(|(name, _)| String::from_utf8_lossy(name).into_owned())
            .collect(),
        start,
        options.width,
        settings,
        reference,
    )?;
    let mut row_packer = RowPacker::default();
    pileup.drain(|column| {
        let rows = row_packer.pack(column)?;
        let original_position = u64::try_from(column.position()).unwrap() + coordinate_origin;
        if column.reference_id() == i32::try_from(reference_id).unwrap()
            && original_position >= start_zero
            && original_position <= query_end
        {
            let position = i64::try_from(original_position).map_err(|_| {
                RsomicsError::InvalidInput(
                    "tview coordinate exceeds signed 64-bit range".to_owned(),
                )
            })?;
            grid.column(own_column(column, position, &rows))?;
        }
        Ok::<_, RsomicsError>(())
    })?;
    grid.finish()
}

fn open_threaded(
    input: &Path,
    index: Option<&Path>,
    reference: Option<&Path>,
    threads: usize,
) -> Result<(hts_bam::IndexedReader, sam::Header)> {
    let mut reader = match index {
        Some(index) => hts_bam::IndexedReader::from_path_and_index(input, index),
        None => hts_bam::IndexedReader::from_path(input),
    }
    .map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "opening indexed alignment {}: {error}",
            input.display()
        ))
    })?;
    if let Some(reference) = reference {
        reader.set_reference(reference).map_err(|error| {
            RsomicsError::ConfigError(format!(
                "attaching reference {} to {}: {error}",
                reference.display(),
                input.display()
            ))
        })?;
    }
    reader.set_threads(threads).map_err(|error| {
        RsomicsError::ConfigError(format!(
            "configuring {threads} tview decoding threads for {}: {error}",
            input.display()
        ))
    })?;
    let header = sam::io::Reader::new(reader.header().as_bytes())
        .read_header()
        .map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "reading alignment header from {}: {error}",
                input.display()
            ))
        })?;
    Ok((reader, header))
}

fn threaded_records(
    reader: &mut hts_bam::IndexedReader,
    input: &Path,
    reference_id: usize,
    start: u64,
    end: u64,
) -> Result<Vec<RecordBuf>> {
    reader
        .fetch((
            u32::try_from(reference_id).unwrap(),
            i64::try_from(start).map_err(|_| {
                RsomicsError::InvalidInput("tview start exceeds signed 64-bit range".to_owned())
            })?,
            i64::try_from(end + 1).map_err(|_| {
                RsomicsError::InvalidInput("tview end exceeds signed 64-bit range".to_owned())
            })?,
        ))
        .map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "querying tview window from {}: {error}",
                input.display()
            ))
        })?;
    reader
        .records()
        .map(|record| {
            let record = record.map_err(|error| {
                RsomicsError::InvalidInput(format!(
                    "reading tview record from {}: {error}",
                    input.display()
                ))
            })?;
            hts_record_buf(&record)
        })
        .collect()
}

fn hts_record_buf(record: &hts_bam::Record) -> Result<RecordBuf> {
    let mut builder = RecordBuf::builder()
        .set_name(record.qname().to_vec())
        .set_flags(Flags::from(record.flags()))
        .set_cigar(hts_cigar(record))
        .set_sequence(Sequence::from(record.seq().as_bytes()))
        .set_data(hts_data(record));
    if record.tid() >= 0 {
        builder = builder.set_reference_sequence_id(usize::try_from(record.tid()).unwrap());
    }
    if record.pos() >= 0 {
        builder = builder.set_alignment_start(hts_position(record.pos())?);
    }
    if let Some(mapping_quality) = MappingQuality::new(record.mapq()) {
        builder = builder.set_mapping_quality(mapping_quality);
    }
    if record.mtid() >= 0 {
        builder = builder.set_mate_reference_sequence_id(usize::try_from(record.mtid()).unwrap());
    }
    if record.mpos() >= 0 {
        builder = builder.set_mate_alignment_start(hts_position(record.mpos())?);
    }
    if !record.qual().iter().all(|quality| *quality == 255) {
        let quality_scores: QualityScores = record.qual().iter().copied().collect();
        builder = builder.set_quality_scores(quality_scores);
    }
    Ok(builder.build())
}

fn hts_position(position: i64) -> Result<Position> {
    let position = position
        .checked_add(1)
        .and_then(|position| usize::try_from(position).ok())
        .ok_or_else(|| {
            RsomicsError::InvalidInput("alignment position exceeds this platform".to_owned())
        })?;
    Position::try_from(position)
        .map_err(|error| RsomicsError::InvalidInput(format!("reading alignment position: {error}")))
}

fn hts_cigar(record: &hts_bam::Record) -> Cigar {
    record
        .cigar()
        .iter()
        .map(|operation| {
            let kind = match operation {
                hts_bam::record::Cigar::Match(_) => Kind::Match,
                hts_bam::record::Cigar::Ins(_) => Kind::Insertion,
                hts_bam::record::Cigar::Del(_) => Kind::Deletion,
                hts_bam::record::Cigar::RefSkip(_) => Kind::Skip,
                hts_bam::record::Cigar::SoftClip(_) => Kind::SoftClip,
                hts_bam::record::Cigar::HardClip(_) => Kind::HardClip,
                hts_bam::record::Cigar::Pad(_) => Kind::Pad,
                hts_bam::record::Cigar::Equal(_) => Kind::SequenceMatch,
                hts_bam::record::Cigar::Diff(_) => Kind::SequenceMismatch,
            };
            Op::new(kind, operation.len() as usize)
        })
        .collect()
}

fn hts_data(record: &hts_bam::Record) -> Data {
    use noodles::sam::alignment::record::data::field::Tag;
    use rust_htslib::bam::record::Aux;

    [
        (b"RG".as_slice(), Tag::READ_GROUP),
        (b"CS".as_slice(), Tag::COLOR_SEQUENCE),
        (b"CQ".as_slice(), Tag::COLOR_QUALITY_SCORES),
    ]
    .into_iter()
    .filter_map(|(source, target)| match record.aux(source) {
        Ok(Aux::String(value)) => Some((target, Value::from(value))),
        _ => None,
    })
    .collect()
}

fn translate_record(record: &mut RecordBuf, origin: u64) -> Result<()> {
    let position = record.alignment_start().ok_or_else(|| {
        RsomicsError::InvalidInput("indexed tview query returned an unmapped record".to_owned())
    })?;
    *record.alignment_start_mut() = Some(local_position(position, origin)?);
    if let Some(position) = record.mate_alignment_start() {
        *record.mate_alignment_start_mut() = local_position(position, origin).ok();
    }
    Ok(())
}

fn local_position(position: Position, origin: u64) -> Result<Position> {
    let position = u64::try_from(usize::from(position) - 1).unwrap();
    let local = position
        .checked_sub(origin)
        .and_then(|position| position.checked_add(1))
        .and_then(|position| usize::try_from(position).ok())
        .ok_or_else(|| {
            RsomicsError::InvalidInput("tview record is outside the translated window".to_owned())
        })?;
    Position::try_from(local).map_err(|error| {
        RsomicsError::InvalidInput(format!("translating tview record coordinate: {error}"))
    })
}

fn resolve_position(
    value: Option<&str>,
    references: &[(Vec<u8>, u64)],
) -> Result<(usize, String, u64)> {
    let first = references.first().ok_or_else(|| {
        RsomicsError::InvalidInput("tview input header has no reference sequences".to_owned())
    })?;
    let (name, start) = match value {
        None => (String::from_utf8_lossy(&first.0).into_owned(), 1),
        Some(value) => {
            let (name, start) = value.split_once(':').unwrap_or((value, "1"));
            if name.is_empty() || start.is_empty() || start.contains('-') {
                return Err(invalid_position(value));
            }
            let digits = start.replace(',', "");
            let start = digits
                .parse::<u64>()
                .ok()
                .filter(|start| *start > 0)
                .ok_or_else(|| invalid_position(value))?;
            (name.to_owned(), start)
        }
    };
    let reference_id = references
        .iter()
        .position(|(candidate, _)| candidate.as_slice() == name.as_bytes())
        .ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "tview reference {name:?} is absent from the alignment header"
            ))
        })?;
    Ok((reference_id, name, start))
}

fn invalid_position(value: &str) -> RsomicsError {
    RsomicsError::ConfigError(format!(
        "invalid tview position {value:?}; expected REFERENCE or REFERENCE:POSITION"
    ))
}

fn position_range(start: u64, end: u64) -> Result<RangeInclusive<Position>> {
    let start = Position::try_from(usize::try_from(start + 1).unwrap())
        .map_err(|error| RsomicsError::ConfigError(format!("invalid tview start: {error}")))?;
    let end = Position::try_from(usize::try_from(end + 1).unwrap())
        .map_err(|error| RsomicsError::ConfigError(format!("invalid tview end: {error}")))?;
    Ok(start..=end)
}

fn sample_groups(header: &noodles::sam::Header, sample: &str) -> Result<HashSet<Vec<u8>>> {
    let requested = sample
        .split(',')
        .filter(|value| !value.is_empty())
        .map(str::as_bytes)
        .collect::<HashSet<_>>();
    let groups = header
        .read_groups()
        .iter()
        .filter_map(|(name, group)| {
            let id: &[u8] = name.as_ref();
            let selected = requested.contains(id)
                || group
                    .other_fields()
                    .get(&tag::SAMPLE)
                    .is_some_and(|value| requested.contains(AsRef::<[u8]>::as_ref(value)));
            selected.then(|| name.to_vec())
        })
        .collect::<HashSet<_>>();
    if groups.is_empty() {
        return Err(RsomicsError::InvalidInput(format!(
            "tview sample {sample:?} is absent from the alignment header"
        )));
    }
    Ok(groups)
}

fn record_in_sample(record: &RawRecord, groups: &HashSet<Vec<u8>>) -> bool {
    record
        .aux_value(*b"RG")
        .and_then(|value| value.strip_suffix(&[0]))
        .is_some_and(|group| groups.contains(group))
}
fn alignment_end(record: &RawRecord) -> Result<i64> {
    let span = record
        .decoded_cigar()?
        .into_iter()
        .filter(|(kind, _)| matches!(kind, 0 | 2 | 3 | 7 | 8))
        .map(|(_, length)| i64::from(length))
        .sum::<i64>();
    Ok(i64::from(record.alignment_start()) + span - 1)
}

fn load_reference(
    path: Option<&Path>,
    name: &str,
    start: u64,
    end: u64,
) -> Result<Option<Vec<u8>>> {
    let Some(path) = path else {
        return Ok(None);
    };
    let mut reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "opening indexed reference {}: {error}",
                path.display()
            ))
        })?;
    let region = Region::new(name.as_bytes(), position_range(start, end)?);
    let record = reader.query(&region).map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "reading reference region {region} from {}: {error}",
            path.display()
        ))
    })?;
    Ok(Some(record.sequence().as_ref().to_vec()))
}
