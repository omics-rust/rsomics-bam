mod key;

use std::collections::{HashMap, VecDeque};
use std::io::Write;
use std::path::{Path, PathBuf};

use noodles::sam;
use noodles::sam::header::record::value::map::header::tag as header_tag;
use rsomics_bamio::raw::{RawRecord, RawRecordEncoder};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use self::key::{PairKey, SingleKey};
use crate::hts_quickcheck::{require_bgzf_eof, require_cram_eof};
use crate::{Program, input, md, output};

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const SECONDARY: u16 = 0x100;
const QC_FAIL: u16 = 0x200;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;
const MIN_SCORE_QUALITY: u8 = 15;
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Mode {
    #[default]
    Template,
    Sequence,
}

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub remove: bool,
    pub clear: bool,
    pub include_fails: bool,
    pub mode: Mode,
    pub max_read_length: u32,
    pub additional_threads: Option<usize>,
    pub reference: Option<&'a Path>,
    pub destination: Option<&'a Path>,
    pub program: Option<Program<'a>>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub input: PathBuf,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output: Option<PathBuf>,
    pub records: u64,
    pub written_records: u64,
    pub excluded_records: u64,
    pub examined_records: u64,
    pub paired_records: u64,
    pub single_records: u64,
    pub duplicate_pair_records: u64,
    pub duplicate_single_records: u64,
    pub additional_threads: usize,
}

#[derive(Default)]
struct Counts {
    records: u64,
    written_records: u64,
    excluded_records: u64,
    examined_records: u64,
    paired_records: u64,
    single_records: u64,
    duplicate_pair_records: u64,
    duplicate_single_records: u64,
}

struct Entry {
    record: RawRecord,
    reference: i32,
    position: i64,
    single_key: Option<SingleKey>,
    pair_key: Option<PairKey>,
}

struct Marker<'a, W>
where
    W: Write + Send + 'static,
{
    writer: &'a mut output::Writer<W>,
    remove: bool,
    clear: bool,
    include_fails: bool,
    mode: Mode,
    max_read_length: i64,
    buffer: VecDeque<Entry>,
    base: usize,
    single: HashMap<SingleKey, usize>,
    pair: HashMap<PairKey, usize>,
    counts: Counts,
    previous_reference: i32,
    previous_position: i32,
}

impl<'a, W> Marker<'a, W>
where
    W: Write + Send + 'static,
{
    fn new(writer: &'a mut output::Writer<W>, options: &Options<'_>) -> Self {
        Self {
            writer,
            remove: options.remove,
            clear: options.clear,
            include_fails: options.include_fails,
            mode: options.mode,
            max_read_length: i64::from(options.max_read_length),
            buffer: VecDeque::new(),
            base: 0,
            single: HashMap::new(),
            pair: HashMap::new(),
            counts: Counts::default(),
            previous_reference: 0,
            previous_position: 0,
        }
    }

    fn process(&mut self, mut record: RawRecord) -> Result<bool> {
        self.counts.records = increment(self.counts.records)?;
        self.require_coordinate_order(&record)?;
        if self.clear && record.flags() & DUPLICATE != 0 {
            record.clear_flag_bits(DUPLICATE);
            record.remove_aux(*b"do");
            record.remove_aux(*b"dt");
        }
        let reference = record.reference_sequence_id();
        let raw_position = record.alignment_start();
        let mut excluded_flags = SECONDARY | SUPPLEMENTARY | UNMAPPED;
        if !self.include_fails {
            excluded_flags |= QC_FAIL;
        }
        let excluded = record.flags() & excluded_flags != 0;
        if excluded {
            self.counts.excluded_records = increment(self.counts.excluded_records)?;
            self.buffer.push_back(Entry {
                record,
                reference,
                position: i64::from(raw_position),
                single_key: None,
                pair_key: None,
            });
            self.flush(i64::from(raw_position), reference)?;
            return Ok(true);
        }

        self.counts.examined_records = increment(self.counts.examined_records)?;
        let single_key = key::single(&record)?;
        let position = key::coordinate(single_key);
        let paired = has_mate(&record);
        let id = self
            .base
            .checked_add(self.buffer.len())
            .ok_or_else(count_overflow)?;
        let mut entry = Entry {
            record,
            reference,
            position,
            single_key: None,
            pair_key: None,
        };
        self.resolve_single(&mut entry, single_key, id, paired)?;
        if paired {
            self.counts.paired_records = increment(self.counts.paired_records)?;
            let pair_key = key::pair(&entry.record, self.mode)?;
            self.resolve_pair(&mut entry, pair_key, id)?;
        } else {
            self.counts.single_records = increment(self.counts.single_records)?;
        }
        self.buffer.push_back(entry);
        self.flush(i64::from(raw_position), reference)?;
        Ok(true)
    }

    fn resolve_single(
        &mut self,
        entry: &mut Entry,
        key: SingleKey,
        id: usize,
        paired: bool,
    ) -> Result<()> {
        let Some(occupant) = self.single.get(&key).copied() else {
            self.single.insert(key, id);
            entry.single_key = Some(key);
            return Ok(());
        };
        let occupant_paired = has_mate(&self.entry(occupant).record);
        if paired {
            if !occupant_paired {
                self.mark_existing_single(occupant)?;
                self.entry_mut(occupant).single_key = None;
                self.single.insert(key, id);
                entry.single_key = Some(key);
            }
        } else if occupant_paired {
            mark(&mut entry.record);
            self.counts.duplicate_single_records = increment(self.counts.duplicate_single_records)?;
        } else if score(&entry.record) > score(&self.entry(occupant).record) {
            self.mark_existing_single(occupant)?;
            self.entry_mut(occupant).single_key = None;
            self.single.insert(key, id);
            entry.single_key = Some(key);
        } else {
            mark(&mut entry.record);
            self.counts.duplicate_single_records = increment(self.counts.duplicate_single_records)?;
        }
        Ok(())
    }

    fn resolve_pair(&mut self, entry: &mut Entry, key: PairKey, id: usize) -> Result<()> {
        let Some(occupant) = self.pair.get(&key).copied() else {
            self.pair.insert(key, id);
            entry.pair_key = Some(key);
            return Ok(());
        };
        let old_qc_fail = self.entry(occupant).record.flags() & QC_FAIL != 0;
        let new_qc_fail = entry.record.flags() & QC_FAIL != 0;
        let (old_score, new_score) = if old_qc_fail != new_qc_fail {
            if old_qc_fail { (0, 1) } else { (1, 0) }
        } else {
            (
                score(&self.entry(occupant).record)
                    .checked_add(mate_score(&self.entry(occupant).record)?)
                    .ok_or_else(score_overflow)?,
                score(&entry.record)
                    .checked_add(mate_score(&entry.record)?)
                    .ok_or_else(score_overflow)?,
            )
        };
        let incoming_wins = new_score > old_score
            || (new_score == old_score && entry.record.name() < self.entry(occupant).record.name());
        if incoming_wins {
            mark(&mut self.entry_mut(occupant).record);
            self.entry_mut(occupant).pair_key = None;
            self.pair.insert(key, id);
            entry.pair_key = Some(key);
        } else {
            mark(&mut entry.record);
        }
        self.counts.duplicate_pair_records = increment(self.counts.duplicate_pair_records)?;
        Ok(())
    }

    fn mark_existing_single(&mut self, id: usize) -> Result<()> {
        mark(&mut self.entry_mut(id).record);
        self.counts.duplicate_single_records = increment(self.counts.duplicate_single_records)?;
        Ok(())
    }

    fn flush(&mut self, current_position: i64, current_reference: i32) -> Result<()> {
        while let Some(front) = self.buffer.front() {
            let boundary = front
                .position
                .checked_add(self.max_read_length)
                .ok_or_else(coordinate_overflow)?;
            if boundary > current_position
                && front.reference == current_reference
                && (current_reference != -1 || current_position != -1)
            {
                break;
            }
            self.emit_front()?;
        }
        Ok(())
    }

    fn emit_front(&mut self) -> Result<()> {
        let entry = self
            .buffer
            .pop_front()
            .ok_or_else(|| RsomicsError::InvalidInput("empty markdup buffer".to_owned()))?;
        self.base = self.base.checked_add(1).ok_or_else(count_overflow)?;
        if let Some(key) = entry.single_key {
            self.single.remove(&key);
        }
        if let Some(key) = entry.pair_key {
            self.pair.remove(&key);
        }
        if !self.remove || entry.record.flags() & DUPLICATE == 0 {
            self.writer.write_owned_raw_record(&entry.record)?;
            self.counts.written_records = increment(self.counts.written_records)?;
        }
        Ok(())
    }

    fn finish(mut self) -> Result<Counts> {
        while !self.buffer.is_empty() {
            self.emit_front()?;
        }
        Ok(self.counts)
    }

    fn entry(&self, id: usize) -> &Entry {
        &self.buffer[id - self.base]
    }

    fn entry_mut(&mut self, id: usize) -> &mut Entry {
        &mut self.buffer[id - self.base]
    }

    fn require_coordinate_order(&mut self, record: &RawRecord) -> Result<()> {
        let reference = record.reference_sequence_id();
        let position = record.alignment_start();
        if reference >= 0
            && (reference < self.previous_reference
                || (reference == self.previous_reference && position < self.previous_position))
        {
            return Err(RsomicsError::InvalidInput(
                "markdup input is not in coordinate order".to_owned(),
            ));
        }
        self.previous_reference = reference;
        self.previous_position = position;
        Ok(())
    }
}

pub fn write<W>(input_path: &Path, options: Options<'_>, output: W) -> Result<Summary>
where
    W: Write + Send + 'static,
{
    let additional_threads = options
        .additional_threads
        .unwrap_or_else(crate::sort::default_additional_threads);
    if additional_threads > 256 {
        return Err(RsomicsError::ConfigError(
            "markdup additional thread count cannot exceed 256".to_owned(),
        ));
    }
    if options.max_read_length == 0 {
        return Err(RsomicsError::ConfigError(
            "markdup maximum read length must be greater than zero".to_owned(),
        ));
    }
    let named_format = validate_named_input(input_path)?;
    let input_threads = match named_format {
        Some(input::Format::Bam) => additional_threads,
        _ => 0,
    };
    let mut reader = input::open(input_path, options.reference, input_threads)?;
    let input_format = reader.format();
    let header = reader.read_header(input_path)?;
    reject_query_name_order(&header, input_path)?;
    let mut header = header;
    if let Some(program) = options.program {
        program.add_to(&mut header)?;
    }
    let mut writer = output::Writer::new(
        output::Format::Bam,
        output::Compression::Default,
        additional_threads,
        output,
    );
    writer.write_header(&header)?;
    let mut marker = Marker::new(&mut writer, &options);
    if input_format == input::Format::Cram {
        let mut reference = options
            .reference
            .map(md::ReferenceCache::open)
            .transpose()?;
        let mut raw_encoder = RawRecordEncoder::new();
        let mut completed_encoder = RawRecordEncoder::new();
        reader.visit_records(&header, input_path, |record| {
            let record = complete_cram_raw(
                &header,
                record,
                reference.as_mut(),
                &mut raw_encoder,
                &mut completed_encoder,
            )?;
            marker.process(record)
        })?;
    } else {
        reader.visit_owned_raw_records(&header, input_path, |record| marker.process(record))?;
    }
    let counts = marker.finish()?;
    writer.finish(&header)?;
    Ok(Summary {
        input: input_path.to_path_buf(),
        output: options.destination.map(Path::to_path_buf),
        records: counts.records,
        written_records: counts.written_records,
        excluded_records: counts.excluded_records,
        examined_records: counts.examined_records,
        paired_records: counts.paired_records,
        single_records: counts.single_records,
        duplicate_pair_records: counts.duplicate_pair_records,
        duplicate_single_records: counts.duplicate_single_records,
        additional_threads,
    })
}

fn complete_cram_raw(
    header: &sam::Header,
    record: &dyn sam::alignment::Record,
    reference: Option<&mut md::ReferenceCache>,
    raw_encoder: &mut RawRecordEncoder,
    completed_encoder: &mut RawRecordEncoder,
) -> Result<RawRecord> {
    let mut raw = raw_encoder.encode(header, record)?;
    if raw.flags() & UNMAPPED != 0 {
        let completed = md::complete(header, record, reference)?;
        return completed_encoder.encode(header, &completed);
    }
    if raw.aux_type(*b"MD").is_some() && raw.aux_type(*b"NM").is_some() {
        return Ok(raw);
    }
    let completed = md::complete(header, record, reference)?;
    let completed = completed_encoder.encode(header, &completed)?;
    for tag in [*b"MD", *b"NM"] {
        if raw.aux_type(tag).is_some() {
            continue;
        }
        if let (Some(type_code), Some(value)) = (completed.aux_type(tag), completed.aux_value(tag))
        {
            raw.append_aux(tag, type_code, value)?;
        }
    }
    Ok(raw)
}

fn validate_named_input(path: &Path) -> Result<Option<input::Format>> {
    if path == Path::new("-") {
        return Ok(None);
    }
    let format = input::detect_format(path)?;
    match format {
        input::Format::Bam | input::Format::Sam if input::is_bgzf(path)? => {
            require_bgzf_eof(path)?;
        }
        input::Format::Cram => require_cram_eof(path)?,
        input::Format::Bam | input::Format::Sam => {}
    }
    Ok(Some(format))
}

fn reject_query_name_order(header: &sam::Header, path: &Path) -> Result<()> {
    let order = header
        .header()
        .and_then(|header| header.other_fields().get(&header_tag::SORT_ORDER));
    if order.is_some_and(|value| value.as_slice() == b"queryname") {
        return Err(RsomicsError::InvalidInput(format!(
            "input {} is query-name sorted; markdup requires coordinate order",
            path.display()
        )));
    }
    Ok(())
}

fn has_mate(record: &RawRecord) -> bool {
    record.flags() & PAIRED != 0
        && record.flags() & MATE_UNMAPPED == 0
        && !(record.mate_reference_sequence_id() == -1 && record.mate_alignment_start() == -1)
}

fn score(record: &RawRecord) -> i64 {
    record
        .quality_scores()
        .iter()
        .filter(|&&quality| quality >= MIN_SCORE_QUALITY)
        .map(|&quality| i64::from(quality))
        .sum()
}

fn mate_score(record: &RawRecord) -> Result<i64> {
    match crate::raw_aux::integer(record, *b"ms") {
        crate::raw_aux::Integer::Value(value) => Some(value),
        crate::raw_aux::Integer::Missing | crate::raw_aux::Integer::Invalid => None,
    }
    .ok_or_else(missing_mate_score)
}

fn mark(record: &mut RawRecord) {
    record.set_flag_bits(DUPLICATE);
}

fn increment(value: u64) -> Result<u64> {
    value.checked_add(1).ok_or_else(count_overflow)
}

fn missing_mate_score() -> RsomicsError {
    RsomicsError::InvalidInput(
        "paired duplicate comparison requires an integer ms tag from fixmate -m".to_owned(),
    )
}

fn count_overflow() -> RsomicsError {
    RsomicsError::InvalidInput("markdup record count overflowed".to_owned())
}

fn score_overflow() -> RsomicsError {
    RsomicsError::InvalidInput("markdup score overflowed".to_owned())
}

fn coordinate_overflow() -> RsomicsError {
    RsomicsError::InvalidInput("markdup coordinate arithmetic overflowed".to_owned())
}
