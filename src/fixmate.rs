use std::io::Write;
use std::path::{Path, PathBuf};

use noodles::sam;
use noodles::sam::header::record::value::map::header::tag as header_tag;
use rsomics_bamio::raw::{RawRecord, RawRecordEncoder};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::hts_quickcheck::{require_bgzf_eof, require_cram_eof};
use crate::{Program, input, md, output};

const PAIRED: u16 = 0x1;
const PROPER_PAIR: u16 = 0x2;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const REVERSE: u16 = 0x10;
const MATE_REVERSE: u16 = 0x20;
const READ2: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const SUPPLEMENTARY: u16 = 0x800;
const MATE_SCORE_MIN_QUALITY: u8 = 15;

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub mate_score: bool,
    pub remove: bool,
    pub proper_pair_check: bool,
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
    pub templates: u64,
    pub paired_templates: u64,
    pub additional_threads: usize,
}

#[derive(Default)]
struct Counts {
    records: u64,
    written_records: u64,
    templates: u64,
    paired_templates: u64,
}

struct MateSnapshot {
    flags: u16,
    reference: i32,
    position: i32,
    mapping_quality: u8,
    end: i64,
    cigar: Vec<u8>,
    score: Option<u32>,
}

impl MateSnapshot {
    fn new(record: &RawRecord, score: bool) -> Result<Self> {
        let (span, cigar) = if record.aux_type(*b"CG") == Some(b'B') {
            cigar_fields(record.decoded_cigar()?)?
        } else {
            cigar_fields(record.cigar_ops())?
        };
        let end = i64::from(record.alignment_start())
            .checked_add(span.max(1))
            .ok_or_else(cigar_overflow)?;
        Ok(Self {
            flags: record.flags(),
            reference: record.reference_sequence_id(),
            position: record.alignment_start(),
            mapping_quality: record.mapping_quality(),
            end,
            cigar,
            score: score.then(|| mate_score(record)),
        })
    }

    fn flag_set(&self, bits: u16) -> bool {
        self.flags & bits != 0
    }

    fn end_position(&self) -> i64 {
        self.end
    }

    fn five_prime(&self) -> i64 {
        if self.flag_set(REVERSE) {
            self.end_position()
        } else {
            i64::from(self.position)
        }
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
            "fixmate additional thread count cannot exceed 256".to_owned(),
        ));
    }

    let named_format = validate_named_input(input_path)?;
    let input_threads = match named_format {
        Some(input::Format::Bam) => additional_threads,
        _ => 0,
    };
    let mut reader = input::open(input_path, options.reference, input_threads)?;
    let input_format = reader.format();
    let mut header = reader.read_header(input_path)?;
    reject_coordinate_order(&header, input_path)?;
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

    let mut counts = Counts::default();
    let mut group = Vec::new();
    let mut ingest = |record: RawRecord| {
        if group
            .first()
            .is_some_and(|first: &RawRecord| first.name() != record.name())
        {
            process_group(&mut group, &options, &mut writer, &mut counts)?;
        }
        group.push(record);
        Ok(true)
    };
    if input_format == input::Format::Cram {
        let mut reference = options
            .reference
            .map(md::ReferenceCache::open)
            .transpose()?;
        let mut encoder = RawRecordEncoder::new();
        reader.visit_records(&header, input_path, |record| {
            let record = md::complete(&header, record, reference.as_mut())?;
            ingest(encoder.encode(&header, &record)?)
        })?;
    } else {
        reader.visit_owned_raw_records(&header, input_path, ingest)?;
    }
    if !group.is_empty() {
        process_group(&mut group, &options, &mut writer, &mut counts)?;
    }
    writer.finish(&header)?;

    Ok(Summary {
        input: input_path.to_path_buf(),
        output: options.destination.map(Path::to_path_buf),
        records: counts.records,
        written_records: counts.written_records,
        templates: counts.templates,
        paired_templates: counts.paired_templates,
        additional_threads,
    })
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

fn reject_coordinate_order(header: &sam::Header, path: &Path) -> Result<()> {
    let order = header
        .header()
        .and_then(|header| header.other_fields().get(&header_tag::SORT_ORDER));
    if order.is_some_and(|value| value.as_slice() == b"coordinate") {
        return Err(RsomicsError::InvalidInput(format!(
            "input {} is coordinate sorted; fixmate requires records grouped by read name",
            path.display()
        )));
    }
    Ok(())
}

fn process_group<W>(
    group: &mut Vec<RawRecord>,
    options: &Options<'_>,
    writer: &mut output::Writer<W>,
    counts: &mut Counts,
) -> Result<()>
where
    W: Write + Send + 'static,
{
    counts.templates = counts.templates.checked_add(1).ok_or_else(count_overflow)?;
    counts.records = counts
        .records
        .checked_add(u64::try_from(group.len()).map_err(|_| count_overflow())?)
        .ok_or_else(count_overflow)?;

    if fix_primary_records(group, options)? {
        counts.paired_templates = counts
            .paired_templates
            .checked_add(1)
            .ok_or_else(count_overflow)?;
    }

    for record in group.iter() {
        if options.remove && flag_set(record, SECONDARY | UNMAPPED) {
            continue;
        }
        writer.write_owned_raw_record(record)?;
        counts.written_records = counts
            .written_records
            .checked_add(1)
            .ok_or_else(count_overflow)?;
    }
    group.clear();
    Ok(())
}

fn fix_primary_records(group: &mut [RawRecord], options: &Options<'_>) -> Result<bool> {
    let primary_indices = group
        .iter()
        .enumerate()
        .filter_map(|(index, record)| {
            (!flag_set(record, SECONDARY | SUPPLEMENTARY)).then_some(index)
        })
        .collect::<Vec<_>>();

    let Some((&first, rest)) = primary_indices.split_first() else {
        return Ok(false);
    };
    if rest.is_empty() {
        fix_orphan(&mut group[first]);
        return Ok(false);
    }

    for &current in rest {
        let (pre, cur) = two_mut(group, first, current);
        fix_pair(pre, cur, options)?;
    }

    let mut primaries = [None, None];
    for &index in &primary_indices {
        primaries[usize::from(flag_set(&group[index], READ2))] = Some(index);
    }
    for index in 0..group.len() {
        if !flag_set(&group[index], SUPPLEMENTARY) || !flag_set(&group[index], PAIRED) {
            continue;
        }
        let mate = primaries[usize::from(!flag_set(&group[index], READ2))];
        if let Some(mate) = mate {
            let source = MateSnapshot::new(&group[mate], false)?;
            sync_mate_inner(&mut group[index], &source);
            sync_mq_mc(&mut group[index], &source)?;
        }
    }
    Ok(true)
}

fn fix_pair(pre: &mut RawRecord, cur: &mut RawRecord, options: &Options<'_>) -> Result<()> {
    pre.set_flag_bits(PAIRED);
    cur.set_flag_bits(PAIRED);

    let pre_unmapped = flag_set(pre, UNMAPPED);
    let cur_unmapped = flag_set(cur, UNMAPPED);
    let mut pre_source = MateSnapshot::new(pre, options.mate_score)?;
    let mut cur_source = MateSnapshot::new(cur, options.mate_score)?;
    let pre_end = (!pre_unmapped).then_some(pre_source.end_position());
    let cur_end = (!cur_unmapped).then_some(cur_source.end_position());

    sync_unmapped_position(cur, &pre_source);
    sync_unmapped_position(pre, &cur_source);

    pre_source.reference = pre.reference_sequence_id();
    pre_source.position = pre.alignment_start();
    cur_source.reference = cur.reference_sequence_id();
    cur_source.position = cur.alignment_start();
    sync_mate_inner(cur, &pre_source);
    sync_mate_inner(pre, &cur_source);
    sync_mq_mc(cur, &pre_source)?;
    sync_mq_mc(pre, &cur_source)?;

    if pre.reference_sequence_id() == cur.reference_sequence_id()
        && !flag_set(pre, UNMAPPED | MATE_UNMAPPED)
        && !flag_set(cur, UNMAPPED | MATE_UNMAPPED)
    {
        let pre_five = if pre_source.flag_set(REVERSE) {
            pre_end.unwrap_or_default()
        } else {
            i64::from(pre_source.position)
        };
        let cur_five = if cur_source.flag_set(REVERSE) {
            cur_end.unwrap_or_default()
        } else {
            i64::from(cur_source.position)
        };
        pre.set_template_length(template_length(cur_five, pre_five)?);
        cur.set_template_length(template_length(pre_five, cur_five)?);
    } else {
        pre.set_template_length(0);
        cur.set_template_length(0);
    }

    if options.proper_pair_check && !plausibly_proper(&pre_source, &cur_source) {
        pre.clear_flag_bits(PROPER_PAIR);
        cur.clear_flag_bits(PROPER_PAIR);
    }
    if options.mate_score {
        set_aux_i(pre, *b"ms", cur_source.score.unwrap_or_default())?;
        set_aux_i(cur, *b"ms", pre_source.score.unwrap_or_default())?;
    }
    if options.remove {
        if pre_unmapped {
            cur.clear_flag_bits(MATE_REVERSE | PROPER_PAIR);
        }
        if cur_unmapped {
            pre.clear_flag_bits(MATE_REVERSE | PROPER_PAIR);
        }
    }
    Ok(())
}

fn sync_unmapped_position(destination: &mut RawRecord, source: &MateSnapshot) {
    if flag_set(destination, UNMAPPED) && !source.flag_set(UNMAPPED) {
        destination.set_reference_sequence_id(source.reference);
        destination.set_alignment_start(source.position);
    }
}

fn sync_mate_inner(destination: &mut RawRecord, source: &MateSnapshot) {
    destination.set_mate_reference_sequence_id(source.reference);
    destination.set_mate_alignment_start(source.position);
    if source.flag_set(REVERSE) {
        destination.set_flag_bits(MATE_REVERSE);
    } else {
        destination.clear_flag_bits(MATE_REVERSE);
    }
    if source.flag_set(UNMAPPED) {
        destination.set_flag_bits(MATE_UNMAPPED);
    }
}

fn sync_mq_mc(destination: &mut RawRecord, source: &MateSnapshot) -> Result<()> {
    let source_mapped = !source.flag_set(UNMAPPED);
    if source_mapped {
        set_aux_i(destination, *b"MQ", u32::from(source.mapping_quality))?;
    }
    if source_mapped || !flag_set(destination, UNMAPPED) {
        destination.set_aux(*b"MC", b'Z', &source.cigar)?;
    }
    Ok(())
}

fn fix_orphan(record: &mut RawRecord) {
    record.set_mate_reference_sequence_id(-1);
    record.set_mate_alignment_start(-1);
    record.set_template_length(0);
    record.clear_flag_bits(MATE_REVERSE | PROPER_PAIR);
}

fn plausibly_proper(left: &MateSnapshot, right: &MateSnapshot) -> bool {
    if left.flag_set(UNMAPPED) || right.flag_set(UNMAPPED) || left.reference != right.reference {
        return false;
    }
    let left_five = left.five_prime();
    let right_five = right.five_prime();
    let (first, second) = if left_five > right_five {
        (right, left)
    } else {
        (left, right)
    };
    !first.flag_set(REVERSE) && second.flag_set(REVERSE)
}

fn template_length(left: i64, right: i64) -> Result<i32> {
    left.checked_sub(right)
        .and_then(|value| i32::try_from(value).ok())
        .ok_or_else(|| {
            RsomicsError::InvalidInput("template length exceeds BAM i32 range".to_owned())
        })
}

fn mate_score(record: &RawRecord) -> u32 {
    let qualities = record.quality_scores();
    if qualities.is_empty() && record.sequence_len() > 0 {
        return u32::try_from(record.sequence_len())
            .unwrap_or(u32::MAX)
            .wrapping_mul(u32::from(u8::MAX));
    }
    qualities
        .iter()
        .filter(|&&quality| quality >= MATE_SCORE_MIN_QUALITY)
        .fold(0u32, |score, &quality| {
            score.wrapping_add(u32::from(quality))
        })
}

fn cigar_fields(operations: impl IntoIterator<Item = (u8, u32)>) -> Result<(i64, Vec<u8>)> {
    let mut span = 0i64;
    let mut cigar = Vec::new();
    let mut number = itoa::Buffer::new();
    for (kind, length) in operations {
        if matches!(kind, 0 | 2 | 3 | 7 | 8) {
            span = span
                .checked_add(i64::from(length))
                .ok_or_else(cigar_overflow)?;
        }
        let operation = *b"MIDNSHP=X".get(usize::from(kind)).ok_or_else(|| {
            RsomicsError::InvalidInput(format!("unsupported BAM CIGAR operation code {kind}"))
        })?;
        cigar.extend_from_slice(number.format(length).as_bytes());
        cigar.push(operation);
    }
    if cigar.is_empty() {
        cigar.push(b'*');
    }
    cigar.push(0);
    Ok((span, cigar))
}

fn set_aux_i(record: &mut RawRecord, tag: [u8; 2], value: u32) -> Result<()> {
    record.set_aux(tag, b'i', &value.to_le_bytes())
}

fn two_mut<T>(slice: &mut [T], left: usize, right: usize) -> (&mut T, &mut T) {
    let (low, high, reverse) = if left < right {
        (left, right, false)
    } else {
        (right, left, true)
    };
    let (before, after) = slice.split_at_mut(high);
    if reverse {
        (&mut after[0], &mut before[low])
    } else {
        (&mut before[low], &mut after[0])
    }
}

fn flag_set(record: &RawRecord, bits: u16) -> bool {
    record.flags() & bits != 0
}

fn cigar_overflow() -> RsomicsError {
    RsomicsError::InvalidInput("CIGAR reference span overflows i64".to_owned())
}

fn count_overflow() -> RsomicsError {
    RsomicsError::InvalidInput("fixmate record count exceeds u64".to_owned())
}
