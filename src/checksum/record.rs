use crc32fast::Hasher;
use rsomics_bamio::raw::RecordRef;
use rsomics_common::Result;

use super::{Options, RecordChecksums, Sanitize, TagSelection};

const FORWARD: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
const REVERSE: &[u8; 16] = b"=TGKCYSBAWRDMHVN";
const FORWARD_PAIRS: [[u8; 2]; 256] = base_pairs(false);
const REVERSE_PAIRS: [[u8; 2]; 256] = base_pairs(true);

const fn base_pairs(reverse: bool) -> [[u8; 2]; 256] {
    let mut pairs = [[0; 2]; 256];
    let mut byte = 0;
    while byte < pairs.len() {
        let high = byte >> 4;
        let low = byte & 0x0f;
        pairs[byte] = if reverse {
            [REVERSE[low], REVERSE[high]]
        } else {
            [FORWARD[high], FORWARD[low]]
        };
        byte += 1;
    }
    pairs
}

#[derive(Default)]
pub(super) struct Scratch {
    sequence: Vec<u8>,
    quality: Vec<u8>,
    auxiliary: Vec<u8>,
    cigar: Vec<(u8, u32)>,
    fields: Vec<AuxField>,
    slots: Vec<Option<AuxField>>,
}

#[derive(Clone, Copy)]
struct AuxField {
    tag: [u8; 2],
    type_code: u8,
    value_start: usize,
    value_end: usize,
}

pub(super) struct ProcessedAlignment<'a> {
    pub checksums: RecordChecksums,
    pub read_group: Option<&'a [u8]>,
    pub qc_fail: bool,
}

pub(super) fn alignment_with_scratch<'a>(
    record: RecordRef<'a>,
    reference_lengths: &[i64],
    options: &Options,
    scratch: &mut Scratch,
) -> Result<Option<ProcessedAlignment<'a>>> {
    let long_cigar = has_long_cigar_encoding(record);
    scratch.cigar.clear();
    if options.check_cigar || !options.sanitize.is_empty() {
        record.decode_cigar_into(&mut scratch.cigar)?;
    }
    let mut fields = SanitizedFields::new(record, &mut scratch.cigar);
    fields.apply(reference_lengths, options.sanitize);
    let flags = fields.flags;
    if flags & options.excluded_flags != 0
        || flags & options.required_flags != options.required_flags
    {
        return Ok(None);
    }

    decode_sequence_and_quality(
        record,
        flags,
        options,
        &mut scratch.sequence,
        &mut scratch.quality,
    );

    let masked_flags = (flags & options.flag_mask) as u8;
    let sequence_crc = crc(&[&[masked_flags], &scratch.sequence]);
    let name_crc = crc(&[record.name(), &[0, masked_flags], &scratch.sequence]);
    let quality_crc = crc_with(sequence_crc, &scratch.quality);
    let (auxiliary_crc, read_group) = auxiliary_crc(
        record,
        &options.tags,
        fields.unmapped_auxiliary_removed,
        long_cigar,
        sequence_crc,
        &mut scratch.auxiliary,
        &mut scratch.fields,
        &mut scratch.slots,
    )?;

    let mut position = sequence_crc;
    if options.check_position {
        position = crc_with(position, &fields.reference_id.to_le_bytes());
        position = crc_with(position, &fields.position.to_le_bytes());
    }
    let mut mate = sequence_crc;
    if options.check_mate {
        mate = crc_with(mate, &record.mate_reference_sequence_id().to_le_bytes());
        mate = crc_with(
            mate,
            &i64::from(record.mate_alignment_start()).to_le_bytes(),
        );
        mate = crc_with(mate, &i64::from(record.template_length()).to_le_bytes());
    }
    let mut cigar = sequence_crc;
    if options.check_cigar {
        cigar = crc_with(cigar, &u32::from(fields.mapping_quality).to_le_bytes());
        for &(kind, length) in fields.cigar.iter() {
            cigar = crc_with(cigar, &((length << 4) | u32::from(kind)).to_le_bytes());
        }
    }

    Ok(Some(ProcessedAlignment {
        checksums: RecordChecksums {
            sequence: sequence_crc,
            name: name_crc,
            quality: quality_crc,
            auxiliary: auxiliary_crc,
            position,
            cigar,
            mate,
        },
        read_group,
        qc_fail: flags & 0x200 != 0,
    }))
}

fn has_long_cigar_encoding(record: RecordRef<'_>) -> bool {
    let mut cigar = record.cigar_ops();
    if cigar.next()
        != u32::try_from(record.sequence_len())
            .ok()
            .map(|length| (4, length))
        || cigar.next().is_none_or(|(kind, _)| kind != 3)
        || cigar.next().is_some()
        || record.aux_type(*b"CG") != Some(b'B')
    {
        return false;
    }
    let Some(value) = record.aux_value(*b"CG") else {
        return false;
    };
    value.len() >= 5
        && matches!(value[0], b'I' | b'i')
        && u32::from_le_bytes(value[1..5].try_into().unwrap()) > u32::from(u16::MAX)
}

fn decode_sequence_and_quality(
    record: RecordRef<'_>,
    flags: u16,
    options: &Options,
    sequence: &mut Vec<u8>,
    quality: &mut Vec<u8>,
) {
    let length = record.sequence_len();
    let packed = record.seq_bytes_packed();
    let bytes = record.payload();
    let quality_start = packed.as_ptr_range().end as usize - bytes.as_ptr() as usize;
    let scores = &bytes[quality_start..quality_start + length];
    sequence.resize(length, 0);
    quality.resize(length, 0);

    if options.reverse_complement && flags & 0x10 != 0 {
        let mut target = 0;
        let mut source = packed.len();
        if !length.is_multiple_of(2) {
            source -= 1;
            sequence[0] = REVERSE[usize::from(packed[source] >> 4)];
            target = 1;
        }
        while source > 0 {
            source -= 1;
            sequence[target..target + 2]
                .copy_from_slice(&REVERSE_PAIRS[usize::from(packed[source])]);
            target += 2;
        }
        for (target, &score) in quality.iter_mut().zip(scores.iter().rev()) {
            *target = score.wrapping_add(33);
        }
    } else {
        for (&source, target) in packed.iter().zip(sequence.chunks_exact_mut(2)) {
            target.copy_from_slice(&FORWARD_PAIRS[usize::from(source)]);
        }
        if !length.is_multiple_of(2) {
            sequence[length - 1] = FORWARD[usize::from(packed[packed.len() - 1] >> 4)];
        }
        for (target, &score) in quality.iter_mut().zip(scores) {
            *target = score.wrapping_add(33);
        }
    }
}

pub(super) fn sequence(
    record: rsomics_seqio::Record<'_>,
    options: &Options,
    scratch: &mut Scratch,
) -> Option<RecordChecksums> {
    let identifier = record
        .id
        .split(|byte| byte.is_ascii_whitespace())
        .next()
        .unwrap_or_default();
    let (name, flags) = if let Some(name) = identifier.strip_suffix(b"/1") {
        (name, 0x4d)
    } else if let Some(name) = identifier.strip_suffix(b"/2") {
        (name, 0x8d)
    } else {
        (identifier, 0x04)
    };
    if flags & options.excluded_flags != 0
        || flags & options.required_flags != options.required_flags
    {
        return None;
    }
    let masked_flags = (flags & options.flag_mask) as u8;
    scratch.sequence.clear();
    scratch
        .sequence
        .extend(record.seq.iter().copied().map(canonical_base));
    let sequence_crc = crc(&[&[masked_flags], &scratch.sequence]);
    let name_crc = crc(&[name, &[0, masked_flags], &scratch.sequence]);
    let quality_crc = if let Some(quality) = record.qual {
        crc_with(sequence_crc, quality)
    } else {
        crc_repeated(sequence_crc, 32, record.seq.len())
    };
    let mut position = sequence_crc;
    if options.check_position {
        position = crc_with(position, &(-1i32).to_le_bytes());
        position = crc_with(position, &(-1i64).to_le_bytes());
    }
    let mut cigar = sequence_crc;
    if options.check_cigar {
        cigar = crc_with(cigar, &0u32.to_le_bytes());
    }
    let mut mate = sequence_crc;
    if options.check_mate {
        mate = crc_with(mate, &(-1i32).to_le_bytes());
        mate = crc_with(mate, &(-1i64).to_le_bytes());
        mate = crc_with(mate, &0i64.to_le_bytes());
    }
    Some(RecordChecksums {
        sequence: sequence_crc,
        name: name_crc,
        quality: quality_crc,
        auxiliary: sequence_crc,
        position,
        cigar,
        mate,
    })
}

#[allow(clippy::too_many_arguments)]
fn auxiliary_crc<'a>(
    record: RecordRef<'a>,
    selection: &TagSelection,
    remove_unmapped: bool,
    remove_long_cigar: bool,
    initial: u32,
    auxiliary: &mut Vec<u8>,
    fields: &mut Vec<AuxField>,
    slots: &mut Vec<Option<AuxField>>,
) -> Result<(u32, Option<&'a [u8]>)> {
    auxiliary.clear();
    fields.clear();
    slots.clear();
    if let TagSelection::Listed(tags) = selection {
        slots.resize(tags.len(), None);
    }
    let mut read_group = None;
    let bytes = record.payload();
    let name_length = usize::from(bytes[8]);
    let cigar_count = usize::from(u16::from_le_bytes([bytes[12], bytes[13]]));
    let mut position = 32
        + name_length
        + 4 * cigar_count
        + record.sequence_len().div_ceil(2)
        + record.sequence_len();
    while position < bytes.len() {
        let tag = [bytes[position], bytes[position + 1]];
        let type_code = bytes[position + 2];
        let value_start = position + 3;
        let value_length = value_length(bytes, value_start, type_code);
        let value_end = value_start + value_length;
        if tag == *b"RG" && type_code == b'Z' {
            read_group = Some(&bytes[value_start..value_end - 1]);
        }
        let removed = remove_long_cigar && tag == *b"CG"
            || remove_unmapped && [*b"NM", *b"MD", *b"CG", *b"SM"].contains(&tag);
        if !removed && tag.iter().all(|byte| (b'0'..=b'z').contains(byte)) {
            let field = AuxField {
                tag,
                type_code,
                value_start,
                value_end,
            };
            match selection {
                TagSelection::Listed(tags) => {
                    if let Some(index) = tags.iter().position(|candidate| *candidate == tag) {
                        slots[index] = Some(field);
                    }
                }
                TagSelection::AllExcept(tags) if !tags.contains(&tag) => fields.push(field),
                TagSelection::AllExcept(_) => {}
            }
        }
        position = value_end;
    }
    match selection {
        TagSelection::Listed(_) => {
            for field in slots.iter().flatten().copied() {
                append_canonical(auxiliary, field, bytes);
            }
        }
        TagSelection::AllExcept(_) => {
            fields.sort_by_key(|field| field.tag);
            for &field in fields.iter() {
                append_canonical(auxiliary, field, bytes);
            }
        }
    }
    Ok((crc_with(initial, auxiliary), read_group))
}

struct SanitizedFields<'a> {
    flags: u16,
    reference_id: i32,
    position: i64,
    mapping_quality: u8,
    cigar: &'a mut Vec<(u8, u32)>,
    unmapped_auxiliary_removed: bool,
}

impl<'a> SanitizedFields<'a> {
    fn new(record: RecordRef<'_>, cigar: &'a mut Vec<(u8, u32)>) -> Self {
        Self {
            flags: record.flags(),
            reference_id: record.reference_sequence_id(),
            position: i64::from(record.alignment_start()),
            mapping_quality: record.mapping_quality(),
            cigar,
            unmapped_auxiliary_removed: false,
        }
    }

    fn apply(&mut self, references: &[i64], sanitize: Sanitize) {
        if sanitize.contains(Sanitize::POSITION) && self.reference_id < 0 {
            self.position = -1;
            if sanitize.contains(Sanitize::UNMAPPED) {
                self.flags |= 0x4;
            }
        }
        if sanitize.contains(Sanitize::CIGAR) && self.flags & 0x4 == 0 {
            if self.position < 0 && sanitize.contains(Sanitize::UNMAPPED) {
                self.flags |= 0x4;
            } else {
                let reference_length = usize::try_from(self.reference_id)
                    .ok()
                    .and_then(|id| references.get(id))
                    .copied()
                    .unwrap_or(-1);
                if self.position >= reference_length && sanitize.contains(Sanitize::UNMAPPED) {
                    self.flags |= 0x4;
                    if sanitize.contains(Sanitize::POSITION) {
                        self.reference_id = -1;
                        self.position = -1;
                    }
                } else if alignment_end(self.position, self.cigar) > reference_length {
                    trim_cigar(self.cigar, self.position, reference_length, &mut self.flags);
                }
            }
        }
        if self.flags & 0x4 != 0 {
            if sanitize.contains(Sanitize::CIGAR) {
                self.cigar.clear();
            }
            if sanitize.contains(Sanitize::MAPPING_QUALITY) {
                self.mapping_quality = 0;
            }
            self.unmapped_auxiliary_removed = sanitize.contains(Sanitize::AUXILIARY);
        }
        if self.flags & 0x4 == 0 && sanitize.contains(Sanitize::CIGAR_EQX) {
            for (kind, _) in self.cigar.iter_mut() {
                if matches!(*kind, 7 | 8) {
                    *kind = 0;
                }
            }
        }
        if self.flags & 0x4 == 0 && sanitize.contains(Sanitize::CIGAR_DUPLICATES) {
            collapse_cigar(self.cigar);
        }
    }
}

fn alignment_end(start: i64, cigar: &[(u8, u32)]) -> i64 {
    cigar
        .iter()
        .filter(|(kind, _)| matches!(*kind, 0 | 2 | 3 | 7 | 8))
        .fold(start, |position, (_, length)| position + i64::from(*length))
}

fn trim_cigar(cigar: &mut Vec<(u8, u32)>, start: i64, end: i64, flags: &mut u16) {
    let mut position = start;
    let mut boundary = None;
    for (index, &(kind, length)) in cigar.iter().enumerate() {
        if !matches!(kind, 0 | 2 | 3 | 7 | 8) {
            continue;
        }
        position += i64::from(length);
        if position > end {
            boundary = Some((index, kind, length, position));
            break;
        }
    }
    let Some((index, kind, length, position)) = boundary else {
        return;
    };
    if position - i64::from(length) >= end {
        *flags |= 0x4;
        *flags &= !0x2;
        return;
    }
    let retained = u32::try_from(end - (position - i64::from(length))).unwrap();
    let mut result = cigar[..index].to_vec();
    if retained > 0 {
        result.push((kind, retained));
    }
    let mut soft_clip = length - retained;
    let mut hard_clip = None;
    for &(trailing_kind, trailing_length) in &cigar[index + 1..] {
        if trailing_kind == 5 {
            hard_clip = Some(hard_clip.unwrap_or(0) + trailing_length);
        } else {
            soft_clip += trailing_length;
        }
    }
    if soft_clip > 0 {
        result.push((4, soft_clip));
    }
    if let Some(length) = hard_clip {
        result.push((5, length));
    }
    *cigar = result;
}

fn collapse_cigar(cigar: &mut Vec<(u8, u32)>) {
    const MAX_LENGTH: u64 = (1u64 << 28) - 1;
    let mut written = 0;
    for read in 0..cigar.len() {
        let (kind, mut length) = cigar[read];
        if length == 0 {
            continue;
        }
        if written > 0 && cigar[written - 1].0 == kind {
            let total = u64::from(cigar[written - 1].1) + u64::from(length);
            cigar[written - 1].1 = total.min(MAX_LENGTH) as u32;
            length = u32::try_from(total.saturating_sub(MAX_LENGTH)).unwrap();
            if length == 0 {
                continue;
            }
        }
        cigar[written] = (kind, length);
        written += 1;
    }
    cigar.truncate(written);
}

fn append_canonical(output: &mut Vec<u8>, field: AuxField, bytes: &[u8]) {
    let value = &bytes[field.value_start..field.value_end];
    output.extend_from_slice(&field.tag);
    if matches!(field.type_code, b'c' | b'C' | b's' | b'S' | b'i' | b'I') {
        let integer = match field.type_code {
            b'c' => i64::from(value[0] as i8),
            b'C' => i64::from(value[0]),
            b's' => i64::from(i16::from_le_bytes(value[..2].try_into().unwrap())),
            b'S' => i64::from(u16::from_le_bytes(value[..2].try_into().unwrap())),
            b'i' => i64::from(i32::from_le_bytes(value[..4].try_into().unwrap())),
            b'I' => i64::from(u32::from_le_bytes(value[..4].try_into().unwrap())),
            _ => unreachable!(),
        };
        if integer >= 0 {
            if integer <= i64::from(u8::MAX) {
                output.push(b'C');
                output.push(integer as u8);
            } else if integer <= i64::from(u16::MAX) {
                output.push(b'S');
                output.extend_from_slice(&(integer as u16).to_le_bytes());
            } else {
                output.push(b'I');
                output.extend_from_slice(&(integer as u32).to_le_bytes());
            }
        } else if integer >= i64::from(i8::MIN) {
            output.push(b'c');
            output.push(integer as i8 as u8);
        } else if integer >= i64::from(i16::MIN) {
            output.push(b's');
            output.extend_from_slice(&(integer as i16).to_le_bytes());
        } else {
            output.push(b'i');
            output.extend_from_slice(&(integer as i32).to_le_bytes());
        }
    } else {
        output.push(field.type_code);
        output.extend_from_slice(value);
    }
}

fn value_length(bytes: &[u8], position: usize, type_code: u8) -> usize {
    match type_code {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'Z' | b'H' => {
            bytes[position..]
                .iter()
                .position(|&byte| byte == 0)
                .unwrap()
                + 1
        }
        b'B' => {
            let width = match bytes[position] {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => unreachable!(),
            };
            let count = u32::from_le_bytes(bytes[position + 1..position + 5].try_into().unwrap());
            5 + usize::try_from(count).unwrap() * width
        }
        _ => unreachable!(),
    }
}

fn canonical_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'=' | b'A' | b'C' | b'M' | b'G' | b'R' | b'S' | b'V' | b'T' | b'W' | b'Y' | b'H'
        | b'K' | b'D' | b'B' | b'N' => base.to_ascii_uppercase(),
        _ => b'N',
    }
}

fn crc(parts: &[&[u8]]) -> u32 {
    let mut hasher = Hasher::new();
    for part in parts {
        hasher.update(part);
    }
    hasher.finalize()
}

fn crc_with(initial: u32, bytes: &[u8]) -> u32 {
    let mut hasher = Hasher::new_with_initial(initial);
    hasher.update(bytes);
    hasher.finalize()
}

fn crc_repeated(initial: u32, byte: u8, length: usize) -> u32 {
    let block = [byte; 1024];
    let mut hasher = Hasher::new_with_initial(initial);
    for _ in 0..length / block.len() {
        hasher.update(&block);
    }
    hasher.update(&block[..length % block.len()]);
    hasher.finalize()
}
