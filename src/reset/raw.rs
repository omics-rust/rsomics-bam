use rsomics_bamio::raw::RecordRef;
use rsomics_common::{Result, RsomicsError};

use super::{FLAG_SECONDARY, FLAG_SUPPLEMENTARY, TagFilter, transformed_flags};

pub(super) fn reset(
    record: &RecordRef<'_>,
    tags: TagFilter<'_>,
    keep_duplicate: bool,
    output: &mut Vec<u8>,
) -> Result<bool> {
    let flags = record.flags();
    if flags & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
        return Ok(false);
    }

    let source = record.payload();
    let name_len = usize::from(source[8]);
    let cigar_len = usize::from(u16::from_le_bytes(source[12..14].try_into().unwrap()))
        .checked_mul(4)
        .ok_or_else(|| invalid_record(record, "CIGAR length overflows"))?;
    let sequence_len = usize::try_from(u32::from_le_bytes(source[16..20].try_into().unwrap()))
        .map_err(|error| invalid_record(record, error))?;
    let name_end = 32usize
        .checked_add(name_len)
        .ok_or_else(|| invalid_record(record, "read name length overflows"))?;
    let sequence_start = name_end
        .checked_add(cigar_len)
        .ok_or_else(|| invalid_record(record, "record layout overflows"))?;
    let quality_start = sequence_start
        .checked_add(sequence_len.div_ceil(2))
        .ok_or_else(|| invalid_record(record, "record layout overflows"))?;
    let aux_start = quality_start
        .checked_add(sequence_len)
        .ok_or_else(|| invalid_record(record, "record layout overflows"))?;

    output.clear();
    output.reserve(source.len().saturating_sub(cigar_len));
    output.extend_from_slice(&source[..32]);
    output[..4].copy_from_slice(&(-1i32).to_le_bytes());
    output[4..8].copy_from_slice(&(-1i32).to_le_bytes());
    output[9] = 0;
    output[10..12].copy_from_slice(&4680u16.to_le_bytes());
    output[12..14].fill(0);
    let (flags, reverse) = transformed_flags(flags, keep_duplicate);
    output[14..16].copy_from_slice(&flags.to_le_bytes());
    output[20..24].copy_from_slice(&(-1i32).to_le_bytes());
    output[24..28].copy_from_slice(&(-1i32).to_le_bytes());
    output[28..32].fill(0);
    output.extend_from_slice(&source[32..name_end]);

    let packed_len = sequence_len.div_ceil(2);
    if reverse {
        reverse_complement(
            &source[sequence_start..sequence_start + packed_len],
            sequence_len,
            output,
        );
        output.extend(source[quality_start..aux_start].iter().rev());
    } else {
        output.extend_from_slice(&source[sequence_start..sequence_start + packed_len]);
        output.extend_from_slice(&source[quality_start..aux_start]);
    }

    let mut position = aux_start;
    while position < source.len() {
        let end = aux_field_end(source, position)
            .ok_or_else(|| invalid_record(record, "auxiliary field is malformed"))?;
        let tag = [source[position], source[position + 1]];
        if !tags.remove(tag) {
            output.extend_from_slice(&source[position..end]);
        }
        position = end;
    }
    Ok(true)
}

fn reverse_complement(source: &[u8], sequence_len: usize, output: &mut Vec<u8>) {
    if sequence_len.is_multiple_of(2) {
        output.extend(
            source
                .iter()
                .rev()
                .map(|byte| (complement_nibble(byte & 0x0f) << 4) | complement_nibble(byte >> 4)),
        );
        return;
    }

    let start = output.len();
    output.resize(start + source.len(), 0);
    for target_index in 0..sequence_len {
        let source_index = sequence_len - target_index - 1;
        let source_byte = source[source_index / 2];
        let base = if source_index.is_multiple_of(2) {
            source_byte >> 4
        } else {
            source_byte & 0x0f
        };
        let base = complement_nibble(base);
        if target_index.is_multiple_of(2) {
            output[start + target_index / 2] = base << 4;
        } else {
            output[start + target_index / 2] |= base;
        }
    }
}

fn aux_field_end(bytes: &[u8], start: usize) -> Option<usize> {
    let value_start = start.checked_add(3)?;
    let type_code = *bytes.get(start + 2)?;
    let value_len = match type_code {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'Z' | b'H' => {
            bytes
                .get(value_start..)?
                .iter()
                .position(|byte| *byte == 0)?
                + 1
        }
        b'B' => {
            let element_size = match *bytes.get(value_start)? {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return None,
            };
            let count = usize::try_from(u32::from_le_bytes(
                bytes
                    .get(value_start + 1..value_start + 5)?
                    .try_into()
                    .ok()?,
            ))
            .ok()?;
            count.checked_mul(element_size)?.checked_add(5)?
        }
        _ => return None,
    };
    value_start
        .checked_add(value_len)
        .filter(|end| *end <= bytes.len())
}

fn invalid_record(record: &RecordRef<'_>, reason: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "read {}: {reason}",
        String::from_utf8_lossy(record.name())
    ))
}

fn complement_nibble(base: u8) -> u8 {
    const COMPLEMENT: [u8; 16] = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];
    COMPLEMENT[usize::from(base)]
}
