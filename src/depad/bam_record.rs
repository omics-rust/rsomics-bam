use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};

use super::cigar::{CigarOp, SKIP, SOFT_CLIP, cigar_overflow, reference_span};

pub(super) fn replace_raw_cigar(
    record: &mut RawRecord,
    cigar: &[CigarOp],
    source_long_cigar: bool,
) -> Result<()> {
    let long_cigar = cigar.len() > usize::from(u16::MAX);
    let mut encoded = Vec::with_capacity(if long_cigar { 2 } else { cigar.len() });
    if long_cigar {
        encoded.push(pack_cigar(
            u32::try_from(record.sequence_len()).map_err(|_| {
                RsomicsError::InvalidInput("BAM sequence length exceeds u32".to_owned())
            })?,
            SOFT_CLIP,
        )?);
        encoded.push(pack_cigar(reference_span(cigar)?, SKIP)?);
    } else {
        for op in cigar {
            encoded.push(pack_cigar(op.len, op.kind)?);
        }
    }

    let source = record.as_bytes();
    let name_len = usize::from(source[8]);
    let old_cigar_count = record.cigar_ops().count();
    let cigar_start = 32usize.checked_add(name_len).ok_or_else(cigar_overflow)?;
    let old_cigar_end = old_cigar_count
        .checked_mul(4)
        .and_then(|length| cigar_start.checked_add(length))
        .ok_or_else(cigar_overflow)?;
    let suffix = source
        .get(old_cigar_end..)
        .ok_or_else(|| RsomicsError::InvalidInput("invalid BAM CIGAR boundary".to_owned()))?;
    let mut bytes = Vec::with_capacity(
        cigar_start
            .checked_add(encoded.len().saturating_mul(4))
            .and_then(|length| length.checked_add(suffix.len()))
            .ok_or_else(cigar_overflow)?,
    );
    bytes.extend_from_slice(&source[..cigar_start]);
    for operation in encoded {
        bytes.extend_from_slice(&operation.to_le_bytes());
    }
    bytes.extend_from_slice(suffix);
    let count = u16::try_from(if long_cigar { 2 } else { cigar.len() }).map_err(|_| {
        RsomicsError::InvalidInput("BAM CIGAR operation count exceeds u16".to_owned())
    })?;
    bytes[12..14].copy_from_slice(&count.to_le_bytes());

    let start = i32::from_le_bytes(bytes[4..8].try_into().expect("BAM position has four bytes"));
    let span = i64::from(reference_span(cigar)?.max(1));
    let end = i64::from(start)
        .checked_add(span)
        .ok_or_else(cigar_overflow)?;
    let bin = reg2bin(i64::from(start), end)?;
    bytes[10..12].copy_from_slice(&bin.to_le_bytes());

    let mut rebuilt = RawRecord::try_from(bytes)?;
    if source_long_cigar || long_cigar {
        rebuilt.remove_aux(*b"CG");
    }
    if long_cigar {
        let mut value = Vec::with_capacity(
            cigar
                .len()
                .checked_mul(4)
                .and_then(|length| length.checked_add(5))
                .ok_or_else(cigar_overflow)?,
        );
        value.push(b'I');
        value.extend_from_slice(
            &u32::try_from(cigar.len())
                .map_err(|_| RsomicsError::InvalidInput("long CIGAR count exceeds u32".to_owned()))?
                .to_le_bytes(),
        );
        for op in cigar {
            value.extend_from_slice(&pack_cigar(op.len, op.kind)?.to_le_bytes());
        }
        rebuilt.append_aux(*b"CG", b'B', &value)?;
    }
    *record = rebuilt;
    Ok(())
}

fn pack_cigar(length: u32, kind: u8) -> Result<u32> {
    if length == 0 || length >= 1 << 28 {
        return Err(RsomicsError::InvalidInput(format!(
            "CIGAR operation length {length} exceeds BAM limits"
        )));
    }
    Ok((length << 4) | u32::from(kind))
}

fn reg2bin(start: i64, end: i64) -> Result<u16> {
    if start < 0 || end <= start {
        return Err(RsomicsError::InvalidInput(
            "cannot calculate BAM bin for invalid coordinates".to_owned(),
        ));
    }
    let end = end - 1;
    let bin = if start >> 14 == end >> 14 {
        ((1 << 15) - 1) / 7 + (start >> 14)
    } else if start >> 17 == end >> 17 {
        ((1 << 12) - 1) / 7 + (start >> 17)
    } else if start >> 20 == end >> 20 {
        ((1 << 9) - 1) / 7 + (start >> 20)
    } else if start >> 23 == end >> 23 {
        ((1 << 6) - 1) / 7 + (start >> 23)
    } else if start >> 26 == end >> 26 {
        1 + (start >> 26)
    } else {
        0
    };
    u16::try_from(bin).map_err(|_| RsomicsError::InvalidInput("BAM bin exceeds u16".to_owned()))
}
