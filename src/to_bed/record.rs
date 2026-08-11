use std::fmt::Write as _;

use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};

use crate::raw_aux::{self, Integer};

pub(super) const UNMAPPED: u16 = 0x4;
pub(super) const PAIRED: u16 = 0x1;
pub(super) const REVERSE: u16 = 0x10;
pub(super) const READ1: u16 = 0x40;
pub(super) const READ2: u16 = 0x80;

pub(super) struct Mapped<'a> {
    pub reference: &'a str,
    pub start: u64,
    pub name: &'a [u8],
    pub flags: u16,
}

pub(super) fn project<'a>(
    record: &'a RawRecord,
    references: &'a [String],
) -> Result<Option<Mapped<'a>>> {
    let flags = record.flags();
    if flags & UNMAPPED != 0 {
        return Ok(None);
    }
    let reference_id = usize::try_from(record.reference_sequence_id()).map_err(|_| {
        RsomicsError::InvalidInput(format!(
            "mapped record {} has no reference sequence",
            String::from_utf8_lossy(record.name())
        ))
    })?;
    let reference = references.get(reference_id).ok_or_else(|| {
        RsomicsError::InvalidInput(format!(
            "record {} has reference ID {} outside the header",
            String::from_utf8_lossy(record.name()),
            reference_id
        ))
    })?;
    let start = u64::try_from(record.alignment_start()).map_err(|_| {
        RsomicsError::InvalidInput(format!(
            "mapped record {} has no alignment start",
            String::from_utf8_lossy(record.name())
        ))
    })?;
    Ok(Some(Mapped {
        reference,
        start,
        name: record.name(),
        flags,
    }))
}

pub(super) fn reference_end(cigar: &[(u8, u32)], start: u64, name: &[u8]) -> Result<u64> {
    let span = cigar.iter().try_fold(0u64, |span, &(kind, length)| {
        if matches!(kind, 0 | 2 | 3 | 7 | 8) {
            span.checked_add(u64::from(length))
        } else {
            Some(span)
        }
    });
    start
        .checked_add(span.ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "record {} reference span overflows",
                String::from_utf8_lossy(name)
            ))
        })?)
        .ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "record {} alignment end overflows",
                String::from_utf8_lossy(name)
            ))
        })
}

pub(super) fn blocks(
    cigar: &[(u8, u32)],
    start: u64,
    split_deletions: bool,
) -> Result<Vec<(u64, u64)>> {
    let mut blocks = Vec::new();
    let mut block_start = start;
    let mut block_length = 0u64;
    for &(kind, length) in cigar {
        let length = u64::from(length);
        match kind {
            0 | 7 | 8 => {
                block_length = block_length.checked_add(length).ok_or_else(|| {
                    RsomicsError::InvalidInput("alignment block end overflows".to_owned())
                })?;
            }
            2 if !split_deletions => {
                block_length = block_length.checked_add(length).ok_or_else(|| {
                    RsomicsError::InvalidInput("alignment block end overflows".to_owned())
                })?;
            }
            2 | 3 => {
                let block_end = block_start.checked_add(block_length).ok_or_else(|| {
                    RsomicsError::InvalidInput("alignment block end overflows".to_owned())
                })?;
                blocks.push((block_start, block_end));
                block_start = block_end.checked_add(length).ok_or_else(|| {
                    RsomicsError::InvalidInput("alignment block start overflows".to_owned())
                })?;
                block_length = 0;
            }
            _ => {}
        }
    }
    let block_end = block_start
        .checked_add(block_length)
        .ok_or_else(|| RsomicsError::InvalidInput("alignment block end overflows".to_owned()))?;
    blocks.push((block_start, block_end));
    Ok(blocks)
}

pub(super) fn score(record: &RawRecord, tag: Option<[u8; 2]>) -> Result<i64> {
    let Some(tag) = tag else {
        return Ok(i64::from(record.mapping_quality()));
    };
    match raw_aux::integer(record, tag) {
        Integer::Value(value) => Ok(value),
        Integer::Missing => Err(RsomicsError::InvalidInput(format!(
            "record {} is missing numeric tag {}",
            String::from_utf8_lossy(record.name()),
            String::from_utf8_lossy(&tag)
        ))),
        Integer::Invalid => Err(RsomicsError::InvalidInput(format!(
            "record {} tag {} is not a valid integer",
            String::from_utf8_lossy(record.name()),
            String::from_utf8_lossy(&tag)
        ))),
    }
}

pub(super) fn cigar_text(cigar_ops: &[(u8, u32)], name: &[u8]) -> Result<String> {
    const KINDS: &[u8] = b"MIDNSHP=X";
    let mut cigar = String::new();
    for &(kind, length) in cigar_ops {
        let kind = KINDS.get(usize::from(kind)).copied().ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "record {} contains an invalid CIGAR operation",
                String::from_utf8_lossy(name)
            ))
        })?;
        write!(cigar, "{length}{}", char::from(kind)).unwrap();
    }
    if cigar.is_empty() {
        cigar.push('*');
    }
    Ok(cigar)
}
