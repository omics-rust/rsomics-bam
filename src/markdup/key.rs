use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};

use super::Mode;

const REVERSE: u16 = 0x10;
const MATE_REVERSE: u16 = 0x20;
const READ1: u16 = 0x40;

const FORWARD_FORWARD: u8 = 2;
const REVERSE_REVERSE: u8 = 3;
const FORWARD_REVERSE: u8 = 5;
const REVERSE_FORWARD: u8 = 7;
const LEFT: u8 = 11;
const RIGHT: u8 = 13;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) struct SingleKey {
    reference: i64,
    coordinate: i64,
    orientation: u8,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) struct PairKey {
    reference: i64,
    coordinate: i64,
    mate_reference: i64,
    mate_coordinate: i64,
    side: u8,
    orientation: u8,
}

pub(crate) fn single(record: &RawRecord) -> Result<SingleKey> {
    let (start, end) = record_coordinates(record)?;
    let (coordinate, orientation) = if flag(record, REVERSE) {
        (end, REVERSE_REVERSE)
    } else {
        (start, FORWARD_FORWARD)
    };
    Ok(SingleKey {
        reference: i64::from(record.reference_sequence_id()) + 1,
        coordinate,
        orientation,
    })
}

pub(crate) fn pair(record: &RawRecord, mode: Mode) -> Result<PairKey> {
    let (start, end) = record_coordinates(record)?;
    let (mate_start, mate_end) = mate_coordinates(record)?;
    let reference = i64::from(record.reference_sequence_id()) + 1;
    let mate_reference = i64::from(record.mate_reference_sequence_id()) + 1;
    let reverse = flag(record, REVERSE);
    let mate_reverse = flag(record, MATE_REVERSE);

    let left = match mode {
        Mode::Template => {
            if reference != mate_reference {
                reference < mate_reference
            } else if reverse == mate_reverse {
                if reverse {
                    end <= mate_end
                } else {
                    start <= mate_start
                }
            } else if reverse {
                end <= mate_start
            } else {
                start <= mate_end
            }
        }
        Mode::Sequence => sequence_left(record, start, end, mate_start, mate_end),
    };

    let (coordinate, mate_coordinate, orientation) = match mode {
        Mode::Template => template_fields(record, left, start, end, mate_start, mate_end),
        Mode::Sequence => (
            if reverse { end } else { start },
            if mate_reverse { mate_end } else { mate_start },
            sequence_orientation(left, reverse, mate_reverse),
        ),
    };

    Ok(PairKey {
        reference,
        coordinate,
        mate_reference,
        mate_coordinate,
        side: if left { LEFT } else { RIGHT },
        orientation,
    })
}

fn sequence_left(record: &RawRecord, start: i64, end: i64, mate_start: i64, mate_end: i64) -> bool {
    let reference = record.reference_sequence_id();
    let mate_reference = record.mate_reference_sequence_id();
    let reverse = flag(record, REVERSE);
    let mate_reverse = flag(record, MATE_REVERSE);
    let ordering = if reference != mate_reference {
        reference.cmp(&mate_reference)
    } else if reverse == mate_reverse {
        if reverse {
            end.cmp(&mate_end)
        } else {
            start.cmp(&mate_start)
        }
    } else if reverse {
        end.cmp(&mate_start)
    } else {
        start.cmp(&mate_end)
    };
    match ordering {
        std::cmp::Ordering::Less => true,
        std::cmp::Ordering::Greater => false,
        std::cmp::Ordering::Equal => {
            if record.alignment_start() == record.mate_alignment_start() {
                flag(record, READ1)
            } else {
                record.alignment_start() < record.mate_alignment_start()
            }
        }
    }
}

fn template_fields(
    record: &RawRecord,
    left: bool,
    start: i64,
    end: i64,
    mate_start: i64,
    mate_end: i64,
) -> (i64, i64, u8) {
    let reverse = flag(record, REVERSE);
    let mate_reverse = flag(record, MATE_REVERSE);
    let read1 = flag(record, READ1);
    if left {
        if reverse == mate_reverse {
            let orientation = if reverse != read1 {
                FORWARD_FORWARD
            } else {
                REVERSE_REVERSE
            };
            (start, mate_end, orientation)
        } else if reverse {
            (end, mate_start, REVERSE_FORWARD)
        } else {
            (start, mate_end, FORWARD_REVERSE)
        }
    } else if reverse == mate_reverse {
        let orientation = if reverse == read1 {
            FORWARD_FORWARD
        } else {
            REVERSE_REVERSE
        };
        (end, mate_start, orientation)
    } else if reverse {
        (end, mate_start, FORWARD_REVERSE)
    } else {
        (start, mate_end, REVERSE_FORWARD)
    }
}

fn sequence_orientation(left: bool, reverse: bool, mate_reverse: bool) -> u8 {
    match (left, reverse, mate_reverse) {
        (true, false, false) => FORWARD_FORWARD,
        (true, true, true) => REVERSE_REVERSE,
        (true, false, true) => FORWARD_REVERSE,
        (true, true, false) => REVERSE_FORWARD,
        (false, false, false) => REVERSE_REVERSE,
        (false, true, true) => FORWARD_FORWARD,
        (false, false, true) => REVERSE_FORWARD,
        (false, true, false) => FORWARD_REVERSE,
    }
}

pub(crate) fn coordinate(key: SingleKey) -> i64 {
    key.coordinate
}

fn record_coordinates(record: &RawRecord) -> Result<(i64, i64)> {
    if record.aux_type(*b"CG") == Some(b'B') {
        return coordinates(
            i64::from(record.alignment_start()),
            record.decoded_cigar()?.into_iter(),
        );
    }
    coordinates(i64::from(record.alignment_start()), record.cigar_ops())
}

fn coordinates(position: i64, cigar: impl Iterator<Item = (u8, u32)>) -> Result<(i64, i64)> {
    let mut leading = 0i64;
    let mut trailing = 0i64;
    let mut span = 0i64;
    let mut before_alignment = true;
    for (operation, length) in cigar {
        let length = i64::from(length);
        if matches!(operation, 4 | 5) {
            if before_alignment {
                leading = leading
                    .checked_add(length)
                    .ok_or_else(coordinate_overflow)?;
            }
            trailing = trailing
                .checked_add(length)
                .ok_or_else(coordinate_overflow)?;
        } else {
            before_alignment = false;
            trailing = 0;
        }
        if matches!(operation, 0 | 2 | 3 | 7 | 8) {
            span = span.checked_add(length).ok_or_else(coordinate_overflow)?;
        }
    }
    let start = position
        .checked_sub(leading)
        .and_then(|value| value.checked_add(1))
        .ok_or_else(coordinate_overflow)?;
    let end = position
        .checked_add(span.max(1))
        .and_then(|value| value.checked_add(trailing))
        .ok_or_else(coordinate_overflow)?;
    Ok((start, end))
}

fn mate_coordinates(record: &RawRecord) -> Result<(i64, i64)> {
    if record.aux_type(*b"MC") != Some(b'Z') {
        return Err(RsomicsError::InvalidInput(
            "paired markdup record requires an MC:Z tag from fixmate -m".to_owned(),
        ));
    }
    let value = record.aux_value(*b"MC").ok_or_else(|| {
        RsomicsError::InvalidInput(
            "paired markdup record requires an MC:Z tag from fixmate -m".to_owned(),
        )
    })?;
    let cigar = value
        .strip_suffix(&[0])
        .ok_or_else(|| RsomicsError::InvalidInput("MC tag is not NUL terminated".to_owned()))?;
    let position = i64::from(record.mate_alignment_start());
    if cigar == b"*" {
        return Ok((position + 1, position));
    }

    let mut index = 0;
    let mut leading = 0i64;
    let mut span = 0i64;
    let mut before_alignment = true;
    while index < cigar.len() {
        let number_start = index;
        while cigar.get(index).is_some_and(u8::is_ascii_digit) {
            index += 1;
        }
        if number_start == index || index == cigar.len() {
            return Err(RsomicsError::InvalidInput(
                "MC tag contains an invalid CIGAR".to_owned(),
            ));
        }
        let length = std::str::from_utf8(&cigar[number_start..index])
            .ok()
            .and_then(|value| value.parse::<i64>().ok())
            .filter(|&value| value > 0)
            .ok_or_else(|| {
                RsomicsError::InvalidInput("MC tag contains an invalid CIGAR".to_owned())
            })?;
        let operation = cigar[index];
        index += 1;
        match operation {
            b'S' | b'H' if before_alignment => {
                leading = leading
                    .checked_add(length)
                    .ok_or_else(coordinate_overflow)?;
            }
            b'M' | b'D' | b'N' | b'=' | b'X' => {
                before_alignment = false;
                span = span.checked_add(length).ok_or_else(coordinate_overflow)?;
            }
            b'S' | b'H' => {
                span = span.checked_add(length).ok_or_else(coordinate_overflow)?;
            }
            b'I' | b'P' => {}
            _ => {
                return Err(RsomicsError::InvalidInput(
                    "MC tag contains an invalid CIGAR".to_owned(),
                ));
            }
        }
    }
    let start = position
        .checked_sub(leading)
        .and_then(|value| value.checked_add(1))
        .ok_or_else(coordinate_overflow)?;
    let end = position.checked_add(span).ok_or_else(coordinate_overflow)?;
    Ok((start, end))
}

fn flag(record: &RawRecord, bits: u16) -> bool {
    record.flags() & bits != 0
}

fn coordinate_overflow() -> RsomicsError {
    RsomicsError::InvalidInput("markdup coordinate arithmetic overflowed".to_owned())
}
