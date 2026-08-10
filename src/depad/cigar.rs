use noodles::core::Position;
use noodles::sam::{
    self,
    alignment::{
        RecordBuf,
        record::cigar::{Op, op::Kind},
        record_buf::Cigar,
    },
};
use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct CigarOp {
    pub(super) kind: u8,
    pub(super) len: u32,
}

impl CigarOp {
    pub(super) fn new(kind: u8, len: u32) -> Result<Self> {
        if len == 0 {
            return Err(RsomicsError::InvalidInput(
                "CIGAR operation length cannot be zero".to_owned(),
            ));
        }
        Ok(Self { kind, len })
    }
}

pub(super) const MATCH: u8 = 0;
pub(super) const INSERTION: u8 = 1;
pub(super) const DELETION: u8 = 2;
pub(super) const SKIP: u8 = 3;
pub(super) const SOFT_CLIP: u8 = 4;
pub(super) const HARD_CLIP: u8 = 5;
pub(super) const PAD: u8 = 6;
pub(super) const EQUAL: u8 = 7;
pub(super) const DIFFERENCE: u8 = 8;

pub(super) fn typed_cigar(cigar: &Cigar) -> Result<Vec<CigarOp>> {
    cigar
        .as_ref()
        .iter()
        .map(|op| {
            let kind = match op.kind() {
                Kind::Match => MATCH,
                Kind::Insertion => INSERTION,
                Kind::Deletion => DELETION,
                Kind::Skip => SKIP,
                Kind::SoftClip => SOFT_CLIP,
                Kind::HardClip => HARD_CLIP,
                Kind::Pad => PAD,
                Kind::SequenceMatch => EQUAL,
                Kind::SequenceMismatch => DIFFERENCE,
            };
            CigarOp::new(
                kind,
                u32::try_from(op.len()).map_err(|_| {
                    RsomicsError::InvalidInput("CIGAR operation exceeds u32".to_owned())
                })?,
            )
        })
        .collect()
}

pub(super) fn typed_projected_cigar(cigar: &[CigarOp]) -> Cigar {
    cigar
        .iter()
        .map(|op| {
            let kind = match op.kind {
                MATCH => Kind::Match,
                INSERTION => Kind::Insertion,
                DELETION => Kind::Deletion,
                SOFT_CLIP => Kind::SoftClip,
                HARD_CLIP => Kind::HardClip,
                PAD => Kind::Pad,
                _ => unreachable!("projected CIGAR kind"),
            };
            Op::new(
                kind,
                usize::try_from(op.len).expect("u32 CIGAR length fits supported targets"),
            )
        })
        .collect()
}

pub(super) fn decode_padded(
    cigar: &[CigarOp],
    query: &[u8],
    name: &str,
    output: &mut Vec<u8>,
) -> Result<()> {
    output.clear();
    let mut query_position = 0usize;
    for op in cigar {
        let len = op.len as usize;
        match op.kind {
            MATCH | EQUAL | DIFFERENCE => {
                let end = query_position.checked_add(len).ok_or_else(cigar_overflow)?;
                let bases = query.get(query_position..end).ok_or_else(|| {
                    RsomicsError::InvalidInput(format!(
                        "record {name} CIGAR consumes beyond its sequence"
                    ))
                })?;
                output.extend(bases.iter().map(|&base| query_code(base)));
                query_position = end;
            }
            DELETION | SKIP => {
                output.resize(output.len().checked_add(len).ok_or_else(cigar_overflow)?, 0)
            }
            SOFT_CLIP => {
                query_position = query_position.checked_add(len).ok_or_else(cigar_overflow)?;
                if query_position > query.len() {
                    return Err(RsomicsError::InvalidInput(format!(
                        "record {name} CIGAR consumes beyond its sequence"
                    )));
                }
            }
            HARD_CLIP => {}
            INSERTION | PAD => {
                return Err(RsomicsError::InvalidInput(format!(
                    "record {name} has unsupported input CIGAR operation {}",
                    cigar_letter(op.kind)
                )));
            }
            _ => {
                return Err(RsomicsError::InvalidInput(format!(
                    "record {name} has unknown CIGAR operation {}",
                    op.kind
                )));
            }
        }
    }
    if query_position != query.len() {
        return Err(RsomicsError::InvalidInput(format!(
            "record {name} CIGAR consumes {query_position} of {} query bases",
            query.len()
        )));
    }
    Ok(())
}

pub(super) fn decode_padded_raw(
    cigar: &[CigarOp],
    record: &RawRecord,
    name: &str,
    output: &mut Vec<u8>,
) -> Result<()> {
    output.clear();
    let mut query_position = 0usize;
    for op in cigar {
        let len = op.len as usize;
        match op.kind {
            MATCH | EQUAL | DIFFERENCE => {
                let end = query_position.checked_add(len).ok_or_else(cigar_overflow)?;
                if end > record.sequence_len() {
                    return Err(RsomicsError::InvalidInput(format!(
                        "record {name} CIGAR consumes beyond its sequence"
                    )));
                }
                output.extend((query_position..end).map(|index| record.seq_nibble(index)));
                query_position = end;
            }
            DELETION | SKIP => {
                output.resize(output.len().checked_add(len).ok_or_else(cigar_overflow)?, 0)
            }
            SOFT_CLIP => {
                query_position = query_position.checked_add(len).ok_or_else(cigar_overflow)?;
                if query_position > record.sequence_len() {
                    return Err(RsomicsError::InvalidInput(format!(
                        "record {name} CIGAR consumes beyond its sequence"
                    )));
                }
            }
            HARD_CLIP => {}
            INSERTION | PAD => {
                return Err(RsomicsError::InvalidInput(format!(
                    "record {name} has unsupported input CIGAR operation {}",
                    cigar_letter(op.kind)
                )));
            }
            _ => {
                return Err(RsomicsError::InvalidInput(format!(
                    "record {name} has unknown CIGAR operation {}",
                    op.kind
                )));
            }
        }
    }
    if query_position != record.sequence_len() {
        return Err(RsomicsError::InvalidInput(format!(
            "record {name} CIGAR consumes {query_position} of {} query bases",
            record.sequence_len()
        )));
    }
    Ok(())
}

pub(super) fn reference_span(cigar: &[CigarOp]) -> Result<u32> {
    cigar
        .iter()
        .filter(|op| matches!(op.kind, MATCH | DELETION | SKIP | EQUAL | DIFFERENCE))
        .try_fold(0u32, |span, op| {
            span.checked_add(op.len).ok_or_else(cigar_overflow)
        })
}

pub(super) fn project_cigar(
    original: &[CigarOp],
    query: &[u8],
    reference: &[u8],
    start: usize,
    name: &str,
    output: &mut Vec<CigarOp>,
) -> Result<()> {
    if query.is_empty() {
        return Err(RsomicsError::InvalidInput(format!(
            "mapped record {name} has no aligned columns"
        )));
    }
    let end = start.checked_add(query.len()).ok_or_else(cigar_overflow)?;
    output.clear();
    append_leading_clips(original, output);

    let reference_segment = reference.get(start..end).ok_or_else(|| {
        RsomicsError::InvalidInput(format!(
            "record {name} alignment exceeds the padded reference"
        ))
    })?;
    let mut kinds = query
        .iter()
        .zip(reference_segment)
        .map(|(&query, &reference)| match (query != 0, reference != 0) {
            (true, true) => MATCH,
            (true, false) => INSERTION,
            (false, true) => DELETION,
            (false, false) => PAD,
        });
    let first = kinds.next().expect("query is nonempty");
    let leading_pads = if matches!(first, INSERTION | PAD) {
        count_preceding_pads(reference, start)?
    } else {
        0
    };
    if first == INSERTION && leading_pads > 0 {
        push_cigar(output, PAD, leading_pads)?;
    }
    let mut current = first;
    let mut length = if first == PAD { leading_pads + 1 } else { 1 };
    for kind in kinds {
        if kind == current {
            length = length.checked_add(1).ok_or_else(cigar_overflow)?;
        } else {
            push_cigar(output, current, length)?;
            current = kind;
            length = 1;
        }
    }
    push_cigar(output, current, length)?;
    remove_redundant_pads(output)?;
    append_trailing_clips(original, output);
    Ok(())
}

fn count_preceding_pads(reference: &[u8], start: usize) -> Result<u32> {
    let mut count = 0usize;
    while count + 1 < start && reference.get(start - count - 1) == Some(&0) {
        count += 1;
    }
    u32::try_from(count)
        .map_err(|_| RsomicsError::InvalidInput("leading pad count exceeds u32".to_owned()))
}

fn append_leading_clips(original: &[CigarOp], output: &mut Vec<CigarOp>) {
    let mut index = 0;
    if original.first().is_some_and(|op| op.kind == HARD_CLIP) {
        output.push(original[0]);
        index = 1;
    }
    if original.get(index).is_some_and(|op| op.kind == SOFT_CLIP) {
        output.push(original[index]);
    }
}

fn append_trailing_clips(original: &[CigarOp], output: &mut Vec<CigarOp>) {
    let mut end = original.len();
    let hard = original.last().filter(|op| op.kind == HARD_CLIP).copied();
    if hard.is_some() {
        end -= 1;
    }
    if let Some(soft) = original
        .get(end.wrapping_sub(1))
        .filter(|op| op.kind == SOFT_CLIP)
    {
        output.push(*soft);
    }
    if let Some(hard) = hard {
        output.push(hard);
    }
}

fn remove_redundant_pads(cigar: &mut Vec<CigarOp>) -> Result<()> {
    let mut index = 2;
    while index < cigar.len() {
        let left = cigar[index - 2];
        let middle = cigar[index - 1];
        let right = cigar[index];
        if middle.kind == PAD
            && matches!(left.kind, MATCH | DELETION)
            && matches!(right.kind, MATCH | DELETION)
        {
            cigar.remove(index - 1);
            if left.kind == right.kind {
                let len = left.len.checked_add(right.len).ok_or_else(cigar_overflow)?;
                cigar[index - 2] = CigarOp::new(left.kind, len)?;
                cigar.remove(index - 1);
            }
            index = index.saturating_sub(1).max(2);
        } else {
            index += 1;
        }
    }
    Ok(())
}

fn push_cigar(output: &mut Vec<CigarOp>, kind: u8, len: u32) -> Result<()> {
    if let Some(last) = output.last_mut().filter(|op| op.kind == kind) {
        last.len = last.len.checked_add(len).ok_or_else(cigar_overflow)?;
    } else {
        output.push(CigarOp::new(kind, len)?);
    }
    Ok(())
}

pub(super) fn build_position_map(sequence: &[u8], output: &mut Vec<usize>) {
    output.clear();
    output.reserve(sequence.len());
    let mut position = 0;
    for &base in sequence {
        output.push(position);
        position += usize::from(base != 0);
    }
}

pub(super) fn reference_name(header: &sam::Header, reference_id: usize) -> Result<&[u8]> {
    header
        .reference_sequences()
        .get_index(reference_id)
        .map(|(name, _)| name.as_ref())
        .ok_or_else(|| unknown_reference(reference_id))
}

pub(super) fn unknown_reference(reference_id: usize) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "record references unknown sequence ID {reference_id}"
    ))
}

pub(super) fn position(value: usize, name: &str) -> Result<Position> {
    value
        .checked_add(1)
        .and_then(Position::new)
        .ok_or_else(|| RsomicsError::InvalidInput(format!("record {name} position overflows")))
}

pub(super) fn record_name(record: &RecordBuf) -> String {
    record
        .name()
        .map(|name| String::from_utf8_lossy(name.as_ref()).into_owned())
        .unwrap_or_else(|| "*".to_owned())
}

fn query_code(base: u8) -> u8 {
    nucleotide_code(base).unwrap_or(15)
}

pub(super) fn reference_code(base: u8) -> Option<u8> {
    if matches!(base, b'*' | b'-') {
        Some(0)
    } else {
        nucleotide_code(base).filter(|&code| code != 0)
    }
}

fn nucleotide_code(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'=' => Some(0),
        b'A' => Some(1),
        b'C' => Some(2),
        b'M' => Some(3),
        b'G' => Some(4),
        b'R' => Some(5),
        b'S' => Some(6),
        b'V' => Some(7),
        b'T' | b'U' => Some(8),
        b'W' => Some(9),
        b'Y' => Some(10),
        b'H' => Some(11),
        b'K' => Some(12),
        b'D' => Some(13),
        b'B' => Some(14),
        b'N' => Some(15),
        _ => None,
    }
}

fn cigar_letter(kind: u8) -> char {
    b"MIDNSHP=X"
        .get(kind as usize)
        .copied()
        .map(char::from)
        .unwrap_or('?')
}

pub(super) fn cigar_overflow() -> RsomicsError {
    RsomicsError::InvalidInput("CIGAR length overflows".to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn position_map_counts_only_reference_bases() {
        let mut map = Vec::new();
        build_position_map(&[1, 2, 0, 8], &mut map);
        assert_eq!(map, [0, 1, 2, 2]);
    }

    #[test]
    fn redundant_internal_pads_are_removed() {
        let mut cigar = vec![
            CigarOp::new(MATCH, 5).unwrap(),
            CigarOp::new(PAD, 2).unwrap(),
            CigarOp::new(MATCH, 10).unwrap(),
        ];
        remove_redundant_pads(&mut cigar).unwrap();
        assert_eq!(cigar, [CigarOp::new(MATCH, 15).unwrap()]);
    }
}
