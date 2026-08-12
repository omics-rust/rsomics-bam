use rsomics_common::{Result, RsomicsError};

use super::builder::EvidenceRecord;

const NIBBLE_BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

pub(super) fn apply(
    record: &impl EvidenceRecord,
    cigar: &[(u8, u32)],
    md: &[u8],
    alignment_start: usize,
    reference: &mut [u8],
    reference_start: usize,
    reference_length: usize,
) -> Result<()> {
    let mut columns = Columns::new(cigar, record.sequence_len());
    let mut cursor = 0;
    while cursor < md.len() {
        match md[cursor] {
            digit if digit.is_ascii_digit() => {
                let mut count = 0usize;
                while cursor < md.len() && md[cursor].is_ascii_digit() {
                    count = count
                        .checked_mul(10)
                        .and_then(|value| value.checked_add(usize::from(md[cursor] - b'0')))
                        .ok_or_else(incompatible_md)?;
                    cursor += 1;
                }
                for _ in 0..count {
                    let column = columns.next_aligned(MdOperation::Match)?;
                    let base = query_base(record, column.query)?;
                    set_base(
                        reference,
                        reference_start,
                        reference_length,
                        alignment_start,
                        column.reference,
                        base,
                    )?;
                }
            }
            b'^' => {
                cursor += 1;
                let start = cursor;
                while cursor < md.len() && md[cursor].is_ascii_alphabetic() {
                    let column = columns.next_aligned(MdOperation::Deletion)?;
                    set_base(
                        reference,
                        reference_start,
                        reference_length,
                        alignment_start,
                        column.reference,
                        md[cursor].to_ascii_uppercase(),
                    )?;
                    cursor += 1;
                }
                if cursor == start {
                    return Err(incompatible_md());
                }
            }
            base if base.is_ascii_alphabetic() => {
                let column = columns.next_aligned(MdOperation::Mismatch)?;
                set_base(
                    reference,
                    reference_start,
                    reference_length,
                    alignment_start,
                    column.reference,
                    base.to_ascii_uppercase(),
                )?;
                cursor += 1;
            }
            _ => return Err(incompatible_md()),
        }
    }
    columns.finish()
}

#[derive(Clone, Copy)]
struct Column {
    query: usize,
    reference: usize,
    operation: u8,
}

enum MdOperation {
    Match,
    Mismatch,
    Deletion,
}

struct Columns<'a> {
    cigar: &'a [(u8, u32)],
    cigar_index: usize,
    operation: u8,
    remaining: usize,
    query: usize,
    reference: usize,
    sequence_len: usize,
}

impl<'a> Columns<'a> {
    fn new(cigar: &'a [(u8, u32)], sequence_len: usize) -> Self {
        Self {
            cigar,
            cigar_index: 0,
            operation: 0,
            remaining: 0,
            query: 0,
            reference: 0,
            sequence_len,
        }
    }

    fn next_aligned(&mut self, expected: MdOperation) -> Result<Column> {
        let column = self.next_column()?.ok_or_else(incompatible_md)?;
        let compatible = match expected {
            MdOperation::Match => matches!(column.operation, 0 | 7),
            MdOperation::Mismatch => matches!(column.operation, 0 | 8),
            MdOperation::Deletion => column.operation == 2,
        };
        if !compatible {
            return Err(incompatible_md());
        }
        Ok(column)
    }

    fn next_column(&mut self) -> Result<Option<Column>> {
        loop {
            if self.remaining == 0 {
                let Some(&(operation, length)) = self.cigar.get(self.cigar_index) else {
                    return Ok(None);
                };
                self.cigar_index += 1;
                self.operation = operation;
                self.remaining = usize::try_from(length).unwrap();
                if self.remaining == 0 {
                    return Err(incompatible_md());
                }
            }
            match self.operation {
                0 | 7 | 8 => {
                    let column = Column {
                        query: self.query,
                        reference: self.reference,
                        operation: self.operation,
                    };
                    self.query = self.query.checked_add(1).ok_or_else(incompatible_md)?;
                    self.reference = self.reference.checked_add(1).ok_or_else(incompatible_md)?;
                    self.remaining -= 1;
                    return Ok(Some(column));
                }
                2 => {
                    let column = Column {
                        query: self.query,
                        reference: self.reference,
                        operation: self.operation,
                    };
                    self.reference = self.reference.checked_add(1).ok_or_else(incompatible_md)?;
                    self.remaining -= 1;
                    return Ok(Some(column));
                }
                1 | 4 => {
                    self.query = self
                        .query
                        .checked_add(self.remaining)
                        .ok_or_else(incompatible_md)?;
                    self.remaining = 0;
                }
                3 => {
                    self.reference = self
                        .reference
                        .checked_add(self.remaining)
                        .ok_or_else(incompatible_md)?;
                    self.remaining = 0;
                }
                5 | 6 => self.remaining = 0,
                _ => return Err(incompatible_md()),
            }
        }
    }

    fn finish(mut self) -> Result<()> {
        if self.next_column()?.is_none() && self.query == self.sequence_len {
            Ok(())
        } else {
            Err(incompatible_md())
        }
    }
}

fn query_base(record: &impl EvidenceRecord, index: usize) -> Result<u8> {
    let code = usize::from(record.seq_nibble(index));
    NIBBLE_BASES
        .get(code)
        .copied()
        .filter(|base| *base != b'=')
        .ok_or_else(incompatible_md)
}

fn set_base(
    reference: &mut [u8],
    reference_start: usize,
    reference_length: usize,
    alignment_start: usize,
    offset: usize,
    base: u8,
) -> Result<()> {
    let absolute = alignment_start
        .checked_add(offset)
        .filter(|position| *position < reference_length)
        .ok_or_else(|| {
            RsomicsError::InvalidInput(
                "alignment evidence extends beyond the declared reference length".to_owned(),
            )
        })?;
    if absolute < reference_start || absolute >= reference_start + reference.len() {
        return Ok(());
    }
    let position = absolute - reference_start;
    match reference[position] {
        b'N' => reference[position] = base,
        existing if existing == base || base == b'N' => {}
        existing => {
            return Err(RsomicsError::InvalidInput(format!(
                "conflicting reference evidence at position {}: {} and {}",
                absolute + 1,
                char::from(existing),
                char::from(base)
            )));
        }
    }
    Ok(())
}

fn incompatible_md() -> RsomicsError {
    RsomicsError::InvalidInput("MD and CIGAR fields are incompatible".to_owned())
}
