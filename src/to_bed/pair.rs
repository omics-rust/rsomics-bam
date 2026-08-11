use std::io::Write;

use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};

use super::record::{PAIRED, READ1, READ2, REVERSE, UNMAPPED, project, reference_end, score};
use super::{PairScore, Score};

pub(super) struct State {
    pending: Option<RawRecord>,
}

impl State {
    pub fn new() -> Self {
        Self { pending: None }
    }

    pub fn push(
        &mut self,
        output: &mut impl Write,
        references: &[String],
        record: &RawRecord,
        score: PairScore,
        mate1_first: bool,
    ) -> Result<bool> {
        let Some(first) = self.pending.take() else {
            self.pending = Some(record.clone());
            return Ok(false);
        };
        write(output, references, &first, record, score, mate1_first)?;
        Ok(true)
    }

    pub fn finish(self) -> Result<()> {
        if let Some(record) = self.pending {
            return Err(RsomicsError::InvalidInput(format!(
                "BEDPE input ends with incomplete pair {}",
                String::from_utf8_lossy(record.name())
            )));
        }
        Ok(())
    }
}

struct End<'a> {
    reference: &'a str,
    start: i64,
    end: i64,
    strand: char,
}

fn write(
    output: &mut impl Write,
    references: &[String],
    first: &RawRecord,
    second: &RawRecord,
    score_mode: PairScore,
    mate1_first: bool,
) -> Result<()> {
    if first.name() != second.name() {
        return Err(RsomicsError::InvalidInput(format!(
            "BEDPE mates are not adjacent: {} followed by {}",
            String::from_utf8_lossy(first.name()),
            String::from_utf8_lossy(second.name())
        )));
    }
    if first.flags() & PAIRED == 0 || second.flags() & PAIRED == 0 {
        return Err(RsomicsError::InvalidInput(format!(
            "BEDPE records must be marked paired: {}",
            String::from_utf8_lossy(first.name())
        )));
    }
    let categories = (
        (first.flags() & READ1 != 0, first.flags() & READ2 != 0),
        (second.flags() & READ1 != 0, second.flags() & READ2 != 0),
    );
    if !matches!(
        categories,
        ((true, false), (false, true)) | ((false, true), (true, false))
    ) {
        return Err(RsomicsError::InvalidInput(format!(
            "BEDPE pair {} must contain one read1 and one read2",
            String::from_utf8_lossy(first.name())
        )));
    }

    let mut left_record = first;
    let mut right_record = second;
    let mut left = end(first, references)?;
    let mut right = end(second, references)?;
    let swap = if mate1_first {
        first.flags() & READ1 == 0
    } else {
        (left.reference, left.start) > (right.reference, right.start)
    };
    if swap {
        std::mem::swap(&mut left_record, &mut right_record);
        std::mem::swap(&mut left, &mut right);
    }

    let pair_score = match score_mode {
        PairScore::EditDistance => [left_record, right_record]
            .into_iter()
            .filter(|record| record.flags() & UNMAPPED == 0)
            .try_fold(0i64, |sum, record| {
                sum.checked_add(score(record, Score::EditDistance)?)
                    .ok_or_else(|| RsomicsError::InvalidInput("BEDPE score overflows".to_owned()))
            })?,
        PairScore::MappingQuality
            if left_record.flags() & UNMAPPED == 0 && right_record.flags() & UNMAPPED == 0 =>
        {
            i64::from(
                left_record
                    .mapping_quality()
                    .min(right_record.mapping_quality()),
            )
        }
        PairScore::MappingQuality => 0,
    };
    write!(
        output,
        "{}\t{}\t{}\t{}\t{}\t{}\t",
        left.reference, left.start, left.end, right.reference, right.start, right.end
    )
    .map_err(RsomicsError::Io)?;
    output.write_all(first.name()).map_err(RsomicsError::Io)?;
    writeln!(output, "\t{pair_score}\t{}\t{}", left.strand, right.strand).map_err(RsomicsError::Io)
}

fn end<'a>(record: &'a RawRecord, references: &'a [String]) -> Result<End<'a>> {
    let Some(mapped) = project(record, references)? else {
        return Ok(End {
            reference: ".",
            start: -1,
            end: -1,
            strand: '.',
        });
    };
    let cigar = record.decoded_cigar()?;
    Ok(End {
        reference: mapped.reference,
        start: i64::try_from(mapped.start)
            .map_err(|_| RsomicsError::InvalidInput("BEDPE start overflows".to_owned()))?,
        end: i64::try_from(reference_end(&cigar, mapped.start, mapped.name)?)
            .map_err(|_| RsomicsError::InvalidInput("BEDPE end overflows".to_owned()))?,
        strand: if mapped.flags & REVERSE != 0 {
            '-'
        } else {
            '+'
        },
    })
}
