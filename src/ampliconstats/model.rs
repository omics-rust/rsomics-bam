use std::io::Write;

use rsomics_common::{Result, RsomicsError};

use crate::amplicon::{Primer, Strand};

#[derive(Clone, Debug)]
pub(super) struct Amplicon {
    pub(super) max_left: i64,
    pub(super) min_right: i64,
    pub(super) min_left: i64,
    pub(super) max_right: i64,
    pub(super) lefts: Vec<i64>,
    pub(super) rights: Vec<i64>,
}

impl Amplicon {
    fn new() -> Self {
        Self {
            max_left: 0,
            min_right: i64::MAX,
            min_left: i64::MAX,
            max_right: 0,
            lefts: Vec::new(),
            rights: Vec::new(),
        }
    }
}

pub(super) fn count(entries: &[Primer]) -> Result<usize> {
    groups(entries).map(|groups| groups.len())
}

pub(super) fn build(
    entries: &[Primer],
    output: &mut impl Write,
    reference: Option<&str>,
    first_number: usize,
    title: bool,
    max_length: usize,
) -> Result<Vec<Amplicon>> {
    let groups = groups(entries)?;
    if title {
        writeln!(output, "# Amplicon locations from BED file.")?;
        writeln!(
            output,
            "# LEFT/RIGHT are <start>-<end> format and comma-separated for alt-primers."
        )?;
        if reference.is_some() {
            writeln!(output, "#\n# AMPLICON\tREF\tNUMBER\tLEFT\tRIGHT")?;
        } else {
            writeln!(output, "#\n# AMPLICON\tNUMBER\tLEFT\tRIGHT")?;
        }
    }

    let mut amplicons = Vec::with_capacity(groups.len());
    for (group_index, (lefts, rights)) in groups.into_iter().enumerate() {
        if let Some(reference) = reference {
            write!(
                output,
                "AMPLICON\t{reference}\t{}",
                first_number + group_index + 1
            )?;
        } else {
            write!(output, "AMPLICON\t{}", first_number + group_index + 1)?;
        }
        for (index, primer) in lefts.iter().enumerate() {
            write!(
                output,
                "{}{}-{}",
                if index == 0 { "\t" } else { "," },
                primer.start + 1,
                primer.end
            )?;
        }
        for (index, primer) in rights.iter().enumerate() {
            write!(
                output,
                "{}{}-{}",
                if index == 0 { "\t" } else { "," },
                primer.start + 1,
                primer.end
            )?;
        }
        writeln!(output)?;

        let left_positions: Vec<_> = lefts.iter().map(|primer| primer.end).collect();
        let right_positions: Vec<_> = rights.iter().map(|primer| primer.start).collect();
        let mut amplicon = Amplicon::new();
        amplicon.max_left = left_positions.iter().copied().max().unwrap() + 1;
        amplicon.min_left = left_positions.iter().copied().min().unwrap() + 1;
        amplicon.min_right = right_positions.iter().copied().min().unwrap() - 1;
        amplicon.max_right = right_positions.iter().copied().max().unwrap() - 1;
        amplicon.lefts = left_positions;
        amplicon.rights = right_positions;
        if amplicon.max_left > amplicon.min_right {
            return Err(RsomicsError::InvalidInput(
                "amplicon primer intervals do not enclose a positive span".to_owned(),
            ));
        }
        let length = amplicon
            .max_right
            .checked_sub(amplicon.min_left)
            .and_then(|length| length.checked_add(2))
            .and_then(|length| usize::try_from(length).ok())
            .ok_or_else(|| {
                RsomicsError::InvalidInput("amplicon length is outside supported bounds".to_owned())
            })?;
        if length > max_length {
            return Err(RsomicsError::ConfigError(format!(
                "amplicon length {length} exceeds --max-amplicon-length {max_length}"
            )));
        }
        amplicons.push(amplicon);
    }
    Ok(amplicons)
}

fn groups(entries: &[Primer]) -> Result<Vec<(Vec<&Primer>, Vec<&Primer>)>> {
    if entries.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "primer reference has no rows".to_owned(),
        ));
    }
    let mut groups = Vec::new();
    let mut index = 0;
    while index < entries.len() {
        let mut lefts = Vec::new();
        while index < entries.len() && entries[index].strand == Some(Strand::Forward) {
            lefts.push(&entries[index]);
            index += 1;
        }
        if lefts.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "each amplicon must start with a forward primer".to_owned(),
            ));
        }
        let mut rights = Vec::new();
        while index < entries.len() && entries[index].strand == Some(Strand::Reverse) {
            rights.push(&entries[index]);
            index += 1;
        }
        if rights.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "each amplicon must end with a reverse primer".to_owned(),
            ));
        }
        groups.push((lefts, rights));
    }
    Ok(groups)
}
