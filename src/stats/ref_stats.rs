use std::ops::Range;

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use super::reference::Reference;
use super::regions::Regions;

#[derive(Clone, Debug, Serialize)]
pub(crate) struct ReferenceStats {
    pub(crate) total_count: usize,
    pub(crate) output_count: usize,
    pub(crate) average_gc: f32,
    pub(crate) minimum_length: u64,
    pub(crate) maximum_length: u64,
    pub(crate) average_length: f32,
    pub(crate) total_length: u64,
    pub(crate) sequences: Vec<ReferenceSequenceStats>,
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct ReferenceSequenceStats {
    pub(crate) name: String,
    pub(crate) length: u64,
    pub(crate) gc: f32,
    pub(crate) unknown: i64,
}

impl ReferenceStats {
    pub(crate) fn collect(
        references: &[(Vec<u8>, u64)],
        selected: Option<&Regions>,
        mut reference: Option<&mut Reference>,
    ) -> Result<Self> {
        let mut sequences = Vec::new();
        if let Some(selected) = selected {
            for (reference_id, range) in selected.iter() {
                let (name, _) = &references[reference_id as usize];
                sequences.push(sequence_stats(name, range, reference.as_deref_mut())?);
            }
        } else {
            for (name, length) in references {
                sequences.push(sequence_stats(
                    name,
                    0..i64::try_from(*length).unwrap(),
                    reference.as_deref_mut(),
                )?);
            }
        }

        let output_count = sequences.len();
        let total_length = sequences.iter().map(|entry| entry.length).sum::<u64>();
        let minimum_length = sequences
            .iter()
            .map(|entry| entry.length)
            .min()
            .unwrap_or(0);
        let maximum_length = sequences
            .iter()
            .map(|entry| entry.length)
            .max()
            .unwrap_or(0);
        let average_length = if output_count == 0 {
            -1.0
        } else {
            total_length as f32 / output_count as f32
        };
        let average_gc = if reference.is_none() || output_count == 0 {
            -1.0
        } else {
            sequences.iter().map(|entry| entry.gc).sum::<f32>() / output_count as f32
        };
        Ok(Self {
            total_count: references.len(),
            output_count,
            average_gc,
            minimum_length,
            maximum_length,
            average_length,
            total_length,
            sequences,
        })
    }
}

fn sequence_stats(
    name: &[u8],
    range: Range<i64>,
    reference: Option<&mut Reference>,
) -> Result<ReferenceSequenceStats> {
    let full_sequence = reference
        .as_ref()
        .and_then(|reference| reference.length(name))
        .is_none_or(|length| range.start == 0 && range.end == length as i64);
    let label = if full_sequence {
        String::from_utf8_lossy(name).into_owned()
    } else {
        format!(
            "{}:{}-{}",
            String::from_utf8_lossy(name),
            range.start + 1,
            range.end
        )
    };
    let length = u64::try_from(range.end - range.start).unwrap();
    let Some(reference) = reference else {
        return Ok(ReferenceSequenceStats {
            name: label,
            length,
            gc: -1.0,
            unknown: -1,
        });
    };
    if !reference.contains(name) {
        return Err(RsomicsError::InvalidInput(format!(
            "reference FASTA is missing sequence {}",
            String::from_utf8_lossy(name)
        )));
    }
    let start = usize::try_from(range.start).unwrap();
    let end = usize::try_from(range.end).unwrap();
    if end > reference.length(name).unwrap() {
        return Err(RsomicsError::InvalidInput(format!(
            "reference FASTA sequence {} is shorter than the selected region",
            String::from_utf8_lossy(name)
        )));
    }
    let (gc, at, unknown) = reference.count_bases(name, start, end)?;
    Ok(ReferenceSequenceStats {
        name: label,
        length,
        gc: if gc + at == 0 {
            0.0
        } else {
            gc as f32 / (gc + at) as f32
        },
        unknown,
    })
}
