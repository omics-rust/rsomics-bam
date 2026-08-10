use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use super::record::{BaseCounts, QualityCycles};

#[derive(Clone, Debug, Serialize)]
pub(crate) struct BarcodeStats {
    pub(crate) sequence_tag: [u8; 2],
    pub(crate) quality_tag: [u8; 2],
    pub(crate) bases: Vec<BaseCounts>,
    pub(crate) qualities: QualityCycles,
    pub(crate) separator: Option<usize>,
    pub(crate) maximum_quality: Option<usize>,
}

impl BarcodeStats {
    pub(crate) fn new(sequence_tag: [u8; 2], quality_tag: [u8; 2]) -> Self {
        Self {
            sequence_tag,
            quality_tag,
            bases: Vec::new(),
            qualities: QualityCycles::default(),
            separator: None,
            maximum_quality: None,
        }
    }

    pub(crate) fn collect(
        &mut self,
        name: &[u8],
        sequence: Option<&[u8]>,
        quality: Option<&[u8]>,
    ) -> Result<()> {
        let Some(sequence) = sequence.filter(|value| !value.is_empty()) else {
            return Ok(());
        };
        if self.bases.is_empty() {
            self.bases.try_reserve(sequence.len()).map_err(|_| {
                RsomicsError::InvalidInput(
                    "barcode cycle count exceeds available memory".to_owned(),
                )
            })?;
            self.bases.resize(sequence.len(), BaseCounts::default());
            self.qualities.ensure_length(sequence.len())?;
        } else if sequence.len() != self.bases.len() {
            return Err(invalid(
                name,
                self.sequence_tag,
                "barcode length differs from earlier records",
            ));
        }

        let mut separator = None;
        for (index, &base) in sequence.iter().enumerate() {
            match base.to_ascii_uppercase() {
                b'A' => self.bases[index].a += 1,
                b'C' => self.bases[index].c += 1,
                b'G' => self.bases[index].g += 1,
                b'T' => self.bases[index].t += 1,
                b'N' => self.bases[index].n += 1,
                _ if separator.is_none() => separator = Some(index),
                _ => {
                    return Err(invalid(
                        name,
                        self.sequence_tag,
                        "barcode contains multiple separators",
                    ));
                }
            }
        }
        if let Some(previous) = self.separator {
            if separator != Some(previous) {
                return Err(invalid(
                    name,
                    self.sequence_tag,
                    "barcode separator is inconsistent",
                ));
            }
        } else if separator.is_some() {
            self.separator = separator;
        }

        let Some(quality) = quality else {
            return Ok(());
        };
        if quality.len() != sequence.len() {
            return Err(invalid(
                name,
                self.quality_tag,
                "barcode sequence and quality lengths differ",
            ));
        }
        for (index, &value) in quality.iter().enumerate() {
            let Some(score) = value.checked_sub(b'!') else {
                continue;
            };
            let score = usize::from(score);
            self.qualities.increment(index, score);
            self.maximum_quality = Some(self.maximum_quality.map_or(score, |old| old.max(score)));
        }
        Ok(())
    }
}

fn invalid(name: &[u8], tag: [u8; 2], reason: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "record {} tag {}: {reason}",
        String::from_utf8_lossy(name),
        String::from_utf8_lossy(&tag)
    ))
}
