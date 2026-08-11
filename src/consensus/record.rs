use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};

use super::call::{BayesianMode, BayesianObservation};

const SOFT_CLIP: u8 = 4;
const HARD_CLIP: u8 = 5;

#[derive(Clone, Copy)]
pub(super) struct RecordOptions {
    pub(super) use_mapping_quality: bool,
    pub(super) adjust_quality: bool,
    pub(super) mismatch_halo: usize,
    pub(super) soft_clip_cost: u32,
    pub(super) homopolymer_fix: f64,
    pub(super) mode: BayesianMode,
}

impl Default for RecordOptions {
    fn default() -> Self {
        Self {
            use_mapping_quality: true,
            adjust_quality: true,
            mismatch_halo: 50,
            soft_clip_cost: 60,
            homopolymer_fix: 0.0,
            mode: BayesianMode::Recall,
        }
    }
}

pub(super) struct RecordState {
    qualities: Vec<u8>,
    local_mismatch_tenths: Vec<u32>,
    homopolymer: Vec<u8>,
    missing_quality: bool,
    has_padding: bool,
}

impl RecordState {
    pub(super) fn new(record: &RawRecord, options: RecordOptions) -> Result<Self> {
        let missing_quality = record.quality_scores().is_empty();
        let has_padding = record.cigar_ops().any(|(kind, _)| kind == 6);
        let mut qualities = if missing_quality {
            vec![u8::MAX; record.sequence_len()]
        } else {
            record.quality_scores().to_vec()
        };
        let mut local_mismatch_tenths = vec![0; record.sequence_len()];
        let mut homopolymer = vec![0; record.sequence_len()];
        if !options.use_mapping_quality || record.sequence_len() == 0 {
            return Ok(Self {
                qualities,
                local_mismatch_tenths,
                homopolymer,
                missing_quality,
                has_padding,
            });
        }

        if options.adjust_quality {
            adjust_local_quality(record, &qualities, &mut local_mismatch_tenths, options);
        }
        if options.homopolymer_fix != 0.0 && !missing_quality {
            redistribute_homopolymer_qualities(record, &mut qualities);
        }
        measure_homopolymers(record, &mut homopolymer);

        let Some(md) = record.aux_value(*b"MD") else {
            return Ok(Self {
                qualities,
                local_mismatch_tenths,
                homopolymer,
                missing_quality,
                has_padding,
            });
        };
        if record.aux_type(*b"MD") != Some(b'Z') {
            return Err(invalid_md(record, "MD is not a string"));
        }
        apply_soft_clip_costs(
            record,
            &mut local_mismatch_tenths,
            options.mismatch_halo,
            options.soft_clip_cost,
        )?;
        apply_md(
            record,
            md.strip_suffix(&[0]).unwrap_or(md),
            &mut local_mismatch_tenths,
            options.mismatch_halo,
        )?;

        Ok(Self {
            qualities,
            local_mismatch_tenths,
            homopolymer,
            missing_quality,
            has_padding,
        })
    }

    pub(super) fn quality(&self, query_position: usize) -> u8 {
        self.qualities
            .get(query_position)
            .copied()
            .unwrap_or(u8::MAX)
    }

    pub(super) fn has_padding(&self) -> bool {
        self.has_padding
    }

    pub(super) fn observation(
        &self,
        record: &RawRecord,
        base: u8,
        quality: u8,
        query_position: usize,
        reference_skip: bool,
    ) -> BayesianObservation {
        let local = self
            .local_mismatch_tenths
            .get(query_position)
            .or_else(|| self.local_mismatch_tenths.last())
            .copied()
            .unwrap_or(0);
        let homopolymer = self.homopolymer.get(query_position).copied().unwrap_or(0);
        BayesianObservation::new(base, quality).with_record_context(
            record.mapping_quality(),
            local,
            homopolymer,
            reference_skip,
            self.missing_quality,
        )
    }
}

fn adjust_local_quality(
    record: &RawRecord,
    qualities: &[u8],
    local: &mut [u32],
    options: RecordOptions,
) {
    let length = qualities.len();
    let mut base = record.seq_nibble(0);
    let mut homopolymer_left = 0usize;
    let mut homopolymer_minimum = qualities[0];
    for (index, &quality) in qualities.iter().enumerate().skip(1) {
        if record.seq_nibble(index) != base {
            break;
        }
        if index < 2 {
            homopolymer_minimum = homopolymer_minimum.min(quality);
        }
    }
    let mut window_minimum = qualities[..length.min(8)].iter().copied().min().unwrap();
    let mut index = length.min(8);
    let homopolymer_adjustment = if options.homopolymer_fix == 0.0 {
        1.0
    } else {
        options.homopolymer_fix
    };

    while index + 8 < length {
        let homopolymer_right =
            if options.homopolymer_fix != 0.0 && record.seq_nibble(index) != base {
                homopolymer_left = index;
                base = record.seq_nibble(index);
                homopolymer_minimum = qualities[index];
                let mut right = index + 1;
                while right < length && record.seq_nibble(right) == base {
                    if index < 2 {
                        homopolymer_minimum = homopolymer_minimum.min(qualities[right]);
                    }
                    right += 1;
                }
                right - 1
            } else {
                homopolymer_left
            };
        let homopolymer_span = homopolymer_right - homopolymer_left;
        let adjusted = if options.mode == BayesianMode::Compatibility116 {
            (f64::from(qualities[index]) + 5.0 * f64::from(window_minimum)) / 4.0
        } else {
            f64::from(qualities[index]) / 3.0
                + (f64::from(homopolymer_minimum) - 2.0 * homopolymer_span as f64)
                    * homopolymer_adjustment
        };
        if adjusted < f64::from(qualities[index]) {
            local[index] += (f64::from(qualities[index]) - adjusted) as u32;
        }

        homopolymer_minimum = qualities[index];
        let halo_start = index.saturating_sub(2).max(homopolymer_left);
        let halo_end = (index + 2).min(homopolymer_right);
        if halo_start <= halo_end {
            for &quality in &qualities[halo_start..=halo_end] {
                homopolymer_minimum = homopolymer_minimum.min(quality);
            }
        }
        if window_minimum > qualities[index + 8] {
            window_minimum = qualities[index + 8];
        } else if window_minimum <= qualities[index - 8] {
            window_minimum = qualities[index - 7..=index + 8]
                .iter()
                .copied()
                .min()
                .unwrap();
        }
        index += 1;
    }

    while index < length {
        let adjusted = if options.mode == BayesianMode::Compatibility116 {
            (f64::from(qualities[index]) + 5.0 * f64::from(window_minimum)) / 4.0
        } else {
            f64::from(qualities[index]) / 3.0
                + f64::from(homopolymer_minimum) * homopolymer_adjustment
        };
        if adjusted < f64::from(qualities[index]) {
            local[index] += (f64::from(qualities[index]) - adjusted) as u32;
        }
        index += 1;
    }
}

fn redistribute_homopolymer_qualities(record: &RawRecord, qualities: &mut [u8]) {
    let mut start = 0usize;
    while start < qualities.len() {
        let base = record.seq_nibble(start);
        let mut end = start + 1;
        while end < qualities.len() && record.seq_nibble(end) == base {
            end += 1;
        }
        let mut left = start;
        let mut right = end.saturating_sub(1);
        while left < right {
            let error = 10.0f64.powf(-f64::from(qualities[left]) / 10.0)
                + 10.0f64.powf(-f64::from(qualities[right]) / 10.0);
            let quality = (-fast_log2(error / 2.0) * 3.0104 + 0.49) as u8;
            qualities[left] = quality;
            qualities[right] = quality;
            left += 1;
            right -= 1;
        }
        start = end;
    }
}

fn measure_homopolymers(record: &RawRecord, output: &mut [u8]) {
    let mut start = 0usize;
    while start < output.len() {
        let base = record.seq_nibble(start);
        let mut end = start + 1;
        while end < output.len() && record.seq_nibble(end) == base {
            end += 1;
        }
        output[start..end].fill((end - start - 1).min(100) as u8);
        start = end;
    }
}

fn apply_soft_clip_costs(
    record: &RawRecord,
    local: &mut [u32],
    halo: usize,
    cost: u32,
) -> Result<()> {
    let cigar = record.decoded_cigar()?;
    let left_clipped = cigar.first().is_some_and(|&(kind, _)| kind == SOFT_CLIP)
        || matches!(cigar.as_slice(), [(HARD_CLIP, _), (SOFT_CLIP, _), ..]);
    let right_clipped = cigar.last().is_some_and(|&(kind, _)| kind == SOFT_CLIP)
        || matches!(cigar.as_slice(), [.., (SOFT_CLIP, _), (HARD_CLIP, _)]);
    if left_clipped {
        for value in local.iter_mut().take(halo) {
            *value += cost;
        }
        for value in local.iter_mut().skip(halo).take(halo) {
            *value += cost >> 1;
        }
    }
    if right_clipped {
        let length = local.len();
        for value in local.iter_mut().skip(length.saturating_sub(halo)) {
            *value += cost;
        }
        for value in local
            .iter_mut()
            .skip(length.saturating_sub(halo * 2))
            .take(halo)
        {
            *value += cost >> 1;
        }
    }
    Ok(())
}

fn apply_md(record: &RawRecord, md: &[u8], local: &mut [u32], halo: usize) -> Result<()> {
    let mut cursor = 0usize;
    let mut position = 0usize;
    while cursor < md.len() {
        if md[cursor].is_ascii_digit() {
            let start = cursor;
            while cursor < md.len() && md[cursor].is_ascii_digit() {
                cursor += 1;
            }
            let count = std::str::from_utf8(&md[start..cursor])
                .ok()
                .and_then(|value| value.parse::<usize>().ok())
                .ok_or_else(|| invalid_md(record, "MD match length is invalid"))?;
            position = position
                .checked_add(count)
                .ok_or_else(|| invalid_md(record, "MD coordinate overflows"))?;
            continue;
        }
        if md[cursor] == b'^' {
            cursor += 1;
            let start = cursor;
            while cursor < md.len() && md[cursor].is_ascii_alphabetic() {
                cursor += 1;
            }
            if cursor == start {
                return Err(invalid_md(record, "MD deletion is empty"));
            }
            continue;
        }
        if !md[cursor].is_ascii_alphabetic() {
            return Err(invalid_md(record, "MD contains an invalid byte"));
        }
        add_mismatch(local, position, halo);
        cursor += 1;
    }
    Ok(())
}

fn add_mismatch(local: &mut [u32], position: usize, halo: usize) {
    let length = local.len();
    let far_left = position.saturating_sub(halo.saturating_mul(2));
    let near_left = position.saturating_sub(halo);
    let near_right = position.saturating_add(halo).min(local.len());
    let far_right = position
        .saturating_add(halo.saturating_mul(2))
        .min(local.len());
    for value in &mut local[far_left.min(length)..near_left.min(length)] {
        *value += 5;
    }
    for value in &mut local[near_left.min(length)..near_right] {
        *value += 10;
    }
    for value in &mut local[near_right..far_right] {
        *value += 5;
    }
}

fn fast_log2(value: f64) -> f64 {
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 2047) as i32 - 1024;
    let mantissa = f64::from_bits((bits & !(2047_u64 << 52)) + (1023_u64 << 52));
    f64::from(exponent) + ((-1.0 / 3.0 * mantissa + 2.0) * mantissa - 2.0 / 3.0)
}

fn invalid_md(record: &RawRecord, reason: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "read {}: {reason}",
        String::from_utf8_lossy(record.name())
    ))
}

#[cfg(test)]
pub(super) fn test_record(md: &[u8]) -> RawRecord {
    test_record_with(b"AAACCCCCGT", &[(0, 10)], md)
}

#[cfg(test)]
pub(super) fn test_record_with(bases: &[u8], cigar: &[(u8, u32)], md: &[u8]) -> RawRecord {
    let name = b"state";
    let mut payload = Vec::new();
    payload.extend_from_slice(&0i32.to_le_bytes());
    payload.extend_from_slice(&0i32.to_le_bytes());
    payload.push((name.len() + 1) as u8);
    payload.push(60);
    payload.extend_from_slice(&0u16.to_le_bytes());
    payload.extend_from_slice(&(cigar.len() as u16).to_le_bytes());
    payload.extend_from_slice(&0u16.to_le_bytes());
    payload.extend_from_slice(&(bases.len() as u32).to_le_bytes());
    payload.extend_from_slice(&(-1i32).to_le_bytes());
    payload.extend_from_slice(&(-1i32).to_le_bytes());
    payload.extend_from_slice(&0i32.to_le_bytes());
    payload.extend_from_slice(name);
    payload.push(0);
    for &(kind, length) in cigar {
        payload.extend_from_slice(&(length << 4 | u32::from(kind)).to_le_bytes());
    }
    for pair in bases.chunks(2) {
        let code = |base| match base {
            b'A' => 1,
            b'C' => 2,
            b'G' => 4,
            b'T' => 8,
            _ => 15,
        };
        payload.push(code(pair[0]) << 4 | pair.get(1).copied().map_or(0, code));
    }
    payload.extend(std::iter::repeat_n(30, bases.len()));
    payload.extend_from_slice(b"MDZ");
    payload.extend_from_slice(md);
    payload.push(0);
    RawRecord::try_from(payload).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn retains_local_mismatch_and_homopolymer_context() {
        let record = test_record(b"2A7");
        let state = RecordState::new(
            &record,
            RecordOptions {
                adjust_quality: false,
                mismatch_halo: 2,
                ..RecordOptions::default()
            },
        )
        .unwrap();

        assert_eq!(
            state.local_mismatch_tenths,
            [10, 10, 10, 10, 5, 5, 0, 0, 0, 0]
        );
        assert_eq!(state.homopolymer, [2, 2, 2, 4, 4, 4, 4, 4, 0, 0]);
        assert_eq!(state.quality(3), 30);
        let observation = state.observation(&record, 2, 30, 3, false);
        assert_eq!(observation.local_mismatch_tenths, 10);
        assert_eq!(observation.homopolymer, 4);
    }

    #[test]
    fn adjusts_long_records_without_homopolymer_redistribution() {
        let bases = b"ACGTACGTACGTACGTACGTACGT";
        let record = test_record_with(bases, &[(0, bases.len() as u32)], b"24");

        let state = RecordState::new(&record, RecordOptions::default()).unwrap();

        assert_eq!(state.local_mismatch_tenths, [0; 24]);
    }

    #[test]
    fn rejects_malformed_md() {
        let record = test_record(b"2^");
        let error = RecordState::new(
            &record,
            RecordOptions {
                adjust_quality: false,
                ..RecordOptions::default()
            },
        )
        .err()
        .unwrap();

        assert!(error.to_string().contains("MD deletion is empty"));
    }
}
