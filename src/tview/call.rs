use rsomics_common::Result;

use crate::errmod::Errmod;

use super::grid::Entry;

pub(super) struct Caller {
    model: Errmod,
    bases: Vec<u16>,
}

pub(super) struct Call {
    pub(super) base: u8,
    pub(super) confidence: u8,
}

impl Caller {
    pub(super) fn new() -> Result<Self> {
        Ok(Self {
            model: Errmod::new(1.0 - 0.83)?,
            bases: Vec::new(),
        })
    }

    pub(super) fn call(&mut self, entries: &[Entry], reference: u8) -> Result<Call> {
        let reference_allele = allele(reference);
        self.bases.clear();
        let mut quality_sums = [0_i32; 4];
        for entry in entries {
            if entry.deletion || entry.reference_skip {
                continue;
            }
            let mut quality = entry.quality;
            if quality < 13 {
                continue;
            }
            let mapping_quality = if entry.mapping_quality == 255 {
                20
            } else {
                entry.mapping_quality.min(60)
            };
            quality = quality.min(mapping_quality).clamp(4, 63);
            let base = if entry.base == 0 {
                reference_allele
            } else {
                allele_nt16(entry.base)
            };
            self.bases
                .push(u16::from(quality) << 5 | u16::from(entry.reverse) << 4 | u16::from(base));
            if base < 4 {
                quality_sums[usize::from(base)] += i32::from(quality);
            }
        }
        let mut likelihoods = [0.0_f32; 25];
        self.model.calculate(&mut self.bases, 5, &mut likelihoods)?;
        let mut order = std::array::from_fn::<_, 4, _>(|index| (quality_sums[index], index));
        order.sort_by(|left, right| right.cmp(left));
        let first = order[0].1;
        let second = order[1].1;
        let prior = 30.0_f32;
        let mut probabilities = [
            likelihoods[first * 5 + first],
            likelihoods[first * 5 + second] + prior,
            likelihoods[second * 5 + second],
        ];
        if first != usize::from(reference_allele) {
            probabilities[0] += prior + 3.0;
        }
        if second != usize::from(reference_allele) {
            probabilities[2] += prior + 3.0;
        }
        let (code, confidence) =
            if probabilities[0] < probabilities[1] && probabilities[0] < probabilities[2] {
                (
                    1 << first,
                    probabilities[1].min(probabilities[2]) - probabilities[0],
                )
            } else if probabilities[2] < probabilities[1] && probabilities[2] < probabilities[0] {
                (
                    1 << second,
                    probabilities[0].min(probabilities[1]) - probabilities[2],
                )
            } else {
                (
                    1 << first | 1 << second,
                    probabilities[0].min(probabilities[2]) - probabilities[1],
                )
            };
        Ok(Call {
            base: b",ACMGRSVTWYHKDBN"[code],
            confidence: (((confidence + 0.499) as i32 / 10) + 1).clamp(1, 4) as u8,
        })
    }
}

fn allele(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

fn allele_nt16(base: u8) -> u8 {
    match base {
        1 => 0,
        2 => 1,
        4 => 2,
        8 => 3,
        _ => 4,
    }
}
