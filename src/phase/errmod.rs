use rsomics_common::Result;

use crate::errmod::Errmod as Model;

pub(super) struct Errmod(Model);

impl Errmod {
    pub(super) fn new() -> Result<Self> {
        Model::new(1.0 - 0.83).map(Self)
    }

    pub(super) fn call(&self, bases: &mut [u16]) -> Result<Option<Call>> {
        if bases.is_empty() {
            return Ok(None);
        }
        let mut likelihoods = [0.0_f32; 16];
        self.0.calculate(bases, 4, &mut likelihoods)?;

        let mut best = f32::INFINITY;
        let mut second = f32::INFINITY;
        let mut genotype = 0;
        for first in 0..4 {
            for second_allele in first..4 {
                let index = first * 4 + second_allele;
                let value = likelihoods[index];
                if value < best {
                    genotype = index;
                    second = best;
                    best = value;
                } else if value < second {
                    second = value;
                }
            }
        }
        let lower = genotype >> 2;
        let upper = genotype & 3;
        Ok((lower != upper).then_some(Call {
            alleles: [upper as u8, lower as u8],
            lod: (second - best + 0.499) as i32,
        }))
    }
}

pub(super) struct Call {
    pub(super) alleles: [u8; 2],
    pub(super) lod: i32,
}
