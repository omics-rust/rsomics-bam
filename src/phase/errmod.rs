#![allow(unsafe_code)]

use rsomics_common::{Result, RsomicsError};
use rust_htslib::htslib;

pub(super) struct Errmod(*mut htslib::errmod_t);

impl Errmod {
    pub(super) fn new() -> Result<Self> {
        let inner = unsafe { htslib::errmod_init(1.0 - 0.83) };
        if inner.is_null() {
            Err(RsomicsError::InvalidInput(
                "initializing the genotype likelihood model failed".to_owned(),
            ))
        } else {
            Ok(Self(inner))
        }
    }

    pub(super) fn call(&self, bases: &mut [u16]) -> Option<Call> {
        if bases.is_empty() {
            return None;
        }
        let mut likelihoods = [0.0_f32; 16];
        unsafe {
            htslib::errmod_cal(
                self.0,
                i32::try_from(bases.len()).ok()?,
                4,
                bases.as_mut_ptr(),
                likelihoods.as_mut_ptr(),
            );
        }

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
        (lower != upper).then_some(Call {
            alleles: [upper as u8, lower as u8],
            lod: (second - best + 0.499) as i32,
        })
    }
}

impl Drop for Errmod {
    fn drop(&mut self) {
        unsafe { htslib::errmod_destroy(self.0) };
    }
}

pub(super) struct Call {
    pub(super) alleles: [u8; 2],
    pub(super) lod: i32,
}
