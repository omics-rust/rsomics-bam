#![allow(unsafe_code)]

use rsomics_common::{Result, RsomicsError};
use rust_htslib::htslib;

pub(crate) struct Errmod(*mut htslib::errmod_t);

impl Errmod {
    pub(crate) fn new(error_dependency: f64) -> Result<Self> {
        let inner = unsafe { htslib::errmod_init(error_dependency) };
        if inner.is_null() {
            Err(RsomicsError::InvalidInput(
                "initializing the genotype likelihood model failed".to_owned(),
            ))
        } else {
            Ok(Self(inner))
        }
    }

    pub(crate) fn calculate(
        &self,
        bases: &mut [u16],
        allele_count: usize,
        likelihoods: &mut [f32],
    ) -> Result<()> {
        if likelihoods.len() != allele_count * allele_count {
            return Err(RsomicsError::ConfigError(
                "genotype likelihood output has the wrong size".to_owned(),
            ));
        }
        likelihoods.fill(0.0);
        if bases.is_empty() {
            return Ok(());
        }
        unsafe {
            htslib::errmod_cal(
                self.0,
                i32::try_from(bases.len()).map_err(|_| {
                    RsomicsError::InvalidInput("pileup depth exceeds errmod limits".to_owned())
                })?,
                i32::try_from(allele_count)
                    .map_err(|_| RsomicsError::ConfigError("too many errmod alleles".to_owned()))?,
                bases.as_mut_ptr(),
                likelihoods.as_mut_ptr(),
            );
        }
        Ok(())
    }
}

impl Drop for Errmod {
    fn drop(&mut self) {
        unsafe { htslib::errmod_destroy(self.0) };
    }
}
