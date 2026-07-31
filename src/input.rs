use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rust_htslib::bam::{Read, Reader};

pub(crate) fn open(
    input: &Path,
    reference: Option<&Path>,
    additional_threads: usize,
) -> Result<Reader> {
    let mut reader = if input == Path::new("-") {
        Reader::from_stdin()
    } else {
        Reader::from_path(input)
    }
    .map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "opening alignment input {}: {error}",
            input.display()
        ))
    })?;

    if let Some(reference) = reference {
        reader.set_reference(reference).map_err(|error| {
            RsomicsError::ConfigError(format!(
                "setting reference {}: {error}",
                reference.display()
            ))
        })?;
    }

    if additional_threads > 0 {
        reader.set_threads(additional_threads).map_err(|error| {
            RsomicsError::ConfigError(format!(
                "configuring {additional_threads} additional input threads: {error}"
            ))
        })?;
    }

    Ok(reader)
}
