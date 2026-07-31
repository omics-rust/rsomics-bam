use std::path::{Path, PathBuf};

use noodles::{bam, cram, csi, fasta, sam};
use noodles_util::alignment;
use rsomics_common::{Result, RsomicsError};

pub(crate) struct ReferenceDictionary {
    pub path: PathBuf,
    pub targets: Vec<(Vec<u8>, u64)>,
}

pub(crate) fn load_reference(path: &Path) -> Result<ReferenceDictionary> {
    let mut index_path = path.as_os_str().to_os_string();
    index_path.push(".fai");
    let index_path = PathBuf::from(index_path);
    let index = fasta::fai::fs::read(&index_path).map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "loading reference index {}: {error}",
            index_path.display()
        ))
    })?;
    let records: &[fasta::fai::Record] = index.as_ref();
    let targets = records
        .iter()
        .map(|record| (record.name().to_vec(), record.length()))
        .collect();

    Ok(ReferenceDictionary {
        path: path.to_path_buf(),
        targets,
    })
}

pub(crate) fn has_index(
    format: crate::input::Format,
    alignment_path: &Path,
    index_path: Option<&Path>,
) -> bool {
    let Some(index_path) = index_path else {
        return alignment::io::indexed_reader::Builder::default()
            .build_from_path(alignment_path)
            .is_ok();
    };

    match format {
        crate::input::Format::Sam => csi::fs::read(index_path).is_ok(),
        crate::input::Format::Bam => {
            bam::bai::fs::read(index_path).is_ok() || csi::fs::read(index_path).is_ok()
        }
        crate::input::Format::Cram => cram::crai::fs::read(index_path).is_ok(),
    }
}

pub(crate) fn header_targets(header: &sam::Header) -> Vec<(Vec<u8>, u64)> {
    header
        .reference_sequences()
        .iter()
        .map(|(name, sequence)| {
            (
                name.to_vec(),
                u64::try_from(usize::from(sequence.length())).unwrap(),
            )
        })
        .collect()
}
