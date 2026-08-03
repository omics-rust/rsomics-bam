use std::path::{Path, PathBuf};
use std::thread;

use noodles::{bam, cram, csi};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::hts_quickcheck::{require_bgzf_eof, require_cram_eof};
use crate::input;
use crate::output::{TransactionalFile, same_target};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum AlignmentFormat {
    Sam,
    Bam,
    Cram,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum IndexKind {
    Bai,
    Csi,
    Crai,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Options {
    pub kind: IndexKind,
    pub min_shift: u8,
    pub additional_threads: Option<usize>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub input: PathBuf,
    pub output: PathBuf,
    pub format: AlignmentFormat,
    pub kind: IndexKind,
    pub additional_threads: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_shift: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub depth: Option<u8>,
}

#[must_use]
pub fn default_output_path(input: &Path, kind: IndexKind) -> PathBuf {
    let mut output = input.as_os_str().to_os_string();
    output.push(match kind {
        IndexKind::Bai => ".bai",
        IndexKind::Csi => ".csi",
        IndexKind::Crai => ".crai",
    });
    PathBuf::from(output)
}

pub fn detect_format(input: &Path) -> Result<AlignmentFormat> {
    input::detect_format(input).map(|format| match format {
        input::Format::Sam => AlignmentFormat::Sam,
        input::Format::Bam => AlignmentFormat::Bam,
        input::Format::Cram => AlignmentFormat::Cram,
    })
}

#[must_use]
pub fn kind_for(format: AlignmentFormat, csi: bool) -> IndexKind {
    match format {
        AlignmentFormat::Sam | AlignmentFormat::Bam if csi => IndexKind::Csi,
        AlignmentFormat::Sam | AlignmentFormat::Bam => IndexKind::Bai,
        AlignmentFormat::Cram => IndexKind::Crai,
    }
}

pub fn create(input: &Path, output: &Path, options: Options) -> Result<Summary> {
    if output == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "index output must be a named file".to_owned(),
        ));
    }
    if same_target(input, output)? {
        return Err(RsomicsError::ConfigError(
            "alignment input and index output must be different files".to_owned(),
        ));
    }
    if options.kind == IndexKind::Csi && !(1..=30).contains(&options.min_shift) {
        return Err(RsomicsError::ConfigError(
            "CSI min-shift must be in 1..=30".to_owned(),
        ));
    }

    let format = detect_format(input)?;
    match (format, options.kind) {
        (AlignmentFormat::Sam | AlignmentFormat::Bam, IndexKind::Bai | IndexKind::Csi) => {
            require_bgzf_eof(input)?;
        }
        (AlignmentFormat::Cram, IndexKind::Crai) => require_cram_eof(input)?,
        _ => {
            return Err(RsomicsError::ConfigError(format!(
                "{} indexes cannot be written for {format:?} input",
                match options.kind {
                    IndexKind::Bai => "BAI",
                    IndexKind::Csi => "CSI",
                    IndexKind::Crai => "CRAI",
                }
            )));
        }
    }

    let additional_threads = options
        .additional_threads
        .unwrap_or_else(default_additional_threads);
    let threads = u32::try_from(additional_threads).map_err(|_| {
        RsomicsError::ConfigError("index thread count exceeds the supported range".to_owned())
    })?;
    let index_type = match options.kind {
        IndexKind::Csi => rust_htslib::bam::index::Type::Csi(u32::from(options.min_shift)),
        IndexKind::Bai | IndexKind::Crai => rust_htslib::bam::index::Type::Bai,
    };
    let transaction = TransactionalFile::new(output)?;
    rust_htslib::bam::index::build(
        input,
        Some(transaction.temporary_path()),
        index_type,
        threads,
    )
    .map_err(|error| {
        RsomicsError::InvalidInput(format!("{}: indexing alignment: {error}", input.display()))
    })?;

    let (min_shift, depth) = validate_index(transaction.temporary_path(), options.kind)?;
    transaction.commit()?;
    Ok(Summary {
        input: input.to_path_buf(),
        output: output.to_path_buf(),
        format,
        kind: options.kind,
        additional_threads,
        min_shift,
        depth,
    })
}

fn default_additional_threads() -> usize {
    thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1)
        .saturating_sub(1)
        .min(4)
}

fn validate_index(output: &Path, kind: IndexKind) -> Result<(Option<u8>, Option<u8>)> {
    match kind {
        IndexKind::Bai => bam::bai::fs::read(output)
            .map(|_| (Some(14), Some(5)))
            .map_err(RsomicsError::Io),
        IndexKind::Csi => csi::fs::read(output)
            .map(|index| {
                use noodles::csi::BinningIndex as _;
                (Some(index.min_shift()), Some(index.depth()))
            })
            .map_err(RsomicsError::Io),
        IndexKind::Crai => cram::crai::fs::read(output)
            .map(|_| (None, None))
            .map_err(RsomicsError::Io),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn automatic_threads_are_bounded() {
        assert!(default_additional_threads() <= 4);
    }
}
