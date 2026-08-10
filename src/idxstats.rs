use std::ffi::OsString;
use std::io::Write;
use std::path::{Path, PathBuf};

use rsomics_common::{Result, RsomicsError};
use rust_htslib::bam::Read as _;
use serde::Serialize;

use crate::input;

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct ReferenceStats {
    pub name: String,
    pub length: u64,
    pub mapped: u64,
    pub unmapped: u64,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Report {
    pub references: Vec<ReferenceStats>,
    pub unplaced_unmapped: u64,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Options<'a> {
    pub reference: Option<&'a Path>,
    pub index: Option<&'a Path>,
    pub additional_threads: usize,
}

pub fn collect(input_path: &Path, options: Options<'_>) -> Result<Report> {
    let mut reader = input::open(input_path, options.reference, options.additional_threads)?;
    let header = reader.read_header(input_path)?;
    let mut report = Report {
        references: header
            .reference_sequences()
            .iter()
            .map(|(name, reference)| ReferenceStats {
                name: String::from_utf8_lossy(name.as_ref()).into_owned(),
                length: u64::try_from(usize::from(reference.length())).unwrap(),
                mapped: 0,
                unmapped: 0,
            })
            .collect(),
        unplaced_unmapped: 0,
    };

    let index = options
        .index
        .map(Path::to_path_buf)
        .or_else(|| find_index(input_path));
    if let Some(index) = index {
        return collect_from_index(input_path, &index, options, report);
    }

    let mut last_reference_id = -2;
    reader.visit_owned_raw_records(&header, input_path, |record| {
        let raw_reference_id = record.reference_sequence_id();
        let reference_id = usize::try_from(raw_reference_id).ok();
        if raw_reference_id != last_reference_id {
            if last_reference_id >= -1 && already_seen(&report, reference_id) {
                return Err(RsomicsError::InvalidInput(format!(
                    "{} is not coordinate sorted: alignments for one reference occur in multiple blocks",
                    input_path.display()
                )));
            }
            last_reference_id = raw_reference_id;
        }
        if record.flags() & 0x04 == 0 {
            let reference_id = reference_id.ok_or_else(|| {
                RsomicsError::InvalidInput("mapped alignment has no reference".to_owned())
            })?;
            let reference = report.references.get_mut(reference_id).ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "alignment reference ID {reference_id} is absent from the header"
                ))
            })?;
            reference.mapped += 1;
        } else if let Some(reference_id) = reference_id {
            let reference = report.references.get_mut(reference_id).ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "alignment reference ID {reference_id} is absent from the header"
                ))
            })?;
            reference.unmapped += 1;
        } else {
            report.unplaced_unmapped += 1;
        }
        Ok(true)
    })?;
    Ok(report)
}

fn collect_from_index(
    input: &Path,
    index: &Path,
    options: Options<'_>,
    mut report: Report,
) -> Result<Report> {
    let mut reader =
        rust_htslib::bam::IndexedReader::from_path_and_index(input, index).map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "loading index {} for {}: {error}",
                index.display(),
                input.display()
            ))
        })?;
    if let Some(reference) = options.reference {
        reader.set_reference(reference).map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "setting reference {} for {}: {error}",
                reference.display(),
                input.display()
            ))
        })?;
    }
    if options.additional_threads > 0 {
        reader
            .set_threads(options.additional_threads)
            .map_err(|error| {
                RsomicsError::ConfigError(format!(
                    "configuring {} decoding threads for {}: {error}",
                    options.additional_threads,
                    input.display()
                ))
            })?;
    }
    let stats = reader.index_stats().map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "reading index statistics from {} for {}: {error}",
            index.display(),
            input.display()
        ))
    })?;
    let mut seen = vec![false; report.references.len()];
    let mut saw_unplaced = false;
    for (reference_id, length, mapped, unmapped) in stats {
        if reference_id == -1 {
            if saw_unplaced || mapped != 0 {
                return Err(invalid_index_stats(input, index));
            }
            saw_unplaced = true;
            report.unplaced_unmapped = unmapped;
            continue;
        }
        let reference_id = usize::try_from(reference_id)
            .ok()
            .filter(|&reference_id| reference_id < report.references.len())
            .ok_or_else(|| invalid_index_stats(input, index))?;
        if seen[reference_id] || report.references[reference_id].length != length {
            return Err(invalid_index_stats(input, index));
        }
        seen[reference_id] = true;
        report.references[reference_id].mapped = mapped;
        report.references[reference_id].unmapped = unmapped;
    }
    if seen.iter().any(|seen| !seen) || !saw_unplaced {
        return Err(invalid_index_stats(input, index));
    }
    Ok(report)
}

fn already_seen(report: &Report, reference_id: Option<usize>) -> bool {
    match reference_id {
        Some(reference_id) => report
            .references
            .get(reference_id)
            .is_some_and(|reference| reference.mapped > 0 || reference.unmapped > 0),
        None => report.unplaced_unmapped > 0,
    }
}

fn find_index(input: &Path) -> Option<PathBuf> {
    let mut candidates = Vec::with_capacity(6);
    for extension in ["bai", "csi", "crai"] {
        let mut appended = OsString::from(input.as_os_str());
        appended.push(format!(".{extension}"));
        candidates.push(PathBuf::from(appended));
        candidates.push(input.with_extension(extension));
    }
    candidates.into_iter().find(|path| path.exists())
}

fn invalid_index_stats(input: &Path, index: &Path) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "index statistics in {} do not match {}",
        index.display(),
        input.display()
    ))
}

impl Report {
    pub fn write(&self, mut output: impl Write) -> Result<()> {
        for reference in &self.references {
            writeln!(
                output,
                "{}\t{}\t{}\t{}",
                reference.name, reference.length, reference.mapped, reference.unmapped
            )
            .map_err(RsomicsError::Io)?;
        }
        writeln!(output, "*\t0\t0\t{}", self.unplaced_unmapped).map_err(RsomicsError::Io)?;
        output.flush().map_err(RsomicsError::Io)
    }
}
