use std::io::Write;
use std::path::{Path, PathBuf};

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::output::same_target;
use crate::{Program, bgzf_rewrite, header_source, hts_quickcheck, input};

#[derive(Clone, Copy, Debug, Default)]
pub struct Options<'a> {
    pub destination: Option<&'a Path>,
    pub program: Option<Program<'a>>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub input: PathBuf,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output: Option<PathBuf>,
    pub reference_sequences: usize,
}

pub fn write<W: Write>(
    header_source: &Path,
    bam: &Path,
    options: Options<'_>,
    mut output: W,
) -> Result<Summary> {
    if bam == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "reheader requires a named BAM input".to_owned(),
        ));
    }
    if input::detect_format(bam)? != input::Format::Bam {
        return Err(RsomicsError::InvalidInput(format!(
            "{}: reheader input must be BAM",
            bam.display()
        )));
    }
    if let Some(destination) = options.destination {
        if same_target(bam, destination)? {
            return Err(RsomicsError::ConfigError(
                "reheader input and output must be different files".to_owned(),
            ));
        }
        if same_target(header_source, destination)? {
            return Err(RsomicsError::ConfigError(
                "reheader header source and output must be different files".to_owned(),
            ));
        }
    }

    hts_quickcheck::require_bgzf_eof(bam)?;
    let (original, _) = bgzf_rewrite::read_header(bam)?;
    let source = header_source::read(header_source)?;
    let mut replacement = source.header;
    let mut text = source.text;
    let expected = original.reference_sequences().len();
    let actual = replacement.reference_sequences().len();
    if actual != expected {
        return Err(RsomicsError::InvalidInput(format!(
            "replacement header has {actual} reference sequences; input has {expected}"
        )));
    }
    if let Some(program) = options.program {
        let id = program.add_to(&mut replacement)?;
        let canonical = bgzf_rewrite::canonical_header_text(&replacement)?;
        let line = header_source::program_line(&canonical, &id).ok_or_else(|| {
            RsomicsError::InvalidInput("serialized program header line is missing".to_owned())
        })?;
        header_source::append_line(&mut text, line);
    }

    bgzf_rewrite::write_header(&mut output, &replacement, &text)?;
    bgzf_rewrite::copy_records(bam, &mut output)?;
    bgzf_rewrite::finish(&mut output)?;
    output.flush().map_err(RsomicsError::Io)?;

    Ok(Summary {
        input: bam.to_path_buf(),
        output: options.destination.map(Path::to_path_buf),
        reference_sequences: actual,
    })
}
