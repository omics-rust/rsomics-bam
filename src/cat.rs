use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use noodles::sam;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::output::same_target;
use crate::{Program, bgzf_rewrite, header_source, hts_quickcheck, input};

const OUTPUT_BUFFER: usize = 2 * 1024 * 1024;

#[derive(Clone, Copy, Debug, Default)]
pub struct Options<'a> {
    pub header: Option<&'a Path>,
    pub destination: Option<&'a Path>,
    pub program: Option<Program<'a>>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub inputs: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output: Option<PathBuf>,
}

pub fn write<W: Write>(inputs: &[PathBuf], options: Options<'_>, output: W) -> Result<Summary> {
    if inputs.is_empty() {
        return Err(RsomicsError::ConfigError(
            "cat requires at least one BAM input".to_owned(),
        ));
    }

    let mut headers = Vec::with_capacity(inputs.len());
    for path in inputs {
        if path == Path::new("-") {
            return Err(RsomicsError::ConfigError(
                "cat requires named BAM inputs".to_owned(),
            ));
        }
        if input::detect_format(path)? != input::Format::Bam {
            return Err(RsomicsError::InvalidInput(format!(
                "{}: cat input must be BAM",
                path.display()
            )));
        }
        if let Some(destination) = options.destination
            && same_target(path, destination)?
        {
            return Err(RsomicsError::ConfigError(format!(
                "cat input and output must be different files: {}",
                path.display()
            )));
        }
        hts_quickcheck::require_bgzf_eof(path)?;
        headers.push(bgzf_rewrite::read_header(path)?);
    }

    let (mut header, mut text) = match options.header {
        Some(path) => {
            if let Some(destination) = options.destination
                && same_target(path, destination)?
            {
                return Err(RsomicsError::ConfigError(
                    "cat header source and output must be different files".to_owned(),
                ));
            }
            let source = header_source::read(path)?;
            (source.header, source.text)
        }
        None => headers[0].clone(),
    };
    for (path, (candidate, candidate_text)) in inputs.iter().zip(&headers) {
        require_same_dictionary(&header, candidate, path)?;
        for (id, read_group) in candidate.read_groups() {
            if !header.read_groups().contains_key(id) {
                let line = header_source::read_group_line(candidate_text, id).ok_or_else(|| {
                    RsomicsError::InvalidInput(format!(
                        "{}: read-group header line is missing",
                        path.display()
                    ))
                })?;
                header_source::append_line(&mut text, line);
                header
                    .read_groups_mut()
                    .insert(id.clone(), read_group.clone());
            }
        }
    }
    if let Some(program) = options.program {
        let id = program.add_to(&mut header)?;
        let canonical = bgzf_rewrite::canonical_header_text(&header)?;
        let line = header_source::program_line(&canonical, &id).ok_or_else(|| {
            RsomicsError::InvalidInput("serialized program header line is missing".to_owned())
        })?;
        header_source::append_line(&mut text, line);
    }

    let mut output = BufWriter::with_capacity(OUTPUT_BUFFER, output);
    bgzf_rewrite::write_header(&mut output, &header, &text)?;
    for path in inputs {
        bgzf_rewrite::copy_records(path, &mut output)?;
    }
    bgzf_rewrite::finish(&mut output)?;
    output.flush().map_err(RsomicsError::Io)?;

    Ok(Summary {
        inputs: inputs.len(),
        output: options.destination.map(Path::to_path_buf),
    })
}

fn require_same_dictionary(base: &sam::Header, candidate: &sam::Header, path: &Path) -> Result<()> {
    let base = base.reference_sequences();
    let candidate = candidate.reference_sequences();
    if base.len() == candidate.len()
        && base.iter().zip(candidate).all(
            |((base_name, base_reference), (candidate_name, candidate_reference))| {
                base_name == candidate_name
                    && base_reference.length() == candidate_reference.length()
            },
        )
    {
        Ok(())
    } else {
        Err(RsomicsError::InvalidInput(format!(
            "{}: reference sequence dictionary differs from the output header",
            path.display()
        )))
    }
}
