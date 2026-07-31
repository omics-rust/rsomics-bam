use std::fs::File;
use std::io::{self, BufRead, BufWriter, IsTerminal, Write};
use std::path::PathBuf;

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::samples;

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Add column names before the results
    #[arg(short = 'H', long)]
    header: bool,

    /// Report whether each alignment has an index
    #[arg(short = 'i', long)]
    index: bool,

    /// Read-group tag to report
    #[arg(short = 'T', long, default_value = "SM")]
    tag: String,

    /// Write tabular output to a file
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Reference FASTA to match against alignment targets
    #[arg(short = 'f', long = "reference", value_name = "FASTA")]
    references: Vec<PathBuf>,

    /// File containing reference FASTA paths
    #[arg(short = 'F', long = "reference-list", value_name = "FILE")]
    reference_lists: Vec<PathBuf>,

    /// Pair each alignment with a custom index
    #[arg(short = 'X', long)]
    custom_index: bool,

    /// Input SAM, BAM, or CRAM files; read paths from stdin when omitted
    #[arg(value_name = "ALIGNMENT")]
    inputs: Vec<PathBuf>,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if json && arguments.output.is_some() {
        return Err(RsomicsError::ConfigError(
            "--json cannot be combined with --output".to_owned(),
        ));
    }

    let tag = samples::Tag::parse(&arguments.tag)?;
    let inputs = inputs(arguments.inputs, arguments.custom_index)?;
    let references = reference_paths(arguments.references, &arguments.reference_lists)?;
    let report = samples::collect(
        &inputs,
        samples::Options {
            tag,
            test_index: arguments.index,
            references: &references,
        },
    )?;

    if !json {
        match arguments.output {
            Some(path) => {
                let file = File::create(&path).map_err(|error| {
                    RsomicsError::Io(io::Error::new(
                        error.kind(),
                        format!("creating {}: {error}", path.display()),
                    ))
                })?;
                let mut output = BufWriter::new(file);
                report.write(arguments.header, &mut output)?;
                output.flush().map_err(RsomicsError::Io)?;
            }
            None => report.write(arguments.header, io::stdout().lock())?,
        }
    }

    Ok(CommandOutput::Samples { report })
}

fn inputs(inputs: Vec<PathBuf>, custom_index: bool) -> Result<Vec<samples::Input>> {
    if inputs.is_empty() {
        if io::stdin().is_terminal() {
            return Err(RsomicsError::ConfigError(
                "provide alignment paths or pipe one path per line to stdin".to_owned(),
            ));
        }
        let lines = io::stdin()
            .lock()
            .lines()
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(RsomicsError::Io)?;
        return lines
            .into_iter()
            .map(|line| {
                if custom_index {
                    let (alignment, index) = line.split_once('\t').ok_or_else(|| {
                        RsomicsError::InvalidInput(format!(
                            "expected ALIGNMENT<TAB>INDEX, got {line:?}"
                        ))
                    })?;
                    if alignment.is_empty() || index.is_empty() {
                        return Err(RsomicsError::InvalidInput(format!(
                            "expected ALIGNMENT<TAB>INDEX, got {line:?}"
                        )));
                    }
                    Ok(samples::Input {
                        path: PathBuf::from(alignment),
                        index: Some(PathBuf::from(index)),
                    })
                } else {
                    Ok(samples::Input {
                        path: PathBuf::from(line),
                        index: None,
                    })
                }
            })
            .collect();
    }

    if !custom_index {
        return Ok(inputs
            .into_iter()
            .map(|path| samples::Input { path, index: None })
            .collect());
    }

    if !inputs.len().is_multiple_of(2) {
        return Err(RsomicsError::ConfigError(
            "custom indexes require one index for every alignment".to_owned(),
        ));
    }
    let split = inputs.len() / 2;
    let (alignments, indexes) = inputs.split_at(split);
    Ok(alignments
        .iter()
        .cloned()
        .zip(indexes.iter().cloned())
        .map(|(path, index)| samples::Input {
            path,
            index: Some(index),
        })
        .collect())
}

fn reference_paths(mut references: Vec<PathBuf>, lists: &[PathBuf]) -> Result<Vec<PathBuf>> {
    for list in lists {
        let file = File::open(list).map_err(|error| {
            RsomicsError::Io(io::Error::new(
                error.kind(),
                format!("opening reference list {}: {error}", list.display()),
            ))
        })?;
        for line in io::BufReader::new(file).lines() {
            references.push(PathBuf::from(line.map_err(RsomicsError::Io)?));
        }
    }
    Ok(references)
}
