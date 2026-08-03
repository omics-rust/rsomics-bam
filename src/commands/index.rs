use std::path::{Path, PathBuf};

use clap::{ArgAction, Args};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::index::{self, AlignmentFormat, IndexKind, Options};
use crate::output::same_target;

#[derive(Debug, Args)]
#[command(after_help = "\
Examples:
  rsomics-bam index sorted.bam
  rsomics-bam index -c sorted.bam
  rsomics-bam index -m 18 sorted.bam
  rsomics-bam index -M first.bam second.cram
  rsomics-bam index -o custom.bai sorted.bam")]
pub(crate) struct Arguments {
    /// Input coordinate-sorted SAM, BAM, or CRAM; a second path is a legacy index output
    #[arg(value_name = "ALIGNMENT_OR_INDEX", required = true, num_args = 1..)]
    paths: Vec<PathBuf>,

    /// Write BAI for BAM or BGZF SAM; this is the default
    #[arg(short = 'b', long, conflicts_with_all = ["csi", "min_shift"])]
    bai: bool,

    /// Write CSI for BAM or BGZF SAM
    #[arg(short = 'c', long, conflicts_with = "bai")]
    csi: bool,

    /// Smallest CSI bin is 2^INT bases; implies --csi
    #[arg(
        short = 'm',
        long,
        value_name = "INT",
        value_parser = clap::value_parser!(u8).range(1..=30),
        conflicts_with = "bai"
    )]
    min_shift: Option<u8>,

    /// Treat every positional path as an input and write default sibling indexes
    #[arg(short = 'M', long, action = ArgAction::SetTrue, conflicts_with = "output")]
    multiple: bool,

    /// Write a single index to FILE
    #[arg(short = 'o', long, value_name = "FILE", conflicts_with = "multiple")]
    output: Option<PathBuf>,

    /// Additional index workers; omit for automatic parallelism, or use 0 for one thread
    #[arg(short = '@', long, value_name = "INT")]
    threads: Option<usize>,
}

struct Job {
    input: PathBuf,
    output: PathBuf,
    format: AlignmentFormat,
    kind: IndexKind,
}

pub(crate) fn execute(arguments: Arguments, _json: bool) -> Result<CommandOutput> {
    let min_shift = arguments.min_shift.unwrap_or(14);
    let csi = !arguments.bai && (arguments.csi || arguments.min_shift.is_some());
    let (inputs, positional_output) = split_paths(&arguments)?;
    if arguments.output.is_some() && positional_output.is_some() {
        return Err(RsomicsError::ConfigError(
            "use either a positional index output or --output, not both".to_owned(),
        ));
    }
    let explicit_output = arguments.output.or(positional_output);
    let mut jobs = Vec::with_capacity(inputs.len());
    for input in inputs {
        if input == Path::new("-") {
            return Err(RsomicsError::ConfigError(
                "index input must be a named file".to_owned(),
            ));
        }
        let format = index::detect_format(&input)?;
        let kind = index::kind_for(format, csi);
        let output = explicit_output
            .clone()
            .unwrap_or_else(|| index::default_output_path(&input, kind));
        jobs.push(Job {
            input,
            output,
            format,
            kind,
        });
    }
    validate_jobs(&jobs)?;

    let mut summaries = Vec::with_capacity(jobs.len());
    for job in jobs {
        let summary = index::create(
            &job.input,
            &job.output,
            Options {
                kind: job.kind,
                min_shift,
                additional_threads: arguments.threads,
            },
        )?;
        debug_assert_eq!(summary.format, job.format);
        summaries.push(summary);
    }
    Ok(CommandOutput::Index { summaries })
}

fn split_paths(arguments: &Arguments) -> Result<(Vec<PathBuf>, Option<PathBuf>)> {
    if arguments.multiple {
        return Ok((arguments.paths.clone(), None));
    }
    match arguments.paths.as_slice() {
        [input] => Ok((vec![input.clone()], None)),
        [input, output] => {
            if existing_alignment(output) {
                return Err(RsomicsError::ConfigError(
                    "multiple alignment inputs require --multiple".to_owned(),
                ));
            }
            Ok((vec![input.clone()], Some(output.clone())))
        }
        _ => Err(RsomicsError::ConfigError(
            "multiple alignment inputs require --multiple".to_owned(),
        )),
    }
}

fn existing_alignment(path: &Path) -> bool {
    if !path.exists() {
        return false;
    }
    let Ok(mut reader) = crate::input::open(path, None, 0) else {
        return false;
    };
    reader.read_header(path).is_ok()
}

fn validate_jobs(jobs: &[Job]) -> Result<()> {
    for (index, job) in jobs.iter().enumerate() {
        if job.output == Path::new("-") {
            return Err(RsomicsError::ConfigError(
                "index output must be a named file".to_owned(),
            ));
        }
        for input in jobs.iter().map(|job| &job.input) {
            if same_target(input, &job.output)? {
                return Err(RsomicsError::ConfigError(format!(
                    "index output cannot overwrite an alignment input: {}",
                    job.output.display()
                )));
            }
        }
        for other in &jobs[..index] {
            if same_target(&other.output, &job.output)? {
                return Err(RsomicsError::ConfigError(format!(
                    "multiple inputs resolve to the same index output: {}",
                    job.output.display()
                )));
            }
        }
    }
    Ok(())
}
