use std::io;
use std::path::{Path, PathBuf};

use clap::{Args, ValueEnum};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::output::{TransactionalFile, same_target};
use crate::{Program, depad};

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum Format {
    Sam,
    Bam,
}

#[derive(Debug, Args)]
#[command(after_help = "\
Examples:
  rsomics-bam depad padded.bam -o unpadded.bam
  rsomics-bam depad -T padded.fa -s padded.bam > unpadded.sam
  rsomics-bam depad -T padded.fa -1 -@ 4 -o unpadded.bam padded.cram")]
pub(crate) struct Arguments {
    /// Write SAM output
    #[arg(short = 's', long = "sam", conflicts_with = "format")]
    sam: bool,

    /// Accept the legacy SAM-input flag; input format is detected automatically
    #[arg(short = 'S')]
    legacy_sam_input: bool,

    /// Write BAM without DEFLATE compression
    #[arg(short = 'u', long = "uncompressed", conflicts_with = "fast")]
    uncompressed: bool,

    /// Write BAM with level-1 compression
    #[arg(
        short = '1',
        long = "fast-compression",
        conflicts_with = "uncompressed"
    )]
    fast: bool,

    /// Select SAM or BAM output
    #[arg(short = 'O', long = "output-fmt", value_enum, conflicts_with = "sam")]
    format: Option<Format>,

    /// Padded FASTA containing * or - gap columns
    #[arg(short = 'T', long = "reference", value_name = "FILE")]
    reference: Option<PathBuf>,

    /// Write output to FILE instead of standard output
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Additional alignment I/O workers
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line
    #[arg(long = "no-pg", visible_alias = "no-PG")]
    suppress_program_record: bool,

    /// Input SAM, BAM, or no-reference CRAM file; use - for standard input
    #[arg(value_name = "ALIGNMENT")]
    input: PathBuf,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let _ = arguments.legacy_sam_input;
    let output = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"));
    if json && output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for depad".to_owned(),
        ));
    }
    if let Some(output) = output {
        if arguments.input != Path::new("-") && same_target(&arguments.input, output)? {
            return Err(RsomicsError::ConfigError(
                "depad input and output must be different files".to_owned(),
            ));
        }
        if let Some(reference) = arguments.reference.as_deref()
            && same_target(reference, output)?
        {
            return Err(RsomicsError::ConfigError(
                "depad reference and output must be different files".to_owned(),
            ));
        }
    }

    let format = output_format(
        arguments.format,
        arguments.sam,
        arguments.uncompressed,
        arguments.fast,
        arguments.output.as_deref(),
    )?;
    let command_line = (!arguments.suppress_program_record).then(crate::program::command_line);
    let program = command_line
        .as_deref()
        .map(|line| Program::new("rsomics-bam", env!("CARGO_PKG_VERSION"), line))
        .transpose()?;
    let options = depad::Options {
        format,
        compression: if arguments.uncompressed {
            depad::Compression::Uncompressed
        } else if arguments.fast {
            depad::Compression::Fast
        } else {
            depad::Compression::Default
        },
        reference: arguments.reference.as_deref(),
        additional_threads: arguments.threads,
        program,
    };

    let summary = match output {
        Some(path) => {
            let transaction = TransactionalFile::new(path)?;
            let summary = depad::write(&arguments.input, options, transaction.reopen()?)?;
            transaction.commit()?;
            summary
        }
        None => depad::write(&arguments.input, options, io::stdout())?,
    };
    if arguments.reference.is_none() {
        eprintln!("warning: reference lengths remain padded without --reference");
    }
    if summary.records_with_reference_skips > 0 {
        eprintln!(
            "warning: CIGAR N was treated as D in {} record(s)",
            summary.records_with_reference_skips
        );
    }
    Ok(CommandOutput::Depad { summary })
}

fn output_format(
    format: Option<Format>,
    sam: bool,
    uncompressed: bool,
    fast: bool,
    output: Option<&Path>,
) -> Result<depad::Format> {
    if (sam || format == Some(Format::Sam)) && (uncompressed || fast) {
        return Err(RsomicsError::ConfigError(
            "BAM compression options cannot be used with SAM output".to_owned(),
        ));
    }
    match format {
        Some(Format::Sam) => Ok(depad::Format::Sam),
        Some(Format::Bam) => Ok(depad::Format::Bam),
        None if sam => Ok(depad::Format::Sam),
        None if output.is_some_and(is_sam_path) && !uncompressed && !fast => Ok(depad::Format::Sam),
        None if output.is_some_and(is_cram_path) => Err(RsomicsError::ConfigError(
            "CRAM output is not available; select -O sam or -O bam".to_owned(),
        )),
        None => Ok(depad::Format::Bam),
    }
}

fn is_sam_path(path: &Path) -> bool {
    path.extension()
        .and_then(std::ffi::OsStr::to_str)
        .is_some_and(|extension| extension.eq_ignore_ascii_case("sam"))
}

fn is_cram_path(path: &Path) -> bool {
    path.extension()
        .and_then(std::ffi::OsStr::to_str)
        .is_some_and(|extension| extension.eq_ignore_ascii_case("cram"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_selection_defaults_to_bam_and_rejects_cram() {
        assert_eq!(
            output_format(None, false, false, false, None).unwrap(),
            depad::Format::Bam
        );
        assert_eq!(
            output_format(None, false, false, false, Some(Path::new("out.sam"))).unwrap(),
            depad::Format::Sam
        );
        assert!(output_format(None, false, false, false, Some(Path::new("out.cram"))).is_err());
        assert!(output_format(Some(Format::Sam), false, true, false, None).is_err());
        assert!(output_format(None, true, false, true, None).is_err());
    }
}
