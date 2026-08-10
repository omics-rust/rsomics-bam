use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};

use clap::{Args, ValueEnum};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::output::{TransactionalFile, same_target};
use crate::{Program, calmd};

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum Format {
    Sam,
    Bam,
}

#[derive(Debug, Args)]
#[command(after_help = "\
Examples:
  rsomics-bam calmd reads.bam reference.fa > reads.sam
  rsomics-bam calmd -b -o reads.calmd.bam reads.bam reference.fa
  rsomics-bam calmd -e -u reads.bam reference.fa > reads.calmd.bam")]
pub(crate) struct Arguments {
    /// Rewrite reference-matching query bases as =
    #[arg(short = 'e', long = "convert-equal")]
    use_equal: bool,

    /// Write compressed BAM output
    #[arg(short = 'b', long = "bam", conflicts_with = "format")]
    bam: bool,

    /// Write BAM without DEFLATE compression
    #[arg(short = 'u', long = "uncompressed")]
    uncompressed: bool,

    /// Select SAM or BAM output
    #[arg(short = 'O', long = "output-fmt", value_enum, conflicts_with = "bam")]
    format: Option<Format>,

    /// Write output to FILE instead of standard output
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Additional alignment I/O workers
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line
    #[arg(long = "no-pg", visible_alias = "no-PG")]
    suppress_program_record: bool,

    /// Input SAM, BAM, or CRAM file; use - for standard input
    #[arg(value_name = "ALIGNMENT")]
    input: PathBuf,

    /// Indexed reference FASTA used to recalculate MD and NM
    #[arg(value_name = "REFERENCE")]
    reference: PathBuf,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let output = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"));
    if json && output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for calmd".to_owned(),
        ));
    }
    if let Some(output) = output {
        for (label, input) in [
            ("input", arguments.input.as_path()),
            ("reference", arguments.reference.as_path()),
        ] {
            if input != Path::new("-") && same_target(input, output)? {
                return Err(RsomicsError::ConfigError(format!(
                    "calmd {label} and output must be different files"
                )));
            }
        }
    }

    let format = output_format(
        arguments.format,
        arguments.bam,
        arguments.uncompressed,
        arguments.output.as_deref(),
    )?;
    let command_line = (!arguments.suppress_program_record).then(crate::program::command_line);
    let program = command_line
        .as_deref()
        .map(|line| Program::new("rsomics-bam", env!("CARGO_PKG_VERSION"), line))
        .transpose()?;
    let options = calmd::Options {
        format,
        compression: if arguments.uncompressed {
            calmd::Compression::Uncompressed
        } else {
            calmd::Compression::Default
        },
        use_equal: arguments.use_equal,
        additional_threads: arguments.threads,
        program,
    };

    let summary = match output {
        Some(path) => {
            let transaction = TransactionalFile::new(path)?;
            let summary = calmd::write(
                &arguments.input,
                &arguments.reference,
                options,
                BufWriter::new(transaction.reopen()?),
            )?;
            transaction.commit()?;
            summary
        }
        None => calmd::write(
            &arguments.input,
            &arguments.reference,
            options,
            io::stdout(),
        )?,
    };
    if summary.records_without_sequence > 0 {
        eprintln!(
            "warning: {} mapped records were preserved because they have no query sequence",
            summary.records_without_sequence
        );
    }
    Ok(CommandOutput::Calmd { summary })
}

fn output_format(
    format: Option<Format>,
    bam: bool,
    uncompressed: bool,
    output: Option<&Path>,
) -> Result<calmd::Format> {
    match format {
        Some(Format::Sam) if uncompressed => Err(RsomicsError::ConfigError(
            "--uncompressed requires BAM output".to_owned(),
        )),
        Some(Format::Sam) => Ok(calmd::Format::Sam),
        Some(Format::Bam) => Ok(calmd::Format::Bam),
        None if bam || uncompressed => Ok(calmd::Format::Bam),
        None if output.is_some_and(is_bam_path) => Ok(calmd::Format::Bam),
        None if output.is_some_and(is_cram_path) => Err(RsomicsError::ConfigError(
            "CRAM output is not available; select -O sam or -O bam".to_owned(),
        )),
        None => Ok(calmd::Format::Sam),
    }
}

fn is_bam_path(path: &Path) -> bool {
    path.extension()
        .and_then(std::ffi::OsStr::to_str)
        .is_some_and(|extension| extension.eq_ignore_ascii_case("bam"))
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
    fn output_selection_rejects_cram_and_sam_compression() {
        assert_eq!(
            output_format(None, false, false, Some(Path::new("reads.BAM"))).unwrap(),
            calmd::Format::Bam
        );
        assert!(output_format(None, false, false, Some(Path::new("reads.cram"))).is_err());
        assert!(output_format(Some(Format::Sam), false, true, None).is_err());
    }
}
