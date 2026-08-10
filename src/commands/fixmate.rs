use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::output::{TransactionalFile, same_target};
use crate::{Program, fixmate, hts_quickcheck};

#[derive(Debug, Args)]
#[command(after_help = "\
Examples:
  rsomics-bam fixmate grouped.bam -o fixed.bam
  rsomics-bam fixmate -m grouped.bam -o fixed.bam")]
pub(crate) struct Arguments {
    /// Input name-grouped SAM, BAM, or CRAM file; omit or use - for standard input
    #[arg(value_name = "ALIGNMENT")]
    input: Option<PathBuf>,

    /// Write BAM output to FILE; omit or use - for standard output
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Add mate-score tags used by markdup
    #[arg(short = 'm', long = "mate-score")]
    mate_score: bool,

    /// Remove unmapped and secondary records
    #[arg(short = 'r', long)]
    remove: bool,

    /// Disable the FR proper-pair check
    #[arg(short = 'p', long = "no-proper-pair")]
    no_proper_pair: bool,

    /// Reference FASTA for CRAM decoding
    #[arg(long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Additional BAM I/O workers; omit for automatic parallelism
    #[arg(short = '@', long, value_name = "INT")]
    threads: Option<usize>,

    /// Do not add an @PG line
    #[arg(long = "no-pg", visible_alias = "no-PG")]
    suppress_program_record: bool,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let input = arguments.input.as_deref().unwrap_or_else(|| Path::new("-"));
    let output = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"));
    if json && output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for fixmate".to_owned(),
        ));
    }
    if let Some(output) = output
        && input != Path::new("-")
        && same_target(input, output)?
    {
        return Err(RsomicsError::ConfigError(
            "fixmate input and output must be different files".to_owned(),
        ));
    }

    let command_line = (!arguments.suppress_program_record).then(crate::program::command_line);
    let program = command_line
        .as_deref()
        .map(|line| Program::new("rsomics-bam", env!("CARGO_PKG_VERSION"), line))
        .transpose()?;
    let options = fixmate::Options {
        mate_score: arguments.mate_score,
        remove: arguments.remove,
        proper_pair_check: !arguments.no_proper_pair,
        additional_threads: arguments.threads,
        reference: arguments.reference.as_deref(),
        destination: output,
        program,
    };

    let summary = match output {
        Some(path) => {
            let transaction = TransactionalFile::new(path)?;
            let summary = fixmate::write(input, options, BufWriter::new(transaction.reopen()?))?;
            hts_quickcheck::require_bgzf_eof(transaction.temporary_path())?;
            transaction.commit()?;
            summary
        }
        None => fixmate::write(input, options, io::stdout())?,
    };
    Ok(CommandOutput::Fixmate { summary })
}
