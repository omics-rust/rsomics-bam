use std::io;
use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::output::TransactionalFile;
use crate::{Program, hts_quickcheck, reheader};

#[derive(Debug, Args)]
#[command(after_help = "\
Examples:
  rsomics-bam reheader replacement.sam input.bam -o output.bam")]
pub(crate) struct Arguments {
    /// Replacement header from a SAM, BAM, or CRAM file
    #[arg(value_name = "HEADER")]
    header: PathBuf,

    /// Input BAM whose alignment records are preserved
    #[arg(value_name = "BAM")]
    input: PathBuf,

    /// Write BAM output to FILE; omit or use - for standard output
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Do not add an @PG line
    #[arg(long = "no-pg", visible_alias = "no-PG")]
    suppress_program_record: bool,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let output = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"));
    if json && output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for reheader".to_owned(),
        ));
    }

    let command_line = (!arguments.suppress_program_record).then(crate::program::command_line);
    let program = command_line
        .as_deref()
        .map(|line| Program::new("rsomics-bam", env!("CARGO_PKG_VERSION"), line))
        .transpose()?;
    let options = reheader::Options {
        destination: output,
        program,
    };
    let summary = match output {
        Some(path) => {
            let transaction = TransactionalFile::new(path)?;
            let summary = reheader::write(
                &arguments.header,
                &arguments.input,
                options,
                transaction.reopen()?,
            )?;
            hts_quickcheck::require_bgzf_eof(transaction.temporary_path())?;
            transaction.commit()?;
            summary
        }
        None => reheader::write(&arguments.header, &arguments.input, options, io::stdout())?,
    };
    Ok(CommandOutput::Reheader { summary })
}
