use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

use clap::{ArgAction, Args};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::cram_size;
use crate::output::{TransactionalFile, same_target};

#[derive(Debug, Args)]
#[command(after_help = "\
Examples:
  rsomics-bam cram-size input.cram
  rsomics-bam cram-size --verbose input.cram
  rsomics-bam cram-size --encodings -o sizes.txt input.cram")]
pub(crate) struct Arguments {
    /// Write the compatibility report to a file
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Report each compression method separately
    #[arg(short = 'v', long, action = ArgAction::Count)]
    verbose: u8,

    /// Include container encoding maps
    #[arg(short = 'e', long, action = ArgAction::Count)]
    encodings: u8,

    /// Input CRAM file, or - for standard input
    #[arg(value_name = "CRAM", default_value = "-")]
    input: PathBuf,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let output = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"));
    if json && output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for cram-size".to_owned(),
        ));
    }
    if let Some(output) = output
        && arguments.input != Path::new("-")
        && same_target(&arguments.input, output)?
    {
        return Err(RsomicsError::ConfigError(
            "cram-size output must differ from its CRAM input".to_owned(),
        ));
    }

    let report = if arguments.input == Path::new("-") {
        cram_size::analyze(io::stdin().lock())?
    } else {
        let file = File::open(&arguments.input).map_err(|error| {
            RsomicsError::Io(io::Error::new(
                error.kind(),
                format!("opening {}: {error}", arguments.input.display()),
            ))
        })?;
        cram_size::analyze(file)?
    };
    let options = cram_size::Options {
        verbose: arguments.verbose != 0,
        encodings: arguments.encodings != 0,
    };
    if let Some(path) = output {
        let transaction = TransactionalFile::new(path)?;
        report.write(options, transaction.reopen()?)?;
        transaction.commit()?;
    } else {
        report.write(options, io::stdout().lock())?;
    }
    Ok(CommandOutput::CramSize { report })
}
