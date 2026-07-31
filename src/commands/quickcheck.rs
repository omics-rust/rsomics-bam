use std::io;
use std::path::PathBuf;

use clap::{ArgAction, Args};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::quickcheck;

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// List failed inputs
    #[arg(short = 'v', action = ArgAction::Count)]
    verbose: u8,

    /// Suppress per-input warnings
    #[arg(short = 'q', long)]
    quiet: bool,

    /// Allow headers without reference targets
    #[arg(short = 'u', long)]
    unmapped: bool,

    /// Input SAM, BAM, or CRAM files
    #[arg(value_name = "ALIGNMENT", required = true)]
    inputs: Vec<PathBuf>,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let report = quickcheck::check_all(
        &arguments.inputs,
        quickcheck::Options {
            allow_no_targets: arguments.unmapped,
        },
    );

    if !report.is_ok() {
        if !json {
            report.write_diagnostics(
                arguments.verbose,
                arguments.quiet,
                io::stdout().lock(),
                io::stderr().lock(),
            )?;
        }
        return Err(RsomicsError::InvalidInput(format!(
            "quickcheck failed for {} of {} inputs",
            report.failed(),
            report.files.len()
        )));
    }

    Ok(CommandOutput::Quickcheck { report })
}
