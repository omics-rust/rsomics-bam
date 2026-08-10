use std::io;
use std::path::PathBuf;

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::{
    idxstats,
    output::{TransactionalFile, same_target},
};

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Input alignment and, with -X, its custom index
    #[arg(value_name = "ALIGNMENT_OR_INDEX", required = true, num_args = 1..=2)]
    paths: Vec<PathBuf>,

    /// Read the second positional path as a custom index
    #[arg(short = 'X')]
    custom_index: bool,

    /// Write tabular output to a file
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Reference FASTA for CRAM decoding
    #[arg(short = 'T', long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Additional alignment decompression threads
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let (input, index) = match (arguments.custom_index, arguments.paths.as_slice()) {
        (false, [input]) => (input.as_path(), None),
        (true, [input, index]) => (input.as_path(), Some(index.as_path())),
        (false, _) => {
            return Err(RsomicsError::ConfigError(
                "a custom index path requires -X".to_owned(),
            ));
        }
        (true, _) => {
            return Err(RsomicsError::ConfigError(
                "-X requires an alignment and custom index path".to_owned(),
            ));
        }
    };
    if let Some(output) = arguments.output.as_deref() {
        if same_target(input, output)? {
            return Err(RsomicsError::ConfigError(
                "idxstats output must differ from the alignment input".to_owned(),
            ));
        }
        if let Some(index) = index
            && same_target(index, output)?
        {
            return Err(RsomicsError::ConfigError(
                "idxstats output must differ from the custom index".to_owned(),
            ));
        }
    }
    let report = idxstats::collect(
        input,
        idxstats::Options {
            reference: arguments.reference.as_deref(),
            index,
            additional_threads: arguments.threads,
        },
    )?;
    if let Some(path) = arguments.output.as_deref() {
        let transaction = TransactionalFile::new(path)?;
        report.write(transaction.reopen()?)?;
        transaction.commit()?;
    } else if !json {
        report.write(io::stdout().lock())?;
    }
    Ok(CommandOutput::Idxstats { report })
}
