use std::io;
use std::path::PathBuf;

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::head;

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Input SAM, BAM, or CRAM file; omit or use - for standard input
    #[arg(value_name = "ALIGNMENT", default_value = "-")]
    input: PathBuf,

    /// Print at most this many header lines
    #[arg(short = 'H', long, value_name = "INT")]
    headers: Option<usize>,

    /// Print at most this many alignment records
    #[arg(short = 'n', long, value_name = "INT", default_value_t = 0)]
    records: usize,

    /// Reference FASTA for CRAM decoding
    #[arg(short = 'T', long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Additional alignment decompression threads
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if json {
        return Err(RsomicsError::ConfigError(
            "--json cannot be combined with SAM stream output".to_owned(),
        ));
    }

    let summary = head::write(
        &arguments.input,
        head::Options {
            header_lines: arguments.headers,
            records: arguments.records,
            reference: arguments.reference.as_deref(),
            additional_threads: arguments.threads,
        },
        io::stdout().lock(),
    )?;

    Ok(CommandOutput::Head { summary })
}
