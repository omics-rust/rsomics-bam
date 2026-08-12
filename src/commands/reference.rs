use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::output::{TransactionalFile, same_target};
use crate::reference::{self, Options};

#[derive(Debug, Args)]
#[command(after_help = "\
Examples:
  rsomics-bam reference reads.bam > reference.fa
  rsomics-bam reference -r chr1:1-1000 -o chr1.fa reads.bam
  rsomics-bam reference --embedded reads.cram > reference.fa")]
pub(crate) struct Arguments {
    /// Extract reference blocks embedded in CRAM instead of using MD tags
    #[arg(short = 'e', long)]
    embedded: bool,

    /// Suppress per-reference coverage diagnostics
    #[arg(short = 'q', long)]
    quiet: bool,

    /// Restrict an indexed input to REGION
    #[arg(short = 'r', long, value_name = "REGION")]
    region: Option<String>,

    /// Write FASTA output to FILE instead of standard output
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Additional alignment I/O workers
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,

    /// Input SAM, BAM, or CRAM file; use - for standard input
    #[arg(value_name = "ALIGNMENT", default_value = "-")]
    input: PathBuf,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let output = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"));
    if json && output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for reference".to_owned(),
        ));
    }
    if let Some(output) = output
        && arguments.input != Path::new("-")
        && same_target(&arguments.input, output)?
    {
        return Err(RsomicsError::ConfigError(
            "reference input and output must be different files".to_owned(),
        ));
    }
    if arguments.region.is_some() && arguments.input == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "--region requires a file-backed indexed input".to_owned(),
        ));
    }

    let options = Options {
        embedded: arguments.embedded,
        region: arguments.region,
        additional_threads: arguments.threads,
    };
    let summary = match output {
        Some(path) => {
            let transaction = TransactionalFile::new(path)?;
            let summary = reference::write(
                &arguments.input,
                options,
                BufWriter::new(transaction.reopen()?),
            )?;
            transaction.commit()?;
            summary
        }
        None => reference::write(&arguments.input, options, io::stdout())?,
    };
    if !arguments.quiet {
        for item in &summary.items {
            eprintln!(
                "reference {}: {} bases, {:.2}% recovered",
                item.name, item.bases, item.coverage
            );
        }
    }
    Ok(CommandOutput::Reference { summary })
}
