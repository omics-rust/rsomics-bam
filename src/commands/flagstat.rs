use std::io::{self, Write};
use std::path::PathBuf;

use clap::{Args, ValueEnum};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::flagstat;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
enum Format {
    #[default]
    Text,
    Json,
    Tsv,
}

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Input SAM, BAM, or CRAM file
    #[arg(value_name = "ALIGNMENT")]
    input: PathBuf,

    /// Output representation
    #[arg(short = 'O', long = "output-fmt", value_enum, default_value_t)]
    format: Format,

    /// Reference FASTA for CRAM decoding
    #[arg(short = 'T', long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Additional alignment decompression threads
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if json && arguments.format != Format::Text {
        return Err(RsomicsError::ConfigError(
            "--json cannot be combined with --output-fmt".to_owned(),
        ));
    }

    let counts = flagstat::count(
        &arguments.input,
        flagstat::Options {
            reference: arguments.reference.as_deref(),
            additional_threads: arguments.threads,
        },
    )?;

    if !json {
        let mut output = io::stdout().lock();
        match arguments.format {
            Format::Text => write!(output, "{counts}"),
            Format::Json => {
                serde_json::to_writer_pretty(&mut output, &counts.to_json())
                    .map_err(|error| io::Error::other(error.to_string()))?;
                writeln!(output)
            }
            Format::Tsv => write!(output, "{}", counts.to_tsv()),
        }
        .map_err(RsomicsError::Io)?;
        output.flush().map_err(RsomicsError::Io)?;
    }

    Ok(CommandOutput::Flagstat {
        counts: Box::new(counts),
    })
}
