use std::io::{self, IsTerminal};
use std::path::PathBuf;

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::{output::TransactionalFile, tview};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Display {
    Terminal,
    Text,
    Html,
}

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Indexed SAM, BAM, or CRAM file
    #[arg(value_name = "ALIGNMENT")]
    input: PathBuf,

    /// Display mode: terminal, text, or html
    #[arg(short = 'd', long, value_name = "MODE", default_value = "terminal", value_parser = parse_display)]
    display: Display,

    /// Initial one-based reference position
    #[arg(short = 'p', long, value_name = "REFERENCE[:POSITION]")]
    position: Option<String>,

    /// Show alignments for one sample
    #[arg(short = 's', long, value_name = "SAMPLE")]
    sample: Option<String>,

    /// Text and HTML viewport width
    #[arg(short = 'w', long, value_name = "INT")]
    width: Option<usize>,

    /// Hide insertion columns
    #[arg(short = 'i', long)]
    hide_insertions: bool,

    /// Custom BAI, CSI, or CRAI index
    #[arg(short = 'X', long, value_name = "FILE")]
    index: Option<PathBuf>,

    /// Indexed reference FASTA
    #[arg(short = 'T', long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Write text or HTML to a file
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Additional alignment decompression threads
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if arguments.display == Display::Terminal {
        if arguments.width.is_some() {
            return Err(RsomicsError::ConfigError(
                "--width is only available for text and HTML tview".to_owned(),
            ));
        }
        if json {
            return Err(RsomicsError::ConfigError(
                "--json cannot be combined with terminal tview".to_owned(),
            ));
        }
        if arguments.output.is_some() {
            return Err(RsomicsError::ConfigError(
                "--output cannot be combined with terminal tview".to_owned(),
            ));
        }
        if !io::stdin().is_terminal() || !io::stdout().is_terminal() {
            return Err(RsomicsError::ConfigError(
                "terminal tview requires terminal input and output".to_owned(),
            ));
        }
        let summary = tview::interactive(
            &arguments.input,
            tview::Options {
                reference: arguments.reference.as_deref(),
                index: arguments.index.as_deref(),
                position: arguments.position.as_deref(),
                sample: arguments.sample.as_deref(),
                width: 80,
                hide_insertions: arguments.hide_insertions,
                additional_threads: arguments.threads,
            },
        )?;
        return Ok(CommandOutput::Tview { summary });
    }
    let format = match arguments.display {
        Display::Text => tview::Format::Text,
        Display::Html => tview::Format::Html,
        Display::Terminal => unreachable!(),
    };
    if json && arguments.output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires --output for tview".to_owned(),
        ));
    }
    let width = arguments.width.unwrap_or_else(default_width);
    if width == 0 {
        return Err(RsomicsError::ConfigError(
            "tview width must be greater than zero".to_owned(),
        ));
    }
    let options = tview::Options {
        reference: arguments.reference.as_deref(),
        index: arguments.index.as_deref(),
        position: arguments.position.as_deref(),
        sample: arguments.sample.as_deref(),
        width,
        hide_insertions: arguments.hide_insertions,
        additional_threads: arguments.threads,
    };
    let summary = if let Some(path) = arguments.output.as_deref() {
        let transaction = TransactionalFile::new(path)?;
        let summary = tview::write(&arguments.input, options, format, transaction.reopen()?)?;
        transaction.commit()?;
        summary
    } else {
        tview::write(&arguments.input, options, format, io::stdout().lock())?
    };
    Ok(CommandOutput::Tview { summary })
}

fn parse_display(value: &str) -> std::result::Result<Display, String> {
    match value.to_ascii_lowercase().as_str() {
        "terminal" | "c" => Ok(Display::Terminal),
        "text" | "t" => Ok(Display::Text),
        "html" | "h" => Ok(Display::Html),
        _ => Err(format!(
            "invalid display mode {value:?}; expected terminal, text, or html"
        )),
    }
}

fn default_width() -> usize {
    std::env::var("COLUMNS")
        .ok()
        .and_then(|value| value.parse().ok())
        .filter(|width| *width > 0)
        .unwrap_or(80)
}
