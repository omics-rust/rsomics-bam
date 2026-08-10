use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::output::TransactionalFile;
use crate::{Program, cat, hts_quickcheck};

#[derive(Debug, Args)]
#[command(after_help = "\
Examples:
  rsomics-bam cat lane1.bam lane2.bam -o combined.bam
  rsomics-bam cat --list shards.txt --header header.sam -o combined.bam")]
pub(crate) struct Arguments {
    /// Input BAM files appended after entries from --list
    #[arg(value_name = "BAM")]
    inputs: Vec<PathBuf>,

    /// Read input BAM paths from FILE, one per line; may be repeated
    #[arg(short = 'b', long = "list", value_name = "FILE")]
    lists: Vec<PathBuf>,

    /// Use the alignment header from FILE instead of the first input
    #[arg(long, value_name = "FILE")]
    header: Option<PathBuf>,

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
            "--json requires a named --output for cat".to_owned(),
        ));
    }

    let mut inputs = read_lists(&arguments.lists)?;
    inputs.extend(arguments.inputs);
    let command_line = (!arguments.suppress_program_record).then(crate::program::command_line);
    let program = command_line
        .as_deref()
        .map(|line| Program::new("rsomics-bam", env!("CARGO_PKG_VERSION"), line))
        .transpose()?;
    let options = cat::Options {
        header: arguments.header.as_deref(),
        destination: output,
        program,
    };

    let summary = match output {
        Some(path) => {
            let transaction = TransactionalFile::new(path)?;
            let summary = cat::write(&inputs, options, BufWriter::new(transaction.reopen()?))?;
            hts_quickcheck::require_bgzf_eof(transaction.temporary_path())?;
            transaction.commit()?;
            summary
        }
        None => cat::write(&inputs, options, io::stdout())?,
    };
    Ok(CommandOutput::Cat { summary })
}

fn read_lists(paths: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let mut inputs = Vec::new();
    for path in paths {
        let text = std::fs::read_to_string(path).map_err(|error| {
            RsomicsError::InvalidInput(format!("reading input list {}: {error}", path.display()))
        })?;
        inputs.extend(
            text.lines()
                .map(str::trim)
                .filter(|line| !line.is_empty())
                .map(PathBuf::from),
        );
    }
    Ok(inputs)
}
