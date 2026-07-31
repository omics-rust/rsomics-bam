use std::io::{self, Write};

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::flags;

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Numeric or comma-separated symbolic flag values
    #[arg(value_name = "FLAGS")]
    values: Vec<String>,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if arguments.values.is_empty() {
        if json {
            let values = flags::definitions()
                .iter()
                .map(|definition| flags::describe(definition.bit))
                .collect();
            return Ok(CommandOutput::Flags { values });
        }

        let mut output = io::stdout().lock();
        for definition in flags::definitions() {
            writeln!(
                output,
                "{:#7x} {:5}  {:<13} {}",
                definition.bit, definition.bit, definition.name, definition.description
            )
            .map_err(RsomicsError::Io)?;
        }
        output.flush().map_err(RsomicsError::Io)?;
        return Ok(CommandOutput::Flags { values: Vec::new() });
    }

    let values = arguments
        .values
        .iter()
        .map(|token| flags::parse(token).map(flags::describe))
        .collect::<Result<Vec<_>>>()?;

    if !json {
        let mut output = io::stdout().lock();
        for value in &values {
            writeln!(output, "{}", flags::render(value)).map_err(RsomicsError::Io)?;
        }
        output.flush().map_err(RsomicsError::Io)?;
    }

    Ok(CommandOutput::Flags { values })
}
