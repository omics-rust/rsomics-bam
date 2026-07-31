use std::process;

use clap::{Parser, Subcommand};
use rsomics_common::{OutputArgs, Result, ToolMeta, run as run_tool};
use serde::Serialize;

use crate::{commands, flags, flagstat, head, quickcheck, samples, view};

const META: ToolMeta = ToolMeta {
    name: "rsomics-bam",
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-bam",
    version,
    about = "SAM, BAM, and CRAM workflows",
    arg_required_else_help = true
)]
struct Cli {
    #[command(flatten)]
    output: OutputArgs,

    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Convert between numeric and symbolic SAM flags
    Flags(commands::flags::Arguments),
    /// Count alignments by SAM flag category
    Flagstat(commands::flagstat::Arguments),
    /// Print header lines and the first alignments as SAM
    Head(commands::head::Arguments),
    /// Check alignment headers and format-specific end markers
    Quickcheck(commands::quickcheck::Arguments),
    /// List samples declared by alignment read groups
    Samples(commands::samples::Arguments),
    /// Filter and convert alignment records
    View(commands::view::Arguments),
}

#[derive(Debug, Serialize)]
#[serde(tag = "command", rename_all = "kebab-case")]
pub(crate) enum CommandOutput {
    Flags { values: Vec<flags::FlagValue> },
    Flagstat { counts: Box<flagstat::Counts> },
    Head { summary: head::Summary },
    Quickcheck { report: quickcheck::Report },
    Samples { report: samples::Report },
    View { summary: view::Summary },
}

#[must_use]
pub(crate) fn run() -> process::ExitCode {
    let cli = rsomics_help::parse::<Cli>();
    let output = cli.output.clone();
    run_tool(&output, META, || execute(cli))
}

fn execute(cli: Cli) -> Result<CommandOutput> {
    match cli.command {
        Command::Flags(arguments) => commands::flags::execute(arguments, cli.output.json),
        Command::Flagstat(arguments) => commands::flagstat::execute(arguments, cli.output.json),
        Command::Head(arguments) => commands::head::execute(arguments, cli.output.json),
        Command::Quickcheck(arguments) => commands::quickcheck::execute(arguments, cli.output.json),
        Command::Samples(arguments) => commands::samples::execute(arguments, cli.output.json),
        Command::View(arguments) => commands::view::execute(arguments, cli.output.json),
    }
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    use super::*;

    #[test]
    fn command_tree_is_valid() {
        Cli::command().debug_assert();
    }

    #[test]
    fn help_uses_one_nested_command_tree() {
        let error =
            rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "flagstat", "--help"])
                .unwrap_err();
        let help = error.to_string();
        assert!(help.contains("Input SAM, BAM, or CRAM file"), "{help}");
        assert!(help.contains("--output-fmt"), "{help}");
        assert!(help.contains("-@, --threads <INT>"), "{help}");
    }

    #[test]
    fn view_exposes_program_provenance_control() {
        let command = Cli::command();
        let view = command
            .get_subcommands()
            .find(|command| command.get_name() == "view")
            .unwrap();
        let argument = view
            .get_arguments()
            .find(|argument| argument.get_long() == Some("no-pg"))
            .unwrap();
        assert!(!argument.is_hide_set());
        assert!(
            rsomics_help::try_parse_from::<Cli, _, _>([
                "rsomics-bam",
                "view",
                "--no-PG",
                "input.sam"
            ])
            .is_ok()
        );
    }
}
