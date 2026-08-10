use std::process;

use clap::{Parser, Subcommand};
use rsomics_common::{OutputArgs, Result, ToolMeta, run as run_tool};
use serde::Serialize;

use crate::{
    cat, collate, commands, depth, fixmate, flags, flagstat, head, index, markdup, merge, mpileup,
    quickcheck, reheader, samples, sort, view,
};

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
    /// Concatenate BAM files without reencoding alignment blocks
    Cat(commands::cat::Arguments),
    /// Group alignments by read name with bounded memory
    Collate(commands::collate::Arguments),
    /// Compute read depth at each position
    Depth(commands::depth::Arguments),
    /// Repair mate fields in name-grouped alignments
    Fixmate(commands::fixmate::Arguments),
    /// Convert between numeric and symbolic SAM flags
    Flags(commands::flags::Arguments),
    /// Count alignments by SAM flag category
    Flagstat(commands::flagstat::Arguments),
    /// Print header lines and the first alignments as SAM
    Head(commands::head::Arguments),
    /// Build BAI, CSI, or CRAI random-access indexes
    Index(commands::index::Arguments),
    /// Merge ordered alignment files
    Merge(commands::merge::Arguments),
    /// Mark duplicate alignments in coordinate order
    Markdup(commands::markdup::Arguments),
    /// Generate per-position text pileup
    Mpileup(commands::mpileup::Arguments),
    /// Check alignment headers and format-specific end markers
    Quickcheck(commands::quickcheck::Arguments),
    /// Replace a BAM header without reencoding alignment blocks
    Reheader(commands::reheader::Arguments),
    /// List samples declared by alignment read groups
    Samples(commands::samples::Arguments),
    /// Sort alignments with bounded memory
    Sort(commands::sort::Arguments),
    /// Filter and convert alignment records
    View(Box<commands::view::Arguments>),
}

#[derive(Debug, Serialize)]
#[serde(tag = "command", rename_all = "kebab-case")]
pub(crate) enum CommandOutput {
    Cat { summary: cat::Summary },
    Collate { summary: collate::Summary },
    Depth { summary: depth::Summary },
    Fixmate { summary: fixmate::Summary },
    Flags { values: Vec<flags::FlagValue> },
    Flagstat { counts: Box<flagstat::Counts> },
    Head { summary: head::Summary },
    Index { summaries: Vec<index::Summary> },
    Merge { summary: merge::Summary },
    Markdup { summary: markdup::Summary },
    Mpileup { summary: mpileup::Summary },
    Quickcheck { report: quickcheck::Report },
    Reheader { summary: reheader::Summary },
    Samples { report: samples::Report },
    Sort { summary: sort::Summary },
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
        Command::Cat(arguments) => commands::cat::execute(arguments, cli.output.json),
        Command::Collate(arguments) => commands::collate::execute(arguments, cli.output.json),
        Command::Depth(arguments) => commands::depth::execute(arguments, cli.output.json),
        Command::Fixmate(arguments) => commands::fixmate::execute(arguments, cli.output.json),
        Command::Flags(arguments) => commands::flags::execute(arguments, cli.output.json),
        Command::Flagstat(arguments) => commands::flagstat::execute(arguments, cli.output.json),
        Command::Head(arguments) => commands::head::execute(arguments, cli.output.json),
        Command::Index(arguments) => commands::index::execute(arguments, cli.output.json),
        Command::Merge(arguments) => commands::merge::execute(arguments, cli.output.json),
        Command::Markdup(arguments) => commands::markdup::execute(arguments, cli.output.json),
        Command::Mpileup(arguments) => commands::mpileup::execute(arguments, cli.output.json),
        Command::Quickcheck(arguments) => commands::quickcheck::execute(arguments, cli.output.json),
        Command::Reheader(arguments) => commands::reheader::execute(arguments, cli.output.json),
        Command::Samples(arguments) => commands::samples::execute(arguments, cli.output.json),
        Command::Sort(arguments) => commands::sort::execute(arguments, cli.output.json),
        Command::View(arguments) => commands::view::execute(*arguments, cli.output.json),
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

    #[test]
    fn markdup_help_exposes_stable_scope() {
        let error = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "markdup", "--help"])
            .unwrap_err();
        let help = error.to_string();
        for option in [
            "-r, --remove",
            "-c, --clear",
            "--include-fails",
            "-m, --mode",
            "-l, --max-read-length",
            "-@, --threads",
            "--reference",
            "--no-pg",
        ] {
            assert!(help.contains(option), "missing {option} in {help}");
        }
        for excluded in ["--barcode-tag", "--duplicate-count", "--read-groups"] {
            assert!(!help.contains(excluded), "unexpected {excluded} in {help}");
        }
    }

    #[test]
    fn cat_and_reheader_help_expose_only_the_stable_slice() {
        let cat = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "cat", "--help"])
            .unwrap_err()
            .to_string();
        for option in ["-b, --list", "--header", "-o, --output", "--no-pg"] {
            assert!(cat.contains(option), "missing {option} in {cat}");
        }
        for excluded in ["--region", "--part", "--fast", "--query"] {
            assert!(!cat.contains(excluded), "unexpected {excluded} in {cat}");
        }

        let reheader =
            rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "reheader", "--help"])
                .unwrap_err()
                .to_string();
        assert!(reheader.contains("-o, --output"), "{reheader}");
        assert!(reheader.contains("--no-pg"), "{reheader}");
        assert!(!reheader.contains("--command"), "{reheader}");
        assert!(!reheader.contains("--in-place"), "{reheader}");
    }
}
