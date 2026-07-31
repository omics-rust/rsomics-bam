use std::fs::File;
use std::io::{self, BufRead, BufWriter, IsTerminal, Write};
use std::path::PathBuf;
use std::process;

use clap::{ArgAction, Args, Parser, Subcommand, ValueEnum};
use rsomics_common::{OutputArgs, Result, RsomicsError, ToolMeta, run as run_tool};
use serde::Serialize;

use crate::{flags, flagstat, head, quickcheck, samples};

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
    Flags(FlagsArgs),
    /// Count alignments by SAM flag category
    Flagstat(FlagstatArgs),
    /// Print header lines and the first alignments as SAM
    Head(HeadArgs),
    /// Check alignment headers and format-specific end markers
    Quickcheck(QuickcheckArgs),
    /// List samples declared by alignment read groups
    Samples(SamplesArgs),
}

#[derive(Debug, Args)]
struct FlagsArgs {
    /// Numeric or comma-separated symbolic flag values
    #[arg(value_name = "FLAGS")]
    values: Vec<String>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
enum FlagstatFormat {
    #[default]
    Text,
    Json,
    Tsv,
}

#[derive(Debug, Args)]
struct FlagstatArgs {
    /// Input SAM, BAM, or CRAM file
    #[arg(value_name = "ALIGNMENT")]
    input: PathBuf,

    /// Output representation
    #[arg(short = 'O', long = "output-fmt", value_enum, default_value_t)]
    format: FlagstatFormat,

    /// Reference FASTA for CRAM decoding
    #[arg(short = 'T', long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Additional alignment decompression threads
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,
}

#[derive(Debug, Args)]
struct HeadArgs {
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

#[derive(Debug, Args)]
struct QuickcheckArgs {
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

#[derive(Debug, Args)]
struct SamplesArgs {
    /// Add column names before the results
    #[arg(short = 'H', long)]
    header: bool,

    /// Report whether each alignment has an index
    #[arg(short = 'i', long)]
    index: bool,

    /// Read-group tag to report
    #[arg(short = 'T', long, default_value = "SM")]
    tag: String,

    /// Write tabular output to a file
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Reference FASTA to match against alignment targets
    #[arg(short = 'f', long = "reference", value_name = "FASTA")]
    references: Vec<PathBuf>,

    /// File containing reference FASTA paths
    #[arg(short = 'F', long = "reference-list", value_name = "FILE")]
    reference_lists: Vec<PathBuf>,

    /// Pair each alignment with a custom index
    #[arg(short = 'X', long)]
    custom_index: bool,

    /// Input SAM, BAM, or CRAM files; read paths from stdin when omitted
    #[arg(value_name = "ALIGNMENT")]
    inputs: Vec<PathBuf>,
}

#[derive(Debug, Serialize)]
#[serde(tag = "command", rename_all = "kebab-case")]
enum CommandOutput {
    Flags { values: Vec<flags::FlagValue> },
    Flagstat { counts: Box<flagstat::Counts> },
    Head { summary: head::Summary },
    Quickcheck { report: quickcheck::Report },
    Samples { report: samples::Report },
}

#[must_use]
pub(crate) fn run() -> process::ExitCode {
    let cli = rsomics_help::parse::<Cli>();
    let output = cli.output.clone();
    run_tool(&output, META, || execute(cli))
}

fn execute(cli: Cli) -> Result<CommandOutput> {
    match cli.command {
        Command::Flags(arguments) => execute_flags(arguments, cli.output.json),
        Command::Flagstat(arguments) => execute_flagstat(arguments, cli.output.json),
        Command::Head(arguments) => execute_head(arguments, cli.output.json),
        Command::Quickcheck(arguments) => execute_quickcheck(arguments, cli.output.json),
        Command::Samples(arguments) => execute_samples(arguments, cli.output.json),
    }
}

fn execute_flags(arguments: FlagsArgs, json: bool) -> Result<CommandOutput> {
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

fn execute_flagstat(arguments: FlagstatArgs, json: bool) -> Result<CommandOutput> {
    if json && arguments.format != FlagstatFormat::Text {
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
            FlagstatFormat::Text => write!(output, "{counts}"),
            FlagstatFormat::Json => {
                serde_json::to_writer_pretty(&mut output, &counts.to_json())
                    .map_err(|error| io::Error::other(error.to_string()))?;
                writeln!(output)
            }
            FlagstatFormat::Tsv => write!(output, "{}", counts.to_tsv()),
        }
        .map_err(RsomicsError::Io)?;
        output.flush().map_err(RsomicsError::Io)?;
    }

    Ok(CommandOutput::Flagstat {
        counts: Box::new(counts),
    })
}

fn execute_head(arguments: HeadArgs, json: bool) -> Result<CommandOutput> {
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

fn execute_quickcheck(arguments: QuickcheckArgs, json: bool) -> Result<CommandOutput> {
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

fn execute_samples(arguments: SamplesArgs, json: bool) -> Result<CommandOutput> {
    if json && arguments.output.is_some() {
        return Err(RsomicsError::ConfigError(
            "--json cannot be combined with --output".to_owned(),
        ));
    }

    let tag = samples::Tag::parse(&arguments.tag)?;
    let inputs = samples_inputs(arguments.inputs, arguments.custom_index)?;
    let references = reference_paths(arguments.references, &arguments.reference_lists)?;
    let report = samples::collect(
        &inputs,
        samples::Options {
            tag,
            test_index: arguments.index,
            references: &references,
        },
    )?;

    if !json {
        match arguments.output {
            Some(path) => {
                let file = File::create(&path).map_err(|error| {
                    RsomicsError::Io(io::Error::new(
                        error.kind(),
                        format!("creating {}: {error}", path.display()),
                    ))
                })?;
                let mut output = BufWriter::new(file);
                report.write(arguments.header, &mut output)?;
                output.flush().map_err(RsomicsError::Io)?;
            }
            None => report.write(arguments.header, io::stdout().lock())?,
        }
    }

    Ok(CommandOutput::Samples { report })
}

fn samples_inputs(inputs: Vec<PathBuf>, custom_index: bool) -> Result<Vec<samples::Input>> {
    if inputs.is_empty() {
        if io::stdin().is_terminal() {
            return Err(RsomicsError::ConfigError(
                "provide alignment paths or pipe one path per line to stdin".to_owned(),
            ));
        }
        let lines = io::stdin()
            .lock()
            .lines()
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(RsomicsError::Io)?;
        return lines
            .into_iter()
            .map(|line| {
                if custom_index {
                    let (alignment, index) = line.split_once('\t').ok_or_else(|| {
                        RsomicsError::InvalidInput(format!(
                            "expected ALIGNMENT<TAB>INDEX, got {line:?}"
                        ))
                    })?;
                    if alignment.is_empty() || index.is_empty() {
                        return Err(RsomicsError::InvalidInput(format!(
                            "expected ALIGNMENT<TAB>INDEX, got {line:?}"
                        )));
                    }
                    Ok(samples::Input {
                        path: PathBuf::from(alignment),
                        index: Some(PathBuf::from(index)),
                    })
                } else {
                    Ok(samples::Input {
                        path: PathBuf::from(line),
                        index: None,
                    })
                }
            })
            .collect();
    }

    if !custom_index {
        return Ok(inputs
            .into_iter()
            .map(|path| samples::Input { path, index: None })
            .collect());
    }

    if !inputs.len().is_multiple_of(2) {
        return Err(RsomicsError::ConfigError(
            "custom indexes require one index for every alignment".to_owned(),
        ));
    }
    let split = inputs.len() / 2;
    let (alignments, indexes) = inputs.split_at(split);
    Ok(alignments
        .iter()
        .cloned()
        .zip(indexes.iter().cloned())
        .map(|(path, index)| samples::Input {
            path,
            index: Some(index),
        })
        .collect())
}

fn reference_paths(mut references: Vec<PathBuf>, lists: &[PathBuf]) -> Result<Vec<PathBuf>> {
    for list in lists {
        let file = File::open(list).map_err(|error| {
            RsomicsError::Io(io::Error::new(
                error.kind(),
                format!("opening reference list {}: {error}", list.display()),
            ))
        })?;
        for line in io::BufReader::new(file).lines() {
            references.push(PathBuf::from(line.map_err(RsomicsError::Io)?));
        }
    }
    Ok(references)
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
}
