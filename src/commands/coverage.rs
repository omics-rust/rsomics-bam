use std::fs;
use std::io;
use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::{
    coverage, flags,
    output::{TransactionalFile, same_target},
};

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Input SAM, BAM, or CRAM files
    #[arg(value_name = "ALIGNMENT", required_unless_present = "input_list")]
    inputs: Vec<PathBuf>,

    /// Read alignment paths from a file
    #[arg(
        short = 'b',
        long = "bam-list",
        value_name = "FILE",
        conflicts_with = "inputs"
    )]
    input_list: Option<PathBuf>,

    /// Write tabular output to a file
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Do not print column names
    #[arg(short = 'H', long = "no-header")]
    no_header: bool,

    /// Ignore reads shorter than this query length
    #[arg(
        short = 'l',
        long = "min-read-len",
        value_name = "INT",
        default_value_t = 0
    )]
    minimum_read_length: usize,

    /// Minimum mapping quality
    #[arg(
        short = 'q',
        long = "min-MQ",
        visible_alias = "min-mq",
        value_name = "INT",
        default_value_t = 0
    )]
    minimum_mapping_quality: u8,

    /// Minimum base quality
    #[arg(
        short = 'Q',
        long = "min-BQ",
        visible_alias = "min-bq",
        value_name = "INT",
        default_value_t = 0
    )]
    minimum_base_quality: u8,

    /// Include records with at least one listed flag
    #[arg(long = "rf", visible_aliases = ["incl-flags", "include-flags"], value_name = "FLAG", default_value = "0", value_parser = parse_flags)]
    required_flags: u16,

    /// Exclude records with any listed flag
    #[arg(long = "ff", visible_aliases = ["excl-flags", "exclude-flags"], value_name = "FLAG", default_value = "UNMAP,SECONDARY,QCFAIL,DUP", value_parser = parse_flags)]
    excluded_flags: u16,

    /// Maximum pileup depth per input; 0 removes the limit
    #[arg(
        short = 'd',
        long = "depth",
        value_name = "INT",
        default_value_t = 1_000_000
    )]
    maximum_depth: usize,

    /// Ignore positions below this combined depth
    #[arg(long = "min-depth", value_name = "INT", default_value_t = 1, value_parser = parse_positive_usize)]
    minimum_depth: usize,

    /// Restrict coverage to an indexed region
    #[arg(short = 'r', long, value_name = "REGION")]
    region: Option<String>,

    /// Reference FASTA for CRAM decoding
    #[arg(short = 'T', long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Additional alignment decompression threads
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let inputs = if let Some(path) = arguments.input_list.as_deref() {
        read_input_list(path)?
    } else {
        arguments.inputs
    };
    if inputs.is_empty() {
        return Err(RsomicsError::ConfigError(
            "coverage input list is empty".to_owned(),
        ));
    }
    if json
        && arguments
            .output
            .as_deref()
            .is_none_or(|path| path == Path::new("-"))
    {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for coverage".to_owned(),
        ));
    }
    if let Some(output) = arguments.output.as_deref()
        && output != Path::new("-")
    {
        for input in &inputs {
            if input != Path::new("-") && same_target(input, output)? {
                return Err(RsomicsError::ConfigError(
                    "coverage output must differ from every alignment input".to_owned(),
                ));
            }
        }
        if let Some(list) = arguments.input_list.as_deref()
            && same_target(list, output)?
        {
            return Err(RsomicsError::ConfigError(
                "coverage output must differ from the input list".to_owned(),
            ));
        }
    }
    let report = coverage::collect(
        &inputs,
        coverage::Options {
            reference: arguments.reference.as_deref(),
            additional_threads: arguments.threads,
            minimum_read_length: arguments.minimum_read_length,
            minimum_mapping_quality: arguments.minimum_mapping_quality,
            minimum_base_quality: arguments.minimum_base_quality,
            required_flags: arguments.required_flags,
            excluded_flags: arguments.excluded_flags,
            maximum_depth: arguments.maximum_depth,
            minimum_depth: arguments.minimum_depth,
            region: arguments.region.as_deref(),
        },
    )?;
    if let Some(path) = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"))
    {
        let transaction = TransactionalFile::new(path)?;
        report.write(!arguments.no_header, transaction.reopen()?)?;
        transaction.commit()?;
    } else {
        report.write(!arguments.no_header, io::stdout().lock())?;
    }
    Ok(CommandOutput::Coverage { report })
}

fn read_input_list(path: &Path) -> Result<Vec<PathBuf>> {
    let content = fs::read_to_string(path).map_err(RsomicsError::Io)?;
    Ok(content
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(PathBuf::from)
        .collect())
}

fn parse_flags(value: &str) -> std::result::Result<u16, String> {
    flags::parse(value).map_err(|error| error.to_string())
}

fn parse_positive_usize(value: &str) -> std::result::Result<usize, String> {
    value
        .parse::<usize>()
        .ok()
        .filter(|&value| value > 0)
        .ok_or_else(|| "value must be greater than zero".to_owned())
}
