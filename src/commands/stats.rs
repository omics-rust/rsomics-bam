use std::collections::HashSet;
use std::io;
use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::output::{TransactionalFile, same_target, target_identity};
use crate::{flags, stats};

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Coverage bins as MIN,MAX,STEP
    #[arg(short = 'c', long = "coverage", value_name = "MIN,MAX,STEP", default_value = "1,1000,1", value_parser = parse_coverage)]
    coverage: stats::CoverageBins,

    /// Exclude duplicate reads
    #[arg(short = 'd', long = "remove-dups")]
    remove_duplicates: bool,

    /// Use a custom index supplied after the alignment
    #[arg(short = 'X', long = "customized-index")]
    customized_index: bool,

    /// Require all listed SAM flags
    #[arg(short = 'f', long = "required-flag", value_name = "FLAG", default_value = "0", value_parser = parse_flags)]
    required_flags: u16,

    /// Exclude records with any listed SAM flag
    #[arg(short = 'F', long = "filtering-flag", value_name = "FLAG", default_value = "0", value_parser = parse_flags)]
    filtered_flags: u16,

    /// GC-depth window size
    #[arg(
        long = "GC-depth",
        visible_alias = "gc-depth",
        value_name = "FLOAT",
        default_value_t = 20_000.0
    )]
    gc_depth: f64,

    /// Maximum insert-size bin
    #[arg(
        short = 'i',
        long = "insert-size",
        value_name = "INT",
        default_value_t = 8000
    )]
    maximum_insert_size: usize,

    /// Include a sample or read-group ID
    #[arg(short = 'I', long = "id", value_name = "ID")]
    id: Option<String>,

    /// Include only reads of this exact sequence length
    #[arg(short = 'l', long = "read-length", value_name = "INT")]
    read_length: Option<usize>,

    /// Fraction of insert observations used for mean and deviation
    #[arg(
        short = 'm',
        long = "most-inserts",
        value_name = "FLOAT",
        default_value_t = 0.99
    )]
    insert_bulk: f64,

    /// Prefix for split report paths
    #[arg(short = 'P', long = "split-prefix", value_name = "PREFIX")]
    split_prefix: Option<PathBuf>,

    /// BWA-style trimming quality
    #[arg(
        short = 'q',
        long = "trim-quality",
        value_name = "INT",
        default_value_t = 0
    )]
    trim_quality: u8,

    /// Reference FASTA
    #[arg(
        short = 'r',
        long = "ref-seq",
        visible_alias = "reference",
        value_name = "FASTA"
    )]
    reference: Option<PathBuf>,

    /// Accepted for samtools command-line compatibility
    #[arg(short = 's', long = "sam")]
    sam: bool,

    /// Split reports by a two-character auxiliary tag (up to 4096 values)
    #[arg(short = 'S', long = "split", value_name = "TAG")]
    split: Option<String>,

    /// Restrict statistics to target regions
    #[arg(short = 't', long = "target-regions", value_name = "FILE")]
    targets: Option<PathBuf>,

    /// Emit only populated insert-size rows
    #[arg(short = 'x', long = "sparse")]
    sparse: bool,

    /// Remove overlap between paired reads from coverage
    #[arg(short = 'p', long = "remove-overlaps")]
    remove_overlaps: bool,

    /// Target coverage threshold
    #[arg(
        short = 'g',
        long = "cov-threshold",
        value_name = "INT",
        default_value_t = 0
    )]
    coverage_threshold: usize,

    /// Include reference sequence statistics
    #[arg(long = "ref-stats")]
    reference_stats: bool,

    /// Reference-statistics read chunk in MiB
    #[arg(
        long = "ref-stats-chunk",
        value_name = "MB",
        default_value_t = 1,
        allow_hyphen_values = true
    )]
    reference_stats_chunk: i64,

    /// Write the compatibility report to a file
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Additional alignment decompression threads
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,

    /// Input alignment, optional custom index, and indexed regions
    #[arg(value_name = "ALIGNMENT [INDEX] [REGION ...]", num_args = 0..)]
    positional: Vec<String>,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let input = arguments
        .positional
        .first()
        .map_or_else(|| PathBuf::from("-"), PathBuf::from);
    let offset = usize::from(arguments.customized_index);
    if arguments.customized_index && arguments.positional.len() < 2 {
        return Err(RsomicsError::ConfigError(
            "--customized-index requires an index after the alignment".to_owned(),
        ));
    }
    let index = arguments
        .customized_index
        .then(|| PathBuf::from(&arguments.positional[1]));
    let regions = arguments
        .positional
        .get(1 + offset..)
        .unwrap_or_default()
        .to_vec();
    let split_tag = arguments.split.as_deref().map(parse_tag).transpose()?;
    if arguments.split_prefix.is_some() && split_tag.is_none() {
        return Err(RsomicsError::ConfigError(
            "--split-prefix requires --split".to_owned(),
        ));
    }
    let named_output = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"));
    if json && named_output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for stats".to_owned(),
        ));
    }
    for path in [
        Some(input.as_path()),
        index.as_deref(),
        arguments.reference.as_deref(),
        arguments.targets.as_deref(),
    ]
    .into_iter()
    .flatten()
    {
        if let Some(output) = named_output
            && path != Path::new("-")
            && same_target(path, output)?
        {
            return Err(RsomicsError::ConfigError(
                "stats output must differ from every input".to_owned(),
            ));
        }
    }
    let _ = arguments.sam;
    let filtered_flags = arguments.filtered_flags
        | if arguments.remove_duplicates {
            0x400
        } else {
            0
        };
    let report = stats::collect(
        &input,
        stats::Options {
            reference: arguments.reference.as_deref(),
            additional_threads: arguments.threads,
            required_flags: arguments.required_flags,
            filtered_flags,
            read_length: arguments.read_length,
            coverage: arguments.coverage,
            maximum_insert_size: arguments.maximum_insert_size,
            insert_bulk: arguments.insert_bulk,
            trim_quality: arguments.trim_quality,
            gc_depth: arguments.gc_depth,
            sparse: arguments.sparse,
            coverage_threshold: arguments.coverage_threshold,
            targets: arguments.targets.as_deref(),
            regions: &regions,
            index: index.as_deref(),
            id: arguments.id.as_deref(),
            split_tag,
            remove_overlaps: arguments.remove_overlaps,
            reference_stats: arguments.reference_stats,
            reference_stats_chunk: usize::try_from(arguments.reference_stats_chunk.max(1)).unwrap(),
        },
    )?;
    if json {
        serde_json::to_writer(io::sink(), &report).map_err(|error| {
            RsomicsError::ConfigError(format!("serializing stats report: {error}"))
        })?;
    }
    let split_paths = report
        .split_values()
        .map(|value| split_path(&input, arguments.split_prefix.as_deref(), value))
        .collect::<Result<Vec<_>>>()?;
    let mut distinct_split_targets = HashSet::with_capacity(split_paths.len());
    for path in &split_paths {
        for input_path in [
            Some(input.as_path()),
            index.as_deref(),
            arguments.reference.as_deref(),
            arguments.targets.as_deref(),
            named_output,
        ]
        .into_iter()
        .flatten()
        {
            if input_path != Path::new("-") && same_target(input_path, path)? {
                return Err(RsomicsError::ConfigError(format!(
                    "split output aliases another input or output: {}",
                    path.display()
                )));
            }
        }
        if !distinct_split_targets.insert(target_identity(path)?) {
            return Err(RsomicsError::ConfigError(format!(
                "duplicate split output: {}",
                path.display()
            )));
        }
    }
    let mut transactions =
        Vec::with_capacity(split_paths.len() + usize::from(named_output.is_some()));
    for (value, path) in report.split_values().zip(&split_paths) {
        let transaction = TransactionalFile::new(path)?;
        report.write_split(value, transaction.reopen()?)?;
        transactions.push(transaction);
    }
    if let Some(path) = named_output {
        let transaction = TransactionalFile::new(path)?;
        report.write(transaction.reopen()?)?;
        transactions.push(transaction);
        TransactionalFile::commit_all(transactions)?;
    } else {
        report.write(io::stdout().lock())?;
        TransactionalFile::commit_all(transactions)?;
    }
    Ok(CommandOutput::Stats {
        report: Box::new(report),
    })
}

fn parse_tag(value: &str) -> Result<[u8; 2]> {
    value.as_bytes().try_into().map_err(|_| {
        RsomicsError::ConfigError(format!(
            "split tag must contain exactly two bytes: {value:?}"
        ))
    })
}

fn split_path(input: &Path, prefix: Option<&Path>, value: &[u8]) -> Result<PathBuf> {
    let value = std::str::from_utf8(value)
        .map_err(|_| RsomicsError::InvalidInput("split tag value is not UTF-8".to_owned()))?;
    let mut path = prefix.unwrap_or(input).as_os_str().to_os_string();
    path.push("_");
    path.push(value);
    path.push(".bamstat");
    Ok(path.into())
}

fn parse_flags(value: &str) -> std::result::Result<u16, String> {
    flags::parse(value).map_err(|error| error.to_string())
}

fn parse_coverage(value: &str) -> std::result::Result<stats::CoverageBins, String> {
    value.parse()
}
