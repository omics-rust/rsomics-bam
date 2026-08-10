use std::io;
use std::path::PathBuf;

use clap::{ArgAction, Args};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::{
    bedcov, flags,
    output::{TransactionalFile, same_target},
};

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Input BED regions
    #[arg(value_name = "BED")]
    bed: PathBuf,

    /// Indexed SAM, BAM, or CRAM files
    #[arg(value_name = "ALIGNMENT", required = true)]
    paths: Vec<PathBuf>,

    /// Treat the second half of alignment paths as custom indices
    #[arg(short = 'X')]
    custom_indices: bool,

    /// Write tabular output to a file
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Minimum mapping quality
    #[arg(
        short = 'Q',
        long = "min-MQ",
        visible_alias = "min-mq",
        value_name = "INT",
        default_value_t = 0
    )]
    minimum_mapping_quality: u8,

    /// Remove flags from the default exclusion set
    #[arg(short = 'g', value_name = "FLAG", action = ArgAction::Append, value_parser = parse_flags)]
    restored_flags: Vec<u16>,

    /// Add flags to the exclusion set
    #[arg(short = 'G', value_name = "FLAG", action = ArgAction::Append, value_parser = parse_flags)]
    excluded_flags: Vec<u16>,

    /// Exclude deletions and reference skips from coverage
    #[arg(short = 'j')]
    skip_deletions_and_skips: bool,

    /// Add the number of bases meeting this depth
    #[arg(short = 'd', value_name = "INT")]
    depth_threshold: Option<usize>,

    /// Maximum pileup depth; 0 removes the limit
    #[arg(long, value_name = "INT", default_value_t = i32::MAX as usize)]
    max_depth: usize,

    /// Add the number of overlapping reads
    #[arg(short = 'c')]
    read_count: bool,

    /// Print column names
    #[arg(short = 'H')]
    header: bool,

    /// Reference FASTA for CRAM decoding
    #[arg(short = 'T', long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Additional alignment decompression threads
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let (inputs, indices) = if arguments.custom_indices {
        if !arguments.paths.len().is_multiple_of(2) {
            return Err(RsomicsError::ConfigError(
                "-X requires one custom index for every alignment input".to_owned(),
            ));
        }
        let split = arguments.paths.len() / 2;
        (
            arguments.paths[..split].to_vec(),
            Some(arguments.paths[split..].to_vec()),
        )
    } else {
        (arguments.paths, None)
    };
    if inputs.is_empty() {
        return Err(RsomicsError::ConfigError(
            "bedcov requires at least one alignment input".to_owned(),
        ));
    }
    if json
        && arguments
            .output
            .as_deref()
            .is_none_or(|path| path == std::path::Path::new("-"))
    {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for bedcov".to_owned(),
        ));
    }
    if let Some(output) = arguments.output.as_deref()
        && output != std::path::Path::new("-")
    {
        if same_target(&arguments.bed, output)? {
            return Err(RsomicsError::ConfigError(
                "bedcov output must differ from the BED input".to_owned(),
            ));
        }
        for input in inputs.iter().chain(indices.iter().flatten()) {
            if same_target(input, output)? {
                return Err(RsomicsError::ConfigError(
                    "bedcov output must differ from alignment and index inputs".to_owned(),
                ));
            }
        }
    }
    let restored_flags = combine_flags(&arguments.restored_flags);
    let options = bedcov::Options {
        reference: arguments.reference.as_deref(),
        indices: indices.as_deref(),
        additional_threads: arguments.threads,
        minimum_mapping_quality: arguments.minimum_mapping_quality,
        excluded_flags: (0x704 & !restored_flags) | combine_flags(&arguments.excluded_flags),
        skip_deletions_and_skips: arguments.skip_deletions_and_skips,
        depth_threshold: arguments.depth_threshold,
        maximum_depth: arguments.max_depth,
        read_count: arguments.read_count,
        header: arguments.header,
    };
    let summary = if let Some(path) = arguments
        .output
        .as_deref()
        .filter(|path| *path != std::path::Path::new("-"))
    {
        let transaction = TransactionalFile::new(path)?;
        let summary = bedcov::write(&arguments.bed, &inputs, options, transaction.reopen()?)?;
        transaction.commit()?;
        summary
    } else {
        bedcov::write(&arguments.bed, &inputs, options, io::stdout().lock())?
    };
    Ok(CommandOutput::Bedcov { summary })
}

fn parse_flags(value: &str) -> std::result::Result<u16, String> {
    flags::parse(value).map_err(|error| error.to_string())
}

fn combine_flags(values: &[u16]) -> u16 {
    values.iter().copied().fold(0, |all, value| all | value)
}
