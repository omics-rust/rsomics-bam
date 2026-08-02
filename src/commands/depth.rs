use std::fs;
use std::io;
use std::path::{Path, PathBuf};

use clap::{ArgAction, Args};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::{depth, flags, output::TransactionalFile};

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Coordinate-sorted SAM, BAM, or CRAM files
    #[arg(value_name = "ALIGNMENT", required_unless_present = "input_list")]
    inputs: Vec<PathBuf>,

    /// Read alignment paths from a file
    #[arg(
        short = 'f',
        long = "input-list",
        value_name = "FILE",
        conflicts_with = "inputs"
    )]
    input_list: Option<PathBuf>,

    /// Write output to a file
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Output zero-depth positions; repeat to include unused references
    #[arg(short = 'a', action = ArgAction::Count)]
    all_positions: u8,

    /// Restrict output to BED intervals
    #[arg(short = 'b', long, value_name = "BED")]
    bed: Option<PathBuf>,

    /// Print column names
    #[arg(short = 'H', long)]
    header: bool,

    /// Ignore reads shorter than this aligned query length
    #[arg(short = 'l', long, value_name = "INT", default_value_t = 0)]
    minimum_read_length: usize,

    /// Minimum base quality
    #[arg(short = 'q', long = "min-BQ", value_name = "INT", default_value_t = 0)]
    minimum_base_quality: u8,

    /// Minimum mapping quality
    #[arg(short = 'Q', long = "min-MQ", value_name = "INT", default_value_t = 0)]
    minimum_mapping_quality: u8,

    /// Indexed region
    #[arg(short = 'r', long, value_name = "REGION")]
    region: Option<String>,

    /// Remove flags from the default exclusion set
    #[arg(short = 'g', value_name = "FLAG", action = ArgAction::Append, value_parser = parse_flags)]
    restored_flags: Vec<u16>,

    /// Add flags to the exclusion set
    #[arg(short = 'G', long = "excl-flags", value_name = "FLAG", action = ArgAction::Append, value_parser = parse_flags)]
    excluded_flags: Vec<u16>,

    /// Include records with at least one listed flag
    #[arg(long = "incl-flags", value_name = "FLAG", action = ArgAction::Append, value_parser = parse_flags)]
    included_flags: Vec<u16>,

    /// Include records with every listed flag
    #[arg(long = "require-flags", value_name = "FLAG", action = ArgAction::Append, value_parser = parse_flags)]
    required_flags: Vec<u16>,

    /// Count deletions
    #[arg(short = 'J', long)]
    include_deletions: bool,

    /// Count only the first read across mate overlaps
    #[arg(short = 's', long)]
    remove_overlaps: bool,

    /// Reference FASTA for CRAM input
    #[arg(long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Additional alignment decompression threads
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if json && arguments.output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires --output for depth".to_owned(),
        ));
    }
    let inputs = if let Some(path) = arguments.input_list.as_deref() {
        read_input_list(path)?
    } else {
        arguments.inputs
    };
    if inputs.is_empty() {
        return Err(RsomicsError::ConfigError(
            "depth input list is empty".to_owned(),
        ));
    }
    let restored = combine_flags(&arguments.restored_flags);
    let excluded = (0x704 & !restored) | combine_flags(&arguments.excluded_flags);
    let options = depth::Options {
        reference: arguments.reference.as_deref(),
        additional_threads: arguments.threads,
        minimum_base_quality: arguments.minimum_base_quality,
        minimum_mapping_quality: arguments.minimum_mapping_quality,
        minimum_read_length: arguments.minimum_read_length,
        excluded_flags: excluded,
        included_flags: combine_flags(&arguments.included_flags),
        required_flags: combine_flags(&arguments.required_flags),
        include_deletions: arguments.include_deletions,
        remove_overlaps: arguments.remove_overlaps,
        positions: match arguments.all_positions {
            0 => depth::PositionMode::Covered,
            1 => depth::PositionMode::UsedReferences,
            _ => depth::PositionMode::AllReferences,
        },
        region: arguments.region.as_deref(),
        bed: arguments.bed.as_deref(),
        header: arguments.header,
    };
    let summary = if let Some(path) = arguments.output.as_deref() {
        let transaction = TransactionalFile::new(path)?;
        let summary = depth::write(&inputs, options, transaction.reopen()?)?;
        transaction.commit()?;
        summary
    } else {
        depth::write(&inputs, options, io::stdout().lock())?
    };
    Ok(CommandOutput::Depth { summary })
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

fn combine_flags(values: &[u16]) -> u16 {
    values.iter().copied().fold(0, |all, value| all | value)
}
