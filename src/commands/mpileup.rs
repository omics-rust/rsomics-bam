use std::io;
use std::path::PathBuf;

use clap::{ArgAction, Args};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::{flags, mpileup, output::TransactionalFile};

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Coordinate-sorted SAM, BAM, or CRAM file
    #[arg(value_name = "ALIGNMENT")]
    input: PathBuf,

    /// Indexed reference FASTA
    #[arg(
        short = 'f',
        long = "fasta-ref",
        visible_alias = "reference",
        value_name = "FASTA"
    )]
    reference: Option<PathBuf>,

    /// Write output to a file
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Minimum base quality
    #[arg(short = 'Q', long = "min-BQ", value_name = "INT", default_value_t = 13)]
    minimum_base_quality: u8,

    /// Minimum mapping quality
    #[arg(short = 'q', long = "min-MQ", value_name = "INT", default_value_t = 0)]
    minimum_mapping_quality: u8,

    /// Maximum per-input depth; 0 disables the limit
    #[arg(short = 'd', long, value_name = "INT", default_value_t = 8000)]
    maximum_depth: usize,

    /// Disable BAQ
    #[arg(short = 'B', long = "no-BAQ", conflicts_with = "recalculate_baq")]
    disable_baq: bool,

    /// Recalculate BAQ and ignore existing adjustments
    #[arg(short = 'E', long = "redo-BAQ", conflicts_with = "disable_baq")]
    recalculate_baq: bool,

    /// Include anomalous read pairs
    #[arg(short = 'A', long = "count-orphans")]
    include_anomalous_pairs: bool,

    /// Disable overlapping-mate quality adjustment
    #[arg(short = 'x', long = "ignore-overlaps-removal")]
    ignore_overlaps: bool,

    /// Exclude records with any FLAG bits
    #[arg(long = "ff", visible_alias = "excl-flags", value_name = "FLAG", default_value = "UNMAP,SECONDARY,QCFAIL,DUP", value_parser = parse_flags)]
    excluded_flags: u16,

    /// Require all FLAG bits
    #[arg(long = "rf", visible_alias = "incl-flags", value_name = "FLAG", default_value = "0", value_parser = parse_flags)]
    required_flags: u16,

    /// Output uncovered positions; repeat for unused references
    #[arg(short = 'a', action = ArgAction::Count)]
    all_positions: u8,

    /// Additional alignment decompression threads
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if json && arguments.output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires --output for mpileup".to_owned(),
        ));
    }
    let options = mpileup::Options {
        reference: arguments.reference.as_deref(),
        additional_threads: arguments.threads,
        minimum_base_quality: arguments.minimum_base_quality,
        minimum_mapping_quality: arguments.minimum_mapping_quality,
        maximum_depth: arguments.maximum_depth,
        adjust_overlaps: !arguments.ignore_overlaps,
        include_anomalous_pairs: arguments.include_anomalous_pairs,
        excluded_flags: arguments.excluded_flags,
        required_flags: arguments.required_flags,
        baq: if arguments.disable_baq {
            mpileup::BaqMode::Disabled
        } else if arguments.recalculate_baq {
            mpileup::BaqMode::Recalculate
        } else {
            mpileup::BaqMode::Calculate
        },
        positions: match arguments.all_positions {
            0 => mpileup::PositionMode::Covered,
            1 => mpileup::PositionMode::UsedReferences,
            _ => mpileup::PositionMode::AllReferences,
        },
    };
    let summary = if let Some(path) = arguments.output.as_deref() {
        let transaction = TransactionalFile::new(path)?;
        let summary = mpileup::write(&arguments.input, options, transaction.reopen()?)?;
        transaction.commit()?;
        summary
    } else {
        mpileup::write(&arguments.input, options, io::stdout().lock())?
    };
    Ok(CommandOutput::Mpileup { summary })
}

fn parse_flags(value: &str) -> std::result::Result<u16, String> {
    flags::parse(value).map_err(|error| error.to_string())
}
