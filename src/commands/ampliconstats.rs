use std::io;
use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::output::{TransactionalFile, same_target};
use crate::{ampliconstats, flags};

#[derive(Debug, Args)]
#[command(after_help = "\
Example:
  rsomics-bam ampliconstats primers.bed clipped.bam -o amplicons.txt")]
pub(crate) struct Arguments {
    /// Include reads containing every listed SAM flag
    #[arg(short = 'f', long = "required-flag", value_name = "FLAG", default_value = "0", value_parser = parse_flags)]
    required_flag: u16,

    /// Exclude reads containing any listed SAM flag
    #[arg(short = 'F', long = "filter-flag", value_name = "FLAG", default_value = "UNMAP,SECONDARY,QCFAIL,SUPPLEMENTARY", value_parser = parse_flags)]
    filter_flag: u16,

    /// Maximum number of amplicons per reference
    #[arg(short = 'a', long, value_name = "INT", default_value_t = 1000)]
    max_amplicons: usize,

    /// Maximum amplicon length
    #[arg(short = 'l', long, value_name = "INT", default_value_t = 1000)]
    max_amplicon_length: usize,

    /// Coverage depth thresholds, comma-separated
    #[arg(short = 'd', long, value_name = "INT[,INT]", default_value = "1")]
    min_depth: String,

    /// Primer-position matching margin
    #[arg(short = 'm', long, value_name = "INT", default_value_t = 30)]
    pos_margin: i64,

    /// Write statistics to FILE instead of standard output
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Use the first read-group sample name for each input
    #[arg(short = 's', long)]
    use_sample_name: bool,

    /// Adjust template length after clipping without fixmate
    #[arg(short = 't', long, value_name = "INT", default_value_t = 0)]
    tlen_adjust: i64,

    /// Bin template coordinates by this width
    #[arg(short = 'b', long, value_name = "INT", default_value_t = 1)]
    tcoord_bin: i64,

    /// Minimum template-coordinate count to report
    #[arg(short = 'c', long, value_name = "INT", default_value_t = 10)]
    tcoord_min_count: u32,

    /// Merge neighbouring depths within this fraction
    #[arg(short = 'D', long, value_name = "FRACTION", default_value_t = 0.01)]
    depth_bin: f64,

    /// Use the legacy single-reference output schema
    #[arg(short = 'S', long)]
    single_ref: bool,

    /// Additional BAM decompression workers
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,

    /// Six-column primer BED
    #[arg(value_name = "PRIMERS.BED")]
    primers: PathBuf,

    /// Coordinate-ordered clipped BAM inputs
    #[arg(value_name = "BAM", required = true)]
    inputs: Vec<PathBuf>,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let output = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"));
    if json && output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for ampliconstats".to_owned(),
        ));
    }
    if arguments.inputs.len() > 1 && arguments.inputs.iter().any(|path| path == Path::new("-")) {
        return Err(RsomicsError::ConfigError(
            "standard input cannot be combined with other BAM inputs".to_owned(),
        ));
    }
    if let Some(output) = output {
        if same_target(&arguments.primers, output)? {
            return Err(RsomicsError::ConfigError(
                "ampliconstats output must differ from the primer BED".to_owned(),
            ));
        }
        for input in &arguments.inputs {
            if input != Path::new("-") && same_target(input, output)? {
                return Err(RsomicsError::ConfigError(
                    "ampliconstats output must differ from every BAM input".to_owned(),
                ));
            }
        }
    }

    let transaction = output.map(TransactionalFile::new).transpose()?;
    let input_paths: Vec<&Path> = arguments.inputs.iter().map(PathBuf::as_path).collect();
    let options = ampliconstats::Options {
        flag_require: arguments.required_flag,
        flag_filter: arguments.filter_flag,
        max_delta: arguments.pos_margin,
        min_depth: parse_depths(&arguments.min_depth)?,
        max_amp: arguments.max_amplicons,
        max_amp_len: arguments.max_amplicon_length,
        tlen_adj: arguments.tlen_adjust,
        depth_bin: arguments.depth_bin,
        tcoord_min_count: arguments.tcoord_min_count,
        tcoord_bin: arguments.tcoord_bin,
        additional_threads: arguments.threads,
        single_ref: arguments.single_ref,
        use_sample_name: arguments.use_sample_name,
    };
    let command_line = crate::program::command_line();
    let summary = match &transaction {
        Some(transaction) => ampliconstats::write(
            &options,
            &arguments.primers,
            &input_paths,
            &command_line,
            &mut io::BufWriter::new(transaction.reopen()?),
        )?,
        None => ampliconstats::write(
            &options,
            &arguments.primers,
            &input_paths,
            &command_line,
            &mut io::BufWriter::new(io::stdout().lock()),
        )?,
    };
    if let Some(transaction) = transaction {
        transaction.commit()?;
    }
    Ok(CommandOutput::Ampliconstats { summary })
}

fn parse_flags(value: &str) -> std::result::Result<u16, String> {
    flags::parse(value).map_err(|error| error.to_string())
}

fn parse_depths(value: &str) -> Result<[u32; ampliconstats::MAX_DEPTH_LEVELS]> {
    let mut depths = [0; ampliconstats::MAX_DEPTH_LEVELS];
    let parts: Vec<_> = value.split(',').collect();
    if parts.is_empty() || parts.len() > depths.len() {
        return Err(RsomicsError::ConfigError(
            "--min-depth accepts one to five values".to_owned(),
        ));
    }
    for (depth, part) in depths.iter_mut().zip(parts) {
        *depth = part.trim().parse().map_err(|_| {
            RsomicsError::ConfigError(format!("invalid --min-depth value: {part:?}"))
        })?;
    }
    Ok(depths)
}
