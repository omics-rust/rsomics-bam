use std::io;
use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::output::{TransactionalFile, same_target};
use crate::{Program, ampliconclip};

#[derive(Debug, Args)]
#[command(after_help = "\
Examples:
  rsomics-bam ampliconclip -b primers.bed input.bam -o clipped.bam
  rsomics-bam ampliconclip -b primers.bed --hard-clip --both-ends input.bam > clipped.bam")]
pub(crate) struct Arguments {
    /// Primer regions to remove
    #[arg(short = 'b', value_name = "BED")]
    bed: PathBuf,

    /// Write BAM to FILE instead of standard output
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Write run statistics to FILE instead of standard error
    #[arg(short = 'f', long = "stats", value_name = "FILE")]
    stats: Option<PathBuf>,

    /// Write BAM without DEFLATE compression
    #[arg(short = 'u', long = "uncompressed")]
    uncompressed: bool,

    /// Soft-clip primer sequence
    #[arg(long, conflicts_with = "hard_clip")]
    soft_clip: bool,

    /// Hard-clip primer sequence and qualities
    #[arg(long)]
    hard_clip: bool,

    /// Clip matching primers at both ends
    #[arg(long)]
    both_ends: bool,

    /// Match the read direction to BED strand
    #[arg(long)]
    strand: bool,

    /// Write only reads that matched a primer
    #[arg(long)]
    clipped: bool,

    /// Mark unmatched mapped reads as QC fail
    #[arg(long)]
    fail: bool,

    /// Reject reads at or below this active query length
    #[arg(long, value_name = "INT", value_parser = parse_length)]
    filter_len: Option<i64>,

    /// Mark reads at or below this active query length as QC fail
    #[arg(long, value_name = "INT", value_parser = parse_length)]
    fail_len: Option<i64>,

    /// Unmap reads at or below this active query length
    #[arg(long, value_name = "INT", default_value_t = 0, value_parser = parse_length)]
    unmap_len: i64,

    /// Reject input reads already marked unmapped or QC fail
    #[arg(long)]
    no_excluded: bool,

    /// Write rejected reads to BAM FILE
    #[arg(long, value_name = "FILE")]
    rejects_file: Option<PathBuf>,

    /// Write per-primer bedGraph counts to FILE
    #[arg(long, value_name = "FILE")]
    primer_counts: Option<PathBuf>,

    /// Add an OA tag containing the original alignment
    #[arg(long)]
    original: bool,

    /// Keep NM and MD tags after clipping
    #[arg(long)]
    keep_tag: bool,

    /// Primer-coordinate matching tolerance
    #[arg(long, value_name = "INT", default_value_t = 5)]
    tolerance: i64,

    /// Do not add an @PG line
    #[arg(long = "no-pg", visible_alias = "no-PG")]
    suppress_program_record: bool,

    /// Additional BAM decompression and compression workers
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,

    /// Coordinate-ordered BAM input
    #[arg(value_name = "BAM")]
    input: PathBuf,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let _ = arguments.soft_clip;
    let output = named(arguments.output.as_deref());
    if json && output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for ampliconclip".to_owned(),
        ));
    }
    let targets = [
        output,
        arguments.stats.as_deref(),
        arguments.rejects_file.as_deref(),
        arguments.primer_counts.as_deref(),
    ];
    validate_targets(&arguments.input, &arguments.bed, &targets)?;

    let output_transaction = output.map(TransactionalFile::new).transpose()?;
    let rejects_transaction = arguments
        .rejects_file
        .as_deref()
        .map(TransactionalFile::new)
        .transpose()?;
    let stats_transaction = arguments
        .stats
        .as_deref()
        .map(TransactionalFile::new)
        .transpose()?;
    let counts_transaction = arguments
        .primer_counts
        .as_deref()
        .map(TransactionalFile::new)
        .transpose()?;

    let output_writer: Box<dyn io::Write + Send> = match &output_transaction {
        Some(transaction) => Box::new(transaction.reopen()?),
        None => Box::new(io::stdout()),
    };
    let reject_writer: Option<Box<dyn io::Write + Send>> = rejects_transaction
        .as_ref()
        .map(|transaction| transaction.reopen().map(|file| Box::new(file) as _))
        .transpose()?;
    let stats_command_line = crate::program::command_line();
    let command_line = (!arguments.suppress_program_record).then(|| stats_command_line.clone());
    let program = command_line
        .as_deref()
        .map(|line| Program::new("rsomics-bam", env!("CARGO_PKG_VERSION"), line))
        .transpose()?;
    let run = ampliconclip::write(
        &arguments.input,
        &arguments.bed,
        ampliconclip::Options {
            mode: if arguments.hard_clip {
                ampliconclip::ClipMode::Hard
            } else {
                ampliconclip::ClipMode::Soft
            },
            both_ends: arguments.both_ends,
            use_strand: arguments.strand,
            tolerance: arguments.tolerance,
            mark_fail: arguments.fail,
            clipped_only: arguments.clipped,
            exclude_flagged: arguments.no_excluded,
            filter_length: arguments.filter_len,
            fail_length: arguments.fail_len,
            unmap_length: (arguments.unmap_len >= 0).then_some(arguments.unmap_len),
            keep_tags: arguments.keep_tag,
            original: arguments.original,
            uncompressed: arguments.uncompressed,
            additional_threads: arguments.threads,
            program,
        },
        output_writer,
        reject_writer,
    )?;

    if let Some(transaction) = &stats_transaction {
        run.summary
            .write(transaction.reopen()?, Some(&stats_command_line))?;
    } else {
        run.summary
            .write(io::stderr().lock(), Some(&stats_command_line))?;
    }
    if let Some(transaction) = &counts_transaction {
        ampliconclip::write_primer_counts(&run.primer_counts, transaction.reopen()?)?;
    }

    if let Some(transaction) = rejects_transaction {
        transaction.commit()?;
    }
    if let Some(transaction) = stats_transaction {
        transaction.commit()?;
    }
    if let Some(transaction) = counts_transaction {
        transaction.commit()?;
    }
    if let Some(transaction) = output_transaction {
        transaction.commit()?;
    }
    Ok(CommandOutput::Ampliconclip { run })
}

fn named(path: Option<&Path>) -> Option<&Path> {
    path.filter(|path| *path != Path::new("-"))
}

fn validate_targets(input: &Path, bed: &Path, targets: &[Option<&Path>]) -> Result<()> {
    for (index, target) in targets.iter().enumerate() {
        let Some(target) = target else { continue };
        if (input != Path::new("-") && same_target(input, target)?) || same_target(bed, target)? {
            return Err(RsomicsError::ConfigError(
                "ampliconclip outputs must differ from its inputs".to_owned(),
            ));
        }
        for other in targets.iter().skip(index + 1).flatten() {
            if same_target(target, other)? {
                return Err(RsomicsError::ConfigError(
                    "ampliconclip output paths must be distinct".to_owned(),
                ));
            }
        }
    }
    Ok(())
}

fn parse_length(value: &str) -> std::result::Result<i64, String> {
    value
        .parse::<i64>()
        .ok()
        .filter(|value| *value >= -1)
        .ok_or_else(|| "length must be -1 or greater".to_owned())
}
