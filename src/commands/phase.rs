use std::io;
use std::path::{Path, PathBuf};

use clap::{Args, ValueEnum};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::output::{TransactionalFile, same_target};
use crate::phase;

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
enum Format {
    Sam,
    #[default]
    Bam,
    Cram,
}

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Coordinate-sorted SAM, BAM, or CRAM; omit or use - for standard input
    #[arg(value_name = "ALIGNMENT")]
    input: Option<PathBuf>,

    /// Write the phase report to FILE
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Prefix for haplotype and chimera alignment outputs
    #[arg(short = 'b', long = "output-prefix", value_name = "PREFIX")]
    output_prefix: Option<PathBuf>,

    /// Alignment output format used with --output-prefix
    #[arg(short = 'O', long = "output-fmt", value_enum, default_value_t)]
    output_format: Format,

    /// Maximum local phasing window
    #[arg(short = 'k', value_name = "INT", default_value_t = 13, value_parser = parse_window)]
    window: usize,

    /// Minimum Phred-scaled LOD for a heterozygous call
    #[arg(short = 'q', value_name = "INT", default_value_t = 37)]
    minimum_lod: i32,

    /// Minimum base quality used in heterozygous calling
    #[arg(
        short = 'Q',
        long = "min-BQ",
        visible_alias = "min-bq",
        value_name = "INT",
        default_value_t = 13
    )]
    minimum_base_quality: u8,

    /// Maximum raw pileup depth
    #[arg(short = 'D', value_name = "INT", default_value_t = 256, value_parser = parse_depth)]
    maximum_depth: usize,

    /// Disable chimeric-read repair
    #[arg(short = 'F')]
    no_chimera_repair: bool,

    /// Route reads with ambiguous phase to the chimera output
    #[arg(short = 'A')]
    drop_ambiguous: bool,

    /// Reference FASTA for CRAM input or output
    #[arg(long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Do not add an @PG line to partition headers
    #[arg(long = "no-pg", visible_alias = "no-PG")]
    suppress_program_record: bool,

    /// Additional alignment I/O workers
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,
}

fn parse_window(value: &str) -> std::result::Result<usize, String> {
    let value = value
        .parse::<usize>()
        .map_err(|_| format!("invalid phase window: {value}"))?;
    (1..=phase::MAX_WINDOW)
        .contains(&value)
        .then_some(value)
        .ok_or_else(|| format!("phase window must be between 1 and {}", phase::MAX_WINDOW))
}

fn parse_depth(value: &str) -> std::result::Result<usize, String> {
    let value = value
        .parse::<usize>()
        .map_err(|_| format!("invalid maximum depth: {value}"))?;
    (1..=phase::MAX_DEPTH)
        .contains(&value)
        .then_some(value)
        .ok_or_else(|| {
            format!(
                "maximum phase depth must be between 1 and {}",
                phase::MAX_DEPTH
            )
        })
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if arguments.drop_ambiguous && arguments.output_prefix.is_none() {
        return Err(RsomicsError::ConfigError(
            "-A requires --output-prefix".to_owned(),
        ));
    }
    if json && arguments.output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires --output for phase".to_owned(),
        ));
    }
    if arguments.output_prefix.is_some()
        && matches!(arguments.output_format, Format::Cram)
        && arguments.reference.is_none()
    {
        return Err(RsomicsError::ConfigError(
            "CRAM partition output requires --reference".to_owned(),
        ));
    }
    let input = arguments.input.as_deref().unwrap_or_else(|| Path::new("-"));
    if let Some(output) = arguments.output.as_deref()
        && input != Path::new("-")
        && same_target(input, output)?
    {
        return Err(RsomicsError::ConfigError(
            "phase output cannot overwrite the alignment input".to_owned(),
        ));
    }
    let options = phase::Options {
        window: arguments.window,
        minimum_lod: arguments.minimum_lod,
        minimum_base_quality: arguments.minimum_base_quality,
        maximum_depth: arguments.maximum_depth,
        fix_chimeras: !arguments.no_chimera_repair,
        reference: arguments.reference.clone(),
        additional_threads: arguments.threads,
    };
    let summary = if let Some(prefix) = arguments.output_prefix.as_deref() {
        let format = match arguments.output_format {
            Format::Sam => phase::PartitionFormat::Sam,
            Format::Bam => phase::PartitionFormat::Bam,
            Format::Cram => phase::PartitionFormat::Cram,
        };
        let targets = phase::partition_paths(prefix, format);
        validate_targets(input, arguments.output.as_deref(), &targets)?;
        let mut transactions = Vec::with_capacity(4);
        let report_index = if let Some(output) = arguments.output.as_deref() {
            transactions.push(TransactionalFile::new(output)?);
            Some(0)
        } else {
            None
        };
        let partition_start = transactions.len();
        for target in &targets {
            transactions.push(TransactionalFile::new(target)?);
        }
        let files: Vec<_> = transactions[partition_start..]
            .iter()
            .map(TransactionalFile::reopen)
            .collect::<Result<_>>()?;
        let files: [std::fs::File; 3] = files.try_into().unwrap();
        let stdout = io::stdout();
        let mut report: Box<dyn io::Write> = if let Some(index) = report_index {
            Box::new(transactions[index].reopen()?)
        } else {
            Box::new(stdout.lock())
        };
        let command_line = (!arguments.suppress_program_record).then(crate::program::command_line);
        let program = command_line
            .as_deref()
            .map(|line| crate::Program::new("rsomics-bam", env!("CARGO_PKG_VERSION"), line))
            .transpose()?;
        let summary = phase::write_partitioned(
            input,
            options,
            &mut report,
            files,
            format,
            arguments.drop_ambiguous,
            program,
        )?;
        drop(report);
        TransactionalFile::commit_all(transactions)?;
        summary
    } else if let Some(output) = arguments.output.as_deref() {
        let transaction = TransactionalFile::new(output)?;
        let summary = phase::write(input, options, transaction.reopen()?)?;
        transaction.commit()?;
        summary
    } else {
        phase::write(input, options, io::stdout().lock())?
    };
    Ok(CommandOutput::Phase { summary })
}

fn validate_targets(input: &Path, report: Option<&Path>, partitions: &[PathBuf; 3]) -> Result<()> {
    for target in partitions {
        if input != Path::new("-") && same_target(input, target)? {
            return Err(RsomicsError::ConfigError(
                "phase partition output cannot overwrite the alignment input".to_owned(),
            ));
        }
        if let Some(report) = report
            && same_target(report, target)?
        {
            return Err(RsomicsError::ConfigError(
                "phase report and partition outputs require different files".to_owned(),
            ));
        }
    }
    for left in 0..partitions.len() {
        for right in left + 1..partitions.len() {
            if same_target(&partitions[left], &partitions[right])? {
                return Err(RsomicsError::ConfigError(
                    "phase partition outputs require different files".to_owned(),
                ));
            }
        }
    }
    Ok(())
}
