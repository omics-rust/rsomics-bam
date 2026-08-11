use std::path::{Path, PathBuf};

use clap::{ArgGroup, Args, ValueEnum};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::{Program, split};

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum Format {
    Sam,
    Bam,
    Cram,
}

#[derive(Debug, Args)]
#[command(
    group(ArgGroup::new("mode").args(["tag", "parts", "genes", "mates"]).multiple(false)),
    after_help = "\
Examples:
  rsomics-bam split reads.bam -b sample
  rsomics-bam split reads.bam --tag NM -b by-edit-distance
  rsomics-bam split reads.bam --parts 4 --seed 7 -b shard
  rsomics-bam split reads.bam --genes genes.bed12 -b genes
  rsomics-bam split reads.bam --mates -b mates"
)]
pub(crate) struct Arguments {
    /// Partition by an auxiliary tag instead of header read groups
    #[arg(long, value_name = "TAG")]
    tag: Option<String>,

    /// Partition every retained record into one of N deterministic parts
    #[arg(long, value_name = "N")]
    parts: Option<usize>,

    /// Partition by the leftmost alignment start and BED12 exons
    #[arg(long, value_name = "BED12")]
    genes: Option<PathBuf>,

    /// Emit R1, R2, and unmapped single-end projections
    #[arg(long)]
    mates: bool,

    /// Prefix for every partition output
    #[arg(short = 'b', long, value_name = "PREFIX", required = true)]
    output_prefix: PathBuf,

    /// Route missing, incompatible, or excess tag values to FILE
    #[arg(short = 'u', long, value_name = "FILE")]
    unaccounted: Option<PathBuf>,

    /// Read the unaccounted output header from ALIGNMENT
    #[arg(long, value_name = "ALIGNMENT", requires = "unaccounted")]
    unaccounted_header: Option<PathBuf>,

    /// Select SAM, BAM, or CRAM output
    #[arg(
        short = 'O',
        long = "output-fmt",
        value_enum,
        ignore_case = true,
        default_value = "bam"
    )]
    format: Format,

    /// Maximum number of partition outputs
    #[arg(
        short = 'M',
        long = "max-outputs",
        value_name = "N",
        default_value_t = 100
    )]
    maximum_outputs: usize,

    /// Minimum digits for part numbers and integer tag values
    #[arg(long, value_name = "N", default_value_t = 0)]
    zero_pad: usize,

    /// Deterministic 64-bit seed for --parts; defaults to 0
    #[arg(long, value_name = "INT", requires = "parts")]
    seed: Option<u64>,

    /// Exclude records carrying the unmapped flag from --parts
    #[arg(long, requires = "parts")]
    skip_unmapped: bool,

    /// Reference FASTA for CRAM input or output
    #[arg(short = 'T', long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Additional alignment I/O workers
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line
    #[arg(long = "no-pg", visible_alias = "no-PG")]
    suppress_program_record: bool,

    /// Input SAM, BAM, or CRAM file; omit or use - for standard input
    #[arg(value_name = "ALIGNMENT")]
    input: Option<PathBuf>,
}

pub(crate) fn execute(arguments: Arguments, _json: bool) -> Result<CommandOutput> {
    let input = arguments.input.as_deref().unwrap_or_else(|| Path::new("-"));
    let mode = if let Some(tag) = arguments.tag.as_deref() {
        split::Mode::Tag(parse_tag(tag)?)
    } else if let Some(count) = arguments.parts {
        split::Mode::Parts {
            count,
            seed: arguments.seed.unwrap_or(0),
            skip_unmapped: arguments.skip_unmapped,
        }
    } else if let Some(path) = arguments.genes.as_deref() {
        split::Mode::Genes(path)
    } else if arguments.mates {
        split::Mode::Mates
    } else {
        split::Mode::ReadGroup
    };
    let command_line = (!arguments.suppress_program_record).then(crate::program::command_line);
    let program = command_line
        .as_deref()
        .map(|line| Program::new("rsomics-bam", env!("CARGO_PKG_VERSION"), line))
        .transpose()?;
    let summary = split::run(
        input,
        split::Options {
            mode,
            output_prefix: &arguments.output_prefix,
            unaccounted: arguments.unaccounted.as_deref(),
            unaccounted_header: arguments.unaccounted_header.as_deref(),
            format: match arguments.format {
                Format::Sam => split::Format::Sam,
                Format::Bam => split::Format::Bam,
                Format::Cram => split::Format::Cram,
            },
            maximum_outputs: arguments.maximum_outputs,
            zero_pad: arguments.zero_pad,
            reference: arguments.reference.as_deref(),
            additional_threads: arguments.threads,
            program,
        },
    )?;
    Ok(CommandOutput::Split { summary })
}

fn parse_tag(value: &str) -> Result<[u8; 2]> {
    let tag: [u8; 2] = value.as_bytes().try_into().map_err(|_| {
        RsomicsError::InvalidInput(format!(
            "split tag must contain exactly two bytes: {value:?}"
        ))
    })?;
    if !tag[0].is_ascii_alphabetic() || !tag[1].is_ascii_alphanumeric() {
        return Err(RsomicsError::InvalidInput(format!(
            "split tag must match [A-Za-z][A-Za-z0-9]: {value:?}"
        )));
    }
    Ok(tag)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tags_are_exactly_two_bytes() {
        assert_eq!(parse_tag("RG").unwrap(), *b"RG");
        for value in ["", "R", "RGA", "é", "1A", "A_"] {
            assert!(parse_tag(value).is_err());
        }
    }
}
