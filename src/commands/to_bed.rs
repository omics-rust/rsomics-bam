use std::io;
use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::output::{TransactionalFile, same_target};
use crate::to_bed::{self, Layout, Options, PairScore, RecordLayout, Score};

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Emit paired alignments as BEDPE
    #[arg(
        long,
        conflicts_with_all = ["bed12", "split", "split_deletions", "tag", "cigar", "color"]
    )]
    bedpe: bool,

    /// Put read1 first in BEDPE output
    #[arg(long, requires = "bedpe")]
    mate1: bool,

    /// Use NM edit distance as the BED score
    #[arg(long = "ed", conflicts_with = "tag")]
    edit_distance: bool,

    /// Use a numeric auxiliary tag as the BED score
    #[arg(long, value_name = "TAG", value_parser = parse_tag)]
    tag: Option<[u8; 2]>,

    /// Append the CIGAR as a seventh BED6 column
    #[arg(long)]
    cigar: bool,

    /// Emit blocked BED12
    #[arg(long, conflicts_with = "cigar")]
    bed12: bool,

    /// Emit blocks separated by CIGAR N as separate BED rows
    #[arg(long, conflicts_with_all = ["cigar", "edit_distance"])]
    split: bool,

    /// Emit blocks separated by CIGAR N or D as separate BED rows
    #[arg(
        long = "split-d",
        conflicts_with_all = ["cigar", "edit_distance"]
    )]
    split_deletions: bool,

    /// BED12 itemRgb value
    #[arg(
        long,
        value_name = "R,G,B",
        requires = "bed12",
        value_parser = parse_color
    )]
    color: Option<[u8; 3]>,

    /// Reference FASTA for CRAM input
    #[arg(short = 'T', long = "reference", value_name = "FILE")]
    reference: Option<PathBuf>,

    /// Additional BAM decompression workers
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,

    /// Write BED output to FILE
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Input SAM, BAM, or CRAM using the bedtools-compatible flag
    #[arg(
        short = 'i',
        long = "input",
        value_name = "ALIGNMENT",
        conflicts_with = "input"
    )]
    flagged_input: Option<PathBuf>,

    /// Input SAM, BAM, or CRAM; use - for standard input
    #[arg(value_name = "ALIGNMENT", required_unless_present = "flagged_input")]
    input: Option<PathBuf>,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let input = arguments
        .flagged_input
        .as_deref()
        .or(arguments.input.as_deref())
        .ok_or_else(|| RsomicsError::ConfigError("to-bed requires an input".to_owned()))?;
    let output = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"));
    if json && output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for to-bed".to_owned(),
        ));
    }
    if let Some(output) = output
        && input != Path::new("-")
        && same_target(input, output)?
    {
        return Err(RsomicsError::ConfigError(
            "to-bed input and output must be different files".to_owned(),
        ));
    }
    if let (Some(reference), Some(output)) = (arguments.reference.as_deref(), output)
        && same_target(reference, output)?
    {
        return Err(RsomicsError::ConfigError(
            "to-bed reference and output must be different files".to_owned(),
        ));
    }
    let score = if arguments.edit_distance {
        Score::EditDistance
    } else if let Some(tag) = arguments.tag {
        Score::Tag(tag)
    } else {
        Score::MappingQuality
    };
    let layout = if arguments.bedpe {
        Layout::Bedpe {
            score: if arguments.edit_distance {
                PairScore::EditDistance
            } else {
                PairScore::MappingQuality
            },
            mate1_first: arguments.mate1,
        }
    } else if arguments.bed12 {
        Layout::Records(RecordLayout::Bed12 {
            score,
            split_deletions: arguments.split_deletions,
            color: arguments.color.unwrap_or([255, 0, 0]),
        })
    } else if arguments.split || arguments.split_deletions {
        Layout::Records(RecordLayout::SplitBed6 {
            score,
            split_deletions: arguments.split_deletions,
        })
    } else {
        Layout::Records(RecordLayout::Bed6 {
            score,
            cigar: arguments.cigar,
        })
    };
    let options = Options {
        reference: arguments.reference.as_deref(),
        additional_threads: arguments.threads,
        layout,
    };
    let summary = if let Some(output) = output {
        let mut transaction = TransactionalFile::new(output)?;
        let summary = to_bed::write(input, options, transaction.file_mut())?;
        transaction.commit()?;
        summary
    } else {
        to_bed::write(input, options, io::stdout().lock())?
    };
    Ok(CommandOutput::ToBed { summary })
}

fn parse_tag(value: &str) -> std::result::Result<[u8; 2], String> {
    let tag: [u8; 2] = value
        .as_bytes()
        .try_into()
        .map_err(|_| format!("tag must contain exactly two bytes: {value:?}"))?;
    if !tag[0].is_ascii_alphabetic() || !tag[1].is_ascii_alphanumeric() {
        return Err(format!("invalid SAM tag: {value:?}"));
    }
    Ok(tag)
}

fn parse_color(value: &str) -> std::result::Result<[u8; 3], String> {
    let channels: Vec<_> = value
        .split(',')
        .map(str::parse::<u8>)
        .collect::<std::result::Result<_, _>>()
        .map_err(|_| format!("color must be an R,G,B triplet from 0 to 255: {value:?}"))?;
    channels
        .try_into()
        .map_err(|_| format!("color must be an R,G,B triplet from 0 to 255: {value:?}"))
}
