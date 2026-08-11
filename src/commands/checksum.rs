use std::io::{self, Write};
use std::path::{Path, PathBuf};

use clap::{ArgAction, Args};
use rsomics_common::{Result, RsomicsError};

use crate::checksum::{self, TagSelection};
use crate::cli::CommandOutput;
use crate::output::{TransactionalFile, same_target};

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Filter records containing any listed SAM flag
    #[arg(short = 'F', long = "exclude-flags", value_name = "FLAG", default_value = "0x900", value_parser = parse_flags)]
    excluded_flags: u16,

    /// Filter records unless all listed SAM flags are present
    #[arg(short = 'f', long = "require-flags", value_name = "FLAG", default_value = "0", value_parser = parse_flags)]
    required_flags: u16,

    /// SAM flag bits included in each checksum
    #[arg(short = 'b', long = "flag-mask", value_name = "FLAG", default_value = "0x0c1", value_parser = parse_flags)]
    flag_mask: u16,

    /// Do not normalize reverse-strand sequence and quality
    #[arg(short = 'c', long = "no-rev-comp")]
    no_reverse_complement: bool,

    /// Ordered tags, or * followed by tags to exclude
    #[arg(short = 't', long, value_name = "TAG[,TAG...]", default_value = "BC,FI,QT,RT,TC", value_parser = parse_tags)]
    tags: TagSelection,

    /// Include record order; repeat for absolute whole-file order
    #[arg(short = 'O', long = "in-order", action = ArgAction::Count)]
    order: u8,

    /// Include reference and position
    #[arg(short = 'P', long = "check-pos")]
    check_position: bool,

    /// Include mapping quality and CIGAR
    #[arg(short = 'C', long = "check-cigar")]
    check_cigar: bool,

    /// Include mate reference, position, and template length
    #[arg(short = 'M', long = "check-mate")]
    check_mate: bool,

    /// Stop after this many accepted records
    #[arg(short = 'N', long = "count", value_name = "INT")]
    maximum_records: Option<u64>,

    /// Normalize selected alignment fields before checksumming
    #[arg(short = 'z', long, value_name = "FLAGS", value_parser = checksum::Sanitize::parse)]
    sanitize: Option<checksum::Sanitize>,

    /// Write the compatibility report to FILE
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Include QC pass and fail rows
    #[arg(short = 'q', long = "show-qc")]
    show_qc: bool,

    /// Include rows whose count is zero
    #[arg(short = 'v', long = "verbose", action = ArgAction::Count)]
    verbosity: u8,

    /// Use tab-delimited output
    #[arg(short = 'T', long)]
    tabs: bool,

    /// Use biobambam2 bamseqchksum output
    #[arg(short = 'B', long)]
    bamseqchksum: bool,

    /// Merge checksum reports instead of reading sequence data
    #[arg(short = 'm', long)]
    merge: bool,

    /// Check every stable content field
    #[arg(
        short = 'a',
        long,
        conflicts_with_all = [
            "excluded_flags",
            "required_flags",
            "flag_mask",
            "no_reverse_complement",
            "tags",
            "order",
            "check_position",
            "check_cigar",
            "check_mate",
            "sanitize"
        ]
    )]
    all: bool,

    /// Additional BAM decompression workers
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,

    /// Input SAM, BAM, CRAM, FASTA, or FASTQ files; use - for standard input
    #[arg(value_name = "INPUT", num_args = 0..)]
    inputs: Vec<PathBuf>,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if arguments.threads > 256 {
        return Err(RsomicsError::ConfigError(
            "checksum additional thread count cannot exceed 256".to_owned(),
        ));
    }
    let inputs = if arguments.inputs.is_empty() {
        if arguments.merge {
            return Err(RsomicsError::ConfigError(
                "--merge requires at least one checksum report".to_owned(),
            ));
        }
        vec![PathBuf::from("-")]
    } else {
        arguments.inputs
    };
    let named_output = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"));
    if json && named_output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for checksum".to_owned(),
        ));
    }
    if let Some(output) = named_output {
        for input in &inputs {
            if input != Path::new("-") && same_target(input, output)? {
                return Err(RsomicsError::ConfigError(
                    "checksum output must differ from every input".to_owned(),
                ));
            }
        }
    }
    if arguments.bamseqchksum
        && (arguments.all
            || arguments.check_position
            || arguments.check_cigar
            || arguments.check_mate)
    {
        return Err(RsomicsError::ConfigError(
            "--bamseqchksum cannot include position, CIGAR, or mate columns".to_owned(),
        ));
    }

    let options = checksum::Options {
        required_flags: if arguments.all {
            0
        } else {
            arguments.required_flags
        },
        excluded_flags: if arguments.all {
            0
        } else {
            arguments.excluded_flags
        },
        flag_mask: if arguments.all {
            0xfff
        } else {
            arguments.flag_mask
        },
        reverse_complement: !arguments.all && !arguments.no_reverse_complement,
        tags: if arguments.all {
            TagSelection::AllExcept(vec![*b"cF", *b"MD", *b"NM"])
        } else {
            arguments.tags
        },
        order: if arguments.all { 1 } else { arguments.order },
        check_position: arguments.all || arguments.check_position,
        check_cigar: arguments.all || arguments.check_cigar,
        check_mate: arguments.all || arguments.check_mate,
        maximum_records: arguments.maximum_records.filter(|&limit| limit != 0),
        show_qc: arguments.show_qc || arguments.bamseqchksum,
        verbose: arguments.verbosity > 0 || arguments.bamseqchksum,
        tabs: arguments.tabs,
        bamseqchksum: arguments.bamseqchksum,
        sanitize: if arguments.all {
            checksum::Sanitize::all()
        } else {
            arguments.sanitize.unwrap_or_default()
        },
        additional_threads: arguments.threads,
    };
    if arguments.merge && arguments.order > 1 {
        return Err(RsomicsError::ConfigError(
            "absolute double-order checksums cannot be merged".to_owned(),
        ));
    }
    let reports = if arguments.merge {
        vec![checksum::merge(&inputs, &options)?]
    } else {
        inputs
            .iter()
            .map(|input| checksum::collect(input, &options))
            .collect::<Result<Vec<_>>>()?
    };

    if let Some(path) = named_output {
        let mut transaction = TransactionalFile::new(path)?;
        for report in &reports {
            report
                .write(transaction.file_mut())
                .map_err(RsomicsError::Io)?;
        }
        transaction.file_mut().flush().map_err(RsomicsError::Io)?;
        transaction.commit()?;
    } else {
        let mut output = io::stdout().lock();
        for report in &reports {
            report.write(&mut output).map_err(RsomicsError::Io)?;
        }
        output.flush().map_err(RsomicsError::Io)?;
    }

    Ok(CommandOutput::Checksum { reports })
}

fn parse_flags(value: &str) -> std::result::Result<u16, String> {
    crate::flags::parse(value).map_err(|error| error.to_string())
}

fn parse_tags(value: &str) -> std::result::Result<TagSelection, String> {
    let mut fields = value.split(',');
    let wildcard = fields.next() == Some("*");
    let values = if wildcard {
        fields.collect::<Vec<_>>()
    } else {
        value.split(',').collect()
    };
    let mut tags = Vec::with_capacity(values.len());
    for tag in values {
        let bytes: [u8; 2] = tag
            .as_bytes()
            .try_into()
            .map_err(|_| format!("tag must contain exactly two bytes: {tag:?}"))?;
        if !(b'0'..=b'z').contains(&bytes[0]) || !(b'0'..=b'z').contains(&bytes[1]) {
            return Err(format!("illegal tag ID {tag:?}"));
        }
        if tags.contains(&bytes) {
            return Err(format!("duplicate tag ID {tag:?}"));
        }
        tags.push(bytes);
    }
    if wildcard {
        Ok(TagSelection::AllExcept(tags))
    } else {
        Ok(TagSelection::Listed(tags))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_ordered_and_wildcard_tag_lists() {
        assert_eq!(
            parse_tags("BC,NM"),
            Ok(TagSelection::Listed(vec![*b"BC", *b"NM"]))
        );
        assert_eq!(
            parse_tags("*,cF,MD,NM"),
            Ok(TagSelection::AllExcept(vec![*b"cF", *b"MD", *b"NM"]))
        );
    }

    #[test]
    fn rejects_malformed_and_duplicate_tags() {
        for value in ["B", "ABC", "BC,,NM", "BC,BC"] {
            assert!(parse_tags(value).is_err(), "{value}");
        }
    }
}
