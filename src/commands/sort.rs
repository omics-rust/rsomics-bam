use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::output::{TransactionalFile, same_target};
use crate::{Program, hts_quickcheck, sort};

#[derive(Debug, Args)]
#[command(after_help = "\
Examples:
  rsomics-bam sort input.bam -o sorted.bam
  rsomics-bam sort -n input.bam -o names.bam
  rsomics-bam sort --template-coordinate input.bam -o templates.bam
  rsomics-bam sort -m 256M -T /scratch/sample input.bam -o sorted.bam")]
pub(crate) struct Arguments {
    /// Input SAM, BAM, or CRAM file; omit or use - for standard input
    #[arg(value_name = "ALIGNMENT")]
    input: Option<PathBuf>,

    /// Write BAM output to FILE; omit or use - for standard output
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Sort read names in natural numeric order
    #[arg(short = 'n', long, conflicts_with_all = ["ascii_name", "template_coordinate"])]
    natural_name: bool,

    /// Sort read names in bytewise lexicographical order
    #[arg(short = 'N', long, conflicts_with_all = ["natural_name", "template_coordinate"])]
    ascii_name: bool,

    /// Sort by unclipped template coordinates
    #[arg(long, conflicts_with_all = ["natural_name", "ascii_name"])]
    template_coordinate: bool,

    /// Total in-memory record budget; suffix K, M, or G is accepted
    #[arg(short = 'm', long, value_name = "SIZE", default_value = "768M", value_parser = parse_memory)]
    memory: u64,

    /// Prefix for temporary BAM runs
    #[arg(short = 'T', long, value_name = "PREFIX")]
    temporary_prefix: Option<PathBuf>,

    /// Reference FASTA for CRAM decoding
    #[arg(long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Additional sort and BAM I/O workers; omit for automatic parallelism
    #[arg(short = '@', long, value_name = "INT")]
    threads: Option<usize>,

    /// Do not add an @PG line
    #[arg(long = "no-pg", visible_alias = "no-PG")]
    suppress_program_record: bool,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let input = arguments.input.as_deref().unwrap_or_else(|| Path::new("-"));
    let output = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"));
    if json && output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for sort".to_owned(),
        ));
    }
    if let Some(output) = output
        && input != Path::new("-")
        && same_target(input, output)?
    {
        return Err(RsomicsError::ConfigError(
            "sort input and output must be different files".to_owned(),
        ));
    }

    let order = if arguments.natural_name {
        sort::Order::QueryNameNatural
    } else if arguments.ascii_name {
        sort::Order::QueryNameLexicographical
    } else if arguments.template_coordinate {
        sort::Order::TemplateCoordinate
    } else {
        sort::Order::Coordinate
    };
    let command_line = (!arguments.suppress_program_record).then(crate::program::command_line);
    let program = command_line
        .as_deref()
        .map(|line| Program::new("rsomics-bam", env!("CARGO_PKG_VERSION"), line))
        .transpose()?;
    let derived_prefix = output.map(|path| {
        let mut value = path.as_os_str().to_os_string();
        value.push(".tmp");
        PathBuf::from(value)
    });
    let options = sort::Options {
        order,
        memory_limit: arguments.memory,
        additional_threads: arguments.threads,
        temporary_prefix: arguments
            .temporary_prefix
            .as_deref()
            .or(derived_prefix.as_deref()),
        reference: arguments.reference.as_deref(),
        destination: output,
        program,
    };

    let summary = match output {
        Some(path) => {
            let transaction = TransactionalFile::new(path)?;
            let file = transaction.reopen()?;
            let summary = sort::write(input, options, BufWriter::new(file))?;
            hts_quickcheck::require_bgzf_eof(transaction.temporary_path())?;
            transaction.commit()?;
            summary
        }
        None => sort::write(input, options, io::stdout())?,
    };
    Ok(CommandOutput::Sort { summary })
}

fn parse_memory(value: &str) -> std::result::Result<u64, String> {
    let (number, multiplier) = match value.as_bytes().last().copied() {
        Some(b'k' | b'K') => (&value[..value.len() - 1], 1u64 << 10),
        Some(b'm' | b'M') => (&value[..value.len() - 1], 1u64 << 20),
        Some(b'g' | b'G') => (&value[..value.len() - 1], 1u64 << 30),
        _ => (value, 1),
    };
    let bytes = number
        .parse::<u64>()
        .map_err(|_| format!("invalid memory size: {value}"))?
        .checked_mul(multiplier)
        .ok_or_else(|| format!("memory size overflows: {value}"))?;
    if bytes < 1 << 20 {
        return Err("sort memory must be at least 1M".to_owned());
    }
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_memory_units() {
        assert_eq!(parse_memory("1M"), Ok(1 << 20));
        assert_eq!(parse_memory("2G"), Ok(2 << 30));
        assert!(parse_memory("1023K").is_err());
        assert!(parse_memory("1.5G").is_err());
    }
}
