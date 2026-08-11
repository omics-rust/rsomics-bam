use std::io::{self, Write as _};
use std::path::{Path, PathBuf};

use clap::{Args, ValueEnum};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::output::same_target;
use crate::{Program, reset};

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum Format {
    Sam,
    Bam,
    Cram,
}

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Write output to FILE instead of standard output
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Select SAM, BAM, or CRAM output
    #[arg(short = 'O', long = "output-fmt", value_enum, ignore_case = true)]
    format: Option<Format>,

    /// Remove comma-separated auxiliary tags
    #[arg(short = 'x', long = "remove-tag", value_name = "TAG[,TAG...]")]
    remove_tags: Vec<String>,

    /// Keep only comma-separated auxiliary tags
    #[arg(long = "keep-tag", value_name = "TAG[,TAG...]")]
    keep_tags: Vec<String>,

    /// Drop @RG lines and RG auxiliary tags
    #[arg(long = "no-RG")]
    no_rg: bool,

    /// Drop the matching @PG line and all following @PG lines
    #[arg(long = "reject-PG", value_name = "ID")]
    reject_pg: Option<String>,

    /// Preserve the duplicate FLAG bit
    #[arg(long)]
    dupflag: bool,

    /// Reference FASTA for CRAM input or output
    #[arg(short = 'T', long = "reference", value_name = "FILE")]
    reference: Option<PathBuf>,

    /// Additional alignment I/O workers
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line
    #[arg(long = "no-pg", visible_alias = "no-PG")]
    suppress_program_record: bool,

    /// Input SAM, BAM, or CRAM file; use - for standard input
    #[arg(value_name = "ALIGNMENT")]
    input: PathBuf,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let named_output = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"));
    if arguments.input != Path::new("-")
        && let Some(output) = named_output
        && same_target(&arguments.input, output)?
    {
        return Err(RsomicsError::ConfigError(
            "reset input and output must be different files".to_owned(),
        ));
    }
    if json
        && arguments
            .output
            .as_deref()
            .is_none_or(|path| path == Path::new("-"))
    {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for reset".to_owned(),
        ));
    }

    let command_line = (!arguments.suppress_program_record).then(crate::program::command_line);
    let program = command_line
        .as_deref()
        .map(|line| Program::new("rsomics-bam", env!("CARGO_PKG_VERSION"), line))
        .transpose()?;
    let format = output_format(arguments.format, arguments.output.as_deref());
    let mut remove_values = Vec::new();
    let mut keep_values = arguments.keep_tags.clone();
    for value in &arguments.remove_tags {
        if let Some(value) = value.strip_prefix('^') {
            keep_values.push(value.to_owned());
        } else {
            remove_values.push(value.clone());
        }
    }
    let remove_tags = parse_tag_lists(&remove_values)?;
    let keep_tags = (!keep_values.is_empty())
        .then(|| parse_tag_lists(&keep_values))
        .transpose()?;
    let options = reset::Options {
        format,
        remove_tags: &remove_tags,
        keep_tags: keep_tags.as_deref(),
        no_rg: arguments.no_rg,
        reject_pg: arguments.reject_pg.as_deref(),
        keep_duplicate: arguments.dupflag,
        reference: arguments.reference.as_deref(),
        additional_threads: arguments.threads,
        program,
    };

    let summary = if format == reset::Format::Cram {
        if let Some(path) = named_output {
            let transaction = crate::output::TransactionalFile::new(path)?;
            let summary =
                reset::write_cram_path(&arguments.input, options, transaction.temporary_path())?;
            transaction.commit()?;
            summary
        } else {
            let temporary = tempfile::NamedTempFile::new().map_err(RsomicsError::Io)?;
            let summary = reset::write_cram_path(&arguments.input, options, temporary.path())?;
            let mut input = temporary.reopen().map_err(RsomicsError::Io)?;
            let mut output = io::stdout().lock();
            io::copy(&mut input, &mut output).map_err(RsomicsError::Io)?;
            output.flush().map_err(RsomicsError::Io)?;
            summary
        }
    } else if let Some(path) = named_output {
        let transaction = crate::output::TransactionalFile::new(path)?;
        let summary = reset::write(&arguments.input, options, transaction.reopen()?)?;
        transaction.commit()?;
        summary
    } else {
        reset::write(&arguments.input, options, io::stdout())?
    };

    Ok(CommandOutput::Reset { summary })
}

fn parse_tag_lists(values: &[String]) -> Result<Vec<[u8; 2]>> {
    let mut tags = Vec::new();
    for value in values {
        let bytes = value.as_bytes();
        let mut position = 0;
        while bytes.len().saturating_sub(position) >= 2 {
            let tag = [bytes[position], bytes[position + 1]];
            position += 2;
            if !tags.contains(&tag) {
                tags.push(tag);
            }
            if position == bytes.len() {
                break;
            }
            if bytes[position] != b',' {
                return Err(invalid_tag_list(value));
            }
            position += 1;
        }
        if position != bytes.len() {
            return Err(invalid_tag_list(value));
        }
    }
    Ok(tags)
}

fn invalid_tag_list(value: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "auxiliary tags must contain exactly two bytes each: {value}"
    ))
}

fn output_format(format: Option<Format>, output: Option<&Path>) -> reset::Format {
    match format {
        Some(Format::Sam) => reset::Format::Sam,
        Some(Format::Bam) => reset::Format::Bam,
        Some(Format::Cram) => reset::Format::Cram,
        None if output.is_some_and(|path| has_extension(path, "bam")) => reset::Format::Bam,
        None if output.is_some_and(|path| has_extension(path, "cram")) => reset::Format::Cram,
        None => reset::Format::Sam,
    }
}

fn has_extension(path: &Path, expected: &str) -> bool {
    path.extension()
        .and_then(std::ffi::OsStr::to_str)
        .is_some_and(|extension| extension.eq_ignore_ascii_case(expected))
}
