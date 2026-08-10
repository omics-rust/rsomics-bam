use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};

use clap::{Args, ValueEnum};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::output::{TransactionalFile, same_target};
use crate::{Program, addreplacerg};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
enum Mode {
    #[value(name = "overwrite_all")]
    #[default]
    OverwriteAll,
    #[value(name = "orphan_only")]
    OrphanOnly,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum Format {
    Sam,
    Bam,
}

#[derive(Debug, Args)]
#[command(after_help = "\
Examples:
  rsomics-bam addreplacerg -r ID:lane1 -r SM:sample reads.bam -o tagged.bam
  rsomics-bam addreplacerg -R lane1 reads.bam > tagged.sam
  rsomics-bam addreplacerg -r 'ID:new\\tSM:sample' -m orphan_only reads.bam")]
pub(crate) struct Arguments {
    /// Append a field to one new @RG record; repeat for multiple fields
    #[arg(
        short = 'r',
        long = "rg-line",
        value_name = "FIELD",
        conflicts_with = "existing_read_group"
    )]
    new_read_group: Vec<String>,

    /// Select an existing @RG ID from the input header
    #[arg(
        short = 'R',
        long = "rg-id",
        value_name = "ID",
        conflicts_with = "new_read_group"
    )]
    existing_read_group: Option<String>,

    /// Choose whether all records or only records without RG are changed
    #[arg(short = 'm', long, value_enum, default_value = "overwrite_all")]
    mode: Mode,

    /// Replace an existing header @RG with the same ID supplied by -r
    #[arg(short = 'w', long = "overwrite-header")]
    overwrite_header: bool,

    /// Write output to FILE instead of standard output
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Select SAM or BAM output
    #[arg(short = 'O', long = "output-fmt", value_enum)]
    format: Option<Format>,

    /// Write BAM without DEFLATE compression
    #[arg(short = 'u', long)]
    uncompressed: bool,

    /// Reference FASTA for CRAM decoding
    #[arg(long, value_name = "FASTA")]
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

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let input = arguments.input.as_deref().unwrap_or_else(|| Path::new("-"));
    let output = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"));
    if json && output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for addreplacerg".to_owned(),
        ));
    }
    if let Some(output) = output
        && input != Path::new("-")
        && same_target(input, output)?
    {
        return Err(RsomicsError::ConfigError(
            "addreplacerg input and output must be different files".to_owned(),
        ));
    }

    let new_read_group =
        (!arguments.new_read_group.is_empty()).then(|| arguments.new_read_group.join("\t"));
    let source = match (
        new_read_group.as_deref(),
        arguments.existing_read_group.as_deref(),
    ) {
        (Some(fields), None) => addreplacerg::Source::New(fields),
        (None, Some(id)) => addreplacerg::Source::Existing(id),
        (None, None) => addreplacerg::Source::First,
        (Some(_), Some(_)) => unreachable!("clap rejects conflicting read-group sources"),
    };
    let command_line = (!arguments.suppress_program_record).then(crate::program::command_line);
    let program = command_line
        .as_deref()
        .map(|line| Program::new("rsomics-bam", env!("CARGO_PKG_VERSION"), line))
        .transpose()?;
    let options = addreplacerg::Options {
        source,
        mode: match arguments.mode {
            Mode::OverwriteAll => addreplacerg::Mode::OverwriteAll,
            Mode::OrphanOnly => addreplacerg::Mode::OrphanOnly,
        },
        overwrite_header: arguments.overwrite_header,
        format: output_format(arguments.format, arguments.output.as_deref())?,
        compression: if arguments.uncompressed {
            addreplacerg::Compression::Uncompressed
        } else {
            addreplacerg::Compression::Default
        },
        reference: arguments.reference.as_deref(),
        additional_threads: arguments.threads,
        program,
    };

    let summary = match output {
        Some(path) => {
            let transaction = TransactionalFile::new(path)?;
            let summary =
                addreplacerg::write(input, options, BufWriter::new(transaction.reopen()?))?;
            transaction.commit()?;
            summary
        }
        None => addreplacerg::write(input, options, io::stdout())?,
    };
    Ok(CommandOutput::Addreplacerg { summary })
}

fn output_format(format: Option<Format>, output: Option<&Path>) -> Result<addreplacerg::Format> {
    match format {
        Some(Format::Sam) => Ok(addreplacerg::Format::Sam),
        Some(Format::Bam) => Ok(addreplacerg::Format::Bam),
        None if output.is_some_and(is_bam_path) => Ok(addreplacerg::Format::Bam),
        None if output.is_some_and(is_cram_path) => Err(RsomicsError::ConfigError(
            "CRAM output is not available; select -O sam or -O bam".to_owned(),
        )),
        None => Ok(addreplacerg::Format::Sam),
    }
}

fn is_bam_path(path: &Path) -> bool {
    path.extension()
        .and_then(std::ffi::OsStr::to_str)
        .is_some_and(|extension| extension.eq_ignore_ascii_case("bam"))
}

fn is_cram_path(path: &Path) -> bool {
    path.extension()
        .and_then(std::ffi::OsStr::to_str)
        .is_some_and(|extension| extension.eq_ignore_ascii_case("cram"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bam_extension_is_case_insensitive() {
        assert_eq!(
            output_format(None, Some(Path::new("reads.BAM"))).unwrap(),
            addreplacerg::Format::Bam,
        );
        assert!(output_format(None, Some(Path::new("reads.cram"))).is_err());
    }
}
