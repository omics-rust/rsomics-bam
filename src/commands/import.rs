use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};

use clap::{Args, ValueEnum};
use rsomics_common::{Result, RsomicsError, reject_output_alias};

use crate::cli::CommandOutput;
use crate::output::TransactionalFile;
use crate::{Program, import};

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum Format {
    Sam,
    Bam,
}

#[derive(Debug, Args)]
#[command(after_help = "\
Examples:
  rsomics-bam import reads.fastq > reads.sam
  rsomics-bam import read1.fastq read2.fastq -o reads.bam
  rsomics-bam import -0 reads.fastq -R sample -O bam -o reads.bam")]
pub(crate) struct Arguments {
    /// Single FASTQ whose /1 and /2 suffixes determine pairing
    #[arg(short = '0', long = "single", value_name = "FILE")]
    single: Option<PathBuf>,

    /// Interleaved FASTQ whose /1 and /2 suffixes determine pairing
    #[arg(short = 's', long, value_name = "FILE")]
    interleaved: Option<PathBuf>,

    /// Read-1 FASTQ of a paired run
    #[arg(short = '1', long = "read1", value_name = "FILE")]
    read1: Option<PathBuf>,

    /// Read-2 FASTQ of a paired run
    #[arg(short = '2', long = "read2", value_name = "FILE")]
    read2: Option<PathBuf>,

    /// One single or two paired FASTQ inputs
    #[arg(value_name = "FASTQ", num_args = 0..=2)]
    inputs: Vec<PathBuf>,

    /// Write output to FILE instead of standard output
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Select SAM or BAM output
    #[arg(short = 'O', long = "output-fmt", value_enum)]
    format: Option<Format>,

    /// Write BAM without DEFLATE compression
    #[arg(short = 'u', long)]
    uncompressed: bool,

    /// Append a field to one @RG line; repeat for multiple fields
    #[arg(short = 'r', long = "rg-line", value_name = "FIELD")]
    read_group_fields: Vec<String>,

    /// Add a read group with ID
    #[arg(short = 'R', long = "rg", value_name = "ID")]
    read_group_id: Option<String>,

    /// Add the zero-based input record number in TAG or fixed-width TAG:WIDTH
    #[arg(long, value_name = "TAG[:WIDTH]", value_parser = parse_order)]
    order: Option<import::OrderTag>,

    /// Additional BAM compression workers
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line
    #[arg(long = "no-pg", visible_alias = "no-PG")]
    suppress_program_record: bool,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let inputs = resolve_inputs(&arguments)?;
    let named_output = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"));
    if json && named_output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for import".to_owned(),
        ));
    }
    if let Some(output) = named_output {
        match inputs {
            import::Inputs::Single(input) if input != Path::new("-") => {
                reject_output_alias(output, [input])?;
            }
            import::Inputs::Paired(read1, read2) => {
                reject_output_alias(output, [read1, read2])?;
            }
            import::Inputs::Single(_) => {}
        }
    }

    let format = output_format(arguments.format, arguments.output.as_deref());
    let read_group = read_group(
        &arguments.read_group_fields,
        arguments.read_group_id.as_deref(),
    )?;
    let command_line = (!arguments.suppress_program_record).then(crate::program::command_line);
    let program = command_line
        .as_deref()
        .map(|line| Program::new("rsomics-bam", env!("CARGO_PKG_VERSION"), line))
        .transpose()?;
    let options = import::Options {
        format,
        compression: if arguments.uncompressed {
            import::Compression::Uncompressed
        } else {
            import::Compression::Default
        },
        additional_threads: arguments.threads,
        read_group: read_group.as_deref(),
        order: arguments.order,
        program,
    };

    let summary = match named_output {
        Some(path) => {
            let transaction = TransactionalFile::new(path)?;
            let file = transaction.reopen()?;
            let summary = import::write(inputs, options, BufWriter::new(file))?;
            transaction.commit()?;
            summary
        }
        None => import::write(inputs, options, io::stdout())?,
    };
    Ok(CommandOutput::Import { summary })
}

fn resolve_inputs(arguments: &Arguments) -> Result<import::Inputs<'_>> {
    let singles =
        usize::from(arguments.single.is_some()) + usize::from(arguments.interleaved.is_some());
    let paired = arguments.read1.is_some() || arguments.read2.is_some();
    let explicit = singles != 0 || paired;

    if explicit && !arguments.inputs.is_empty() {
        return Err(RsomicsError::ConfigError(
            "positional FASTQ inputs cannot be combined with -0, -s, -1, or -2".to_owned(),
        ));
    }
    if singles > 1 || (singles != 0 && paired) {
        return Err(RsomicsError::ConfigError(
            "select exactly one of -0, -s, or the -1/-2 pair".to_owned(),
        ));
    }
    if singles == 1 {
        return Ok(import::Inputs::Single(
            arguments
                .single
                .as_deref()
                .or(arguments.interleaved.as_deref())
                .expect("one single input is present"),
        ));
    }
    if paired {
        return match (arguments.read1.as_deref(), arguments.read2.as_deref()) {
            (Some(read1), Some(read2)) => Ok(import::Inputs::Paired(read1, read2)),
            _ => Err(RsomicsError::ConfigError(
                "paired import requires both -1 and -2".to_owned(),
            )),
        };
    }
    match arguments.inputs.as_slice() {
        [single] => Ok(import::Inputs::Single(single)),
        [read1, read2] => Ok(import::Inputs::Paired(read1, read2)),
        [] => Err(RsomicsError::ConfigError(
            "import requires one or two FASTQ inputs".to_owned(),
        )),
        _ => unreachable!("clap limits positional inputs to two"),
    }
}

fn output_format(format: Option<Format>, output: Option<&Path>) -> import::Format {
    match format {
        Some(Format::Sam) => import::Format::Sam,
        Some(Format::Bam) => import::Format::Bam,
        None if output.is_some_and(is_bam_path) => import::Format::Bam,
        None => import::Format::Sam,
    }
}

fn is_bam_path(path: &Path) -> bool {
    path.extension()
        .and_then(std::ffi::OsStr::to_str)
        .is_some_and(|extension| extension.eq_ignore_ascii_case("bam"))
}

fn read_group(fields: &[String], id: Option<&str>) -> Result<Option<String>> {
    if !fields.is_empty() {
        return Ok(Some(fields.join("\t")));
    }
    match id {
        Some(id) if id.is_empty() || id.contains(['\t', '\r', '\n']) => Err(
            RsomicsError::InvalidInput("read-group ID must be one nonempty field".to_owned()),
        ),
        Some(id) => Ok(Some(format!("ID:{id}"))),
        None => Ok(None),
    }
}

fn parse_order(value: &str) -> std::result::Result<import::OrderTag, String> {
    import::OrderTag::parse(value).map_err(|error| error.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bam_extension_is_case_insensitive() {
        assert_eq!(
            output_format(None, Some(Path::new("reads.BAM"))),
            import::Format::Bam
        );
        assert_eq!(
            output_format(None, Some(Path::new("reads.data"))),
            import::Format::Sam
        );
    }

    #[test]
    fn repeated_read_group_fields_form_one_line() {
        assert_eq!(
            read_group(&["ID:lib1".to_owned(), "SM:sample".to_owned()], None).unwrap(),
            Some("ID:lib1\tSM:sample".to_owned())
        );
    }
}
