use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{Result, RsomicsError, reject_output_alias, write_output};
use rsomics_seqio::Compression;

use crate::cli::CommandOutput;
use crate::{fastx, flags};

#[derive(Debug, Args)]
pub(crate) struct FastaArguments {
    #[command(flatten)]
    common: CommonArguments,
}

#[derive(Debug, Args)]
pub(crate) struct FastqArguments {
    #[command(flatten)]
    common: CommonArguments,

    /// Use original qualities from the OQ tag when present
    #[arg(short = 'O', long)]
    original_quality: bool,

    /// Phred quality used when a record has no qualities
    #[arg(short = 'v', long, default_value_t = 1, value_parser = clap::value_parser!(u8).range(0..=93))]
    default_quality: u8,
}

#[derive(Debug, Args)]
struct CommonArguments {
    /// Input SAM, BAM, or CRAM file; omit or use - for standard input
    #[arg(value_name = "ALIGNMENT")]
    input: Option<PathBuf>,

    /// Write the unified FASTA or FASTQ stream to a file
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Do not append /1 and /2 to paired read names
    #[arg(short = 'n', long)]
    no_mate_suffix: bool,

    /// Require all FLAG bits
    #[arg(short = 'f', long, value_name = "FLAG", value_parser = parse_flags)]
    require_flags: Option<u16>,

    /// Exclude records with any FLAG bits instead of the default 0x900
    #[arg(short = 'F', long, value_name = "FLAG", value_parser = parse_flags)]
    exclude_flags: Option<u16>,

    /// Require at least one FLAG bit
    #[arg(long = "rf", aliases = ["incl-flags", "include-flags"], value_name = "FLAG", value_parser = parse_flags)]
    include_flags: Option<u16>,

    /// Exclude records with all FLAG bits
    #[arg(short = 'G', value_name = "FLAG", value_parser = parse_flags)]
    exclude_all_flags: Option<u16>,

    /// BGZF compression level for .gz, .bgz, and .bgzf output
    #[arg(short = 'c', long, default_value_t = 1, value_parser = clap::value_parser!(u8).range(0..=9))]
    compression_level: u8,

    /// Reference FASTA for CRAM decoding
    #[arg(long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Additional alignment decompression threads
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,
}

pub(crate) fn execute_fasta(arguments: FastaArguments, json: bool) -> Result<CommandOutput> {
    execute(arguments.common, json, fastx::Format::Fasta, false, 1)
}

pub(crate) fn execute_fastq(arguments: FastqArguments, json: bool) -> Result<CommandOutput> {
    execute(
        arguments.common,
        json,
        fastx::Format::Fastq,
        arguments.original_quality,
        arguments.default_quality,
    )
}

fn execute(
    arguments: CommonArguments,
    json: bool,
    format: fastx::Format,
    use_original_quality: bool,
    default_quality: u8,
) -> Result<CommandOutput> {
    if json
        && arguments
            .output
            .as_deref()
            .is_none_or(|path| path == Path::new("-"))
    {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for fasta and fastq".to_owned(),
        ));
    }

    let input = arguments.input.as_deref().unwrap_or_else(|| Path::new("-"));
    if let Some(output) = arguments.output.as_deref() {
        reject_output_alias(output, [input])?;
    }
    let options = fastx::Options {
        format,
        compression: compression(arguments.output.as_deref(), arguments.compression_level),
        append_mate_suffix: !arguments.no_mate_suffix,
        use_original_quality,
        default_quality,
        require_flags: arguments.require_flags.unwrap_or_default(),
        exclude_flags: arguments.exclude_flags.unwrap_or(0x900),
        include_flags: arguments.include_flags.unwrap_or_default(),
        exclude_all_flags: arguments.exclude_all_flags.unwrap_or_default(),
        reference: arguments.reference.as_deref(),
        additional_threads: arguments.threads,
    };
    let summary = write_output(arguments.output.as_deref(), |output| {
        fastx::write(input, options, output)
    })?;
    Ok(match format {
        fastx::Format::Fasta => CommandOutput::Fasta { summary },
        fastx::Format::Fastq => CommandOutput::Fastq { summary },
    })
}

fn compression(path: Option<&Path>, level: u8) -> Compression {
    let extension = path
        .and_then(Path::extension)
        .and_then(std::ffi::OsStr::to_str)
        .map(str::to_ascii_lowercase);
    if matches!(extension.as_deref(), Some("gz" | "bgz" | "bgzf")) {
        Compression::Bgzf { level }
    } else {
        Compression::Plain
    }
}

fn parse_flags(value: &str) -> std::result::Result<u16, String> {
    flags::parse(value).map_err(|error| error.to_string())
}
