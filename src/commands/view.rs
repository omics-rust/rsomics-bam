use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use clap::{ArgAction, Args};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::{flags, view};

#[derive(Debug, Args)]
#[command(disable_help_flag = true)]
pub(crate) struct Arguments {
    /// Input SAM, BAM, or CRAM file; omit or use - for standard input
    #[arg(value_name = "ALIGNMENT", default_value = "-")]
    input: PathBuf,

    /// Include the SAM header
    #[arg(short = 'h', long, conflicts_with = "header_only")]
    with_header: bool,

    /// Print only the SAM header
    #[arg(short = 'H', long, conflicts_with = "with_header")]
    header_only: bool,

    /// Explicitly omit the SAM header
    #[arg(long, conflicts_with_all = ["with_header", "header_only"])]
    no_header: bool,

    /// Print only the count of selected records
    #[arg(short = 'c', long, conflicts_with = "header_only")]
    count: bool,

    /// Write output to a file
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Require all FLAG bits
    #[arg(short = 'f', long, value_name = "FLAG", value_parser = parse_flags)]
    require_flags: Option<u16>,

    /// Exclude records with any FLAG bits
    #[arg(short = 'F', long, value_name = "FLAG", value_parser = parse_flags)]
    exclude_flags: Option<u16>,

    /// Require at least one FLAG bit
    #[arg(long = "rf", aliases = ["incl-flags", "include-flags"], value_name = "FLAG", value_parser = parse_flags)]
    include_flags: Option<u16>,

    /// Exclude records with all FLAG bits
    #[arg(short = 'G', value_name = "FLAG", value_parser = parse_flags)]
    exclude_all_flags: Option<u16>,

    /// Minimum mapping quality
    #[arg(short = 'q', long = "min-MQ", value_name = "INT", default_value_t = 0)]
    minimum_mapping_quality: u8,

    /// Reference FASTA for CRAM decoding
    #[arg(short = 'T', long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Additional alignment decompression threads
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,

    /// Print command help
    #[arg(short = '?', long, action = ArgAction::Help)]
    help: Option<bool>,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if json && !arguments.count {
        return Err(RsomicsError::ConfigError(
            "--json requires --count for view".to_owned(),
        ));
    }
    if json && arguments.output.is_some() {
        return Err(RsomicsError::ConfigError(
            "--json cannot be combined with --output".to_owned(),
        ));
    }

    let options = view::Options {
        with_header: arguments.with_header,
        header_only: arguments.header_only,
        count_only: arguments.count,
        require_flags: arguments.require_flags.unwrap_or_default(),
        exclude_flags: arguments.exclude_flags.unwrap_or_default(),
        include_flags: arguments.include_flags.unwrap_or_default(),
        exclude_all_flags: arguments.exclude_all_flags.unwrap_or_default(),
        minimum_mapping_quality: arguments.minimum_mapping_quality,
        reference: arguments.reference.as_deref(),
        additional_threads: arguments.threads,
    };

    let summary = if json {
        view::write(&arguments.input, options, io::sink())?
    } else if let Some(path) = arguments.output {
        run_to_path(&arguments.input, options, &path)?
    } else {
        run_to(&arguments.input, options, io::stdout().lock())?
    };

    Ok(CommandOutput::View { summary })
}

fn run_to(
    input: &std::path::Path,
    options: view::Options<'_>,
    mut output: impl Write,
) -> Result<view::Summary> {
    let summary = view::write(input, options, &mut output)?;
    if options.count_only {
        writeln!(output, "{}", summary.selected).map_err(RsomicsError::Io)?;
        output.flush().map_err(RsomicsError::Io)?;
    }
    Ok(summary)
}

fn run_to_path(input: &Path, options: view::Options<'_>, output: &Path) -> Result<view::Summary> {
    let parent = output
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let permissions = std::fs::metadata(output)
        .ok()
        .map(|metadata| metadata.permissions());
    let mut temporary = tempfile::NamedTempFile::new_in(parent).map_err(|error| {
        RsomicsError::Io(io::Error::new(
            error.kind(),
            format!(
                "creating temporary output beside {}: {error}",
                output.display()
            ),
        ))
    })?;
    let summary = run_to(input, options, BufWriter::new(temporary.as_file_mut()))?;
    temporary
        .as_file_mut()
        .sync_all()
        .map_err(RsomicsError::Io)?;
    if let Some(permissions) = permissions {
        temporary
            .as_file_mut()
            .set_permissions(permissions)
            .map_err(RsomicsError::Io)?;
    }
    temporary.persist(output).map_err(|error| {
        RsomicsError::Io(io::Error::new(
            error.error.kind(),
            format!("committing output {}: {}", output.display(), error.error),
        ))
    })?;
    Ok(summary)
}

fn parse_flags(value: &str) -> std::result::Result<u16, String> {
    flags::parse(value).map_err(|error| error.to_string())
}
