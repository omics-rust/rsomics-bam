use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use clap::{ArgAction, Args, ValueEnum};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::{flags, view};

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
enum Format {
    Sam,
    Bam,
}

#[derive(Debug, Args)]
#[command(disable_help_flag = true)]
pub(crate) struct Arguments {
    /// Input SAM, BAM, or CRAM file; omit or use - for standard input
    #[arg(value_name = "ALIGNMENT")]
    input: Option<PathBuf>,

    /// Indexed region such as chr1:1-100 or * for unmapped records
    #[arg(index = 2, value_name = "REGION", value_parser = parse_region)]
    regions: Vec<view::Region>,

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

    /// Write BAM output
    #[arg(short = 'b', long, conflicts_with = "format")]
    bam: bool,

    /// Select the output alignment format
    #[arg(short = 'O', long = "output-fmt", value_enum, conflicts_with = "bam")]
    format: Option<Format>,

    /// Use fast BAM compression
    #[arg(short = '1', long, conflicts_with = "uncompressed")]
    fast: bool,

    /// Write uncompressed BAM
    #[arg(short = 'u', long, conflicts_with = "fast")]
    uncompressed: bool,

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

    /// Additional alignment I/O thread budget
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
        output_format: output_format(&arguments)?,
        compression: compression(&arguments),
        reference: arguments.reference.as_deref(),
        additional_threads: arguments.threads,
        regions: &arguments.regions,
    };
    let input = arguments.input.as_deref().unwrap_or_else(|| Path::new("-"));

    let summary = if json {
        view::write(input, options, io::sink())?
    } else if let Some(path) = arguments.output.as_deref() {
        run_to_path(input, options, path)?
    } else {
        run_to(input, options, io::stdout())?
    };

    Ok(CommandOutput::View { summary })
}

fn run_to(
    input: &std::path::Path,
    options: view::Options<'_>,
    mut output: impl Write + Send + 'static,
) -> Result<view::Summary> {
    if options.count_only {
        let summary = view::write(input, options, io::sink())?;
        writeln!(output, "{}", summary.selected).map_err(RsomicsError::Io)?;
        output.flush().map_err(RsomicsError::Io)?;
        Ok(summary)
    } else {
        view::write(input, options, output)
    }
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
    let file = temporary.reopen().map_err(|error| {
        RsomicsError::Io(io::Error::new(
            error.kind(),
            format!(
                "opening temporary output beside {}: {error}",
                output.display()
            ),
        ))
    })?;
    let summary = run_to(input, options, BufWriter::new(file))?;
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

fn parse_region(value: &str) -> std::result::Result<view::Region, String> {
    value.parse()
}

fn output_format(arguments: &Arguments) -> Result<view::Format> {
    if (arguments.fast || arguments.uncompressed) && arguments.format == Some(Format::Sam) {
        return Err(RsomicsError::ConfigError(
            "BAM compression options cannot be combined with SAM output".to_owned(),
        ));
    }
    if arguments.bam {
        return Ok(view::Format::Bam);
    }
    if let Some(format) = arguments.format {
        return Ok(match format {
            Format::Sam => view::Format::Sam,
            Format::Bam => view::Format::Bam,
        });
    }
    if arguments.fast || arguments.uncompressed {
        return Ok(view::Format::Bam);
    }

    Ok(
        match arguments
            .output
            .as_deref()
            .and_then(Path::extension)
            .and_then(|extension| extension.to_str())
        {
            Some("bam") => view::Format::Bam,
            _ => view::Format::Sam,
        },
    )
}

fn compression(arguments: &Arguments) -> view::Compression {
    if arguments.fast {
        view::Compression::Fast
    } else if arguments.uncompressed {
        view::Compression::Uncompressed
    } else {
        view::Compression::Default
    }
}
