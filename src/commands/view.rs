use std::ffi::OsStr;
use std::fs;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use clap::{ArgAction, Args, ValueEnum};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

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

    /// Do not add an @PG line to the output header
    #[arg(long = "no-pg", visible_alias = "no-PG")]
    suppress_program_record: bool,

    /// Print only the count of selected records
    #[arg(short = 'c', long, conflicts_with = "header_only")]
    count: bool,

    /// Write output to a file
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Write processed, accepted, and rejected record counts as JSON
    #[arg(long, value_name = "FILE")]
    save_counts: Option<PathBuf>,

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

    /// Select records in this read group or with no read group
    #[arg(short = 'r', long = "read-group", value_name = "STR")]
    read_groups: Vec<String>,

    /// Select read names listed in [^]FILE; ^ negates the selection
    #[arg(short = 'N', long = "qname-file", value_name = "[^]FILE")]
    qname_files: Vec<PathBuf>,

    /// Minimum mapping quality
    #[arg(short = 'q', long = "min-MQ", value_name = "INT", default_value_t = 0)]
    minimum_mapping_quality: u8,

    /// Minimum query length measured from the CIGAR
    #[arg(
        short = 'm',
        long = "min-qlen",
        value_name = "INT",
        default_value_t = 0
    )]
    minimum_query_length: u64,

    /// Add FLAG bits to selected output records
    #[arg(long, value_name = "FLAG", value_parser = parse_flags)]
    add_flags: Vec<u16>,

    /// Remove FLAG bits from selected output records
    #[arg(long, value_name = "FLAG", value_parser = parse_flags)]
    remove_flags: Vec<u16>,

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

    let command_line = (!arguments.suppress_program_record).then(command_line);
    let program = command_line
        .as_deref()
        .map(|command_line| {
            view::Program::new("rsomics-bam", env!("CARGO_PKG_VERSION"), command_line)
        })
        .transpose()?;
    let options = view::Options {
        with_header: arguments.with_header,
        header_only: arguments.header_only,
        count_only: arguments.count,
        require_flags: arguments.require_flags.unwrap_or_default(),
        exclude_flags: arguments.exclude_flags.unwrap_or_default(),
        include_flags: arguments.include_flags.unwrap_or_default(),
        exclude_all_flags: arguments.exclude_all_flags.unwrap_or_default(),
        read_groups: &arguments.read_groups,
        qname_files: &arguments.qname_files,
        minimum_mapping_quality: arguments.minimum_mapping_quality,
        minimum_query_length: arguments.minimum_query_length,
        add_flags: combine_flags(&arguments.add_flags),
        remove_flags: combine_flags(&arguments.remove_flags),
        output_format: output_format(&arguments)?,
        compression: compression(&arguments),
        reference: arguments.reference.as_deref(),
        additional_threads: arguments.threads,
        regions: &arguments.regions,
        program,
    };
    let input = arguments.input.as_deref().unwrap_or_else(|| Path::new("-"));
    validate_count_output(
        input,
        arguments.output.as_deref(),
        arguments.save_counts.as_deref(),
    )?;

    let summary = if json {
        view::write(input, options, io::sink())?
    } else if let Some(path) = arguments.output.as_deref() {
        run_to_path(input, options, path)?
    } else {
        run_to(input, options, io::stdout())?
    };
    if let Some(path) = arguments.save_counts.as_deref() {
        write_counts(path, summary)?;
    }

    Ok(CommandOutput::View { summary })
}

fn command_line() -> String {
    std::env::args_os()
        .map(|argument| sanitize_argument(&argument))
        .collect::<Vec<_>>()
        .join(" ")
}

fn sanitize_argument(argument: &OsStr) -> String {
    argument.to_string_lossy().replace(['\t', '\r', '\n'], " ")
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
    let transaction = TransactionalFile::new(output)?;
    let file = transaction.reopen()?;
    let summary = run_to(input, options, BufWriter::new(file))?;
    transaction.commit()?;
    Ok(summary)
}

#[derive(Serialize)]
struct SavedCounts {
    records_processed: u64,
    records_filter_accepted: u64,
    records_filter_rejected: u64,
}

fn write_counts(path: &Path, summary: view::Summary) -> Result<()> {
    let counts = SavedCounts {
        records_processed: summary
            .selected
            .checked_add(summary.rejected)
            .ok_or_else(|| {
                RsomicsError::InvalidInput("processed alignment count exceeds u64".to_owned())
            })?,
        records_filter_accepted: summary.selected,
        records_filter_rejected: summary.rejected,
    };
    let mut data = serde_json::to_vec_pretty(&counts)
        .map_err(|error| RsomicsError::InvalidInput(format!("serializing counts: {error}")))?;
    data.push(b'\n');

    let mut transaction = TransactionalFile::new(path)?;
    transaction
        .temporary
        .as_file_mut()
        .write_all(&data)
        .map_err(RsomicsError::Io)?;
    transaction.commit()
}

fn validate_count_output(input: &Path, output: Option<&Path>, counts: Option<&Path>) -> Result<()> {
    let Some(counts) = counts else {
        return Ok(());
    };
    if counts == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "--save-counts requires a named file".to_owned(),
        ));
    }
    if input != Path::new("-") && same_target(input, counts)? {
        return Err(RsomicsError::ConfigError(
            "--save-counts cannot overwrite the alignment input".to_owned(),
        ));
    }
    if let Some(output) = output
        && same_target(output, counts)?
    {
        return Err(RsomicsError::ConfigError(
            "--save-counts and --output require different files".to_owned(),
        ));
    }
    Ok(())
}

fn same_target(left: &Path, right: &Path) -> Result<bool> {
    Ok(target_identity(left)? == target_identity(right)?)
}

fn target_identity(path: &Path) -> Result<PathBuf> {
    match fs::canonicalize(path) {
        Ok(path) => return Ok(path),
        Err(error) if error.kind() == io::ErrorKind::NotFound => {}
        Err(error) => return Err(RsomicsError::Io(error)),
    }
    let name = path.file_name().ok_or_else(|| {
        RsomicsError::ConfigError(format!("output path has no file name: {}", path.display()))
    })?;
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    fs::canonicalize(parent)
        .map(|parent| parent.join(name))
        .map_err(RsomicsError::Io)
}

struct TransactionalFile<'a> {
    target: &'a Path,
    temporary: tempfile::NamedTempFile,
    permissions: Option<fs::Permissions>,
}

impl<'a> TransactionalFile<'a> {
    fn new(target: &'a Path) -> Result<Self> {
        let parent = target
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        let temporary = tempfile::NamedTempFile::new_in(parent).map_err(|error| {
            RsomicsError::Io(io::Error::new(
                error.kind(),
                format!(
                    "creating temporary output beside {}: {error}",
                    target.display()
                ),
            ))
        })?;
        let permissions = fs::metadata(target)
            .ok()
            .map(|metadata| metadata.permissions());
        Ok(Self {
            target,
            temporary,
            permissions,
        })
    }

    fn reopen(&self) -> Result<fs::File> {
        self.temporary.reopen().map_err(|error| {
            RsomicsError::Io(io::Error::new(
                error.kind(),
                format!(
                    "opening temporary output beside {}: {error}",
                    self.target.display()
                ),
            ))
        })
    }

    fn commit(mut self) -> Result<()> {
        if let Some(permissions) = self.permissions {
            self.temporary
                .as_file_mut()
                .set_permissions(permissions)
                .map_err(RsomicsError::Io)?;
        }
        self.temporary
            .as_file_mut()
            .sync_all()
            .map_err(RsomicsError::Io)?;
        self.temporary.persist(self.target).map_err(|error| {
            RsomicsError::Io(io::Error::new(
                error.error.kind(),
                format!(
                    "committing output {}: {}",
                    self.target.display(),
                    error.error
                ),
            ))
        })?;
        Ok(())
    }
}

fn parse_flags(value: &str) -> std::result::Result<u16, String> {
    flags::parse(value).map_err(|error| error.to_string())
}

fn parse_region(value: &str) -> std::result::Result<view::Region, String> {
    value.parse()
}

fn combine_flags(flags: &[u16]) -> u16 {
    flags
        .iter()
        .copied()
        .fold(0, |combined, flag| combined | flag)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn command_line_arguments_cannot_create_header_fields_or_lines() {
        assert_eq!(
            sanitize_argument(OsStr::new("input\tname\r\n.sam")),
            "input name  .sam"
        );
    }
}
