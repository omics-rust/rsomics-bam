use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::output::{TransactionalFile, same_target};
use crate::{Program, hts_quickcheck, merge};

#[derive(Debug, Args)]
#[command(after_help = "\
Examples:
  rsomics-bam merge part1.bam part2.bam -o merged.bam
  rsomics-bam merge -n lane1.bam lane2.bam -o names.bam
  rsomics-bam merge --template-coordinate a.bam b.bam -o templates.bam")]
pub(crate) struct Arguments {
    /// Ordered input SAM, BAM, or CRAM files
    #[arg(required = true, value_name = "ALIGNMENT", num_args = 1..)]
    inputs: Vec<PathBuf>,

    /// Write BAM output to FILE; omit or use - for standard output
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Inputs use natural numeric read-name order
    #[arg(short = 'n', long, conflicts_with_all = ["ascii_name", "template_coordinate"])]
    natural_name: bool,

    /// Inputs use bytewise lexicographical read-name order
    #[arg(short = 'N', long, conflicts_with_all = ["natural_name", "template_coordinate"])]
    ascii_name: bool,

    /// Inputs use unclipped template-coordinate order
    #[arg(long, conflicts_with_all = ["natural_name", "ascii_name"])]
    template_coordinate: bool,

    /// Keep the first @RG record for each colliding ID
    #[arg(short = 'c', long)]
    combine_read_groups: bool,

    /// Keep the first @PG record for each colliding ID
    #[arg(short = 'p', long)]
    combine_programs: bool,

    /// Reference FASTA for CRAM decoding
    #[arg(long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Additional BAM output workers; omit for automatic parallelism
    #[arg(short = '@', long, value_name = "INT")]
    threads: Option<usize>,

    /// Do not add an @PG line
    #[arg(long = "no-pg", visible_alias = "no-PG")]
    suppress_program_record: bool,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let output = arguments
        .output
        .as_deref()
        .filter(|path| *path != Path::new("-"));
    if json && output.is_none() {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for merge".to_owned(),
        ));
    }
    if let Some(output) = output {
        for input in &arguments.inputs {
            if same_target(input, output)? {
                return Err(RsomicsError::ConfigError(format!(
                    "merge input and output must be different files: {}",
                    input.display()
                )));
            }
        }
    }

    let order = if arguments.natural_name {
        merge::Order::QueryNameNatural
    } else if arguments.ascii_name {
        merge::Order::QueryNameLexicographical
    } else if arguments.template_coordinate {
        merge::Order::TemplateCoordinate
    } else {
        merge::Order::Coordinate
    };
    let command_line = (!arguments.suppress_program_record).then(crate::program::command_line);
    let program = command_line
        .as_deref()
        .map(|line| Program::new("rsomics-bam", env!("CARGO_PKG_VERSION"), line))
        .transpose()?;
    let options = merge::Options {
        order,
        additional_threads: arguments.threads,
        reference: arguments.reference.as_deref(),
        destination: output,
        combine_read_groups: arguments.combine_read_groups,
        combine_programs: arguments.combine_programs,
        program,
    };

    let summary = match output {
        Some(path) => {
            let transaction = TransactionalFile::new(path)?;
            let file = transaction.reopen()?;
            let summary = merge::write(&arguments.inputs, options, BufWriter::new(file))?;
            hts_quickcheck::require_bgzf_eof(transaction.temporary_path())?;
            transaction.commit()?;
            summary
        }
        None => merge::write(&arguments.inputs, options, io::stdout())?,
    };
    Ok(CommandOutput::Merge { summary })
}
