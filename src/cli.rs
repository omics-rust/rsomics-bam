use std::process;

use clap::{Parser, Subcommand};
use rsomics_common::{OutputArgs, Result, ToolMeta, run as run_tool};
use serde::Serialize;

use crate::{
    addreplacerg, ampliconclip, ampliconstats, bedcov, calmd, cat, collate, commands, coverage,
    depad, depth, fixmate, flags, flagstat, head, idxstats, index, markdup, merge, mpileup,
    quickcheck, reheader, reset, samples, sort, split, stats, to_bed, view,
};

const META: ToolMeta = ToolMeta {
    name: "rsomics-bam",
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-bam",
    version,
    about = "SAM, BAM, and CRAM workflows",
    arg_required_else_help = true
)]
struct Cli {
    #[command(flatten)]
    output: OutputArgs,

    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Add or replace alignment read groups
    Addreplacerg(commands::addreplacerg::Arguments),
    /// Clip amplicon primer sequence from alignments
    Ampliconclip(commands::ampliconclip::Arguments),
    /// Report amplicon sequencing statistics
    Ampliconstats(commands::ampliconstats::Arguments),
    /// Report coverage totals over BED regions
    Bedcov(commands::bedcov::Arguments),
    /// Recalculate alignment MD and NM tags
    Calmd(commands::calmd::Arguments),
    /// Concatenate BAM files without reencoding alignment blocks
    Cat(commands::cat::Arguments),
    /// Compute content checksums or merge prior reports
    Checksum(commands::checksum::Arguments),
    /// Group alignments by read name with bounded memory
    Collate(commands::collate::Arguments),
    /// Call a consensus sequence from coordinate-sorted alignments
    Consensus(commands::consensus::Arguments),
    /// Summarize coverage by reference sequence
    Coverage(commands::coverage::Arguments),
    /// Report CRAM data-series storage by block
    CramSize(commands::cram_size::Arguments),
    /// Compute read depth at each position
    Depth(commands::depth::Arguments),
    /// Convert padded-reference alignments to unpadded coordinates
    Depad(commands::depad::Arguments),
    /// Repair mate fields in name-grouped alignments
    Fixmate(commands::fixmate::Arguments),
    /// Convert between numeric and symbolic SAM flags
    Flags(commands::flags::Arguments),
    /// Count alignments by SAM flag category
    Flagstat(commands::flagstat::Arguments),
    /// Convert name-grouped alignments to FASTA
    Fasta(commands::fastx::FastaArguments),
    /// Convert name-grouped alignments to FASTQ
    Fastq(commands::fastx::FastqArguments),
    /// Print header lines and the first alignments as SAM
    Head(commands::head::Arguments),
    /// Build BAI, CSI, or CRAI random-access indexes
    Index(commands::index::Arguments),
    /// Report mapped and unmapped segments by reference
    Idxstats(commands::idxstats::Arguments),
    /// Convert FASTQ reads to unmapped SAM or BAM
    Import(commands::import::Arguments),
    /// Merge ordered alignment files
    Merge(commands::merge::Arguments),
    /// Mark duplicate alignments in coordinate order
    Markdup(commands::markdup::Arguments),
    /// Generate per-position text pileup
    Mpileup(commands::mpileup::Arguments),
    /// Call and phase heterozygous SNPs
    Phase(commands::phase::Arguments),
    /// Check alignment headers and format-specific end markers
    Quickcheck(commands::quickcheck::Arguments),
    /// Replace a BAM header without reencoding alignment blocks
    Reheader(commands::reheader::Arguments),
    /// Restore primary alignments to unaligned reads
    Reset(commands::reset::Arguments),
    /// List samples declared by alignment read groups
    Samples(commands::samples::Arguments),
    /// Sort alignments with bounded memory
    Sort(commands::sort::Arguments),
    /// Partition alignments into a complete transactional output set
    Split(commands::split::Arguments),
    /// Produce comprehensive alignment statistics
    Stats(Box<commands::stats::Arguments>),
    /// Convert alignments to BED intervals
    ToBed(commands::to_bed::Arguments),
    /// Filter and convert alignment records
    View(Box<commands::view::Arguments>),
}

#[derive(Debug, Serialize)]
#[serde(tag = "command", rename_all = "kebab-case")]
pub(crate) enum CommandOutput {
    Addreplacerg {
        summary: addreplacerg::Summary,
    },
    Ampliconclip {
        run: ampliconclip::Run,
    },
    Ampliconstats {
        summary: ampliconstats::Summary,
    },
    Bedcov {
        summary: bedcov::Summary,
    },
    Calmd {
        summary: calmd::Summary,
    },
    Cat {
        summary: cat::Summary,
    },
    Checksum {
        reports: Vec<crate::checksum::Report>,
    },
    Collate {
        summary: collate::Summary,
    },
    Consensus {
        summary: crate::consensus::Summary,
    },
    Coverage {
        report: coverage::Report,
    },
    CramSize {
        report: crate::cram_size::Report,
    },
    Depth {
        summary: depth::Summary,
    },
    Depad {
        summary: depad::Summary,
    },
    Fixmate {
        summary: fixmate::Summary,
    },
    Flags {
        values: Vec<flags::FlagValue>,
    },
    Flagstat {
        counts: Box<flagstat::Counts>,
    },
    Fasta {
        summary: crate::fastx::Summary,
    },
    Fastq {
        summary: crate::fastx::Summary,
    },
    Head {
        summary: head::Summary,
    },
    Index {
        summaries: Vec<index::Summary>,
    },
    Idxstats {
        report: idxstats::Report,
    },
    Import {
        summary: crate::import::Summary,
    },
    Merge {
        summary: merge::Summary,
    },
    Markdup {
        summary: markdup::Summary,
    },
    Mpileup {
        summary: mpileup::Summary,
    },
    Phase {
        summary: crate::phase::Summary,
    },
    Quickcheck {
        report: quickcheck::Report,
    },
    Reheader {
        summary: reheader::Summary,
    },
    Reset {
        summary: reset::Summary,
    },
    Samples {
        report: samples::Report,
    },
    Sort {
        summary: sort::Summary,
    },
    Split {
        summary: split::Summary,
    },
    Stats {
        report: Box<stats::Report>,
    },
    ToBed {
        summary: to_bed::Summary,
    },
    View {
        summary: view::Summary,
    },
}

#[must_use]
pub(crate) fn run() -> process::ExitCode {
    let cli = rsomics_help::parse::<Cli>();
    let output = cli.output.clone();
    run_tool(&output, META, || execute(cli))
}

fn execute(cli: Cli) -> Result<CommandOutput> {
    match cli.command {
        Command::Addreplacerg(arguments) => {
            commands::addreplacerg::execute(arguments, cli.output.json)
        }
        Command::Ampliconclip(arguments) => {
            commands::ampliconclip::execute(arguments, cli.output.json)
        }
        Command::Ampliconstats(arguments) => {
            commands::ampliconstats::execute(arguments, cli.output.json)
        }
        Command::Bedcov(arguments) => commands::bedcov::execute(arguments, cli.output.json),
        Command::Calmd(arguments) => commands::calmd::execute(arguments, cli.output.json),
        Command::Cat(arguments) => commands::cat::execute(arguments, cli.output.json),
        Command::Checksum(arguments) => commands::checksum::execute(arguments, cli.output.json),
        Command::Collate(arguments) => commands::collate::execute(arguments, cli.output.json),
        Command::Consensus(arguments) => commands::consensus::execute(arguments, cli.output.json),
        Command::Coverage(arguments) => commands::coverage::execute(arguments, cli.output.json),
        Command::CramSize(arguments) => commands::cram_size::execute(arguments, cli.output.json),
        Command::Depth(arguments) => commands::depth::execute(arguments, cli.output.json),
        Command::Depad(arguments) => commands::depad::execute(arguments, cli.output.json),
        Command::Fixmate(arguments) => commands::fixmate::execute(arguments, cli.output.json),
        Command::Flags(arguments) => commands::flags::execute(arguments, cli.output.json),
        Command::Flagstat(arguments) => commands::flagstat::execute(arguments, cli.output.json),
        Command::Fasta(arguments) => commands::fastx::execute_fasta(arguments, cli.output.json),
        Command::Fastq(arguments) => commands::fastx::execute_fastq(arguments, cli.output.json),
        Command::Head(arguments) => commands::head::execute(arguments, cli.output.json),
        Command::Index(arguments) => commands::index::execute(arguments, cli.output.json),
        Command::Idxstats(arguments) => commands::idxstats::execute(arguments, cli.output.json),
        Command::Import(arguments) => commands::import::execute(arguments, cli.output.json),
        Command::Merge(arguments) => commands::merge::execute(arguments, cli.output.json),
        Command::Markdup(arguments) => commands::markdup::execute(arguments, cli.output.json),
        Command::Mpileup(arguments) => commands::mpileup::execute(arguments, cli.output.json),
        Command::Phase(arguments) => commands::phase::execute(arguments, cli.output.json),
        Command::Quickcheck(arguments) => commands::quickcheck::execute(arguments, cli.output.json),
        Command::Reheader(arguments) => commands::reheader::execute(arguments, cli.output.json),
        Command::Reset(arguments) => commands::reset::execute(arguments, cli.output.json),
        Command::Samples(arguments) => commands::samples::execute(arguments, cli.output.json),
        Command::Sort(arguments) => commands::sort::execute(arguments, cli.output.json),
        Command::Split(arguments) => commands::split::execute(arguments, cli.output.json),
        Command::Stats(arguments) => commands::stats::execute(*arguments, cli.output.json),
        Command::ToBed(arguments) => commands::to_bed::execute(arguments, cli.output.json),
        Command::View(arguments) => commands::view::execute(*arguments, cli.output.json),
    }
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    use super::*;

    #[test]
    fn command_tree_is_valid() {
        Cli::command().debug_assert();
    }

    #[test]
    fn consensus_help_exposes_complete_stable_surface() {
        let help =
            rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "consensus", "--help"])
                .unwrap_err()
                .to_string();
        for option in [
            "-m, --mode",
            "--regions-file",
            "-f, --format",
            "--show-del",
            "--show-ins",
            "--mark-ins",
            "--min-MQ",
            "--min-BQ",
            "-c, --call-fract",
            "-C, --cutoff",
            "--no-adj-qual",
            "--no-use-MQ",
            "--no-adj-MQ",
            "--NM-halo",
            "--scale-MQ",
            "--low-MQ",
            "--high-MQ",
            "--P-het",
            "--P-indel",
            "--het-scale",
            "-t, --qual-calibration",
            "-X, --config",
            "-T, --reference",
            "-@, --threads",
        ] {
            assert!(help.contains(option), "missing {option} in {help}");
        }
        for excluded in [
            "--SC-cost",
            "--default-qual",
            "--het-only",
            "--homopoly-redux",
            "--block-size",
        ] {
            assert!(!help.contains(excluded), "unexpected {excluded} in {help}");
        }
    }

    #[test]
    fn phase_help_exposes_the_stable_slice() {
        let help = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "phase", "--help"])
            .unwrap_err()
            .to_string();
        for option in [
            "-o, --output",
            "-b, --output-prefix",
            "-O, --output-fmt",
            "-k <INT>",
            "-q <INT>",
            "-Q, --min-BQ",
            "-D <INT>",
            "-F",
            "-A",
            "--reference",
            "--no-pg",
            "-@, --threads",
        ] {
            assert!(help.contains(option), "missing {option} in {help}");
        }
        for excluded in ["-l", "-e", "--input-fmt-option", "--output-fmt-option"] {
            assert!(!help.contains(excluded), "unexpected {excluded} in {help}");
        }
    }

    #[test]
    fn split_help_exposes_one_unified_partition_surface() {
        let help = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "split", "--help"])
            .unwrap_err()
            .to_string();
        for option in [
            "--tag",
            "--parts",
            "--genes",
            "--mates",
            "-b, --output-prefix",
            "-u, --unaccounted",
            "--unaccounted-header",
            "-O, --output-fmt",
            "-M, --max-outputs",
            "--zero-pad",
            "--seed",
            "--skip-unmapped",
            "-T, --reference",
            "-@, --threads",
            "--no-pg",
        ] {
            assert!(help.contains(option), "missing {option} in {help}");
        }
        assert!(!help.contains("--filename-format"));
    }

    #[test]
    fn split_mode_selectors_are_mutually_exclusive() {
        assert!(
            rsomics_help::try_parse_from::<Cli, _, _>([
                "rsomics-bam",
                "split",
                "--tag",
                "RG",
                "--parts",
                "2",
                "--output-prefix",
                "out",
            ])
            .is_err()
        );
        assert!(
            rsomics_help::try_parse_from::<Cli, _, _>([
                "rsomics-bam",
                "split",
                "--skip-unmapped",
                "--output-prefix",
                "out",
            ])
            .is_err()
        );
    }

    #[test]
    fn consensus_rejects_experimental_options() {
        for arguments in [
            ["--SC-cost", "20"],
            ["--default-qual", "20"],
            ["--het-only", "yes"],
            ["--homopoly-redux", "0.1"],
            ["--block-size", "1000"],
        ] {
            assert!(
                rsomics_help::try_parse_from::<Cli, _, _>([
                    "rsomics-bam",
                    "consensus",
                    arguments[0],
                    arguments[1],
                ])
                .is_err(),
                "accepted {}",
                arguments[0]
            );
        }
    }

    #[test]
    fn consensus_boolean_options_require_yes_or_no_values() {
        assert!(
            rsomics_help::try_parse_from::<Cli, _, _>([
                "rsomics-bam",
                "consensus",
                "--show-del",
                "yes",
                "--show-ins",
                "no",
            ])
            .is_ok()
        );
        assert!(
            rsomics_help::try_parse_from::<Cli, _, _>([
                "rsomics-bam",
                "consensus",
                "--show-del",
                "true",
            ])
            .is_err()
        );
    }

    #[test]
    fn consensus_command_writes_the_samtools_oracle_and_summary() {
        let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/upstream/samtools-consensus");
        let directory = tempfile::tempdir().unwrap();
        let output = directory.path().join("consensus.fastq");
        let cli = rsomics_help::try_parse_from::<Cli, _, _>([
            "rsomics-bam".into(),
            "consensus".into(),
            "--mode".into(),
            "simple".into(),
            "--call-fract".into(),
            "0.6".into(),
            "--format".into(),
            "fastq".into(),
            "--output".into(),
            output.as_os_str().to_owned(),
            root.join("consen1.sam").into_os_string(),
        ])
        .unwrap();

        let result = execute(cli).unwrap();

        let CommandOutput::Consensus { summary } = result else {
            panic!("unexpected command output")
        };
        assert_eq!(summary.sequence_records, 1);
        assert_eq!(summary.sequence_symbols, 16);
        assert_eq!(summary.pileup_rows, 0);
        assert_eq!(
            std::fs::read(output).unwrap(),
            std::fs::read(root.join("expected/1q.out")).unwrap()
        );
    }

    #[test]
    fn cram_size_help_exposes_the_complete_diagnostic_surface() {
        let help =
            rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "cram-size", "--help"])
                .unwrap_err()
                .to_string();
        for option in ["-o, --output", "-v, --verbose", "-e, --encodings"] {
            assert!(help.contains(option), "missing {option} in {help}");
        }
    }

    #[test]
    fn help_uses_one_nested_command_tree() {
        let error =
            rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "flagstat", "--help"])
                .unwrap_err();
        let help = error.to_string();
        assert!(help.contains("Input SAM, BAM, or CRAM file"), "{help}");
        assert!(help.contains("--output-fmt"), "{help}");
        assert!(help.contains("-@, --threads <INT>"), "{help}");
    }

    #[test]
    fn view_exposes_program_provenance_control() {
        let command = Cli::command();
        let view = command
            .get_subcommands()
            .find(|command| command.get_name() == "view")
            .unwrap();
        let argument = view
            .get_arguments()
            .find(|argument| argument.get_long() == Some("no-pg"))
            .unwrap();
        assert!(!argument.is_hide_set());
        assert!(
            rsomics_help::try_parse_from::<Cli, _, _>([
                "rsomics-bam",
                "view",
                "--no-PG",
                "input.sam"
            ])
            .is_ok()
        );
    }

    #[test]
    fn markdup_help_exposes_stable_scope() {
        let error = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "markdup", "--help"])
            .unwrap_err();
        let help = error.to_string();
        for option in [
            "-r, --remove",
            "-c, --clear",
            "--include-fails",
            "-m, --mode",
            "-l, --max-read-length",
            "-@, --threads",
            "--reference",
            "--no-pg",
        ] {
            assert!(help.contains(option), "missing {option} in {help}");
        }
        for excluded in ["--barcode-tag", "--duplicate-count", "--read-groups"] {
            assert!(!help.contains(excluded), "unexpected {excluded} in {help}");
        }
    }

    #[test]
    fn cat_and_reheader_help_expose_only_the_stable_slice() {
        let cat = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "cat", "--help"])
            .unwrap_err()
            .to_string();
        for option in ["-b, --list", "--header", "-o, --output", "--no-pg"] {
            assert!(cat.contains(option), "missing {option} in {cat}");
        }
        for excluded in ["--region", "--part", "--fast", "--query"] {
            assert!(!cat.contains(excluded), "unexpected {excluded} in {cat}");
        }

        let reheader =
            rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "reheader", "--help"])
                .unwrap_err()
                .to_string();
        assert!(reheader.contains("-o, --output"), "{reheader}");
        assert!(reheader.contains("--no-pg"), "{reheader}");
        assert!(!reheader.contains("--command"), "{reheader}");
        assert!(!reheader.contains("--in-place"), "{reheader}");
    }

    #[test]
    fn reset_help_exposes_the_stable_slice() {
        let help = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "reset", "--help"])
            .unwrap_err()
            .to_string();
        for option in [
            "-o, --output",
            "-O, --output-fmt",
            "-x, --remove-tag",
            "--keep-tag",
            "--no-RG",
            "--reject-PG",
            "--dupflag",
            "-T, --reference",
            "-@, --threads",
            "--no-pg",
        ] {
            assert!(help.contains(option), "missing {option} in {help}");
        }
        for excluded in ["--input-fmt-option", "--output-fmt-option"] {
            assert!(!help.contains(excluded), "unexpected {excluded} in {help}");
        }
    }

    #[test]
    fn addreplacerg_help_exposes_the_stable_slice() {
        let help =
            rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "addreplacerg", "--help"])
                .unwrap_err()
                .to_string();
        for option in [
            "-r, --rg-line",
            "-R, --rg-id",
            "-m, --mode",
            "-w, --overwrite-header",
            "-o, --output",
            "-O, --output-fmt",
            "-u, --uncompressed",
            "--reference",
            "-@, --threads",
            "--no-pg",
        ] {
            assert!(help.contains(option), "missing {option} in {help}");
        }
        for excluded in ["--write-index", "--output-fmt-option", "--verbosity"] {
            assert!(!help.contains(excluded), "unexpected {excluded} in {help}");
        }
    }

    #[test]
    fn coverage_summary_help_exposes_only_the_stable_table_slice() {
        let bedcov = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "bedcov", "--help"])
            .unwrap_err()
            .to_string();
        for option in [
            "-Q, --min-MQ",
            "-X",
            "-g <FLAG>",
            "-G <FLAG>",
            "-j",
            "-d <INT>",
            "--max-depth",
            "-c",
            "-H",
            "-o, --output",
            "-@, --threads",
            "--reference",
        ] {
            assert!(bedcov.contains(option), "missing {option} in {bedcov}");
        }

        let coverage =
            rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "coverage", "--help"])
                .unwrap_err()
                .to_string();
        for option in [
            "-b, --bam-list",
            "-l, --min-read-len",
            "-q, --min-MQ",
            "-Q, --min-BQ",
            "--rf",
            "--ff",
            "-d, --depth",
            "--min-depth",
            "-r, --region",
            "-H, --no-header",
            "-o, --output",
            "-@, --threads",
            "--reference",
        ] {
            assert!(coverage.contains(option), "missing {option} in {coverage}");
        }
        for excluded in ["--histogram", "--plot-depth", "--ascii", "--n-bins"] {
            assert!(
                !coverage.contains(excluded),
                "unexpected {excluded} in {coverage}"
            );
        }

        let idxstats =
            rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "idxstats", "--help"])
                .unwrap_err()
                .to_string();
        for option in ["-X", "-o, --output", "-@, --threads", "--reference"] {
            assert!(idxstats.contains(option), "missing {option} in {idxstats}");
        }
    }

    #[test]
    fn amplicon_help_exposes_the_paired_workflow() {
        let clip =
            rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "ampliconclip", "--help"])
                .unwrap_err()
                .to_string();
        for option in [
            "-b <BED>",
            "-o, --output",
            "-f, --stats",
            "--soft-clip",
            "--hard-clip",
            "--both-ends",
            "--strand",
            "--rejects-file",
            "--primer-counts",
            "--original",
            "--keep-tag",
            "--no-pg",
            "-@, --threads",
        ] {
            assert!(clip.contains(option), "missing {option} in {clip}");
        }

        let stats =
            rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "ampliconstats", "--help"])
                .unwrap_err()
                .to_string();
        for option in [
            "-f, --required-flag",
            "-F, --filter-flag",
            "-a, --max-amplicons",
            "-l, --max-amplicon-length",
            "-d, --min-depth",
            "-m, --pos-margin",
            "-o, --output",
            "-s, --use-sample-name",
            "-t, --tlen-adjust",
            "-b, --tcoord-bin",
            "-c, --tcoord-min-count",
            "-D, --depth-bin",
            "-S, --single-ref",
            "-@, --threads",
        ] {
            assert!(stats.contains(option), "missing {option} in {stats}");
        }
    }

    #[test]
    fn calmd_help_exposes_only_the_stable_recalculation_slice() {
        let help = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "calmd", "--help"])
            .unwrap_err()
            .to_string();
        for option in [
            "-e, --convert-equal",
            "-b, --bam",
            "-u, --uncompressed",
            "-O, --output-fmt",
            "-o, --output",
            "-@, --threads",
            "--no-pg",
        ] {
            assert!(help.contains(option), "missing {option} in {help}");
        }
        for excluded in ["--baq", "--extended-baq", "--cap-mapq", "--output-cram"] {
            assert!(!help.contains(excluded), "unexpected {excluded} in {help}");
        }
    }

    #[test]
    fn depad_help_exposes_only_the_stable_projection_slice() {
        let help = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "depad", "--help"])
            .unwrap_err()
            .to_string();
        for option in [
            "-s, --sam",
            "-S",
            "-u, --uncompressed",
            "-1, --fast-compression",
            "-O, --output-fmt",
            "-T, --reference",
            "-o, --output",
            "-@, --threads",
            "--no-pg",
        ] {
            assert!(help.contains(option), "missing {option} in {help}");
        }
        for excluded in ["--write-index", "--output-fmt-option", "--output-cram"] {
            assert!(!help.contains(excluded), "unexpected {excluded} in {help}");
        }
    }

    #[test]
    fn fasta_and_fastq_help_expose_the_unified_output_slice() {
        let fasta = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "fasta", "--help"])
            .unwrap_err()
            .to_string();
        for option in [
            "-o, --output",
            "-n, --no-mate-suffix",
            "-f, --require-flags",
            "-F, --exclude-flags",
            "--rf",
            "-G",
            "-c, --compression-level",
            "--reference",
            "-@, --threads",
        ] {
            assert!(fasta.contains(option), "missing {option} in {fasta}");
        }
        for excluded in ["-0", "-1", "-2", "--singleton", "--tag", "--UMI"] {
            assert!(
                !fasta.contains(excluded),
                "unexpected {excluded} in {fasta}"
            );
        }

        let fastq = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "fastq", "--help"])
            .unwrap_err()
            .to_string();
        assert!(fastq.contains("-O, --original-quality"), "{fastq}");
        assert!(fastq.contains("-v, --default-quality"), "{fastq}");
    }

    #[test]
    fn import_help_exposes_the_stable_fastq_slice() {
        let help = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "import", "--help"])
            .unwrap_err()
            .to_string();
        for option in [
            "-0, --single",
            "-s, --interleaved",
            "-1, --read1",
            "-2, --read2",
            "-o, --output",
            "-O, --output-fmt",
            "-r, --rg-line",
            "-R, --rg",
            "--order",
            "-@, --threads",
            "--no-pg",
        ] {
            assert!(help.contains(option), "missing {option} in {help}");
        }
        for excluded in ["--i1", "--CASAVA", "--UMI", "--name2", "--barcode-tag"] {
            assert!(!help.contains(excluded), "unexpected {excluded} in {help}");
        }
    }

    #[test]
    fn to_bed_help_exposes_the_complete_conversion_surface() {
        let help = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-bam", "to-bed", "--help"])
            .unwrap_err()
            .to_string();
        for option in [
            "--bedpe",
            "--mate1",
            "--ed",
            "--tag",
            "--cigar",
            "--bed12",
            "--split",
            "--split-d",
            "--color",
            "-T, --reference",
            "-@, --threads",
            "-o, --output",
            "-i, --input",
        ] {
            assert!(help.contains(option), "missing {option} in {help}");
        }
        assert!(help.contains("Input SAM, BAM, or CRAM"), "{help}");
        assert!(help.contains("Global options:"), "{help}");
    }
}
