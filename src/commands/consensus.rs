use std::path::{Path, PathBuf};

use clap::{ArgAction, Args, ValueEnum};
use rsomics_common::{Result, RsomicsError, reject_output_alias, write_output};

use crate::cli::CommandOutput;
use crate::consensus::{
    BayesianOverrides, CalibrationPreset, Format, Options, Profile, write_pileup,
};
use crate::flags;

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Mode {
    Simple,
    Bayesian,
    #[value(name = "bayesian-116", alias = "bayesian_116")]
    Bayesian116,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum OutputFormat {
    Fasta,
    Fastq,
    Pileup,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum MachineProfile {
    Hifi,
    Hiseq,
    #[value(name = "r10.4-sup", alias = "r10.4_sup")]
    R10_4Sup,
    #[value(name = "r10.4-dup", alias = "r10.4_dup")]
    R10_4Dup,
    Ultima,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum YesNo {
    Yes,
    No,
}

impl YesNo {
    fn enabled(self) -> bool {
        matches!(self, Self::Yes)
    }
}

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Coordinate-sorted SAM, BAM, or CRAM; omit or use - for standard input
    #[arg(value_name = "ALIGNMENT")]
    input: Option<PathBuf>,

    /// Write consensus output to FILE
    #[arg(short = 'o', long, value_name = "FILE")]
    output: Option<PathBuf>,

    /// Output encoding
    #[arg(short = 'f', long, value_enum, default_value_t = OutputFormat::Fasta)]
    format: OutputFormat,

    /// Consensus model
    #[arg(short = 'm', long, value_enum, default_value_t = Mode::Bayesian)]
    mode: Mode,

    /// Restrict an indexed input to REGION
    #[arg(
        short = 'r',
        long,
        value_name = "REGION",
        conflicts_with = "regions_file"
    )]
    region: Option<String>,

    /// Process BED intervals independently without merging overlaps
    #[arg(long, value_name = "BED", conflicts_with = "region")]
    regions_file: Option<PathBuf>,

    /// Indexed reference FASTA for CRAM decoding and uncovered bases
    #[arg(short = 'T', long, value_name = "FASTA")]
    reference: Option<PathBuf>,

    /// Wrap FASTA and FASTQ records at this width
    #[arg(short = 'l', long = "line-len", value_name = "INT", default_value_t = 70, value_parser = parse_positive_usize)]
    line_width: usize,

    /// Output uncovered positions; repeat for unused references
    #[arg(short = 'a', action = ArgAction::Count)]
    all_positions: u8,

    /// Include deletion calls
    #[arg(long = "show-del", visible_alias = "show-deletions", value_name = "BOOL", value_enum, default_value_t = YesNo::No)]
    show_deletions: YesNo,

    /// Include insertion calls
    #[arg(long = "show-ins", visible_alias = "show-insertions", value_name = "BOOL", value_enum, default_value_t = YesNo::Yes)]
    show_insertions: YesNo,

    /// Prefix inserted sequence symbols with an underscore
    #[arg(long = "mark-ins", visible_alias = "mark-insertions")]
    mark_insertions: bool,

    /// Minimum alignment mapping quality
    #[arg(
        long = "min-MQ",
        visible_alias = "min-mq",
        value_name = "INT",
        default_value_t = 0
    )]
    minimum_mapping_quality: u8,

    /// Minimum observed base quality
    #[arg(
        long = "min-BQ",
        visible_alias = "min-bq",
        value_name = "INT",
        default_value_t = 0
    )]
    minimum_base_quality: u8,

    /// Minimum observations required to call a base
    #[arg(short = 'd', long = "min-depth", value_name = "INT", default_value_t = 1, value_parser = parse_positive_usize)]
    minimum_depth: usize,

    /// Enable IUPAC ambiguity calls
    #[arg(short = 'A', long = "ambig")]
    ambiguous: bool,

    /// Require reads with any listed FLAG bit
    #[arg(long = "rf", visible_alias = "incl-flags", value_name = "FLAG", default_value = "0", value_parser = parse_flags)]
    required_flags: u16,

    /// Exclude reads with any listed FLAG bit
    #[arg(long = "ff", visible_alias = "excl-flags", value_name = "FLAG", default_value = "UNMAP,SECONDARY,QCFAIL,DUP", value_parser = parse_flags)]
    excluded_flags: u16,

    /// Additional alignment decompression threads
    #[arg(short = '@', long, value_name = "INT", default_value_t = 0)]
    threads: usize,

    /// Quality assigned to uncovered reference bases
    #[arg(long = "ref-qual", value_name = "INT", default_value_t = 0, value_parser = clap::value_parser!(u8).range(0..=93))]
    reference_quality: u8,

    /// Weight simple-mode votes by base quality
    #[arg(short = 'q', long = "use-qual", conflicts_with = "no_use_quality")]
    use_quality: bool,

    /// Ignore base quality in simple mode
    #[arg(long = "no-use-qual", conflicts_with = "use_quality")]
    no_use_quality: bool,

    /// Required fraction supporting the simple-mode call
    #[arg(short = 'c', long = "call-fract", value_name = "FLOAT", value_parser = parse_fraction)]
    call_fraction: Option<f64>,

    /// Required second-to-first support ratio for ambiguity [default: 0.5]
    #[arg(short = 'H', long = "het-fract", value_name = "FLOAT", value_parser = parse_fraction)]
    heterozygous_fraction: Option<f64>,

    /// Bayesian call-quality cutoff
    #[arg(short = 'C', long, value_name = "INT")]
    cutoff: Option<i32>,

    /// Disable local-minimum base-quality adjustment
    #[arg(long = "no-adj-qual", conflicts_with = "adjust_quality")]
    no_adjust_quality: bool,

    /// Enable local-minimum base-quality adjustment
    #[arg(long = "adj-qual", conflicts_with = "no_adjust_quality")]
    adjust_quality: bool,

    /// Ignore mapping quality in the Bayesian model
    #[arg(
        long = "no-use-MQ",
        visible_alias = "no-use-mq",
        conflicts_with = "use_mapping_quality"
    )]
    no_use_mapping_quality: bool,

    /// Use mapping quality in the Bayesian model
    #[arg(
        long = "use-MQ",
        visible_alias = "use-mq",
        conflicts_with = "no_use_mapping_quality"
    )]
    use_mapping_quality: bool,

    /// Disable local-NM mapping-quality adjustment
    #[arg(
        long = "no-adj-MQ",
        visible_alias = "no-adj-mq",
        conflicts_with = "adjust_mapping_quality"
    )]
    no_adjust_mapping_quality: bool,

    /// Enable local-NM mapping-quality adjustment
    #[arg(
        long = "adj-MQ",
        visible_alias = "adj-mq",
        conflicts_with = "no_adjust_mapping_quality"
    )]
    adjust_mapping_quality: bool,

    /// Window radius for local NM adjustment
    #[arg(long = "NM-halo", visible_alias = "nm-halo", value_name = "INT")]
    mismatch_halo: Option<usize>,

    /// Scale mapping quality inside the Bayesian model
    #[arg(long = "scale-MQ", visible_alias = "scale-mq", value_name = "FLOAT", value_parser = parse_positive_f64)]
    mapping_quality_scale: Option<f64>,

    /// Lower mapping-quality cap inside the Bayesian model
    #[arg(long = "low-MQ", visible_alias = "low-mq", value_name = "INT")]
    low_mapping_quality: Option<u8>,

    /// Upper mapping-quality cap inside the Bayesian model
    #[arg(long = "high-MQ", visible_alias = "high-mq", value_name = "INT")]
    high_mapping_quality: Option<u8>,

    /// Prior probability of a heterozygous site
    #[arg(long = "P-het", visible_alias = "p-het", value_name = "FLOAT", value_parser = parse_probability)]
    heterozygous_probability: Option<f64>,

    /// Prior probability of an indel site
    #[arg(long = "P-indel", visible_alias = "p-indel", value_name = "FLOAT", value_parser = parse_probability)]
    indel_probability: Option<f64>,

    /// Heterozygous SNP probability multiplier
    #[arg(long = "het-scale", value_name = "FLOAT", value_parser = parse_positive_f64)]
    heterozygous_scale: Option<f64>,

    /// Redistribute low qualities across homopolymer ends
    #[arg(short = 'p', long = "homopoly-fix", visible_alias = "homopoly-score", value_name = "FRACTION", num_args = 0..=1, default_missing_value = "0.5", value_parser = parse_fraction)]
    homopolymer_fix: Option<f64>,

    /// Quality calibration file or :flat, :hifi, :hiseq, :r10.4_sup, :r10.4_dup, :ultima
    #[arg(short = 't', long = "qual-calibration", value_name = "FILE|PRESET")]
    calibration: Option<String>,

    /// Sequencing-platform parameter profile
    #[arg(short = 'X', long = "config", value_enum)]
    profile: Option<MachineProfile>,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let input = arguments
        .input
        .clone()
        .unwrap_or_else(|| PathBuf::from("-"));
    if json
        && arguments
            .output
            .as_deref()
            .is_none_or(|path| path == Path::new("-"))
    {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output for consensus".to_owned(),
        ));
    }
    let mut sources = vec![input.as_path()];
    sources.extend(arguments.reference.as_deref());
    sources.extend(arguments.regions_file.as_deref());
    if let Some(calibration) = arguments
        .calibration
        .as_deref()
        .filter(|value| !value.starts_with(':'))
    {
        sources.push(Path::new(calibration));
    }
    if let Some(output) = arguments.output.as_deref() {
        reject_output_alias(output, sources)?;
    }
    if input == Path::new("-") && (arguments.region.is_some() || arguments.regions_file.is_some()) {
        return Err(RsomicsError::ConfigError(
            "indexed region selection requires a named alignment input".to_owned(),
        ));
    }
    if arguments.reference_quality > 0 && arguments.reference.is_none() {
        return Err(RsomicsError::ConfigError(
            "--ref-qual requires --reference".to_owned(),
        ));
    }
    if arguments.mark_insertions && matches!(arguments.format, OutputFormat::Pileup) {
        return Err(RsomicsError::ConfigError(
            "--mark-ins is only valid for FASTA and FASTQ output".to_owned(),
        ));
    }
    if arguments.all_positions > 2 {
        return Err(RsomicsError::ConfigError(
            "-a can be supplied at most twice".to_owned(),
        ));
    }

    let mut options = match arguments.mode {
        Mode::Simple => build_simple(&arguments)?,
        Mode::Bayesian | Mode::Bayesian116 => build_bayesian(&arguments)?,
    };
    options.minimum_mapping_quality = arguments.minimum_mapping_quality;
    options.excluded_flags = arguments.excluded_flags;
    options.required_flags = arguments.required_flags;
    options.format = match arguments.format {
        OutputFormat::Fasta => Format::Fasta,
        OutputFormat::Fastq => Format::Fastq,
        OutputFormat::Pileup => Format::Pileup,
    };
    options.show_deletions = arguments.show_deletions.enabled();
    options.show_insertions = arguments.show_insertions.enabled();
    options.mark_insertions = arguments.mark_insertions;
    options.all_positions = arguments.all_positions;
    options.reference = arguments.reference;
    options.reference_quality = arguments.reference_quality;
    options.region = arguments.region;
    options.regions_file = arguments.regions_file;
    options.additional_threads = arguments.threads;
    options.line_width = arguments.line_width;

    let summary = write_output(arguments.output.as_deref(), |output| {
        write_pileup(&input, options, output)
    })?;
    Ok(CommandOutput::Consensus { summary })
}

fn build_simple(arguments: &Arguments) -> Result<Options> {
    if has_bayesian_arguments(arguments) {
        return Err(RsomicsError::ConfigError(
            "Bayesian-only options cannot be used with --mode simple".to_owned(),
        ));
    }
    let mut options = Options::simple(arguments.call_fraction.unwrap_or(0.75));
    options.configure_simple(
        arguments.use_quality && !arguments.no_use_quality,
        arguments.minimum_base_quality,
        arguments.minimum_depth,
        arguments.heterozygous_fraction.unwrap_or(0.5),
        arguments.ambiguous,
    )?;
    Ok(options)
}

fn build_bayesian(arguments: &Arguments) -> Result<Options> {
    if arguments.use_quality
        || arguments.no_use_quality
        || arguments.call_fraction.is_some()
        || arguments.heterozygous_fraction.is_some()
    {
        return Err(RsomicsError::ConfigError(
            "simple-only options cannot be used with a Bayesian mode".to_owned(),
        ));
    }
    let compatibility_116 = matches!(arguments.mode, Mode::Bayesian116);
    if compatibility_116 && arguments.profile.is_some() {
        return Err(RsomicsError::ConfigError(
            "machine profiles cannot be combined with bayesian-116".to_owned(),
        ));
    }
    let mut options = Options::bayesian(arguments.cutoff.unwrap_or(10), arguments.ambiguous);
    if let Some(profile) = arguments.profile {
        options.apply_profile(match profile {
            MachineProfile::Hifi => Profile::Hifi,
            MachineProfile::Hiseq => Profile::Hiseq,
            MachineProfile::R10_4Sup => Profile::R10_4Sup,
            MachineProfile::R10_4Dup => Profile::R10_4Dup,
            MachineProfile::Ultima => Profile::Ultima,
        })?;
    }
    options.configure_bayesian(
        arguments.minimum_base_quality,
        arguments.minimum_depth,
        arguments.ambiguous,
        compatibility_116,
    )?;
    options.apply_bayesian_overrides(BayesianOverrides {
        adjust_quality: arguments.adjust_quality || !arguments.no_adjust_quality,
        use_mapping_quality: arguments.use_mapping_quality || !arguments.no_use_mapping_quality,
        adjust_mapping_quality: arguments.adjust_mapping_quality
            || !arguments.no_adjust_mapping_quality,
        mismatch_halo: arguments.mismatch_halo.unwrap_or(50),
        soft_clip_cost: 60,
        mapping_quality_scale: arguments.mapping_quality_scale,
        minimum_mapping_quality: arguments.low_mapping_quality,
        maximum_mapping_quality: arguments.high_mapping_quality.unwrap_or(60),
        default_quality: 10,
        heterozygous_probability: arguments.heterozygous_probability.unwrap_or(1e-3),
        indel_probability: arguments.indel_probability.unwrap_or(2e-4),
        heterozygous_scale: arguments.heterozygous_scale,
        homopolymer_fix: arguments.homopolymer_fix,
    })?;
    if let Some(value) = arguments.calibration.as_deref() {
        match value {
            ":flat" => options.set_calibration_preset(CalibrationPreset::Flat)?,
            ":hifi" => options.set_calibration_preset(CalibrationPreset::Hifi)?,
            ":hiseq" => options.set_calibration_preset(CalibrationPreset::Hiseq)?,
            ":r10.4_sup" => options.set_calibration_preset(CalibrationPreset::R10_4Sup)?,
            ":r10.4_dup" => options.set_calibration_preset(CalibrationPreset::R10_4Dup)?,
            ":ultima" => options.set_calibration_preset(CalibrationPreset::Ultima)?,
            value if value.starts_with(':') => {
                return Err(RsomicsError::ConfigError(format!(
                    "unknown quality calibration preset: {value}"
                )));
            }
            path => options.set_calibration_file(Path::new(path))?,
        }
    }
    Ok(options)
}

fn has_bayesian_arguments(arguments: &Arguments) -> bool {
    arguments.cutoff.is_some()
        || arguments.no_adjust_quality
        || arguments.adjust_quality
        || arguments.no_use_mapping_quality
        || arguments.use_mapping_quality
        || arguments.no_adjust_mapping_quality
        || arguments.adjust_mapping_quality
        || arguments.mismatch_halo.is_some()
        || arguments.mapping_quality_scale.is_some()
        || arguments.low_mapping_quality.is_some()
        || arguments.high_mapping_quality.is_some()
        || arguments.heterozygous_probability.is_some()
        || arguments.indel_probability.is_some()
        || arguments.heterozygous_scale.is_some()
        || arguments.homopolymer_fix.is_some()
        || arguments.calibration.is_some()
        || arguments.profile.is_some()
}

fn parse_flags(value: &str) -> std::result::Result<u16, String> {
    flags::parse(value).map_err(|error| error.to_string())
}

fn parse_positive_usize(value: &str) -> std::result::Result<usize, String> {
    value
        .parse::<usize>()
        .ok()
        .filter(|&value| value > 0)
        .ok_or_else(|| "value must be greater than zero".to_owned())
}

fn parse_fraction(value: &str) -> std::result::Result<f64, String> {
    value
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite() && (0.0..=1.0).contains(value))
        .ok_or_else(|| "value must be a finite number from 0 to 1".to_owned())
}

fn parse_probability(value: &str) -> std::result::Result<f64, String> {
    value
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite() && *value > 0.0 && *value < 1.0)
        .ok_or_else(|| "value must be a finite number between 0 and 1".to_owned())
}

fn parse_positive_f64(value: &str) -> std::result::Result<f64, String> {
    value
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite() && *value > 0.0)
        .ok_or_else(|| "value must be a finite positive number".to_owned())
}
