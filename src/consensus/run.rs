use std::io::Write;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::Arc;

use noodles::core::Region;
use noodles::fasta::record::Sequence;
use rsomics_common::{Result, RsomicsError};
use rsomics_pileup::{FlagFilter, PileupEngine, PileupOptions, RecordFilter};

use super::{
    call::{
        BayesianCaller, BayesianMode, BayesianOptions, CalibrationPreset, Caller, SimpleOptions,
    },
    output::{Configuration as OutputConfiguration, Format, Output, Reference, Summary},
    record::{RecordOptions, RecordState},
    regions,
    walker::Walker,
};
use crate::input;

const DEFAULT_EXCLUDED_FLAGS: u16 = 0x704;
#[derive(Clone)]
pub(super) enum Model {
    Simple(SimpleOptions),
    Bayesian {
        caller: Box<BayesianOptions>,
        record: RecordOptions,
    },
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Profile {
    Hifi,
    Hiseq,
    R10_4Sup,
    R10_4Dup,
    Ultima,
}

pub(crate) struct BayesianOverrides {
    pub(crate) adjust_quality: bool,
    pub(crate) use_mapping_quality: bool,
    pub(crate) adjust_mapping_quality: bool,
    pub(crate) mismatch_halo: usize,
    pub(crate) soft_clip_cost: u32,
    pub(crate) mapping_quality_scale: Option<f64>,
    pub(crate) minimum_mapping_quality: Option<u8>,
    pub(crate) maximum_mapping_quality: u8,
    pub(crate) default_quality: u8,
    pub(crate) heterozygous_probability: f64,
    pub(crate) indel_probability: f64,
    pub(crate) heterozygous_scale: Option<f64>,
    pub(crate) homopolymer_fix: Option<f64>,
}

#[derive(Clone)]
pub(crate) struct Options {
    model: Model,
    pub(crate) minimum_mapping_quality: u8,
    pub(crate) excluded_flags: u16,
    pub(crate) required_flags: u16,
    pub(crate) format: Format,
    pub(crate) show_deletions: bool,
    pub(crate) show_insertions: bool,
    pub(crate) mark_insertions: bool,
    pub(crate) all_positions: u8,
    pub(crate) reference: Option<PathBuf>,
    pub(crate) reference_quality: u8,
    pub(crate) region: Option<String>,
    pub(crate) regions_file: Option<PathBuf>,
    pub(crate) additional_threads: usize,
    pub(crate) line_width: usize,
}

impl Options {
    pub(crate) fn simple(call_fraction: f64) -> Self {
        Self {
            model: Model::Simple(SimpleOptions {
                use_quality: false,
                minimum_quality: 0,
                minimum_depth: 1,
                call_fraction,
                heterozygous_fraction: 0.5,
                ambiguous: false,
            }),
            minimum_mapping_quality: 0,
            excluded_flags: DEFAULT_EXCLUDED_FLAGS,
            required_flags: 0,
            format: Format::Pileup,
            show_deletions: false,
            show_insertions: true,
            mark_insertions: false,
            all_positions: 0,
            reference: None,
            reference_quality: 0,
            region: None,
            regions_file: None,
            additional_threads: 0,
            line_width: 70,
        }
    }

    pub(crate) fn bayesian(cutoff: i32, ambiguous: bool) -> Self {
        Self::with_bayesian(
            BayesianOptions::with_call_options(cutoff, ambiguous),
            RecordOptions::default(),
        )
    }

    #[cfg(test)]
    pub(super) fn bayesian_without_mapping_quality(cutoff: i32, ambiguous: bool) -> Self {
        Self::with_bayesian(
            BayesianOptions::without_mapping_quality(cutoff, ambiguous),
            RecordOptions {
                use_mapping_quality: false,
                ..RecordOptions::default()
            },
        )
    }

    fn with_bayesian(caller: BayesianOptions, record: RecordOptions) -> Self {
        Self {
            model: Model::Bayesian {
                caller: Box::new(caller),
                record,
            },
            minimum_mapping_quality: 0,
            excluded_flags: DEFAULT_EXCLUDED_FLAGS,
            required_flags: 0,
            format: Format::Pileup,
            show_deletions: false,
            show_insertions: true,
            mark_insertions: false,
            all_positions: 0,
            reference: None,
            reference_quality: 0,
            region: None,
            regions_file: None,
            additional_threads: 0,
            line_width: 70,
        }
    }

    pub(crate) fn apply_profile(&mut self, profile: Profile) -> Result<()> {
        let Model::Bayesian { caller, record } = &mut self.model else {
            return Err(RsomicsError::ConfigError(
                "machine profiles require Bayesian consensus".to_owned(),
            ));
        };
        caller.mode = BayesianMode::Recall;
        record.mode = BayesianMode::Recall;
        let calibration = match profile {
            Profile::Hifi => CalibrationPreset::Hifi,
            Profile::Hiseq => CalibrationPreset::Hiseq,
            Profile::R10_4Sup => CalibrationPreset::R10_4Sup,
            Profile::R10_4Dup => CalibrationPreset::R10_4Dup,
            Profile::Ultima => CalibrationPreset::Ultima,
        };
        caller.set_calibration_preset(calibration);
        if profile != Profile::Hiseq {
            record.homopolymer_fix = 0.3;
            caller.homopolymer_reduction = 0.01;
            caller.heterozygous_scale = 0.37;
            if profile == Profile::Ultima {
                caller.minimum_mapping_quality = 10;
                caller.mapping_quality_scale = 2.0;
            } else {
                caller.minimum_mapping_quality = 5;
                caller.mapping_quality_scale = 1.5;
            }
        }
        Ok(())
    }

    pub(crate) fn configure_simple(
        &mut self,
        use_quality: bool,
        minimum_base_quality: u8,
        minimum_depth: usize,
        heterozygous_fraction: f64,
        ambiguous: bool,
    ) -> Result<()> {
        let Model::Simple(options) = &mut self.model else {
            return Err(RsomicsError::ConfigError(
                "simple consensus options require simple mode".to_owned(),
            ));
        };
        options.use_quality = use_quality;
        options.minimum_quality = minimum_base_quality;
        options.minimum_depth = minimum_depth;
        options.heterozygous_fraction = heterozygous_fraction;
        options.ambiguous = ambiguous;
        Ok(())
    }

    pub(crate) fn configure_bayesian(
        &mut self,
        minimum_base_quality: u8,
        minimum_depth: usize,
        ambiguous: bool,
        compatibility_116: bool,
    ) -> Result<()> {
        let Model::Bayesian { caller, record } = &mut self.model else {
            return Err(RsomicsError::ConfigError(
                "Bayesian consensus options require Bayesian mode".to_owned(),
            ));
        };
        caller.minimum_base_quality = minimum_base_quality;
        caller.minimum_depth = minimum_depth;
        caller.ambiguous = ambiguous;
        if compatibility_116 {
            caller.mode = BayesianMode::Compatibility116;
            record.mode = BayesianMode::Compatibility116;
        }
        Ok(())
    }

    pub(crate) fn apply_bayesian_overrides(&mut self, overrides: BayesianOverrides) -> Result<()> {
        let Model::Bayesian { caller, record } = &mut self.model else {
            return Err(RsomicsError::ConfigError(
                "Bayesian consensus options require Bayesian mode".to_owned(),
            ));
        };
        let minimum_mapping_quality = overrides
            .minimum_mapping_quality
            .unwrap_or(caller.minimum_mapping_quality);
        if overrides.maximum_mapping_quality < minimum_mapping_quality {
            return Err(RsomicsError::ConfigError(format!(
                "Bayesian high-MQ {} is smaller than low-MQ {minimum_mapping_quality}",
                overrides.maximum_mapping_quality
            )));
        }
        for (name, value) in [
            (
                "heterozygous probability",
                overrides.heterozygous_probability,
            ),
            ("indel probability", overrides.indel_probability),
        ] {
            if !value.is_finite() || value <= 0.0 || value >= 1.0 {
                return Err(RsomicsError::ConfigError(format!(
                    "{name} must be a finite number between 0 and 1"
                )));
            }
        }
        for (name, value) in [
            ("mapping-quality scale", overrides.mapping_quality_scale),
            ("heterozygous scale", overrides.heterozygous_scale),
        ] {
            if value.is_some_and(|value| !value.is_finite() || value <= 0.0) {
                return Err(RsomicsError::ConfigError(format!(
                    "{name} must be a finite positive number"
                )));
            }
        }
        if overrides
            .homopolymer_fix
            .is_some_and(|value| !value.is_finite() || !(0.0..=1.0).contains(&value))
        {
            return Err(RsomicsError::ConfigError(
                "homopolymer adjustment must be a finite number from 0 to 1".to_owned(),
            ));
        }
        record.adjust_quality = overrides.adjust_quality;
        record.use_mapping_quality = overrides.use_mapping_quality;
        caller.use_mapping_quality = overrides.use_mapping_quality;
        caller.adjust_mapping_quality = overrides.adjust_mapping_quality;
        record.mismatch_halo = overrides.mismatch_halo;
        record.soft_clip_cost = overrides.soft_clip_cost;
        caller.maximum_mapping_quality = overrides.maximum_mapping_quality;
        caller.default_quality = overrides.default_quality;
        caller.heterozygous_probability = overrides.heterozygous_probability;
        caller.indel_probability = overrides.indel_probability;
        if let Some(value) = overrides.mapping_quality_scale {
            caller.mapping_quality_scale = value;
        }
        if let Some(value) = overrides.minimum_mapping_quality {
            caller.minimum_mapping_quality = value;
        }
        if let Some(value) = overrides.heterozygous_scale {
            caller.heterozygous_scale = value;
        }
        if let Some(value) = overrides.homopolymer_fix {
            record.homopolymer_fix = value;
        }
        Ok(())
    }

    pub(crate) fn set_calibration_preset(&mut self, preset: CalibrationPreset) -> Result<()> {
        let Model::Bayesian { caller, .. } = &mut self.model else {
            return Err(RsomicsError::ConfigError(
                "quality calibration requires Bayesian mode".to_owned(),
            ));
        };
        caller.set_calibration_preset(preset);
        Ok(())
    }

    pub(crate) fn set_calibration_file(&mut self, path: &Path) -> Result<()> {
        let Model::Bayesian { caller, .. } = &mut self.model else {
            return Err(RsomicsError::ConfigError(
                "quality calibration requires Bayesian mode".to_owned(),
            ));
        };
        caller.set_calibration_file(path)
    }
}

pub(crate) fn write_pileup(
    input_path: &Path,
    options: Options,
    output: impl Write,
) -> Result<Summary> {
    if options.region.is_some() && options.regions_file.is_some() {
        return Err(RsomicsError::ConfigError(
            "region and regions file are mutually exclusive".to_owned(),
        ));
    }
    let region = options
        .region
        .as_deref()
        .map(Region::from_str)
        .transpose()
        .map_err(|error| RsomicsError::ConfigError(format!("invalid region: {error}")))?;
    let mut reader = if region.is_some() || options.regions_file.is_some() {
        input::open_indexed(input_path, options.reference.as_deref())?
    } else {
        input::open(
            input_path,
            options.reference.as_deref(),
            options.additional_threads,
        )?
    };
    let header = reader.read_header(input_path)?;
    let references = header
        .reference_sequences()
        .iter()
        .map(|(name, reference)| {
            let name = name.to_vec().into_boxed_slice();
            let length = usize::from(reference.length()) as u64;
            Reference {
                label: name.clone(),
                name,
                length,
                start: 0,
                end: length as i64,
                enabled: true,
            }
        })
        .collect::<Vec<_>>();
    let reference_sequences = options
        .reference
        .as_deref()
        .map(|path| load_reference_sequences(path, &references))
        .transpose()?;
    let pileup_options = PileupOptions {
        filter: RecordFilter {
            flags: FlagFilter {
                skip_any_set: options.excluded_flags,
                skip_any_unset: options.required_flags,
                ..FlagFilter::default()
            },
            minimum_mapping_quality: options.minimum_mapping_quality,
            include_anomalous_pairs: true,
        },
        adjust_overlaps: false,
        maximum_depth_per_source: None,
    };
    let mut output = output;
    let mut summary = Summary::default();
    if let Some(path) = &options.regions_file {
        for interval in regions::read(path)? {
            let Some(region) = interval_region(&references, &interval)? else {
                continue;
            };
            summary.add(write_selection(
                &mut reader,
                &header,
                input_path,
                &references,
                reference_sequences.as_ref(),
                pileup_options,
                &options,
                Some(&region),
                &mut output,
            )?);
        }
    } else {
        summary.add(write_selection(
            &mut reader,
            &header,
            input_path,
            &references,
            reference_sequences.as_ref(),
            pileup_options,
            &options,
            region.as_ref(),
            &mut output,
        )?);
    }
    Ok(summary)
}

#[allow(clippy::too_many_arguments)]
fn write_selection(
    reader: &mut input::Reader,
    header: &noodles::sam::Header,
    input_path: &Path,
    base_references: &[Reference],
    reference_sequences: Option<&Vec<Arc<Sequence>>>,
    pileup_options: PileupOptions,
    options: &Options,
    region: Option<&Region>,
    output: &mut impl Write,
) -> Result<Summary> {
    let mut references = base_references.to_vec();
    if let Some(region) = region {
        select_region(&mut references, region)?;
    }
    let mut pileup = PileupEngine::with_record_state(
        references.iter().map(|reference| reference.length),
        pileup_options,
    );
    let (caller, record_options) = match options.model.clone() {
        Model::Simple(options) => (
            Caller::Simple(options),
            RecordOptions {
                use_mapping_quality: false,
                ..RecordOptions::default()
            },
        ),
        Model::Bayesian { caller, record } => (
            Caller::Bayesian(Box::new(BayesianCaller::new(*caller))),
            record,
        ),
    };
    let mut walker = Walker::new(caller);
    let mut output = Output::new(
        output,
        OutputConfiguration {
            format: options.format,
            show_deletions: options.show_deletions,
            show_insertions: options.show_insertions,
            mark_insertions: options.mark_insertions,
            all_positions: options.all_positions,
            reference_sequences: reference_sequences.cloned(),
            reference_quality: options.reference_quality,
            line_width: options.line_width,
        },
    );
    if region.is_some() {
        output.begin_selected_reference(&references);
    }

    let mut visit = |record| {
        let state = RecordState::new(&record, record_options)?;
        pileup
            .push_with_state(record, state)
            .map_err(|error| pileup_error(input_path, error))?;
        drain(&mut pileup, &mut walker, &references, &mut output)
    };
    if let Some(region) = &region {
        reader.visit_owned_raw_region(header, input_path, region, &mut visit)?;
    } else {
        reader.visit_owned_raw_records(header, input_path, &mut visit)?;
    }
    pileup
        .finish()
        .map_err(|error| pileup_error(input_path, error))?;
    drain(&mut pileup, &mut walker, &references, &mut output)?;
    output.finish(&references)
}

fn select_region(references: &mut [Reference], region: &Region) -> Result<()> {
    let reference_id = references
        .iter()
        .position(|reference| reference.name.as_ref() == region.name())
        .ok_or_else(|| {
            RsomicsError::ConfigError(format!(
                "region reference is absent from the alignment header: {}",
                String::from_utf8_lossy(region.name())
            ))
        })?;
    for reference in references.iter_mut() {
        reference.enabled = false;
    }
    let reference = &mut references[reference_id];
    let interval = region.interval();
    let start = interval
        .start()
        .map_or(0, |position| usize::from(position) - 1);
    let end = interval
        .end()
        .map_or(reference.length as usize, usize::from);
    if end > reference.length as usize {
        return Err(RsomicsError::ConfigError(format!(
            "region end {end} exceeds reference {} length {}",
            String::from_utf8_lossy(&reference.name),
            reference.length
        )));
    }
    reference.start = start as i64;
    reference.end = end as i64;
    reference.label = if start == 0 && end == reference.length as usize {
        reference.name.clone()
    } else {
        region.to_string().into_bytes().into_boxed_slice()
    };
    reference.enabled = true;
    Ok(())
}

fn interval_region(
    references: &[Reference],
    interval: &regions::Interval,
) -> Result<Option<Region>> {
    let reference = references
        .iter()
        .find(|reference| reference.name == interval.name)
        .ok_or_else(|| {
            RsomicsError::ConfigError(format!(
                "region-file reference is absent from the alignment header: {}",
                String::from_utf8_lossy(&interval.name)
            ))
        })?;
    if interval.start > reference.length {
        return Err(RsomicsError::ConfigError(format!(
            "region-file start {} exceeds reference {} length {}",
            interval.start,
            String::from_utf8_lossy(&interval.name),
            reference.length
        )));
    }
    let end = interval.end.min(reference.length);
    if interval.start == end {
        return Ok(None);
    }
    let start =
        noodles::core::Position::try_from(interval.start as usize + 1).map_err(|error| {
            RsomicsError::ConfigError(format!("region-file start is out of range: {error}"))
        })?;
    let end = noodles::core::Position::try_from(end as usize).map_err(|error| {
        RsomicsError::ConfigError(format!("region-file end is out of range: {error}"))
    })?;
    Ok(Some(Region::new(interval.name.to_vec(), start..=end)))
}

fn load_reference_sequences(path: &Path, references: &[Reference]) -> Result<Vec<Arc<Sequence>>> {
    let repository = input::reference_repository(path)?;
    references
        .iter()
        .map(|reference| {
            let sequence = repository
                .get(&reference.name)
                .ok_or_else(|| {
                    RsomicsError::ConfigError(format!(
                        "reference {} is absent from {}",
                        String::from_utf8_lossy(&reference.name),
                        path.display()
                    ))
                })?
                .map_err(|error| {
                    RsomicsError::ConfigError(format!(
                        "reading reference {} from {}: {error}",
                        String::from_utf8_lossy(&reference.name),
                        path.display()
                    ))
                })?;
            if sequence.len() as u64 != reference.length {
                return Err(RsomicsError::ConfigError(format!(
                    "reference {} has length {} in {} but length {} in the alignment header",
                    String::from_utf8_lossy(&reference.name),
                    sequence.len(),
                    path.display(),
                    reference.length
                )));
            }
            Ok(sequence)
        })
        .collect()
}

fn drain<W: Write>(
    pileup: &mut PileupEngine<RecordState>,
    walker: &mut Walker,
    references: &[Reference],
    output: &mut Output<W>,
) -> Result<bool> {
    pileup.drain(|column| {
        walker.visit(column, |call, observations| {
            output.write(references, call, observations)
        })
    })?;
    Ok(true)
}

fn pileup_error(input: &Path, error: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!("processing {}: {error}", input.display()))
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use std::process::Command;

    use super::*;

    fn root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/upstream/samtools-consensus")
    }

    fn assert_input(input: &str, expected: &str, options: Options) {
        let root = root();
        let mut output = Vec::new();

        write_pileup(&root.join(input), options, &mut output).unwrap();

        assert_eq!(
            output,
            std::fs::read(root.join("expected").join(expected)).unwrap()
        );
    }

    fn assert_output(expected: &str, options: Options) {
        assert_input("consen1.sam", expected, options);
    }

    fn simple(options: &mut Options) -> &mut SimpleOptions {
        let Model::Simple(options) = &mut options.model else {
            unreachable!()
        };
        options
    }

    fn fastq(mut options: Options) -> Options {
        options.format = Format::Fastq;
        options
    }

    #[test]
    fn simple_pileup_matches_samtools_1_24() {
        assert_output("1p.out", Options::simple(0.6));
    }

    #[test]
    fn simple_pileup_option_matrix_matches_samtools_1_24() {
        let options = Options::simple(0.601);
        assert_output("5p.out", options);

        let mut options = Options::simple(0.6);
        simple(&mut options).ambiguous = true;
        simple(&mut options).heterozygous_fraction = 0.66;
        assert_output("6p.out", options);

        let mut options = Options::simple(0.6);
        simple(&mut options).ambiguous = true;
        simple(&mut options).heterozygous_fraction = 0.33;
        assert_output("7p.out", options);

        let mut options = Options::simple(0.6);
        simple(&mut options).ambiguous = true;
        simple(&mut options).heterozygous_fraction = 0.25;
        assert_output("8p.out", options);

        let mut options = Options::simple(0.666);
        simple(&mut options).use_quality = true;
        assert_output("9p.out", options);

        let mut options = Options::simple(0.667);
        simple(&mut options).use_quality = true;
        assert_output("9.1p.out", options);

        let mut options = Options::simple(0.6);
        options.show_deletions = true;
        simple(&mut options).minimum_depth = 5;
        assert_output("11.3p.out", options);

        let mut options = Options::simple(0.6);
        options.show_deletions = true;
        simple(&mut options).minimum_depth = 6;
        assert_output("11.6p.out", options);

        let mut options = Options::simple(0.6);
        options.minimum_mapping_quality = 40;
        simple(&mut options).minimum_depth = 5;
        assert_output("22p.out", options);

        let mut options = Options::simple(0.6);
        options.minimum_mapping_quality = 41;
        simple(&mut options).minimum_depth = 4;
        assert_output("23p.out", options);

        let mut options = Options::simple(0.6);
        options.minimum_mapping_quality = 41;
        simple(&mut options).minimum_depth = 5;
        assert_output("24p.out", options);

        let mut options = Options::simple(0.6);
        simple(&mut options).minimum_quality = 5;
        simple(&mut options).minimum_depth = 2;
        assert_output("25p.out", options);

        let mut options = Options::simple(0.6);
        simple(&mut options).minimum_quality = 6;
        simple(&mut options).minimum_depth = 2;
        assert_output("26p.out", options);
    }

    #[test]
    fn bayesian_pileup_option_matrix_matches_samtools_1_24() {
        assert_output(
            "18p.out",
            Options::bayesian_without_mapping_quality(0, false),
        );
        assert_output(
            "19p.out",
            Options::bayesian_without_mapping_quality(19, false),
        );
        assert_output(
            "20p.out",
            Options::bayesian_without_mapping_quality(30, true),
        );
        assert_output(
            "21p.out",
            Options::bayesian_without_mapping_quality(31, true),
        );
    }

    #[test]
    fn fasta_and_fastq_match_samtools_1_24() {
        let mut fasta = Options::simple(0.6);
        fasta.format = Format::Fasta;
        assert_output("1.out", fasta);

        let mut fastq = Options::simple(0.6);
        fastq.format = Format::Fastq;
        assert_output("1q.out", fastq);

        let mut bayesian = Options::bayesian_without_mapping_quality(0, false);
        bayesian.format = Format::Fastq;
        assert_output("18q.out", bayesian);
    }

    #[test]
    fn fastq_option_matrix_matches_samtools_1_24() {
        let mut options = fastq(Options::simple(0.6));
        options.show_deletions = true;
        assert_output("2q.out", options);

        let mut options = fastq(Options::simple(0.6));
        options.show_insertions = false;
        assert_output("3q.out", options);

        let mut options = fastq(Options::simple(0.6));
        options.show_deletions = true;
        options.show_insertions = false;
        assert_output("4q.out", options);

        assert_output("5q.out", fastq(Options::simple(0.601)));

        for (expected, fraction) in [("6q.out", 0.66), ("7q.out", 0.33), ("8q.out", 0.25)] {
            let mut options = fastq(Options::simple(0.6));
            simple(&mut options).ambiguous = true;
            simple(&mut options).heterozygous_fraction = fraction;
            assert_output(expected, options);
        }

        for (expected, fraction) in [("9q.out", 0.666), ("9.1q.out", 0.667)] {
            let mut options = fastq(Options::simple(fraction));
            simple(&mut options).use_quality = true;
            assert_output(expected, options);
        }

        for (expected, fraction) in [("10q.out", 0.375), ("10.1q.out", 0.376)] {
            let mut options = fastq(Options::simple(0.75));
            simple(&mut options).use_quality = true;
            simple(&mut options).ambiguous = true;
            simple(&mut options).heterozygous_fraction = fraction;
            assert_output(expected, options);
        }

        let mut options = fastq(Options::simple(0.6));
        options.show_deletions = true;
        simple(&mut options).minimum_depth = 5;
        assert_output("11.3q.out", options);

        let mut options = fastq(Options::simple(0.6));
        options.show_deletions = true;
        simple(&mut options).minimum_depth = 6;
        assert_output("11.6q.out", options);
    }

    #[test]
    fn reference_extent_and_line_wrapping_match_samtools_1_24() {
        assert_input("consen2.sam", "12q.out", fastq(Options::simple(0.75)));
        assert_input("consen2.sam", "12p.out", Options::simple(0.75));

        let mut fastq_all = fastq(Options::simple(0.75));
        fastq_all.all_positions = 1;
        assert_input("consen2.sam", "13q.out", fastq_all);

        let mut pileup_all = Options::simple(0.75);
        pileup_all.all_positions = 1;
        assert_input("consen2.sam", "13p.out", pileup_all);

        let mut wrapped = fastq(Options::simple(0.75));
        wrapped.all_positions = 1;
        wrapped.line_width = 7;
        assert_input("consen2.sam", "17q.out", wrapped);
    }

    #[test]
    fn all_reference_modes_match_samtools_1_24() {
        for (expected, all_positions) in [("30.out", 0), ("31.out", 1), ("32.out", 2)] {
            let mut options = fastq(Options::bayesian(0, false));
            options.show_deletions = true;
            options.show_insertions = false;
            options.all_positions = all_positions;
            assert_input("consen1c.sam", expected, options);
        }
    }

    #[test]
    fn reference_backfill_matches_samtools_1_24() {
        for (expected, all_positions) in [("30T.out", 0), ("31T.out", 1), ("32T.out", 2)] {
            let mut options = fastq(Options::bayesian(0, false));
            options.show_deletions = true;
            options.show_insertions = false;
            options.all_positions = all_positions;
            options.reference = Some(root().join("consen1c.fa"));
            options.reference_quality = 20;
            assert_input("consen1c.sam", expected, options);
        }
    }

    #[test]
    fn pileup_reference_extent_matches_samtools_1_24() {
        for (expected, all_positions) in [("40.out", 0), ("41.out", 1), ("42.out", 2)] {
            let mut options = Options::bayesian(0, false);
            options.show_deletions = true;
            options.show_insertions = false;
            options.all_positions = all_positions;
            assert_input("consen1c.sam", expected, options);
        }

        for (expected, all_positions) in [("40T.out", 0), ("41T.out", 1), ("42T.out", 2)] {
            let mut options = Options::bayesian(0, false);
            options.show_deletions = true;
            options.show_insertions = false;
            options.all_positions = all_positions;
            options.reference = Some(root().join("consen1c.fa"));
            options.reference_quality = 20;
            assert_input("consen1c.sam", expected, options);
        }
    }

    #[test]
    fn insertion_markers_match_samtools_1_24() {
        let mut options = fastq(Options::simple(0.6));
        options.mark_insertions = true;
        let mut output = Vec::new();

        write_pileup(&root().join("consen1.sam"), options, &mut output).unwrap();

        assert_eq!(
            output,
            b"@c2\nCC_T_T_TAAGGAA_T_T_TCC\n+\n~~_]_]_]~~]q~~_]_]_]~~\n"
        );
    }

    #[test]
    fn machine_profiles_match_samtools_1_24_parameter_sets() {
        for (profile, minimum_mapping_quality, mapping_quality_scale) in [
            (Profile::Hifi, 5, 1.5),
            (Profile::R10_4Sup, 5, 1.5),
            (Profile::R10_4Dup, 5, 1.5),
            (Profile::Ultima, 10, 2.0),
        ] {
            let mut options = Options::bayesian(10, false);
            options.apply_profile(profile).unwrap();
            let Model::Bayesian { caller, record } = options.model else {
                unreachable!()
            };
            assert_eq!(caller.minimum_mapping_quality, minimum_mapping_quality);
            assert_eq!(caller.mapping_quality_scale, mapping_quality_scale);
            assert_eq!(caller.heterozygous_scale, 0.37);
            assert_eq!(caller.homopolymer_reduction, 0.01);
            assert_eq!(record.homopolymer_fix, 0.3);
            assert_eq!(caller.mode, BayesianMode::Recall);
            assert_eq!(record.mode, BayesianMode::Recall);
        }

        let mut hiseq = Options::bayesian(10, false);
        hiseq.apply_profile(Profile::Hiseq).unwrap();
        let Model::Bayesian { caller, record } = hiseq.model else {
            unreachable!()
        };
        assert_eq!(caller.minimum_mapping_quality, 1);
        assert_eq!(caller.mapping_quality_scale, 1.0);
        assert_eq!(caller.heterozygous_scale, 1.0);
        assert_eq!(record.homopolymer_fix, 0.0);

        assert!(Options::simple(0.6).apply_profile(Profile::Hifi).is_err());
    }

    #[test]
    fn hiseq_insertion_context_matches_samtools_1_24() {
        let mut options = Options::bayesian(10, false);
        options.apply_profile(Profile::Hiseq).unwrap();
        options.show_deletions = true;
        let mut output = Vec::new();

        write_pileup(&root().join("consen1c.sam"), options, &mut output).unwrap();

        let output = String::from_utf8(output).unwrap();
        let line = output
            .lines()
            .find(|line| line.starts_with("c2\t4\t3\t"))
            .unwrap();
        assert_eq!(line, "c2\t4\t3\t5\t*\t2\t**TTT\tIIIII");
    }

    #[test]
    fn region_selection_sets_checked_output_bounds() {
        let mut references = vec![Reference {
            name: Box::from(&b"c1"[..]),
            label: Box::from(&b"c1"[..]),
            length: 20,
            start: 0,
            end: 20,
            enabled: true,
        }];

        select_region(&mut references, &Region::from_str("c1:3-8").unwrap()).unwrap();

        assert_eq!(references[0].start, 2);
        assert_eq!(references[0].end, 8);
        assert_eq!(references[0].label.as_ref(), b"c1:3-8");
        assert!(select_region(&mut references, &Region::from_str("c2:1-2").unwrap()).is_err());
        assert!(select_region(&mut references, &Region::from_str("c1:1-21").unwrap()).is_err());
    }

    #[test]
    #[ignore = "release oracle: requires samtools 1.24"]
    fn indexed_region_matches_samtools_1_24() {
        let version = Command::new("samtools").arg("--version").output().unwrap();
        assert!(String::from_utf8_lossy(&version.stdout).starts_with("samtools 1.24"));
        let directory = tempfile::tempdir().unwrap();
        let bam = directory.path().join("consen2.bam");
        let status = Command::new("samtools")
            .args(["view", "--write-index"])
            .arg(root().join("consen2.sam"))
            .arg("-o")
            .arg(&bam)
            .status()
            .unwrap();
        assert!(status.success());

        let mut pileup = Options::simple(0.75);
        pileup.region = Some("c2:2-13".to_owned());
        pileup.all_positions = 1;
        let mut output = Vec::new();
        write_pileup(&bam, pileup, &mut output).unwrap();
        assert_eq!(
            output,
            std::fs::read(root().join("expected/16p.out")).unwrap()
        );

        let mut fastq = fastq(Options::simple(0.75));
        fastq.region = Some("c2:2-13".to_owned());
        fastq.all_positions = 1;
        output.clear();
        write_pileup(&bam, fastq, &mut output).unwrap();
        assert_eq!(
            output,
            std::fs::read(root().join("expected/16q.out")).unwrap()
        );
    }

    #[test]
    #[ignore = "release oracle: requires samtools 1.24"]
    fn indexed_bed_regions_match_samtools_1_24() {
        let version = Command::new("samtools").arg("--version").output().unwrap();
        assert!(String::from_utf8_lossy(&version.stdout).starts_with("samtools 1.24"));
        let directory = tempfile::tempdir().unwrap();
        let bam = directory.path().join("consen4.bam");
        let status = Command::new("samtools")
            .args(["view", "--write-index"])
            .arg(root().join("consen4.sam"))
            .arg("-o")
            .arg(&bam)
            .status()
            .unwrap();
        assert!(status.success());

        for (expected, format, all_positions) in [
            ("bf1.out", Format::Fasta, 0),
            ("bf2.out", Format::Fasta, 1),
            ("bp1.out", Format::Pileup, 0),
            ("bp2.out", Format::Pileup, 1),
        ] {
            let mut options = Options::bayesian(0, false);
            options.format = format;
            options.all_positions = all_positions;
            options.regions_file = Some(root().join("consen4.bed"));
            let mut output = Vec::new();

            write_pileup(&bam, options, &mut output).unwrap();

            assert_eq!(
                output,
                std::fs::read(root().join("expected").join(expected)).unwrap(),
                "{expected}"
            );
        }
    }

    #[test]
    #[ignore = "release oracle: requires samtools 1.24"]
    fn cram_with_reference_matches_samtools_1_24() {
        let version = Command::new("samtools").arg("--version").output().unwrap();
        assert!(String::from_utf8_lossy(&version.stdout).starts_with("samtools 1.24"));
        let directory = tempfile::tempdir().unwrap();
        let cram = directory.path().join("consen1c.cram");
        let reference = root().join("consen1c.fa");
        let status = Command::new("samtools")
            .args(["view", "-C", "-T"])
            .arg(&reference)
            .arg("-o")
            .arg(&cram)
            .arg(root().join("consen1c.sam"))
            .status()
            .unwrap();
        assert!(status.success());

        let mut options = fastq(Options::bayesian(0, false));
        options.show_deletions = true;
        options.show_insertions = false;
        options.reference = Some(reference);
        options.reference_quality = 20;
        let mut output = Vec::new();

        write_pileup(&cram, options, &mut output).unwrap();

        assert_eq!(
            output,
            std::fs::read(root().join("expected/30T.out")).unwrap()
        );
    }
}
