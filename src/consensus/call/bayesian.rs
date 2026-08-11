use std::sync::LazyLock;

use super::Call;

mod calibration;
pub(crate) use calibration::CalibrationPreset;

const HYPOTHESIS_TO_BASE: [usize; 15] = [0, 5, 5, 5, 5, 1, 5, 5, 5, 2, 5, 5, 3, 5, 4];
const HYPOTHESIS_TO_PAIR: [usize; 15] = [0, 1, 2, 3, 4, 6, 7, 8, 9, 12, 13, 14, 18, 19, 24];
const PURE_HYPOTHESES: [usize; 5] = [0, 5, 9, 12, 14];
const TEN_LOG2_OVER_LOG10: f64 = 3.0103;
const CACHED_DEPTHS: usize = 31;
const CACHED_MAPPING_QUALITIES: usize = 256;
const CACHED_BASE_QUALITIES: usize = 101;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(in crate::consensus) enum BayesianMode {
    Recall,
    Compatibility116,
}

#[derive(Clone)]
pub(in crate::consensus) struct QualityCalibration {
    substitution: [u8; 101],
    undercall: [u8; 101],
    overcall: [u8; 101],
}

impl Default for QualityCalibration {
    fn default() -> Self {
        let identity = std::array::from_fn(|quality| if quality < 100 { quality as u8 } else { 0 });
        Self {
            substitution: identity,
            undercall: identity,
            overcall: identity,
        }
    }
}

impl QualityCalibration {
    pub(in crate::consensus) fn preset(preset: CalibrationPreset) -> Self {
        calibration::preset(preset)
    }

    pub(in crate::consensus) fn from_path(path: &std::path::Path) -> rsomics_common::Result<Self> {
        calibration::read(path)
    }
}

#[derive(Clone)]
pub(in crate::consensus) struct BayesianOptions {
    pub(in crate::consensus) use_mapping_quality: bool,
    pub(in crate::consensus) adjust_mapping_quality: bool,
    pub(in crate::consensus) mapping_quality_scale: f64,
    pub(in crate::consensus) minimum_mapping_quality: u8,
    pub(in crate::consensus) maximum_mapping_quality: u8,
    pub(in crate::consensus) minimum_base_quality: u8,
    pub(in crate::consensus) minimum_depth: usize,
    pub(in crate::consensus) default_quality: u8,
    pub(super) cutoff: i32,
    pub(in crate::consensus) ambiguous: bool,
    pub(in crate::consensus) heterozygous_probability: f64,
    pub(in crate::consensus) indel_probability: f64,
    pub(in crate::consensus) heterozygous_scale: f64,
    pub(in crate::consensus) homopolymer_reduction: f64,
    pub(super) calibration: QualityCalibration,
}

impl Default for BayesianOptions {
    fn default() -> Self {
        Self {
            use_mapping_quality: true,
            adjust_mapping_quality: true,
            mapping_quality_scale: 1.0,
            minimum_mapping_quality: 1,
            maximum_mapping_quality: 60,
            minimum_base_quality: 0,
            minimum_depth: 1,
            default_quality: 10,
            cutoff: 10,
            ambiguous: false,
            heterozygous_probability: 1e-3,
            indel_probability: 2e-4,
            heterozygous_scale: 1.0,
            homopolymer_reduction: 0.01,
            calibration: QualityCalibration::default(),
        }
    }
}

impl BayesianOptions {
    pub(in crate::consensus) fn with_call_options(cutoff: i32, ambiguous: bool) -> Self {
        Self {
            cutoff,
            ambiguous,
            ..Self::default()
        }
    }

    #[cfg(test)]
    pub(in crate::consensus) fn without_mapping_quality(cutoff: i32, ambiguous: bool) -> Self {
        Self {
            use_mapping_quality: false,
            ..Self::with_call_options(cutoff, ambiguous)
        }
    }

    pub(in crate::consensus) fn set_calibration_preset(&mut self, preset: CalibrationPreset) {
        self.calibration = QualityCalibration::preset(preset);
    }

    pub(in crate::consensus) fn set_calibration_file(
        &mut self,
        path: &std::path::Path,
    ) -> rsomics_common::Result<()> {
        self.calibration = QualityCalibration::from_path(path)?;
        Ok(())
    }
}

#[derive(Clone, Copy)]
pub(in crate::consensus) struct BayesianObservation {
    pub(in crate::consensus) base: u8,
    pub(in crate::consensus) quality: u8,
    pub(in crate::consensus) mapping_quality: u8,
    pub(in crate::consensus) local_mismatch_tenths: u32,
    pub(in crate::consensus) homopolymer: u8,
    pub(in crate::consensus) reference_skip: bool,
    pub(in crate::consensus) missing_quality: bool,
}

impl BayesianObservation {
    pub(in crate::consensus) fn new(base: u8, quality: u8) -> Self {
        Self {
            base,
            quality,
            mapping_quality: 0,
            local_mismatch_tenths: 0,
            homopolymer: 0,
            reference_skip: false,
            missing_quality: false,
        }
    }

    pub(in crate::consensus) fn with_record_context(
        mut self,
        mapping_quality: u8,
        local_mismatch_tenths: u32,
        homopolymer: u8,
        reference_skip: bool,
        missing_quality: bool,
    ) -> Self {
        self.mapping_quality = mapping_quality;
        self.local_mismatch_tenths = local_mismatch_tenths;
        self.homopolymer = homopolymer;
        self.reference_skip = reference_skip;
        self.missing_quality = missing_quality;
        self
    }
}

pub(in crate::consensus) struct BayesianCaller {
    options: BayesianOptions,
    probabilities: Probabilities,
    zero_local_quality: Option<Box<[u8]>>,
}

impl BayesianCaller {
    pub(in crate::consensus) fn new(options: BayesianOptions) -> Self {
        let probabilities = Probabilities::new(&options);
        let zero_local_quality = build_zero_local_quality_cache(&options);
        Self {
            options,
            probabilities,
            zero_local_quality,
        }
    }

    pub(in crate::consensus) fn call(&self, observations: &[BayesianObservation]) -> Call {
        let consensus = self.calculate(observations);
        let (mut base, mut quality) =
            if consensus.depth < self.options.minimum_depth && consensus.call != 4 {
                (b'N', 0)
            } else if consensus.heterozygous_log_odds > 0 && self.options.ambiguous {
                (
                    b"AMRWaMCSYcRSGKgWYKTtacgt*"[consensus.heterozygous_call],
                    consensus.heterozygous_log_odds,
                )
            } else {
                (b"ACGT*"[consensus.call], consensus.phred)
            };
        if quality < self.options.cutoff
            && base != b'*'
            && consensus.heterozygous_call % 5 != 4
            && consensus.heterozygous_call / 5 != 4
        {
            base = b'N';
            quality = 0;
        }
        Call { base, quality }
    }

    fn calculate(&self, observations: &[BayesianObservation]) -> Consensus {
        let mut scores = [0.0; 15];
        let mut counts = [0usize; 6];
        let original_depth = observations.len();
        let mut depth = 0usize;

        for observation in observations {
            if observation.quality < self.options.minimum_base_quality || observation.reference_skip
            {
                continue;
            }
            let mut quality = if observation.missing_quality {
                self.options.default_quality
            } else {
                observation.quality
            };
            if self.options.use_mapping_quality {
                quality = self.cached_adjusted_quality(quality, observation, original_depth);
            }
            quality = quality.clamp(1, 100);
            let homopolymer_quality = ((f64::from(quality)
                - (f64::from(observation.homopolymer) - 2.0)
                    * self.probabilities.homopolymer_multiplier)
                .max(1.0) as usize)
                .min(100);
            let quality = usize::from(quality);
            let base = map_base(observation.base);
            let mismatch = self.probabilities.mismatch[quality];
            let both_match = self.probabilities.both_match[quality] - mismatch;
            let one_match = self.probabilities.one_match[quality] - mismatch;
            let both_overcall = self.probabilities.both_overcall[homopolymer_quality] - mismatch;
            let one_overcall = self.probabilities.one_overcall[homopolymer_quality] - mismatch;
            let overcall_mismatch =
                self.probabilities.overcall_mismatch[homopolymer_quality] - mismatch;
            let both_undercall = self.probabilities.both_undercall[homopolymer_quality] - mismatch;
            let one_undercall = self.probabilities.one_undercall[homopolymer_quality] - mismatch;
            let both_gap_match = self.probabilities.both_gap_match[homopolymer_quality] - mismatch;

            counts[base] += 1;
            accumulate(
                &mut scores,
                base,
                both_match,
                one_match,
                both_overcall,
                one_overcall,
                overcall_mismatch,
                both_undercall,
                one_undercall,
                both_gap_match,
            );
            depth += 1;
        }

        if depth == 0 || depth == counts[5] {
            return Consensus::default();
        }

        let mut shift = f64::NEG_INFINITY;
        let mut pure_score = f64::NEG_INFINITY;
        let mut heterozygous_score = f64::NEG_INFINITY;
        let mut pure = 0usize;
        let mut heterozygous = 0usize;
        for (index, score) in scores.iter_mut().enumerate() {
            *score += self.probabilities.log_prior[index];
            shift = shift.max(*score);
            if PURE_HYPOTHESES.contains(&index) {
                if *score > pure_score {
                    pure_score = *score;
                    pure = index;
                }
            } else if *score > heterozygous_score {
                heterozygous_score = *score;
                heterozygous = index;
            }
        }

        let minimum_exponent = -1021.0 * std::f64::consts::LN_2 + 1.0;
        let mut normalized = [0.0; 15];
        let mut alternatives = [0.0; 15];
        for (score, normalized) in scores.iter_mut().zip(&mut normalized) {
            *score -= shift;
            *normalized = if *score > minimum_exponent {
                fast_exp(*score)
            } else {
                f64::MIN_POSITIVE
            };
        }

        let mut left = 0.0;
        let mut right = 0.0;
        for index in 0..15 {
            alternatives[index] += left;
            alternatives[14 - index] += right;
            left += normalized[index];
            right += normalized[14 - index];
        }

        let pure_alternatives = alternatives[pure].max(f64::MIN_POSITIVE);
        let pure_probability = normalized[pure];
        let phred = if pure_probability == 1.0 && pure_alternatives < 0.01 {
            (phred_log(pure_alternatives) + 0.5) as i32
        } else {
            (phred_log(1.0 - pure_probability / (pure_alternatives + pure_probability)) + 0.5)
                as i32
        }
        .max(0);
        let heterozygous_alternatives = alternatives[heterozygous].max(f64::MIN_POSITIVE);
        let heterozygous_log_odds = (TEN_LOG2_OVER_LOG10
            * (fast_log2(normalized[heterozygous]) - fast_log2(heterozygous_alternatives))
            + 0.5) as i32;

        Consensus {
            call: HYPOTHESIS_TO_BASE[pure],
            heterozygous_call: HYPOTHESIS_TO_PAIR[heterozygous],
            heterozygous_log_odds,
            phred,
            depth,
        }
    }

    fn cached_adjusted_quality(
        &self,
        quality: u8,
        observation: &BayesianObservation,
        depth: usize,
    ) -> u8 {
        if quality <= 100
            && observation.local_mismatch_tenths == 0
            && let Some(cache) = &self.zero_local_quality
        {
            return cache[quality_cache_index(depth, observation.mapping_quality, quality)];
        }
        adjusted_quality(quality, observation, depth, &self.options)
    }
}

#[derive(Clone, Copy, Default)]
struct Consensus {
    call: usize,
    heterozygous_call: usize,
    heterozygous_log_odds: i32,
    phred: i32,
    depth: usize,
}

struct Probabilities {
    log_prior: [f64; 15],
    both_match: [f64; 101],
    mismatch: [f64; 101],
    one_match: [f64; 101],
    overcall_mismatch: [f64; 101],
    one_overcall: [f64; 101],
    both_overcall: [f64; 101],
    both_undercall: [f64; 101],
    one_undercall: [f64; 101],
    both_gap_match: [f64; 101],
    homopolymer_multiplier: f64,
}

impl Probabilities {
    fn new(options: &BayesianOptions) -> Self {
        let mut prior = [options.heterozygous_probability / 6.0; 25];
        for index in [0, 6, 12, 18, 24] {
            prior[index] = 1.0;
        }
        for index in (4..24).step_by(5) {
            prior[index] = options.indel_probability / 6.0;
        }
        for value in &mut prior[20..24] {
            *value = options.indel_probability / 6.0;
        }
        let log_prior = HYPOTHESIS_TO_PAIR.map(|index| prior[index].ln());
        let mut probabilities = Self {
            log_prior,
            both_match: [0.0; 101],
            mismatch: [0.0; 101],
            one_match: [0.0; 101],
            overcall_mismatch: [0.0; 101],
            one_overcall: [0.0; 101],
            both_overcall: [0.0; 101],
            both_undercall: [0.0; 101],
            one_undercall: [0.0; 101],
            both_gap_match: [0.0; 101],
            homopolymer_multiplier: options.homopolymer_reduction,
        };

        for quality in 1..=100 {
            let substitution = calibrated_probability(options.calibration.substitution[quality]);
            probabilities.both_match[quality] = substitution.ln();
            probabilities.mismatch[quality] = ((1.0 - substitution) / 3.0).ln();
            probabilities.one_match[quality] = ((probabilities.both_match[quality].exp()
                + probabilities.mismatch[quality].exp())
                / 2.0)
                .ln()
                + options.heterozygous_scale.ln();

            let overcall = calibrated_probability(options.calibration.overcall[quality]);
            probabilities.both_overcall[quality] = ((1.0 - overcall) / 3.0).ln();
            if probabilities.both_overcall[quality] > probabilities.both_match[quality] - 0.5 {
                probabilities.both_overcall[quality] = probabilities.both_match[quality] - 0.5;
            }
            probabilities.overcall_mismatch[quality] = ((probabilities.both_overcall[quality]
                .exp()
                + probabilities.mismatch[quality].exp())
                / 2.0)
                .ln();
            probabilities.one_overcall[quality] = ((probabilities.both_overcall[quality].exp()
                + probabilities.both_match[quality].exp())
                / 2.0)
                .ln();
            if probabilities.one_overcall[quality] > probabilities.one_match[quality] + 0.5 {
                probabilities.one_overcall[quality] = probabilities.one_match[quality] + 0.5;
            }

            let undercall = calibrated_probability(options.calibration.undercall[quality]);
            probabilities.both_gap_match[quality] = undercall.ln();
            probabilities.both_undercall[quality] = ((1.0 - undercall) / 3.0).ln();
            if probabilities.both_undercall[quality] > probabilities.both_match[quality] - 0.5 {
                probabilities.both_undercall[quality] = probabilities.both_match[quality] - 0.5;
            }
            probabilities.one_undercall[quality] = ((probabilities.both_undercall[quality].exp()
                + probabilities.both_gap_match[quality].exp())
                / 2.0)
                .ln();
        }

        for table in [
            &mut probabilities.both_match,
            &mut probabilities.mismatch,
            &mut probabilities.one_match,
            &mut probabilities.overcall_mismatch,
            &mut probabilities.one_overcall,
            &mut probabilities.both_overcall,
            &mut probabilities.both_undercall,
            &mut probabilities.one_undercall,
            &mut probabilities.both_gap_match,
        ] {
            table[0] = table[1];
        }
        probabilities
    }
}

fn calibrated_probability(quality: u8) -> f64 {
    1.0 - 10.0f64.powf(-f64::from(quality) / 10.0)
}

fn adjusted_quality(
    quality: u8,
    observation: &BayesianObservation,
    depth: usize,
    options: &BayesianOptions,
) -> u8 {
    static BASE_ERROR: LazyLock<[f64; 256]> =
        LazyLock::new(|| std::array::from_fn(|quality| 10.0f64.powf(-(quality as f64) / 10.0)));
    static MAPPING_ERROR: LazyLock<[f64; 256]> = LazyLock::new(|| {
        std::array::from_fn(|quality| 10.0f64.powf(-((quality as f64) * 0.9) / 10.0))
    });
    let mut mapping_quality = f64::from(observation.mapping_quality);
    if options.adjust_mapping_quality {
        mapping_quality /= f64::from(observation.local_mismatch_tenths) / 10.0 + 1.0;
        mapping_quality *= 2.0 - depth.min(30) as f64 / 30.0;
    }
    mapping_quality = (mapping_quality * options.mapping_quality_scale).clamp(
        f64::from(options.minimum_mapping_quality),
        f64::from(options.maximum_mapping_quality),
    );
    let base_error = BASE_ERROR[usize::from(quality)];
    let mapping_error = if observation.mapping_quality == 255 {
        MAPPING_ERROR[10]
    } else {
        MAPPING_ERROR[mapping_quality.trunc() as usize]
    };
    phred_log(base_error + 0.75 * mapping_error - base_error * mapping_error) as u8
}

fn build_zero_local_quality_cache(options: &BayesianOptions) -> Option<Box<[u8]>> {
    if !options.use_mapping_quality {
        return None;
    }
    let mut cache = vec![0; CACHED_DEPTHS * CACHED_MAPPING_QUALITIES * CACHED_BASE_QUALITIES];
    for depth in 0..CACHED_DEPTHS {
        for mapping_quality in 0..CACHED_MAPPING_QUALITIES {
            let observation = BayesianObservation::new(0, 0).with_record_context(
                mapping_quality as u8,
                0,
                0,
                false,
                false,
            );
            for quality in 0..CACHED_BASE_QUALITIES {
                cache[quality_cache_index(depth, mapping_quality as u8, quality as u8)] =
                    adjusted_quality(quality as u8, &observation, depth, options);
            }
        }
    }
    Some(cache.into_boxed_slice())
}

fn quality_cache_index(depth: usize, mapping_quality: u8, quality: u8) -> usize {
    (depth.min(CACHED_DEPTHS - 1) * CACHED_MAPPING_QUALITIES + usize::from(mapping_quality))
        * CACHED_BASE_QUALITIES
        + usize::from(quality)
}

#[allow(clippy::too_many_arguments)]
fn accumulate(
    scores: &mut [f64; 15],
    base: usize,
    both_match: f64,
    one_match: f64,
    both_overcall: f64,
    one_overcall: f64,
    overcall_mismatch: f64,
    both_undercall: f64,
    one_undercall: f64,
    both_gap_match: f64,
) {
    match base {
        0 => {
            scores[0] += both_match;
            scores[1] += one_match;
            scores[2] += one_match;
            scores[3] += one_match;
            scores[4] += one_overcall;
            scores[8] += overcall_mismatch;
            scores[11] += overcall_mismatch;
            scores[13] += overcall_mismatch;
            scores[14] += both_overcall;
        }
        1 => {
            scores[1] += one_match;
            scores[5] += both_match;
            scores[6] += one_match;
            scores[7] += one_match;
            scores[8] += one_overcall;
            scores[4] += overcall_mismatch;
            scores[11] += overcall_mismatch;
            scores[13] += overcall_mismatch;
            scores[14] += both_overcall;
        }
        2 => {
            scores[2] += one_match;
            scores[6] += one_match;
            scores[9] += both_match;
            scores[10] += one_match;
            scores[11] += one_overcall;
            scores[4] += overcall_mismatch;
            scores[8] += overcall_mismatch;
            scores[13] += overcall_mismatch;
            scores[14] += both_overcall;
        }
        3 => {
            scores[3] += one_match;
            scores[7] += one_match;
            scores[10] += one_match;
            scores[12] += both_match;
            scores[13] += one_overcall;
            scores[4] += overcall_mismatch;
            scores[8] += overcall_mismatch;
            scores[11] += overcall_mismatch;
            scores[14] += both_overcall;
        }
        4 => {
            scores[0] += both_undercall;
            scores[1] += both_undercall;
            scores[2] += both_undercall;
            scores[3] += both_undercall;
            scores[4] += one_undercall;
            scores[5] += both_undercall;
            scores[6] += both_undercall;
            scores[7] += both_undercall;
            scores[8] += one_undercall;
            scores[9] += both_undercall;
            scores[10] += both_undercall;
            scores[11] += one_undercall;
            scores[12] += both_undercall;
            scores[13] += one_undercall;
            scores[14] += both_gap_match;
        }
        5 => {
            scores[0] += both_match;
            scores[1] += both_match;
            scores[2] += both_match;
            scores[3] += both_match;
            scores[4] += one_overcall;
            scores[5] += both_match;
            scores[6] += both_match;
            scores[7] += both_match;
            scores[8] += one_overcall;
            scores[9] += both_match;
            scores[10] += both_match;
            scores[11] += one_overcall;
            scores[12] += both_match;
            scores[13] += one_overcall;
            scores[14] += both_overcall;
        }
        _ => unreachable!(),
    }
}

fn map_base(base: u8) -> usize {
    const MAP: [usize; 32] = [
        5, 0, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4,
    ];
    MAP[usize::from(base.min(31))]
}

fn fast_exp(value: f64) -> f64 {
    static TENTHS: LazyLock<[f64; 1001]> =
        LazyLock::new(|| std::array::from_fn(|index| (f64::from(index as i32 - 500) / 10.0).exp()));
    static INTEGERS: LazyLock<[f64; 1001]> =
        LazyLock::new(|| std::array::from_fn(|index| f64::from(index as i32 - 500).exp()));
    if (-50.0..=50.0).contains(&value) {
        TENTHS[((value * 10.0) as i32 + 500) as usize]
    } else {
        INTEGERS[(value.clamp(-500.0, 500.0) as i32 + 500) as usize]
    }
}

fn fast_log2(value: f64) -> f64 {
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 2047) as i32 - 1024;
    let mantissa = f64::from_bits((bits & !(2047_u64 << 52)) + (1023_u64 << 52));
    f64::from(exponent) + ((-1.0 / 3.0 * mantissa + 2.0) * mantissa - 2.0 / 3.0)
}

fn phred_log(value: f64) -> f64 {
    -TEN_LOG2_OVER_LOG10 * fast_log2(value)
}
