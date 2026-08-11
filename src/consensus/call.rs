mod bayesian;

pub(crate) use bayesian::CalibrationPreset;
pub(super) use bayesian::{BayesianCaller, BayesianMode, BayesianObservation, BayesianOptions};

#[derive(Clone, Copy)]
pub(super) struct SimpleOptions {
    pub(super) use_quality: bool,
    pub(super) minimum_quality: u8,
    pub(super) minimum_depth: usize,
    pub(super) call_fraction: f64,
    pub(super) heterozygous_fraction: f64,
    pub(super) ambiguous: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct Call {
    pub(super) base: u8,
    pub(super) quality: i32,
}

pub(super) enum Caller {
    Simple(SimpleOptions),
    Bayesian(Box<BayesianCaller>),
}

impl Caller {
    pub(super) fn call(&self, observations: &[BayesianObservation]) -> Call {
        match self {
            Self::Simple(options) => simple(observations, *options),
            Self::Bayesian(caller) => caller.call(observations),
        }
    }
}

pub(super) fn simple(observations: &[BayesianObservation], options: SimpleOptions) -> Call {
    const A: [u8; 16] = [0, 8, 0, 4, 0, 4, 0, 2, 0, 4, 0, 2, 0, 2, 0, 1];
    const C: [u8; 16] = [0, 0, 8, 4, 0, 0, 4, 2, 0, 0, 4, 2, 0, 0, 2, 1];
    const G: [u8; 16] = [0, 0, 0, 0, 8, 4, 4, 1, 0, 0, 0, 0, 4, 2, 2, 1];
    const T: [u8; 16] = [0, 0, 0, 0, 0, 0, 0, 0, 8, 4, 4, 2, 8, 2, 2, 1];
    const HETEROZYGOUS: &[u8; 32] = b"NACMGRSVTWYHKDBN*ac?g???t???????";

    let mut scores = [0u64; 5];
    let mut depth = 0usize;
    for observation in observations {
        if observation.quality < options.minimum_quality {
            continue;
        }
        let quality = if options.use_quality {
            u64::from(observation.quality)
        } else {
            1
        };
        if observation.base < 16 {
            let base = usize::from(observation.base);
            for (score, compatibility) in scores[..4]
                .iter_mut()
                .zip([A[base], C[base], G[base], T[base]])
            {
                *score += u64::from(compatibility) * quality;
            }
        } else {
            scores[4] += 8 * quality;
        }
        depth += 1;
    }

    let total = scores.iter().sum::<u64>();
    let mut first = 15usize;
    let mut second = 15usize;
    let mut first_score = 0u64;
    let mut second_score = 0u64;
    for (index, &score) in scores.iter().enumerate() {
        let base = 1usize << index;
        if score > first_score {
            second = first;
            second_score = first_score;
            first = base;
            first_score = score;
        } else if score > second_score {
            second = base;
            second_score = score;
        }
    }

    let mut called = first;
    let mut called_score = first_score;
    if options.ambiguous
        && first_score != 0
        && second_score as f64 >= options.heterozygous_fraction * first_score as f64
    {
        called |= second;
        called_score += second_score;
    }
    if depth < options.minimum_depth || (called_score as f64) < options.call_fraction * total as f64
    {
        called = if first == 16 { 16 } else { 0 };
    }

    let quality = if called == 0 || total == 0 {
        0
    } else {
        ((u128::from(called_score) * 100) / u128::from(total)).min(100) as i32
    };
    Call {
        base: HETEROZYGOUS[called],
        quality,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nt16(base: u8) -> u8 {
        match base.to_ascii_uppercase() {
            b'A' => 1,
            b'C' => 2,
            b'G' => 4,
            b'T' => 8,
            b'*' | b'#' => 16,
            _ => 15,
        }
    }

    fn options() -> SimpleOptions {
        SimpleOptions {
            use_quality: false,
            minimum_quality: 0,
            minimum_depth: 1,
            call_fraction: 0.75,
            heterozygous_fraction: 0.5,
            ambiguous: false,
        }
    }

    fn observation(base: u8, quality: u8) -> BayesianObservation {
        BayesianObservation::new(base, quality)
    }

    #[test]
    fn calls_dominant_base_and_confidence() {
        let observations = [
            observation(1, 30),
            observation(1, 20),
            observation(1, 10),
            observation(2, 40),
        ];

        assert_eq!(
            simple(&observations, options()),
            Call {
                base: b'A',
                quality: 75
            }
        );
    }

    #[test]
    fn quality_weighting_changes_the_call() {
        let observations = [
            observation(1, 40),
            observation(2, 5),
            observation(2, 5),
            observation(2, 5),
        ];
        let mut weighted = options();
        weighted.use_quality = true;
        weighted.call_fraction = 0.7;

        assert_eq!(simple(&observations, options()).base, b'C');
        assert_eq!(simple(&observations, weighted).base, b'A');
    }

    #[test]
    fn ambiguity_and_gap_use_upstream_codes() {
        let mut ambiguous = options();
        ambiguous.ambiguous = true;
        ambiguous.call_fraction = 1.0;
        let snp = [observation(1, 20), observation(4, 20)];
        let deletion = [observation(4, 20), observation(16, 20)];

        assert_eq!(simple(&snp, ambiguous).base, b'R');
        assert_eq!(simple(&deletion, ambiguous).base, b'g');
    }

    #[test]
    fn depth_and_quality_filters_preserve_a_gap_call() {
        let observations = [observation(16, 5)];
        let mut shallow = options();
        shallow.minimum_depth = 2;
        let mut filtered = options();
        filtered.minimum_quality = 6;

        assert_eq!(simple(&observations, shallow).base, b'*');
        assert_eq!(simple(&observations, filtered).base, b'N');
    }

    #[test]
    fn nt16_ambiguity_contributes_to_each_compatible_base() {
        let observations = [observation(3, 30)];
        let mut ambiguous = options();
        ambiguous.ambiguous = true;

        assert_eq!(simple(&observations, ambiguous).base, b'M');
        assert_eq!(simple(&observations, ambiguous).quality, 100);
    }

    #[test]
    fn caller_dispatches_both_models() {
        let observations = [observation(1, 30)];
        let simple = Caller::Simple(options());
        let bayesian_options = BayesianOptions {
            use_mapping_quality: false,
            cutoff: 0,
            ..BayesianOptions::default()
        };
        let bayesian = Caller::Bayesian(Box::new(BayesianCaller::new(bayesian_options)));

        assert_eq!(simple.call(&observations).base, b'A');
        assert_eq!(bayesian.call(&observations).base, b'A');
    }

    fn assert_bayesian_oracle(oracle: &str, options: BayesianOptions) {
        let caller = BayesianCaller::new(options);
        for line in oracle.lines() {
            let fields = line.split('\t').collect::<Vec<_>>();
            let bases = fields[6].as_bytes();
            let qualities = fields[7].as_bytes();
            let observations = bases
                .iter()
                .zip(qualities)
                .map(|(&base, &quality)| BayesianObservation::new(nt16(base), quality - b'!'))
                .collect::<Vec<_>>();
            let expected = Call {
                base: fields[4].as_bytes()[0],
                quality: fields[5].parse().unwrap(),
            };

            assert_eq!(
                caller.call(&observations),
                expected,
                "{}:{}:{}",
                fields[0],
                fields[1],
                fields[2]
            );
        }
    }

    #[test]
    fn bayesian_without_mapping_quality_matches_samtools_1_24_columns() {
        for (oracle, cutoff, ambiguous) in [
            (
                include_str!("../../tests/upstream/samtools-consensus/expected/18p.out"),
                0,
                false,
            ),
            (
                include_str!("../../tests/upstream/samtools-consensus/expected/19p.out"),
                19,
                false,
            ),
            (
                include_str!("../../tests/upstream/samtools-consensus/expected/20p.out"),
                30,
                true,
            ),
            (
                include_str!("../../tests/upstream/samtools-consensus/expected/21p.out"),
                31,
                true,
            ),
        ] {
            assert_bayesian_oracle(
                oracle,
                BayesianOptions {
                    use_mapping_quality: false,
                    cutoff,
                    ambiguous,
                    ..BayesianOptions::default()
                },
            );
        }
    }
}
