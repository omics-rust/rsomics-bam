use rsomics_pileup::Column;

use super::{
    call::{BayesianObservation, Caller},
    columns::visit_observation_columns,
    record::RecordState,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct CalledColumn {
    pub(super) reference_id: i32,
    pub(super) position: i64,
    pub(super) offset: usize,
    pub(super) depth: usize,
    pub(super) base: u8,
    pub(super) quality: i32,
}

pub(super) struct Walker {
    caller: Caller,
    observations: Vec<BayesianObservation>,
}

impl Walker {
    pub(super) fn new(caller: Caller) -> Self {
        Self {
            caller,
            observations: Vec::new(),
        }
    }

    pub(super) fn visit<E>(
        &mut self,
        column: &Column<'_, RecordState>,
        mut emit: impl FnMut(CalledColumn, &[BayesianObservation]) -> Result<(), E>,
    ) -> Result<(), E> {
        let reference_id = column.reference_id();
        let position = column.position();
        let depth = column.len();
        let caller = &self.caller;
        visit_observation_columns(column, &mut self.observations, |offset, observations| {
            let call = caller.call(observations);
            emit(
                CalledColumn {
                    reference_id,
                    position,
                    offset,
                    depth,
                    base: call.base,
                    quality: call.quality,
                },
                observations,
            )
        })
    }
}

#[cfg(test)]
mod tests {
    use rsomics_pileup::{PileupEngine, PileupOptions};

    use super::*;
    use crate::consensus::{
        call::{Caller, SimpleOptions},
        record::{RecordOptions, RecordState, test_record_with},
    };

    #[test]
    fn calls_reference_and_insertion_columns() {
        let options = RecordOptions {
            adjust_quality: false,
            ..RecordOptions::default()
        };
        let mut pileup = PileupEngine::with_record_state([100], PileupOptions::default());
        for _ in 0..3 {
            let record = test_record_with(b"ACGGGTACGTACG", &[(0, 2), (1, 3), (0, 8)], b"10");
            let state = RecordState::new(&record, options).unwrap();
            pileup.push_with_state(record, state).unwrap();
        }
        for _ in 0..2 {
            let record = test_record_with(b"ACTACGTACG", &[(0, 10)], b"10");
            let state = RecordState::new(&record, options).unwrap();
            pileup.push_with_state(record, state).unwrap();
        }
        pileup.finish().unwrap();
        let mut walker = Walker::new(Caller::Simple(SimpleOptions {
            use_quality: false,
            minimum_quality: 0,
            minimum_depth: 1,
            call_fraction: 0.6,
            heterozygous_fraction: 0.5,
            ambiguous: false,
        }));
        let mut calls = Vec::new();

        pileup
            .drain(|column| {
                if column.position() == 1 {
                    walker.visit(column, |call, _| {
                        calls.push(call);
                        Ok::<_, ()>(())
                    })?;
                }
                Ok::<_, ()>(())
            })
            .unwrap();

        assert_eq!(
            calls,
            [
                CalledColumn {
                    reference_id: 0,
                    position: 1,
                    offset: 0,
                    depth: 5,
                    base: b'C',
                    quality: 100,
                },
                CalledColumn {
                    reference_id: 0,
                    position: 1,
                    offset: 1,
                    depth: 5,
                    base: b'G',
                    quality: 60,
                },
                CalledColumn {
                    reference_id: 0,
                    position: 1,
                    offset: 2,
                    depth: 5,
                    base: b'G',
                    quality: 60,
                },
                CalledColumn {
                    reference_id: 0,
                    position: 1,
                    offset: 3,
                    depth: 5,
                    base: b'G',
                    quality: 60,
                },
            ]
        );
    }
}
