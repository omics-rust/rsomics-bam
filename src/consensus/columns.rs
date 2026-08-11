use rsomics_pileup::{Column, ColumnEntry};

use super::{call::BayesianObservation, record::RecordState};

pub(super) fn reference_observations(
    column: &Column<'_, RecordState>,
    observations: &mut Vec<BayesianObservation>,
) {
    observations.clear();
    observations.reserve(column.len());
    for entry in column.entries() {
        let projection = entry.projection();
        let (base, quality, context_position) = if projection.is_reference_skip {
            (0, 0, projection.qpos)
        } else if projection.is_deletion {
            let previous = entry.state().quality(projection.qpos.saturating_sub(1));
            let next = entry.state().quality(projection.qpos);
            (16, previous.min(next), projection.qpos)
        } else {
            (
                entry.record().seq_nibble(projection.qpos),
                entry.state().quality(projection.qpos),
                projection.qpos + 1,
            )
        };
        observations.push(entry.state().observation(
            entry.record(),
            base,
            quality,
            context_position,
            projection.is_reference_skip,
        ));
    }
}

pub(super) fn insertion_width(column: &Column<'_, RecordState>) -> usize {
    column
        .entries()
        .map(|entry| {
            let projection = entry.projection();
            if projection.indel <= 0 && !entry.state().has_padding() {
                return 0;
            }
            let padding_only = projection.indel == 0
                && entry
                    .cigar()
                    .nth(projection.cigar_index + 1)
                    .is_some_and(|(kind, _)| kind == 6);
            if projection.indel <= 0 && (!padding_only || !at_cigar_end(entry, column.position())) {
                return 0;
            }
            entry
                .cigar()
                .skip(projection.cigar_index + 1)
                .take_while(|(kind, _)| matches!(kind, 1 | 6))
                .map(|(_, length)| length as usize)
                .sum()
        })
        .max()
        .unwrap_or(0)
}

pub(super) fn insertion_observations(
    column: &Column<'_, RecordState>,
    offset: usize,
    observations: &mut Vec<BayesianObservation>,
) {
    observations.clear();
    observations.reserve(column.len());
    for entry in column.entries() {
        observations.push(insertion_observation(entry, column.position(), offset));
    }
}

pub(super) fn visit_observation_columns<E>(
    column: &Column<'_, RecordState>,
    observations: &mut Vec<BayesianObservation>,
    mut visit: impl FnMut(usize, &[BayesianObservation]) -> Result<(), E>,
) -> Result<(), E> {
    reference_observations(column, observations);
    visit(0, observations)?;
    for offset in 0..insertion_width(column) {
        insertion_observations(column, offset, observations);
        visit(offset + 1, observations)?;
    }
    Ok(())
}

fn insertion_observation(
    entry: ColumnEntry<'_, RecordState>,
    column_position: i64,
    offset: usize,
) -> BayesianObservation {
    let projection = entry.projection();
    let query_start = projection.qpos + usize::from(!projection.is_deletion);
    let mut query_position = query_start;
    let mut insertion_position = 0usize;
    let mut inserted_bases = 0usize;
    let mut observed = None;
    if at_cigar_end(entry, column_position) {
        for (kind, length) in entry.cigar().skip(projection.cigar_index + 1) {
            if !matches!(kind, 1 | 6) {
                break;
            }
            let length = length as usize;
            if observed.is_none()
                && kind == 1
                && offset >= insertion_position
                && offset < insertion_position + length
            {
                let position = query_position + offset - insertion_position;
                observed = Some((
                    entry.record().seq_nibble(position),
                    entry.state().quality(position),
                    position + 1,
                ));
            }
            insertion_position += length;
            if kind == 1 {
                query_position += length;
                inserted_bases += length;
            }
        }
    }
    let (base, quality, query_position) = observed.unwrap_or_else(|| {
        let previous = entry.state().quality(query_start.saturating_sub(1));
        let next = entry.state().quality(query_start + inserted_bases);
        (16, previous.min(next), query_start)
    });
    entry.state().observation(
        entry.record(),
        base,
        quality,
        query_position,
        projection.is_reference_skip,
    )
}

fn at_cigar_end(entry: ColumnEntry<'_, RecordState>, column_position: i64) -> bool {
    let reference_length = entry
        .cigar()
        .take(entry.projection().cigar_index + 1)
        .filter(|(kind, _)| matches!(kind, 0 | 2 | 3 | 7 | 8))
        .map(|(_, length)| i64::from(length))
        .sum::<i64>();
    i64::from(entry.record().alignment_start()) + reference_length - 1 == column_position
}

#[cfg(test)]
mod tests {
    use rsomics_pileup::{PileupEngine, PileupOptions};

    use super::*;
    use crate::consensus::{
        call::{BayesianCaller, BayesianOptions},
        record::{RecordOptions, test_record, test_record_with},
    };

    #[test]
    fn retained_record_state_drives_bayesian_columns() {
        let record = test_record(b"10");
        let state = RecordState::new(
            &record,
            RecordOptions {
                adjust_quality: false,
                ..RecordOptions::default()
            },
        )
        .unwrap();
        let mut pileup = PileupEngine::with_record_state([100], PileupOptions::default());
        pileup.push_with_state(record, state).unwrap();
        pileup.finish().unwrap();
        let caller = BayesianCaller::new(BayesianOptions::default());
        let mut observations = Vec::new();
        let mut sequence = Vec::new();

        pileup
            .drain(|column| {
                reference_observations(column, &mut observations);
                sequence.push(caller.call(&observations).base);
                Ok::<_, ()>(())
            })
            .unwrap();

        assert_eq!(sequence, b"AAACCCCCGT");
    }

    #[test]
    fn expands_insertions_and_pads_shorter_records() {
        let inserted = test_record_with(b"ACGGGTACGTACG", &[(0, 2), (1, 3), (0, 8)], b"10");
        let plain = test_record(b"10");
        let options = RecordOptions {
            adjust_quality: false,
            ..RecordOptions::default()
        };
        let inserted_state = RecordState::new(&inserted, options).unwrap();
        let plain_state = RecordState::new(&plain, options).unwrap();
        let mut pileup = PileupEngine::with_record_state([100], PileupOptions::default());
        pileup.push_with_state(inserted, inserted_state).unwrap();
        pileup.push_with_state(plain, plain_state).unwrap();
        pileup.finish().unwrap();
        let mut observations = Vec::new();

        pileup
            .drain(|column| {
                if column.position() == 0 {
                    assert_eq!(insertion_width(column), 0);
                }
                if column.position() == 1 {
                    assert_eq!(insertion_width(column), 3);
                    for offset in 0..3 {
                        insertion_observations(column, offset, &mut observations);
                        assert_eq!(
                            observations
                                .iter()
                                .map(|entry| entry.base)
                                .collect::<Vec<_>>(),
                            [4, 16]
                        );
                        assert_eq!(
                            observations
                                .iter()
                                .map(|entry| entry.quality)
                                .collect::<Vec<_>>(),
                            [30, 30]
                        );
                    }
                }
                Ok::<_, ()>(())
            })
            .unwrap();
    }

    #[test]
    fn preserves_cigar_padding_as_an_insertion_gap() {
        let padded = test_record_with(b"ACGGGTACGTACG", &[(0, 2), (6, 1), (1, 3), (0, 8)], b"10");
        let plain = test_record(b"10");
        let options = RecordOptions {
            adjust_quality: false,
            ..RecordOptions::default()
        };
        let padded_state = RecordState::new(&padded, options).unwrap();
        let plain_state = RecordState::new(&plain, options).unwrap();
        let mut pileup = PileupEngine::with_record_state([100], PileupOptions::default());
        pileup.push_with_state(padded, padded_state).unwrap();
        pileup.push_with_state(plain, plain_state).unwrap();
        pileup.finish().unwrap();
        let mut observations = Vec::new();

        pileup
            .drain(|column| {
                if column.position() == 1 {
                    assert_eq!(insertion_width(column), 4);
                    for (offset, expected) in [[16, 16], [4, 16], [4, 16], [4, 16]]
                        .into_iter()
                        .enumerate()
                    {
                        insertion_observations(column, offset, &mut observations);
                        assert_eq!(
                            observations
                                .iter()
                                .map(|entry| entry.base)
                                .collect::<Vec<_>>(),
                            expected
                        );
                    }
                }
                Ok::<_, ()>(())
            })
            .unwrap();
    }

    #[test]
    fn reads_insertions_after_a_deletion_from_the_next_query_base() {
        let inserted =
            test_record_with(b"ACTTTACGTACG", &[(0, 2), (2, 1), (1, 3), (0, 7)], b"2^A7");
        let plain = test_record(b"10");
        let options = RecordOptions {
            adjust_quality: false,
            ..RecordOptions::default()
        };
        let inserted_state = RecordState::new(&inserted, options).unwrap();
        let plain_state = RecordState::new(&plain, options).unwrap();
        let mut pileup = PileupEngine::with_record_state([100], PileupOptions::default());
        pileup.push_with_state(inserted, inserted_state).unwrap();
        pileup.push_with_state(plain, plain_state).unwrap();
        pileup.finish().unwrap();
        let mut observations = Vec::new();

        pileup
            .drain(|column| {
                if column.position() == 2 {
                    assert_eq!(insertion_width(column), 3);
                    for offset in 0..3 {
                        insertion_observations(column, offset, &mut observations);
                        assert_eq!(
                            observations
                                .iter()
                                .map(|entry| entry.base)
                                .collect::<Vec<_>>(),
                            [8, 16]
                        );
                    }
                }
                Ok::<_, ()>(())
            })
            .unwrap();
    }

    #[test]
    fn visits_reference_and_insertion_columns_in_output_order() {
        let inserted = test_record_with(b"ACGGGTACGTACG", &[(0, 2), (1, 3), (0, 8)], b"10");
        let state = RecordState::new(
            &inserted,
            RecordOptions {
                adjust_quality: false,
                ..RecordOptions::default()
            },
        )
        .unwrap();
        let mut pileup = PileupEngine::with_record_state([100], PileupOptions::default());
        pileup.push_with_state(inserted, state).unwrap();
        pileup.finish().unwrap();
        let mut observations = Vec::new();
        let mut visited = Vec::new();

        pileup
            .drain(|column| {
                if column.position() == 1 {
                    visit_observation_columns(column, &mut observations, |offset, entries| {
                        visited.push((offset, entries[0].base));
                        Ok::<_, ()>(())
                    })?;
                }
                Ok::<_, ()>(())
            })
            .unwrap();

        assert_eq!(visited, [(0, 2), (1, 4), (2, 4), (3, 4)]);
    }
}
