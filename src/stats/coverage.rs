use std::collections::BTreeMap;
use std::ops::Bound::{Excluded, Unbounded};

use serde::Serialize;

use super::CoverageBins;

#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize)]
pub(crate) struct Coverage {
    reference: Option<i32>,
    position: i64,
    depth: i64,
    changes: BTreeMap<i64, i64>,
    lengths: BTreeMap<usize, u64>,
}

impl Coverage {
    pub(crate) fn advance(&mut self, reference: i32, position: i64) {
        if self.reference != Some(reference) {
            self.finish_reference();
            self.reference = Some(reference);
            self.position = position;
            return;
        }
        self.flush_to(position);
    }

    pub(crate) fn insert(&mut self, start: i64, end: i64) {
        if end <= start {
            return;
        }
        if start == self.position {
            self.depth += 1;
        } else {
            *self.changes.entry(start).or_default() += 1;
        }
        *self.changes.entry(end).or_default() -= 1;
    }

    pub(crate) fn histogram(&self, bins: CoverageBins) -> Vec<u64> {
        let lengths = self.finalized_lengths();
        let mut histogram = vec![0; bins.count()];
        for (depth, length) in lengths {
            histogram[bins.index(depth)] += length;
        }
        histogram
    }

    pub(crate) fn bases_above(&self, threshold: usize) -> u64 {
        self.finalized_lengths()
            .range((Excluded(threshold), Unbounded))
            .map(|(_, length)| length)
            .sum()
    }

    fn flush_to(&mut self, end: i64) {
        while self
            .changes
            .first_key_value()
            .is_some_and(|(&position, _)| position <= end)
        {
            let (position, change) = self.changes.pop_first().unwrap();
            self.record_length(position);
            self.depth += change;
        }
        self.record_length(end);
    }

    fn record_length(&mut self, end: i64) {
        if self.depth > 0 && end > self.position {
            let depth = self.depth as usize;
            *self.lengths.entry(depth).or_default() += u64::try_from(end - self.position).unwrap();
        }
        self.position = end;
    }

    fn finish_reference(&mut self) {
        while let Some((position, change)) = self.changes.pop_first() {
            self.record_length(position);
            self.depth += change;
        }
        self.reference = None;
        self.position = 0;
        self.depth = 0;
    }

    fn finalized_lengths(&self) -> BTreeMap<usize, u64> {
        let mut coverage = self.clone();
        coverage.finish_reference();
        coverage.lengths
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn streams_changes_into_depth_lengths() {
        let mut coverage = Coverage::default();
        coverage.advance(0, 10);
        coverage.insert(10, 20);
        coverage.advance(0, 15);
        coverage.insert(15, 25);
        coverage.advance(0, 30);

        assert_eq!(coverage.histogram("1,10,1".parse().unwrap())[1..3], [10, 5]);
        assert_eq!(coverage.bases_above(1), 5);
    }

    #[test]
    fn flushes_between_references() {
        let mut coverage = Coverage::default();
        coverage.advance(0, 0);
        coverage.insert(0, 5);
        coverage.advance(1, 0);
        coverage.insert(0, 7);

        assert_eq!(coverage.histogram("1,10,1".parse().unwrap())[1], 12);
    }

    #[test]
    fn stores_extreme_depth_sparsely() {
        let mut coverage = Coverage {
            position: 0,
            depth: i64::MAX,
            ..Coverage::default()
        };
        coverage.record_length(1);

        let histogram = coverage.histogram("1,1000,1".parse().unwrap());
        assert_eq!(histogram.last(), Some(&1));
        assert_eq!(coverage.bases_above(1000), 1);
        assert_eq!(coverage.lengths.len(), 1);
    }
}
