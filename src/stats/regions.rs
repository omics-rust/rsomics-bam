use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::Path;

use noodles::core::Region;
use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Debug)]
pub(crate) struct Regions {
    intervals: HashMap<i32, Vec<Range<i64>>>,
    bases: u64,
}

impl Regions {
    pub(crate) fn from_targets(path: &Path, references: &[(Vec<u8>, u64)]) -> Result<Self> {
        let names = references
            .iter()
            .enumerate()
            .map(|(index, (name, _))| (name.clone(), index as i32))
            .collect::<HashMap<_, _>>();
        let input = File::open(path).map_err(|error| {
            RsomicsError::Io(std::io::Error::new(
                error.kind(),
                format!("opening target regions {}: {error}", path.display()),
            ))
        })?;
        let mut intervals = HashMap::<i32, Vec<Range<i64>>>::new();
        for (line_number, result) in BufReader::new(input).lines().enumerate() {
            let line = result.map_err(RsomicsError::Io)?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields = line.split_ascii_whitespace().collect::<Vec<_>>();
            if fields.len() < 3 {
                return Err(invalid_target(
                    path,
                    line_number + 1,
                    "expected name, start, and end",
                ));
            }
            let Some(&reference) = names.get(fields[0].as_bytes()) else {
                continue;
            };
            let start = fields[1]
                .parse::<u64>()
                .map_err(|_| invalid_target(path, line_number + 1, "invalid start"))?;
            let end = fields[2]
                .parse::<u64>()
                .map_err(|_| invalid_target(path, line_number + 1, "invalid end"))?;
            if start == 0 || start > end {
                return Err(invalid_target(
                    path,
                    line_number + 1,
                    "coordinates must be 1-based and inclusive",
                ));
            }
            let reference_length = references[reference as usize].1;
            if start > reference_length {
                continue;
            }
            intervals.entry(reference).or_default().push(
                i64::try_from(start - 1).unwrap()
                    ..i64::try_from(end.min(reference_length)).unwrap(),
            );
        }
        Self::finish(intervals)
    }

    pub(crate) fn from_cli(values: &[String], references: &[(Vec<u8>, u64)]) -> Result<Self> {
        let names = references
            .iter()
            .enumerate()
            .map(|(index, (name, _))| (name.clone(), index as i32))
            .collect::<HashMap<_, _>>();
        let mut intervals = HashMap::<i32, Vec<Range<i64>>>::new();
        for value in values {
            let region = value.parse::<Region>().map_err(|error| {
                RsomicsError::ConfigError(format!("invalid region {value}: {error}"))
            })?;
            let reference = *names.get(region.name()).ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "region reference is absent from the alignment header: {}",
                    String::from_utf8_lossy(region.name())
                ))
            })?;
            let reference_length = references[reference as usize].1;
            let start = region
                .interval()
                .start()
                .map(|position| usize::from(position) as u64 - 1)
                .unwrap_or(0);
            let end = region
                .interval()
                .end()
                .map(|position| usize::from(position) as u64)
                .unwrap_or(reference_length)
                .min(reference_length);
            if start >= end {
                return Err(RsomicsError::InvalidInput(format!(
                    "region {value} starts outside its reference"
                )));
            }
            intervals
                .entry(reference)
                .or_default()
                .push(i64::try_from(start).unwrap()..i64::try_from(end).unwrap());
        }
        Self::finish(intervals)
    }

    pub(crate) fn intersect(self, other: Self) -> Result<Self> {
        let mut intervals = HashMap::new();
        for (reference, left) in self.intervals {
            let Some(right) = other.intervals.get(&reference) else {
                continue;
            };
            let mut output = Vec::new();
            for a in &left {
                for b in right {
                    let start = a.start.max(b.start);
                    let end = a.end.min(b.end);
                    if start < end {
                        output.push(start..end);
                    }
                }
            }
            if !output.is_empty() {
                intervals.insert(reference, output);
            }
        }
        Self::finish(intervals)
    }

    pub(crate) fn bases(&self) -> u64 {
        self.bases
    }

    pub(crate) fn overlaps(&self, reference: i32, start: i64, end: i64) -> bool {
        self.intervals.get(&reference).is_some_and(|ranges| {
            let index = ranges.partition_point(|range| range.end <= start);
            ranges.get(index).is_some_and(|range| range.start < end)
        })
    }

    pub(crate) fn intersections(
        &self,
        reference: i32,
        start: i64,
        end: i64,
    ) -> impl Iterator<Item = Range<i64>> + '_ {
        self.intervals
            .get(&reference)
            .into_iter()
            .flatten()
            .skip_while(move |range| range.end <= start)
            .take_while(move |range| range.start < end)
            .map(move |range| range.start.max(start)..range.end.min(end))
    }

    pub(crate) fn first_overlap(&self, reference: i32, start: i64, end: i64) -> Option<Range<i64>> {
        let ranges = self.intervals.get(&reference)?;
        let index = ranges.partition_point(|range| range.end <= start);
        ranges.get(index).filter(|range| range.start < end).cloned()
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = (i32, Range<i64>)> + '_ {
        let mut references = self.intervals.iter().collect::<Vec<_>>();
        references.sort_unstable_by_key(|(reference, _)| **reference);
        references.into_iter().flat_map(|(&reference, ranges)| {
            ranges.iter().cloned().map(move |range| (reference, range))
        })
    }

    fn finish(mut intervals: HashMap<i32, Vec<Range<i64>>>) -> Result<Self> {
        for ranges in intervals.values_mut() {
            ranges.sort_unstable_by_key(|range| (range.start, range.end));
            let mut merged = Vec::<Range<i64>>::with_capacity(ranges.len());
            for range in ranges.drain(..) {
                if let Some(previous) = merged.last_mut()
                    && range.start <= previous.end
                {
                    previous.end = previous.end.max(range.end);
                } else {
                    merged.push(range);
                }
            }
            *ranges = merged;
        }
        let bases = intervals
            .values()
            .flatten()
            .try_fold(0u64, |total, range| {
                let length = u64::try_from(range.end - range.start).map_err(|_| {
                    RsomicsError::InvalidInput("target region length overflows".to_owned())
                })?;
                total.checked_add(length).ok_or_else(|| {
                    RsomicsError::InvalidInput("total target length overflows".to_owned())
                })
            })?;
        if bases == 0 {
            return Err(RsomicsError::InvalidInput(
                "no selected region maps to the alignment header".to_owned(),
            ));
        }
        Ok(Self { intervals, bases })
    }
}

fn invalid_target(path: &Path, line: usize, reason: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{}:{line}: {reason}", path.display()))
}
