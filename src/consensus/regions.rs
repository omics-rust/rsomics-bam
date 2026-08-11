use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use rsomics_common::{Result, RsomicsError};

#[derive(Debug, Eq, PartialEq)]
pub(super) struct Interval {
    pub(super) name: Box<[u8]>,
    pub(super) start: u64,
    pub(super) end: u64,
}

struct Entry {
    name: Box<[u8]>,
    intervals: Vec<(u64, u64)>,
}

#[derive(Default)]
struct SamtoolsTable {
    entries: Vec<Entry>,
    buckets: Vec<Option<usize>>,
    upper_bound: usize,
}

impl SamtoolsTable {
    fn push(&mut self, name: &[u8], start: u64, end: u64) {
        if let Some(entry) = self
            .entries
            .iter_mut()
            .find(|entry| entry.name.as_ref() == name)
        {
            entry.intervals.push((start, end));
            return;
        }
        if self.entries.len() >= self.upper_bound {
            self.resize(self.buckets.len() + 1);
        }
        let entry = self.entries.len();
        self.entries.push(Entry {
            name: name.into(),
            intervals: vec![(start, end)],
        });
        let bucket = find_empty(&self.buckets, hash(name));
        self.buckets[bucket] = Some(entry);
    }

    fn resize(&mut self, requested: usize) {
        let count = requested.next_power_of_two().max(4);
        let mut old = std::mem::take(&mut self.buckets);
        let mut buckets = vec![None; count];
        for index in 0..old.len() {
            let Some(mut entry) = old[index].take() else {
                continue;
            };
            loop {
                let bucket = find_empty(&buckets, hash(&self.entries[entry].name));
                buckets[bucket] = Some(entry);
                let Some(displaced) = old.get_mut(bucket).and_then(Option::take) else {
                    break;
                };
                entry = displaced;
            }
        }
        self.buckets = buckets;
        self.upper_bound = ((count as f64 * 0.77) + 0.5) as usize;
    }

    fn finish(mut self) -> Vec<Interval> {
        let mut output = Vec::new();
        for entry in self.buckets.into_iter().flatten() {
            let entry = &mut self.entries[entry];
            entry.intervals.sort_unstable();
            output.extend(entry.intervals.drain(..).map(|(start, end)| Interval {
                name: entry.name.clone(),
                start,
                end,
            }));
        }
        output
    }
}

pub(super) fn read(path: &Path) -> Result<Vec<Interval>> {
    let mut source = File::open(path).map_err(|error| io_error(path, "opening", error))?;
    let mut magic = [0; 2];
    let bytes = source
        .read(&mut magic)
        .map_err(|error| io_error(path, "reading", error))?;
    drop(source);
    let source = File::open(path).map_err(|error| io_error(path, "opening", error))?;
    let reader: Box<dyn BufRead> = if bytes == 2 && magic == [0x1f, 0x8b] {
        Box::new(BufReader::new(MultiGzDecoder::new(source)))
    } else {
        Box::new(BufReader::new(source))
    };
    let mut table = SamtoolsTable::default();
    for (index, result) in reader.split(b'\n').enumerate() {
        let line = result.map_err(|error| io_error(path, "reading", error))?;
        let line = trim_ascii(&line);
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        let mut fields = line.split(|byte| byte.is_ascii_whitespace());
        let name = fields.next().unwrap();
        if name == b"track" || name == b"browser" {
            continue;
        }
        if name.contains(&0) {
            return Err(invalid(
                path,
                index + 1,
                "reference name contains a NUL byte",
            ));
        }
        let start = parse_coordinate(path, index + 1, fields.next(), "start")?;
        let (start, end) = match fields.next() {
            Some(value) if !value.is_empty() => {
                let end = parse_coordinate(path, index + 1, Some(value), "end")?;
                if end < start {
                    return Err(invalid(path, index + 1, "end precedes start"));
                }
                (start, end)
            }
            _ if start > 0 => (start - 1, start),
            _ => {
                return Err(invalid(
                    path,
                    index + 1,
                    "a one-coordinate region must be 1-based",
                ));
            }
        };
        table.push(name, start, end);
    }
    let intervals = table.finish();
    if intervals.is_empty() {
        return Err(RsomicsError::InvalidInput(format!(
            "region file {} contains no intervals",
            path.display()
        )));
    }
    Ok(intervals)
}

fn find_empty(buckets: &[Option<usize>], hash: u32) -> usize {
    let mask = buckets.len() - 1;
    let mut bucket = hash as usize & mask;
    let mut step = 0usize;
    while buckets[bucket].is_some() {
        step += 1;
        bucket = (bucket + step) & mask;
    }
    bucket
}

fn hash(value: &[u8]) -> u32 {
    value.iter().fold(2_166_136_261, |hash, byte| {
        (hash ^ u32::from(*byte)).wrapping_mul(16_777_619)
    })
}

fn trim_ascii(mut value: &[u8]) -> &[u8] {
    while value.first().is_some_and(u8::is_ascii_whitespace) {
        value = &value[1..];
    }
    while value.last().is_some_and(u8::is_ascii_whitespace) {
        value = &value[..value.len() - 1];
    }
    value
}

fn parse_coordinate(path: &Path, line: usize, value: Option<&[u8]>, field: &str) -> Result<u64> {
    let value = value.ok_or_else(|| invalid(path, line, &format!("missing {field}")))?;
    std::str::from_utf8(value)
        .ok()
        .and_then(|value| value.parse().ok())
        .ok_or_else(|| invalid(path, line, &format!("invalid {field}")))
}

fn invalid(path: &Path, line: usize, reason: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{}:{line}: {reason}", path.display()))
}

fn io_error(path: &Path, action: &str, error: std::io::Error) -> RsomicsError {
    RsomicsError::Io(std::io::Error::new(
        error.kind(),
        format!("{action} region file {}: {error}", path.display()),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn root() -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/upstream/samtools-consensus")
    }

    #[test]
    fn matches_samtools_1_24_bed_order_and_interval_sorting() {
        let intervals = read(&root().join("consen4.bed")).unwrap();
        let actual = intervals
            .iter()
            .map(|interval| {
                (
                    String::from_utf8_lossy(&interval.name).into_owned(),
                    interval.start,
                    interval.end,
                )
            })
            .collect::<Vec<_>>();

        assert_eq!(
            actual,
            [
                ("c2".to_owned(), 0, 2),
                ("c2".to_owned(), 8, 10),
                ("c1".to_owned(), 0, 2),
                ("c1".to_owned(), 0, 10),
                ("c1".to_owned(), 8, 15),
                ("c3".to_owned(), 0, 2),
                ("c3".to_owned(), 8, 10),
                ("c5".to_owned(), 0, 5),
                ("c5".to_owned(), 3, 7),
                ("c5".to_owned(), 5, 10),
            ]
        );
    }
}
