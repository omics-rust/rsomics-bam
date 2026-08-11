use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::Path;

use noodles::sam;
use rsomics_common::{Result, RsomicsError};

pub(super) struct ExonIndex {
    by_reference: Vec<Vec<Range<i32>>>,
}

impl ExonIndex {
    pub(super) fn read(path: &Path, header: &sam::Header) -> Result<Self> {
        let references = header
            .reference_sequences()
            .keys()
            .enumerate()
            .map(|(index, name)| (name.to_vec(), index))
            .collect::<HashMap<_, _>>();
        let mut by_reference = vec![Vec::new(); references.len()];
        let file = File::open(path).map_err(RsomicsError::Io)?;
        let mut reader = BufReader::new(file);
        let mut line = Vec::new();
        let mut line_number = 0usize;

        loop {
            line.clear();
            if reader
                .read_until(b'\n', &mut line)
                .map_err(RsomicsError::Io)?
                == 0
            {
                break;
            }
            line_number += 1;
            if line.last() == Some(&b'\n') {
                line.pop();
            }
            if line.last() == Some(&b'\r') {
                line.pop();
            }
            if line.is_empty()
                || line.starts_with(b"#")
                || metadata_line(&line, b"track")
                || metadata_line(&line, b"browser")
            {
                continue;
            }

            let fields = line.splitn(13, |byte| *byte == b'\t').collect::<Vec<_>>();
            if fields.len() != 12 {
                return Err(invalid(path, line_number, "expected exactly 12 fields"));
            }
            let reference = references.get(fields[0]).copied().ok_or_else(|| {
                invalid(
                    path,
                    line_number,
                    format!("unknown reference {}", String::from_utf8_lossy(fields[0])),
                )
            })?;
            let transcript_start = nonnegative(fields[1], path, line_number, "start")?;
            let transcript_end = nonnegative(fields[2], path, line_number, "end")?;
            if transcript_start > transcript_end {
                return Err(invalid(path, line_number, "start exceeds end"));
            }
            let block_count = positive(fields[9], path, line_number, "block count")?;
            let block_count = usize::try_from(block_count)
                .map_err(|_| invalid(path, line_number, "block count exceeds usize"))?;
            let block_sizes = list(fields[10], path, line_number, "block size")?;
            let block_starts = list(fields[11], path, line_number, "block start")?;
            if block_sizes.len() != block_count || block_starts.len() != block_count {
                return Err(invalid(
                    path,
                    line_number,
                    "block count does not match the block lists",
                ));
            }

            for (&size, &start) in block_sizes.iter().zip(&block_starts) {
                if size <= 0 {
                    return Err(invalid(path, line_number, "block size must be positive"));
                }
                if start < 0 {
                    return Err(invalid(path, line_number, "block start cannot be negative"));
                }
                let exon_start = transcript_start
                    .checked_add(start)
                    .ok_or_else(|| invalid(path, line_number, "block start overflows"))?;
                let exon_end = exon_start
                    .checked_add(size)
                    .ok_or_else(|| invalid(path, line_number, "block end overflows"))?;
                if exon_end > transcript_end {
                    return Err(invalid(
                        path,
                        line_number,
                        "block extends beyond the transcript",
                    ));
                }
                by_reference[reference].push(
                    i32::try_from(exon_start)
                        .map_err(|_| invalid(path, line_number, "block start exceeds i32"))?
                        ..i32::try_from(exon_end)
                            .map_err(|_| invalid(path, line_number, "block end exceeds i32"))?,
                );
            }
        }

        for ranges in &mut by_reference {
            ranges.sort_unstable_by_key(|range| (range.start, range.end));
            let mut merged: Vec<Range<i32>> = Vec::with_capacity(ranges.len());
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
        Ok(Self { by_reference })
    }

    pub(super) fn contains(&self, reference: i32, position: i32) -> Result<bool> {
        if reference < 0 || position < 0 {
            return Err(RsomicsError::InvalidInput(
                "gene split record has a negative reference or position".to_owned(),
            ));
        }
        let reference = usize::try_from(reference).map_err(|_| {
            RsomicsError::InvalidInput("gene split reference exceeds usize".to_owned())
        })?;
        let ranges = self.by_reference.get(reference).ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "gene split record reference {reference} is outside the header"
            ))
        })?;
        let index = ranges.partition_point(|range| range.end <= position);
        Ok(ranges
            .get(index)
            .is_some_and(|range| range.start <= position))
    }
}

fn metadata_line(line: &[u8], name: &[u8]) -> bool {
    line.strip_prefix(name)
        .is_some_and(|suffix| suffix.is_empty() || matches!(suffix[0], b' ' | b'\t'))
}

fn nonnegative(field: &[u8], path: &Path, line: usize, name: &str) -> Result<i64> {
    let value = integer(field, path, line, name)?;
    if value < 0 {
        return Err(invalid(path, line, format!("{name} cannot be negative")));
    }
    Ok(value)
}

fn positive(field: &[u8], path: &Path, line: usize, name: &str) -> Result<i64> {
    let value = integer(field, path, line, name)?;
    if value <= 0 {
        return Err(invalid(path, line, format!("{name} must be positive")));
    }
    Ok(value)
}

fn list(field: &[u8], path: &Path, line: usize, name: &str) -> Result<Vec<i64>> {
    let field = field.strip_suffix(b",").unwrap_or(field);
    if field.is_empty() {
        return Ok(Vec::new());
    }
    field
        .split(|byte| *byte == b',')
        .map(|value| integer(value, path, line, name))
        .collect()
}

fn integer(field: &[u8], path: &Path, line: usize, name: &str) -> Result<i64> {
    std::str::from_utf8(field)
        .ok()
        .and_then(|value| value.parse().ok())
        .ok_or_else(|| invalid(path, line, format!("invalid {name}")))
}

fn invalid(path: &Path, line: usize, message: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{}:{line}: {message}", path.display()))
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::num::NonZero;

    use noodles::sam;
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};

    use super::*;

    fn header() -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZero::new(100).unwrap()),
            )
            .add_reference_sequence(
                "Chr1",
                Map::<ReferenceSequence>::new(NonZero::new(100).unwrap()),
            )
            .build()
    }

    #[test]
    fn rows_merge_by_reference_and_keep_half_open_boundaries() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("genes.bed");
        fs::write(
            &path,
            b"# comment\ntrack name=genes\nbrowser position chr1:1-100\nchr1\t10\t40\ta\t0\t+\t10\t40\t0\t2\t5,10,\t0,20,\nchr1\t15\t30\tb\t0\t+\t15\t30\t0\t1\t15,\t0,\nChr1\t50\t60\tc\t0\t+\t50\t60\t0\t1\t10,\t0,\n",
        )
        .unwrap();
        let index = ExonIndex::read(&path, &header()).unwrap();

        for position in [10, 14, 15, 29, 30, 39] {
            assert!(index.contains(0, position).unwrap(), "position {position}");
        }
        for position in [9, 40, 50] {
            assert!(!index.contains(0, position).unwrap(), "position {position}");
        }
        assert!(index.contains(1, 50).unwrap());
        assert!(!index.contains(1, 49).unwrap());
        assert!(!index.contains(1, 60).unwrap());
    }

    #[test]
    fn malformed_rows_fail_loudly() {
        let cases = [
            ("short", "chr1\t0\t10"),
            ("negative", "chr1\t-1\t10\ta\t0\t+\t0\t10\t0\t1\t5,\t0,"),
            ("reversed", "chr1\t10\t9\ta\t0\t+\t9\t10\t0\t1\t1,\t0,"),
            ("zero blocks", "chr1\t0\t10\ta\t0\t+\t0\t10\t0\t0\t\t"),
            (
                "count mismatch",
                "chr1\t0\t10\ta\t0\t+\t0\t10\t0\t2\t5,\t0,",
            ),
            (
                "invalid integer",
                "chr1\tx\t10\ta\t0\t+\t0\t10\t0\t1\t5,\t0,",
            ),
            ("zero block", "chr1\t0\t10\ta\t0\t+\t0\t10\t0\t1\t0,\t0,"),
            (
                "overflow",
                "chr1\t0\t2147483647\ta\t0\t+\t0\t10\t0\t1\t10,\t2147483647,",
            ),
            ("outside", "chr1\t10\t20\ta\t0\t+\t10\t20\t0\t1\t11,\t0,"),
            (
                "unknown reference",
                "chr2\t0\t10\ta\t0\t+\t0\t10\t0\t1\t5,\t0,",
            ),
        ];
        let directory = tempfile::tempdir().unwrap();
        for (name, row) in cases {
            let path = directory.path().join(name);
            fs::write(&path, format!("{row}\n")).unwrap();
            assert!(ExonIndex::read(&path, &header()).is_err(), "{name}");
        }
    }
}
