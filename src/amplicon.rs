use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Strand {
    Forward,
    Reverse,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Primer {
    pub(crate) start: i64,
    pub(crate) end: i64,
    pub(crate) name: String,
    pub(crate) score: String,
    pub(crate) strand: Option<Strand>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct ReferencePrimers {
    pub(crate) name: String,
    pub(crate) primers: Vec<Primer>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct PrimerBed {
    references: Vec<ReferencePrimers>,
    by_name: HashMap<String, usize>,
}

impl PrimerBed {
    pub(crate) fn read(path: &Path) -> Result<Self> {
        let file = File::open(path).map_err(RsomicsError::Io)?;
        let mut references: Vec<ReferencePrimers> = Vec::new();
        let mut by_name = HashMap::new();

        for (line_index, result) in BufReader::new(file).lines().enumerate() {
            let line = result.map_err(RsomicsError::Io)?;
            let line = line.trim();
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track ")
                || line.starts_with("browser ")
            {
                continue;
            }
            let columns: Vec<_> = line.split_whitespace().collect();
            if columns.len() < 3 {
                return Err(invalid_bed(
                    path,
                    line_index,
                    "expected at least three columns",
                ));
            }
            let start = parse_coordinate(path, line_index, columns[1], "start")?;
            let end = parse_coordinate(path, line_index, columns[2], "end")?;
            if start >= end {
                return Err(invalid_bed(
                    path,
                    line_index,
                    "start must be smaller than end",
                ));
            }
            let strand = match columns.get(5).copied() {
                Some("+") => Some(Strand::Forward),
                Some("-") => Some(Strand::Reverse),
                Some(value) => {
                    return Err(invalid_bed(
                        path,
                        line_index,
                        format!("invalid strand {value:?}"),
                    ));
                }
                None => None,
            };
            let reference_index = match by_name.get(columns[0]).copied() {
                Some(index) => index,
                None => {
                    let index = references.len();
                    by_name.insert(columns[0].to_owned(), index);
                    references.push(ReferencePrimers {
                        name: columns[0].to_owned(),
                        primers: Vec::new(),
                    });
                    index
                }
            };
            references[reference_index].primers.push(Primer {
                start,
                end,
                name: columns.get(3).copied().unwrap_or(".").to_owned(),
                score: columns.get(4).copied().unwrap_or("0").to_owned(),
                strand,
            });
        }

        if references.is_empty() {
            return Err(RsomicsError::InvalidInput(format!(
                "primer BED is empty: {}",
                path.display()
            )));
        }
        Ok(Self {
            references,
            by_name,
        })
    }

    pub(crate) fn references(&self) -> &[ReferencePrimers] {
        &self.references
    }

    pub(crate) fn get(&self, name: &str) -> Option<&ReferencePrimers> {
        self.by_name.get(name).map(|&index| &self.references[index])
    }
}

fn parse_coordinate(path: &Path, line_index: usize, value: &str, label: &str) -> Result<i64> {
    let coordinate = value
        .parse::<i64>()
        .map_err(|_| invalid_bed(path, line_index, format!("invalid {label} coordinate")))?;
    if coordinate < 0 {
        return Err(invalid_bed(
            path,
            line_index,
            format!("{label} coordinate must not be negative"),
        ));
    }
    Ok(coordinate)
}

fn invalid_bed(path: &Path, line_index: usize, message: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{}:{}: {message}", path.display(), line_index + 1))
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;

    #[test]
    fn preserves_reference_and_primer_order() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("primers.bed");
        fs::write(
            &path,
            b"track name=x\nchr2\t20\t30\tb\t0\t-\nchr1\t1\t4\ta\t0\t+\nchr2\t10\t15\tc\t0\t+\n",
        )
        .unwrap();

        let bed = PrimerBed::read(&path).unwrap();
        assert_eq!(bed.references()[0].name, "chr2");
        assert_eq!(bed.references()[1].name, "chr1");
        assert_eq!(bed.get("chr2").unwrap().primers[0].start, 20);
        assert_eq!(bed.get("chr2").unwrap().primers[1].start, 10);
    }

    #[test]
    fn rejects_invalid_coordinates_and_strands() {
        let directory = tempfile::tempdir().unwrap();
        for (name, row) in [
            ("negative", "chr1 -1 5"),
            ("empty", "chr1 5 5"),
            ("reverse", "chr1 6 5"),
            ("strand", "chr1 1 5 name 0 x"),
        ] {
            let path = directory.path().join(name);
            fs::write(&path, row).unwrap();
            let error = PrimerBed::read(&path).unwrap_err().to_string();
            assert!(error.contains(":1:"), "{error}");
        }
    }
}
