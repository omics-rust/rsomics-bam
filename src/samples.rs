use std::io::Write;
use std::path::{Path, PathBuf};

use rsomics_common::{Result, RsomicsError};
use rust_htslib::bam::Read;
use serde::Serialize;

use crate::{hts_metadata, input};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Tag([u8; 2]);

impl Tag {
    pub fn parse(value: &str) -> Result<Self> {
        value.as_bytes().try_into().map(Self).map_err(|_| {
            RsomicsError::ConfigError(format!(
                "read-group tag must contain exactly two bytes: {value:?}"
            ))
        })
    }

    #[must_use]
    pub fn as_str(self) -> String {
        String::from_utf8_lossy(&self.0).into_owned()
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Input {
    pub path: PathBuf,
    pub index: Option<PathBuf>,
}

#[derive(Clone, Debug)]
pub struct Options<'a> {
    pub tag: Tag,
    pub test_index: bool,
    pub references: &'a [PathBuf],
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Entry {
    pub value: String,
    pub path: PathBuf,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub indexed: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reference: Option<PathBuf>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Report {
    pub tag: String,
    pub index_column: bool,
    pub reference_column: bool,
    pub entries: Vec<Entry>,
}

impl Report {
    pub fn write(&self, header: bool, mut output: impl Write) -> Result<()> {
        if header {
            write!(output, "#{}\tPATH", self.tag).map_err(RsomicsError::Io)?;
            if self.index_column {
                write!(output, "\tINDEX").map_err(RsomicsError::Io)?;
            }
            if self.reference_column {
                write!(output, "\tREFERENCE").map_err(RsomicsError::Io)?;
            }
            writeln!(output).map_err(RsomicsError::Io)?;
        }

        for entry in &self.entries {
            write!(output, "{}\t{}", entry.value, entry.path.display())
                .map_err(RsomicsError::Io)?;
            if let Some(indexed) = entry.indexed {
                write!(output, "\t{}", if indexed { "Y" } else { "N" })
                    .map_err(RsomicsError::Io)?;
            }
            if self.reference_column {
                write!(output, "\t").map_err(RsomicsError::Io)?;
                if let Some(reference) = &entry.reference {
                    write!(output, "{}", reference.display()).map_err(RsomicsError::Io)?;
                } else {
                    write!(output, ".").map_err(RsomicsError::Io)?;
                }
            }
            writeln!(output).map_err(RsomicsError::Io)?;
        }

        output.flush().map_err(RsomicsError::Io)
    }
}

pub fn collect(inputs: &[Input], options: Options<'_>) -> Result<Report> {
    let references = options
        .references
        .iter()
        .map(|path| hts_metadata::load_reference(path))
        .collect::<Result<Vec<_>>>()?;
    let mut entries = Vec::new();

    for input_spec in inputs {
        let reader = input::open(&input_spec.path, None, 0)?;
        let indexed = options
            .test_index
            .then(|| {
                hts_metadata::has_index(&reader, &input_spec.path, input_spec.index.as_deref())
            })
            .transpose()?;
        let reference = matching_reference(reader.header(), &references);
        let mut values = sample_values(reader.header().as_bytes(), options.tag, &input_spec.path)?;
        if values.is_empty() {
            values.push(".".to_owned());
        }
        entries.extend(values.into_iter().map(|value| Entry {
            value,
            path: input_spec.path.clone(),
            indexed,
            reference: reference.cloned(),
        }));
    }

    Ok(Report {
        tag: options.tag.as_str(),
        index_column: options.test_index,
        reference_column: !references.is_empty(),
        entries,
    })
}

fn sample_values(header: &[u8], tag: Tag, path: &Path) -> Result<Vec<String>> {
    let mut values = SamtoolsStringSet::new();
    for line in header.split(|byte| *byte == b'\n') {
        let Some(fields) = line.strip_prefix(b"@RG\t") else {
            continue;
        };
        for field in fields.split(|byte| *byte == b'\t') {
            if field.len() >= 3 && field[..2] == tag.0 && field[2] == b':' {
                let value = std::str::from_utf8(&field[3..]).map_err(|error| {
                    RsomicsError::InvalidInput(format!(
                        "reading tag {} from {}: {error}",
                        tag.as_str(),
                        path.display()
                    ))
                })?;
                values.insert(value.to_owned());
                break;
            }
        }
    }
    Ok(values.into_values())
}

fn matching_reference<'a>(
    header: &rust_htslib::bam::HeaderView,
    references: &'a [hts_metadata::ReferenceDictionary],
) -> Option<&'a PathBuf> {
    references.iter().rev().find_map(|reference| {
        if reference.targets.len() != header.target_count() as usize {
            return None;
        }
        reference
            .targets
            .iter()
            .enumerate()
            .all(|(position, (name, length))| {
                header.tid2name(position as u32) == name
                    && header.target_len(position as u32) == Some(*length)
            })
            .then_some(&reference.path)
    })
}

struct SamtoolsStringSet {
    buckets: Vec<Option<String>>,
    len: usize,
    upper_bound: usize,
}

impl SamtoolsStringSet {
    fn new() -> Self {
        Self {
            buckets: Vec::new(),
            len: 0,
            upper_bound: 0,
        }
    }

    fn insert(&mut self, value: String) {
        if self.len >= self.upper_bound {
            self.resize(self.buckets.len() + 1);
        }
        let slot = probe(&self.buckets, &value);
        if self.buckets[slot].is_none() {
            self.buckets[slot] = Some(value);
            self.len += 1;
        }
    }

    fn resize(&mut self, requested: usize) {
        let size = requested.next_power_of_two().max(4);
        let mut old = std::mem::take(&mut self.buckets);
        let mut resized = std::iter::repeat_with(|| None)
            .take(size)
            .collect::<Vec<Option<String>>>();

        for position in 0..old.len() {
            let Some(mut value) = old[position].take() else {
                continue;
            };
            loop {
                let slot = probe(&resized, &value);
                if slot < old.len()
                    && let Some(displaced) = old[slot].take()
                {
                    resized[slot] = Some(value);
                    value = displaced;
                } else {
                    resized[slot] = Some(value);
                    break;
                }
            }
        }

        self.buckets = resized;
        self.upper_bound = (size as f64 * 0.77 + 0.5) as usize;
    }

    fn into_values(self) -> Vec<String> {
        self.buckets.into_iter().flatten().collect()
    }
}

fn probe(buckets: &[Option<String>], value: &str) -> usize {
    let mask = buckets.len() - 1;
    let mut position = fnv1a(value.as_bytes()) as usize & mask;
    let mut step = 0;
    while buckets[position]
        .as_ref()
        .is_some_and(|existing| existing != value)
    {
        step += 1;
        position = (position + step) & mask;
    }
    position
}

fn fnv1a(value: &[u8]) -> u32 {
    value.iter().fold(2_166_136_261, |hash, byte| {
        (hash ^ u32::from(*byte)).wrapping_mul(16_777_619)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tags_are_exactly_two_bytes() {
        assert_eq!(Tag::parse("SM").unwrap(), Tag(*b"SM"));
        assert!(Tag::parse("S").is_err());
        assert!(Tag::parse("SAM").is_err());
        assert!(Tag::parse("样").is_err());
    }

    #[test]
    fn sample_order_matches_samtools_hash_iteration() {
        let header = b"@RG\tID:r1\tSM:zeta\tLB:L1\n\
                       @RG\tID:r2\tSM:alpha\tLB:L2\n\
                       @RG\tID:r3\tSM:zeta\tLB:L3\n\
                       @RG\tID:r4\tLB:L4\n\
                       @RG\tID:r5\tSM:beta\tLB:L5\n";
        assert_eq!(
            sample_values(header, Tag(*b"SM"), Path::new("input")).unwrap(),
            ["alpha", "beta", "zeta"]
        );
        assert_eq!(
            sample_values(header, Tag(*b"LB"), Path::new("input")).unwrap(),
            ["L3", "L5", "L2", "L4", "L1"]
        );
    }
}
