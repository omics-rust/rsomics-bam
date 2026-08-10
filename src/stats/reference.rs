use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use noodles::core::{Position, Region};
use noodles::fasta;
use rsomics_common::{Result, RsomicsError};

pub(crate) struct Reference {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    lengths: HashMap<Vec<u8>, usize>,
    cache: Cache,
    chunk_size: usize,
}

#[derive(Default)]
struct Cache {
    name: Vec<u8>,
    start: usize,
    sequence: Vec<u8>,
}

#[derive(Clone, Copy)]
pub(crate) struct Slice<'a> {
    pub(crate) start: usize,
    pub(crate) sequence: &'a [u8],
}

impl Reference {
    pub(crate) fn open(path: &Path, chunk_size: usize) -> Result<Self> {
        let index = load_index(path)?;
        let lengths = index
            .as_ref()
            .iter()
            .map(|record| {
                usize::try_from(record.length())
                    .map(|length| (record.name().to_vec(), length))
                    .map_err(|_| {
                        RsomicsError::InvalidInput(format!(
                            "reference sequence {} is too long",
                            String::from_utf8_lossy(record.name())
                        ))
                    })
            })
            .collect::<Result<HashMap<_, _>>>()?;
        let reader = fasta::io::indexed_reader::Builder::default()
            .set_index(index)
            .build_from_path(path)
            .map_err(|error| reference_io(path, "opening", error))?;
        Ok(Self {
            reader,
            lengths,
            cache: Cache::default(),
            chunk_size: chunk_size.max(1),
        })
    }

    pub(crate) fn contains(&self, name: &[u8]) -> bool {
        self.lengths.contains_key(name)
    }

    pub(crate) fn length(&self, name: &[u8]) -> Option<usize> {
        self.lengths.get(name).copied()
    }

    pub(crate) fn get(&mut self, name: &[u8], start: usize, end: usize) -> Result<Slice<'_>> {
        let length = self.lengths.get(name).copied().ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "reference FASTA is missing sequence {}",
                String::from_utf8_lossy(name)
            ))
        })?;
        if start > end || end > length {
            return Err(RsomicsError::InvalidInput(format!(
                "reference FASTA sequence {} is shorter than the requested range",
                String::from_utf8_lossy(name)
            )));
        }
        let cached_end = self.cache.start + self.cache.sequence.len();
        if self.cache.name != name || start < self.cache.start || end > cached_end {
            let chunk_end = start.saturating_add(self.chunk_size).max(end).min(length);
            self.load(name, start, chunk_end)?;
        }
        let local_start = start - self.cache.start;
        let local_end = end - self.cache.start;
        Ok(Slice {
            start,
            sequence: &self.cache.sequence[local_start..local_end],
        })
    }

    pub(crate) fn count_bases(
        &mut self,
        name: &[u8],
        start: usize,
        end: usize,
    ) -> Result<(u64, u64, i64)> {
        let mut gc = 0u64;
        let mut at = 0u64;
        let mut unknown = 0i64;
        let mut position = start;
        while position < end {
            let chunk_end = position.saturating_add(self.chunk_size).min(end);
            let slice = self.get(name, position, chunk_end)?.sequence;
            for base in slice {
                match base.to_ascii_uppercase() {
                    b'G' | b'C' => gc += 1,
                    b'A' | b'T' => at += 1,
                    b'N' => unknown += 1,
                    _ => {}
                }
            }
            position = chunk_end;
        }
        Ok((gc, at, unknown))
    }

    fn load(&mut self, name: &[u8], start: usize, end: usize) -> Result<()> {
        if start == end {
            self.cache.name.clear();
            self.cache.name.extend_from_slice(name);
            self.cache.start = start;
            self.cache.sequence.clear();
            return Ok(());
        }
        let start_position = start
            .checked_add(1)
            .and_then(|value| Position::try_from(value).ok())
            .ok_or_else(|| {
                RsomicsError::InvalidInput("reference range start overflows".to_owned())
            })?;
        let end_position = Position::try_from(end)
            .map_err(|_| RsomicsError::InvalidInput("reference range end overflows".to_owned()))?;
        let region = Region::new(name, start_position..=end_position);
        let record = self.reader.query(&region).map_err(|error| {
            RsomicsError::InvalidInput(format!("reading reference range {region}: {error}"))
        })?;
        self.cache.name.clear();
        self.cache.name.extend_from_slice(name);
        self.cache.start = start;
        self.cache.sequence.clear();
        self.cache.sequence.extend(
            record
                .sequence()
                .as_ref()
                .iter()
                .map(u8::to_ascii_uppercase),
        );
        Ok(())
    }
}

fn load_index(path: &Path) -> Result<fasta::fai::Index> {
    let index_path = index_path(path);
    if index_path.exists() {
        return fasta::fai::fs::read(&index_path)
            .map_err(|error| reference_io(&index_path, "reading index for", error));
    }
    let file = File::open(path).map_err(|error| reference_io(path, "opening", error))?;
    let mut indexer = fasta::io::Indexer::new(BufReader::new(file));
    let mut records = Vec::new();
    while let Some(record) = indexer.index_record().map_err(|error| {
        RsomicsError::InvalidInput(format!("indexing reference {}: {error}", path.display()))
    })? {
        records.push(record);
    }
    if records.is_empty() {
        return Err(RsomicsError::InvalidInput(format!(
            "reference {} contains no sequences",
            path.display()
        )));
    }
    Ok(records.into())
}

fn index_path(path: &Path) -> PathBuf {
    let mut value = OsString::from(path);
    value.push(".fai");
    value.into()
}

fn reference_io(path: &Path, action: &str, error: std::io::Error) -> RsomicsError {
    RsomicsError::Io(std::io::Error::new(
        error.kind(),
        format!("{action} reference {}: {error}", path.display()),
    ))
}
