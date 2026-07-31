use std::fs::File;
use std::path::{Path, PathBuf};

use noodles::core::Region;
use noodles::fasta;
use noodles::sam::{
    self,
    alignment::{
        Record, RecordBuf,
        record::MappingQuality,
        record::{cigar::op::Kind, data::field::Tag},
        record_buf::data::field::Value,
    },
};
use rsomics_common::{Result, RsomicsError};

pub(crate) struct ReferenceCache {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    path: PathBuf,
    reference_id: Option<usize>,
    sequence: Vec<u8>,
}

impl ReferenceCache {
    pub(crate) fn open(path: &Path) -> Result<Self> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .map_err(|error| {
                RsomicsError::ConfigError(format!(
                    "opening indexed reference {}: {error}",
                    path.display()
                ))
            })?;
        Ok(Self {
            reader,
            path: path.to_path_buf(),
            reference_id: None,
            sequence: Vec::new(),
        })
    }

    fn sequence<'a>(&'a mut self, header: &sam::Header, reference_id: usize) -> Result<&'a [u8]> {
        if self.reference_id != Some(reference_id) {
            let (name, _) = header
                .reference_sequences()
                .get_index(reference_id)
                .ok_or_else(|| {
                    RsomicsError::InvalidInput(format!(
                        "record references unknown sequence ID {reference_id}"
                    ))
                })?;
            let record = self
                .reader
                .query(&Region::new(name.clone(), ..))
                .map_err(|error| {
                    RsomicsError::InvalidInput(format!(
                        "reading {} from {}: {error}",
                        String::from_utf8_lossy(name),
                        self.path.display()
                    ))
                })?;
            self.sequence.clear();
            self.sequence.extend_from_slice(record.sequence().as_ref());
            self.reference_id = Some(reference_id);
        }
        Ok(&self.sequence)
    }
}

pub(crate) fn complete(
    header: &sam::Header,
    record: &dyn Record,
    cache: Option<&mut ReferenceCache>,
) -> Result<RecordBuf> {
    let mut record =
        RecordBuf::try_from_alignment_record(header, record).map_err(RsomicsError::Io)?;
    if record.flags().is_unmapped() {
        if record.mapping_quality().is_none() {
            *record.mapping_quality_mut() = Some(MappingQuality::MIN);
        }
        return Ok(record);
    }
    if record.data().get(&Tag::MISMATCHED_POSITIONS).is_some()
        && record.data().get(&Tag::EDIT_DISTANCE).is_some()
    {
        return Ok(record);
    }

    let cache = cache.ok_or_else(|| {
        RsomicsError::ConfigError(
            "a reference FASTA is required to restore CRAM MD and NM tags".to_owned(),
        )
    })?;
    let reference_id = record.reference_sequence_id().ok_or_else(|| {
        RsomicsError::InvalidInput("mapped CRAM record has no reference sequence".to_owned())
    })?;
    let start = record.alignment_start().ok_or_else(|| {
        RsomicsError::InvalidInput("mapped CRAM record has no alignment start".to_owned())
    })?;
    let reference = cache.sequence(header, reference_id)?;
    let (md, nm) = calculate(
        usize::from(start) - 1,
        record.cigar().as_ref(),
        record.sequence().as_ref(),
        reference,
    )?;

    let fields = record
        .data()
        .iter()
        .filter(|(tag, _)| !matches!(*tag, Tag::MISMATCHED_POSITIONS | Tag::EDIT_DISTANCE))
        .map(|(tag, value)| (tag, value.clone()))
        .collect::<Vec<_>>();
    let data = record.data_mut();
    data.clear();
    data.insert(Tag::MISMATCHED_POSITIONS, Value::String(md.into()));
    data.insert(Tag::EDIT_DISTANCE, Value::from(nm));
    for (tag, value) in fields {
        data.insert(tag, value);
    }

    Ok(record)
}

fn calculate(
    mut reference_position: usize,
    cigar: &[sam::alignment::record::cigar::Op],
    read: &[u8],
    reference: &[u8],
) -> Result<(String, u32)> {
    let mut read_position = 0usize;
    let mut matched = 0;
    let mut edit_distance = 0u32;
    let mut md = String::new();

    for operation in cigar {
        let length = operation.len();
        match operation.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let read_end = read_position
                    .checked_add(length)
                    .ok_or_else(record_overflow)?;
                let reference_end = reference_position
                    .checked_add(length)
                    .ok_or_else(record_overflow)?;
                let read_bases = read
                    .get(read_position..read_end)
                    .ok_or_else(record_overflow)?;
                let reference_bases = reference
                    .get(reference_position..reference_end)
                    .ok_or_else(record_overflow)?;

                for (&read_base, &reference_base) in read_bases.iter().zip(reference_bases) {
                    if bases_match(read_base, reference_base) {
                        matched += 1;
                    } else {
                        md.push_str(&matched.to_string());
                        md.push(char::from(reference_base.to_ascii_uppercase()));
                        matched = 0;
                        edit_distance = edit_distance.checked_add(1).ok_or_else(record_overflow)?;
                    }
                }
                read_position = read_end;
                reference_position = reference_end;
            }
            Kind::Deletion => {
                let reference_end = reference_position
                    .checked_add(length)
                    .ok_or_else(record_overflow)?;
                let deleted = reference
                    .get(reference_position..reference_end)
                    .ok_or_else(record_overflow)?;
                md.push_str(&matched.to_string());
                md.push('^');
                for base in deleted {
                    md.push(char::from(base.to_ascii_uppercase()));
                }
                matched = 0;
                reference_position = reference_end;
                edit_distance = edit_distance
                    .checked_add(u32::try_from(length).map_err(|_| record_overflow())?)
                    .ok_or_else(record_overflow)?;
            }
            Kind::Insertion => {
                read_position = read_position
                    .checked_add(length)
                    .filter(|end| *end <= read.len())
                    .ok_or_else(record_overflow)?;
                edit_distance = edit_distance
                    .checked_add(u32::try_from(length).map_err(|_| record_overflow())?)
                    .ok_or_else(record_overflow)?;
            }
            Kind::SoftClip => {
                read_position = read_position
                    .checked_add(length)
                    .filter(|end| *end <= read.len())
                    .ok_or_else(record_overflow)?;
            }
            Kind::Skip => {
                reference_position = reference_position
                    .checked_add(length)
                    .filter(|end| *end <= reference.len())
                    .ok_or_else(record_overflow)?;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    if read_position != read.len() {
        return Err(record_overflow());
    }
    md.push_str(&matched.to_string());
    Ok((md, edit_distance))
}

fn bases_match(read: u8, reference: u8) -> bool {
    let read = base_code(read);
    let reference = base_code(reference);
    read == reference && read != 15
}

fn base_code(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'=' => 0,
        b'A' => 1,
        b'C' => 2,
        b'M' => 3,
        b'G' => 4,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'T' => 8,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        _ => 15,
    }
}

fn record_overflow() -> RsomicsError {
    RsomicsError::InvalidInput(
        "alignment CIGAR, sequence, and reference lengths are inconsistent".to_owned(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record::cigar::{Op, op::Kind};

    #[test]
    fn md_and_nm_cover_mismatches_indels_and_skips() {
        let cigar = [
            Op::new(Kind::Match, 4),
            Op::new(Kind::Insertion, 1),
            Op::new(Kind::Deletion, 2),
            Op::new(Kind::Skip, 1),
            Op::new(Kind::Match, 2),
        ];
        let (md, nm) = calculate(0, &cigar, b"ACTGACC", b"ACGTTAACC").unwrap();
        assert_eq!(md, "2G0T0^TA2");
        assert_eq!(nm, 5);
    }

    #[test]
    fn equal_bases_are_mismatches() {
        let cigar = [Op::new(Kind::Match, 4)];
        let (md, nm) = calculate(0, &cigar, b"====", b"ACGT").unwrap();
        assert_eq!(md, "0A0C0G0T0");
        assert_eq!(nm, 4);
    }
}
