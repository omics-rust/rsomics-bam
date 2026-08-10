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
use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};

const UNMAPPED: u16 = 0x04;
const NM: [u8; 2] = *b"NM";
const MD: [u8; 2] = *b"MD";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Recalculation {
    Unchanged,
    MissingSequence,
    Updated { corrected_existing: bool },
}

#[derive(Clone, Copy)]
enum EqualBase {
    Literal,
    ReferenceMatch,
}

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
    if record.sequence().is_empty() {
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
    let cigar = record.cigar().as_ref().to_vec();
    let (md, nm) = calculate(
        usize::from(start) - 1,
        &cigar,
        record.sequence_mut().as_mut(),
        reference,
        false,
        EqualBase::Literal,
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

pub(crate) fn recalculate_record(
    header: &sam::Header,
    record: &mut RecordBuf,
    cache: &mut ReferenceCache,
    use_equal: bool,
) -> Result<Recalculation> {
    if record.flags().is_unmapped() {
        return Ok(Recalculation::Unchanged);
    }
    if record.sequence().is_empty() {
        return Ok(Recalculation::MissingSequence);
    }
    let reference_id = record.reference_sequence_id().ok_or_else(|| {
        RsomicsError::InvalidInput("mapped record has no reference sequence".to_owned())
    })?;
    let start = record.alignment_start().ok_or_else(|| {
        RsomicsError::InvalidInput("mapped record has no alignment start".to_owned())
    })?;
    let reference = cache.sequence(header, reference_id)?;
    let cigar = record.cigar().as_ref().to_vec();
    let (md, nm) = calculate(
        usize::from(start) - 1,
        &cigar,
        record.sequence_mut().as_mut(),
        reference,
        use_equal,
        EqualBase::ReferenceMatch,
    )?;
    let corrected_existing = apply_decoded_tags(record, nm, &md)?;
    Ok(Recalculation::Updated { corrected_existing })
}

pub(crate) fn recalculate_raw(
    header: &sam::Header,
    record: &mut RawRecord,
    cache: &mut ReferenceCache,
    use_equal: bool,
    cigar: &mut Vec<(u8, u32)>,
    md: &mut Vec<u8>,
) -> Result<Recalculation> {
    if record.flags() & UNMAPPED != 0 {
        return Ok(Recalculation::Unchanged);
    }
    if record.sequence_len() == 0 {
        return Ok(Recalculation::MissingSequence);
    }
    let reference_id = usize::try_from(record.reference_sequence_id()).map_err(|_| {
        RsomicsError::InvalidInput("mapped BAM record has no reference sequence".to_owned())
    })?;
    let start = usize::try_from(record.alignment_start()).map_err(|_| {
        RsomicsError::InvalidInput("mapped BAM record has no alignment start".to_owned())
    })?;
    let reference = cache.sequence(header, reference_id)?;
    record.decode_cigar_into(cigar)?;
    let sequence_len = record.sequence_len();
    let edit_distance = calculate_raw(
        start,
        cigar,
        record.seq_bytes_mut(),
        sequence_len,
        reference,
        use_equal,
        md,
    )?;
    let corrected_existing = apply_raw_tags(record, edit_distance, md)?;
    Ok(Recalculation::Updated { corrected_existing })
}

fn calculate(
    mut reference_position: usize,
    cigar: &[sam::alignment::record::cigar::Op],
    read: &mut [u8],
    reference: &[u8],
    use_equal: bool,
    equal_base: EqualBase,
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
                    .get_mut(read_position..read_end)
                    .ok_or_else(record_overflow)?;
                let reference_bases = reference
                    .get(reference_position..reference_end)
                    .ok_or_else(record_overflow)?;

                for (read_base, &reference_base) in read_bases.iter_mut().zip(reference_bases) {
                    if bases_match(*read_base, reference_base, equal_base) {
                        if use_equal {
                            *read_base = b'=';
                        }
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

fn bases_match(read: u8, reference: u8, equal_base: EqualBase) -> bool {
    let read = base_code(read);
    let reference = base_code(reference);
    (matches!(equal_base, EqualBase::ReferenceMatch) && read == 0)
        || (read == reference && read != 15)
}

fn calculate_raw(
    mut reference_position: usize,
    cigar: &[(u8, u32)],
    read: &mut [u8],
    read_len: usize,
    reference: &[u8],
    use_equal: bool,
    md: &mut Vec<u8>,
) -> Result<u32> {
    let mut read_position = 0usize;
    let mut matched = 0u64;
    let mut edit_distance = 0u32;
    md.clear();

    for &(kind, length) in cigar {
        let length = usize::try_from(length).map_err(|_| record_overflow())?;
        match kind {
            0 | 7 | 8 => {
                let read_end = read_position
                    .checked_add(length)
                    .filter(|end| *end <= read_len)
                    .ok_or_else(record_overflow)?;
                let reference_end = reference_position
                    .checked_add(length)
                    .filter(|end| *end <= reference.len())
                    .ok_or_else(record_overflow)?;
                for offset in 0..length {
                    let query = read_position + offset;
                    let byte = query / 2;
                    let read_code = if query.is_multiple_of(2) {
                        read[byte] >> 4
                    } else {
                        read[byte] & 0x0f
                    };
                    let reference_base = reference[reference_position + offset];
                    let reference_code = base_code(reference_base);
                    if read_code == 0 || (read_code == reference_code && read_code != 15) {
                        if use_equal {
                            if query.is_multiple_of(2) {
                                read[byte] &= 0x0f;
                            } else {
                                read[byte] &= 0xf0;
                            }
                        }
                        matched += 1;
                    } else {
                        append_number(md, matched);
                        md.push(reference_base.to_ascii_uppercase());
                        matched = 0;
                        edit_distance = edit_distance.checked_add(1).ok_or_else(record_overflow)?;
                    }
                }
                read_position = read_end;
                reference_position = reference_end;
            }
            1 => {
                read_position = read_position
                    .checked_add(length)
                    .filter(|end| *end <= read_len)
                    .ok_or_else(record_overflow)?;
                edit_distance = edit_distance
                    .checked_add(u32::try_from(length).map_err(|_| record_overflow())?)
                    .ok_or_else(record_overflow)?;
            }
            2 => {
                let reference_end = reference_position
                    .checked_add(length)
                    .filter(|end| *end <= reference.len())
                    .ok_or_else(record_overflow)?;
                append_number(md, matched);
                md.push(b'^');
                md.extend(
                    reference[reference_position..reference_end]
                        .iter()
                        .map(u8::to_ascii_uppercase),
                );
                matched = 0;
                reference_position = reference_end;
                edit_distance = edit_distance
                    .checked_add(u32::try_from(length).map_err(|_| record_overflow())?)
                    .ok_or_else(record_overflow)?;
            }
            3 => {
                reference_position = reference_position
                    .checked_add(length)
                    .filter(|end| *end <= reference.len())
                    .ok_or_else(record_overflow)?;
            }
            4 => {
                read_position = read_position
                    .checked_add(length)
                    .filter(|end| *end <= read_len)
                    .ok_or_else(record_overflow)?;
            }
            5 | 6 => {}
            _ => {
                return Err(RsomicsError::InvalidInput(format!(
                    "unsupported BAM CIGAR operation code {kind}"
                )));
            }
        }
    }
    if read_position != read_len {
        return Err(record_overflow());
    }
    append_number(md, matched);
    Ok(edit_distance)
}

fn apply_decoded_tags(record: &mut RecordBuf, nm: u32, md: &str) -> Result<bool> {
    let data = record.data_mut();
    let old_nm = data.get(&Tag::EDIT_DISTANCE);
    let old_md = data.get(&Tag::MISMATCHED_POSITIONS);
    let nm_matches = old_nm.and_then(Value::as_int) == Some(i64::from(nm));
    let md_matches = old_md.is_some_and(|value| match value {
        Value::String(value) => value.eq_ignore_ascii_case(md.as_bytes()),
        _ => false,
    });
    let corrected = old_nm.is_some_and(|_| !nm_matches) || old_md.is_some_and(|_| !md_matches);
    if !nm_matches {
        let nm = i32::try_from(nm).map_err(|_| record_overflow())?;
        replace_decoded_tag(data, Tag::EDIT_DISTANCE, Value::Int32(nm));
    }
    if !md_matches {
        replace_decoded_tag(data, Tag::MISMATCHED_POSITIONS, Value::String(md.into()));
    }
    Ok(corrected)
}

fn replace_decoded_tag(data: &mut sam::alignment::record_buf::Data, tag: Tag, value: Value) {
    let fields = data
        .iter()
        .filter(|(candidate, _)| *candidate != tag)
        .map(|(candidate, value)| (candidate, value.clone()))
        .chain(std::iter::once((tag, value)))
        .collect();
    *data = fields;
}

fn apply_raw_tags(record: &mut RawRecord, edit_distance: u32, md: &mut Vec<u8>) -> Result<bool> {
    let nm_type = record.aux_type(NM);
    let old_edit_distance = nm_type.and_then(|type_code| {
        record
            .aux_value(NM)
            .and_then(|value| raw_edit_distance(type_code, value))
    });
    let nm_matches = old_edit_distance == Some(i64::from(edit_distance));
    let md_type = record.aux_type(MD);
    let md_matches = md_type == Some(b'Z')
        && record
            .aux_value(MD)
            .and_then(|value| value.strip_suffix(&[0]))
            .is_some_and(|value| value.eq_ignore_ascii_case(md));
    let corrected = nm_type.is_some() && !nm_matches || md_type.is_some() && !md_matches;
    if !nm_matches {
        let value = i32::try_from(edit_distance).map_err(|_| record_overflow())?;
        if nm_type.is_some() {
            record.set_aux(NM, b'i', &value.to_le_bytes())?;
        } else {
            record.append_aux(NM, b'i', &value.to_le_bytes())?;
        }
    }
    if !md_matches {
        md.push(0);
        let result = if md_type.is_some() {
            record.set_aux(MD, b'Z', md)
        } else {
            record.append_aux(MD, b'Z', md)
        };
        md.pop();
        result?;
    }
    Ok(corrected)
}

fn raw_edit_distance(type_code: u8, value: &[u8]) -> Option<i64> {
    match type_code {
        b'c' => value
            .first()
            .map(|value| i64::from(i8::from_le_bytes([*value]))),
        b'C' => value.first().map(|value| i64::from(*value)),
        b's' => value
            .get(..2)
            .map(|value| i64::from(i16::from_le_bytes([value[0], value[1]]))),
        b'S' => value
            .get(..2)
            .map(|value| i64::from(u16::from_le_bytes([value[0], value[1]]))),
        b'i' => value
            .get(..4)
            .map(|value| i64::from(i32::from_le_bytes([value[0], value[1], value[2], value[3]]))),
        b'I' => value
            .get(..4)
            .map(|value| i64::from(u32::from_le_bytes([value[0], value[1], value[2], value[3]]))),
        _ => None,
    }
}

fn append_number(output: &mut Vec<u8>, value: u64) {
    if value == 0 {
        output.push(b'0');
        return;
    }
    let mut digits = [0; 20];
    let mut start = digits.len();
    let mut value = value;
    while value > 0 {
        start -= 1;
        digits[start] = b'0' + (value % 10) as u8;
        value /= 10;
    }
    output.extend_from_slice(&digits[start..]);
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
    use rsomics_bamio::raw::RawRecordEncoder;

    #[test]
    fn md_and_nm_cover_mismatches_indels_and_skips() {
        let cigar = [
            Op::new(Kind::Match, 4),
            Op::new(Kind::Insertion, 1),
            Op::new(Kind::Deletion, 2),
            Op::new(Kind::Skip, 1),
            Op::new(Kind::Match, 2),
        ];
        let mut read = b"ACTGACC".to_vec();
        let (md, nm) = calculate(
            0,
            &cigar,
            &mut read,
            b"ACGTTAACC",
            false,
            EqualBase::ReferenceMatch,
        )
        .unwrap();
        assert_eq!(md, "2G0T0^TA2");
        assert_eq!(nm, 5);
    }

    #[test]
    fn equal_bases_match_any_reference_base() {
        let cigar = [Op::new(Kind::Match, 4)];
        let mut read = b"====".to_vec();
        let (md, nm) = calculate(
            0,
            &cigar,
            &mut read,
            b"ACGT",
            false,
            EqualBase::ReferenceMatch,
        )
        .unwrap();
        assert_eq!(md, "4");
        assert_eq!(nm, 0);
    }

    #[test]
    fn equal_mode_changes_only_reference_matches() {
        let cigar = [Op::new(Kind::Match, 4)];
        let mut read = b"ATGT".to_vec();
        let (md, nm) = calculate(
            0,
            &cigar,
            &mut read,
            b"ACGT",
            true,
            EqualBase::ReferenceMatch,
        )
        .unwrap();
        assert_eq!(read, b"=T==");
        assert_eq!(md, "1C2");
        assert_eq!(nm, 1);
    }

    #[test]
    fn inconsistent_cigar_is_rejected() {
        let cigar = [Op::new(Kind::Match, 5)];
        let mut read = b"ACGT".to_vec();
        assert!(
            calculate(
                0,
                &cigar,
                &mut read,
                b"ACGT",
                false,
                EqualBase::ReferenceMatch,
            )
            .is_err()
        );
    }

    #[test]
    fn cram_completion_treats_equal_as_a_literal_query_base() {
        let cigar = [Op::new(Kind::Match, 4)];
        let mut read = b"====".to_vec();
        let (md, nm) = calculate(0, &cigar, &mut read, b"ACGT", false, EqualBase::Literal).unwrap();
        assert_eq!(md, "0A0C0G0T0");
        assert_eq!(nm, 4);
    }

    #[test]
    fn raw_tags_preserve_a_correct_numeric_subtype_and_replace_a_wrong_value() {
        let input = b"@SQ\tSN:chr1\tLN:4\nread\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        let mut reader = sam::io::Reader::new(&input[..]);
        let header = reader.read_header().unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let mut record = RawRecordEncoder::new().encode(&header, &record).unwrap();
        record.set_aux(NM, b'C', &[0]).unwrap();
        record.set_aux(MD, b'Z', b"4\0").unwrap();

        assert!(!apply_raw_tags(&mut record, 0, &mut b"4".to_vec()).unwrap());
        assert_eq!(record.aux_type(NM), Some(b'C'));

        assert!(apply_raw_tags(&mut record, 1, &mut b"3A0".to_vec()).unwrap());
        assert_eq!(record.aux_type(NM), Some(b'i'));
        assert_eq!(
            raw_edit_distance(record.aux_type(NM).unwrap(), record.aux_value(NM).unwrap()),
            Some(1)
        );
        assert_eq!(record.aux_value(MD), Some(b"3A0\0".as_slice()));
    }
}
