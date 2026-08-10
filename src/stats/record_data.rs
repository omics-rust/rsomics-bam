use noodles::sam;
use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use rsomics_bamio::raw::RecordRef;
use rsomics_common::{Result, RsomicsError};
use rust_htslib::bam::record::{Aux, Cigar};

#[derive(Default)]
pub(crate) struct RecordData {
    pub(crate) name: Vec<u8>,
    pub(crate) flags: u16,
    pub(crate) reference: i32,
    pub(crate) position: i64,
    pub(crate) mapping_quality: u8,
    pub(crate) mate_reference: i32,
    pub(crate) mate_position: i64,
    pub(crate) template_length: i64,
    pub(crate) cigar: Vec<(u8, u32)>,
    pub(crate) sequence: Vec<u8>,
    pub(crate) packed_sequence: Vec<u8>,
    pub(crate) qualities: Vec<u8>,
    pub(crate) edit_distance: Option<u64>,
    pub(crate) barcodes: [Option<Vec<u8>>; 4],
    pub(crate) barcode_qualities: [Option<Vec<u8>>; 4],
    pub(crate) read_group: Option<Vec<u8>>,
    pub(crate) split_value: Option<Vec<u8>>,
}

impl RecordData {
    pub(crate) fn decode(
        header: &sam::Header,
        source: &dyn Record,
        split_tag: Option<[u8; 2]>,
    ) -> Result<Self> {
        let record = sam::alignment::RecordBuf::try_from_alignment_record(header, source)
            .map_err(RsomicsError::Io)?;
        let reference = option_index(record.reference_sequence_id())?;
        let mate_reference = option_index(record.mate_reference_sequence_id())?;
        let position = record
            .alignment_start()
            .map(|position| i64::try_from(usize::from(position) - 1))
            .transpose()
            .map_err(|_| RsomicsError::InvalidInput("alignment position overflows".to_owned()))?
            .unwrap_or(-1);
        let mate_position = record
            .mate_alignment_start()
            .map(|position| i64::try_from(usize::from(position) - 1))
            .transpose()
            .map_err(|_| RsomicsError::InvalidInput("mate position overflows".to_owned()))?
            .unwrap_or(-1);
        let cigar = record
            .cigar()
            .as_ref()
            .iter()
            .map(|&operation| (cigar_code(operation.kind()), operation.len() as u32))
            .collect();
        let sequence: Vec<_> = record
            .sequence()
            .as_ref()
            .iter()
            .copied()
            .map(sequence_code)
            .collect();
        let packed_sequence = pack_sequence(&sequence);
        let edit_distance = record
            .data()
            .get(&Tag::EDIT_DISTANCE)
            .map(integer)
            .transpose()?;
        let barcodes = [
            string_field(record.data().get(&Tag::SAMPLE_BARCODE_SEQUENCE), b"BC")?,
            string_field(record.data().get(&Tag::CELL_BARCODE_SEQUENCE), b"CR")?,
            string_field(
                record.data().get(&Tag::ORIGINAL_UMI_BARCODE_SEQUENCE),
                b"OX",
            )?,
            string_field(record.data().get(&Tag::UMI_SEQUENCE), b"RX")?,
        ];
        let barcode_qualities = [
            string_field(
                record.data().get(&Tag::SAMPLE_BARCODE_QUALITY_SCORES),
                b"QT",
            )?,
            string_field(record.data().get(&Tag::CELL_BARCODE_QUALITY_SCORES), b"CY")?,
            string_field(record.data().get(&Tag::ORIGINAL_UMI_QUALITY_SCORES), b"BZ")?,
            string_field(record.data().get(&Tag::UMI_QUALITY_SCORES), b"QX")?,
        ];
        let read_group = string_field(record.data().get(&Tag::READ_GROUP), b"RG")?;
        let split_value = split_tag
            .map(|tag| {
                string_field(record.data().get(&Tag::from(tag)), &tag)?.ok_or_else(|| {
                    RsomicsError::InvalidInput(format!(
                        "record {} is missing split tag {}",
                        String::from_utf8_lossy(record.name().map_or(b"*", |name| name.as_ref())),
                        String::from_utf8_lossy(&tag)
                    ))
                })
            })
            .transpose()?;
        Ok(Self {
            name: record.name().map_or_else(Vec::new, |name| name.to_vec()),
            flags: u16::from(record.flags()),
            reference,
            position,
            mapping_quality: record.mapping_quality().map_or(255, u8::from),
            mate_reference,
            mate_position,
            template_length: i64::from(record.template_length()),
            cigar,
            sequence,
            packed_sequence,
            qualities: record.quality_scores().as_ref().to_vec(),
            edit_distance,
            barcodes,
            barcode_qualities,
            read_group,
            split_value,
        })
    }

    pub(crate) fn decode_raw(
        &mut self,
        source: &RecordRef<'_>,
        split_tag: Option<[u8; 2]>,
    ) -> Result<()> {
        self.name.clear();
        self.name.extend_from_slice(source.name());
        self.flags = source.flags();
        self.reference = source.reference_sequence_id();
        self.position = i64::from(source.alignment_start());
        self.mapping_quality = source.mapping_quality();
        self.mate_reference = source.mate_reference_sequence_id();
        self.mate_position = i64::from(source.mate_alignment_start());
        self.template_length = i64::from(source.template_length());
        source.decode_cigar_into(&mut self.cigar)?;
        self.sequence.clear();
        self.sequence
            .extend((0..source.sequence_len()).map(|index| source.seq_nibble(index)));
        self.packed_sequence.clear();
        self.packed_sequence
            .extend_from_slice(source.seq_bytes_packed());
        self.qualities.clear();
        self.qualities.extend_from_slice(source.quality_scores());
        self.edit_distance = raw_integer(source, *b"NM")?;
        for (index, tag) in [*b"BC", *b"CR", *b"OX", *b"RX"].into_iter().enumerate() {
            raw_string_into(source, tag, &mut self.barcodes[index])?;
        }
        for (index, tag) in [*b"QT", *b"CY", *b"BZ", *b"QX"].into_iter().enumerate() {
            raw_string_into(source, tag, &mut self.barcode_qualities[index])?;
        }
        raw_string_into(source, *b"RG", &mut self.read_group)?;
        if let Some(tag) = split_tag {
            raw_string_into(source, tag, &mut self.split_value)?;
            if self.split_value.is_none() {
                return Err(RsomicsError::InvalidInput(format!(
                    "record {} is missing split tag {}",
                    String::from_utf8_lossy(&self.name),
                    String::from_utf8_lossy(&tag)
                )));
            }
        } else {
            self.split_value = None;
        }
        Ok(())
    }

    pub(crate) fn decode_hts(
        &mut self,
        source: &rust_htslib::bam::Record,
        split_tag: Option<[u8; 2]>,
    ) -> Result<()> {
        self.name.clear();
        self.name.extend_from_slice(source.qname());
        self.flags = source.flags();
        self.reference = source.tid();
        self.position = source.pos();
        self.mapping_quality = source.mapq();
        self.mate_reference = source.mtid();
        self.mate_position = source.mpos();
        self.template_length = source.insert_size();
        self.cigar.clear();
        self.cigar.extend(source.cigar().iter().map(|operation| {
            let kind = match operation {
                Cigar::Match(_) => 0,
                Cigar::Ins(_) => 1,
                Cigar::Del(_) => 2,
                Cigar::RefSkip(_) => 3,
                Cigar::SoftClip(_) => 4,
                Cigar::HardClip(_) => 5,
                Cigar::Pad(_) => 6,
                Cigar::Equal(_) => 7,
                Cigar::Diff(_) => 8,
            };
            (kind, operation.len())
        }));
        let sequence = source.seq();
        self.sequence.clear();
        self.sequence
            .extend((0..sequence.len()).map(|index| sequence.encoded_base(index)));
        self.packed_sequence.clear();
        self.packed_sequence.extend(
            self.sequence
                .chunks(2)
                .map(|bases| (bases[0] << 4) | bases.get(1).copied().unwrap_or(0)),
        );
        self.qualities.clear();
        if !source.qual().iter().all(|&quality| quality == 0xff) {
            self.qualities.extend_from_slice(source.qual());
        }
        self.edit_distance = hts_integer(source, *b"NM")?;
        for (index, tag) in [*b"BC", *b"CR", *b"OX", *b"RX"].into_iter().enumerate() {
            hts_string_into(source, tag, &mut self.barcodes[index])?;
        }
        for (index, tag) in [*b"QT", *b"CY", *b"BZ", *b"QX"].into_iter().enumerate() {
            hts_string_into(source, tag, &mut self.barcode_qualities[index])?;
        }
        hts_string_into(source, *b"RG", &mut self.read_group)?;
        if let Some(tag) = split_tag {
            hts_string_into(source, tag, &mut self.split_value)?;
            if self.split_value.is_none() {
                return Err(RsomicsError::InvalidInput(format!(
                    "record {} is missing split tag {}",
                    String::from_utf8_lossy(&self.name),
                    String::from_utf8_lossy(&tag)
                )));
            }
        } else {
            self.split_value = None;
        }
        Ok(())
    }

    pub(crate) fn reference_end(&self) -> Result<i64> {
        let span = self.cigar.iter().try_fold(0i64, |span, &(kind, count)| {
            if matches!(kind, 0 | 2 | 3 | 7 | 8) {
                span.checked_add(i64::from(count))
            } else {
                Some(span)
            }
        });
        let span = span.ok_or_else(|| {
            RsomicsError::InvalidInput("alignment reference span overflows".to_owned())
        })?;
        self.position.checked_add(span.max(1)).ok_or_else(|| {
            RsomicsError::InvalidInput("alignment end position overflows".to_owned())
        })
    }
}

fn pack_sequence(sequence: &[u8]) -> Vec<u8> {
    sequence
        .chunks(2)
        .map(|bases| (bases[0] << 4) | bases.get(1).copied().unwrap_or(0))
        .collect()
}

fn raw_string_into(
    source: &RecordRef<'_>,
    tag: [u8; 2],
    destination: &mut Option<Vec<u8>>,
) -> Result<()> {
    let Some(value) = source.aux_value(tag) else {
        *destination = None;
        return Ok(());
    };
    if source.aux_type(tag) != Some(b'Z') {
        return Err(RsomicsError::InvalidInput(format!(
            "tag {} must be a string",
            String::from_utf8_lossy(&tag)
        )));
    }
    let value = value.strip_suffix(&[0]).ok_or_else(|| {
        RsomicsError::InvalidInput(format!(
            "tag {} has no string terminator",
            String::from_utf8_lossy(&tag)
        ))
    })?;
    let buffer = destination.get_or_insert_with(Vec::new);
    buffer.clear();
    buffer.extend_from_slice(value);
    Ok(())
}

fn raw_integer(source: &RecordRef<'_>, tag: [u8; 2]) -> Result<Option<u64>> {
    let Some(value) = source.aux_value(tag) else {
        return Ok(None);
    };
    let signed = match source.aux_type(tag) {
        Some(b'c') if value.len() == 1 => i64::from(value[0] as i8),
        Some(b'C') if value.len() == 1 => i64::from(value[0]),
        Some(b's') if value.len() == 2 => i64::from(i16::from_le_bytes(value.try_into().unwrap())),
        Some(b'S') if value.len() == 2 => i64::from(u16::from_le_bytes(value.try_into().unwrap())),
        Some(b'i') if value.len() == 4 => i64::from(i32::from_le_bytes(value.try_into().unwrap())),
        Some(b'I') if value.len() == 4 => i64::from(u32::from_le_bytes(value.try_into().unwrap())),
        _ => {
            return Err(RsomicsError::InvalidInput(
                "NM auxiliary tag is not an integer".to_owned(),
            ));
        }
    };
    u64::try_from(signed)
        .map(Some)
        .map_err(|_| RsomicsError::InvalidInput("negative NM auxiliary tag".to_owned()))
}

fn hts_string_into(
    source: &rust_htslib::bam::Record,
    tag: [u8; 2],
    destination: &mut Option<Vec<u8>>,
) -> Result<()> {
    let value = match source.aux(&tag) {
        Ok(Aux::String(value)) => value.as_bytes(),
        Err(rust_htslib::errors::Error::BamAuxTagNotFound) => {
            *destination = None;
            return Ok(());
        }
        Ok(_) => {
            return Err(RsomicsError::InvalidInput(format!(
                "tag {} must be a string",
                String::from_utf8_lossy(&tag)
            )));
        }
        Err(error) => {
            return Err(RsomicsError::InvalidInput(format!(
                "reading tag {}: {error}",
                String::from_utf8_lossy(&tag)
            )));
        }
    };
    let buffer = destination.get_or_insert_with(Vec::new);
    buffer.clear();
    buffer.extend_from_slice(value);
    Ok(())
}

fn hts_integer(source: &rust_htslib::bam::Record, tag: [u8; 2]) -> Result<Option<u64>> {
    let value = match source.aux(&tag) {
        Ok(Aux::I8(value)) => i64::from(value),
        Ok(Aux::U8(value)) => i64::from(value),
        Ok(Aux::I16(value)) => i64::from(value),
        Ok(Aux::U16(value)) => i64::from(value),
        Ok(Aux::I32(value)) => i64::from(value),
        Ok(Aux::U32(value)) => i64::from(value),
        Err(rust_htslib::errors::Error::BamAuxTagNotFound) => return Ok(None),
        Ok(_) => {
            return Err(RsomicsError::InvalidInput(
                "NM auxiliary tag is not an integer".to_owned(),
            ));
        }
        Err(error) => {
            return Err(RsomicsError::InvalidInput(format!(
                "reading NM auxiliary tag: {error}"
            )));
        }
    };
    u64::try_from(value)
        .map(Some)
        .map_err(|_| RsomicsError::InvalidInput("negative NM auxiliary tag".to_owned()))
}

fn string_field(value: Option<&Value>, tag: &[u8; 2]) -> Result<Option<Vec<u8>>> {
    match value {
        None => Ok(None),
        Some(Value::String(value)) => Ok(Some(value.to_vec())),
        Some(_) => Err(RsomicsError::InvalidInput(format!(
            "tag {} must be a string",
            String::from_utf8_lossy(tag)
        ))),
    }
}

fn option_index(value: Option<usize>) -> Result<i32> {
    value.map_or(Ok(-1), |value| {
        i32::try_from(value)
            .map_err(|_| RsomicsError::InvalidInput("reference index overflows".to_owned()))
    })
}

fn cigar_code(kind: Kind) -> u8 {
    match kind {
        Kind::Match => 0,
        Kind::Insertion => 1,
        Kind::Deletion => 2,
        Kind::Skip => 3,
        Kind::SoftClip => 4,
        Kind::HardClip => 5,
        Kind::Pad => 6,
        Kind::SequenceMatch => 7,
        Kind::SequenceMismatch => 8,
    }
}

fn sequence_code(base: u8) -> u8 {
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

fn integer(value: &Value) -> Result<u64> {
    let value = match value {
        Value::Int8(value) => i64::from(*value),
        Value::UInt8(value) => i64::from(*value),
        Value::Int16(value) => i64::from(*value),
        Value::UInt16(value) => i64::from(*value),
        Value::Int32(value) => i64::from(*value),
        Value::UInt32(value) => i64::from(*value),
        _ => {
            return Err(RsomicsError::InvalidInput(
                "NM auxiliary tag is not an integer".to_owned(),
            ));
        }
    };
    u64::try_from(value)
        .map_err(|_| RsomicsError::InvalidInput("negative NM auxiliary tag".to_owned()))
}
