use std::io::Write;

use noodles::core::Region;
use rsomics_bamio::raw::{RawRecord, RecordRef};
use rsomics_common::{Result, RsomicsError};
use rust_htslib::bam;

use super::md;
use super::{ReferenceSummary, Summary};

const UNMAPPED: u16 = 0x04;

#[derive(Clone)]
pub(super) struct Reference {
    pub(super) name: Box<[u8]>,
    pub(super) length: usize,
}

#[derive(Clone)]
pub(super) struct Selection {
    pub(super) reference_id: usize,
    pub(super) start: usize,
    pub(super) end: usize,
    label: Box<[u8]>,
}

impl Selection {
    pub(super) fn new(references: &[Reference], region: &Region) -> Result<Self> {
        let reference_id = references
            .iter()
            .position(|reference| reference.name.as_ref() == region.name())
            .ok_or_else(|| {
                RsomicsError::ConfigError(format!(
                    "region reference is absent from the alignment header: {}",
                    String::from_utf8_lossy(region.name())
                ))
            })?;
        let length = references[reference_id].length;
        let start = region
            .interval()
            .start()
            .map_or(0, |position| usize::from(position) - 1);
        let end = region.interval().end().map_or(length, usize::from);
        if start >= end || end > length {
            return Err(RsomicsError::ConfigError(format!(
                "region interval {start}..{end} is outside reference length {length}"
            )));
        }
        let label = if start == 0 && end == length {
            region.name().to_vec().into_boxed_slice()
        } else {
            format!(
                "{}:{}-{end}",
                String::from_utf8_lossy(region.name()),
                start + 1
            )
            .into_bytes()
            .into_boxed_slice()
        };
        Ok(Self {
            reference_id,
            start,
            end,
            label,
        })
    }
}

struct WorkingReference {
    id: usize,
    name: Box<[u8]>,
    sequence: Vec<u8>,
    coordinate_start: usize,
    reference_length: usize,
    last_alignment_start: Option<usize>,
}

pub(super) struct Builder<W: Write> {
    references: Vec<Reference>,
    pub(super) selection: Option<Selection>,
    current: Option<WorkingReference>,
    output: W,
    summary: Summary,
    cigar: Vec<(u8, u32)>,
}

impl<W: Write> Builder<W> {
    pub(super) fn new(references: Vec<Reference>, selection: Option<Selection>, output: W) -> Self {
        Self {
            references,
            selection,
            current: None,
            output,
            summary: Summary::default(),
            cigar: Vec::new(),
        }
    }

    pub(super) fn add(&mut self, record: &impl EvidenceRecord) -> Result<bool> {
        if record.flags() & UNMAPPED != 0 || record.reference_sequence_id() < 0 {
            return Ok(true);
        }
        let reference_id = usize::try_from(record.reference_sequence_id()).map_err(|_| {
            RsomicsError::InvalidInput("alignment reference ID is out of range".to_owned())
        })?;
        if reference_id >= self.references.len() {
            return Err(RsomicsError::InvalidInput(format!(
                "alignment reference ID {reference_id} is absent from the header"
            )));
        }
        if self
            .selection
            .as_ref()
            .is_some_and(|selection| selection.reference_id != reference_id)
        {
            return Ok(true);
        }
        if self
            .current
            .as_ref()
            .is_some_and(|current| reference_id < current.id)
        {
            return Err(RsomicsError::InvalidInput(
                "reference reconstruction requires coordinate-sorted alignments".to_owned(),
            ));
        }
        if self
            .current
            .as_ref()
            .is_none_or(|current| current.id != reference_id)
        {
            self.flush()?;
            self.current = Some(self.start(reference_id)?);
        }

        record.decode_cigar_into(&mut self.cigar)?;
        let start = usize::try_from(record.alignment_start()).map_err(|_| {
            RsomicsError::InvalidInput("mapped alignment has a negative position".to_owned())
        })?;
        let current = self.current.as_mut().unwrap();
        if current
            .last_alignment_start
            .is_some_and(|previous| start < previous)
        {
            return Err(RsomicsError::InvalidInput(
                "reference reconstruction requires coordinate-sorted alignments".to_owned(),
            ));
        }
        current.last_alignment_start = Some(start);
        let Some(aux_type) = record.aux_type(*b"MD") else {
            return Ok(true);
        };
        if aux_type != b'Z' {
            return Err(RsomicsError::InvalidInput(
                "MD auxiliary field must have type Z".to_owned(),
            ));
        }
        let md = record.aux_value(*b"MD").ok_or_else(|| {
            RsomicsError::InvalidInput("MD auxiliary field cannot be decoded".to_owned())
        })?;
        let md = md.strip_suffix(&[0]).unwrap_or(md);
        md::apply(
            record,
            &self.cigar,
            md,
            start,
            &mut current.sequence,
            current.coordinate_start,
            current.reference_length,
        )?;
        Ok(true)
    }

    fn start(&self, reference_id: usize) -> Result<WorkingReference> {
        let reference = &self.references[reference_id];
        let (name, coordinate_start, coordinate_end) = match &self.selection {
            Some(selection) => (selection.label.clone(), selection.start, selection.end),
            None => (reference.name.clone(), 0, reference.length),
        };
        let sequence = vec![b'N'; coordinate_end - coordinate_start];
        Ok(WorkingReference {
            id: reference_id,
            name,
            sequence,
            coordinate_start,
            reference_length: reference.length,
            last_alignment_start: None,
        })
    }

    pub(super) fn add_embedded(
        &mut self,
        reference_id: i32,
        start: i64,
        data: &[u8],
    ) -> Result<()> {
        let reference_id = usize::try_from(reference_id).map_err(|_| {
            RsomicsError::InvalidInput(
                "embedded reference slice has no single reference ID".to_owned(),
            )
        })?;
        if reference_id >= self.references.len() {
            return Err(RsomicsError::InvalidInput(format!(
                "embedded reference ID {reference_id} is absent from the header"
            )));
        }
        if self
            .selection
            .as_ref()
            .is_some_and(|selection| selection.reference_id != reference_id)
        {
            return Ok(());
        }
        if self
            .current
            .as_ref()
            .is_some_and(|current| reference_id < current.id)
        {
            return Err(RsomicsError::InvalidInput(
                "embedded reference slices are not coordinate sorted".to_owned(),
            ));
        }
        if self
            .current
            .as_ref()
            .is_none_or(|current| current.id != reference_id)
        {
            self.flush()?;
            self.current = Some(self.start(reference_id)?);
        }
        let start = usize::try_from(start)
            .ok()
            .and_then(|start| start.checked_sub(1))
            .ok_or_else(|| {
                RsomicsError::InvalidInput(
                    "embedded reference slice has an invalid start".to_owned(),
                )
            })?;
        let current = self.current.as_mut().unwrap();
        let data_end = start.checked_add(data.len()).ok_or_else(|| {
            RsomicsError::InvalidInput("embedded reference slice length overflows".to_owned())
        })?;
        if start >= current.reference_length {
            return Err(RsomicsError::InvalidInput(
                "embedded reference slice starts beyond the declared reference length".to_owned(),
            ));
        }
        let overlap_start = start.max(current.coordinate_start);
        let overlap_end = data_end
            .min(current.reference_length)
            .min(current.coordinate_start + current.sequence.len());
        for position in overlap_start..overlap_end {
            let base = data[position - start].to_ascii_uppercase();
            set_selected_base(current, position, base)?;
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        let Some(reference) = self.current.take() else {
            return Ok(());
        };
        let sequence = &reference.sequence;
        self.output.write_all(b">")?;
        self.output.write_all(&reference.name)?;
        self.output.write_all(b"\n")?;
        for line in sequence.chunks(60) {
            self.output.write_all(line)?;
            self.output.write_all(b"\n")?;
        }
        let known = sequence.iter().filter(|base| **base != b'N').count() as u64;
        let bases = sequence.len() as u64;
        let coverage = if bases == 0 {
            0.0
        } else {
            known as f64 * 100.0 / bases as f64
        };
        self.summary.references += 1;
        self.summary.bases += bases;
        self.summary.known_bases += known;
        self.summary.items.push(ReferenceSummary {
            name: String::from_utf8_lossy(&reference.name).into_owned(),
            bases,
            known_bases: known,
            coverage,
        });
        Ok(())
    }

    pub(super) fn finish(mut self) -> Result<Summary> {
        if self.current.is_none()
            && let Some(selection) = &self.selection
        {
            self.current = Some(self.start(selection.reference_id)?);
        }
        self.flush()?;
        self.output.flush()?;
        Ok(self.summary)
    }
}

pub(super) trait EvidenceRecord {
    fn flags(&self) -> u16;
    fn reference_sequence_id(&self) -> i32;
    fn alignment_start(&self) -> i32;
    fn decode_cigar_into(&self, cigar: &mut Vec<(u8, u32)>) -> Result<()>;
    fn sequence_len(&self) -> usize;
    fn seq_nibble(&self, index: usize) -> u8;
    fn aux_value(&self, tag: [u8; 2]) -> Option<&[u8]>;
    fn aux_type(&self, tag: [u8; 2]) -> Option<u8>;
}

macro_rules! evidence_record {
    ($record:ty) => {
        impl EvidenceRecord for $record {
            fn flags(&self) -> u16 {
                self.flags()
            }
            fn reference_sequence_id(&self) -> i32 {
                self.reference_sequence_id()
            }
            fn alignment_start(&self) -> i32 {
                self.alignment_start()
            }
            fn decode_cigar_into(&self, cigar: &mut Vec<(u8, u32)>) -> Result<()> {
                self.decode_cigar_into(cigar)
            }
            fn sequence_len(&self) -> usize {
                self.sequence_len()
            }
            fn seq_nibble(&self, index: usize) -> u8 {
                self.seq_nibble(index)
            }
            fn aux_value(&self, tag: [u8; 2]) -> Option<&[u8]> {
                self.aux_value(tag)
            }
            fn aux_type(&self, tag: [u8; 2]) -> Option<u8> {
                self.aux_type(tag)
            }
        }
    };
}

evidence_record!(RawRecord);
evidence_record!(RecordRef<'_>);

impl EvidenceRecord for bam::Record {
    fn flags(&self) -> u16 {
        self.flags()
    }

    fn reference_sequence_id(&self) -> i32 {
        self.tid()
    }

    fn alignment_start(&self) -> i32 {
        i32::try_from(self.pos()).unwrap_or(i32::MIN)
    }

    fn decode_cigar_into(&self, cigar: &mut Vec<(u8, u32)>) -> Result<()> {
        use rust_htslib::bam::record::Cigar;

        cigar.clear();
        for operation in self.cigar().iter() {
            let pair = match operation {
                Cigar::Match(length) => (0, *length),
                Cigar::Ins(length) => (1, *length),
                Cigar::Del(length) => (2, *length),
                Cigar::RefSkip(length) => (3, *length),
                Cigar::SoftClip(length) => (4, *length),
                Cigar::HardClip(length) => (5, *length),
                Cigar::Pad(length) => (6, *length),
                Cigar::Equal(length) => (7, *length),
                Cigar::Diff(length) => (8, *length),
            };
            cigar.push(pair);
        }
        Ok(())
    }

    fn sequence_len(&self) -> usize {
        self.seq_len()
    }

    fn seq_nibble(&self, index: usize) -> u8 {
        self.seq().encoded_base(index)
    }

    fn aux_value(&self, tag: [u8; 2]) -> Option<&[u8]> {
        match self.aux(&tag).ok()? {
            bam::record::Aux::String(value) => Some(value.as_bytes()),
            _ => None,
        }
    }

    fn aux_type(&self, tag: [u8; 2]) -> Option<u8> {
        self.aux(&tag).ok().map(|value| match value {
            bam::record::Aux::String(_) => b'Z',
            bam::record::Aux::Char(_) => b'A',
            bam::record::Aux::I8(_) => b'c',
            bam::record::Aux::U8(_) => b'C',
            bam::record::Aux::I16(_) => b's',
            bam::record::Aux::U16(_) => b'S',
            bam::record::Aux::I32(_) => b'i',
            bam::record::Aux::U32(_) => b'I',
            bam::record::Aux::Float(_) => b'f',
            bam::record::Aux::Double(_) => b'd',
            bam::record::Aux::HexByteArray(_) => b'H',
            bam::record::Aux::ArrayI8(_)
            | bam::record::Aux::ArrayU8(_)
            | bam::record::Aux::ArrayI16(_)
            | bam::record::Aux::ArrayU16(_)
            | bam::record::Aux::ArrayI32(_)
            | bam::record::Aux::ArrayU32(_)
            | bam::record::Aux::ArrayFloat(_) => b'B',
        })
    }
}

fn set_selected_base(reference: &mut WorkingReference, absolute: usize, base: u8) -> Result<()> {
    let position = absolute - reference.coordinate_start;
    match reference.sequence[position] {
        b'N' => reference.sequence[position] = base,
        existing if existing == base || base == b'N' => {}
        existing => {
            return Err(RsomicsError::InvalidInput(format!(
                "conflicting embedded reference evidence at position {}: {} and {}",
                absolute + 1,
                char::from(existing),
                char::from(base)
            )));
        }
    }
    Ok(())
}
