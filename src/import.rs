use std::io::{BufRead, Write};
use std::path::Path;

use noodles::sam;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::{Flags, MappingQuality};
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::{Data, QualityScores, Sequence};
use rsomics_bamio::raw::RecordRef;
use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::{Format as SequenceFormat, PathReader, Reader, Record};
use serde::Serialize;

use crate::{Program, output};

const PAIRED: u16 = 0x01;
const UNMAPPED: u16 = 0x04;
const MATE_UNMAPPED: u16 = 0x08;
const READ1: u16 = 0x40;
const READ2: u16 = 0x80;
const UNMAPPED_BIN: u16 = 4680;

#[derive(Clone, Copy, Debug)]
pub enum Inputs<'a> {
    Single(&'a Path),
    Paired(&'a Path, &'a Path),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Format {
    Sam,
    Bam,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Compression {
    Default,
    Uncompressed,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OrderTag {
    Integer([u8; 2]),
    String { tag: [u8; 2], width: usize },
}

impl OrderTag {
    pub fn parse(value: &str) -> Result<Self> {
        let (tag, width) = value
            .split_once(':')
            .map_or((value, None), |(tag, width)| (tag, Some(width)));
        let tag = parse_tag(tag)?;
        match width {
            None => Ok(Self::Integer(tag)),
            Some("") => Err(RsomicsError::InvalidInput(
                "order tag width cannot be empty".to_owned(),
            )),
            Some(width) => {
                let width = width.parse::<usize>().map_err(|error| {
                    RsomicsError::InvalidInput(format!("invalid order tag width: {error}"))
                })?;
                if !(1..=20).contains(&width) {
                    return Err(RsomicsError::InvalidInput(
                        "order tag width must be between 1 and 20".to_owned(),
                    ));
                }
                Ok(Self::String { tag, width })
            }
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub format: Format,
    pub compression: Compression,
    pub additional_threads: usize,
    pub read_group: Option<&'a str>,
    pub order: Option<OrderTag>,
    pub program: Option<Program<'a>>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub records: u64,
    pub paired_records: u64,
    pub output_format: Format,
    pub additional_threads: usize,
}

pub fn write<W>(inputs: Inputs<'_>, options: Options<'_>, output: W) -> Result<Summary>
where
    W: Write + Send + 'static,
{
    if options.additional_threads > 256 {
        return Err(RsomicsError::ConfigError(
            "import additional thread count cannot exceed 256".to_owned(),
        ));
    }

    let (mut header, read_group_id) = build_header(options.read_group)?;
    if let Some(program) = options.program {
        program.add_to(&mut header)?;
    }

    let format = match options.format {
        Format::Sam => output::Format::Sam,
        Format::Bam => output::Format::Bam,
    };
    let compression = match options.compression {
        Compression::Default => output::Compression::Default,
        Compression::Uncompressed => output::Compression::Uncompressed,
    };
    let mut writer = output::Writer::new(format, compression, options.additional_threads, output);
    writer.write_header(&header)?;
    let mut emitter = Emitter {
        writer,
        header: &header,
        format: options.format,
        read_group_id: read_group_id.as_deref(),
        order: options.order,
        payload: Vec::with_capacity(512),
        records: 0,
        paired_records: 0,
    };

    match inputs {
        Inputs::Single(path) => write_single(path, &mut emitter)?,
        Inputs::Paired(read1, read2) => write_paired(read1, read2, &mut emitter)?,
    }
    let records = emitter.records;
    let paired_records = emitter.paired_records;
    emitter.writer.finish(&header)?;

    Ok(Summary {
        records,
        paired_records,
        output_format: options.format,
        additional_threads: options.additional_threads,
    })
}

trait FastqReader {
    fn read_record(&mut self) -> Result<Option<Record<'_>>>;
}

impl FastqReader for PathReader {
    fn read_record(&mut self) -> Result<Option<Record<'_>>> {
        self.read_record()
    }
}

impl<R: BufRead> FastqReader for Reader<R> {
    fn read_record(&mut self) -> Result<Option<Record<'_>>> {
        self.read_record()
    }
}

fn open_fastq(path: &Path) -> Result<Box<dyn FastqReader>> {
    if path == Path::new("-") {
        let reader = rsomics_seqio::open_reader(std::io::stdin())?;
        if reader.format() != SequenceFormat::Fastq {
            return Err(RsomicsError::InvalidInput(
                "import input is not FASTQ".to_owned(),
            ));
        }
        Ok(Box::new(reader))
    } else {
        let reader = rsomics_seqio::open_path(path)?;
        if reader.format() != SequenceFormat::Fastq {
            return Err(RsomicsError::InvalidInput(format!(
                "import input {} is not FASTQ",
                path.display()
            )));
        }
        Ok(Box::new(reader))
    }
}

fn write_single(path: &Path, emitter: &mut Emitter<'_, impl Write + Send + 'static>) -> Result<()> {
    let mut reader = open_fastq(path)?;
    while let Some(record) = reader.read_record()? {
        let end = read_end(record.id);
        emitter.write(record, flags_for_single(end))?;
    }
    Ok(())
}

fn write_paired(
    read1: &Path,
    read2: &Path,
    emitter: &mut Emitter<'_, impl Write + Send + 'static>,
) -> Result<()> {
    if read1 == Path::new("-") || read2 == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "paired import requires named FASTQ inputs".to_owned(),
        ));
    }
    let mut first = open_fastq(read1)?;
    let mut second = open_fastq(read2)?;
    loop {
        match (first.read_record()?, second.read_record()?) {
            (Some(r1), Some(r2)) => {
                emitter.write(r1, flags_for_read1(read_end(r1.id)))?;
                emitter.write(r2, flags_for_read2(read_end(r2.id)))?;
            }
            (None, None) => return Ok(()),
            (Some(_), None) => {
                return Err(RsomicsError::InvalidInput(
                    "read-2 FASTQ has fewer records than read-1".to_owned(),
                ));
            }
            (None, Some(_)) => {
                return Err(RsomicsError::InvalidInput(
                    "read-2 FASTQ has more records than read-1".to_owned(),
                ));
            }
        }
    }
}

struct Emitter<'a, W>
where
    W: Write + Send + 'static,
{
    writer: output::Writer<W>,
    header: &'a sam::Header,
    format: Format,
    read_group_id: Option<&'a [u8]>,
    order: Option<OrderTag>,
    payload: Vec<u8>,
    records: u64,
    paired_records: u64,
}

impl<W> Emitter<'_, W>
where
    W: Write + Send + 'static,
{
    fn write(&mut self, record: Record<'_>, flags: u16) -> Result<()> {
        let quality = record.qual.ok_or_else(|| {
            RsomicsError::InvalidInput("import record has no FASTQ qualities".to_owned())
        })?;
        let name = normalized_name(record.id);
        if name.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "FASTQ read name is empty after normalization".to_owned(),
            ));
        }
        let order_value = self
            .order
            .map(|order| order_value(order, self.records))
            .transpose()?;

        match self.format {
            Format::Bam => {
                encode_bam_payload(
                    &mut self.payload,
                    name,
                    record.seq,
                    quality,
                    flags,
                    self.read_group_id,
                    order_value.as_ref(),
                )?;
                let record = RecordRef::from_bytes(&self.payload)?;
                self.writer.write_raw_record(&record)?;
            }
            Format::Sam => {
                let record = build_sam_record(
                    name,
                    record.seq,
                    quality,
                    flags,
                    self.read_group_id,
                    order_value,
                )?;
                self.writer.write_record(self.header, &record)?;
            }
        }
        self.records = self.records.checked_add(1).ok_or_else(|| {
            RsomicsError::InvalidInput("import record count exceeds u64".to_owned())
        })?;
        if flags & PAIRED != 0 {
            self.paired_records = self.paired_records.checked_add(1).ok_or_else(|| {
                RsomicsError::InvalidInput("paired import record count exceeds u64".to_owned())
            })?;
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
enum OrderValue {
    Integer { tag: [u8; 2], value: u32 },
    String { tag: [u8; 2], value: String },
}

fn order_value(order: OrderTag, record: u64) -> Result<OrderValue> {
    match order {
        OrderTag::Integer(tag) => Ok(OrderValue::Integer {
            tag,
            value: record.try_into().map_err(|_| {
                RsomicsError::InvalidInput(
                    "integer order tag exceeds the 32-bit SAM limit; use TAG:WIDTH".to_owned(),
                )
            })?,
        }),
        OrderTag::String { tag, width } => {
            let value = format!("{record:0width$}");
            if value.len() > width {
                return Err(RsomicsError::InvalidInput(format!(
                    "order value {record} exceeds configured width {width}"
                )));
            }
            Ok(OrderValue::String { tag, value })
        }
    }
}

fn build_sam_record(
    name: &[u8],
    sequence: &[u8],
    quality: &[u8],
    flags: u16,
    read_group_id: Option<&[u8]>,
    order: Option<OrderValue>,
) -> Result<sam::alignment::RecordBuf> {
    let mut data = Data::default();
    if let Some(id) = read_group_id {
        data.insert(Tag::READ_GROUP, Value::String(id.to_vec().into()));
    }
    if let Some(order) = order {
        match order {
            OrderValue::Integer { tag, value } => {
                data.insert(Tag::from(tag), Value::from(value));
            }
            OrderValue::String { tag, value } => {
                data.insert(Tag::from(tag), Value::String(value.into()));
            }
        }
    }

    let scores = quality.iter().map(|score| score - 33).collect::<Vec<_>>();
    Ok(sam::alignment::RecordBuf::builder()
        .set_name(name.to_vec())
        .set_flags(Flags::from(flags))
        .set_mapping_quality(MappingQuality::new(0).expect("zero is a valid mapping quality"))
        .set_sequence(Sequence::from(sequence.to_vec()))
        .set_quality_scores(QualityScores::from(scores))
        .set_data(data)
        .build())
}

fn encode_bam_payload(
    output: &mut Vec<u8>,
    name: &[u8],
    sequence: &[u8],
    quality: &[u8],
    flags: u16,
    read_group_id: Option<&[u8]>,
    order: Option<&OrderValue>,
) -> Result<()> {
    let name_length: u8 = name
        .len()
        .checked_add(1)
        .and_then(|length| length.try_into().ok())
        .ok_or_else(|| {
            RsomicsError::InvalidInput(
                "FASTQ read name exceeds the BAM limit of 254 bytes".to_owned(),
            )
        })?;
    let sequence_length = u32::try_from(sequence.len())
        .ok()
        .filter(|length| *length <= i32::MAX as u32)
        .ok_or_else(|| {
            RsomicsError::InvalidInput("FASTQ sequence exceeds the BAM length limit".to_owned())
        })?;
    if sequence.len() != quality.len() {
        return Err(RsomicsError::InvalidInput(format!(
            "FASTQ sequence and quality lengths differ: {} and {}",
            sequence.len(),
            quality.len()
        )));
    }

    output.clear();
    output.extend_from_slice(&(-1i32).to_le_bytes());
    output.extend_from_slice(&(-1i32).to_le_bytes());
    output.push(name_length);
    output.push(0);
    output.extend_from_slice(&UNMAPPED_BIN.to_le_bytes());
    output.extend_from_slice(&0u16.to_le_bytes());
    output.extend_from_slice(&flags.to_le_bytes());
    output.extend_from_slice(&sequence_length.to_le_bytes());
    output.extend_from_slice(&(-1i32).to_le_bytes());
    output.extend_from_slice(&(-1i32).to_le_bytes());
    output.extend_from_slice(&0i32.to_le_bytes());
    output.extend_from_slice(name);
    output.push(0);

    for pair in sequence.chunks(2) {
        let high = nt16(pair[0]) << 4;
        let low = pair.get(1).copied().map_or(0, nt16);
        output.push(high | low);
    }
    output.extend(quality.iter().map(|score| score - 33));
    if let Some(id) = read_group_id {
        append_string_tag(output, *b"RG", id);
    }
    if let Some(order) = order {
        match order {
            OrderValue::Integer { tag, value } => {
                output.extend_from_slice(tag);
                output.push(b'I');
                output.extend_from_slice(&value.to_le_bytes());
            }
            OrderValue::String { tag, value } => append_string_tag(output, *tag, value.as_bytes()),
        }
    }
    Ok(())
}

fn append_string_tag(output: &mut Vec<u8>, tag: [u8; 2], value: &[u8]) {
    output.extend_from_slice(&tag);
    output.push(b'Z');
    output.extend_from_slice(value);
    output.push(0);
}

fn nt16(base: u8) -> u8 {
    match base {
        b'=' => 0,
        b'A' | b'a' => 1,
        b'C' | b'c' => 2,
        b'M' | b'm' => 3,
        b'G' | b'g' => 4,
        b'R' | b'r' => 5,
        b'S' | b's' => 6,
        b'V' | b'v' => 7,
        b'T' | b't' => 8,
        b'W' | b'w' => 9,
        b'Y' | b'y' => 10,
        b'H' | b'h' => 11,
        b'K' | b'k' => 12,
        b'D' | b'd' => 13,
        b'B' | b'b' => 14,
        _ => 15,
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ReadEnd {
    First,
    Last,
    Unknown,
}

fn read_end(id: &[u8]) -> ReadEnd {
    let name = first_field(id);
    if name.ends_with(b"/1") {
        ReadEnd::First
    } else if name.ends_with(b"/2") {
        ReadEnd::Last
    } else {
        ReadEnd::Unknown
    }
}

fn normalized_name(id: &[u8]) -> &[u8] {
    let name = first_field(id);
    if name.ends_with(b"/1") || name.ends_with(b"/2") {
        &name[..name.len() - 2]
    } else {
        name
    }
}

fn first_field(id: &[u8]) -> &[u8] {
    id.split(|byte| byte.is_ascii_whitespace())
        .next()
        .unwrap_or_default()
}

fn flags_for_single(end: ReadEnd) -> u16 {
    match end {
        ReadEnd::First => PAIRED | UNMAPPED | MATE_UNMAPPED | READ1,
        ReadEnd::Last => PAIRED | UNMAPPED | MATE_UNMAPPED | READ2,
        ReadEnd::Unknown => UNMAPPED,
    }
}

fn flags_for_read1(end: ReadEnd) -> u16 {
    let end_flag = match end {
        ReadEnd::Last => READ2,
        ReadEnd::First | ReadEnd::Unknown => READ1,
    };
    PAIRED | UNMAPPED | MATE_UNMAPPED | end_flag
}

fn flags_for_read2(end: ReadEnd) -> u16 {
    let suffix_flag = match end {
        ReadEnd::First => READ1,
        ReadEnd::Last | ReadEnd::Unknown => 0,
    };
    PAIRED | UNMAPPED | MATE_UNMAPPED | READ2 | suffix_flag
}

fn build_header(read_group: Option<&str>) -> Result<(sam::Header, Option<Vec<u8>>)> {
    let mut text = "@HD\tVN:1.6\tSO:unsorted\tGO:query\n".to_owned();
    if let Some(line) = read_group {
        if line.contains(['\r', '\n']) {
            return Err(RsomicsError::InvalidInput(
                "read-group line cannot contain line breaks".to_owned(),
            ));
        }
        if line.starts_with("@RG\t") {
            text.push_str(line);
        } else if line.starts_with('@') {
            return Err(RsomicsError::InvalidInput(
                "read-group line must start with @RG followed by a tab".to_owned(),
            ));
        } else {
            text.push_str("@RG\t");
            text.push_str(line);
        }
        text.push('\n');
    }
    let header = text.parse::<sam::Header>().map_err(|error| {
        RsomicsError::InvalidInput(format!("invalid read-group header: {error}"))
    })?;
    let read_group_id = match read_group {
        Some(_) => {
            if header.read_groups().len() != 1 {
                return Err(RsomicsError::InvalidInput(
                    "read-group line must contain exactly one ID".to_owned(),
                ));
            }
            Some(
                header
                    .read_groups()
                    .get_index(0)
                    .expect("one read group is present")
                    .0
                    .to_vec(),
            )
        }
        None => None,
    };
    Ok((header, read_group_id))
}

fn parse_tag(value: &str) -> Result<[u8; 2]> {
    let bytes = value.as_bytes();
    if bytes.len() != 2 || !bytes[0].is_ascii_alphabetic() || !bytes[1].is_ascii_alphanumeric() {
        return Err(RsomicsError::InvalidInput(format!(
            "SAM tag must match [A-Za-z][A-Za-z0-9], got {value:?}"
        )));
    }
    Ok([bytes[0], bytes[1]])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn name_suffixes_drive_single_input_flags() {
        assert_eq!(normalized_name(b"read/1 comment"), b"read");
        assert_eq!(flags_for_single(read_end(b"read/1")), 77);
        assert_eq!(flags_for_single(read_end(b"read/2")), 141);
        assert_eq!(flags_for_single(read_end(b"read")), 4);
    }

    #[test]
    fn bam_payload_is_validated_and_packs_iupac_bases() {
        let mut payload = Vec::new();
        encode_bam_payload(
            &mut payload,
            b"read",
            b"ACGTN",
            b"IIIII",
            UNMAPPED,
            Some(b"lib1"),
            None,
        )
        .unwrap();
        let record = RecordRef::from_bytes(&payload).unwrap();
        assert_eq!(record.name(), b"read");
        assert_eq!(record.flags(), UNMAPPED);
        assert_eq!(record.sequence_len(), 5);
        assert_eq!(record.aux_type(*b"RG"), Some(b'Z'));
        assert_eq!(record.aux_value(*b"RG"), Some(b"lib1\0".as_slice()));
    }

    #[test]
    fn read_group_requires_an_id() {
        assert!(build_header(Some("SM:sample")).is_err());
        let (header, id) = build_header(Some("ID:lib1\tSM:sample")).unwrap();
        assert_eq!(header.read_groups().len(), 1);
        assert_eq!(id.as_deref(), Some(b"lib1".as_slice()));
    }

    #[test]
    fn order_tag_parser_is_strict() {
        assert_eq!(OrderTag::parse("ro").unwrap(), OrderTag::Integer(*b"ro"));
        assert_eq!(
            OrderTag::parse("ro:12").unwrap(),
            OrderTag::String {
                tag: *b"ro",
                width: 12
            }
        );
        assert!(OrderTag::parse("1o").is_err());
        assert!(OrderTag::parse("ro:0").is_err());
    }
}
