use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

use noodles::{bam, bgzf, sam};
use rsomics_common::{Result, RsomicsError};

pub(crate) const EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

const READ_BUFFER: usize = 1024 * 1024;
const MAX_FRAME_SIZE: usize = 1 << 16;

#[derive(Debug)]
struct Frame {
    data: Vec<u8>,
    payload_start: usize,
    uncompressed_size: usize,
}

struct FrameLayout {
    payload_start: usize,
    uncompressed_size: usize,
}

struct OpenedBam {
    header: sam::Header,
    text: Vec<u8>,
    reader: BufReader<File>,
    boundary: Vec<u8>,
}

pub(crate) fn read_header(path: &Path) -> Result<(sam::Header, Vec<u8>)> {
    open_at_records(path).map(|bam| (bam.header, bam.text))
}

pub(crate) fn write_header<W: Write>(
    output: &mut W,
    header: &sam::Header,
    text: &[u8],
) -> Result<()> {
    let mut raw = Vec::new();
    raw.extend_from_slice(b"BAM\x01");
    raw.extend_from_slice(
        &i32::try_from(text.len())
            .map_err(|_| RsomicsError::InvalidInput("BAM header text is too large".to_owned()))?
            .to_le_bytes(),
    );
    raw.extend_from_slice(text);
    raw.extend_from_slice(
        &i32::try_from(header.reference_sequences().len())
            .map_err(|_| RsomicsError::InvalidInput("BAM reference count is too large".to_owned()))?
            .to_le_bytes(),
    );
    for (name, reference) in header.reference_sequences() {
        let name: &[u8] = name.as_ref();
        let name_len = name.len().checked_add(1).ok_or_else(|| {
            RsomicsError::InvalidInput("BAM reference name length overflows".to_owned())
        })?;
        raw.extend_from_slice(
            &i32::try_from(name_len)
                .map_err(|_| {
                    RsomicsError::InvalidInput("BAM reference name is too large".to_owned())
                })?
                .to_le_bytes(),
        );
        raw.extend_from_slice(name);
        raw.push(0);
        raw.extend_from_slice(
            &i32::try_from(usize::from(reference.length()))
                .map_err(|_| {
                    RsomicsError::InvalidInput("BAM reference length exceeds i32".to_owned())
                })?
                .to_le_bytes(),
        );
    }
    write_encoded(output, &raw)
}

pub(crate) fn canonical_header_text(header: &sam::Header) -> Result<Vec<u8>> {
    let mut writer = bam::io::Writer::from(Vec::new());
    writer.write_header(header).map_err(RsomicsError::Io)?;
    header_text(&writer.into_inner())
}

pub(crate) fn copy_records<W: Write>(path: &Path, output: &mut W) -> Result<()> {
    let OpenedBam {
        mut reader,
        boundary,
        ..
    } = open_at_records(path)?;
    write_encoded(output, &boundary)?;
    copy_frames(path, &mut reader, output)
}

pub(crate) fn finish<W: Write>(output: &mut W) -> Result<()> {
    output.write_all(&EOF).map_err(RsomicsError::Io)
}

fn open_at_records(path: &Path) -> Result<OpenedBam> {
    let file = File::open(path).map_err(|error| {
        RsomicsError::InvalidInput(format!("opening {}: {error}", path.display()))
    })?;
    let mut stream = HeaderStream::new(BufReader::with_capacity(READ_BUFFER, file));
    let header = {
        let mut reader = bam::io::Reader::from(&mut stream);
        reader.read_header().map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "reading BAM header from {}: {error}",
                path.display()
            ))
        })?
    };
    let (reader, boundary, consumed) = stream.into_records();
    let text = header_text(&consumed)?;
    Ok(OpenedBam {
        header,
        text,
        reader,
        boundary,
    })
}

fn header_text(raw: &[u8]) -> Result<Vec<u8>> {
    if raw.len() < 8 || &raw[..4] != b"BAM\x01" {
        return Err(RsomicsError::InvalidInput(
            "invalid BAM header bytes".to_owned(),
        ));
    }
    let len = i32::from_le_bytes(raw[4..8].try_into().unwrap());
    let len = usize::try_from(len)
        .map_err(|_| RsomicsError::InvalidInput("negative BAM header length".to_owned()))?;
    let end = 8usize
        .checked_add(len)
        .filter(|end| *end <= raw.len())
        .ok_or_else(|| RsomicsError::InvalidInput("truncated BAM header text".to_owned()))?;
    Ok(raw[8..end].to_vec())
}

fn write_encoded<W: Write>(output: &mut W, data: &[u8]) -> Result<()> {
    if data.is_empty() {
        return Ok(());
    }
    let mut writer = bgzf::io::Writer::new(Vec::new());
    writer.write_all(data).map_err(RsomicsError::Io)?;
    writer.try_finish().map_err(RsomicsError::Io)?;
    let mut encoded = writer.into_inner();
    if !encoded.ends_with(&EOF) {
        return Err(RsomicsError::InvalidInput(
            "BGZF encoder did not produce an end-of-file marker".to_owned(),
        ));
    }
    encoded.truncate(encoded.len() - EOF.len());
    output.write_all(&encoded).map_err(RsomicsError::Io)
}

struct HeaderStream<R> {
    inner: R,
    decoded: Vec<u8>,
    position: usize,
    consumed: Vec<u8>,
}

impl<R> HeaderStream<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            decoded: Vec::new(),
            position: 0,
            consumed: Vec::new(),
        }
    }

    fn into_records(self) -> (R, Vec<u8>, Vec<u8>) {
        (
            self.inner,
            self.decoded[self.position..].to_vec(),
            self.consumed,
        )
    }
}

impl<R: Read> Read for HeaderStream<R> {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        while self.position == self.decoded.len() {
            let Some(frame) = read_frame(&mut self.inner).map_err(invalid_data)? else {
                return Ok(0);
            };
            self.decoded.clear();
            self.position = 0;
            if frame.data == EOF {
                return Ok(0);
            }
            decode_frame(&frame, &mut self.decoded).map_err(invalid_data)?;
        }

        let available = &self.decoded[self.position..];
        let len = available.len().min(buffer.len());
        buffer[..len].copy_from_slice(&available[..len]);
        self.consumed.extend_from_slice(&available[..len]);
        self.position += len;
        Ok(len)
    }
}

fn read_frame<R: Read>(reader: &mut R) -> Result<Option<Frame>> {
    let mut data = Vec::with_capacity(MAX_FRAME_SIZE);
    let Some(layout) = read_frame_into(reader, &mut data)? else {
        return Ok(None);
    };
    Ok(Some(Frame {
        data,
        payload_start: layout.payload_start,
        uncompressed_size: layout.uncompressed_size,
    }))
}

fn read_frame_into<R: Read>(reader: &mut R, data: &mut Vec<u8>) -> Result<Option<FrameLayout>> {
    data.clear();
    data.resize(12, 0);
    loop {
        match reader.read(&mut data[..1]) {
            Ok(0) => {
                data.clear();
                return Ok(None);
            }
            Ok(_) => break,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(RsomicsError::Io(error)),
        }
    }
    reader.read_exact(&mut data[1..]).map_err(truncated_frame)?;
    if data[..3] != [0x1f, 0x8b, 0x08] || data[3] != 0x04 {
        return Err(RsomicsError::InvalidInput(
            "invalid BGZF gzip header".to_owned(),
        ));
    }

    let extra_len = usize::from(u16::from_le_bytes([data[10], data[11]]));
    let extra_end = 12 + extra_len;
    data.resize(extra_end, 0);
    reader
        .read_exact(&mut data[12..extra_end])
        .map_err(truncated_frame)?;
    let block_size = parse_block_size(&data[12..extra_end])?;
    let header_size = extra_end;
    if !(header_size + 8..=MAX_FRAME_SIZE).contains(&block_size) {
        return Err(RsomicsError::InvalidInput(
            "invalid BGZF frame size".to_owned(),
        ));
    }

    data.resize(block_size, 0);
    reader
        .read_exact(&mut data[header_size..])
        .map_err(truncated_frame)?;
    let uncompressed_size = u32::from_le_bytes(data[block_size - 4..].try_into().unwrap());
    let uncompressed_size = usize::try_from(uncompressed_size).unwrap();
    if uncompressed_size > MAX_FRAME_SIZE {
        return Err(RsomicsError::InvalidInput(
            "invalid BGZF uncompressed size".to_owned(),
        ));
    }
    Ok(Some(FrameLayout {
        payload_start: header_size,
        uncompressed_size,
    }))
}

fn copy_frames<W: Write>(path: &Path, reader: &mut BufReader<File>, output: &mut W) -> Result<()> {
    let mut partial = Vec::with_capacity(MAX_FRAME_SIZE);
    loop {
        let buffer = reader.fill_buf().map_err(RsomicsError::Io)?;
        if buffer.is_empty() {
            return Err(RsomicsError::InvalidInput(format!(
                "{}: BGZF end-of-file marker is missing",
                path.display()
            )));
        }

        let mut offset = 0;
        let mut eof = None;
        while offset < buffer.len() {
            let Some(size) = complete_frame_size(&buffer[offset..])? else {
                break;
            };
            if buffer[offset..offset + size] == EOF {
                eof = Some(offset + size);
                break;
            }
            offset += size;
        }

        if let Some(end) = eof {
            output
                .write_all(&buffer[..end - EOF.len()])
                .map_err(RsomicsError::Io)?;
            reader.consume(end);
            let mut trailing = [0; 1];
            if reader.read(&mut trailing).map_err(RsomicsError::Io)? != 0 {
                return Err(RsomicsError::InvalidInput(format!(
                    "{}: data follows the BGZF end-of-file marker",
                    path.display()
                )));
            }
            return Ok(());
        }

        if offset > 0 {
            output
                .write_all(&buffer[..offset])
                .map_err(RsomicsError::Io)?;
            reader.consume(offset);
            continue;
        }

        partial.extend_from_slice(buffer);
        let consumed = buffer.len();
        reader.consume(consumed);
        complete_partial_frame(reader, &mut partial)?;
        if partial == EOF {
            let mut trailing = [0; 1];
            if reader.read(&mut trailing).map_err(RsomicsError::Io)? != 0 {
                return Err(RsomicsError::InvalidInput(format!(
                    "{}: data follows the BGZF end-of-file marker",
                    path.display()
                )));
            }
            return Ok(());
        }
        output.write_all(&partial).map_err(RsomicsError::Io)?;
        partial.clear();
    }
}

fn complete_frame_size(data: &[u8]) -> Result<Option<usize>> {
    if data.len() < 12 {
        return Ok(None);
    }
    if data[..3] != [0x1f, 0x8b, 0x08] || data[3] != 0x04 {
        return Err(RsomicsError::InvalidInput(
            "invalid BGZF gzip header".to_owned(),
        ));
    }
    let extra_len = usize::from(u16::from_le_bytes([data[10], data[11]]));
    let header_size = 12 + extra_len;
    if data.len() < header_size {
        return Ok(None);
    }
    let block_size = parse_block_size(&data[12..header_size])?;
    if !(header_size + 8..=MAX_FRAME_SIZE).contains(&block_size) {
        return Err(RsomicsError::InvalidInput(
            "invalid BGZF frame size".to_owned(),
        ));
    }
    if data.len() < block_size {
        return Ok(None);
    }
    let uncompressed_size =
        u32::from_le_bytes(data[block_size - 4..block_size].try_into().unwrap());
    if usize::try_from(uncompressed_size).unwrap() > MAX_FRAME_SIZE {
        return Err(RsomicsError::InvalidInput(
            "invalid BGZF uncompressed size".to_owned(),
        ));
    }
    Ok(Some(block_size))
}

fn complete_partial_frame<R: Read>(reader: &mut R, data: &mut Vec<u8>) -> Result<()> {
    fill_to(reader, data, 12)?;
    if data[..3] != [0x1f, 0x8b, 0x08] || data[3] != 0x04 {
        return Err(RsomicsError::InvalidInput(
            "invalid BGZF gzip header".to_owned(),
        ));
    }
    let extra_len = usize::from(u16::from_le_bytes([data[10], data[11]]));
    let header_size = 12 + extra_len;
    fill_to(reader, data, header_size)?;
    let block_size = parse_block_size(&data[12..header_size])?;
    if !(header_size + 8..=MAX_FRAME_SIZE).contains(&block_size) {
        return Err(RsomicsError::InvalidInput(
            "invalid BGZF frame size".to_owned(),
        ));
    }
    fill_to(reader, data, block_size)?;
    let uncompressed_size = u32::from_le_bytes(data[block_size - 4..].try_into().unwrap());
    let uncompressed_size = usize::try_from(uncompressed_size).unwrap();
    if uncompressed_size > MAX_FRAME_SIZE {
        return Err(RsomicsError::InvalidInput(
            "invalid BGZF uncompressed size".to_owned(),
        ));
    }
    Ok(())
}

fn fill_to<R: Read>(reader: &mut R, data: &mut Vec<u8>, len: usize) -> Result<()> {
    if data.len() >= len {
        return Ok(());
    }
    let start = data.len();
    data.resize(len, 0);
    reader
        .read_exact(&mut data[start..])
        .map_err(truncated_frame)
}

fn decode_frame(frame: &Frame, output: &mut Vec<u8>) -> Result<()> {
    let trailer = frame.data.len() - 8;
    output.clear();
    flate2::read::DeflateDecoder::new(&frame.data[frame.payload_start..trailer])
        .read_to_end(output)
        .map_err(RsomicsError::Io)?;
    if output.len() != frame.uncompressed_size {
        return Err(RsomicsError::InvalidInput(
            "BGZF uncompressed size does not match its trailer".to_owned(),
        ));
    }
    let expected_crc = u32::from_le_bytes(frame.data[trailer..trailer + 4].try_into().unwrap());
    let mut crc = flate2::Crc::new();
    crc.update(output);
    if crc.sum() != expected_crc {
        return Err(RsomicsError::InvalidInput(
            "BGZF checksum does not match its trailer".to_owned(),
        ));
    }
    Ok(())
}

fn parse_block_size(extra: &[u8]) -> Result<usize> {
    let mut position = 0;
    let mut block_size = None;
    while position < extra.len() {
        if extra.len() - position < 4 {
            return Err(RsomicsError::InvalidInput(
                "truncated BGZF extra subfield".to_owned(),
            ));
        }
        let id = &extra[position..position + 2];
        let len = usize::from(u16::from_le_bytes([
            extra[position + 2],
            extra[position + 3],
        ]));
        position += 4;
        let end = position.checked_add(len).ok_or_else(|| {
            RsomicsError::InvalidInput("BGZF extra subfield size overflows".to_owned())
        })?;
        if end > extra.len() {
            return Err(RsomicsError::InvalidInput(
                "truncated BGZF extra subfield".to_owned(),
            ));
        }
        if id == b"BC" {
            if len != 2 || block_size.is_some() {
                return Err(RsomicsError::InvalidInput(
                    "invalid BGZF BC subfield".to_owned(),
                ));
            }
            block_size =
                Some(usize::from(u16::from_le_bytes([extra[position], extra[position + 1]])) + 1);
        }
        position = end;
    }
    block_size.ok_or_else(|| RsomicsError::InvalidInput("BGZF BC subfield is missing".to_owned()))
}

fn truncated_frame(error: io::Error) -> RsomicsError {
    if error.kind() == io::ErrorKind::UnexpectedEof {
        RsomicsError::InvalidInput("truncated BGZF frame".to_owned())
    } else {
        RsomicsError::Io(error)
    }
}

fn invalid_data(error: impl std::fmt::Display) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, error.to_string())
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn variable_extra_subfields_are_accepted() {
        let mut encoded = Vec::new();
        write_encoded(&mut encoded, b"payload").unwrap();
        let old_size = usize::from(u16::from_le_bytes([encoded[16], encoded[17]])) + 1;
        encoded.splice(12..12, [b'X', b'Y', 0, 0]);
        encoded[10..12].copy_from_slice(&10u16.to_le_bytes());
        encoded[20..22].copy_from_slice(&u16::try_from(old_size + 3).unwrap().to_le_bytes());
        encoded.extend_from_slice(&EOF);

        let mut stream = HeaderStream::new(Cursor::new(encoded));
        let mut decoded = Vec::new();
        stream.read_to_end(&mut decoded).unwrap();
        assert_eq!(decoded, b"payload");
    }

    #[test]
    fn missing_bc_subfield_is_rejected() {
        let mut frame = Vec::new();
        write_encoded(&mut frame, b"payload").unwrap();
        frame[12..14].copy_from_slice(b"XY");
        let error = read_frame(&mut Cursor::new(frame)).unwrap_err();
        assert!(error.to_string().contains("BC subfield is missing"));
    }

    #[test]
    fn truncated_frame_is_rejected() {
        let mut frame = Vec::new();
        write_encoded(&mut frame, b"payload").unwrap();
        frame.truncate(frame.len() - 1);
        let error = read_frame(&mut Cursor::new(frame)).unwrap_err();
        assert!(error.to_string().contains("truncated BGZF frame"));
    }

    #[test]
    fn frame_copy_handles_reader_boundaries() {
        let mut expected = Vec::new();
        write_encoded(&mut expected, &vec![b'a'; 48_000]).unwrap();
        write_encoded(&mut expected, &vec![b'b'; 48_000]).unwrap();
        let mut encoded = expected.clone();
        encoded.extend_from_slice(&EOF);
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(&encoded).unwrap();
        let mut reader = BufReader::with_capacity(37, file.reopen().unwrap());
        let mut actual = Vec::new();
        copy_frames(Path::new("input.bam"), &mut reader, &mut actual).unwrap();
        assert_eq!(actual, expected);
    }
}
