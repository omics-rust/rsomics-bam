use std::io::{self, Read};

use flate2::{CrcReader, read::GzDecoder};

use super::encoding;
use super::varint::{
    read_itf8, read_nonnegative_itf8, read_nonnegative_ltf8, read_u8, read_u32_le,
};
use super::{Accumulator, BlockSummary, Method, MethodSummary, Report, invalid};

const MAX_CONTAINER_SIZE: usize = 2 * 1024 * 1024 * 1024;
const MAX_HEADER_BLOCK_SIZE: usize = 64 * 1024 * 1024;

#[derive(Clone, Copy)]
struct Version {
    major: u8,
    minor: u8,
}

#[derive(Debug)]
struct ContainerHeader {
    length: usize,
    reference_id: i32,
    alignment_start: i32,
    records: u64,
    bases: u64,
    blocks: usize,
    landmarks: usize,
    crc: Option<u32>,
}

struct Block {
    content_type: u8,
    content_id: i32,
    compressed_size: usize,
    uncompressed_size: usize,
    method: Method,
    compression_code: u8,
    data: Option<Vec<u8>>,
}

struct CountingReader<R> {
    inner: R,
    position: u64,
}

impl<R> CountingReader<R> {
    fn new(inner: R) -> Self {
        Self { inner, position: 0 }
    }
}

impl<R: Read> Read for CountingReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let count = self.inner.read(buffer)?;
        self.position = self
            .position
            .checked_add(count as u64)
            .ok_or_else(|| invalid("CRAM input position overflow"))?;
        Ok(count)
    }
}

pub(super) fn parse(reader: impl Read) -> io::Result<Report> {
    let mut reader = CountingReader::new(reader);
    let version = read_file_definition(&mut reader)?;
    read_file_header(&mut reader, version)?;

    let mut accumulator = Accumulator::default();
    loop {
        let header = read_container_header(&mut reader, version)?;
        let eof = is_eof(&header, version);
        let mut body = (&mut reader).take(header.length as u64);
        if eof {
            for _ in 0..header.blocks {
                let _ = read_block(&mut body, version, false)?;
            }
            discard(&mut body)?;
            if body.limit() != 0 {
                return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
            }
            let mut trailing = [0];
            if reader.read(&mut trailing)? != 0 {
                return Err(invalid("trailing bytes after the CRAM EOF container"));
            }
            break;
        }

        parse_data_container(&mut body, version, &header, &mut accumulator)?;
        if body.limit() != 0 {
            return Err(invalid(format!(
                "{} unread bytes in CRAM container",
                body.limit()
            )));
        }
    }

    build_report(version, accumulator, reader.position)
}

fn read_file_definition(reader: &mut impl Read) -> io::Result<Version> {
    let mut magic = [0; 4];
    reader.read_exact(&mut magic)?;
    if magic != *b"CRAM" {
        return Err(invalid("missing CRAM magic number"));
    }
    let version = Version {
        major: read_u8(reader)?,
        minor: read_u8(reader)?,
    };
    if !matches!((version.major, version.minor), (2, 1) | (3, 0) | (3, 1)) {
        return Err(invalid(format!(
            "unsupported CRAM version {}.{}",
            version.major, version.minor
        )));
    }
    let mut file_id = [0; 20];
    reader.read_exact(&mut file_id)?;
    Ok(version)
}

fn read_file_header(reader: &mut impl Read, version: Version) -> io::Result<()> {
    let header = read_container_header(reader, version)?;
    if is_eof(&header, version) || header.blocks == 0 {
        return Err(invalid("missing CRAM file header block"));
    }
    let mut body = reader.take(header.length as u64);
    let block = read_block(&mut body, version, false)?;
    if block.content_type != 0 {
        return Err(invalid(format!(
            "expected CRAM file header block, found content type {}",
            block.content_type
        )));
    }
    for _ in 1..header.blocks {
        let _ = read_block(&mut body, version, false)?;
    }
    discard(&mut body)?;
    if body.limit() == 0 {
        Ok(())
    } else {
        Err(io::Error::from(io::ErrorKind::UnexpectedEof))
    }
}

fn read_container_header(reader: &mut impl Read, version: Version) -> io::Result<ContainerHeader> {
    if version.major >= 3 {
        let mut digest = CrcReader::new(reader);
        let mut header = read_container_fields(&mut digest, version)?;
        let actual = digest.crc().sum();
        let expected = read_u32_le(digest.get_mut())?;
        if actual != expected {
            return Err(invalid(format!(
                "CRAM container header checksum mismatch: expected {expected:08x}, got {actual:08x}"
            )));
        }
        header.crc = Some(actual);
        Ok(header)
    } else {
        read_container_fields(reader, version)
    }
}

fn read_container_fields(reader: &mut impl Read, version: Version) -> io::Result<ContainerHeader> {
    let mut length = [0; 4];
    reader.read_exact(&mut length)?;
    let length = usize::try_from(i32::from_le_bytes(length))
        .map_err(|_| invalid("negative CRAM container length"))?;
    if length > MAX_CONTAINER_SIZE {
        return Err(invalid(format!("oversized CRAM container: {length}")));
    }
    let reference_id = read_itf8(reader)?;
    let alignment_start = read_itf8(reader)?;
    let _alignment_span = read_nonnegative_itf8(reader, "alignment span")?;
    let records = read_nonnegative_itf8(reader, "record count")? as u64;
    if version.major == 2 {
        let _ = read_nonnegative_itf8(reader, "record counter")?;
    } else {
        let _ = read_nonnegative_ltf8(reader, "record counter")?;
    }
    let bases = read_nonnegative_ltf8(reader, "base count")?;
    let blocks = read_nonnegative_itf8(reader, "block count")?;
    let landmarks = read_nonnegative_itf8(reader, "landmark count")?;
    if blocks > 1_000_000 || landmarks > 100_000 {
        return Err(invalid("oversized CRAM container index"));
    }
    for _ in 0..landmarks {
        let _ = read_nonnegative_itf8(reader, "landmark")?;
    }
    Ok(ContainerHeader {
        length,
        reference_id,
        alignment_start,
        records,
        bases,
        blocks,
        landmarks,
        crc: None,
    })
}

fn is_eof(header: &ContainerHeader, version: Version) -> bool {
    header.length == if version.major == 2 { 11 } else { 15 }
        && header.reference_id == -1
        && header.alignment_start == 4_542_278
        && header.records == 0
        && header.blocks == 1
        && (header.crc.is_none() || header.crc == Some(0x4fd9_bd05))
}

fn parse_data_container(
    reader: &mut impl Read,
    version: Version,
    header: &ContainerHeader,
    accumulator: &mut Accumulator,
) -> io::Result<()> {
    if header.landmarks == 0 {
        return Err(invalid("data container has no slices"));
    }
    let compression_header = read_block(reader, version, true)?;
    if compression_header.content_type != 1 {
        return Err(invalid(format!(
            "expected compression header block, found content type {}",
            compression_header.content_type
        )));
    }
    let decoded = decode_header_block(&compression_header)?;
    let (encodings, mappings) = encoding::parse(&decoded)?;
    accumulator.encodings.push(encodings);
    for (content_id, data_series) in mappings {
        let series = accumulator.data_series.entry(content_id).or_default();
        if !series.contains(&data_series) {
            if series.is_empty() {
                series.push(data_series);
            } else {
                series.insert(0, data_series);
            }
        }
    }

    accumulator.containers += 1;
    accumulator.slices = accumulator
        .slices
        .checked_add(header.landmarks as u64)
        .ok_or_else(|| invalid("slice count overflow"))?;
    accumulator.sequences = accumulator
        .sequences
        .checked_add(header.records)
        .ok_or_else(|| invalid("sequence count overflow"))?;
    accumulator.bases = accumulator
        .bases
        .checked_add(header.bases)
        .ok_or_else(|| invalid("base count overflow"))?;

    for _ in 0..header.landmarks {
        let slice_header = read_block(reader, version, true)?;
        if !matches!(slice_header.content_type, 2 | 3) {
            return Err(invalid(format!(
                "expected slice header block, found content type {}",
                slice_header.content_type
            )));
        }
        let decoded = decode_header_block(&slice_header)?;
        let (block_count, embedded_reference) =
            parse_slice_header(&decoded, slice_header.content_type, version)?;
        if let Some(content_id) = embedded_reference {
            match accumulator.embedded_reference {
                Some(previous) if previous != content_id => {
                    return Err(invalid("inconsistent embedded-reference content IDs"));
                }
                None => accumulator.embedded_reference = Some(content_id),
                _ => {}
            }
        }
        for _ in 0..block_count {
            let block = read_block(reader, version, false)?;
            let content_id = match block.content_type {
                4 => block.content_id,
                5 => -1,
                other => {
                    return Err(invalid(format!(
                        "unexpected slice data block content type {other}"
                    )));
                }
            };
            let methods = accumulator.methods.entry(content_id).or_default();
            if let Some((_, uncompressed, compressed)) = methods
                .iter_mut()
                .find(|(method, _, _)| *method == block.method)
            {
                *uncompressed = uncompressed
                    .checked_add(block.uncompressed_size as u64)
                    .ok_or_else(|| invalid("uncompressed block size overflow"))?;
                *compressed = compressed
                    .checked_add(block.compressed_size as u64)
                    .ok_or_else(|| invalid("compressed block size overflow"))?;
            } else {
                methods.push((
                    block.method,
                    block.uncompressed_size as u64,
                    block.compressed_size as u64,
                ));
            }
        }
    }
    Ok(())
}

fn read_block(reader: &mut impl Read, version: Version, keep_data: bool) -> io::Result<Block> {
    let mut digest = CrcReader::new(reader);
    let compression_code = read_u8(&mut digest)?;
    let content_type = read_u8(&mut digest)?;
    let content_id = read_itf8(&mut digest)?;
    let compressed_size = read_nonnegative_itf8(&mut digest, "compressed block size")?;
    let uncompressed_size = read_nonnegative_itf8(&mut digest, "uncompressed block size")?;
    if compression_code == 0 && compressed_size != uncompressed_size {
        return Err(invalid("raw CRAM block sizes differ"));
    }
    if keep_data && compressed_size > MAX_HEADER_BLOCK_SIZE {
        return Err(invalid(format!(
            "oversized CRAM header block: {compressed_size}"
        )));
    }

    let mut data = keep_data.then(|| Vec::with_capacity(compressed_size));
    let mut prefix = Vec::with_capacity(compressed_size.min(9));
    let mut remaining = compressed_size;
    let mut buffer = [0; 64 * 1024];
    while remaining > 0 {
        let count = remaining.min(buffer.len());
        digest.read_exact(&mut buffer[..count])?;
        if prefix.len() < 9 {
            let take = (9 - prefix.len()).min(count);
            prefix.extend_from_slice(&buffer[..take]);
        }
        if let Some(data) = &mut data {
            data.extend_from_slice(&buffer[..count]);
        }
        remaining -= count;
    }
    let actual = digest.crc().sum();
    if version.major >= 3 {
        let expected = read_u32_le(digest.get_mut())?;
        if actual != expected {
            return Err(invalid(format!(
                "CRAM block checksum mismatch: expected {expected:08x}, got {actual:08x}"
            )));
        }
    }
    let method = if uncompressed_size == 0 {
        Method::Raw
    } else {
        Method::classify(compression_code, &prefix)?
    };
    Ok(Block {
        content_type,
        content_id,
        compressed_size,
        uncompressed_size,
        method,
        compression_code,
        data,
    })
}

fn decode_header_block(block: &Block) -> io::Result<Vec<u8>> {
    let data = block
        .data
        .as_deref()
        .ok_or_else(|| invalid("missing retained CRAM header block"))?;
    let decoded = match block.compression_code {
        0 => data.to_vec(),
        1 => {
            let mut decoded = Vec::with_capacity(block.uncompressed_size);
            GzDecoder::new(data).read_to_end(&mut decoded)?;
            decoded
        }
        method => {
            return Err(invalid(format!(
                "unsupported compression method {method} on CRAM header block"
            )));
        }
    };
    if decoded.len() != block.uncompressed_size {
        return Err(invalid(format!(
            "decoded CRAM header size mismatch: expected {}, got {}",
            block.uncompressed_size,
            decoded.len()
        )));
    }
    Ok(decoded)
}

fn parse_slice_header(
    data: &[u8],
    content_type: u8,
    version: Version,
) -> io::Result<(usize, Option<i32>)> {
    let mut reader = data;
    if content_type == 2 {
        let _ = read_itf8(&mut reader)?;
        let _ = read_nonnegative_itf8(&mut reader, "slice alignment start")?;
        let _ = read_nonnegative_itf8(&mut reader, "slice alignment span")?;
    }
    let _ = read_nonnegative_itf8(&mut reader, "slice record count")?;
    if version.major == 2 {
        let _ = read_nonnegative_itf8(&mut reader, "slice record counter")?;
    } else {
        let _ = read_nonnegative_ltf8(&mut reader, "slice record counter")?;
    }
    let block_count = read_nonnegative_itf8(&mut reader, "slice block count")?;
    let content_ids = read_nonnegative_itf8(&mut reader, "slice content ID count")?;
    if content_ids == 0 || content_ids > 100_000 || block_count > 1_000_000 {
        return Err(invalid("invalid slice block index"));
    }
    for _ in 0..content_ids {
        let _ = read_itf8(&mut reader)?;
    }
    let embedded_reference = if content_type == 2 {
        match read_itf8(&mut reader)? {
            -1 => None,
            id => Some(id),
        }
    } else {
        None
    };
    let mut md5 = [0; 16];
    reader.read_exact(&mut md5)?;
    Ok((block_count, embedded_reference))
}

fn build_report(version: Version, accumulator: Accumulator, file_size: u64) -> io::Result<Report> {
    let mut blocks = Vec::with_capacity(accumulator.methods.len());
    let mut compressed_total = 0u64;
    for (content_id, methods) in accumulator.methods {
        let mut methods = methods
            .into_iter()
            .map(|(method, uncompressed_size, compressed_size)| {
                compressed_total = compressed_total
                    .checked_add(compressed_size)
                    .ok_or_else(|| invalid("compressed total overflow"))?;
                Ok((
                    method,
                    MethodSummary {
                        method: method.name().to_owned(),
                        short: method.short().to_owned(),
                        uncompressed_size,
                        compressed_size,
                    },
                ))
            })
            .collect::<io::Result<Vec<_>>>()?;
        methods.sort_by(|(left_method, left), (right_method, right)| {
            right
                .compressed_size
                .cmp(&left.compressed_size)
                .then_with(|| left_method.rank().cmp(&right_method.rank()))
        });
        blocks.push(BlockSummary {
            content_id: (content_id >= 0).then_some(content_id),
            methods: methods.into_iter().map(|(_, summary)| summary).collect(),
            data_series: accumulator
                .data_series
                .get(&content_id)
                .cloned()
                .unwrap_or_default(),
            embedded_reference: accumulator.embedded_reference == Some(content_id),
        });
    }
    let format_overhead_size = file_size
        .checked_sub(compressed_total)
        .ok_or_else(|| invalid("compressed blocks exceed CRAM file size"))?;
    Ok(Report {
        version: format!("{}.{}", version.major, version.minor),
        blocks,
        encodings: accumulator.encodings,
        containers: accumulator.containers,
        slices: accumulator.slices,
        sequences: accumulator.sequences,
        bases: accumulator.bases,
        file_size,
        format_overhead_size,
    })
}

fn discard(reader: &mut impl Read) -> io::Result<()> {
    io::copy(reader, &mut io::sink()).map(|_| ())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_oversized_container_indexes() {
        let mut data = vec![0; 4];
        data.extend([0, 0, 0, 0, 0, 0]);
        put_itf8(&mut data, 1_000_001);
        data.push(0);
        let error = read_container_fields(&mut data.as_slice(), Version { major: 3, minor: 0 })
            .unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn rejects_raw_blocks_with_inconsistent_sizes() {
        let error = read_block(
            &mut &[0, 4, 1, 1, 2][..],
            Version { major: 2, minor: 1 },
            false,
        )
        .err()
        .unwrap();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    }

    fn put_itf8(dst: &mut Vec<u8>, value: u32) {
        if value < 1 << 7 {
            dst.push(value as u8);
        } else if value < 1 << 14 {
            dst.extend([(value >> 8) as u8 | 0x80, value as u8]);
        } else if value < 1 << 21 {
            dst.extend([(value >> 16) as u8 | 0xc0, (value >> 8) as u8, value as u8]);
        } else {
            unreachable!();
        }
    }
}
