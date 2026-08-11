use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read};
use std::num::NonZero;
use std::path::Path;

use flate2::read::MultiGzDecoder;
use noodles::{bam, bgzf, core::Region, cram, fasta, sam};
use noodles_util::alignment;
use rsomics_bamio::raw::{self, RawRecord, RawRecordEncoder, RecordReader, RecordRef};
use rsomics_bamio::{IndexedAlignmentReader, open_indexed_alignment};
use rsomics_common::{Context, Result, RsomicsError};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Format {
    Sam,
    Bam,
    Cram,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Compression {
    None,
    Bgzf,
    Gzip,
}

enum Inner {
    General(alignment::io::Reader<Box<dyn Read>>),
    Sam(sam::io::Reader<BufReader<File>>),
    SamGz(sam::io::Reader<bgzf::io::Reader<BufReader<File>>>),
    SamGzip(Box<sam::io::Reader<BufReader<MultiGzDecoder<BufReader<File>>>>>),
    Bam(ParallelBamReader),
    BamRaw(bam::io::Reader<BufReader<File>>),
    Cram(cram::io::Reader<BufReader<File>>),
    Indexed(IndexedAlignmentReader),
    IndexedDirect(alignment::io::IndexedReader<File>),
}

type ParallelBamReader = bam::io::Reader<Box<dyn BufRead + Send>>;

const READ_BUFFER: usize = 256 * 1024;

pub(crate) struct Reader {
    inner: Inner,
    format: Format,
}

impl Reader {
    pub(crate) fn format(&self) -> Format {
        self.format
    }

    pub(crate) fn has_reusable_raw_bam_path(&self) -> bool {
        matches!(self.inner, Inner::Bam(_) | Inner::BamRaw(_))
    }

    pub(crate) fn read_header(&mut self, input: &Path) -> Result<sam::Header> {
        let result = match &mut self.inner {
            Inner::General(reader) => reader.read_header(),
            Inner::Sam(reader) => reader.read_header(),
            Inner::SamGz(reader) => reader.read_header(),
            Inner::SamGzip(reader) => reader.read_header(),
            Inner::Bam(reader) => reader.read_header(),
            Inner::BamRaw(reader) => reader.read_header(),
            Inner::Cram(reader) => reader.read_header(),
            Inner::Indexed(reader) => reader.read_header(),
            Inner::IndexedDirect(reader) => reader.read_header(),
        };
        result.map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "reading alignment header from {}: {error}",
                input.display()
            ))
        })
    }

    pub(crate) fn visit_records(
        &mut self,
        header: &sam::Header,
        input: &Path,
        mut visit: impl FnMut(&dyn sam::alignment::Record) -> Result<bool>,
    ) -> Result<()> {
        match &mut self.inner {
            Inner::General(reader) => {
                for result in reader.records(header) {
                    let record = result.map_err(|error| record_error(input, error))?;
                    if !visit(record.as_ref())? {
                        break;
                    }
                }
            }
            Inner::Sam(reader) => {
                for result in reader.records() {
                    let record = result.map_err(|error| record_error(input, error))?;
                    if !visit(&record)? {
                        break;
                    }
                }
            }
            Inner::SamGz(reader) => {
                for result in reader.records() {
                    let record = result.map_err(|error| record_error(input, error))?;
                    if !visit(&record)? {
                        break;
                    }
                }
            }
            Inner::SamGzip(reader) => {
                for result in reader.records() {
                    let record = result.map_err(|error| record_error(input, error))?;
                    if !visit(&record)? {
                        break;
                    }
                }
            }
            Inner::Bam(reader) => {
                for result in reader.records() {
                    let record = result.map_err(|error| record_error(input, error))?;
                    if !visit(&record)? {
                        break;
                    }
                }
            }
            Inner::BamRaw(reader) => {
                for result in reader.records() {
                    let record = result.map_err(|error| record_error(input, error))?;
                    if !visit(&record)? {
                        break;
                    }
                }
            }
            Inner::Cram(reader) => {
                for result in reader.records(header) {
                    let record = result.map_err(|error| record_error(input, error))?;
                    if !visit(&record)? {
                        break;
                    }
                }
            }
            Inner::Indexed(reader) => {
                for result in reader.records(header) {
                    let record = result.map_err(|error| record_error(input, error))?;
                    if !visit(record.as_ref())? {
                        break;
                    }
                }
            }
            Inner::IndexedDirect(reader) => {
                for result in reader.records(header) {
                    let record = result.map_err(|error| record_error(input, error))?;
                    if !visit(record.as_ref())? {
                        break;
                    }
                }
            }
        }
        Ok(())
    }

    pub(crate) fn visit_raw_bam_records(
        &mut self,
        input: &Path,
        visit: impl FnMut(RecordRef<'_>) -> Result<bool>,
    ) -> Result<()> {
        match &mut self.inner {
            Inner::Bam(reader) => visit_raw_records(reader.get_mut(), input, visit),
            Inner::BamRaw(reader) => visit_raw_records(reader.get_mut(), input, visit),
            _ => Err(RsomicsError::ConfigError(
                "raw BAM records require sequential BAM input".to_owned(),
            )),
        }
    }

    pub(crate) fn visit_owned_raw_records(
        &mut self,
        header: &sam::Header,
        input: &Path,
        mut visit: impl FnMut(RawRecord) -> Result<bool>,
    ) -> Result<()> {
        match &mut self.inner {
            Inner::Bam(reader) => {
                return visit_owned_raw_records(reader.get_mut(), input, visit);
            }
            Inner::BamRaw(reader) => {
                return visit_owned_raw_records(reader.get_mut(), input, visit);
            }
            _ => {}
        }

        let mut encoder = RawRecordEncoder::new();
        self.visit_records(header, input, |record| {
            visit(encoder.encode(header, record)?)
        })
    }

    pub(crate) fn visit_mut_raw_bam_records(
        &mut self,
        input: &Path,
        visit: impl FnMut(&mut RawRecord) -> Result<bool>,
    ) -> Result<()> {
        match &mut self.inner {
            Inner::Bam(reader) => visit_mut_raw_records(reader.get_mut(), input, visit),
            Inner::BamRaw(reader) => visit_mut_raw_records(reader.get_mut(), input, visit),
            _ => Err(RsomicsError::ConfigError(
                "mutable raw records require BAM input".to_owned(),
            )),
        }
    }

    pub(crate) fn visit_region(
        &mut self,
        header: &sam::Header,
        input: &Path,
        region: Option<&Region>,
        mut visit: impl FnMut(&dyn sam::alignment::Record) -> Result<bool>,
    ) -> Result<()> {
        let records: Box<dyn Iterator<Item = io::Result<Box<dyn sam::alignment::Record>>> + '_> =
            match &mut self.inner {
                Inner::Indexed(reader) => match region {
                    Some(region) => reader
                        .query(header, region)
                        .map(Box::new)
                        .map_err(|error| query_error(input, region.to_string(), error))?,
                    None => reader
                        .query_unmapped(header)
                        .map(Box::new)
                        .map_err(|error| query_error(input, "*", error))?,
                },
                Inner::IndexedDirect(reader) => match region {
                    Some(region) => reader
                        .query(header, region)
                        .map(Box::new)
                        .map_err(|error| query_error(input, region.to_string(), error))?,
                    None => reader
                        .query_unmapped(header)
                        .map(Box::new)
                        .map_err(|error| query_error(input, "*", error))?,
                },
                _ => {
                    return Err(RsomicsError::ConfigError(
                        "region query requires an indexed alignment reader".to_owned(),
                    ));
                }
            };

        for result in records {
            let record = result.map_err(|error| record_error(input, error))?;
            if !visit(record.as_ref())? {
                break;
            }
        }
        Ok(())
    }

    pub(crate) fn visit_owned_raw_region(
        &mut self,
        header: &sam::Header,
        input: &Path,
        region: &Region,
        mut visit: impl FnMut(RawRecord) -> Result<bool>,
    ) -> Result<()> {
        let mut encoder = RawRecordEncoder::new();
        self.visit_region(header, input, Some(region), |record| {
            visit(encoder.encode(header, record)?)
        })
    }
}

fn visit_raw_records(
    reader: &mut impl BufRead,
    input: &Path,
    mut visit: impl FnMut(RecordRef<'_>) -> Result<bool>,
) -> Result<()> {
    let mut records = RecordReader::new(reader);
    while let Some(record) = records
        .next()
        .rs_with_context(|| format!("reading alignment record from {}", input.display()))?
    {
        if !visit(record)? {
            break;
        }
    }
    Ok(())
}

fn visit_owned_raw_records(
    reader: &mut impl BufRead,
    input: &Path,
    mut visit: impl FnMut(RawRecord) -> Result<bool>,
) -> Result<()> {
    let mut records = RecordReader::new(reader);
    while let Some(record) = records
        .next()
        .rs_with_context(|| format!("reading alignment record from {}", input.display()))?
    {
        if !visit(RawRecord::try_from(record.payload().to_vec())?)? {
            break;
        }
    }
    Ok(())
}

fn visit_mut_raw_records(
    reader: &mut impl Read,
    input: &Path,
    mut visit: impl FnMut(&mut RawRecord) -> Result<bool>,
) -> Result<()> {
    let mut record = RawRecord::default();
    while raw::read_record(reader, &mut record)
        .rs_with_context(|| format!("reading alignment record from {}", input.display()))?
        != 0
    {
        if !visit(&mut record)? {
            break;
        }
    }
    Ok(())
}

pub(crate) fn open(
    input: &Path,
    reference: Option<&Path>,
    additional_threads: usize,
) -> Result<Reader> {
    if input == Path::new("-") {
        if additional_threads > 0 {
            return Err(RsomicsError::ConfigError(
                "additional decoding threads require a file-backed BAM input".to_owned(),
            ));
        }
        let mut stdin = io::stdin();
        let mut prefix = Vec::new();
        let (mut format, compression) = detect_stream(&mut stdin, &mut prefix)?;
        let source = Cursor::new(prefix).chain(stdin);
        let source: Box<dyn Read> = if compression == Compression::Gzip {
            let mut decoder = MultiGzDecoder::new(BufReader::new(source));
            let mut magic = [0; 4];
            decoder
                .read_exact(&mut magic)
                .map_err(|error| open_error(input, error))?;
            format = match magic {
                value if value == *b"CRAM" => Format::Cram,
                value if value == *b"BAM\x01" => Format::Bam,
                _ => Format::Sam,
            };
            Box::new(Cursor::new(magic).chain(decoder))
        } else {
            Box::new(source)
        };
        let mut builder = alignment::io::reader::Builder::default();
        if let Some(reference) = reference {
            builder = builder.set_reference_sequence_repository(reference_repository(reference)?);
        }
        let inner = builder
            .build_from_reader(source)
            .map_err(|error| open_error(input, error))?;
        return Ok(Reader {
            inner: Inner::General(inner),
            format,
        });
    }

    let (format, compression) = detect_source(input)?;
    if format == Format::Cram && additional_threads > 0 {
        return Err(RsomicsError::ConfigError(
            "additional CRAM decoding threads are not available yet".to_owned(),
        ));
    }

    let file = File::open(input).map_err(|error| open_error(input, error))?;
    let inner = match (format, compression) {
        (Format::Sam, Compression::None) => Inner::Sam(sam::io::Reader::new(BufReader::new(file))),
        (Format::Sam, Compression::Bgzf) => Inner::SamGz(sam::io::Reader::new(
            bgzf::io::Reader::new(BufReader::new(file)),
        )),
        (Format::Sam, Compression::Gzip) => Inner::SamGzip(Box::new(sam::io::Reader::new(
            BufReader::new(MultiGzDecoder::new(BufReader::new(file))),
        ))),
        (Format::Bam, Compression::Bgzf) => {
            let inner: Box<dyn BufRead + Send> =
                if let Some(workers) = NonZero::new(additional_threads) {
                    Box::new(bgzf::io::MultithreadedReader::with_worker_count(
                        workers, file,
                    ))
                } else {
                    Box::new(bgzf::io::Reader::new(BufReader::with_capacity(
                        READ_BUFFER,
                        file,
                    )))
                };
            Inner::Bam(bam::io::Reader::from(inner))
        }
        (Format::Bam, Compression::None) => {
            Inner::BamRaw(bam::io::Reader::from(BufReader::new(file)))
        }
        (Format::Bam, Compression::Gzip) => {
            return Err(RsomicsError::InvalidInput(format!(
                "BAM input must use BGZF compression: {}",
                input.display()
            )));
        }
        (Format::Cram, Compression::None) => {
            let mut builder = cram::io::reader::Builder::default();
            if let Some(reference) = reference {
                builder =
                    builder.set_reference_sequence_repository(reference_repository(reference)?);
            }
            Inner::Cram(builder.build_from_reader(BufReader::new(file)))
        }
        (Format::Cram, Compression::Bgzf) => {
            return Err(RsomicsError::InvalidInput(format!(
                "CRAM input cannot be BGZF-compressed: {}",
                input.display()
            )));
        }
        (Format::Cram, Compression::Gzip) => {
            return Err(RsomicsError::InvalidInput(format!(
                "CRAM input cannot be gzip-compressed: {}",
                input.display()
            )));
        }
    };

    Ok(Reader { inner, format })
}

pub(crate) fn open_indexed(input: &Path, reference: Option<&Path>) -> Result<Reader> {
    let format = detect_format(input)?;
    let inner = open_indexed_alignment(input, reference).map(Inner::Indexed)?;
    Ok(Reader { inner, format })
}

pub(crate) fn open_indexed_with_index(
    input: &Path,
    index: &Path,
    reference: Option<&Path>,
) -> Result<Reader> {
    let format = detect_format(input)?;
    let mut builder = alignment::io::indexed_reader::Builder::default();
    if let Some(reference) = reference {
        builder = builder.set_reference_sequence_repository(reference_repository(reference)?);
    }
    builder = match format {
        Format::Sam | Format::Bam => match bam::bai::fs::read(index) {
            Ok(index) => builder.set_index(index),
            Err(bai_error) => match noodles::csi::fs::read(index) {
                Ok(index) => builder.set_index(index),
                Err(csi_error) => {
                    return Err(RsomicsError::InvalidInput(format!(
                        "reading custom alignment index {}: BAI: {bai_error}; CSI: {csi_error}",
                        index.display()
                    )));
                }
            },
        },
        Format::Cram => {
            let index = cram::crai::fs::read(index).map_err(|error| {
                RsomicsError::InvalidInput(format!(
                    "reading custom CRAM index {}: {error}",
                    index.display()
                ))
            })?;
            builder.set_index(index)
        }
    };
    let inner = builder
        .build_from_path(input)
        .map(Inner::IndexedDirect)
        .map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "opening {} with custom index {}: {error}",
                input.display(),
                index.display()
            ))
        })?;
    Ok(Reader { inner, format })
}

pub(crate) fn detect_format(input: &Path) -> Result<Format> {
    detect_source(input).map(|(format, _)| format)
}

pub(crate) fn is_bgzf(input: &Path) -> Result<bool> {
    detect_source(input).map(|(_, compression)| compression == Compression::Bgzf)
}

fn detect_source(input: &Path) -> Result<(Format, Compression)> {
    let mut source = BufReader::new(File::open(input).map_err(|error| open_error(input, error))?);
    let mut magic = [0; 4];
    source.read_exact(&mut magic).map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "detecting alignment format for {}: {error}",
            input.display()
        ))
    })?;

    if magic == *b"CRAM" {
        return Ok((Format::Cram, Compression::None));
    }
    if magic == *b"BAM\x01" {
        return Ok((Format::Bam, Compression::None));
    }
    if magic[..2] != [0x1f, 0x8b] {
        return Ok((Format::Sam, Compression::None));
    }

    let file = File::open(input).map_err(|error| open_error(input, error))?;
    let mut reader = bgzf::io::Reader::new(file);
    if reader.read_exact(&mut magic).is_ok() {
        return Ok((
            if magic == *b"BAM\x01" {
                Format::Bam
            } else {
                Format::Sam
            },
            Compression::Bgzf,
        ));
    }

    let file = File::open(input).map_err(|error| open_error(input, error))?;
    let mut reader = MultiGzDecoder::new(BufReader::new(file));
    reader.read_exact(&mut magic).map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "detecting gzip alignment format for {}: {error}",
            input.display()
        ))
    })?;
    Ok((
        if magic == *b"BAM\x01" {
            Format::Bam
        } else if magic == *b"CRAM" {
            Format::Cram
        } else {
            Format::Sam
        },
        Compression::Gzip,
    ))
}

fn detect_stream(source: &mut impl Read, prefix: &mut Vec<u8>) -> Result<(Format, Compression)> {
    let mut magic = [0; 4];
    source
        .read_exact(&mut magic)
        .map_err(|error| open_error(Path::new("-"), error))?;
    prefix.extend_from_slice(&magic);

    if magic == *b"CRAM" {
        return Ok((Format::Cram, Compression::None));
    }
    if magic == *b"BAM\x01" {
        return Ok((Format::Bam, Compression::None));
    }
    if magic[..2] != [0x1f, 0x8b] {
        return Ok((Format::Sam, Compression::None));
    }

    let mut gzip_header = [0; 8];
    source
        .read_exact(&mut gzip_header)
        .map_err(|error| open_error(Path::new("-"), error))?;
    prefix.extend_from_slice(&gzip_header);
    if prefix[3] & 0x04 == 0 {
        return Ok((Format::Sam, Compression::Gzip));
    }

    let extra_length = usize::from(u16::from_le_bytes([prefix[10], prefix[11]]));
    let mut extra = vec![0; extra_length];
    source
        .read_exact(&mut extra)
        .map_err(|error| open_error(Path::new("-"), error))?;
    prefix.extend_from_slice(&extra);

    let mut offset = 0;
    let mut block_size = None;
    while offset + 4 <= extra.len() {
        let length = usize::from(u16::from_le_bytes([extra[offset + 2], extra[offset + 3]]));
        let end = offset
            .checked_add(4)
            .and_then(|value| value.checked_add(length))
            .filter(|end| *end <= extra.len())
            .ok_or_else(|| {
                RsomicsError::InvalidInput("invalid BGZF extra field on standard input".to_owned())
            })?;
        if extra[offset..offset + 2] == *b"BC" && length == 2 {
            block_size =
                Some(usize::from(u16::from_le_bytes([extra[offset + 4], extra[offset + 5]])) + 1);
            break;
        }
        offset = end;
    }
    let Some(block_size) = block_size else {
        return Ok((Format::Sam, Compression::Gzip));
    };
    let remaining = block_size.checked_sub(prefix.len()).ok_or_else(|| {
        RsomicsError::InvalidInput("invalid BGZF block size on standard input".to_owned())
    })?;
    let start = prefix.len();
    prefix.resize(block_size, 0);
    source
        .read_exact(&mut prefix[start..start + remaining])
        .map_err(|error| open_error(Path::new("-"), error))?;

    let mut reader = bgzf::io::Reader::new(Cursor::new(prefix.as_slice()));
    reader.read_exact(&mut magic).map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "detecting BGZF alignment format on standard input: {error}"
        ))
    })?;
    Ok((
        if magic == *b"BAM\x01" {
            Format::Bam
        } else {
            Format::Sam
        },
        Compression::Bgzf,
    ))
}

pub(crate) fn reference_repository(path: &Path) -> Result<fasta::Repository> {
    fasta::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .map(fasta::repository::adapters::IndexedReader::new)
        .map(fasta::Repository::new)
        .map_err(|error| {
            RsomicsError::ConfigError(format!(
                "opening indexed reference {}: {error}",
                path.display()
            ))
        })
}

fn open_error(input: &Path, error: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "opening alignment input {}: {error}",
        input.display()
    ))
}

fn record_error(input: &Path, error: io::Error) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "reading alignment record from {}: {error}",
        input.display()
    ))
}

fn query_error(input: &Path, region: impl std::fmt::Display, error: io::Error) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "querying region {region} from {}: {error}",
        input.display()
    ))
}
