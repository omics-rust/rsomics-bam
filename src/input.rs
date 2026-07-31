use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read};
use std::num::NonZero;
use std::path::Path;

use noodles::{bam, bgzf, cram, fasta, sam};
use noodles_util::alignment;
use rsomics_common::{Result, RsomicsError};

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
}

enum Inner {
    General(alignment::io::Reader<Box<dyn Read>>),
    Sam(sam::io::Reader<BufReader<File>>),
    SamGz(sam::io::Reader<bgzf::io::Reader<BufReader<File>>>),
    Bam(ParallelBamReader),
    BamRaw(bam::io::Reader<BufReader<File>>),
    Cram(cram::io::Reader<BufReader<File>>),
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

    pub(crate) fn read_header(&mut self, input: &Path) -> Result<sam::Header> {
        let result = match &mut self.inner {
            Inner::General(reader) => reader.read_header(),
            Inner::Sam(reader) => reader.read_header(),
            Inner::SamGz(reader) => reader.read_header(),
            Inner::Bam(reader) => reader.read_header(),
            Inner::BamRaw(reader) => reader.read_header(),
            Inner::Cram(reader) => reader.read_header(),
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
        }
        Ok(())
    }
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
        let format = detect_stream(&mut stdin, &mut prefix)?;
        let source = Cursor::new(prefix).chain(stdin);
        let mut builder = alignment::io::reader::Builder::default();
        if let Some(reference) = reference {
            builder = builder.set_reference_sequence_repository(reference_repository(reference)?);
        }
        let inner = builder
            .build_from_reader(Box::new(source) as Box<dyn Read>)
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
        (Format::Bam, Compression::Bgzf) => {
            let workers = additional_threads
                .checked_add(1)
                .and_then(NonZero::new)
                .ok_or_else(|| {
                    RsomicsError::ConfigError(
                        "alignment thread count exceeds the supported range".to_owned(),
                    )
                })?;
            let inner: Box<dyn BufRead + Send> = if workers.get() == 1 {
                Box::new(bgzf::io::Reader::new(BufReader::with_capacity(
                    READ_BUFFER,
                    file,
                )))
            } else {
                Box::new(bgzf::io::MultithreadedReader::with_worker_count(
                    workers, file,
                ))
            };
            Inner::Bam(bam::io::Reader::from(inner))
        }
        (Format::Bam, Compression::None) => {
            Inner::BamRaw(bam::io::Reader::from(BufReader::new(file)))
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
    };

    Ok(Reader { inner, format })
}

pub(crate) fn detect_format(input: &Path) -> Result<Format> {
    detect_source(input).map(|(format, _)| format)
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
    reader.read_exact(&mut magic).map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "detecting BGZF alignment format for {}: {error}",
            input.display()
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

fn detect_stream(source: &mut impl Read, prefix: &mut Vec<u8>) -> Result<Format> {
    let mut magic = [0; 4];
    source
        .read_exact(&mut magic)
        .map_err(|error| open_error(Path::new("-"), error))?;
    prefix.extend_from_slice(&magic);

    if magic == *b"CRAM" {
        return Ok(Format::Cram);
    }
    if magic == *b"BAM\x01" {
        return Ok(Format::Bam);
    }
    if magic[..2] != [0x1f, 0x8b] {
        return Ok(Format::Sam);
    }

    let mut gzip_header = [0; 8];
    source
        .read_exact(&mut gzip_header)
        .map_err(|error| open_error(Path::new("-"), error))?;
    prefix.extend_from_slice(&gzip_header);
    if prefix[3] & 0x04 == 0 {
        return Err(RsomicsError::InvalidInput(
            "gzip-compressed alignment input is not BGZF".to_owned(),
        ));
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
    let block_size = block_size.ok_or_else(|| {
        RsomicsError::InvalidInput("missing BGZF block size on standard input".to_owned())
    })?;
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
    Ok(if magic == *b"BAM\x01" {
        Format::Bam
    } else {
        Format::Sam
    })
}

fn reference_repository(path: &Path) -> Result<fasta::Repository> {
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
