use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read};
use std::num::NonZero;
use std::path::{Path, PathBuf};

use noodles::{bam, bgzf, core::Region, cram, csi, fasta, sam};
use noodles_util::alignment;
use rsomics_bamio::raw::{RecordReader, RecordRef};
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
}

enum Inner {
    General(alignment::io::Reader<Box<dyn Read>>),
    Sam(sam::io::Reader<BufReader<File>>),
    SamGz(sam::io::Reader<bgzf::io::Reader<BufReader<File>>>),
    Bam(ParallelBamReader),
    BamRaw(bam::io::Reader<BufReader<File>>),
    Cram(cram::io::Reader<BufReader<File>>),
    Indexed(alignment::io::IndexedReader<File>),
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
            Inner::Indexed(reader) => reader.read_header(),
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
            Inner::Indexed(reader) => {
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

    pub(crate) fn visit_region(
        &mut self,
        header: &sam::Header,
        input: &Path,
        region: Option<&Region>,
        mut visit: impl FnMut(&dyn sam::alignment::Record) -> Result<bool>,
    ) -> Result<()> {
        let Inner::Indexed(reader) = &mut self.inner else {
            return Err(RsomicsError::ConfigError(
                "region query requires an indexed alignment reader".to_owned(),
            ));
        };

        let records: Box<dyn Iterator<Item = io::Result<Box<dyn sam::alignment::Record>>> + '_> =
            match region {
                Some(region) => reader
                    .query(header, region)
                    .map(Box::new)
                    .map_err(|error| query_error(input, region.to_string(), error))?,
                None => reader
                    .query_unmapped(header)
                    .map(Box::new)
                    .map_err(|error| query_error(input, "*", error))?,
            };

        for result in records {
            let record = result.map_err(|error| record_error(input, error))?;
            if !visit(record.as_ref())? {
                break;
            }
        }
        Ok(())
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
            let inner: Box<dyn BufRead + Send> =
                if let Some(workers) = additional_threads.checked_sub(1).and_then(NonZero::new) {
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

pub(crate) fn open_indexed(input: &Path, reference: Option<&Path>) -> Result<Reader> {
    if input == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "region queries require a file-backed alignment input".to_owned(),
        ));
    }

    let (format, compression) = detect_source(input)?;
    if matches!(
        (format, compression),
        (Format::Sam | Format::Bam, Compression::None)
    ) {
        return Err(RsomicsError::InvalidInput(format!(
            "region queries require BGZF SAM, BAM, or CRAM input: {}",
            input.display()
        )));
    }

    let mut builder = alignment::io::indexed_reader::Builder::default();
    if let Some(reference) = reference {
        builder = builder.set_reference_sequence_repository(reference_repository(reference)?);
    }
    builder = set_alternative_index(builder, input, format)?;

    let inner = builder
        .build_from_path(input)
        .map(Inner::Indexed)
        .map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "opening indexed alignment input {}: {error}",
                input.display()
            ))
        })?;
    Ok(Reader { inner, format })
}

pub(crate) fn detect_format(input: &Path) -> Result<Format> {
    detect_source(input).map(|(format, _)| format)
}

fn set_alternative_index(
    mut builder: alignment::io::indexed_reader::Builder,
    input: &Path,
    format: Format,
) -> Result<alignment::io::indexed_reader::Builder> {
    match format {
        Format::Sam => {
            let appended = append_extension(input, "csi");
            let alternative = input.with_extension("csi");
            if !index_exists(input, &appended)? && index_exists(input, &alternative)? {
                let index = csi::fs::read(&alternative)
                    .map_err(|error| index_error(input, &alternative, error))?;
                builder = builder.set_index(index);
            }
        }
        Format::Bam => {
            let appended_bai = append_extension(input, "bai");
            let appended_csi = append_extension(input, "csi");
            if !index_exists(input, &appended_bai)? && !index_exists(input, &appended_csi)? {
                let alternative_bai = input.with_extension("bai");
                let alternative_csi = input.with_extension("csi");
                if index_exists(input, &alternative_bai)? {
                    let index = bam::bai::fs::read(&alternative_bai)
                        .map_err(|error| index_error(input, &alternative_bai, error))?;
                    builder = builder.set_index(index);
                } else if index_exists(input, &alternative_csi)? {
                    let index = csi::fs::read(&alternative_csi)
                        .map_err(|error| index_error(input, &alternative_csi, error))?;
                    builder = builder.set_index(index);
                }
            }
        }
        Format::Cram => {
            let appended = append_extension(input, "crai");
            let alternative = input.with_extension("crai");
            if !index_exists(input, &appended)? && index_exists(input, &alternative)? {
                let index = cram::crai::fs::read(&alternative)
                    .map_err(|error| index_error(input, &alternative, error))?;
                builder = builder.set_index(index);
            }
        }
    }
    Ok(builder)
}

fn append_extension(path: &Path, extension: &str) -> PathBuf {
    let mut path = path.as_os_str().to_owned();
    path.push(".");
    path.push(extension);
    PathBuf::from(path)
}

fn index_exists(input: &Path, index: &Path) -> Result<bool> {
    index
        .try_exists()
        .map_err(|error| index_error(input, index, error))
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

fn query_error(input: &Path, region: impl std::fmt::Display, error: io::Error) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "querying region {region} from {}: {error}",
        input.display()
    ))
}

fn index_error(input: &Path, index: &Path, error: io::Error) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "reading alignment index {} for {}: {error}",
        index.display(),
        input.display()
    ))
}
