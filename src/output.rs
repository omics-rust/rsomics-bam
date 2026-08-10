use std::fs;
use std::io::{self, BufWriter, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};

use noodles::sam::alignment::io::Write as _;
use noodles::{bam, bgzf, sam};
use rsomics_bamio::RingBgzfWriter;
use rsomics_bamio::raw::{self, RawRecord, RecordRef};
use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum Format {
    #[default]
    Sam,
    Bam,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum Compression {
    #[default]
    Default,
    Fast,
    Uncompressed,
}

pub(crate) struct Writer<W>
where
    W: Write + Send + 'static,
{
    inner: Inner<W>,
}

enum Inner<W>
where
    W: Write + Send + 'static,
{
    Sam(sam::io::Writer<BufWriter<W>>),
    BamSingle(bam::io::Writer<bgzf::io::Writer<W>>),
    BamParallel(bam::io::Writer<RingBgzfWriter<W>>),
    BamParallelLevel(bam::io::Writer<bgzf::io::MultithreadedWriter<W>>),
}

pub(crate) struct TransactionalFile<'a> {
    target: &'a Path,
    temporary: tempfile::NamedTempFile,
    permissions: Option<fs::Permissions>,
}

impl<'a> TransactionalFile<'a> {
    pub(crate) fn new(target: &'a Path) -> Result<Self> {
        let parent = target
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        let temporary = tempfile::NamedTempFile::new_in(parent).map_err(|error| {
            RsomicsError::Io(io::Error::new(
                error.kind(),
                format!(
                    "creating temporary output beside {}: {error}",
                    target.display()
                ),
            ))
        })?;
        let permissions = fs::metadata(target)
            .ok()
            .map(|metadata| metadata.permissions());
        Ok(Self {
            target,
            temporary,
            permissions,
        })
    }

    pub(crate) fn reopen(&self) -> Result<fs::File> {
        self.temporary.reopen().map_err(|error| {
            RsomicsError::Io(io::Error::new(
                error.kind(),
                format!(
                    "opening temporary output beside {}: {error}",
                    self.target.display()
                ),
            ))
        })
    }

    pub(crate) fn file_mut(&mut self) -> &mut fs::File {
        self.temporary.as_file_mut()
    }

    pub(crate) fn temporary_path(&self) -> &Path {
        self.temporary.path()
    }

    pub(crate) fn commit(mut self) -> Result<()> {
        if let Some(permissions) = self.permissions {
            self.temporary
                .as_file_mut()
                .set_permissions(permissions)
                .map_err(RsomicsError::Io)?;
        }
        self.temporary.persist(self.target).map_err(|error| {
            RsomicsError::Io(io::Error::new(
                error.error.kind(),
                format!(
                    "committing output {}: {}",
                    self.target.display(),
                    error.error
                ),
            ))
        })?;
        Ok(())
    }
}

pub(crate) fn same_target(left: &Path, right: &Path) -> Result<bool> {
    if left.exists() && right.exists() {
        return same_file::is_same_file(left, right).map_err(RsomicsError::Io);
    }
    Ok(target_identity(left)? == target_identity(right)?)
}

fn target_identity(path: &Path) -> Result<PathBuf> {
    match fs::canonicalize(path) {
        Ok(path) => return Ok(path),
        Err(error) if error.kind() == io::ErrorKind::NotFound => {}
        Err(error) => return Err(RsomicsError::Io(error)),
    }
    let name = path.file_name().ok_or_else(|| {
        RsomicsError::ConfigError(format!("output path has no file name: {}", path.display()))
    })?;
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    fs::canonicalize(parent)
        .map(|parent| parent.join(name))
        .map_err(RsomicsError::Io)
}

impl<W> Writer<W>
where
    W: Write + Send + 'static,
{
    pub(crate) fn new(
        format: Format,
        compression: Compression,
        additional_threads: usize,
        output: W,
    ) -> Self {
        let inner = match format {
            Format::Sam => Inner::Sam(sam::io::Writer::new(BufWriter::new(output))),
            Format::Bam => {
                let level = match compression {
                    Compression::Default => bgzf::io::writer::CompressionLevel::default(),
                    Compression::Fast => bgzf::io::writer::CompressionLevel::FAST,
                    Compression::Uncompressed => bgzf::io::writer::CompressionLevel::NONE,
                };
                if let Some(workers) = NonZero::new(additional_threads) {
                    if compression == Compression::Default {
                        Inner::BamParallel(bam::io::Writer::from(RingBgzfWriter::new(
                            output, workers,
                        )))
                    } else {
                        let writer = bgzf::io::multithreaded_writer::Builder::default()
                            .set_compression_level(level)
                            .set_worker_count(workers)
                            .build_from_writer(output);
                        Inner::BamParallelLevel(bam::io::Writer::from(writer))
                    }
                } else {
                    let writer = bgzf::io::writer::Builder::default()
                        .set_compression_level(level)
                        .build_from_writer(output);
                    Inner::BamSingle(bam::io::Writer::from(writer))
                }
            }
        };
        Self { inner }
    }

    pub(crate) fn write_header(&mut self, header: &sam::Header) -> Result<()> {
        match &mut self.inner {
            Inner::Sam(writer) => writer.write_header(header),
            Inner::BamSingle(writer) => writer.write_header(header),
            Inner::BamParallel(writer) => writer.write_header(header),
            Inner::BamParallelLevel(writer) => writer.write_header(header),
        }
        .map_err(RsomicsError::Io)
    }

    pub(crate) fn write_record(
        &mut self,
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
    ) -> Result<()> {
        match &mut self.inner {
            Inner::Sam(writer) => writer.write_alignment_record(header, record),
            Inner::BamSingle(writer) => writer.write_alignment_record(header, record),
            Inner::BamParallel(writer) => writer.write_alignment_record(header, record),
            Inner::BamParallelLevel(writer) => writer.write_alignment_record(header, record),
        }
        .map_err(RsomicsError::Io)
    }

    pub(crate) fn write_raw_record(&mut self, record: &RecordRef<'_>) -> Result<()> {
        match &mut self.inner {
            Inner::Sam(_) => Err(RsomicsError::ConfigError(
                "raw BAM records require BAM output".to_owned(),
            )),
            Inner::BamSingle(writer) => raw::write_record_ref(writer.get_mut(), record),
            Inner::BamParallel(writer) => raw::write_record_ref(writer.get_mut(), record),
            Inner::BamParallelLevel(writer) => raw::write_record_ref(writer.get_mut(), record),
        }
    }

    pub(crate) fn write_owned_raw_record(&mut self, record: &RawRecord) -> Result<()> {
        match &mut self.inner {
            Inner::Sam(_) => Err(RsomicsError::ConfigError(
                "raw BAM records require BAM output".to_owned(),
            )),
            Inner::BamSingle(writer) => raw::write_record(writer.get_mut(), record),
            Inner::BamParallel(writer) => raw::write_record(writer.get_mut(), record),
            Inner::BamParallelLevel(writer) => raw::write_record(writer.get_mut(), record),
        }
    }

    pub(crate) fn finish(self, header: &sam::Header) -> Result<()> {
        match self.inner {
            Inner::Sam(mut writer) => {
                writer.finish(header)?;
                writer.get_mut().flush()
            }
            Inner::BamSingle(writer) => writer.into_inner().finish().map(drop),
            Inner::BamParallel(writer) => writer.into_inner().finish().map(drop),
            Inner::BamParallelLevel(mut writer) => writer.get_mut().finish().map(drop),
        }
        .map_err(RsomicsError::Io)
    }
}
