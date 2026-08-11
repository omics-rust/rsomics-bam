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

pub(crate) struct TransactionalFile {
    target: PathBuf,
    temporary: tempfile::NamedTempFile,
    permissions: Option<fs::Permissions>,
}

impl TransactionalFile {
    pub(crate) fn new(target: &Path) -> Result<Self> {
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
        let metadata = fs::metadata(target).ok();
        if metadata
            .as_ref()
            .is_some_and(|metadata| !metadata.is_file())
        {
            return Err(RsomicsError::ConfigError(format!(
                "output target is not a regular file: {}",
                target.display()
            )));
        }
        let permissions = metadata.map(|metadata| metadata.permissions());
        Ok(Self {
            target: target.to_owned(),
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
        self.apply_permissions()?;
        self.temporary.persist(&self.target).map_err(|error| {
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

    pub(crate) fn commit_all(mut transactions: Vec<Self>) -> Result<()> {
        for transaction in &mut transactions {
            transaction.apply_permissions()?;
        }

        let mut backups = Vec::new();
        for transaction in &transactions {
            if transaction.target.exists() {
                let backup = reserve_backup_path(&transaction.target)?;
                if let Err(error) = fs::rename(&transaction.target, &backup) {
                    return restore_backups(backups, error);
                }
                backups.push((transaction.target.clone(), backup));
            }
        }

        let mut committed = Vec::new();
        for transaction in transactions {
            let target = transaction.target;
            if let Err(error) = transaction.temporary.persist(&target) {
                return rollback_outputs(committed, backups, error.error);
            }
            committed.push(target);
        }
        drop(backups);
        Ok(())
    }

    fn apply_permissions(&mut self) -> Result<()> {
        if let Some(permissions) = self.permissions.take() {
            self.temporary
                .as_file_mut()
                .set_permissions(permissions)
                .map_err(RsomicsError::Io)?;
        }
        Ok(())
    }
}

fn reserve_backup_path(target: &Path) -> Result<tempfile::TempPath> {
    let parent = target
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let backup = tempfile::NamedTempFile::new_in(parent)
        .map_err(RsomicsError::Io)?
        .into_temp_path();
    fs::remove_file(&backup).map_err(RsomicsError::Io)?;
    Ok(backup)
}

fn rollback_outputs(
    committed: Vec<PathBuf>,
    backups: Vec<(PathBuf, tempfile::TempPath)>,
    cause: io::Error,
) -> Result<()> {
    let mut rollback_error = None;
    for target in committed {
        if let Err(error) = fs::remove_file(&target)
            && error.kind() != io::ErrorKind::NotFound
        {
            rollback_error.get_or_insert(error);
        }
    }
    for (target, backup) in backups {
        if let Err(error) = fs::rename(&backup, &target) {
            rollback_error.get_or_insert(error);
        }
    }
    output_failure(cause, rollback_error)
}

fn restore_backups(backups: Vec<(PathBuf, tempfile::TempPath)>, cause: io::Error) -> Result<()> {
    let mut rollback_error = None;
    for (target, backup) in backups {
        if let Err(error) = fs::rename(&backup, &target) {
            rollback_error.get_or_insert(error);
        }
    }
    output_failure(cause, rollback_error)
}

fn output_failure(cause: io::Error, rollback_error: Option<io::Error>) -> Result<()> {
    let error = if let Some(rollback_error) = rollback_error {
        io::Error::new(
            cause.kind(),
            format!(
                "committing output failed: {cause}; restoring prior outputs failed: {rollback_error}"
            ),
        )
    } else {
        cause
    };
    Err(RsomicsError::Io(error))
}

pub(crate) fn same_target(left: &Path, right: &Path) -> Result<bool> {
    if left.exists() && right.exists() {
        return same_file::is_same_file(left, right).map_err(RsomicsError::Io);
    }
    Ok(target_identity(left)? == target_identity(right)?)
}

pub(crate) fn target_identity(path: &Path) -> Result<PathBuf> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transaction_owns_its_target() {
        let directory = tempfile::tempdir().unwrap();
        let mut transaction = {
            let target = directory.path().join("owned");
            TransactionalFile::new(&target).unwrap()
        };
        transaction.file_mut().write_all(b"complete").unwrap();
        transaction.commit().unwrap();
        assert_eq!(
            fs::read(directory.path().join("owned")).unwrap(),
            b"complete"
        );
    }

    #[test]
    fn grouped_commit_restores_every_prior_target() {
        let directory = tempfile::tempdir().unwrap();
        let first_path = directory.path().join("first");
        let second_path = directory.path().join("second");
        fs::write(&first_path, b"old first").unwrap();
        fs::write(&second_path, b"old second").unwrap();

        let mut first = TransactionalFile::new(&first_path).unwrap();
        first.file_mut().write_all(b"new first").unwrap();
        let mut second = TransactionalFile::new(&second_path).unwrap();
        second.file_mut().write_all(b"new second").unwrap();
        fs::remove_file(second.temporary_path()).unwrap();

        assert!(TransactionalFile::commit_all(vec![first, second]).is_err());
        assert_eq!(fs::read(first_path).unwrap(), b"old first");
        assert_eq!(fs::read(second_path).unwrap(), b"old second");
    }

    #[test]
    fn grouped_commit_removes_new_outputs_during_rollback() {
        let directory = tempfile::tempdir().unwrap();
        let first_path = directory.path().join("first");
        let second_path = directory.path().join("second");

        let mut first = TransactionalFile::new(&first_path).unwrap();
        first.file_mut().write_all(b"new first").unwrap();
        let second = TransactionalFile::new(&second_path).unwrap();
        fs::remove_file(second.temporary_path()).unwrap();

        assert!(TransactionalFile::commit_all(vec![first, second]).is_err());
        assert!(!first_path.exists());
        assert!(!second_path.exists());
    }
}
