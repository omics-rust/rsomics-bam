use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::path::{Path, PathBuf};

use noodles::{fasta, sam};
use noodles_util::alignment;
use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};

use super::{Format, Mode, Options, OutputSummary, Summary, label, tag};
use crate::output::{self, TransactionalFile};

pub(super) struct Router {
    protected: Vec<PathBuf>,
    header: sam::Header,
    format: Format,
    prefix: PathBuf,
    maximum_outputs: usize,
    fixed_labels: bool,
    destinations: HashMap<Vec<u8>, usize>,
    sinks: Vec<Sink>,
    unaccounted: Option<usize>,
    repository: Option<fasta::Repository>,
}

struct Sink {
    label: String,
    path: PathBuf,
    header: sam::Header,
    writer: Writer,
    transaction: TransactionalFile,
    records: u64,
}

enum Writer {
    SamBam(output::Writer<File>),
    Cram(alignment::io::Writer<File>),
}

impl Writer {
    fn write_raw(&mut self, record: &RawRecord) -> Result<()> {
        match self {
            Self::SamBam(writer) => writer.write_owned_raw_record(record),
            Self::Cram(_) => Err(RsomicsError::ConfigError(
                "raw records cannot be written directly to CRAM".to_owned(),
            )),
        }
    }

    fn write_record(
        &mut self,
        header: &sam::Header,
        record: &sam::alignment::RecordBuf,
    ) -> Result<()> {
        match self {
            Self::SamBam(writer) => writer.write_record(header, record),
            Self::Cram(writer) => writer
                .write_record(header, record)
                .map_err(RsomicsError::Io),
        }
    }

    fn finish(self, header: &sam::Header) -> Result<()> {
        match self {
            Self::SamBam(writer) => writer.finish(header),
            Self::Cram(mut writer) => writer.finish(header).map_err(RsomicsError::Io),
        }
    }
}

impl Router {
    pub(super) fn new(
        input: &Path,
        header: &sam::Header,
        unaccounted_header: Option<sam::Header>,
        options: Options<'_>,
    ) -> Result<Self> {
        let repository = if options.format == Format::Cram {
            Some(crate::input::reference_repository(
                options.reference.ok_or_else(|| {
                    RsomicsError::ConfigError(
                        "CRAM split output requires an indexed reference".to_owned(),
                    )
                })?,
            )?)
        } else {
            None
        };
        let mut protected = vec![input.to_owned()];
        if let Some(path) = options.reference {
            protected.push(path.to_owned());
            protected.push(with_suffix(path, ".fai"));
        }
        if let Some(path) = options.unaccounted_header {
            protected.push(path.to_owned());
        }
        let mut router = Self {
            protected,
            header: header.clone(),
            format: options.format,
            prefix: options.output_prefix.to_owned(),
            maximum_outputs: options.maximum_outputs,
            fixed_labels: options.mode == Mode::ReadGroup || options.mode == Mode::Tag(*b"RG"),
            destinations: HashMap::new(),
            sinks: Vec::new(),
            unaccounted: None,
            repository,
        };

        if options.mode == Mode::ReadGroup || options.mode == Mode::Tag(*b"RG") {
            for id in header.read_groups().keys() {
                router.add_destination(id.as_ref())?;
            }
        }
        if let Some(path) = options.unaccounted {
            let index = router.add_sink(
                "unaccounted".to_owned(),
                path.to_owned(),
                unaccounted_header.unwrap_or_else(|| header.clone()),
            )?;
            router.unaccounted = Some(index);
        }
        Ok(router)
    }

    pub(super) fn write_raw(&mut self, record: &RawRecord, outcome: tag::Outcome) -> Result<()> {
        let destination = self.destination(record.name(), outcome)?;
        let sink = &mut self.sinks[destination];
        sink.writer.write_raw(record)?;
        sink.records = increment(sink.records)?;
        Ok(())
    }

    pub(super) fn write_record(
        &mut self,
        record: &sam::alignment::RecordBuf,
        outcome: tag::Outcome,
    ) -> Result<()> {
        let name = record.name().map_or(b"".as_slice(), |name| name.as_ref());
        let destination = self.destination(name, outcome)?;
        let sink = &mut self.sinks[destination];
        sink.writer.write_record(&sink.header, record)?;
        sink.records = increment(sink.records)?;
        Ok(())
    }

    fn destination(&mut self, name: &[u8], outcome: tag::Outcome) -> Result<usize> {
        let destination = match outcome {
            tag::Outcome::Present(value) => {
                if let Some(&index) = self.destinations.get(value.as_slice()) {
                    index
                } else if self.fixed_labels {
                    self.unaccounted.ok_or_else(|| {
                        RsomicsError::InvalidInput(format!(
                            "read {} has unknown read group {}",
                            String::from_utf8_lossy(name),
                            String::from_utf8_lossy(&value)
                        ))
                    })?
                } else if self.destinations.len() < self.maximum_outputs {
                    self.add_destination(&value)?
                } else if let Some(index) = self.unaccounted {
                    index
                } else {
                    return Err(RsomicsError::InvalidInput(format!(
                        "split output count exceeds {} at label {}",
                        self.maximum_outputs,
                        String::from_utf8_lossy(&value)
                    )));
                }
            }
            tag::Outcome::Missing | tag::Outcome::Invalid => self.unaccounted.ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "read {} has a missing or incompatible split tag",
                    String::from_utf8_lossy(name)
                ))
            })?,
        };
        Ok(destination)
    }

    pub(super) fn finish(self) -> Result<Summary> {
        let mut transactions = Vec::with_capacity(self.sinks.len());
        let mut outputs = Vec::with_capacity(self.sinks.len());
        for sink in self.sinks {
            sink.writer.finish(&sink.header)?;
            transactions.push(sink.transaction);
            outputs.push(OutputSummary {
                label: sink.label,
                path: sink.path,
                records: sink.records,
            });
        }
        TransactionalFile::commit_all(transactions)?;
        Ok(Summary {
            records: 0,
            outputs,
            skipped: 0,
        })
    }

    fn add_destination(&mut self, value: &[u8]) -> Result<usize> {
        if value.is_empty() {
            return self.unaccounted.ok_or_else(|| {
                RsomicsError::InvalidInput("split tag value cannot be empty".to_owned())
            });
        }
        if let Some(&index) = self.destinations.get(value) {
            return Ok(index);
        }
        let encoded = label::encode(value)?;
        let path = output_path(&self.prefix, &encoded, self.format);
        let mut header = self.header.clone();
        if self.fixed_labels && self.header.read_groups().contains_key(value) {
            header
                .read_groups_mut()
                .retain(|id, _| id.as_slice() == value);
        } else if self.fixed_labels {
            header.read_groups_mut().clear();
        }
        let index = self.add_sink(encoded, path, header)?;
        self.destinations.insert(value.to_vec(), index);
        Ok(index)
    }

    fn add_sink(&mut self, label: String, path: PathBuf, header: sam::Header) -> Result<usize> {
        for input in &self.protected {
            if output::same_target(input, &path)? {
                return Err(RsomicsError::ConfigError(format!(
                    "split output aliases an input: {}",
                    path.display()
                )));
            }
        }
        for sink in &self.sinks {
            if output::same_target(&sink.path, &path)? {
                return Err(RsomicsError::ConfigError(format!(
                    "split outputs alias each other: {}",
                    path.display()
                )));
            }
        }
        let transaction = TransactionalFile::new(&path)?;
        let file = transaction.reopen()?;
        let writer = match self.format {
            Format::Sam | Format::Bam => {
                let format = match self.format {
                    Format::Sam => output::Format::Sam,
                    Format::Bam => output::Format::Bam,
                    Format::Cram => unreachable!(),
                };
                let mut writer = output::Writer::new(format, output::Compression::Default, 0, file);
                writer.write_header(&header)?;
                Writer::SamBam(writer)
            }
            Format::Cram => {
                let repository = self.repository.clone().ok_or_else(|| {
                    RsomicsError::ConfigError(
                        "CRAM split output requires an indexed reference".to_owned(),
                    )
                })?;
                let mut writer = alignment::io::writer::Builder::default()
                    .set_format(alignment::io::Format::Cram)
                    .set_reference_sequence_repository(repository)
                    .build_from_writer(file)
                    .map_err(RsomicsError::Io)?;
                writer.write_header(&header).map_err(RsomicsError::Io)?;
                Writer::Cram(writer)
            }
        };
        let index = self.sinks.len();
        self.sinks.push(Sink {
            label,
            path,
            header,
            writer,
            transaction,
            records: 0,
        });
        Ok(index)
    }
}

fn increment(value: u64) -> Result<u64> {
    value
        .checked_add(1)
        .ok_or_else(|| RsomicsError::InvalidInput("split output count exceeds u64".to_owned()))
}

fn output_path(prefix: &Path, label: &str, format: Format) -> PathBuf {
    let mut path = OsString::from(prefix.as_os_str());
    path.push(".");
    path.push(label);
    path.push(".");
    path.push(format.extension());
    path.into()
}

fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut output = OsString::from(path.as_os_str());
    output.push(suffix);
    output.into()
}
