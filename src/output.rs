use std::io::{BufWriter, Write};
use std::num::NonZero;

use noodles::sam::alignment::io::Write as _;
use noodles::{bam, bgzf, sam};
use rsomics_bamio::RingBgzfWriter;
use rsomics_bamio::raw::{self, RecordRef};
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

    pub(crate) fn finish(self, header: &sam::Header) -> Result<()> {
        match self.inner {
            Inner::Sam(mut writer) => {
                writer.finish(header)?;
                writer.get_mut().flush()
            }
            Inner::BamSingle(mut writer) => writer.try_finish(),
            Inner::BamParallel(writer) => writer.into_inner().finish().map(drop),
            Inner::BamParallelLevel(mut writer) => writer.get_mut().finish().map(drop),
        }
        .map_err(RsomicsError::Io)
    }
}
