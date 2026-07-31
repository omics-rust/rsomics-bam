use std::io::{BufWriter, Write};

use noodles::sam::alignment::io::Write as _;
use noodles::{bam, bgzf, sam};
use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum Format {
    #[default]
    Sam,
    Bam,
}

pub(crate) struct Writer<W>
where
    W: Write,
{
    inner: Inner<W>,
}

enum Inner<W>
where
    W: Write,
{
    Sam(sam::io::Writer<BufWriter<W>>),
    Bam(bam::io::Writer<bgzf::io::Writer<W>>),
}

impl<W> Writer<W>
where
    W: Write,
{
    pub(crate) fn new(format: Format, output: W) -> Self {
        let inner = match format {
            Format::Sam => Inner::Sam(sam::io::Writer::new(BufWriter::new(output))),
            Format::Bam => Inner::Bam(bam::io::Writer::new(output)),
        };
        Self { inner }
    }

    pub(crate) fn write_header(&mut self, header: &sam::Header) -> Result<()> {
        match &mut self.inner {
            Inner::Sam(writer) => writer.write_header(header),
            Inner::Bam(writer) => writer.write_header(header),
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
            Inner::Bam(writer) => writer.write_alignment_record(header, record),
        }
        .map_err(RsomicsError::Io)
    }

    pub(crate) fn finish(&mut self, header: &sam::Header) -> Result<()> {
        match &mut self.inner {
            Inner::Sam(writer) => {
                writer.finish(header)?;
                writer.get_mut().flush()
            }
            Inner::Bam(writer) => writer.try_finish(),
        }
        .map_err(RsomicsError::Io)
    }
}
