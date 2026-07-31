use std::io::Write;
use std::path::Path;

use noodles::sam::{self, alignment::io::Write as _};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::{input, md};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Options<'a> {
    pub header_lines: Option<usize>,
    pub records: usize,
    pub reference: Option<&'a Path>,
    pub additional_threads: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub header_lines: usize,
    pub records: usize,
}

pub fn write(input_path: &Path, options: Options<'_>, mut output: impl Write) -> Result<Summary> {
    let mut reader = input::open(input_path, options.reference, options.additional_threads)?;
    let header = reader.read_header(input_path)?;
    let header_bytes = format_header(&header)?;
    let header_lines = write_header(&mut output, &header_bytes, options.header_lines)?;
    let format = reader.format();
    let mut reference = if format == input::Format::Cram {
        options
            .reference
            .map(md::ReferenceCache::open)
            .transpose()?
    } else {
        None
    };

    let mut records = 0;
    if options.records > 0 {
        let mut writer = sam::io::Writer::new(&mut output);
        reader.visit_records(&header, input_path, |record| {
            if format == input::Format::Cram {
                let record = md::complete(&header, record, reference.as_mut())?;
                writer
                    .write_alignment_record(&header, &record)
                    .map_err(RsomicsError::Io)?;
            } else {
                writer
                    .write_alignment_record(&header, record)
                    .map_err(RsomicsError::Io)?;
            }
            records += 1;
            Ok(records < options.records)
        })?;
    }

    output.flush().map_err(RsomicsError::Io)?;
    Ok(Summary {
        header_lines,
        records,
    })
}

fn format_header(header: &sam::Header) -> Result<Vec<u8>> {
    let mut writer = sam::io::Writer::new(Vec::new());
    writer.write_header(header).map_err(RsomicsError::Io)?;
    Ok(writer.into_inner())
}

fn write_header(output: &mut impl Write, header: &[u8], limit: Option<usize>) -> Result<usize> {
    let mut written = 0;
    for line in header.split_inclusive(|byte| *byte == b'\n') {
        if limit.is_some_and(|limit| written == limit) {
            break;
        }
        output.write_all(line).map_err(RsomicsError::Io)?;
        if !line.ends_with(b"\n") {
            output.write_all(b"\n").map_err(RsomicsError::Io)?;
        }
        written += 1;
    }
    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_limit_counts_complete_lines() {
        let mut output = Vec::new();
        assert_eq!(
            write_header(&mut output, b"@HD\tVN:1.6\n@SQ\tSN:chr1\n", Some(1)).unwrap(),
            1
        );
        assert_eq!(output, b"@HD\tVN:1.6\n");
    }

    #[test]
    fn unterminated_header_line_gets_a_newline() {
        let mut output = Vec::new();
        assert_eq!(write_header(&mut output, b"@HD\tVN:1.6", None).unwrap(), 1);
        assert_eq!(output, b"@HD\tVN:1.6\n");
    }
}
