use std::io::Write;
use std::path::Path;

use noodles::sam;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::Data;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::header::record::value::{Map, map::ReadGroup};
use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::{Program, input, md, output};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Mode {
    #[default]
    OverwriteAll,
    OrphanOnly,
}

#[derive(Clone, Copy, Debug)]
pub enum Source<'a> {
    New(&'a str),
    Existing(&'a str),
    First,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Format {
    Sam,
    Bam,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Compression {
    #[default]
    Default,
    Uncompressed,
}

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub source: Source<'a>,
    pub mode: Mode,
    pub overwrite_header: bool,
    pub format: Format,
    pub compression: Compression,
    pub reference: Option<&'a Path>,
    pub additional_threads: usize,
    pub program: Option<Program<'a>>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub records_read: u64,
    pub records_modified: u64,
    pub records_preserved: u64,
}

pub fn write<W>(input_path: &Path, options: Options<'_>, output: W) -> Result<Summary>
where
    W: Write + Send + 'static,
{
    if options.additional_threads > 256 {
        return Err(RsomicsError::ConfigError(
            "addreplacerg additional thread count cannot exceed 256".to_owned(),
        ));
    }

    let input_threads = match (input_path == Path::new("-"), options.format) {
        (true, _) | (_, Format::Bam) => 0,
        (false, Format::Sam) if input::detect_format(input_path)? == input::Format::Bam => {
            options.additional_threads
        }
        _ => 0,
    };
    let mut reader = input::open(input_path, options.reference, input_threads)?;
    let input_format = reader.format();
    let input_header = reader.read_header(input_path)?;
    let mut output_header = input_header.clone();
    let read_group_id = configure_header(
        &mut output_header,
        options.source,
        options.mode,
        options.overwrite_header,
    )?;
    if let Some(program) = options.program {
        program.add_to(&mut output_header)?;
    }

    let output_format = match options.format {
        Format::Sam => output::Format::Sam,
        Format::Bam => output::Format::Bam,
    };
    let compression = match options.compression {
        Compression::Default => output::Compression::Default,
        Compression::Uncompressed => output::Compression::Uncompressed,
    };
    let mut writer = output::Writer::new(
        output_format,
        compression,
        options.additional_threads,
        output,
    );
    writer.write_header(&output_header)?;

    let mut summary = Summary::default();
    if input_format == input::Format::Bam && options.format == Format::Bam {
        let mut encoded_id = read_group_id.clone();
        encoded_id.push(0);
        reader.visit_owned_raw_records(&input_header, input_path, |mut record| {
            let modified = stamp_raw(&mut record, &encoded_id, options.mode)?;
            update_summary(&mut summary, modified)?;
            writer.write_owned_raw_record(&record)?;
            Ok(true)
        })?;
    } else {
        let mut reference = if input_format == input::Format::Cram {
            options
                .reference
                .map(md::ReferenceCache::open)
                .transpose()?
        } else {
            None
        };
        reader.visit_records(&input_header, input_path, |record| {
            let mut record = if input_format == input::Format::Cram {
                md::complete(&input_header, record, reference.as_mut())?
            } else {
                sam::alignment::RecordBuf::try_from_alignment_record(&input_header, record)
                    .map_err(RsomicsError::Io)?
            };
            let modified = stamp_record(&mut record, &read_group_id, options.mode);
            update_summary(&mut summary, modified)?;
            writer.write_record(&output_header, &record)?;
            Ok(true)
        })?;
    }
    writer.finish(&output_header)?;
    Ok(summary)
}

fn configure_header(
    header: &mut sam::Header,
    source: Source<'_>,
    mode: Mode,
    overwrite_header: bool,
) -> Result<Vec<u8>> {
    match source {
        Source::New(fields) => {
            let (id, read_group) = parse_read_group(fields)?;
            if header.read_groups().contains_key(id.as_slice()) && !overwrite_header {
                return Err(RsomicsError::InvalidInput(format!(
                    "read group {} already exists; use --overwrite-header to replace it",
                    String::from_utf8_lossy(&id)
                )));
            }
            header.read_groups_mut().shift_remove(id.as_slice());
            header
                .read_groups_mut()
                .insert(id.clone().into(), read_group);
            if mode == Mode::OverwriteAll {
                header
                    .read_groups_mut()
                    .retain(|candidate, _| candidate.as_slice() == id);
            }
            Ok(id)
        }
        Source::Existing(id) => {
            validate_id(id)?;
            header
                .read_groups()
                .get_key_value(id.as_bytes())
                .map(|(id, _)| id.to_vec())
                .ok_or_else(|| {
                    RsomicsError::InvalidInput(format!(
                        "read group {id} does not exist in the input header"
                    ))
                })
        }
        Source::First => header
            .read_groups()
            .first()
            .map(|(id, _)| id.to_vec())
            .ok_or_else(|| {
                RsomicsError::InvalidInput(
                    "no read group was supplied and the input header has none".to_owned(),
                )
            }),
    }
}

fn parse_read_group(fields: &str) -> Result<(Vec<u8>, Map<ReadGroup>)> {
    let fields = unescape(fields)?;
    let fields = fields.strip_prefix("@RG\t").unwrap_or(&fields);
    if fields.is_empty() || fields.starts_with('@') {
        return Err(RsomicsError::InvalidInput(
            "new read group must contain @RG fields".to_owned(),
        ));
    }
    let text = format!("@RG\t{fields}\n");
    let mut reader = sam::io::Reader::new(text.as_bytes());
    let mut parsed = reader.read_header().map_err(|error| {
        RsomicsError::InvalidInput(format!("invalid read-group fields: {error}"))
    })?;
    if parsed.read_groups().len() != 1 {
        return Err(RsomicsError::InvalidInput(
            "new read group must contain exactly one ID field".to_owned(),
        ));
    }
    let (id, read_group) = parsed
        .read_groups_mut()
        .shift_remove_index(0)
        .expect("one read group was parsed");
    validate_id_bytes(id.as_slice())?;
    Ok((id.to_vec(), read_group))
}

fn unescape(value: &str) -> Result<String> {
    let mut output = String::with_capacity(value.len());
    let mut characters = value.chars();
    while let Some(character) = characters.next() {
        if character != '\\' {
            if matches!(character, '\0' | '\r' | '\n') {
                return Err(RsomicsError::InvalidInput(
                    "read-group fields cannot contain NUL or line breaks".to_owned(),
                ));
            }
            output.push(character);
            continue;
        }
        match characters.next() {
            Some('\\') => output.push('\\'),
            Some('t') => output.push('\t'),
            Some('n') => {
                return Err(RsomicsError::InvalidInput(
                    "\\n is not supported in read-group fields".to_owned(),
                ));
            }
            Some(_) => {
                return Err(RsomicsError::InvalidInput(
                    "unsupported escape in read-group fields".to_owned(),
                ));
            }
            None => {
                return Err(RsomicsError::InvalidInput(
                    "unterminated escape in read-group fields".to_owned(),
                ));
            }
        }
    }
    Ok(output)
}

fn validate_id(id: &str) -> Result<()> {
    validate_id_bytes(id.as_bytes())
}

fn validate_id_bytes(id: &[u8]) -> Result<()> {
    if id.is_empty() || id.iter().any(|byte| !matches!(byte, 0x21..=0x7e)) {
        return Err(RsomicsError::InvalidInput(
            "read-group ID must be nonempty printable ASCII without whitespace".to_owned(),
        ));
    }
    Ok(())
}

fn stamp_raw(record: &mut RawRecord, encoded_id: &[u8], mode: Mode) -> Result<bool> {
    if mode == Mode::OrphanOnly && record.aux_value(*b"RG").is_some() {
        return Ok(false);
    }
    record.set_aux(*b"RG", b'Z', encoded_id)?;
    Ok(true)
}

fn stamp_record(record: &mut sam::alignment::RecordBuf, id: &[u8], mode: Mode) -> bool {
    if mode == Mode::OrphanOnly && record.data().get(&Tag::READ_GROUP).is_some() {
        return false;
    }
    let mut data = record
        .data()
        .iter()
        .filter(|(tag, _)| *tag != Tag::READ_GROUP)
        .map(|(tag, value)| (tag, value.clone()))
        .collect::<Data>();
    data.insert(Tag::READ_GROUP, Value::String(id.to_vec().into()));
    *record.data_mut() = data;
    true
}

fn update_summary(summary: &mut Summary, modified: bool) -> Result<()> {
    summary.records_read = checked_increment(summary.records_read)?;
    if modified {
        summary.records_modified = checked_increment(summary.records_modified)?;
    } else {
        summary.records_preserved = checked_increment(summary.records_preserved)?;
    }
    Ok(())
}

fn checked_increment(value: u64) -> Result<u64> {
    value
        .checked_add(1)
        .ok_or_else(|| RsomicsError::InvalidInput("alignment count exceeds u64".to_owned()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_repeated_fields_and_supported_escapes() {
        let (id, read_group) = parse_read_group("ID:rg1\\tSM:sample\\\\one").unwrap();
        assert_eq!(id, b"rg1");
        assert_eq!(
            read_group.other_fields().get(b"SM"),
            Some(&b"sample\\one"[..].into())
        );
    }

    #[test]
    fn rejects_unsupported_or_missing_fields() {
        assert!(parse_read_group("SM:sample").is_err());
        assert!(parse_read_group("ID:rg1\\nSM:sample").is_err());
        assert!(parse_read_group("ID:rg1\\qSM:sample").is_err());
        assert!(parse_read_group("ID:rg1\\").is_err());
    }
}
