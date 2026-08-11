use std::fs::OpenOptions;
use std::io;
use std::path::Path;

use noodles::sam;
use rsomics_common::{Result, RsomicsError};
use rust_htslib::bam::{self, Read as _};

use super::{
    FLAG_SECONDARY, FLAG_SUPPLEMENTARY, Options, Summary, TagFilter, complement,
    increment_record_count, reset_header, transformed_flags,
};

pub(super) fn write(
    input_path: &Path,
    options: Options<'_>,
    output_path: &Path,
) -> Result<Summary> {
    let mut reader = if input_path == Path::new("-") {
        bam::Reader::from_stdin()
    } else {
        bam::Reader::from_path(input_path)
    }
    .map_err(|error| input_error("opening", input_path, error))?;
    if let Some(reference) = options.reference {
        reader
            .set_reference(reference)
            .map_err(|error| input_error("setting the reference for", input_path, error))?;
    }
    if options.additional_threads > 0 {
        reader
            .set_threads(options.additional_threads)
            .map_err(|error| input_error("configuring threads for", input_path, error))?;
    }

    let input_header = parse_header(reader.header(), input_path)?;
    let mut output_header = reset_header(&input_header, options.no_rg, options.reject_pg);
    if let Some(program) = options.program {
        program.add_to(&mut output_header)?;
    }
    let hts_header = build_header(&output_header)?;
    let mut writer = bam::Writer::from_path(output_path, &hts_header, bam::Format::Cram)
        .map_err(|error| output_error("opening", output_path, error))?;
    if let Some(reference) = options.reference {
        writer
            .set_reference(reference)
            .map_err(|error| output_error("setting the reference for", output_path, error))?;
    }
    if options.additional_threads > 0 {
        writer
            .set_threads(options.additional_threads)
            .map_err(|error| output_error("configuring threads for", output_path, error))?;
    }

    let tags = TagFilter::new(options.remove_tags, options.keep_tags, options.no_rg);
    let mut record = bam::Record::new();
    let mut written = 0u64;
    while let Some(result) = reader.read(&mut record) {
        result.map_err(|error| input_error("reading", input_path, error))?;
        if reset_record(&mut record, tags, options.keep_duplicate)? {
            writer
                .write(&record)
                .map_err(|error| output_error("writing", output_path, error))?;
            written = increment_record_count(written)?;
        }
    }
    drop(writer);

    OpenOptions::new()
        .write(true)
        .open(output_path)
        .and_then(|file| file.sync_all())
        .map_err(|error| {
            RsomicsError::Io(io::Error::new(
                error.kind(),
                format!("syncing CRAM output {}: {error}", output_path.display()),
            ))
        })?;
    validate(output_path, options.reference, written)?;
    Ok(Summary { written })
}

fn parse_header(header: &bam::HeaderView, input_path: &Path) -> Result<sam::Header> {
    sam::io::Reader::new(header.as_bytes())
        .read_header()
        .map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "reading the alignment header from {}: {error}",
                input_path.display()
            ))
        })
}

fn build_header(header: &sam::Header) -> Result<bam::Header> {
    let mut text = Vec::new();
    sam::io::Writer::new(&mut text)
        .write_header(header)
        .map_err(RsomicsError::Io)?;
    let mut output = bam::Header::new();
    for line in text
        .split(|byte| *byte == b'\n')
        .filter(|line| !line.is_empty())
    {
        if line.len() < 3 || line[0] != b'@' {
            return Err(RsomicsError::InvalidInput(
                "reset produced an invalid alignment header".to_owned(),
            ));
        }
        let mut fields = line.split(|byte| *byte == b'\t');
        let record_type = &fields.next().expect("header line is nonempty")[1..];
        let mut record = bam::header::HeaderRecord::new(record_type);
        for field in fields {
            if field.len() < 4 || field[2] != b':' {
                return Err(RsomicsError::InvalidInput(
                    "reset produced an invalid alignment header field".to_owned(),
                ));
            }
            record.push_tag(&field[..2], String::from_utf8_lossy(&field[3..]));
        }
        output.push_record(&record);
    }
    Ok(output)
}

fn reset_record(
    record: &mut bam::Record,
    tags: TagFilter<'_>,
    keep_duplicate: bool,
) -> Result<bool> {
    let flags = record.flags();
    if flags & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
        return Ok(false);
    }

    let (flags, reverse) = transformed_flags(flags, keep_duplicate);
    let name = record.qname().to_vec();
    let mut sequence = record.seq().as_bytes();
    let mut quality = record.qual().to_vec();
    if reverse {
        sequence.reverse();
        sequence
            .iter_mut()
            .for_each(|base| *base = complement(*base));
        quality.reverse();
    }
    record.set(&name, None, &sequence, &quality);
    record.set_flags(flags);
    record.set_tid(-1);
    record.set_pos(-1);
    record.set_mapq(0);
    record.set_mtid(-1);
    record.set_mpos(-1);
    record.set_insert_size(0);

    let remove = record
        .aux_iter()
        .map(|field| {
            field
                .map(|(tag, _)| <[u8; 2]>::try_from(tag).expect("HTSlib returned a two-byte tag"))
                .map_err(|error| {
                    input_error("reading auxiliary data from", Path::new("record"), error)
                })
        })
        .collect::<Result<Vec<_>>>()?;
    for tag in remove.into_iter().filter(|tag| tags.remove(*tag)) {
        record.remove_aux(&tag).map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "removing auxiliary tag {}{}: {error}",
                char::from(tag[0]),
                char::from(tag[1])
            ))
        })?;
    }
    Ok(true)
}

fn validate(path: &Path, reference: Option<&Path>, expected_records: u64) -> Result<()> {
    let mut reader =
        bam::Reader::from_path(path).map_err(|error| output_error("validating", path, error))?;
    if let Some(reference) = reference {
        reader
            .set_reference(reference)
            .map_err(|error| output_error("validating the reference for", path, error))?;
    }
    let mut record = bam::Record::new();
    let mut records = 0u64;
    while let Some(result) = reader.read(&mut record) {
        result.map_err(|error| output_error("decoding", path, error))?;
        records = increment_record_count(records)?;
    }
    if records != expected_records {
        return Err(RsomicsError::InvalidInput(format!(
            "CRAM output {} contains {records} records; expected {expected_records}",
            path.display()
        )));
    }
    Ok(())
}

fn input_error(action: &str, path: &Path, error: rust_htslib::errors::Error) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{action} {}: {error}", path.display()))
}

fn output_error(action: &str, path: &Path, error: rust_htslib::errors::Error) -> RsomicsError {
    RsomicsError::Io(io::Error::other(format!(
        "{action} {}: {error}",
        path.display()
    )))
}
