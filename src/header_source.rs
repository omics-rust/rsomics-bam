use std::io::{BufRead, BufReader};
use std::path::Path;

use noodles::sam;
use rsomics_common::Result;

use crate::{bgzf_rewrite, input};

pub(crate) struct Data {
    pub(crate) header: sam::Header,
    pub(crate) text: Vec<u8>,
}

pub(crate) fn read(path: &Path) -> Result<Data> {
    if input::detect_format(path)? == input::Format::Bam {
        let (header, text) = bgzf_rewrite::read_header(path)?;
        return Ok(Data { header, text });
    }

    let mut reader = input::open(path, None, 0)?;
    let header = reader.read_header(path)?;
    let text = match raw_sam_header(path) {
        Some(text) => text,
        None => bgzf_rewrite::canonical_header_text(&header)?,
    };
    Ok(Data { header, text })
}

pub(crate) fn append_line(text: &mut Vec<u8>, line: &[u8]) {
    if !text.is_empty() && !text.ends_with(b"\n") {
        text.push(b'\n');
    }
    text.extend_from_slice(line);
    if !text.ends_with(b"\n") {
        text.push(b'\n');
    }
}

pub(crate) fn read_group_line<'a>(text: &'a [u8], id: &[u8]) -> Option<&'a [u8]> {
    text.split_inclusive(|byte| *byte == b'\n').find(|line| {
        line.starts_with(b"@RG\t")
            && line.split(|byte| *byte == b'\t').any(|field| {
                field
                    .strip_prefix(b"ID:")
                    .is_some_and(|value| value.trim_ascii_end() == id)
            })
    })
}

pub(crate) fn program_line<'a>(text: &'a [u8], id: &[u8]) -> Option<&'a [u8]> {
    text.split_inclusive(|byte| *byte == b'\n').find(|line| {
        line.starts_with(b"@PG\t")
            && line.split(|byte| *byte == b'\t').any(|field| {
                field
                    .strip_prefix(b"ID:")
                    .is_some_and(|value| value.trim_ascii_end() == id)
            })
    })
}

fn raw_sam_header(path: &Path) -> Option<Vec<u8>> {
    let file = std::fs::File::open(path).ok()?;
    let mut reader = BufReader::new(file);
    if reader.fill_buf().ok()?.first() != Some(&b'@') {
        return None;
    }
    let mut text = Vec::new();
    loop {
        let mut line = Vec::new();
        if reader.read_until(b'\n', &mut line).ok()? == 0 || !line.starts_with(b"@") {
            break;
        }
        text.extend_from_slice(&line);
    }
    Some(text)
}
