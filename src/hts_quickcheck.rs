use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use noodles::sam::alignment::RecordBuf;

use crate::input::{self, Format};
use crate::quickcheck::{FileReport, Problem};

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

const CRAM_EOF: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

pub(crate) fn check(path: &Path, allow_no_targets: bool) -> FileReport {
    let mut problems = Vec::new();
    let format = match input::detect_format(path) {
        Ok(format) => format,
        Err(_) => {
            problems.push(Problem::Open);
            return FileReport {
                path: path.to_path_buf(),
                problems,
            };
        }
    };

    match input::open(path, None, 0).and_then(|mut reader| {
        let header = reader.read_header(path)?;
        if !allow_no_targets && header.reference_sequences().is_empty() {
            problems.push(Problem::NoTargets);
        }
        if header.reference_sequences().is_empty()
            && reader
                .visit_records(&header, path, |record| {
                    RecordBuf::try_from_alignment_record(&header, record)
                        .map_err(rsomics_common::RsomicsError::Io)?;
                    record.flags().map_err(rsomics_common::RsomicsError::Io)?;
                    Ok(false)
                })
                .is_err()
        {
            problems.push(Problem::NotSequence);
        }
        Ok(())
    }) {
        Ok(()) => {}
        Err(_) => problems.push(Problem::Header),
    }

    match format {
        Format::Sam => {}
        Format::Bam => check_eof(path, &BGZF_EOF, &mut problems),
        Format::Cram => check_eof(path, &CRAM_EOF, &mut problems),
    }

    FileReport {
        path: path.to_path_buf(),
        problems,
    }
}

fn check_eof<const N: usize>(path: &Path, expected: &[u8; N], problems: &mut Vec<Problem>) {
    let result = File::open(path).and_then(|mut file| {
        let offset = i64::try_from(N).unwrap();
        file.seek(SeekFrom::End(-offset))?;
        let mut actual = [0; N];
        file.read_exact(&mut actual)?;
        Ok(actual == *expected)
    });

    match result {
        Ok(true) => {}
        Ok(false) => problems.push(Problem::MissingEof),
        Err(error) if error.kind() == std::io::ErrorKind::InvalidInput => {
            problems.push(Problem::MissingEof);
        }
        Err(_) => problems.push(Problem::EofCheck),
    }
}
