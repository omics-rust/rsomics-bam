use std::io::BufReader;
use std::path::{Path, PathBuf};

use noodles::sam::alignment::io::Write as _;
use noodles::{bam, sam};
use rsomics_bam::mpileup::{self, BaqMode, Options};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn default_text_matches_golden() {
    let mut output = Vec::new();
    let summary = mpileup::write(&fixture("records.sam"), Options::default(), &mut output).unwrap();
    assert_eq!(summary.positions, 16);
    assert_eq!(output, include_bytes!("golden/mpileup-default.txt"));
}

#[test]
fn reference_baq_matches_golden() {
    let reference = fixture("reference.fa");
    let mut output = Vec::new();
    let summary = mpileup::write(
        &fixture("records.sam"),
        Options {
            reference: Some(&reference),
            ..Options::default()
        },
        &mut output,
    )
    .unwrap();
    assert_eq!(summary.positions, 16);
    assert_eq!(output, include_bytes!("golden/mpileup-reference.txt"));
}

#[test]
fn raw_bam_and_sam_emit_the_same_pileup() {
    let input = fixture("records.sam");
    let mut reader = sam::io::Reader::new(BufReader::new(std::fs::File::open(&input).unwrap()));
    let header = reader.read_header().unwrap();
    let records = reader
        .record_bufs(&header)
        .collect::<std::io::Result<Vec<_>>>()
        .unwrap();
    let file = tempfile::NamedTempFile::new().unwrap();
    let mut writer = bam::io::Writer::from(file.reopen().unwrap());
    writer.write_header(&header).unwrap();
    for record in records {
        writer.write_alignment_record(&header, &record).unwrap();
    }

    let mut output = Vec::new();
    mpileup::write(file.path(), Options::default(), &mut output).unwrap();
    assert_eq!(output, include_bytes!("golden/mpileup-default.txt"));
}

#[test]
fn indels_and_overlapping_mates_use_samtools_text_encoding() {
    let input = fixture("mpileup-records.sam");
    let mut output = Vec::new();
    mpileup::write(&input, Options::default(), &mut output).unwrap();
    let text = String::from_utf8(output).unwrap();
    assert!(text.contains("chr1\t17\tN\t2\tAT+2TT\tII\n"));
    assert!(text.contains("chr1\t44\tN\t1\tA-3NNN\tI\n"));
    assert!(text.contains("chr1\t45\tN\t1\t*\tI\n"));

    let reference = fixture("mpileup-reference.fa");
    let mut output = Vec::new();
    mpileup::write(
        &input,
        Options {
            reference: Some(&reference),
            baq: BaqMode::Disabled,
            ..Options::default()
        },
        &mut output,
    )
    .unwrap();
    let text = String::from_utf8(output).unwrap();
    assert!(text.contains("chr1\t44\tT\t1\tA-3ACG\tI\n"));
}
