use std::collections::{BTreeMap, HashSet};
use std::fs::{self, File};
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use noodles::bam;
use noodles::sam::header::record::value::map::{header::tag, tag::Other};
use rsomics_bamio::raw::RecordReader;

type HeaderTag = Other<tag::Standard>;
type RecordSnapshot = (Vec<u8>, u16, Vec<u8>);

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/sort-records.sam")
}

fn run(arguments: &[&str]) -> Output {
    let output = bin().args(arguments).output().unwrap();
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn fail(arguments: &[&str]) -> Output {
    let output = bin().args(arguments).output().unwrap();
    assert!(!output.status.success());
    output
}

fn snapshot(path: &Path) -> (Vec<RecordSnapshot>, noodles::sam::Header) {
    let mut reader = bam::io::reader::Builder.build_from_path(path).unwrap();
    let header = reader.read_header().unwrap();
    let mut records = RecordReader::new(reader.get_mut());
    let mut values = Vec::new();
    while let Some(record) = records.next().unwrap() {
        values.push((
            record.name().to_vec(),
            record.flags(),
            record.as_bytes().to_vec(),
        ));
    }
    (values, header)
}

fn header_value(header: &noodles::sam::Header, key: HeaderTag) -> Option<&[u8]> {
    header
        .header()?
        .other_fields()
        .get(&key)
        .map(|value| value.as_slice())
}

fn assert_grouped(records: &[RecordSnapshot]) {
    let mut completed = HashSet::new();
    let mut previous = None::<&[u8]>;
    for (name, _, _) in records {
        if previous != Some(name) {
            if let Some(previous) = previous {
                completed.insert(previous.to_vec());
            }
            assert!(!completed.contains(name), "non-contiguous QNAME group");
            previous = Some(name);
        }
    }
}

#[test]
fn groups_names_and_sets_collate_header() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("collated.bam");
    run(&[
        "collate",
        "--no-PG",
        "-@",
        "0",
        fixture().to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    let (records, header) = snapshot(&output);
    assert_eq!(records.len(), 9);
    assert_grouped(&records);
    assert_eq!(
        header_value(&header, tag::SORT_ORDER),
        Some(b"unsorted".as_slice())
    );
    assert_eq!(
        header_value(&header, tag::GROUP_ORDER),
        Some(b"query".as_slice())
    );
    assert_eq!(
        header_value(&header, tag::SUBSORT_ORDER),
        Some(b"unsorted:old".as_slice())
    );

    for group in records.chunk_by(|left, right| left.0 == right.0) {
        if group.len() == 2 {
            assert_eq!((group[0].1 >> 6) & 3, 1);
            assert_eq!((group[1].1 >> 6) & 3, 2);
        }
    }
}

#[test]
fn external_runs_are_bounded_merged_and_removed() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("large.sam");
    write_large_sam(&input, 1_800, 16_000);
    let output = directory.path().join("collated.bam");
    let result = run(&[
        "--json",
        "collate",
        "-m",
        "1M",
        "-@",
        "0",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    let value: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(value["result"]["summary"]["records"], 1_800);
    assert!(
        value["result"]["summary"]["temporary_runs"]
            .as_u64()
            .unwrap()
            > 32
    );
    assert!(value["result"]["summary"]["merge_passes"].as_u64().unwrap() >= 2);

    let (records, _) = snapshot(&output);
    assert_eq!(records.len(), 1_800);
    assert_grouped(&records);
    let prefix = output.file_name().unwrap().to_string_lossy();
    assert!(!fs::read_dir(directory.path()).unwrap().any(|entry| {
        entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .starts_with(&format!("{prefix}.tmp."))
    }));
}

#[test]
fn failures_preserve_outputs_and_aliases_are_rejected() {
    let directory = tempfile::tempdir().unwrap();
    let invalid = directory.path().join("invalid.sam");
    fs::write(
        &invalid,
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nread\t0\tchr1\t1\t60\t10M\t*\t0\t0\tA\tF\n",
    )
    .unwrap();
    let output = directory.path().join("output.bam");
    fs::write(&output, b"sentinel").unwrap();
    fail(&[
        "collate",
        invalid.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    assert_eq!(fs::read(&output).unwrap(), b"sentinel");

    let alias = directory.path().join("alias.sam");
    fs::hard_link(&invalid, &alias).unwrap();
    let failed = fail(&[
        "collate",
        invalid.to_str().unwrap(),
        "-o",
        alias.to_str().unwrap(),
    ]);
    assert!(String::from_utf8_lossy(&failed.stderr).contains("different files"));
}

#[test]
fn stdout_is_a_complete_bam() {
    let output = run(&["collate", "--no-PG", "-@", "0", fixture().to_str().unwrap()]);
    let mut reader = bam::io::Reader::new(output.stdout.as_slice());
    let header = reader.read_header().unwrap();
    assert_eq!(
        header_value(&header, tag::GROUP_ORDER),
        Some(b"query".as_slice())
    );
    let mut records = RecordReader::new(reader.get_mut());
    let mut count = 0;
    while records.next().unwrap().is_some() {
        count += 1;
    }
    assert_eq!(count, 9);
}

#[test]
#[ignore = "requires samtools 1.24"]
fn sam_bam_cram_and_external_runs_match_samtools_1_24() {
    assert_eq!(samtools_version(), "samtools 1.24");
    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    let mut fasta = File::create(&reference).unwrap();
    writeln!(
        fasta,
        ">chr1\n{}\n>chr2\n{}",
        "A".repeat(1000),
        "C".repeat(1000)
    )
    .unwrap();
    samtools(&["faidx", reference.to_str().unwrap()]);

    let bam = directory.path().join("input.bam");
    samtools(&[
        "view",
        "-b",
        "-o",
        bam.to_str().unwrap(),
        fixture().to_str().unwrap(),
    ]);
    let cram = directory.path().join("input.cram");
    samtools(&[
        "view",
        "-C",
        "-T",
        reference.to_str().unwrap(),
        "-o",
        cram.to_str().unwrap(),
        fixture().to_str().unwrap(),
    ]);

    for input in [fixture(), bam, cram] {
        let extension = input.extension().unwrap().to_string_lossy();
        let ours = directory.path().join(format!("ours-{extension}.bam"));
        let oracle = directory.path().join(format!("oracle-{extension}.bam"));
        let mut arguments = vec!["collate", "--no-PG", "-@", "0"];
        if extension == "cram" {
            arguments.extend(["--reference", reference.to_str().unwrap()]);
        }
        arguments.push(input.to_str().unwrap());
        arguments.extend(["-o", ours.to_str().unwrap()]);
        run(&arguments);
        samtools(&[
            "collate",
            "--no-PG",
            "-@",
            "0",
            "-o",
            oracle.to_str().unwrap(),
            input.to_str().unwrap(),
        ]);

        let (ours_records, ours_header) = snapshot(&ours);
        let (oracle_records, oracle_header) = snapshot(&oracle);
        assert_grouped(&ours_records);
        assert_grouped(&oracle_records);
        assert_eq!(record_multiset(&ours), record_multiset(&oracle));
        assert_eq!(
            header_value(&ours_header, tag::SORT_ORDER),
            Some(b"unsorted".as_slice())
        );
        assert_eq!(
            header_value(&ours_header, tag::GROUP_ORDER),
            Some(b"query".as_slice())
        );
        assert_eq!(
            header_value(&oracle_header, tag::SORT_ORDER),
            Some(b"unsorted".as_slice())
        );
        assert_eq!(
            header_value(&oracle_header, tag::GROUP_ORDER),
            Some(b"query".as_slice())
        );
    }

    let large = directory.path().join("oracle-large.sam");
    write_large_sam(&large, 4_000, 200);
    let ours = directory.path().join("ours-large.bam");
    let oracle = directory.path().join("samtools-large.bam");
    run(&[
        "collate",
        "--no-PG",
        "-m",
        "1M",
        "-@",
        "0",
        large.to_str().unwrap(),
        "-o",
        ours.to_str().unwrap(),
    ]);
    samtools(&[
        "collate",
        "--no-PG",
        "-@",
        "0",
        "-o",
        oracle.to_str().unwrap(),
        large.to_str().unwrap(),
    ]);
    let ours_records = snapshot(&ours).0;
    assert_grouped(&ours_records);
    assert_eq!(record_multiset(&ours), record_multiset(&oracle));
}

fn record_multiset(path: &Path) -> BTreeMap<Vec<u8>, usize> {
    let mut counts = BTreeMap::new();
    let output = samtools(&["view", path.to_str().unwrap()]);
    for line in output.stdout.split(|byte| *byte == b'\n') {
        if line.is_empty() {
            continue;
        }
        let mut fields = line
            .split(|byte| *byte == b'\t')
            .map(<[u8]>::to_vec)
            .collect::<Vec<_>>();
        fields[11..].sort_unstable();
        *counts.entry(fields.join(&b'\t')).or_default() += 1;
    }
    counts
}

fn write_large_sam(path: &Path, records: usize, read_length: usize) {
    let mut file = File::create(path).unwrap();
    let sequence = "A".repeat(read_length);
    let quality = "F".repeat(read_length);
    writeln!(file, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000").unwrap();
    for i in (0..records).rev() {
        writeln!(
            file,
            "read{:08}\t{}\tchr1\t{}\t60\t{read_length}M\t*\t0\t0\t{}\t{}",
            i / 2,
            if i % 2 == 0 { 65 } else { 129 },
            i + 1,
            sequence,
            quality
        )
        .unwrap();
    }
}

fn samtools(arguments: &[&str]) -> Output {
    let output = Command::new("samtools").args(arguments).output().unwrap();
    assert!(
        output.status.success(),
        "samtools stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn samtools_version() -> String {
    let output = Command::new("samtools").arg("--version").output().unwrap();
    String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .to_owned()
}
