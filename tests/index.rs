use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use noodles::csi::BinningIndex as _;
use noodles::sam::alignment::io::Write as _;
use noodles::{bgzf, cram, fasta, sam};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn appended_extension(path: &Path, extension: &str) -> PathBuf {
    let mut path = path.as_os_str().to_owned();
    path.push(".");
    path.push(extension);
    PathBuf::from(path)
}

fn run_ours(arguments: &[&str]) -> Output {
    let output = binary().args(arguments).output().expect("spawn command");
    assert!(
        output.status.success(),
        "status={:?}\nstdout={}\nstderr={}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn fail_ours(arguments: &[&str]) -> Output {
    let output = binary().args(arguments).output().expect("spawn command");
    assert!(
        !output.status.success(),
        "stdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn make_bam(source: &Path, target: &Path) {
    run_ours(&[
        "view",
        "--no-pg",
        "-b",
        "-o",
        target.to_str().unwrap(),
        source.to_str().unwrap(),
    ]);
}

fn make_bgzf_sam(target: &Path) {
    let mut writer = bgzf::io::Writer::new(File::create(target).unwrap());
    writer
        .write_all(&fs::read(golden("records.sam")).unwrap())
        .unwrap();
    writer.try_finish().unwrap();
}

fn make_cram(directory: &Path) -> (PathBuf, PathBuf) {
    let reference = directory.join("reference.fa");
    fs::copy(golden("reference.fa"), &reference).unwrap();
    fs::copy(
        golden("reference.fa.fai"),
        appended_extension(&reference, "fai"),
    )
    .unwrap();

    let bytes = fs::read(golden("records.sam")).unwrap();
    let mut source = sam::io::Reader::new(bytes.as_slice());
    let header = source.read_header().unwrap();
    let records = source
        .record_bufs(&header)
        .collect::<std::io::Result<Vec<_>>>()
        .unwrap();
    let repository = fasta::io::indexed_reader::Builder::default()
        .build_from_path(&reference)
        .map(fasta::repository::adapters::IndexedReader::new)
        .map(fasta::Repository::new)
        .unwrap();
    let input = directory.join("records.cram");
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(&input)
        .unwrap();
    writer.write_header(&header).unwrap();
    for record in records {
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish(&header).unwrap();
    (input, reference)
}

#[test]
fn builds_queryable_bai_and_csi_indexes() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("records.bam");
    make_bam(&golden("records.sam"), &input);

    run_ours(&["index", input.to_str().unwrap()]);
    let bai = appended_extension(&input, "bai");
    assert_eq!(
        noodles::bam::bai::fs::read(&bai)
            .unwrap()
            .reference_sequences()
            .len(),
        1
    );
    assert_eq!(
        run_ours(&["view", "-c", input.to_str().unwrap(), "chr1:1-8"]).stdout,
        b"1\n"
    );

    fs::remove_file(&bai).unwrap();
    run_ours(&["index", "-c", input.to_str().unwrap()]);
    let csi_path = appended_extension(&input, "csi");
    let csi = noodles::csi::fs::read(&csi_path).unwrap();
    assert_eq!((csi.min_shift(), csi.depth()), (14, 0));
    assert_eq!(
        run_ours(&["view", "-c", input.to_str().unwrap(), "chr1:1-8"]).stdout,
        b"1\n"
    );

    let custom = directory.path().join("custom.csi");
    let output = run_ours(&[
        "--json",
        "index",
        "-m",
        "1",
        "-o",
        custom.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["result"]["command"], "index");
    assert_eq!(value["result"]["summaries"][0]["kind"], "csi");
    assert_eq!(value["result"]["summaries"][0]["min_shift"], 1);
    assert!(value["result"]["summaries"][0]["additional_threads"].is_number());
    let csi = noodles::csi::fs::read(custom).unwrap();
    assert_eq!((csi.min_shift(), csi.depth()), (1, 3));

    let single_threaded = directory.path().join("single-threaded.bai");
    let output = run_ours(&[
        "--json",
        "index",
        "-@",
        "0",
        "-o",
        single_threaded.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["result"]["summaries"][0]["additional_threads"], 0);
}

#[test]
fn indexes_bgzf_sam_with_bai_or_csi() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("records.sam.gz");
    make_bgzf_sam(&input);

    run_ours(&["index", input.to_str().unwrap()]);
    let bai = appended_extension(&input, "bai");
    assert!(bai.is_file());
    assert_eq!(
        run_ours(&["view", "-c", input.to_str().unwrap(), "chr1:1-8"]).stdout,
        b"1\n"
    );

    fs::remove_file(bai).unwrap();
    run_ours(&["index", "-c", input.to_str().unwrap()]);
    assert_eq!(
        run_ours(&["view", "-c", input.to_str().unwrap(), "chr1:1-8"]).stdout,
        b"1\n"
    );
}

#[test]
fn indexes_cram_as_crai() {
    let directory = tempfile::tempdir().unwrap();
    let (input, reference) = make_cram(directory.path());

    run_ours(&["index", "-c", input.to_str().unwrap()]);
    let index = appended_extension(&input, "crai");
    assert!(!cram::crai::fs::read(&index).unwrap().is_empty());
    assert_eq!(
        run_ours(&[
            "view",
            "-c",
            "-T",
            reference.to_str().unwrap(),
            input.to_str().unwrap(),
            "chr1:1-8",
        ])
        .stdout,
        b"1\n"
    );

    let threaded = directory.path().join("threaded.crai");
    run_ours(&[
        "index",
        "-@",
        "1",
        "-o",
        threaded.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert!(!cram::crai::fs::read(threaded).unwrap().is_empty());
}

#[test]
fn indexing_failures_preserve_existing_outputs() {
    let directory = tempfile::tempdir().unwrap();
    let valid = directory.path().join("valid.bam");
    make_bam(&golden("records.sam"), &valid);
    let truncated = directory.path().join("truncated.bam");
    let mut bytes = fs::read(&valid).unwrap();
    bytes.truncate(bytes.len() / 2);
    fs::write(&truncated, bytes).unwrap();
    let output = directory.path().join("existing.bai");
    fs::write(&output, b"keep\n").unwrap();

    fail_ours(&[
        "index",
        "-o",
        output.to_str().unwrap(),
        truncated.to_str().unwrap(),
    ]);
    assert_eq!(fs::read(&output).unwrap(), b"keep\n");

    let unsorted = directory.path().join("unsorted.bam");
    make_bam(&golden("depth-unsorted.sam"), &unsorted);
    let unsorted_index = appended_extension(&unsorted, "bai");
    let failed = fail_ours(&["index", unsorted.to_str().unwrap()]);
    assert!(String::from_utf8_lossy(&failed.stderr).contains("indexing alignment"));
    assert!(!unsorted_index.exists());

    let alias = directory.path().join("alias.bai");
    fs::hard_link(&valid, &alias).unwrap();
    let failed = fail_ours(&[
        "index",
        "-o",
        alias.to_str().unwrap(),
        valid.to_str().unwrap(),
    ]);
    assert!(
        String::from_utf8_lossy(&failed.stderr).contains("cannot overwrite an alignment input")
    );
    assert_eq!(fs::read(&alias).unwrap(), fs::read(&valid).unwrap());
}

#[test]
fn supports_legacy_and_multiple_input_forms() {
    let directory = tempfile::tempdir().unwrap();
    let first = directory.path().join("first.bam");
    let second = directory.path().join("second.bam");
    make_bam(&golden("records.sam"), &first);
    fs::copy(&first, &second).unwrap();

    run_ours(&[
        "index",
        "-M",
        first.to_str().unwrap(),
        second.to_str().unwrap(),
    ]);
    assert!(appended_extension(&first, "bai").is_file());
    assert!(appended_extension(&second, "bai").is_file());

    let legacy = directory.path().join("legacy.bai");
    run_ours(&["index", first.to_str().unwrap(), legacy.to_str().unwrap()]);
    assert!(legacy.is_file());

    let failed = fail_ours(&["index", first.to_str().unwrap(), second.to_str().unwrap()]);
    assert!(String::from_utf8_lossy(&failed.stderr).contains("require --multiple"));
}
