use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn run(arguments: &[&str]) -> Output {
    binary().args(arguments).output().unwrap()
}

fn run_stdin(arguments: &[&str], input: &[u8]) -> Output {
    let mut command = binary();
    command
        .args(arguments)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let mut child = command.spawn().unwrap();
    child.stdin.take().unwrap().write_all(input).unwrap();
    child.wait_with_output().unwrap()
}

fn sam_records(bytes: &[u8]) -> Vec<Vec<&str>> {
    std::str::from_utf8(bytes)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('@'))
        .map(|line| line.split('\t').collect())
        .collect()
}

#[test]
fn single_fastq_writes_unmapped_sam_to_stdout() {
    let input = fixture("import-se.fastq");
    let output = run(&["import", "-0", input.to_str().unwrap(), "--no-PG"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::str::from_utf8(&output.stdout).unwrap();
    assert!(
        text.starts_with("@HD\tVN:1.6\tSO:unsorted\tGO:query\n"),
        "{text}"
    );
    let records = sam_records(&output.stdout);
    assert_eq!(records.len(), 3);
    assert_eq!(
        &records[0][..11],
        &[
            "read1",
            "4",
            "*",
            "0",
            "0",
            "*",
            "*",
            "0",
            "0",
            "ACGTACGTACGT",
            "IIIIIIIIIIII"
        ]
    );
}

#[test]
fn single_input_detects_mate_suffixes_per_record() {
    let input = fixture("import-interleaved.fastq");
    let output = run(&["import", input.to_str().unwrap(), "--no-PG"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let records = sam_records(&output.stdout);
    assert_eq!(
        records.iter().map(|record| record[0]).collect::<Vec<_>>(),
        ["pair1", "pair1", "pair2", "pair2"]
    );
    assert_eq!(
        records.iter().map(|record| record[1]).collect::<Vec<_>>(),
        ["77", "141", "77", "141"]
    );
}

#[test]
fn paired_fastq_writes_bam_in_alternating_order() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("reads.bam");
    let r1 = fixture("import-r1.fastq");
    let r2 = fixture("import-r2.fastq");
    let output = run(&[
        "import",
        "-1",
        r1.to_str().unwrap(),
        "-2",
        r2.to_str().unwrap(),
        "--no-PG",
        "-o",
        bam.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let decoded = run(&["view", bam.to_str().unwrap()]);
    assert!(
        decoded.status.success(),
        "{}",
        String::from_utf8_lossy(&decoded.stderr)
    );
    let records = sam_records(&decoded.stdout);
    assert_eq!(
        records.iter().map(|record| record[0]).collect::<Vec<_>>(),
        ["read1", "read1", "read2", "read2"]
    );
    assert_eq!(
        records.iter().map(|record| record[1]).collect::<Vec<_>>(),
        ["77", "141", "77", "141"]
    );
}

#[test]
fn read_group_and_fixed_width_order_tags_are_written() {
    let input = fixture("import-se.fastq");
    let output = run(&[
        "import",
        "-0",
        input.to_str().unwrap(),
        "-r",
        "ID:lib1",
        "-r",
        "SM:sample1",
        "--order",
        "ro:6",
        "--no-PG",
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::str::from_utf8(&output.stdout).unwrap();
    assert!(text.contains("@RG\tID:lib1\tSM:sample1\n"), "{text}");
    let records = sam_records(&output.stdout);
    assert_eq!(&records[0][11..], &["RG:Z:lib1", "ro:Z:000000"]);
    assert_eq!(&records[2][11..], &["RG:Z:lib1", "ro:Z:000002"]);
}

#[test]
fn failed_pair_does_not_replace_named_output() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("reads.bam");
    fs::write(&bam, b"keep me").unwrap();
    let r1 = fixture("import-r1.fastq");
    let longer = fixture("import-se.fastq");
    let output = run(&[
        "import",
        "-1",
        r1.to_str().unwrap(),
        "-2",
        longer.to_str().unwrap(),
        "-o",
        bam.to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert_eq!(fs::read(bam).unwrap(), b"keep me");
}

#[test]
fn named_sam_is_inferred_from_extension() {
    let directory = tempfile::tempdir().unwrap();
    let sam = directory.path().join("reads.sam");
    let input = fixture("import-se.fastq");
    let output = run(&[
        "import",
        "-0",
        input.to_str().unwrap(),
        "--no-PG",
        "-o",
        sam.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(fs::read(sam).unwrap().starts_with(b"@HD\t"));
}

#[test]
fn gzip_and_standard_input_follow_the_same_reader_contract() {
    let directory = tempfile::tempdir().unwrap();
    let input = fs::read(fixture("import-se.fastq")).unwrap();
    let gzip = directory.path().join("reads.fastq.gz");
    let mut encoder = flate2::write::GzEncoder::new(
        fs::File::create(&gzip).unwrap(),
        flate2::Compression::default(),
    );
    encoder.write_all(&input).unwrap();
    encoder.finish().unwrap();

    let compressed = run(&["import", gzip.to_str().unwrap(), "--no-PG"]);
    let standard_input = run_stdin(&["import", "-", "--no-PG"], &input);
    assert!(compressed.status.success());
    assert!(standard_input.status.success());
    assert_eq!(compressed.stdout, standard_input.stdout);
}

#[test]
fn named_output_cannot_alias_an_input() {
    let input = fixture("import-se.fastq");
    let output = run(&[
        "import",
        input.to_str().unwrap(),
        "-o",
        input.to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("also an input path"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn default_header_records_program_provenance() {
    let input = fixture("import-se.fastq");
    let output = run(&["import", input.to_str().unwrap()]);
    assert!(output.status.success());
    let text = std::str::from_utf8(&output.stdout).unwrap();
    assert!(
        text.contains("@PG\tID:rsomics-bam\tPN:rsomics-bam\t"),
        "{text}"
    );
}

#[test]
fn mixed_explicit_and_positional_inputs_are_rejected() {
    let input = fixture("import-se.fastq");
    let output = run(&[
        "import",
        "-0",
        input.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert!(!output.status.success());
}
