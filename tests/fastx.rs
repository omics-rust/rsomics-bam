use std::fs;
use std::path::Path;
use std::process::Command;

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn fixture(path: &Path) {
    fs::write(
        path,
        b"@HD\tVN:1.6\tSO:queryname\n\
q1\t69\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n\
q1\t149\t*\t0\t0\t*\t*\t0\t0\tAGTC\tABCD\n\
q2\t69\t*\t0\t0\t*\t*\t0\t0\tAAAA\t*\n\
q2\t69\t*\t0\t0\t*\t*\t0\t0\tCCCC\tJJJJ\n\
q3\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n\
q4\t2052\t*\t0\t0\t*\t*\t0\t0\tTTTT\tIIII\n",
    )
    .unwrap();
}

fn run(args: &[&str]) -> std::process::Output {
    let output = binary().args(args).output().unwrap();
    assert!(
        output.status.success(),
        "stdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

#[test]
fn fasta_and_fastq_match_name_group_selection_and_orientation() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("records.sam");
    fixture(&input);

    let fasta = run(&["fasta", input.to_str().unwrap()]);
    assert_eq!(
        fasta.stdout,
        b">q1/1\nACGT\n>q1/2\nGACT\n>q2/1\nCCCC\n>q3\n\n"
    );

    let fastq = run(&["fastq", input.to_str().unwrap()]);
    assert_eq!(
        fastq.stdout,
        b"@q1/1\nACGT\n+\nIIII\n@q1/2\nGACT\n+\nDCBA\n@q2/1\nCCCC\n+\nJJJJ\n@q3\n\n+\n\n"
    );
}

#[test]
fn named_bgzf_output_is_transactional_and_keeps_json_separate() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("records.sam");
    let output = directory.path().join("records.fq.bgzf");
    fixture(&input);

    let result = run(&[
        "--json",
        "fastq",
        "-o",
        output.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    let envelope: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(envelope["result"]["summary"]["records_written"], 4);

    let mut reader = rsomics_seqio::open_path(&output).unwrap();
    let mut names = Vec::new();
    while let Some(record) = reader.read_record().unwrap() {
        names.push(record.id.to_vec());
    }
    assert_eq!(
        names,
        vec![
            b"q1/1".to_vec(),
            b"q1/2".to_vec(),
            b"q2/1".to_vec(),
            b"q3".to_vec(),
        ]
    );

    fs::write(&output, b"existing\n").unwrap();
    fs::write(&input, b"@HD\tVN:1.6\nbroken\n").unwrap();
    let failed = binary()
        .args(["fastq", "-o"])
        .arg(&output)
        .arg(&input)
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert_eq!(fs::read(output).unwrap(), b"existing\n");
}

#[test]
fn fastq_original_default_quality_and_filter_controls_are_exact() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("quality.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\tSO:queryname\n\
missing\t4\t*\t0\t0\t*\t*\t0\t0\tAC\t*\n\
oq\t20\t*\t0\t0\t*\t*\t0\t0\tAGTC\t!!!!\tOQ:Z:JKLM\n\
secondary\t260\t*\t0\t0\t*\t*\t0\t0\tGG\tHH\n",
    )
    .unwrap();

    let default = run(&["fastq", "-n", "-O", "-v", "7", input.to_str().unwrap()]);
    assert_eq!(default.stdout, b"@missing\nAC\n+\n((\n@oq\nGACT\n+\nMLKJ\n");

    let unfiltered = run(&["fastq", "-n", "-F", "0", input.to_str().unwrap()]);
    assert!(unfiltered.stdout.ends_with(b"@secondary\nGG\n+\nHH\n"));
}

#[test]
fn malformed_original_quality_fails_without_replacing_output() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("bad-oq.sam");
    let output = directory.path().join("output.fq");
    fs::write(
        &input,
        b"@HD\tVN:1.6\nbad\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tOQ:Z:ABC\n",
    )
    .unwrap();
    fs::write(&output, b"existing\n").unwrap();

    let failed = binary()
        .args(["fastq", "-O", "-o"])
        .arg(&output)
        .arg(&input)
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert!(String::from_utf8_lossy(&failed.stderr).contains("OQ length"));
    assert_eq!(fs::read(output).unwrap(), b"existing\n");
}

#[test]
fn json_requires_a_named_sequence_output() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("records.sam");
    fixture(&input);

    let failed = binary()
        .args(["--json", "fasta"])
        .arg(&input)
        .output()
        .unwrap();
    assert_eq!(failed.status.code(), Some(2));
    assert!(failed.stdout.is_empty());
    assert!(String::from_utf8_lossy(&failed.stderr).contains("named --output"));
}
