use std::io::Write;
use std::process::{Command, Stdio};

fn root() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/upstream/samtools-consensus")
}

#[test]
fn consensus_reads_sam_from_stdin() {
    let mut child = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args([
            "consensus",
            "--mode",
            "simple",
            "--call-fract",
            "0.6",
            "--format",
            "fastq",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&std::fs::read(root().join("consen1.sam")).unwrap())
        .unwrap();

    let output = child.wait_with_output().unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        std::fs::read(root().join("expected/1q.out")).unwrap()
    );
}

#[test]
fn consensus_failure_preserves_existing_output() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("broken.sam");
    let output = directory.path().join("consensus.fa");
    std::fs::write(&input, b"@HD\tVN:1.6\nnot-a-record\n").unwrap();
    std::fs::write(&output, b"preserve me\n").unwrap();

    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["consensus", "--output"])
        .arg(&output)
        .arg(&input)
        .output()
        .unwrap();

    assert!(!result.status.success());
    assert_eq!(std::fs::read(output).unwrap(), b"preserve me\n");
}
