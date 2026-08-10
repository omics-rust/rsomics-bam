use std::fs;
use std::io::Write;
use std::path::Path;
use std::process::{Command, Output, Stdio};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn run(mut command: Command) -> Output {
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "stdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn run_with_stdin(mut command: Command, input: &[u8]) -> Output {
    command
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let mut child = command.spawn().unwrap();
    child.stdin.take().unwrap().write_all(input).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "stdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

const PADDED: &[u8] = b"@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:10\n\
chr1\t0\tchr1\t1\t0\t3M1D1M1D4M\t*\t0\t0\tACGTACGT\t~~~~~~~~\n\
r001\t0\tchr1\t1\t60\t3M1D1M1D4M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
r002\t0\tchr1\t1\t60\t3M\t*\t0\t0\tACG\tIII\n\
r003\t0\tchr1\t6\t60\t1D4M\t*\t0\t0\tACGT\tIIII\n\
r004\t0\tchr1\t5\t60\t1M1D4M\t*\t0\t0\tTACGT\tIIIII\n\
r005\t0\tchr1\t5\t60\t1M\t*\t0\t0\tT\tI\n";

const UNPADDED: &[u8] = b"@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:10\n\
chr1\t0\tchr1\t1\t0\t8M\t*\t0\t0\tACGTACGT\t~~~~~~~~\n\
r001\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
r002\t0\tchr1\t1\t60\t3M\t*\t0\t0\tACG\tIII\n\
r003\t0\tchr1\t5\t60\t1P4M\t*\t0\t0\tACGT\tIIII\n\
r004\t0\tchr1\t4\t60\t5M\t*\t0\t0\tTACGT\tIIIII\n\
r005\t0\tchr1\t4\t60\t1M\t*\t0\t0\tT\tI\n";

fn write_fixture(directory: &Path) -> (std::path::PathBuf, std::path::PathBuf) {
    let input = directory.join("padded.sam");
    let reference = directory.join("padded.fa");
    fs::write(&input, PADDED).unwrap();
    fs::write(&reference, b">chr1\nACG*T*ACGT\n").unwrap();
    (input, reference)
}

#[test]
fn embedded_reference_and_fasta_header_projection_match_the_stable_contract() {
    let directory = tempfile::tempdir().unwrap();
    let (input, reference) = write_fixture(directory.path());

    let embedded = run({
        let mut command = binary();
        command.args(["depad", "--no-pg", "-s"]).arg(&input);
        command
    });
    assert_eq!(embedded.stdout, UNPADDED);
    assert_eq!(
        embedded.stderr,
        b"warning: reference lengths remain padded without --reference\n"
    );

    let projected = run({
        let mut command = binary();
        command
            .args(["depad", "--no-pg", "-S", "-O", "sam", "-T"])
            .arg(&reference)
            .arg(&input);
        command
    });
    let expected = String::from_utf8(UNPADDED.to_vec())
        .unwrap()
        .replace("SN:chr1\tLN:10", "SN:chr1\tLN:8");
    assert_eq!(projected.stdout, expected.as_bytes());
    assert!(projected.stderr.is_empty());
    assert!(!directory.path().join("padded.fa.fai").exists());
}

#[test]
fn named_bam_stdin_and_json_use_product_io_contracts() {
    let directory = tempfile::tempdir().unwrap();
    let (_, reference) = write_fixture(directory.path());
    let output = directory.path().join("depad.bam");

    let result = run_with_stdin(
        {
            let mut command = binary();
            command
                .arg("--json")
                .args(["depad", "--no-pg", "-@", "1", "-1", "-T"])
                .arg(&reference)
                .args(["-o"])
                .arg(&output)
                .arg("-");
            command
        },
        PADDED,
    );
    let value: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(value["status"], "ok");
    assert_eq!(value["result"]["command"], "depad");
    assert_eq!(value["result"]["summary"]["records_read"], 6);
    assert_eq!(value["result"]["summary"]["records_projected"], 6);
    assert_eq!(value["result"]["summary"]["embedded_references"], 1);

    let decoded = run({
        let mut command = binary();
        command.args(["view", "-h", "--no-pg"]).arg(&output);
        command
    });
    let expected = String::from_utf8(UNPADDED.to_vec())
        .unwrap()
        .replace("SN:chr1\tLN:10", "SN:chr1\tLN:8");
    assert_eq!(decoded.stdout, expected.as_bytes());
}

#[test]
fn failures_preserve_named_output_and_reject_aliases() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("missing-reference.sam");
    let output = directory.path().join("out.bam");
    fs::write(
        &input,
        b"@SQ\tSN:chr1\tLN:10\nread\t0\tchr1\t1\t60\t3M\t*\t0\t0\tACG\tIII\n",
    )
    .unwrap();
    fs::write(&output, b"sentinel").unwrap();

    let failed = binary()
        .args(["depad", "-o"])
        .arg(&output)
        .arg(&input)
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert!(String::from_utf8_lossy(&failed.stderr).contains("embedded reference"));
    assert_eq!(fs::read(&output).unwrap(), b"sentinel");

    let alias = binary()
        .args(["depad", "-o"])
        .arg(&input)
        .arg(&input)
        .output()
        .unwrap();
    assert!(!alias.status.success());
    assert!(String::from_utf8_lossy(&alias.stderr).contains("different files"));
}

#[test]
fn unmapped_records_are_not_projected() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("unmapped.sam");
    let source = b"@SQ\tSN:chr1\tLN:10\n\
chr1\t0\tchr1\t1\t0\t3M1D1M1D4M\t*\t0\t0\tACGTACGT\t~~~~~~~~\n\
unmapped\t4\tchr1\t6\t0\t3M\tchr1\t5\t0\tACG\tIII\n";
    fs::write(&input, source).unwrap();

    let result = run({
        let mut command = binary();
        command.args(["depad", "--no-pg", "-s"]).arg(&input);
        command
    });
    let text = String::from_utf8(result.stdout).unwrap();
    assert!(text.contains("unmapped\t4\tchr1\t6\t0\t3M\t=\t5\t0\tACG\tIII\n"));
}

#[test]
fn invalid_cigar_and_reference_data_fail_loudly() {
    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("padded.fa");
    fs::write(&reference, b">chr1\nACG*T*ACGT\n").unwrap();

    for (name, cigar, sequence) in [("insertion", "3M1I", "ACGT"), ("pad", "3M1P", "ACG")] {
        let input = directory.path().join(format!("{name}.sam"));
        fs::write(
            &input,
            format!(
                "@SQ\tSN:chr1\tLN:10\n{name}\t0\tchr1\t1\t60\t{cigar}\t*\t0\t0\t{sequence}\t{}\n",
                "I".repeat(sequence.len())
            ),
        )
        .unwrap();
        let result = binary()
            .args(["depad", "-T"])
            .arg(&reference)
            .arg(&input)
            .output()
            .unwrap();
        assert!(!result.status.success(), "{name}");
        assert!(
            String::from_utf8_lossy(&result.stderr).contains("unsupported input CIGAR"),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }

    let (_, reference) = write_fixture(directory.path());
    fs::write(&reference, b">chr1\nACG?T*ACGT\n").unwrap();
    let result = binary()
        .args(["depad", "-T"])
        .arg(&reference)
        .arg(directory.path().join("padded.sam"))
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("invalid base '?'"));
}

#[test]
fn compression_modes_cannot_select_sam_output() {
    let directory = tempfile::tempdir().unwrap();
    let (input, _) = write_fixture(directory.path());
    for option in ["-u", "-1"] {
        let result = binary()
            .args(["depad", "-s", option])
            .arg(&input)
            .output()
            .unwrap();
        assert!(!result.status.success());
        assert!(String::from_utf8_lossy(&result.stderr).contains("compression"));
    }
}
