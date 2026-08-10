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

const DEFAULT: &[u8] = b"@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:40\n\
@RG\tID:rg1\tSM:sample-a\tLB:lib-a\n\
read1\t99\tchr1\t1\t60\t8M\t=\t17\t24\tACGTACGT\tIIIIIIII\tRG:Z:rg1\tNM:i:0\tMD:Z:8\n\
read1\t147\tchr1\t17\t60\t8M\t=\t1\t-24\tACGTACGT\tIIIIIIII\tRG:Z:rg1\tNM:i:0\tMD:Z:8\n\
unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tN\t!\n";

const COMPLEX_REFERENCE: &[u8] =
    b">chr1\nACGTACGTACGTACGTACGTAACCGGTTACGTACGTACGTACGTNNNNACGTACGTACGTACGT\n";

const COMPLEX_INPUT: &[u8] = b"@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:64\n\
perfect\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
wrong\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACGAAC\tIIIIIIIIII\tAA:Z:first\tMD:Z:10\tBB:i:2\tNM:i:0\tCC:Z:last\n\
deletion\t0\tchr1\t1\t60\t5M2D5M\t*\t0\t0\tACGTATACGT\tIIIIIIIIII\n\
insertion\t0\tchr1\t1\t60\t5M3I5M\t*\t0\t0\tACGTAGGGCGTAC\tIIIIIIIIIIIII\n\
clipped\t0\tchr1\t1\t60\t3S10M\t*\t0\t0\tTTTACGTACGTAC\tIIIIIIIIIIIII\n\
skipped\t0\tchr1\t1\t60\t5M3N5M\t*\t0\t0\tACGTAACGTA\tIIIIIIIIII\n\
correct\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tMD:Z:10\tNM:i:0\n\
unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n\
missing\t256\tchr1\t1\t60\t10M\t*\t0\t0\t*\t*\n";

const COMPLEX_EXPECTED: &[u8] = b"@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:64\n\
perfect\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\tMD:Z:10\n\
wrong\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACGAAC\tIIIIIIIIII\tAA:Z:first\tBB:i:2\tCC:Z:last\tNM:i:2\tMD:Z:2G4T2\n\
deletion\t0\tchr1\t1\t60\t5M2D5M\t*\t0\t0\tACGTATACGT\tIIIIIIIIII\tNM:i:2\tMD:Z:5^CG5\n\
insertion\t0\tchr1\t1\t60\t5M3I5M\t*\t0\t0\tACGTAGGGCGTAC\tIIIIIIIIIIIII\tNM:i:3\tMD:Z:10\n\
clipped\t0\tchr1\t1\t60\t3S10M\t*\t0\t0\tTTTACGTACGTAC\tIIIIIIIIIIIII\tNM:i:0\tMD:Z:10\n\
skipped\t0\tchr1\t1\t60\t5M3N5M\t*\t0\t0\tACGTAACGTA\tIIIIIIIIII\tNM:i:0\tMD:Z:10\n\
correct\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tMD:Z:10\tNM:i:0\n\
unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n\
missing\t256\tchr1\t1\t60\t10M\t*\t0\t0\t*\t*\n";

fn write_complex_fixture(directory: &Path) -> (PathBuf, PathBuf) {
    let input = directory.join("complex.sam");
    let reference = directory.join("complex.fa");
    fs::write(&input, COMPLEX_INPUT).unwrap();
    fs::write(&reference, COMPLEX_REFERENCE).unwrap();
    fs::write(directory.join("complex.fa.fai"), b"chr1\t64\t6\t64\t65\n").unwrap();
    (input, reference)
}

#[test]
fn default_and_equal_modes_match_the_stable_sam_contract() {
    let input = fixture("records.sam");
    let reference = fixture("reference.fa");
    let default = run({
        let mut command = binary();
        command
            .args(["calmd", "--no-pg"])
            .arg(&input)
            .arg(&reference);
        command
    });
    assert_eq!(default.stdout, DEFAULT);

    let equal = run({
        let mut command = binary();
        command
            .args(["calmd", "--no-pg", "-e"])
            .arg(&input)
            .arg(&reference);
        command
    });
    let expected = String::from_utf8(DEFAULT.to_vec())
        .unwrap()
        .replace("ACGTACGT\tIIIIIIII", "========\tIIIIIIII");
    assert_eq!(equal.stdout, expected.as_bytes());
}

#[test]
fn named_bam_stdin_and_json_use_product_io_contracts() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("calmd.bam");
    let input = fixture("records.sam");
    let reference = fixture("reference.fa");
    let json = run({
        let mut command = binary();
        command
            .arg("--json")
            .args(["calmd", "--no-pg", "-@", "1", "-o"])
            .arg(&output)
            .arg(&input)
            .arg(&reference);
        command
    });
    let value: serde_json::Value = serde_json::from_slice(&json.stdout).unwrap();
    assert_eq!(value["result"]["command"], "calmd");
    assert_eq!(value["result"]["summary"]["records_read"], 3);
    assert_eq!(value["result"]["summary"]["records_recalculated"], 2);

    let decoded = run({
        let mut command = binary();
        command.args(["view", "--no-pg", "-h"]).arg(&output);
        command
    });
    assert_eq!(decoded.stdout, DEFAULT);

    let piped = run_with_stdin(
        {
            let mut command = binary();
            command.args(["calmd", "--no-pg", "-"]).arg(&reference);
            command
        },
        &fs::read(input).unwrap(),
    );
    assert_eq!(piped.stdout, DEFAULT);
}

#[test]
fn failures_preserve_named_outputs_and_reject_aliases() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("output.sam");
    fs::write(&output, b"keep me").unwrap();
    let input = directory.path().join("missing-reference.sam");
    fs::write(
        &input,
        b"@SQ\tSN:missing\tLN:4\nread\t0\tmissing\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
    )
    .unwrap();
    let failed = binary()
        .args(["calmd", "-o"])
        .arg(&output)
        .arg(&input)
        .arg(fixture("reference.fa"))
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert_eq!(fs::read(&output).unwrap(), b"keep me");

    let source = fixture("records.sam");
    let alias = binary()
        .args(["calmd", "-o"])
        .arg(&source)
        .arg(&source)
        .arg(fixture("reference.fa"))
        .output()
        .unwrap();
    assert!(!alias.status.success());

    let json_without_output = binary()
        .arg("--json")
        .arg("calmd")
        .arg(&source)
        .arg(fixture("reference.fa"))
        .output()
        .unwrap();
    assert!(!json_without_output.status.success());
}

#[test]
fn complex_cigars_tag_replacement_and_missing_sequence_are_stable() {
    let directory = tempfile::tempdir().unwrap();
    let (input, reference) = write_complex_fixture(directory.path());
    let sam = run({
        let mut command = binary();
        command
            .args(["calmd", "--no-pg"])
            .arg(&input)
            .arg(&reference);
        command
    });
    assert_eq!(sam.stdout, COMPLEX_EXPECTED);
    assert!(
        String::from_utf8_lossy(&sam.stderr)
            .contains("1 mapped records were preserved because they have no query sequence")
    );

    let bam_input = directory.path().join("complex.bam");
    run({
        let mut command = binary();
        command
            .args(["view", "--no-pg", "-b", "-o"])
            .arg(&bam_input)
            .arg(&input);
        command
    });
    let bam_output = directory.path().join("complex.calmd.bam");
    run({
        let mut command = binary();
        command
            .args(["calmd", "--no-pg", "-u", "-o"])
            .arg(&bam_output)
            .arg(&bam_input)
            .arg(&reference);
        command
    });
    let decoded = run({
        let mut command = binary();
        command.args(["view", "--no-pg", "-h"]).arg(&bam_output);
        command
    });
    assert_eq!(decoded.stdout, COMPLEX_EXPECTED);
}
