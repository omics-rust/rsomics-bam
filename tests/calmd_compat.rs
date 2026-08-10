use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn ours() -> Command {
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

fn write_fixture(directory: &Path) -> (PathBuf, PathBuf) {
    let input = directory.join("input.sam");
    let reference = directory.join("reference.fa");
    fs::write(
        &input,
        b"@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:64\n\
perfect\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n\
wrong\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACGAAC\tIIIIIIIIII\tAA:Z:first\tMD:Z:10\tBB:i:2\tNM:i:0\tCC:Z:last\n\
deletion\t0\tchr1\t1\t60\t5M2D5M\t*\t0\t0\tACGTATACGT\tIIIIIIIIII\n\
insertion\t0\tchr1\t1\t60\t5M3I5M\t*\t0\t0\tACGTAGGGCGTAC\tIIIIIIIIIIIII\n\
clipped\t0\tchr1\t1\t60\t3S10M\t*\t0\t0\tTTTACGTACGTAC\tIIIIIIIIIIIII\n\
skipped\t0\tchr1\t1\t60\t5M3N5M\t*\t0\t0\tACGTAACGTA\tIIIIIIIIII\n\
ambiguous\t0\tchr1\t41\t60\t10M\t*\t0\t0\tACGTNNNNAC\tIIIIIIIIII\n\
cigar-equal\t0\tchr1\t1\t60\t5=2X3=\t*\t0\t0\tACGTAAATAC\tIIIIIIIIII\n\
literal-equal\t0\tchr1\t1\t60\t4M\t*\t0\t0\t====\tIIII\n\
correct\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tMD:Z:10\tNM:i:0\n\
missing\t256\tchr1\t1\t60\t10M\t*\t0\t0\t*\t*\n\
unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n",
    )
    .unwrap();
    fs::write(
        &reference,
        b">chr1\nACGTACGTACGTACGTACGTAACCGGTTACGTACGTACGTACGTNNNNACGTACGTACGTACGT\n",
    )
    .unwrap();
    fs::write(directory.join("reference.fa.fai"), b"chr1\t64\t6\t64\t65\n").unwrap();
    (input, reference)
}

fn samtools_calmd(input: &Path, reference: &Path, use_equal: bool) -> Vec<u8> {
    let mut command = Command::new("samtools");
    command.args(["calmd", "--no-PG"]);
    if use_equal {
        command.arg("-e");
    }
    command.arg(input).arg(reference);
    run(command).stdout
}

fn rsomics_calmd(input: &Path, reference: &Path, use_equal: bool) -> Vec<u8> {
    let mut command = ours();
    command.args(["calmd", "--no-pg"]);
    if use_equal {
        command.arg("-e");
    }
    command.arg(input).arg(reference);
    run(command).stdout
}

#[test]
#[ignore = "requires samtools 1.24"]
fn sam_bam_and_cram_match_samtools_1_24() {
    let version = run({
        let mut command = Command::new("samtools");
        command.arg("--version");
        command
    });
    assert_eq!(
        String::from_utf8_lossy(&version.stdout).lines().next(),
        Some("samtools 1.24")
    );

    let directory = tempfile::tempdir().unwrap();
    let (sam, reference) = write_fixture(directory.path());
    let bam = directory.path().join("input.bam");
    run({
        let mut command = Command::new("samtools");
        command
            .args(["view", "--no-PG", "-b", "-o"])
            .arg(&bam)
            .arg(&sam);
        command
    });
    let cram = directory.path().join("input.cram");
    run({
        let mut command = Command::new("samtools");
        command
            .args(["view", "--no-PG", "-C", "-T"])
            .arg(&reference)
            .args(["-o"])
            .arg(&cram)
            .arg(&sam);
        command
    });

    for input in [&sam, &bam, &cram] {
        for use_equal in [false, true] {
            assert_eq!(
                rsomics_calmd(input, &reference, use_equal),
                samtools_calmd(input, &reference, use_equal),
                "input={} equal={use_equal}",
                input.display()
            );
        }
    }
}
