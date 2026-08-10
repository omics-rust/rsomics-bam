use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

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

fn run_stdin(mut command: Command, input: &[u8]) -> Output {
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

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn samtools() -> Command {
    Command::new("samtools")
}

fn build_inputs(directory: &Path) -> (PathBuf, PathBuf, PathBuf, PathBuf) {
    let reference = directory.join("reference.fa");
    fs::write(
        &reference,
        b">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
    )
    .unwrap();
    run({
        let mut command = samtools();
        command.args(["faidx"]).arg(&reference);
        command
    });

    let sam = directory.join("records.sam");
    fs::write(
        &sam,
        b"@HD\tVN:1.6\tSO:queryname\n\
@SQ\tSN:chr1\tLN:60\n\
q1\t65\tchr1\t1\t60\t4M\t=\t9\t12\tACGT\tIIII\n\
q1\t145\tchr1\t9\t60\t4M\t=\t1\t-12\tAGTC\tABCD\n\
q2\t69\t*\t0\t0\t*\t*\t0\t0\tAAAA\t*\n\
q2\t69\t*\t0\t0\t*\t*\t0\t0\tCCCC\tJJJJ\tOQ:Z:KLMN\n\
q3\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n\
q4\t2052\t*\t0\t0\t*\t*\t0\t0\tTTTT\tIIII\n",
    )
    .unwrap();
    let bam = directory.join("records.bam");
    run({
        let mut command = samtools();
        command.args(["view", "-b", "-o"]).arg(&bam).arg(&sam);
        command
    });
    let cram = directory.join("records.cram");
    run({
        let mut command = samtools();
        command
            .args(["view", "-C", "-T"])
            .arg(&reference)
            .args(["-o"])
            .arg(&cram)
            .arg(&sam);
        command
    });
    (sam, bam, cram, reference)
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn fasta_and_fastq_match_samtools_1_24_for_sam_bam_cram_and_stdin() {
    let version = run({
        let mut command = samtools();
        command.arg("--version");
        command
    });
    assert!(version.stdout.starts_with(b"samtools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let (sam, bam, cram, reference) = build_inputs(directory.path());
    for input in [&sam, &bam, &cram] {
        for options in [
            vec!["fasta"],
            vec!["fasta", "-n"],
            vec!["fastq"],
            vec!["fastq", "-n", "-O", "-v", "7"],
            vec!["fastq", "-F", "0"],
        ] {
            let mut ours = binary();
            ours.args(&options);
            let mut oracle = samtools();
            oracle.args(&options);
            if input == &cram {
                ours.args(["--reference"]).arg(&reference);
                oracle.args(["--reference"]).arg(&reference);
            }
            ours.arg(input);
            oracle.arg(input);
            assert_eq!(
                run(ours).stdout,
                run(oracle).stdout,
                "{} {options:?}",
                input.display()
            );
        }
    }

    let bytes = fs::read(&sam).unwrap();
    for command in ["fasta", "fastq"] {
        assert_eq!(
            run_stdin(
                {
                    let mut ours = binary();
                    ours.arg(command);
                    ours
                },
                &bytes
            )
            .stdout,
            run_stdin(
                {
                    let mut oracle = samtools();
                    oracle.arg(command);
                    oracle
                },
                &bytes
            )
            .stdout
        );
    }
}
