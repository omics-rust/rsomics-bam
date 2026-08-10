use std::fs;
use std::path::Path;
use std::process::{Command, Output};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn samtools(arguments: &[&str]) -> Output {
    let output = Command::new("samtools").args(arguments).output().unwrap();
    assert!(
        output.status.success(),
        "samtools stdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn ours(arguments: &[&str]) -> Output {
    let output = binary().args(arguments).output().unwrap();
    assert!(
        output.status.success(),
        "rsomics stdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

const PADDED: &[u8] = b"@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:10\n\
@SQ\tSN:chr2\tLN:6\n\
chr1\t0\tchr1\t1\t0\t3M1D1M1D4M\t*\t0\t0\tACGTACGT\t~~~~~~~~\n\
clipped\t0\tchr1\t1\t60\t2S3M1D1M1D4M2S\t*\t0\t0\tTTACGTACGTGG\tIIIIIIIIIIII\n\
leading-pad\t0\tchr1\t6\t60\t1D4M\t*\t0\t0\tACGT\tIIII\n\
skip\t0\tchr1\t1\t60\t3M1N1M1D4M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
same-mate\t1\tchr1\t1\t60\t3M\t=\t5\t0\tACG\tIII\n\
cross-mate\t1\tchr1\t1\t60\t3M\tchr2\t4\t0\tACG\tIII\n\
unmapped\t4\tchr1\t6\t0\t3M\t=\t5\t0\tACG\tIII\n\
chr2\t0\tchr2\t1\t0\t1M2D3M\t*\t0\t0\tACGT\t~~~~\n\
second\t0\tchr2\t1\t60\t1M2D3M\t*\t0\t0\tACGT\tIIII\n";

fn write_fixture(directory: &Path) -> (std::path::PathBuf, std::path::PathBuf) {
    let input = directory.join("padded.sam");
    let reference = directory.join("padded.fa");
    fs::write(&input, PADDED).unwrap();
    fs::write(&reference, b">chr1\nACG*T*ACGT\n>chr2\nA**CGT\n").unwrap();
    (input, reference)
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn sam_bam_cram_outputs_and_projection_edges_match_samtools_1_24() {
    let version = samtools(&["--version"]);
    assert!(version.stdout.starts_with(b"samtools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let (sam, reference) = write_fixture(directory.path());
    let bam = directory.path().join("padded.bam");
    let cram = directory.path().join("padded.cram");
    samtools(&[
        "view",
        "-b",
        "-o",
        bam.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]);
    samtools(&[
        "view",
        "-C",
        "--output-fmt-option",
        "no_ref=1",
        "-o",
        cram.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]);

    for input in [&sam, &bam, &cram] {
        let expected = samtools(&[
            "depad",
            "-s",
            "--no-PG",
            "-T",
            reference.to_str().unwrap(),
            input.to_str().unwrap(),
        ]);
        let actual = ours(&[
            "depad",
            "-s",
            "--no-pg",
            "-T",
            reference.to_str().unwrap(),
            input.to_str().unwrap(),
        ]);
        assert_eq!(actual.stdout, expected.stdout, "input={}", input.display());
        assert!(String::from_utf8_lossy(&actual.stderr).contains("CIGAR N"));
    }

    for compression in [None, Some("-u"), Some("-1")] {
        let suffix = compression.unwrap_or("default").trim_start_matches('-');
        let expected = directory.path().join(format!("samtools-{suffix}.bam"));
        let actual = directory.path().join(format!("rsomics-{suffix}.bam"));
        let mut upstream = vec![
            "depad",
            "--no-PG",
            "-T",
            reference.to_str().unwrap(),
            "-o",
            expected.to_str().unwrap(),
        ];
        let mut rsomics = vec![
            "depad",
            "--no-pg",
            "-S",
            "-T",
            reference.to_str().unwrap(),
            "-@",
            "2",
            "-o",
            actual.to_str().unwrap(),
        ];
        if let Some(option) = compression {
            upstream.push(option);
            rsomics.push(option);
        }
        upstream.push(bam.to_str().unwrap());
        rsomics.push(bam.to_str().unwrap());
        samtools(&upstream);
        ours(&rsomics);
        let expected = samtools(&["view", "-h", "--no-PG", expected.to_str().unwrap()]);
        let actual = samtools(&["view", "-h", "--no-PG", actual.to_str().unwrap()]);
        assert_eq!(
            actual.stdout, expected.stdout,
            "compression={compression:?}"
        );
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn long_cigar_bam_output_matches_samtools_1_24() {
    let version = samtools(&["--version"]);
    assert!(version.stdout.starts_with(b"samtools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("long.fa");
    let input_sam = directory.path().join("long.sam");
    let input_bam = directory.path().join("long.bam");
    let expected = directory.path().join("samtools.bam");
    let actual = directory.path().join("rsomics.bam");
    let padded = "A*".repeat(35_000);
    let query = "A".repeat(70_000);
    let qualities = "I".repeat(70_000);
    fs::write(&reference, format!(">chr1\n{padded}\n")).unwrap();
    fs::write(
        &input_sam,
        format!(
            "@SQ\tSN:chr1\tLN:70000\nread\t0\tchr1\t1\t60\t70000M\t*\t0\t0\t{query}\t{qualities}\n"
        ),
    )
    .unwrap();
    samtools(&[
        "view",
        "-b",
        "-o",
        input_bam.to_str().unwrap(),
        input_sam.to_str().unwrap(),
    ]);
    samtools(&[
        "depad",
        "--no-PG",
        "-T",
        reference.to_str().unwrap(),
        "-o",
        expected.to_str().unwrap(),
        input_bam.to_str().unwrap(),
    ]);
    ours(&[
        "depad",
        "--no-pg",
        "-T",
        reference.to_str().unwrap(),
        "-o",
        actual.to_str().unwrap(),
        input_bam.to_str().unwrap(),
    ]);
    let expected = samtools(&["view", "--no-PG", expected.to_str().unwrap()]);
    let actual = samtools(&["view", "--no-PG", actual.to_str().unwrap()]);
    assert_eq!(actual.stdout, expected.stdout);
}
