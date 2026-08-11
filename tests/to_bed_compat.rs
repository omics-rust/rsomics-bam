use std::fs;
use std::path::Path;
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

fn encode_bam(sam: &Path, bam: &Path) {
    run({
        let mut command = ours();
        command
            .args(["view", "-b", "--no-pg", "-o"])
            .arg(bam)
            .arg(sam);
        command
    });
}

fn assert_mode(bam: &Path, ours_args: &[&str], bedtools_args: &[&str]) {
    let ours = run({
        let mut command = ours();
        command.arg("to-bed").args(ours_args).arg(bam);
        command
    });
    let oracle = run({
        let mut command = Command::new("bedtools");
        command
            .arg("bamtobed")
            .arg("-i")
            .arg(bam)
            .args(bedtools_args);
        command
    });
    assert_eq!(ours.stdout, oracle.stdout, "mode {ours_args:?}");
}

#[test]
#[ignore = "requires bedtools 2.31.1"]
fn documented_bed6_and_bed12_modes_match_bedtools_2311() {
    let version = run({
        let mut command = Command::new("bedtools");
        command.arg("--version");
        command
    });
    assert_eq!(
        String::from_utf8_lossy(&version.stdout).trim(),
        "bedtools v2.31.1"
    );

    let directory = tempfile::tempdir().unwrap();
    let sam = directory.path().join("records.sam");
    let bam = directory.path().join("records.bam");
    fs::write(
        &sam,
        b"@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:1000\n\
forward\t65\tchr1\t5\t60\t2S3M1I2=1X4D5N4M2S\t=\t50\t0\tAACCGGTTAAACCCC\tIIIIIIIIIIIIIII\tNM:i:6\tXI:i:9\n\
reverse\t145\tchr1\t50\t37\t3M2D4M\t=\t5\t0\tAAAAAAA\tIIIIIII\tNM:i:2\tXI:i:4\n\
unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tAA\tII\n",
    )
    .unwrap();
    encode_bam(&sam, &bam);

    for (ours_args, bedtools_args) in [
        (&[][..], &[][..]),
        (&["--split"][..], &["-split"][..]),
        (&["--split-d"][..], &["-splitD"][..]),
        (&["--bed12"][..], &["-bed12"][..]),
        (&["--bed12", "--split-d"][..], &["-bed12", "-splitD"][..]),
        (
            &["--bed12", "--color", "1,2,3"][..],
            &["-bed12", "-color", "1,2,3"][..],
        ),
        (&["--ed"][..], &["-ed"][..]),
        (&["--tag", "XI"][..], &["-tag", "XI"][..]),
        (&["--cigar"][..], &["-cigar"][..]),
    ] {
        assert_mode(&bam, ours_args, bedtools_args);
    }
}

#[test]
#[ignore = "requires bedtools 2.31.1"]
fn documented_bedpe_modes_match_bedtools_2311() {
    let directory = tempfile::tempdir().unwrap();
    let sam = directory.path().join("pairs.sam");
    let bam = directory.path().join("pairs.bam");
    fs::write(
        &sam,
        b"@HD\tVN:1.6\tSO:queryname\n\
@SQ\tSN:chr1\tLN:1000\n\
@SQ\tSN:chr2\tLN:1000\n\
mapped\t145\tchr1\t10\t30\t5M\tchr2\t50\t0\tAAAAA\tIIIII\tNM:i:2\n\
mapped\t65\tchr2\t50\t40\t5M\tchr1\t10\t0\tCCCCC\tIIIII\tNM:i:3\n\
half\t73\tchr1\t80\t55\t5M\t*\t0\t0\tGGGGG\tIIIII\tNM:i:4\n\
half\t133\t*\t0\t0\t*\tchr1\t80\t0\tTTTTT\tIIIII\n",
    )
    .unwrap();
    encode_bam(&sam, &bam);

    for (ours_args, bedtools_args) in [
        (&["--bedpe"][..], &["-bedpe"][..]),
        (&["--bedpe", "--mate1"][..], &["-bedpe", "-mate1"][..]),
        (&["--bedpe", "--ed"][..], &["-bedpe", "-ed"][..]),
    ] {
        assert_mode(&bam, ours_args, bedtools_args);
    }
}
