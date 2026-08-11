use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/upstream/samtools-consensus")
}

fn run(program: &Path, arguments: &[String]) -> Output {
    Command::new(program).args(arguments).output().unwrap()
}

fn assert_case(input: &Path, arguments: &[&str]) {
    let ours = Path::new(env!("CARGO_BIN_EXE_rsomics-bam"));
    let samtools = Path::new("samtools");
    let mut ours_arguments = vec!["consensus".to_owned()];
    ours_arguments.extend(arguments.iter().map(|value| (*value).to_owned()));
    ours_arguments.push(input.display().to_string());
    let samtools_arguments = ours_arguments.clone();
    let ours_output = run(ours, &ours_arguments);
    let samtools_output = run(samtools, &samtools_arguments);

    assert!(
        ours_output.status.success(),
        "{}",
        String::from_utf8_lossy(&ours_output.stderr)
    );
    assert!(
        samtools_output.status.success(),
        "{}",
        String::from_utf8_lossy(&samtools_output.stderr)
    );
    assert_eq!(ours_output.stdout, samtools_output.stdout, "{arguments:?}");
}

fn assert_samtools_1_24() {
    let output = Command::new("samtools").arg("--version").output().unwrap();
    assert!(output.status.success());
    assert!(String::from_utf8_lossy(&output.stdout).starts_with("samtools 1.24"));
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn modes_profiles_and_calibrations_match_samtools_1_24() {
    assert_samtools_1_24();
    let input = root().join("consen1c.sam");

    for arguments in [
        &[][..],
        &["--mode", "simple", "--ambig", "--format", "pileup"][..],
        &[
            "--mode",
            "bayesian_116",
            "--format",
            "fastq",
            "--show-del",
            "yes",
        ][..],
        &[
            "--mode",
            "simple",
            "--call-fract",
            "0.6",
            "--format",
            "fastq",
            "--mark-ins",
        ][..],
    ] {
        assert_case(&input, arguments);
    }
    for profile in ["hifi", "hiseq", "r10.4_sup", "r10.4_dup", "ultima"] {
        assert_case(
            &input,
            &[
                "--config",
                profile,
                "--format",
                "pileup",
                "--show-del",
                "yes",
            ],
        );
    }
    for calibration in [
        ":flat",
        ":hifi",
        ":hiseq",
        ":r10.4_sup",
        ":r10.4_dup",
        ":ultima",
    ] {
        assert_case(
            &input,
            &["--qual-calibration", calibration, "--format", "pileup"],
        );
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn indexed_regions_and_reference_fill_match_samtools_1_24() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("consen4.bam");
    let status = Command::new("samtools")
        .args(["view", "--write-index"])
        .arg(root().join("consen4.sam"))
        .arg("-o")
        .arg(&bam)
        .status()
        .unwrap();
    assert!(status.success());

    assert_case(&bam, &["--region", "c1:2-9", "--format", "fastq", "-a"]);
    let bed = root().join("consen4.bed");
    assert_case(
        &bam,
        &[
            "--regions-file",
            bed.to_str().unwrap(),
            "--format",
            "pileup",
            "-a",
        ],
    );

    let input = root().join("consen1c.sam");
    let reference = root().join("consen1c.fa");
    assert_case(
        &input,
        &[
            "--reference",
            reference.to_str().unwrap(),
            "--ref-qual",
            "20",
            "--format",
            "fastq",
            "--show-del",
            "yes",
            "--show-ins",
            "no",
            "-aa",
        ],
    );
}
