use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn samtools() -> Command {
    Command::new("samtools")
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

fn stable_sam(bytes: &[u8]) -> Vec<&[u8]> {
    bytes
        .split(|byte| *byte == b'\n')
        .filter(|line| !line.is_empty())
        .filter(|line| !line.starts_with(b"@CO\t") && !line.starts_with(b"@PG\t"))
        .collect()
}

fn decoded(path: &Path) -> Output {
    run({
        let mut command = samtools();
        command.args(["view", "-h"]).arg(path);
        command
    })
}

fn assert_stdout_matches(our_arguments: &[&str], oracle_arguments: &[&str]) {
    let ours = run({
        let mut command = binary();
        command.arg("import").args(our_arguments).arg("--no-PG");
        command
    });
    let oracle = run({
        let mut command = samtools();
        command.arg("import").args(oracle_arguments).arg("--no-PG");
        command
    });
    assert_eq!(stable_sam(&ours.stdout), stable_sam(&oracle.stdout));
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn import_matches_samtools_1_24_for_core_input_modes_and_tags() {
    let version = run({
        let mut command = samtools();
        command.arg("--version");
        command
    });
    assert!(version.stdout.starts_with(b"samtools 1.24\n"));

    let single = fixture("import-se.fastq");
    let iupac = fixture("import-iupac.fastq");
    let interleaved = fixture("import-interleaved.fastq");
    let read1 = fixture("import-r1.fastq");
    let read2 = fixture("import-r2.fastq");
    assert_stdout_matches(&[single.to_str().unwrap()], &[single.to_str().unwrap()]);
    assert_stdout_matches(&[iupac.to_str().unwrap()], &[iupac.to_str().unwrap()]);
    assert_stdout_matches(
        &["-0", interleaved.to_str().unwrap(), "--order", "ro"],
        &["-0", interleaved.to_str().unwrap(), "--order", "ro"],
    );
    assert_stdout_matches(
        &["-s", interleaved.to_str().unwrap()],
        &["-s", interleaved.to_str().unwrap()],
    );
    assert_stdout_matches(
        &[
            read1.to_str().unwrap(),
            read2.to_str().unwrap(),
            "-r",
            "ID:lib1",
            "-r",
            "SM:sample",
            "--order",
            "ro:6",
        ],
        &[
            read1.to_str().unwrap(),
            read2.to_str().unwrap(),
            "-r",
            "ID:lib1",
            "-r",
            "SM:sample",
            "--order",
            "ro:6",
        ],
    );
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn import_bam_gzip_and_stdin_match_samtools_1_24() {
    let version = run({
        let mut command = samtools();
        command.arg("--version");
        command
    });
    assert!(version.stdout.starts_with(b"samtools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let read1 = fixture("import-r1.fastq");
    let read2 = fixture("import-r2.fastq");
    let ours_bam = directory.path().join("ours.bam");
    let oracle_bam = directory.path().join("oracle.bam");
    run({
        let mut command = binary();
        command
            .args(["import", "-1"])
            .arg(&read1)
            .args(["-2"])
            .arg(&read2)
            .args(["-R", "lib1", "--order", "ro", "--no-PG", "-o"])
            .arg(&ours_bam);
        command
    });
    run({
        let mut command = samtools();
        command
            .args(["import", "-1"])
            .arg(&read1)
            .args(["-2"])
            .arg(&read2)
            .args(["-R", "lib1", "--order", "ro", "--no-PG", "-o"])
            .arg(&oracle_bam);
        command
    });
    assert_eq!(
        stable_sam(&decoded(&ours_bam).stdout),
        stable_sam(&decoded(&oracle_bam).stdout)
    );

    let plain = fs::read(fixture("import-se.fastq")).unwrap();
    let gzip = directory.path().join("reads.fastq.gz");
    let mut encoder = flate2::write::GzEncoder::new(
        fs::File::create(&gzip).unwrap(),
        flate2::Compression::default(),
    );
    encoder.write_all(&plain).unwrap();
    encoder.finish().unwrap();
    assert_stdout_matches(&[gzip.to_str().unwrap()], &[gzip.to_str().unwrap()]);

    let ours = run_stdin(
        {
            let mut command = binary();
            command.args(["import", "-", "--no-PG"]);
            command
        },
        &plain,
    );
    let oracle = run_stdin(
        {
            let mut command = samtools();
            command.args(["import", "-", "--no-PG"]);
            command
        },
        &plain,
    );
    assert_eq!(stable_sam(&ours.stdout), stable_sam(&oracle.stdout));
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn import_rejects_non_iupac_fastq_like_samtools_1_24() {
    for name in ["import-unknown.fastq", "import-equals.fastq"] {
        let input = fixture(name);
        let ours = binary()
            .args(["import", input.to_str().unwrap(), "--no-PG"])
            .output()
            .unwrap();
        let oracle = samtools()
            .args(["import", input.to_str().unwrap(), "--no-PG"])
            .output()
            .unwrap();
        assert!(!ours.status.success(), "rsomics accepted {name}");
        assert!(!oracle.status.success(), "samtools accepted {name}");
    }
}
