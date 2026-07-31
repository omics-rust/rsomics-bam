use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn run(mut command: Command) -> Output {
    let output = command.output().expect("spawn command");
    assert!(
        output.status.success(),
        "status={:?}\nstdout={}\nstderr={}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn run_ours(arguments: &[&str]) -> Output {
    let mut command = binary();
    command.args(arguments);
    run(command)
}

fn run_samtools(arguments: &[&str]) -> Output {
    let mut command = Command::new("samtools");
    command.args(arguments);
    run(command)
}

fn assert_samtools_1_24() {
    let output = run_samtools(&["--version"]);
    let version = String::from_utf8(output.stdout).unwrap();
    assert!(version.starts_with("samtools 1.24\n"), "{version}");
}

#[test]
fn flagstat_matches_committed_samtools_text() {
    let input = golden("flagstat-small.bam");
    let output = run_ours(&["flagstat", input.to_str().unwrap()]);
    assert_eq!(
        output.stdout,
        fs::read(golden("flagstat-small.txt")).unwrap()
    );
}

#[test]
fn flagstat_accepts_sam() {
    let input = golden("records.sam");
    let output = run_ours(&["flagstat", input.to_str().unwrap()]);
    let text = String::from_utf8(output.stdout).unwrap();
    assert!(text.starts_with("3 + 0 in total"), "{text}");
    assert!(text.contains("2 + 0 mapped (66.67% : N/A)"), "{text}");
}

#[test]
fn truncated_alignment_fails_loudly() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("truncated.bam");
    let mut bytes = fs::read(golden("flagstat-small.bam")).unwrap();
    bytes.truncate(bytes.len() / 2);
    fs::write(&input, bytes).unwrap();

    let output = binary().args(["flagstat"]).arg(input).output().unwrap();
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("reading alignment record"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn machine_output_uses_the_shared_envelope() {
    let input = golden("records.sam");
    let output = run_ours(&["--json", "flagstat", input.to_str().unwrap()]);
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["status"], "ok");
    assert_eq!(value["tool"], "rsomics-bam");
    assert_eq!(value["result"]["command"], "flagstat");
    assert_eq!(
        value["result"]["counts"]["total"],
        serde_json::json!([3, 0])
    );
}

#[test]
fn domain_and_envelope_json_are_not_mixed() {
    let input = golden("records.sam");
    let output = binary()
        .args(["--json", "flagstat", "--output-fmt", "json"])
        .arg(input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stderr).unwrap();
    assert_eq!(value["status"], "error");
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn flags_match_samtools_1_24() {
    assert_samtools_1_24();
    for value in [
        "0",
        "16",
        "020",
        "0x10",
        "paired",
        "paired,read1",
        "SECONDARY,SUPPLEMENTARY",
        "4096",
    ] {
        let ours = run_ours(&["flags", value]);
        let oracle = run_samtools(&["flags", value]);
        assert_eq!(ours.stdout, oracle.stdout, "{value}");
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn flagstat_matches_samtools_1_24_for_sam_bam_and_cram() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let sam = golden("records.sam");
    let reference = directory.path().join("reference.fa");
    let bam = directory.path().join("records.bam");
    let cram = directory.path().join("records.cram");
    fs::copy(golden("reference.fa"), &reference).unwrap();

    run_samtools(&["faidx", reference.to_str().unwrap()]);
    run_samtools(&[
        "view",
        "-b",
        "-o",
        bam.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]);
    run_samtools(&[
        "view",
        "-C",
        "-T",
        reference.to_str().unwrap(),
        "-o",
        cram.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]);

    for input in [&sam, &bam, &cram] {
        for format in ["text", "tsv", "json"] {
            let mut our_arguments =
                vec!["flagstat", "--output-fmt", format, input.to_str().unwrap()];
            if input == &cram {
                our_arguments.splice(1..1, ["--reference", reference.to_str().unwrap()]);
            }
            let ours = run_ours(&our_arguments);
            let oracle =
                run_samtools(&["flagstat", "--output-fmt", format, input.to_str().unwrap()]);
            if format == "json" {
                let ours: serde_json::Value = serde_json::from_slice(&ours.stdout).unwrap();
                let oracle: serde_json::Value = serde_json::from_slice(&oracle.stdout).unwrap();
                assert_eq!(ours, oracle, "{} {format}", input.display());
            } else {
                assert_eq!(ours.stdout, oracle.stdout, "{} {format}", input.display());
            }
        }
    }
}
