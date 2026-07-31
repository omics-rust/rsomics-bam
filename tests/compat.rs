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

struct AlignmentSet {
    sam: PathBuf,
    bam: PathBuf,
    cram: PathBuf,
    reference: PathBuf,
}

fn build_alignment_set(directory: &Path) -> AlignmentSet {
    let sam = golden("records.sam");
    let reference = directory.join("reference.fa");
    let bam = directory.join("records.bam");
    let cram = directory.join("records.cram");
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

    AlignmentSet {
        sam,
        bam,
        cram,
        reference,
    }
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
fn head_defaults_to_header_only() {
    let input = golden("records.sam");
    let output = run_ours(&["head", input.to_str().unwrap()]);
    let text = String::from_utf8(output.stdout).unwrap();
    assert_eq!(
        text,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:40\n@RG\tID:rg1\tSM:sample-a\n"
    );
}

#[test]
fn head_limits_headers_and_records() {
    let input = golden("records.sam");
    let output = run_ours(&["head", "-H", "1", "-n", "2", input.to_str().unwrap()]);
    let text = String::from_utf8(output.stdout).unwrap();
    assert_eq!(text.lines().count(), 3, "{text}");
    assert!(text.starts_with("@HD\tVN:1.6\tSO:coordinate\n"), "{text}");
    assert!(text.contains("read1\t99\tchr1\t1\t60\t8M"), "{text}");
    assert!(text.contains("read1\t147\tchr1\t17\t60\t8M"), "{text}");
}

#[test]
fn head_rejects_json_stream_mixing() {
    let input = golden("records.sam");
    let output = binary()
        .args(["--json", "head"])
        .arg(input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stderr).unwrap();
    assert_eq!(value["status"], "error");
}

#[test]
fn quickcheck_accepts_sam_and_bam() {
    for input in [golden("records.sam"), golden("flagstat-small.bam")] {
        let output = run_ours(&["quickcheck", input.to_str().unwrap()]);
        assert!(output.stdout.is_empty());
        assert!(output.stderr.is_empty());
    }
}

#[test]
fn quickcheck_detects_missing_bam_eof() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("missing-eof.bam");
    let mut bytes = fs::read(golden("flagstat-small.bam")).unwrap();
    bytes.truncate(bytes.len() - 28);
    fs::write(&input, bytes).unwrap();

    let output = binary()
        .args(["quickcheck", "-v"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        format!("{}\n", input.display())
    );
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("was missing EOF block when one should be present"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn quickcheck_unmapped_allows_header_without_targets() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("unmapped.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\nread1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n",
    )
    .unwrap();

    let rejected = binary().args(["quickcheck"]).arg(&input).output().unwrap();
    assert!(!rejected.status.success());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("had no targets in header"));

    let accepted = run_ours(&["quickcheck", "-u", input.to_str().unwrap()]);
    assert!(accepted.stdout.is_empty());
    assert!(accepted.stderr.is_empty());
}

#[test]
fn quickcheck_quiet_suppresses_per_input_warning() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("empty");
    fs::write(&input, []).unwrap();

    let output = binary()
        .args(["quickcheck", "-q"])
        .arg(input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(
        !stderr.contains("could not be opened for reading"),
        "{stderr}"
    );
    assert!(
        stderr.contains("quickcheck failed for 1 of 1 inputs"),
        "{stderr}"
    );
}

#[test]
fn quickcheck_uses_the_shared_json_envelope() {
    let valid = run_ours(&[
        "--json",
        "quickcheck",
        golden("records.sam").to_str().unwrap(),
    ]);
    let value: serde_json::Value = serde_json::from_slice(&valid.stdout).unwrap();
    assert_eq!(value["status"], "ok");
    assert_eq!(value["result"]["command"], "quickcheck");
    assert_eq!(
        value["result"]["report"]["files"][0]["problems"],
        serde_json::json!([])
    );

    let directory = tempfile::tempdir().unwrap();
    let invalid = directory.path().join("invalid");
    fs::write(&invalid, b"not alignment data").unwrap();
    let failed = binary()
        .args(["--json", "quickcheck"])
        .arg(invalid)
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert!(failed.stdout.is_empty());
    let value: serde_json::Value = serde_json::from_slice(&failed.stderr).unwrap();
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
    let inputs = build_alignment_set(directory.path());

    for input in [&inputs.sam, &inputs.bam, &inputs.cram] {
        for format in ["text", "tsv", "json"] {
            let mut our_arguments =
                vec!["flagstat", "--output-fmt", format, input.to_str().unwrap()];
            if input == &inputs.cram {
                our_arguments.splice(1..1, ["--reference", inputs.reference.to_str().unwrap()]);
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

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn head_matches_samtools_1_24_for_sam_bam_and_cram() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());

    for input in [&inputs.sam, &inputs.bam, &inputs.cram] {
        for (our_options, oracle_options) in [
            (vec![], vec![]),
            (vec!["-H", "1"], vec!["-h", "1"]),
            (vec!["-H", "2", "-n", "2"], vec!["-h", "2", "-n", "2"]),
            (vec!["-H", "0", "-n", "3"], vec!["-h", "0", "-n", "3"]),
        ] {
            let mut our_arguments = vec!["head"];
            our_arguments.extend(our_options);
            if input == &inputs.cram {
                our_arguments.extend(["--reference", inputs.reference.to_str().unwrap()]);
            }
            our_arguments.push(input.to_str().unwrap());

            let mut oracle_arguments = vec!["head"];
            oracle_arguments.extend(oracle_options);
            if input == &inputs.cram {
                oracle_arguments.extend(["-T", inputs.reference.to_str().unwrap()]);
            }
            oracle_arguments.push(input.to_str().unwrap());

            let ours = run_ours(&our_arguments);
            let oracle = run_samtools(&oracle_arguments);
            assert_eq!(ours.stdout, oracle.stdout, "{}", input.display());
        }
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn quickcheck_matches_samtools_1_24_decisions() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());
    let no_targets = directory.path().join("unmapped.sam");
    fs::write(
        &no_targets,
        b"@HD\tVN:1.6\nread1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n",
    )
    .unwrap();
    let missing_eof = directory.path().join("missing-eof.bam");
    let mut bytes = fs::read(&inputs.bam).unwrap();
    bytes.truncate(bytes.len() - 28);
    fs::write(&missing_eof, bytes).unwrap();
    let invalid = directory.path().join("invalid");
    fs::write(&invalid, b"not alignment data").unwrap();

    for input in [
        &inputs.sam,
        &inputs.bam,
        &inputs.cram,
        &no_targets,
        &missing_eof,
        &invalid,
    ] {
        for options in [Vec::new(), vec!["-u"], vec!["-q"], vec!["-q", "-v"]] {
            let mut ours = binary();
            ours.arg("quickcheck").args(&options).arg(input);
            let ours = ours.output().unwrap();

            let mut oracle = Command::new("samtools");
            oracle.arg("quickcheck").args(&options).arg(input);
            let oracle = oracle.output().unwrap();

            assert_eq!(
                ours.status.success(),
                oracle.status.success(),
                "{} {options:?}\nours={}\noracle={}",
                input.display(),
                String::from_utf8_lossy(&ours.stderr),
                String::from_utf8_lossy(&oracle.stderr)
            );
            if options.contains(&"-v") {
                assert_eq!(ours.stdout, oracle.stdout, "{}", input.display());
            }
        }
    }
}
