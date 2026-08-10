use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/records.sam")
}

fn run(arguments: &[&str]) -> Output {
    binary().args(arguments).output().unwrap()
}

fn run_stdin(arguments: &[&str], input: &[u8]) -> Output {
    let mut command = binary();
    command
        .args(arguments)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let mut child = command.spawn().unwrap();
    child.stdin.take().unwrap().write_all(input).unwrap();
    child.wait_with_output().unwrap()
}

fn lines(output: &Output) -> Vec<&str> {
    std::str::from_utf8(&output.stdout)
        .unwrap()
        .lines()
        .collect()
}

fn tags(output: &Output) -> Vec<Option<&str>> {
    lines(output)
        .into_iter()
        .filter(|line| !line.starts_with('@'))
        .map(|line| {
            line.split('\t')
                .skip(11)
                .find_map(|field| field.strip_prefix("RG:Z:"))
        })
        .collect()
}

#[test]
fn overwrite_all_replaces_header_and_every_record_tag() {
    let input = fixture();
    let output = run(&[
        "addreplacerg",
        "-r",
        "ID:new",
        "-r",
        "SM:sample-b",
        "--no-PG",
        input.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::str::from_utf8(&output.stdout).unwrap();
    assert!(text.contains("@RG\tID:new\tSM:sample-b\n"), "{text}");
    assert!(!text.contains("@RG\tID:rg1\t"), "{text}");
    assert_eq!(tags(&output), [Some("new"), Some("new"), Some("new")]);
}

#[test]
fn orphan_only_preserves_existing_tags_and_stamps_missing_records() {
    let input = fixture();
    let output = run(&[
        "addreplacerg",
        "-r",
        "ID:new\\tSM:sample-b",
        "-m",
        "orphan_only",
        "--no-PG",
        input.to_str().unwrap(),
    ]);
    assert!(output.status.success());
    let text = std::str::from_utf8(&output.stdout).unwrap();
    assert!(text.contains("@RG\tID:rg1\tSM:sample-a\tLB:lib-a\n"));
    assert!(text.contains("@RG\tID:new\tSM:sample-b\n"));
    assert_eq!(tags(&output), [Some("rg1"), Some("rg1"), Some("new")]);
}

#[test]
fn omitted_source_uses_the_first_header_read_group() {
    let input = fixture();
    let output = run(&["addreplacerg", "--no-PG", input.to_str().unwrap()]);
    assert!(output.status.success());
    assert_eq!(tags(&output), [Some("rg1"), Some("rg1"), Some("rg1")]);
}

#[test]
fn same_id_requires_overwrite_and_failed_output_is_transactional() {
    let directory = tempfile::tempdir().unwrap();
    let destination = directory.path().join("output.sam");
    fs::write(&destination, b"keep me").unwrap();
    let input = fixture();
    let failed = run(&[
        "addreplacerg",
        "-r",
        "ID:rg1\\tSM:replacement",
        "-o",
        destination.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert!(!failed.status.success());
    assert_eq!(fs::read(&destination).unwrap(), b"keep me");

    let replaced = run(&[
        "addreplacerg",
        "-w",
        "-r",
        "ID:rg1\\tSM:replacement",
        "--no-PG",
        "-o",
        destination.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert!(replaced.status.success());
    let text = fs::read_to_string(destination).unwrap();
    assert!(text.contains("@RG\tID:rg1\tSM:replacement\n"));
    assert!(!text.contains("sample-a"));
}

#[test]
fn named_bam_and_standard_input_use_the_same_contract() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("output.bam");
    let input = fixture();
    let written = run(&[
        "addreplacerg",
        "-R",
        "rg1",
        "--no-PG",
        "-u",
        "-o",
        bam.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert!(written.status.success());
    assert!(fs::read(&bam).unwrap().starts_with(&[0x1f, 0x8b]));

    let decoded = run(&["view", bam.to_str().unwrap()]);
    assert!(decoded.status.success());
    assert_eq!(tags(&decoded), [Some("rg1"), Some("rg1"), Some("rg1")]);

    let piped = run_stdin(
        &["addreplacerg", "-R", "rg1", "--no-PG", "-"],
        &fs::read(input).unwrap(),
    );
    assert!(piped.status.success());
    assert_eq!(tags(&piped), [Some("rg1"), Some("rg1"), Some("rg1")]);
}

#[test]
fn invalid_source_and_output_alias_fail_loudly() {
    let input = fixture();
    let missing = run(&["addreplacerg", "-R", "missing", input.to_str().unwrap()]);
    assert!(!missing.status.success());

    let malformed = run(&["addreplacerg", "-r", "SM:no-id", input.to_str().unwrap()]);
    assert!(!malformed.status.success());

    let alias = run(&[
        "addreplacerg",
        "-o",
        input.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert!(!alias.status.success());
}

#[test]
fn json_summary_uses_the_shared_output_envelope() {
    let directory = tempfile::tempdir().unwrap();
    let destination = directory.path().join("output.bam");
    let input = fixture();
    let output = run(&[
        "--json",
        "addreplacerg",
        "-r",
        "ID:new",
        "--no-PG",
        "-o",
        destination.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["result"]["command"], "addreplacerg");
    assert_eq!(value["result"]["summary"]["records_read"], 3);
    assert_eq!(value["result"]["summary"]["records_modified"], 3);
    assert_eq!(value["result"]["summary"]["records_preserved"], 0);
}

#[test]
fn record_failure_does_not_replace_named_output() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("malformed.sam");
    let destination = directory.path().join("output.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\n@RG\tID:old\nbroken\t0\t*\t0\t0\t*\t*\t0\t0\tA\n",
    )
    .unwrap();
    fs::write(&destination, b"keep me").unwrap();
    let output = run(&[
        "addreplacerg",
        "-R",
        "old",
        "-o",
        destination.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert_eq!(fs::read(destination).unwrap(), b"keep me");
}
