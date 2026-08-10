use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden/markdup")
        .join(name)
}

fn samtools_available() -> bool {
    Command::new("samtools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|status| status.success())
}

fn require_both_fail(name: &str) {
    if !samtools_available() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();
    let input = fixture(name);
    let ours = binary()
        .args(["markdup", "--no-PG"])
        .arg(&input)
        .arg("-o")
        .arg(directory.path().join("ours.bam"))
        .output()
        .unwrap();
    let oracle = Command::new("samtools")
        .args(["markdup", "--no-PG"])
        .arg(input)
        .arg(directory.path().join("oracle.bam"))
        .output()
        .unwrap();
    assert!(!ours.status.success(), "rsomics-bam accepted {name}");
    assert!(!oracle.status.success(), "samtools accepted {name}");
}

fn require_equal(name: &str, remove: bool) {
    let mut arguments = Vec::new();
    if remove {
        arguments.push("-r");
    }
    require_equal_args(name, &arguments);
}

fn require_equal_args(name: &str, arguments: &[&str]) {
    if !samtools_available() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();
    let input = fixture(name);
    require_equal_path(&input, arguments, directory.path());
}

fn require_equal_path(input: &Path, arguments: &[&str], directory: &Path) {
    let ours_path = directory.join("ours.bam");
    let oracle_path = directory.join("oracle.bam");
    let mut ours = binary();
    ours.args(["markdup", "--no-PG"]);
    ours.args(arguments);
    let ours = ours.arg(input).arg("-o").arg(&ours_path).output().unwrap();
    assert!(
        ours.status.success(),
        "{}",
        String::from_utf8_lossy(&ours.stderr)
    );

    let mut oracle = Command::new("samtools");
    oracle.args(["markdup", "--no-PG"]);
    oracle.args(arguments);
    let oracle = oracle.arg(input).arg(&oracle_path).output().unwrap();
    assert!(
        oracle.status.success(),
        "{}",
        String::from_utf8_lossy(&oracle.stderr)
    );

    let ours = sam_text(&ours_path);
    let oracle = sam_text(&oracle_path);
    assert_sam_eq(&ours, &oracle);
}

fn assert_sam_eq(ours: &[u8], oracle: &[u8]) {
    let ours = std::str::from_utf8(ours).unwrap();
    let oracle = std::str::from_utf8(oracle).unwrap();
    let ours_lines: Vec<_> = ours.lines().collect();
    let oracle_lines: Vec<_> = oracle.lines().collect();
    assert_eq!(ours_lines.len(), oracle_lines.len(), "SAM line count");
    for (line_index, (ours, oracle)) in ours_lines.iter().zip(&oracle_lines).enumerate() {
        if ours == oracle {
            continue;
        }
        let ours_fields: Vec<_> = ours.split('\t').collect();
        let oracle_fields: Vec<_> = oracle.split('\t').collect();
        if ours_fields.len() != oracle_fields.len() {
            let ours_tags: Vec<_> = ours_fields
                .iter()
                .skip(11)
                .map(|field| &field[..2])
                .collect();
            let oracle_tags: Vec<_> = oracle_fields
                .iter()
                .skip(11)
                .map(|field| &field[..2])
                .collect();
            panic!(
                "SAM line {} tags differ: ours={ours_tags:?} oracle={oracle_tags:?}",
                line_index + 1
            );
        }
        for (field_index, (ours, oracle)) in ours_fields.iter().zip(&oracle_fields).enumerate() {
            assert_eq!(
                ours,
                oracle,
                "SAM line {} field {}",
                line_index + 1,
                field_index + 1
            );
        }
    }
}

fn sam_text(path: &Path) -> Vec<u8> {
    let output = Command::new("samtools")
        .args(["view", "--no-PG", "-h"])
        .arg(path)
        .output()
        .unwrap();
    assert!(output.status.success());
    output.stdout
}

#[test]
fn query_name_order_failure_matches_samtools() {
    require_both_fail("1_name_sort.sam");
}

#[test]
fn observed_bad_order_failure_matches_samtools() {
    require_both_fail("2_bad_order.sam");
}

#[test]
fn missing_mc_failure_matches_samtools() {
    require_both_fail("3_missing_mc.sam");
}

#[test]
fn missing_ms_failure_matches_samtools() {
    require_both_fail("4_missing_ms.sam");
}

#[test]
fn default_output_matches_samtools() {
    require_equal("5_markdup.sam", false);
}

#[test]
fn removal_output_matches_samtools() {
    require_equal("6_remove_dups.sam", true);
}

#[test]
fn sequence_mode_matches_samtools() {
    require_equal_args("8_sequence_mode.sam", &["--mode", "s"]);
}

#[test]
fn include_fails_matches_samtools() {
    require_equal_args("9_include_fails.sam", &["--mode", "s", "--include-fails"]);
}

#[test]
fn bam_input_with_threads_matches_samtools() {
    if !samtools_available() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.bam");
    let converted = Command::new("samtools")
        .args(["view", "-b", "-o"])
        .arg(&input)
        .arg(fixture("5_markdup.sam"))
        .status()
        .unwrap();
    assert!(converted.success());
    require_equal_path(&input, &["-@", "2"], directory.path());
}

#[test]
fn cram_input_matches_samtools() {
    if !samtools_available() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    std::fs::write(
        &reference,
        format!(">contig_000000000\n{}\n", "N".repeat(11_391)),
    )
    .unwrap();
    assert!(
        Command::new("samtools")
            .arg("faidx")
            .arg(&reference)
            .status()
            .unwrap()
            .success()
    );
    let input = directory.path().join("input.cram");
    let converted = Command::new("samtools")
        .args(["view", "-C", "-T"])
        .arg(&reference)
        .args(["-o", input.to_str().unwrap()])
        .arg(fixture("5_markdup.sam"))
        .status()
        .unwrap();
    assert!(converted.success());
    require_equal_path(
        &input,
        &["--reference", reference.to_str().unwrap()],
        directory.path(),
    );
}
