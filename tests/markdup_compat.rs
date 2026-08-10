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
    if !samtools_available() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();
    let input = fixture(name);
    let ours_path = directory.path().join("ours.bam");
    let oracle_path = directory.path().join("oracle.bam");
    let mut ours = binary();
    ours.args(["markdup", "--no-PG"]);
    if remove {
        ours.arg("-r");
    }
    let ours = ours.arg(&input).arg("-o").arg(&ours_path).output().unwrap();
    assert!(
        ours.status.success(),
        "{}",
        String::from_utf8_lossy(&ours.stderr)
    );

    let mut oracle = Command::new("samtools");
    oracle.args(["markdup", "--no-PG"]);
    if remove {
        oracle.arg("-r");
    }
    let oracle = oracle.arg(input).arg(&oracle_path).output().unwrap();
    assert!(
        oracle.status.success(),
        "{}",
        String::from_utf8_lossy(&oracle.stderr)
    );

    let ours = sam_text(&ours_path);
    let oracle = sam_text(&oracle_path);
    assert_eq!(ours, oracle);
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
