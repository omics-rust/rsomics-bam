use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/markdup.sam")
}

fn run_text(text: &str) -> (tempfile::TempDir, PathBuf, Output) {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.sam");
    let output = directory.path().join("output.bam");
    std::fs::write(&input, text).unwrap();
    let result = binary()
        .args(["markdup", input.to_str().unwrap(), "-o"])
        .arg(&output)
        .output()
        .unwrap();
    (directory, output, result)
}

#[test]
fn lower_scoring_pair_is_marked_duplicate() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("marked.bam");
    let result = binary()
        .args(["markdup", fixture().to_str().unwrap(), "-o"])
        .arg(&output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );

    let duplicates = Command::new("samtools")
        .args(["view", "-f", "1024"])
        .arg(output)
        .output()
        .unwrap();
    assert!(duplicates.status.success());
    let names = String::from_utf8(duplicates.stdout)
        .unwrap()
        .lines()
        .map(|line| line.split('\t').next().unwrap().to_owned())
        .collect::<Vec<_>>();
    assert_eq!(names, ["low", "low"]);
}

#[test]
fn removal_omits_duplicate_records() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("deduplicated.bam");
    let result = binary()
        .args(["markdup", "-r", fixture().to_str().unwrap(), "-o"])
        .arg(&output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );

    let records = Command::new("samtools")
        .args(["view"])
        .arg(output)
        .output()
        .unwrap();
    assert!(records.status.success());
    let names = String::from_utf8(records.stdout)
        .unwrap()
        .lines()
        .map(|line| line.split('\t').next().unwrap().to_owned())
        .collect::<Vec<_>>();
    assert_eq!(names, ["high", "high"]);
}

#[test]
fn query_name_sort_declaration_is_rejected() {
    let input = include_str!("golden/markdup.sam").replace("SO:coordinate", "SO:queryname");
    let (_directory, output, result) = run_text(&input);
    assert!(!result.status.success());
    assert!(!output.exists());
}

#[test]
fn observed_coordinate_decrease_is_rejected() {
    let input = "@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:1000\n\
later\t0\tchr1\t200\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n\
earlier\t0\tchr1\t100\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n";
    let (_directory, output, result) = run_text(input);
    assert!(!result.status.success());
    assert!(!output.exists());
}

#[test]
fn paired_record_without_mc_is_rejected() {
    let input = include_str!("golden/markdup.sam").replace("\tMC:Z:10M", "");
    let (_directory, output, result) = run_text(&input);
    assert!(!result.status.success());
    assert!(!output.exists());
}

#[test]
fn paired_collision_without_ms_is_rejected() {
    let input = include_str!("golden/markdup.sam")
        .replace("\tms:i:400", "")
        .replace("\tms:i:0", "");
    let (_directory, output, result) = run_text(&input);
    assert!(!result.status.success());
    assert!(!output.exists());
}
