use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/markdup.sam")
}

fn run_text(text: &str) -> (tempfile::TempDir, PathBuf, Output) {
    run_text_args(text, &[])
}

fn run_text_args(text: &str, arguments: &[&str]) -> (tempfile::TempDir, PathBuf, Output) {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.sam");
    let output = directory.path().join("output.bam");
    std::fs::write(&input, text).unwrap();
    let mut command = binary();
    command.arg("markdup").args(arguments).arg(&input).arg("-o");
    let result = command.arg(&output).output().unwrap();
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

#[test]
fn clear_removes_prior_duplicate_state_before_marking() {
    let input = include_str!("golden/markdup.sam")
        .replacen("high\t99", "high\t1123", 1)
        .replacen("high\t147", "high\t1171", 1)
        .replace(
            "\tMC:Z:10M\tms:i:400",
            "\tMC:Z:10M\tms:i:400\tdo:Z:old\tdt:Z:LB",
        );
    let (_directory, output, result) = run_text_args(&input, &["-c"]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let records = Command::new("samtools")
        .args(["view", "-f", "1024"])
        .arg(&output)
        .output()
        .unwrap();
    assert!(records.status.success());
    let text = String::from_utf8(records.stdout).unwrap();
    assert_eq!(
        text.lines()
            .map(|line| line.split('\t').next().unwrap())
            .collect::<Vec<_>>(),
        ["low", "low"]
    );
    let all = Command::new("samtools")
        .arg("view")
        .arg(output)
        .output()
        .unwrap();
    let all = String::from_utf8(all.stdout).unwrap();
    assert!(!all.contains("\tdo:Z:"));
    assert!(!all.contains("\tdt:Z:"));
}

#[test]
fn max_read_length_retains_clipped_duplicate_candidate() {
    let sequence_a = "A".repeat(100);
    let quality_a = "I".repeat(100);
    let sequence_b = "A".repeat(401);
    let quality_b = "!".repeat(401);
    let input = format!(
        "@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:1000\n\
original\t0\tchr1\t100\t60\t100M\t*\t0\t0\t{sequence_a}\t{quality_a}\n\
duplicate\t0\tchr1\t401\t60\t301S100M\t*\t0\t0\t{sequence_b}\t{quality_b}\n"
    );
    let (_directory, output, result) = run_text_args(&input, &["-l", "301"]);
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
    let text = String::from_utf8(duplicates.stdout).unwrap();
    assert_eq!(text.split('\t').next(), Some("duplicate"));
    assert_eq!(text.lines().count(), 1);
}

#[test]
fn paired_score_tie_uses_lexical_name() {
    let input = "@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:1000\n\
zeta\t99\tchr1\t100\t60\t10M\t=\t200\t110\tAAAAAAAAAA\tIIIIIIIIII\tMC:Z:10M\tms:i:400\n\
alpha\t99\tchr1\t100\t60\t10M\t=\t200\t110\tAAAAAAAAAA\tIIIIIIIIII\tMC:Z:10M\tms:i:400\n\
zeta\t147\tchr1\t200\t60\t10M\t=\t100\t-110\tTTTTTTTTTT\tIIIIIIIIII\tMC:Z:10M\tms:i:400\n\
alpha\t147\tchr1\t200\t60\t10M\t=\t100\t-110\tTTTTTTTTTT\tIIIIIIIIII\tMC:Z:10M\tms:i:400\n";
    let (_directory, output, result) = run_text(input);
    assert!(result.status.success());
    assert_eq!(duplicate_names(&output), ["zeta", "zeta"]);
}

#[test]
fn paired_record_outranks_single_record() {
    let input = "@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:1000\n\
single\t0\tchr1\t100\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n\
pair\t99\tchr1\t100\t60\t10M\t=\t200\t110\tAAAAAAAAAA\t!!!!!!!!!!\tMC:Z:10M\tms:i:0\n\
pair\t147\tchr1\t200\t60\t10M\t=\t100\t-110\tTTTTTTTTTT\t!!!!!!!!!!\tMC:Z:10M\tms:i:0\n";
    let (_directory, output, result) = run_text(input);
    assert!(result.status.success());
    assert_eq!(duplicate_names(&output), ["single"]);
}

#[test]
fn prior_duplicate_flags_remain_without_clear() {
    let input = include_str!("golden/markdup.sam")
        .replacen("high\t99", "high\t1123", 1)
        .replacen("high\t147", "high\t1171", 1);
    let (_directory, output, result) = run_text(&input);
    assert!(result.status.success());
    assert_eq!(duplicate_names(&output).len(), 4);
}

#[test]
fn wrong_ms_type_is_rejected() {
    let input = include_str!("golden/markdup.sam")
        .replace("ms:i:400", "ms:Z:400")
        .replace("ms:i:0", "ms:Z:0");
    let (_directory, output, result) = run_text(&input);
    assert!(!result.status.success());
    assert!(!output.exists());
}

#[test]
fn hard_clipping_contributes_to_unclipped_start() {
    let sequence = "A".repeat(100);
    let high_quality = "I".repeat(100);
    let low_quality = "!".repeat(100);
    let input = format!(
        "@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:1000\n\
original\t0\tchr1\t100\t60\t100M\t*\t0\t0\t{sequence}\t{high_quality}\n\
duplicate\t0\tchr1\t401\t60\t301H100M\t*\t0\t0\t{sequence}\t{low_quality}\n"
    );
    let (_directory, output, result) = run_text_args(&input, &["-l", "301"]);
    assert!(result.status.success());
    assert_eq!(duplicate_names(&output), ["duplicate"]);
}

#[test]
fn long_cigar_uses_decoded_operations() {
    let long_cigar = "1M1I".repeat(32_768);
    let original_sequence = "A".repeat(32_768);
    let duplicate_sequence = "A".repeat(65_536);
    let original_quality = "I".repeat(32_768);
    let duplicate_quality = "!".repeat(65_536);
    let input = format!(
        "@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:100000\n\
original\t0\tchr1\t100\t60\t32768M\t*\t0\t0\t{original_sequence}\t{original_quality}\n\
duplicate\t0\tchr1\t100\t60\t{long_cigar}\t*\t0\t0\t{duplicate_sequence}\t{duplicate_quality}\n"
    );
    let (_directory, output, result) = run_text(&input);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(duplicate_names(&output), ["duplicate"]);
}

fn duplicate_names(path: &Path) -> Vec<String> {
    let output = Command::new("samtools")
        .args(["view", "-f", "1024"])
        .arg(path)
        .output()
        .unwrap();
    assert!(output.status.success());
    String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|line| line.split('\t').next().unwrap().to_owned())
        .collect()
}
