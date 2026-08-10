use std::fs::OpenOptions;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/markdup.sam")
}

fn samtools_available() -> bool {
    Command::new("samtools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|status| status.success())
}

fn view(path: &Path, arguments: &[&str]) -> Output {
    let mut command = binary();
    command.arg("view").args(arguments).arg("--no-PG").arg(path);
    command.output().unwrap()
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

    let duplicates = view(&output, &["-f", "1024"]);
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

    let records = view(&output, &[]);
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
    let records = view(&output, &["-f", "1024"]);
    assert!(records.status.success());
    let text = String::from_utf8(records.stdout).unwrap();
    assert_eq!(
        text.lines()
            .map(|line| line.split('\t').next().unwrap())
            .collect::<Vec<_>>(),
        ["low", "low"]
    );
    let all = view(&output, &[]);
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
    let duplicates = view(&output, &["-f", "1024"]);
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

#[test]
fn named_bam_without_eof_is_rejected() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("truncated.bam");
    let converted = binary()
        .args(["view", "-b", "--no-PG"])
        .arg(fixture())
        .arg("-o")
        .arg(&input)
        .output()
        .unwrap();
    assert!(converted.status.success());
    let file = OpenOptions::new().write(true).open(&input).unwrap();
    let length = file.metadata().unwrap().len();
    file.set_len(length - 28).unwrap();
    let quickcheck = binary().arg("quickcheck").arg(&input).output().unwrap();
    assert!(!quickcheck.status.success());
    let output = directory.path().join("output.bam");
    let result = binary()
        .args(["markdup", input.to_str().unwrap(), "-o"])
        .arg(&output)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(!output.exists());
}

#[test]
fn named_cram_without_eof_is_rejected() {
    if !samtools_available() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    std::fs::write(&reference, format!(">chr1\n{}\n", "N".repeat(1000))).unwrap();
    assert!(
        Command::new("samtools")
            .arg("faidx")
            .arg(&reference)
            .status()
            .unwrap()
            .success()
    );
    let input = directory.path().join("truncated.cram");
    let converted = Command::new("samtools")
        .args(["view", "-C", "-T"])
        .arg(&reference)
        .args(["-o", input.to_str().unwrap()])
        .arg(fixture())
        .output()
        .unwrap();
    assert!(converted.status.success());
    let file = OpenOptions::new().write(true).open(&input).unwrap();
    let length = file.metadata().unwrap().len();
    file.set_len(length - 38).unwrap();
    let output = directory.path().join("output.bam");
    let result = binary()
        .args(["markdup", "--reference"])
        .arg(&reference)
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(!output.exists());
}

#[test]
fn standard_input_writes_named_bam() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("output.bam");
    let mut child = binary()
        .args(["markdup", "--no-PG", "-o"])
        .arg(&output)
        .stdin(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(include_bytes!("golden/markdup.sam"))
        .unwrap();
    assert!(child.wait().unwrap().success());
    assert_eq!(duplicate_names(&output), ["low", "low"]);
}

#[test]
fn standard_output_is_complete_bam() {
    let result = binary()
        .args(["markdup", "--no-PG"])
        .arg(fixture())
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("stdout.bam");
    std::fs::write(&output, result.stdout).unwrap();
    assert!(
        binary()
            .arg("quickcheck")
            .arg(&output)
            .status()
            .unwrap()
            .success()
    );
    assert_eq!(duplicate_names(&output), ["low", "low"]);
}

#[test]
fn failed_input_preserves_existing_output() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("bad.sam");
    let output = directory.path().join("output.bam");
    std::fs::write(
        &input,
        include_str!("golden/markdup.sam").replace("SO:coordinate", "SO:queryname"),
    )
    .unwrap();
    std::fs::write(&output, b"existing output").unwrap();
    let result = binary()
        .args(["markdup", input.to_str().unwrap(), "-o"])
        .arg(&output)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(output).unwrap(), b"existing output");
}

#[test]
fn same_input_and_output_is_rejected() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.sam");
    std::fs::copy(fixture(), &input).unwrap();
    let before = std::fs::read(&input).unwrap();
    let result = binary()
        .args(["markdup", input.to_str().unwrap(), "-o"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(input).unwrap(), before);
}

#[test]
fn json_reports_markdup_summary() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("output.bam");
    let result = binary()
        .args(["--json", "markdup", "--no-PG"])
        .arg(fixture())
        .arg("-o")
        .arg(&output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let value: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(value["result"]["command"], "markdup");
    assert_eq!(value["result"]["summary"]["records"], 4);
    assert_eq!(value["result"]["summary"]["written_records"], 4);
}

#[test]
fn json_requires_named_output() {
    let result = binary()
        .args(["--json", "markdup"])
        .arg(fixture())
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(result.stdout.is_empty());
}

#[test]
fn zero_max_read_length_is_rejected() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("output.bam");
    let result = binary()
        .args(["markdup", "-l", "0"])
        .arg(fixture())
        .arg("-o")
        .arg(&output)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(!output.exists());
}

#[test]
fn default_header_records_program() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("output.bam");
    let result = binary()
        .args(["markdup", fixture().to_str().unwrap(), "-o"])
        .arg(&output)
        .output()
        .unwrap();
    assert!(result.status.success());
    let header = view(&output, &["-H"]);
    assert!(header.status.success());
    assert!(
        String::from_utf8(header.stdout)
            .unwrap()
            .contains("@PG\tID:rsomics-bam")
    );
}

#[test]
fn output_failure_is_propagated() {
    let result = rsomics_bam::markdup::write(
        &fixture(),
        rsomics_bam::markdup::Options {
            remove: false,
            clear: false,
            include_fails: false,
            mode: rsomics_bam::markdup::Mode::Template,
            max_read_length: 300,
            additional_threads: Some(0),
            reference: None,
            destination: None,
            program: None,
        },
        FailingWriter,
    );
    assert!(result.is_err());
}

fn duplicate_names(path: &Path) -> Vec<String> {
    let output = view(path, &["-f", "1024"]);
    assert!(output.status.success());
    String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|line| line.split('\t').next().unwrap().to_owned())
        .collect()
}

struct FailingWriter;

impl Write for FailingWriter {
    fn write(&mut self, _buffer: &[u8]) -> io::Result<usize> {
        Err(io::Error::new(io::ErrorKind::BrokenPipe, "closed"))
    }

    fn flush(&mut self) -> io::Result<()> {
        Err(io::Error::new(io::ErrorKind::BrokenPipe, "closed"))
    }
}
