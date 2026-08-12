use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn write_alignment(root: &Path, records: &str) -> PathBuf {
    let path = root.join("input.sam");
    fs::write(
        &path,
        format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:14\n{records}"),
    )
    .unwrap();
    path
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

fn complete_records() -> &'static str {
    "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tMD:Z:4\n\
r2\t0\tchr1\t5\t60\t4M\t*\t0\t0\tTCGT\tIIII\tMD:Z:1A2\n\
r3\t0\tchr1\t9\t60\t2M2D2M\t*\t0\t0\tAATT\tIIII\tMD:Z:2^CG2\n"
}

#[test]
fn reconstructs_matches_substitutions_and_deletions() {
    let directory = tempfile::tempdir().unwrap();
    let input = write_alignment(directory.path(), complete_records());
    let output = run({
        let mut command = binary();
        command.args(["reference", "--quiet"]).arg(input);
        command
    });
    assert_eq!(output.stdout, b">chr1\nACGTTAGTAACGTT\n");
    assert!(output.stderr.is_empty());
}

#[test]
fn standard_input_matches_file_input() {
    let directory = tempfile::tempdir().unwrap();
    let input = write_alignment(directory.path(), complete_records());
    let expected = run({
        let mut command = binary();
        command.args(["reference", "--quiet"]).arg(&input);
        command
    });
    let mut child = binary()
        .args(["reference", "--quiet", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&fs::read(input).unwrap())
        .unwrap();
    let actual = child.wait_with_output().unwrap();
    assert!(actual.status.success());
    assert_eq!(actual.stdout, expected.stdout);
}

#[test]
fn named_output_and_json_use_product_contract() {
    let directory = tempfile::tempdir().unwrap();
    let input = write_alignment(directory.path(), complete_records());
    let output_path = directory.path().join("reference.fa");
    let output = run({
        let mut command = binary();
        command
            .args(["--json", "reference", "--quiet", "--output"])
            .arg(&output_path)
            .arg(input);
        command
    });
    assert_eq!(fs::read(output_path).unwrap(), b">chr1\nACGTTAGTAACGTT\n");
    let json: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(json["result"]["command"], "reference");
    assert_eq!(json["result"]["summary"]["references"], 1);
    assert_eq!(json["result"]["summary"]["bases"], 14);
    assert_eq!(json["result"]["summary"]["known_bases"], 14);
}

#[test]
fn invalid_evidence_fails_without_replacing_output() {
    let directory = tempfile::tempdir().unwrap();
    let output_path = directory.path().join("reference.fa");
    fs::write(&output_path, b"sentinel").unwrap();
    let input = write_alignment(
        directory.path(),
        "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tMD:Z:4\n\
r2\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTCGT\tIIII\tMD:Z:4\n",
    );
    let output = binary()
        .args(["reference", "--quiet", "--output"])
        .arg(&output_path)
        .arg(input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_eq!(fs::read(output_path).unwrap(), b"sentinel");
}

#[test]
fn json_and_regions_require_named_or_indexed_input() {
    let directory = tempfile::tempdir().unwrap();
    let input = write_alignment(directory.path(), complete_records());
    let json = binary()
        .args(["--json", "reference"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(!json.status.success());

    let region = binary()
        .args(["reference", "--region", "chr1:2-4"])
        .arg(input)
        .output()
        .unwrap();
    assert!(!region.status.success());
}

#[test]
fn extracts_embedded_cram_reference_by_indexed_region() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.cram");
    fs::copy(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/checksum/chk2.cram"),
        &input,
    )
    .unwrap();
    run({
        let mut command = binary();
        command.arg("index").arg(&input);
        command
    });
    let output = run({
        let mut command = binary();
        command
            .args(["reference", "--embedded", "--quiet", "--region", "17:1-500"])
            .arg(input);
        command
    });
    assert_eq!(
        output.stdout,
        b">17:1-500\n\
AAGCTTCTCACCCTGTTCCTGCATAGATAATTGCATGACAATTGCCTTGTCCCTGCTGAA\n\
TGTGCTCTGGGGTCTCTGGGGTCTCACCCACGACCAACTCCCTGGGCCTGGCACCAGGGA\n\
GCTTAACAAACATCTGTCCAGCGAATACCTGCATCCCTAGAAGTGAAGCCACCGCCCAAA\n\
GACACGCCCATGTCCAGCTTAACCTGCATCCCTAGAAGTGAAGGCACCGCCCAAAGACAC\n\
GCCCATGTCCAGCTTATTCTGCCCAGTTCCTCTCCAGAAAGGCTGCATGGTTGACACACA\n\
GTGCCTGCGACAAAGCTGAATGCTATCATTTAAAAACTCCTTGCTGGTTTGAGAGGCAGA\n\
AAATGATATCTCATAGTTGCTTTACTTTGCATATTTTAAAATTGTGACTTTCATGGCATA\n\
AATAATACTGGTTTATTACAGAAGCACTAGAAAATGCATGTGGACAAAAGTTGGGATTAG\n\
GAGAGAGAAATGAAGACATA\n"
    );
}

#[test]
fn embedded_mode_rejects_non_cram_and_missing_blocks() {
    let directory = tempfile::tempdir().unwrap();
    let sam = write_alignment(directory.path(), complete_records());
    let non_cram = binary()
        .args(["reference", "--embedded"])
        .arg(sam)
        .output()
        .unwrap();
    assert!(!non_cram.status.success());

    let cram =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/cram-size/version-3.1.cram");
    let missing = binary()
        .args(["reference", "--embedded"])
        .arg(cram)
        .output()
        .unwrap();
    assert!(!missing.status.success());
    assert!(String::from_utf8_lossy(&missing.stderr).contains("no embedded reference"));
}

#[test]
fn truncated_embedded_cram_preserves_named_output() {
    let directory = tempfile::tempdir().unwrap();
    let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/checksum/chk2.cram");
    let mut bytes = fs::read(source).unwrap();
    bytes.truncate(bytes.len() - 100);
    let input = directory.path().join("truncated.cram");
    fs::write(&input, bytes).unwrap();
    let output_path = directory.path().join("reference.fa");
    fs::write(&output_path, b"sentinel").unwrap();
    let output = binary()
        .args(["reference", "--embedded", "--quiet", "--output"])
        .arg(&output_path)
        .arg(input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_eq!(fs::read(output_path).unwrap(), b"sentinel");
}

#[test]
fn unsorted_records_fail() {
    let directory = tempfile::tempdir().unwrap();
    let input = write_alignment(
        directory.path(),
        "r2\t0\tchr1\t5\t60\t4M\t*\t0\t0\tTCGT\tIIII\tMD:Z:1A2\n\
r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tMD:Z:4\n",
    );
    let output = binary()
        .args(["reference", "--quiet"])
        .arg(input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("coordinate-sorted"));
}

#[test]
fn handles_clipping_insertions_and_reference_skips() {
    let directory = tempfile::tempdir().unwrap();
    let input = write_alignment(
        directory.path(),
        "r1\t0\tchr1\t1\t60\t2S2M1I2M2N2M1D1M1S\t*\t0\t0\tNNACGTTAACN\tIIIIIIIIIII\tMD:Z:6^G1\n",
    );
    let output = run({
        let mut command = binary();
        command.args(["reference", "--quiet"]).arg(input);
        command
    });
    assert_eq!(output.stdout, b">chr1\nACTTNNAAGCNNNN\n");
}

#[test]
fn rejects_md_operations_that_disagree_with_explicit_cigar_operators() {
    let directory = tempfile::tempdir().unwrap();
    for (cigar, md) in [("4X", "4"), ("4=", "A3")] {
        let input = write_alignment(
            directory.path(),
            &format!("r1\t0\tchr1\t1\t60\t{cigar}\t*\t0\t0\tACGT\tIIII\tMD:Z:{md}\n"),
        );
        let output = binary()
            .args(["reference", "--quiet"])
            .arg(input)
            .output()
            .unwrap();
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains("incompatible"));
    }
}

#[test]
fn bam_streaming_and_indexed_region_paths_match() {
    let directory = tempfile::tempdir().unwrap();
    let sam = write_alignment(directory.path(), complete_records());
    let bam = directory.path().join("input.bam");
    run({
        let mut command = binary();
        command
            .args(["view", "--no-pg", "-b", "-o"])
            .arg(&bam)
            .arg(sam);
        command
    });
    let full = run({
        let mut command = binary();
        command
            .args(["reference", "--quiet", "--threads", "2"])
            .arg(&bam);
        command
    });
    assert_eq!(full.stdout, b">chr1\nACGTTAGTAACGTT\n");

    run({
        let mut command = binary();
        command.arg("index").arg(&bam);
        command
    });
    let region = run({
        let mut command = binary();
        command
            .args(["reference", "--quiet", "--region", "chr1:2-10"])
            .arg(bam);
        command
    });
    assert_eq!(region.stdout, b">chr1:2-10\nCGTTAGTAA\n");
}

#[test]
#[ignore = "requires the pinned samtools 1.24 compatibility oracle"]
fn md_and_embedded_regions_match_samtools_1_24() {
    let version = Command::new("samtools").arg("--version").output().unwrap();
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("samtools 1.24\n"));
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.cram");
    fs::copy(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/checksum/chk2.cram"),
        &input,
    )
    .unwrap();
    run({
        let mut command = binary();
        command.arg("index").arg(&input);
        command
    });

    for extra in [None, Some("--embedded")] {
        let mut oracle = Command::new("samtools");
        oracle.args(["reference", "--quiet", "--region", "17:1-500"]);
        if let Some(extra) = extra {
            oracle.arg(extra);
        }
        let expected = run({
            oracle.arg(&input);
            oracle
        });

        let actual = run({
            let mut command = binary();
            command.args(["reference", "--quiet", "--region", "17:1-500"]);
            if let Some(extra) = extra {
                command.arg(extra);
            }
            command.arg(&input);
            command
        });
        assert_eq!(actual.stdout, expected.stdout);
    }
}
