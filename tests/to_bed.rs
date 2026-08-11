use std::fs;
use std::io::Write;
use std::process::{Command, Stdio};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn from_stdin(arguments: &[&str], input: &[u8]) -> std::process::Output {
    let mut child = binary()
        .args(arguments)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let write_result = child.stdin.take().unwrap().write_all(input);
    if let Err(error) = write_result {
        assert_eq!(error.kind(), std::io::ErrorKind::BrokenPipe);
    }
    child.wait_with_output().unwrap()
}

#[test]
fn default_bed6_projects_mapped_records_and_mate_suffixes() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("records.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:1000\n\
q1\t65\tchr1\t10\t60\t5M2D3M\t=\t30\t0\tACGTACGT\tIIIIIIII\n\
q2\t145\tchr1\t21\t42\t2S4=1X3M\t=\t10\t0\tAACCGGTTAA\tJJJJJJJJJJ\n\
q3\t4\t*\t0\t0\t*\t*\t0\t0\tAC\tII\n",
    )
    .unwrap();

    let output = binary().arg("to-bed").arg(&input).output().unwrap();

    assert!(
        output.status.success(),
        "stdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b"chr1\t9\t19\tq1/1\t60\t+\nchr1\t20\t28\tq2/2\t42\t-\n"
    );
}

#[test]
fn split_modes_distinguish_reference_skips_and_deletions() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("split.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:1000\n\
split\t0\tchr1\t1\t30\t10M5D5M10N4M2D6M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIII\n",
    )
    .unwrap();

    let split = binary()
        .args(["to-bed", "--split"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        split.status.success(),
        "{}",
        String::from_utf8_lossy(&split.stderr)
    );
    assert_eq!(
        split.stdout,
        b"chr1\t0\t20\tsplit\t30\t+\nchr1\t30\t42\tsplit\t30\t+\n"
    );

    let split_d = binary()
        .args(["to-bed", "--split-d"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        split_d.status.success(),
        "{}",
        String::from_utf8_lossy(&split_d.stderr)
    );
    assert_eq!(
        split_d.stdout,
        b"chr1\t0\t10\tsplit\t30\t+\n\
chr1\t15\t20\tsplit\t30\t+\n\
chr1\t30\t34\tsplit\t30\t+\n\
chr1\t36\t42\tsplit\t30\t+\n"
    );
}

#[test]
fn bed12_renders_blocks_and_color() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("blocks.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:1000\n\
blocks\t16\tchr1\t1\t30\t10M5D5M10N4M2D6M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIII\n",
    )
    .unwrap();

    let output = binary()
        .args(["to-bed", "--bed12", "--color", "1,2,3"])
        .arg(&input)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b"chr1\t0\t42\tblocks\t30\t-\t0\t42\t1,2,3\t2\t20,12\t0,30\n"
    );
}

#[test]
fn integer_scores_and_cigar_preserve_their_text_contracts() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("scores.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:1000\n\
score\t0\tchr1\t5\t60\t2S3M1I2=1X4D\t*\t0\t0\tAACCGGTTT\tIIIIIIIII\tNM:i:3\tXI:i:-7\n",
    )
    .unwrap();

    let edit_distance = binary()
        .args(["to-bed", "--ed"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(edit_distance.status.success());
    assert_eq!(edit_distance.stdout, b"chr1\t4\t14\tscore\t3\t+\n");

    let tag = binary()
        .args(["to-bed", "--tag", "XI"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(tag.status.success());
    assert_eq!(tag.stdout, b"chr1\t4\t14\tscore\t-7\t+\n");

    let cigar = binary()
        .args(["to-bed", "--cigar"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(cigar.status.success());
    assert_eq!(cigar.stdout, b"chr1\t4\t14\tscore\t60\t+\t2S3M1I2=1X4D\n");
}

#[test]
fn bedpe_orders_real_pairs_and_mate1_can_override_coordinate_order() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("pairs.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\tSO:queryname\n\
@SQ\tSN:chr1\tLN:1000\n\
@SQ\tSN:chr2\tLN:1000\n\
pair\t145\tchr1\t10\t30\t5M\tchr2\t50\t0\tAAAAA\tIIIII\tNM:i:2\n\
pair\t65\tchr2\t50\t40\t5M\tchr1\t10\t0\tCCCCC\tIIIII\tNM:i:3\n",
    )
    .unwrap();

    let default = binary()
        .args(["to-bed", "--bedpe"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        default.status.success(),
        "{}",
        String::from_utf8_lossy(&default.stderr)
    );
    assert_eq!(
        default.stdout,
        b"chr1\t9\t14\tchr2\t49\t54\tpair\t30\t-\t+\n"
    );

    let mate1 = binary()
        .args(["to-bed", "--bedpe", "--mate1"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(mate1.status.success());
    assert_eq!(mate1.stdout, b"chr2\t49\t54\tchr1\t9\t14\tpair\t30\t+\t-\n");

    let edit_distance = binary()
        .args(["to-bed", "--bedpe", "--ed"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(edit_distance.status.success());
    assert_eq!(
        edit_distance.stdout,
        b"chr1\t9\t14\tchr2\t49\t54\tpair\t5\t-\t+\n"
    );
}

#[test]
fn named_output_json_and_upstream_input_alias_share_product_contracts() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("records.sam");
    let output = directory.path().join("records.bed");
    fs::write(
        &input,
        b"@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:1000\n\
q1\t0\tchr1\t2\t20\t3M\t*\t0\t0\tAAA\tIII\n\
q2\t4\t*\t0\t0\t*\t*\t0\t0\tCC\tII\n",
    )
    .unwrap();

    let result = binary()
        .args(["--json", "to-bed", "-i"])
        .arg(&input)
        .args(["-o"])
        .arg(&output)
        .output()
        .unwrap();

    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(fs::read(&output).unwrap(), b"chr1\t1\t4\tq1\t20\t+\n");
    let value: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(value["status"], "ok");
    assert_eq!(value["result"]["command"], "to-bed");
    assert_eq!(value["result"]["summary"]["format"], "bed6");
    assert_eq!(value["result"]["summary"]["records_read"], 2);
    assert_eq!(value["result"]["summary"]["records_mapped"], 1);
    assert_eq!(value["result"]["summary"]["records_skipped"], 1);
    assert_eq!(value["result"]["summary"]["rows_written"], 1);
}

#[test]
fn incompatible_output_modes_fail_before_reading_records() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("records.sam");
    fs::write(&input, b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n").unwrap();

    for options in [
        &["--color", "1,2,3"][..],
        &["--mate1"][..],
        &["--bedpe", "--bed12"][..],
        &["--bedpe", "--split"][..],
        &["--bedpe", "--tag", "NM"][..],
        &["--bedpe", "--cigar"][..],
        &["--cigar", "--split"][..],
        &["--cigar", "--bed12"][..],
        &["--ed", "--split"][..],
    ] {
        let result = binary()
            .arg("to-bed")
            .args(options)
            .arg(&input)
            .output()
            .unwrap();
        assert!(
            !result.status.success(),
            "options {options:?} unexpectedly succeeded"
        );
        assert!(result.stdout.is_empty(), "options {options:?}");
    }
}

#[test]
fn failures_preserve_outputs_and_reject_reference_aliases_and_excess_workers() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("records.sam");
    let output = directory.path().join("records.bed");
    let reference = directory.path().join("reference.fa");
    fs::write(
        &input,
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\nq\t0\tchr1\t1\t20\t3M\t*\t0\t0\tAAA\tIII\n",
    )
    .unwrap();
    fs::write(&output, b"existing\n").unwrap();
    fs::write(&reference, b">chr1\nAAAAAAAAAA\n").unwrap();

    let missing_tag = binary()
        .args(["to-bed", "--ed", "-o"])
        .arg(&output)
        .arg(&input)
        .output()
        .unwrap();
    assert!(!missing_tag.status.success());
    assert_eq!(fs::read(&output).unwrap(), b"existing\n");

    let reference_alias = binary()
        .args(["to-bed", "-T"])
        .arg(&reference)
        .args(["-o"])
        .arg(&reference)
        .arg(&input)
        .output()
        .unwrap();
    assert!(!reference_alias.status.success());
    assert_eq!(fs::read(&reference).unwrap(), b">chr1\nAAAAAAAAAA\n");

    let workers = binary()
        .args(["to-bed", "-@", "257"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(!workers.status.success());
    assert!(workers.stdout.is_empty());

    let json = binary()
        .args(["--json", "to-bed"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(!json.status.success());
    assert!(json.stdout.is_empty());
}

#[test]
fn upstream_gzip_sam_with_an_extra_long_header_is_accepted() {
    let input = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/upstream/bedtools-bamtobed/extra-long-header.sam.gz");

    let output = binary()
        .args(["to-bed", "--tag", "NM"])
        .arg(input)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stdout.is_empty());
}

#[test]
fn upstream_2311_regression_corpus_is_exact() {
    let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/upstream/bedtools-bamtobed");
    let cases: &[(&str, &[&str], &str)] = &[
        ("one_block.sam", &[], "chr1\t0\t30\tone_blocks\t40\t-\n"),
        (
            "one_block.sam",
            &["--split"],
            "chr1\t0\t30\tone_blocks\t40\t-\n",
        ),
        ("two_blocks.sam", &[], "chr1\t0\t40\ttwo_blocks\t40\t-\n"),
        (
            "two_blocks.sam",
            &["--split"],
            "chr1\t0\t15\ttwo_blocks\t40\t-\nchr1\t25\t40\ttwo_blocks\t40\t-\n",
        ),
        (
            "three_blocks.sam",
            &[],
            "chr1\t0\t50\tthree_blocks\t40\t-\n",
        ),
        (
            "three_blocks.sam",
            &["--split"],
            "chr1\t0\t10\tthree_blocks\t40\t-\nchr1\t20\t30\tthree_blocks\t40\t-\nchr1\t40\t50\tthree_blocks\t40\t-\n",
        ),
        (
            "three_blocks.sam",
            &["--bed12"],
            "chr1\t0\t50\tthree_blocks\t40\t-\t0\t50\t255,0,0\t3\t10,10,10\t0,20,40\n",
        ),
        (
            "two_blocks_w_D.sam",
            &["--split"],
            "chr1\t0\t15\ttwo_blocks_1_1/2\t40\t+\nchr1\t25\t40\ttwo_blocks_1_1/2\t40\t+\nchr1\t99\t129\ttwo_blocks_1_2/1\t40\t+\nchr1\t0\t15\ttwo_blocks_2_1/2\t40\t+\nchr1\t25\t42\ttwo_blocks_2_1/2\t40\t+\nchr1\t99\t129\ttwo_blocks_2_2/1\t40\t+\n",
        ),
        (
            "two_blocks_w_D.sam",
            &["--split-d"],
            "chr1\t0\t15\ttwo_blocks_1_1/2\t40\t+\nchr1\t25\t40\ttwo_blocks_1_1/2\t40\t+\nchr1\t99\t129\ttwo_blocks_1_2/1\t40\t+\nchr1\t0\t15\ttwo_blocks_2_1/2\t40\t+\nchr1\t25\t35\ttwo_blocks_2_1/2\t40\t+\nchr1\t37\t42\ttwo_blocks_2_1/2\t40\t+\nchr1\t99\t129\ttwo_blocks_2_2/1\t40\t+\n",
        ),
        (
            "two_blocks_w_D.sam",
            &["--bed12"],
            "chr1\t0\t40\ttwo_blocks_1_1/2\t40\t+\t0\t40\t255,0,0\t2\t15,15\t0,25\nchr1\t99\t129\ttwo_blocks_1_2/1\t40\t+\t99\t129\t255,0,0\t1\t30\t0\nchr1\t0\t42\ttwo_blocks_2_1/2\t40\t+\t0\t42\t255,0,0\t2\t15,17\t0,25\nchr1\t99\t129\ttwo_blocks_2_2/1\t40\t+\t99\t129\t255,0,0\t1\t30\t0\n",
        ),
        (
            "two_blocks_w_D.sam",
            &["--bed12", "--split-d"],
            "chr1\t0\t40\ttwo_blocks_1_1/2\t40\t+\t0\t40\t255,0,0\t2\t15,15\t0,25\nchr1\t99\t129\ttwo_blocks_1_2/1\t40\t+\t99\t129\t255,0,0\t1\t30\t0\nchr1\t0\t42\ttwo_blocks_2_1/2\t40\t+\t0\t42\t255,0,0\t3\t15,10,5\t0,25,37\nchr1\t99\t129\ttwo_blocks_2_2/1\t40\t+\t99\t129\t255,0,0\t1\t30\t0\n",
        ),
        (
            "numeric_tag.sam",
            &["--tag", "NM"],
            "1\t9998\t10056\tHISEQ1:18:H8VC6ADXX:1:1201:3360:80789/2\t2\t+\n1\t9998\t10064\tHISEQ1:18:H8VC6ADXX:1:1208:16920:47717/2\t2\t+\n",
        ),
    ];

    for &(name, options, expected) in cases {
        let output = binary()
            .arg("to-bed")
            .args(options)
            .arg(root.join(name))
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{name} {options:?}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(output.stdout, expected.as_bytes(), "{name} {options:?}");
    }
}

#[test]
fn ordinary_gzip_sam_is_detected_on_standard_input() {
    let input = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/upstream/bedtools-bamtobed/extra-long-header.sam.gz");
    let mut child = binary()
        .args(["to-bed", "--tag", "NM", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let write_result = child
        .stdin
        .take()
        .unwrap()
        .write_all(&fs::read(input).unwrap());
    if let Err(error) = write_result {
        assert_eq!(error.kind(), std::io::ErrorKind::BrokenPipe);
    }
    let output = child.wait_with_output().unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stdout.is_empty());
}

#[test]
fn sam_bam_cram_and_standard_input_have_one_projection_contract() {
    let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let sam = root.join("tests/golden/stats/1_map_cigar.sam");
    let cram = root.join("tests/golden/stats/1_map_cigar.cram");
    let reference = root.join("tests/golden/stats/test.fa");
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("records.bam");
    let encoded = binary()
        .args(["view", "-b", "--no-pg", "-o"])
        .arg(&bam)
        .arg(&sam)
        .output()
        .unwrap();
    assert!(
        encoded.status.success(),
        "{}",
        String::from_utf8_lossy(&encoded.stderr)
    );
    let expected = b"alpha\t0\t35\tr1/1\t40\t+\nalpha\t65\t100\tr1/2\t40\t-\n";

    for input in [&sam, &bam, &cram] {
        let file = binary()
            .args(["to-bed", "-T"])
            .arg(&reference)
            .arg(input)
            .output()
            .unwrap();
        assert!(
            file.status.success(),
            "{}: {}",
            input.display(),
            String::from_utf8_lossy(&file.stderr)
        );
        assert_eq!(file.stdout, expected, "{}", input.display());

        let reference = reference.to_str().unwrap();
        let stdin = from_stdin(&["to-bed", "-T", reference, "-"], &fs::read(input).unwrap());
        assert!(
            stdin.status.success(),
            "{} stdin: {}",
            input.display(),
            String::from_utf8_lossy(&stdin.stderr)
        );
        assert_eq!(stdin.stdout, expected, "{} stdin", input.display());
    }
}

#[test]
fn split_uses_the_decoded_long_cigar_in_bam() {
    let directory = tempfile::tempdir().unwrap();
    let sam = directory.path().join("long.sam");
    let bam = directory.path().join("long.bam");
    let cigar = "1M1I".repeat(32_768);
    let sequence = "A".repeat(65_536);
    fs::write(
        &sam,
        format!(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40000\nlong\t0\tchr1\t1\t60\t{cigar}\t*\t0\t0\t{sequence}\t*\n"
        ),
    )
    .unwrap();
    let encoded = binary()
        .args(["view", "-b", "--no-pg", "-o"])
        .arg(&bam)
        .arg(&sam)
        .output()
        .unwrap();
    assert!(
        encoded.status.success(),
        "{}",
        String::from_utf8_lossy(&encoded.stderr)
    );

    let output = binary()
        .args(["to-bed", "--split"])
        .arg(&bam)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b"chr1\t0\t32768\tlong\t60\t+\n");
}

#[test]
fn bedpe_renders_unmapped_ends_and_sums_only_present_edit_distance() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("unmapped.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\tSO:queryname\n\
@SQ\tSN:chr1\tLN:1000\n\
pair\t73\tchr1\t10\t55\t5M\t*\t0\t0\tAAAAA\tIIIII\tNM:i:3\n\
pair\t133\t*\t0\t0\t*\tchr1\t10\t0\tCCCCC\tIIIII\n",
    )
    .unwrap();

    let default = binary()
        .args(["to-bed", "--bedpe"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        default.status.success(),
        "{}",
        String::from_utf8_lossy(&default.stderr)
    );
    assert_eq!(default.stdout, b".\t-1\t-1\tchr1\t9\t14\tpair\t0\t.\t+\n");

    let mate1 = binary()
        .args(["to-bed", "--bedpe", "--mate1", "--ed"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        mate1.status.success(),
        "{}",
        String::from_utf8_lossy(&mate1.stderr)
    );
    assert_eq!(mate1.stdout, b"chr1\t9\t14\t.\t-1\t-1\tpair\t3\t+\t.\n");
}

#[test]
fn bedpe_pair_failures_roll_back_named_output() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("misgrouped.sam");
    let output = directory.path().join("pairs.bedpe");
    fs::write(
        &input,
        b"@HD\tVN:1.6\tSO:queryname\n\
@SQ\tSN:chr1\tLN:1000\n\
good\t65\tchr1\t1\t30\t2M\t=\t10\t0\tAA\tII\n\
good\t129\tchr1\t10\t30\t2M\t=\t1\t0\tCC\tII\n\
orphan\t65\tchr1\t20\t30\t2M\t=\t30\t0\tGG\tII\n\
other\t129\tchr1\t30\t30\t2M\t=\t20\t0\tTT\tII\n",
    )
    .unwrap();
    fs::write(&output, b"existing\n").unwrap();

    let result = binary()
        .args(["to-bed", "--bedpe", "-o"])
        .arg(&output)
        .arg(&input)
        .output()
        .unwrap();

    assert!(!result.status.success());
    assert!(result.stdout.is_empty());
    assert_eq!(fs::read(&output).unwrap(), b"existing\n");
}

#[test]
fn numeric_tag_stays_in_the_score_column_for_every_split_block() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("tagged.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:1000\n\
signed\t0\tchr1\t5\t60\t3M2N4M\t*\t0\t0\tAAAAAAA\tIIIIIII\tXI:i:-7\n",
    )
    .unwrap();

    let result = binary()
        .args(["to-bed", "--tag", "XI", "--split"])
        .arg(&input)
        .output()
        .unwrap();

    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        result.stdout,
        b"chr1\t4\t7\tsigned\t-7\t+\nchr1\t9\t13\tsigned\t-7\t+\n"
    );
    assert!(
        result
            .stdout
            .split(|byte| *byte == b'\n')
            .filter(|line| !line.is_empty())
            .all(|line| line.split(|byte| *byte == b'\t').count() == 6)
    );
}

#[test]
fn malformed_score_requests_and_cli_values_fail_without_output() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("tags.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:1000\n\
tagged\t0\tchr1\t1\t60\t3M\t*\t0\t0\tAAA\tIII\tXF:f:1.5\n",
    )
    .unwrap();

    for options in [
        &["--tag", "NM"][..],
        &["--tag", "XF"][..],
        &["--tag", "A"][..],
        &["--tag", "1A"][..],
        &["--bed12", "--color", "256,0,0"][..],
        &["--bed12", "--color", "1,2"][..],
    ] {
        let result = binary()
            .arg("to-bed")
            .args(options)
            .arg(&input)
            .output()
            .unwrap();
        assert!(!result.status.success(), "options {options:?}");
        assert!(result.stdout.is_empty(), "options {options:?}");
    }
}

#[test]
fn integer_score_boundaries_preserve_signed_and_unsigned_values() {
    let directory = tempfile::tempdir().unwrap();
    let sam = directory.path().join("integer-types.sam");
    let bam = directory.path().join("integer-types.bam");
    fs::write(
        &sam,
        b"@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:1000\n\
limits\t0\tchr1\t1\t60\t3M\t*\t0\t0\tAAA\tIII\tTC:i:-128\tTU:i:255\tTS:i:-32768\tTV:i:65535\tTI:i:-2147483648\tTJ:i:4294967295\n",
    )
    .unwrap();
    let encoded = binary()
        .args(["view", "-b", "--no-pg", "-o"])
        .arg(&bam)
        .arg(&sam)
        .output()
        .unwrap();
    assert!(
        encoded.status.success(),
        "{}",
        String::from_utf8_lossy(&encoded.stderr)
    );

    for (tag, score) in [
        ("TC", "-128"),
        ("TU", "255"),
        ("TS", "-32768"),
        ("TV", "65535"),
        ("TI", "-2147483648"),
        ("TJ", "4294967295"),
    ] {
        let output = binary()
            .args(["to-bed", "--tag", tag])
            .arg(&bam)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{tag}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(
            String::from_utf8(output.stdout).unwrap(),
            format!("chr1\t0\t3\tlimits\t{score}\t+\n")
        );
    }
}

#[test]
fn split_d_retains_bedtools_zero_length_boundary_blocks() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("adjacent-boundaries.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:1000\n\
boundary\t0\tchr1\t5\t60\t3M4D5N4M\t*\t0\t0\tAAAAAAA\tIIIIIII\n",
    )
    .unwrap();

    let output = binary()
        .args(["to-bed", "--split-d"])
        .arg(&input)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b"chr1\t4\t7\tboundary\t60\t+\n\
chr1\t11\t11\tboundary\t60\t+\n\
chr1\t16\t20\tboundary\t60\t+\n"
    );
}

#[test]
fn closed_standard_output_is_a_command_failure() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("many.sam");
    let mut input = std::io::BufWriter::new(fs::File::create(&path).unwrap());
    writeln!(input, "@HD\tVN:1.6").unwrap();
    writeln!(input, "@SQ\tSN:chr1\tLN:1").unwrap();
    for index in 0..100_000 {
        writeln!(input, "q{index}\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI").unwrap();
    }
    drop(input);

    let mut child = binary()
        .arg("to-bed")
        .arg(path)
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    drop(child.stdout.take());
    assert!(!child.wait().unwrap().success());
}
