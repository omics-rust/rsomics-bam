use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn fixture(name: impl AsRef<Path>) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden/stats")
        .join(name)
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

fn stable_body(output: &[u8]) -> &[u8] {
    let mut start = 0;
    for _ in 0..3 {
        start += output[start..]
            .iter()
            .position(|&byte| byte == b'\n')
            .unwrap()
            + 1;
    }
    &output[start..]
}

fn assert_bytes(actual: &[u8], expected: &[u8], label: &str) {
    if actual == expected {
        return;
    }
    let offset = actual
        .iter()
        .zip(expected)
        .position(|(left, right)| left != right)
        .unwrap_or(actual.len().min(expected.len()));
    let start = offset.saturating_sub(80);
    let actual_end = (offset + 160).min(actual.len());
    let expected_end = (offset + 160).min(expected.len());
    panic!(
        "{label} differs at byte {offset}\nactual: {:?}\nexpected: {:?}",
        String::from_utf8_lossy(&actual[start..actual_end]),
        String::from_utf8_lossy(&expected[start..expected_end])
    );
}

fn assert_report(output: &Output, expected: &str, label: &str) {
    assert_bytes(
        stable_body(&output.stdout),
        &fs::read(fixture(expected)).unwrap(),
        label,
    );
}

fn rfs_lines(output: &[u8]) -> Vec<u8> {
    output
        .split_inclusive(|&byte| byte == b'\n')
        .filter(|line| line.starts_with(b"RFS\t"))
        .flatten()
        .copied()
        .collect()
}

fn line_with_prefix(output: &[u8], prefix: &[u8]) -> Vec<u8> {
    output
        .split_inclusive(|&byte| byte == b'\n')
        .find(|line| line.starts_with(prefix))
        .unwrap()
        .to_vec()
}

#[test]
fn default_reference_report_matches_samtools_1_24() {
    let output = run({
        let mut command = binary();
        command
            .args(["stats", "-r"])
            .arg(fixture("test.fa"))
            .arg(fixture("1_map_cigar.sam"));
        command
    });
    assert_report(&output, "1.stats.expected", "default reference report");
}

#[test]
fn multithreaded_cram_matches_samtools_1_24() {
    let output = run({
        let mut command = binary();
        command
            .args(["stats", "-@", "2", "-r"])
            .arg(fixture("test.fa"))
            .arg(fixture("1_map_cigar.cram"));
        command
    });
    assert_report(&output, "1.stats.expected", "multithreaded CRAM");
}

#[test]
fn sam_stdin_matches_file_input() {
    let output = run({
        let mut command = binary();
        command
            .args(["stats", "-r"])
            .arg(fixture("test.fa"))
            .stdin(Stdio::from(
                fs::File::open(fixture("1_map_cigar.sam")).unwrap(),
            ));
        command
    });
    assert_report(&output, "1.stats.expected", "SAM stdin");
}

#[test]
fn cigar_and_alignment_classification_matrix_matches_samtools_1_24() {
    for (input, expected, reference, extra) in [
        ("1_map_cigar.sam", "1.stats.expected", true, &[][..]),
        (
            "1_map_cigar_large.sam",
            "1.stats.large.expected",
            false,
            &[][..],
        ),
        (
            "2_equal_cigar_full_seq.sam",
            "2.stats.expected",
            true,
            &[][..],
        ),
        (
            "2_equal_cigar_full_seq_large.sam",
            "2.stats.large.expected",
            false,
            &[][..],
        ),
        (
            "3_map_cigar_equal_seq.sam",
            "3.stats.expected",
            true,
            &[][..],
        ),
        (
            "3_map_cigar_equal_seq_large.sam",
            "3.stats.large.expected",
            false,
            &[][..],
        ),
        ("4_X_cigar_full_seq.sam", "4.stats.expected", true, &[][..]),
        (
            "4_X_cigar_full_seq_large.sam",
            "4.stats.large.expected",
            false,
            &[][..],
        ),
        ("5_insert_cigar.sam", "5.stats.expected", true, &[][..]),
        (
            "5_insert_cigar_large.sam",
            "5.stats.large.expected",
            false,
            &[][..],
        ),
        ("5_insert_cigar.sam", "6.stats.expected", true, &["-i", "0"]),
        ("7_supp.sam", "7.stats.expected", true, &[][..]),
        ("7_supp_large.sam", "7.stats.large.expected", false, &[][..]),
        ("8_secondary.sam", "8.stats.expected", true, &[][..]),
        (
            "8_secondary_large.sam",
            "8.stats.large.expected",
            false,
            &[][..],
        ),
    ] {
        let output = run({
            let mut command = binary();
            command.arg("stats").args(extra);
            if reference {
                command.args(["-r"]).arg(fixture("test.fa"));
            }
            command.arg(fixture(input));
            command
        });
        assert_report(&output, expected, input);
    }
}

#[test]
fn unsorted_input_stops_streaming_coverage_like_samtools_1_24() {
    let output = run({
        let mut command = binary();
        command.arg("stats").arg(fixture("unsorted.sam"));
        command
    });
    assert_eq!(
        line_with_prefix(&output.stdout, b"SN\tis sorted:"),
        b"SN\tis sorted:\t0\t# not sorted by coordinate\n"
    );
    assert!(
        !output
            .stdout
            .split(|&byte| byte == b'\n')
            .any(|line| line.starts_with(b"COV\t"))
    );
}

#[test]
fn indexed_multi_region_query_counts_each_physical_record_once() {
    let output = run({
        let mut command = binary();
        command
            .arg("stats")
            .arg(fixture("11_target.bam"))
            .args(["ref1:1-3", "ref1:5-7"]);
        command
    });
    assert_eq!(
        line_with_prefix(&output.stdout, b"CHK\t"),
        b"CHK\tc03ce880\t86839efe\t807ea85e\n"
    );
    assert_eq!(
        line_with_prefix(&output.stdout, b"SN\traw total sequences:"),
        b"SN\traw total sequences:\t4\t# excluding supplementary and secondary reads\n"
    );
}

#[test]
fn targets_regions_overlap_and_id_selection_match_samtools_1_24() {
    for (label, expected, arguments) in [
        (
            "targets",
            "11.stats.expected",
            vec![
                "-t".into(),
                fixture("11.stats.targets").into_os_string(),
                fixture("11_target.sam").into_os_string(),
            ],
        ),
        (
            "targets threshold",
            "11.stats.g4.expected",
            vec![
                "-g".into(),
                "4".into(),
                "-t".into(),
                fixture("11.stats.targets").into_os_string(),
                fixture("11_target.sam").into_os_string(),
            ],
        ),
        (
            "indexed regions",
            "11.stats.expected",
            vec![
                fixture("11_target.bam").into_os_string(),
                "ref1:10-24".into(),
                "ref1:30-46".into(),
                "ref1:39-56".into(),
            ],
        ),
        (
            "custom index",
            "11.stats.expected",
            vec![
                "-X".into(),
                fixture("11_target.bam").into_os_string(),
                fixture("11_target.bam.bai").into_os_string(),
                "ref1:10-24".into(),
                "ref1:30-46".into(),
                "ref1:39-56".into(),
            ],
        ),
        (
            "remove three-read overlap",
            "12.3reads.nooverlap.expected",
            vec![
                "-p".into(),
                "-t".into(),
                fixture("12_3reads.bed").into_os_string(),
                fixture("12_overlaps.bam").into_os_string(),
            ],
        ),
        (
            "remove two-read overlap",
            "12.2reads.nooverlap.expected",
            vec![
                "-p".into(),
                "-t".into(),
                fixture("12_2reads.bed").into_os_string(),
                fixture("12_overlaps.bam").into_os_string(),
            ],
        ),
        (
            "read group",
            "14.rg.s1.expected",
            vec![
                "-I".into(),
                "s1".into(),
                fixture("11_target.bam").into_os_string(),
            ],
        ),
        (
            "sample",
            "14.rg.Sample.expected",
            vec![
                "-I".into(),
                "Sample".into(),
                fixture("11_target.bam").into_os_string(),
            ],
        ),
    ] {
        let output = run({
            let mut command = binary();
            command.arg("stats").args(arguments);
            command
        });
        assert_report(&output, expected, label);
    }
}

#[test]
fn barcode_reports_and_large_deletion_match_samtools_1_24() {
    for (input, expected) in [
        ("13_barcodes_ok.sam", "13.barcodes.bc.ok.expected"),
        ("13_barcodes_ok_ox_bz.sam", "13.barcodes.ox.ok.expected"),
    ] {
        let output = run({
            let mut command = binary();
            command.arg("stats").arg(fixture(input));
            command
        });
        assert_report(&output, expected, input);
    }

    let output = run({
        let mut command = binary();
        command
            .args(["stats", "-r"])
            .arg(fixture("ce.fa"))
            .arg(fixture("15.big_del.sam"));
        command
    });
    assert_report(&output, "15.stats.expected", "large deletion");

    for input in [
        "13_barcodes_fail_bc_length.sam",
        "13_barcodes_fail_hyphen.sam",
        "13_barcodes_fail_qt_length.sam",
    ] {
        let output = binary().arg("stats").arg(fixture(input)).output().unwrap();
        assert!(!output.status.success(), "{input} unexpectedly succeeded");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("barcode"),
            "{input}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

#[test]
fn split_and_reference_statistics_match_samtools_1_24() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.sam");
    fs::copy(fixture("10_map_cigar.sam"), &input).unwrap();
    let prefix = directory.path().join("split");
    let output = run({
        let mut command = binary();
        command
            .args(["stats", "-S", "RG", "-P"])
            .arg(&prefix)
            .args(["-r"])
            .arg(fixture("test.fa"))
            .arg(&input);
        command
    });
    assert_report(&output, "10.stats.expected", "split main report");
    for value in ["s1_a_1", "s1_b_1"] {
        let actual = fs::read(directory.path().join(format!("split_{value}.bamstat"))).unwrap();
        assert_bytes(
            stable_body(&actual),
            &fs::read(fixture(format!(
                "10_map_cigar.sam_{value}.expected.bamstat"
            )))
            .unwrap(),
            value,
        );
    }

    for (label, expected, arguments) in [
        (
            "reference headers",
            "16.stats.expected",
            vec![fixture("11_target.sam").into_os_string()],
        ),
        (
            "reference fasta",
            "17.stats.expected",
            vec![
                "-r".into(),
                fixture("test1.fa").into_os_string(),
                fixture("11_target.sam").into_os_string(),
            ],
        ),
        (
            "reference chunks",
            "17.stats.expected",
            vec![
                "--ref-stats-chunk".into(),
                "-1".into(),
                "-r".into(),
                fixture("test1.fa").into_os_string(),
                fixture("11_target.sam").into_os_string(),
            ],
        ),
        (
            "reference region",
            "18.stats.expected",
            vec![
                "-r".into(),
                fixture("test1.fa").into_os_string(),
                fixture("11_target.bam").into_os_string(),
                "alpha:10-20".into(),
            ],
        ),
        (
            "reference targets",
            "19.stats.expected",
            vec![
                "-r".into(),
                fixture("test1.fa").into_os_string(),
                "-t".into(),
                fixture("11.stats.targets").into_os_string(),
                fixture("11_target.sam").into_os_string(),
            ],
        ),
    ] {
        let output = run({
            let mut command = binary();
            command.arg("stats").arg("--ref-stats").args(arguments);
            command
        });
        assert_bytes(
            &rfs_lines(&output.stdout),
            &fs::read(fixture(expected)).unwrap(),
            label,
        );
    }
}

#[test]
fn custom_index_and_output_alias_fail_without_overwriting_input() {
    let missing_index = binary()
        .args(["stats", "-X"])
        .arg(fixture("11_target.bam"))
        .arg(fixture("missing.bai"))
        .arg("ref1:1-10")
        .output()
        .unwrap();
    assert!(!missing_index.status.success());

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.sam");
    fs::copy(fixture("1_map_cigar.sam"), &input).unwrap();
    let before = fs::read(&input).unwrap();
    let aliased = binary()
        .arg("stats")
        .args(["-o"])
        .arg(&input)
        .arg(&input)
        .output()
        .unwrap();
    assert!(!aliased.status.success());
    assert_eq!(fs::read(input).unwrap(), before);
}

#[test]
fn json_is_separate_from_text_and_serializes_split_keys() {
    let directory = tempfile::tempdir().unwrap();
    let report = directory.path().join("report.txt");
    let prefix = directory.path().join("split");
    let output = run({
        let mut command = binary();
        command
            .arg("--json")
            .args(["stats", "-S", "RG", "-P"])
            .arg(&prefix)
            .args(["-o"])
            .arg(&report)
            .args(["-r"])
            .arg(fixture("test.fa"))
            .arg(fixture("10_map_cigar.sam"));
        command
    });
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["result"]["command"], "stats");
    assert!(value["result"]["report"]["splits"]["s1_a_1"].is_object());
    let stats = &value["result"]["report"]["stats"];
    assert!(stats["coverage"].is_object());
    assert_eq!(
        stats["coverage"]["histogram"].is_array(),
        stats["sorted"] == true
    );
    assert!(stats["coverage"].get("changes").is_none());
    assert!(stats.get("previous_coordinate").is_none());
    assert_bytes(
        stable_body(&fs::read(report).unwrap()),
        &fs::read(fixture("10.stats.expected")).unwrap(),
        "named stats report",
    );
}

#[test]
fn malformed_and_mismatched_inputs_fail_loudly() {
    let directory = tempfile::tempdir().unwrap();
    let malformed = directory.path().join("malformed.sam");
    fs::write(&malformed, b"not a SAM header\n").unwrap();

    let truncated = directory.path().join("truncated.bam");
    let mut bytes = fs::read(fixture("11_target.bam")).unwrap();
    bytes.truncate(bytes.len() / 2);
    fs::write(&truncated, bytes).unwrap();

    let targets = directory.path().join("targets.bed");
    fs::write(&targets, b"unknown\t0\t10\n").unwrap();

    for (label, arguments) in [
        ("malformed header", vec![malformed.into_os_string()]),
        ("truncated BAM", vec![truncated.into_os_string()]),
        (
            "missing reference",
            vec![
                "-r".into(),
                fixture("ce.fa").into_os_string(),
                fixture("1_map_cigar.sam").into_os_string(),
            ],
        ),
        (
            "unknown target",
            vec![
                "-t".into(),
                targets.into_os_string(),
                fixture("1_map_cigar.sam").into_os_string(),
            ],
        ),
        (
            "missing split tag",
            vec![
                "-S".into(),
                "RG".into(),
                fixture("unsorted.sam").into_os_string(),
            ],
        ),
    ] {
        let output = binary().arg("stats").args(arguments).output().unwrap();
        assert!(!output.status.success(), "{label} unexpectedly succeeded");
        assert!(!output.stderr.is_empty(), "{label} had no diagnostic");
    }
}

#[test]
fn output_failures_and_adversarial_insert_limit_are_bounded() {
    let directory = tempfile::tempdir().unwrap();
    let directory_output = binary()
        .arg("stats")
        .args(["-o"])
        .arg(directory.path())
        .arg(fixture("1_map_cigar.sam"))
        .output()
        .unwrap();
    assert!(!directory_output.status.success());

    let mut child = binary()
        .arg("stats")
        .arg(fixture("1_map_cigar.sam"))
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    drop(child.stdout.take());
    assert!(!child.wait().unwrap().success());

    let output = run({
        let mut command = binary();
        command
            .args(["stats", "-i", "1000000000", "-r"])
            .arg(fixture("test.fa"))
            .arg(fixture("1_map_cigar.sam"));
        command
    });
    assert_report(&output, "1.stats.expected", "large insert-size limit");
}

#[test]
fn adversarial_coverage_and_split_cardinality_fail_loudly() {
    for bins in [
        "0,1000000,1".to_owned(),
        format!("0,{},{}", usize::MAX, usize::MAX),
    ] {
        let output = binary()
            .args(["stats", "-c", &bins])
            .arg(fixture("1_map_cigar.sam"))
            .output()
            .unwrap();
        assert!(
            !output.status.success(),
            "coverage {bins} unexpectedly succeeded"
        );
        assert!(!output.stderr.is_empty());
    }

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("many-splits.sam");
    let mut sam = String::from("@HD\tVN:1.6\n");
    for index in 0..=4096 {
        writeln!(sam, "r{index}\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\tZX:Z:v{index}").unwrap();
    }
    fs::write(&input, sam).unwrap();
    let output = binary()
        .args(["stats", "-S", "ZX"])
        .arg(input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("more than 4096"));
}

#[test]
#[ignore = "requires samtools 1.24"]
fn live_samtools_1_24_oracle_matches_representative_surface() {
    let version = Command::new("samtools").arg("--version").output().unwrap();
    assert!(version.status.success());
    assert!(version.stdout.starts_with(b"samtools 1.24\n"));

    for arguments in [
        vec![
            "-r".into(),
            fixture("test.fa").into_os_string(),
            fixture("1_map_cigar.sam").into_os_string(),
        ],
        vec![
            "-r".into(),
            fixture("test.fa").into_os_string(),
            fixture("7_supp.sam").into_os_string(),
        ],
        vec![fixture("13_barcodes_ok.sam").into_os_string()],
        vec![
            "-t".into(),
            fixture("11.stats.targets").into_os_string(),
            fixture("11_target.sam").into_os_string(),
        ],
        vec![
            fixture("11_target.bam").into_os_string(),
            "ref1:10-24".into(),
            "ref1:30-46".into(),
            "ref1:39-56".into(),
        ],
        vec![
            "--ref-stats".into(),
            "-r".into(),
            fixture("test1.fa").into_os_string(),
            fixture("11_target.sam").into_os_string(),
        ],
        vec![
            "-@".into(),
            "2".into(),
            "-r".into(),
            fixture("test.fa").into_os_string(),
            fixture("1_map_cigar.cram").into_os_string(),
        ],
    ] {
        let ours = run({
            let mut command = binary();
            command.arg("stats").args(&arguments);
            command
        });
        let oracle = Command::new("samtools")
            .arg("stats")
            .args(&arguments)
            .output()
            .unwrap();
        assert!(
            oracle.status.success(),
            "arguments={arguments:?}\nstderr={}",
            String::from_utf8_lossy(&oracle.stderr)
        );
        assert_bytes(
            stable_body(&ours.stdout),
            stable_body(&oracle.stdout),
            &format!("{arguments:?}"),
        );
    }
}
