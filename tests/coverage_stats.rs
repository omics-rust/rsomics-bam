use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use flate2::Compression;
use flate2::write::GzEncoder;

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
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

fn run_failure(mut command: Command) -> Output {
    let output = command.output().unwrap();
    assert!(
        !output.status.success(),
        "command unexpectedly succeeded\nstdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn run_with_stdin(mut command: Command, input: &[u8]) -> Output {
    command
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let mut child = command.spawn().unwrap();
    child.stdin.take().unwrap().write_all(input).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "stdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn indexed_bam(directory: &Path) -> PathBuf {
    let bam = directory.join("records.bam");
    run({
        let mut command = binary();
        command
            .args(["view", "--no-pg", "-b", "-o"])
            .arg(&bam)
            .arg(fixture("records.sam"));
        command
    });
    run({
        let mut command = binary();
        command.arg("index").arg(&bam);
        command
    });
    bam
}

#[test]
fn bedcov_preserves_bed_rows_and_appends_coverage() {
    let directory = tempfile::tempdir().unwrap();
    let bam = indexed_bam(directory.path());
    let bed = directory.path().join("regions.bed");
    fs::write(
        &bed,
        b"chr1\t0\t10\tfirst\nchr1\t15\t25\tsecond\nchr1\t0\t25\tall\n",
    )
    .unwrap();

    let output = run({
        let mut command = binary();
        command.arg("bedcov").arg(&bed).arg(&bam);
        command
    });
    assert_eq!(
        output.stdout,
        b"chr1\t0\t10\tfirst\t8\nchr1\t15\t25\tsecond\t8\nchr1\t0\t25\tall\t16\n"
    );
}

#[test]
fn bedcov_matches_gap_filter_threshold_flag_and_depth_semantics() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("bedcov.bam");
    run({
        let mut command = binary();
        command
            .args(["view", "--no-pg", "-b", "-o"])
            .arg(&bam)
            .arg(fixture("bedcov.sam"));
        command
    });
    run({
        let mut command = binary();
        command.arg("index").arg(&bam);
        command
    });

    let cases: &[(&[&str], &[u8])] = &[
        (&[], b"chr1\t0\t13\tall\t31\nchr1\t0\t5\tleft\t15\n"),
        (&["-j"], b"chr1\t0\t13\tall\t25\nchr1\t0\t5\tleft\t15\n"),
        (
            &["-Q", "20"],
            b"chr1\t0\t13\tall\t26\nchr1\t0\t5\tleft\t10\n",
        ),
        (
            &["-d", "2"],
            b"chr1\t0\t13\tall\t25\t10\nchr1\t0\t5\tleft\t15\t5\n",
        ),
        (
            &["-g", "DUP"],
            b"chr1\t0\t13\tall\t36\nchr1\t0\t5\tleft\t20\n",
        ),
        (
            &["-c"],
            b"chr1\t0\t13\tall\t31\t3\nchr1\t0\t5\tleft\t15\t3\n",
        ),
        (
            &["--max-depth", "1", "-c"],
            b"chr1\t0\t13\tall\t13\t1\nchr1\t0\t5\tleft\t5\t1\n",
        ),
        (
            &["-@", "1"],
            b"chr1\t0\t13\tall\t31\nchr1\t0\t5\tleft\t15\n",
        ),
    ];
    for (options, expected) in cases {
        let output = run({
            let mut command = binary();
            command
                .arg("bedcov")
                .args(*options)
                .arg(fixture("bedcov.bed"))
                .arg(&bam);
            command
        });
        assert_eq!(&output.stdout, expected, "{options:?}");
    }
}

#[test]
fn bedcov_dense_sweep_preserves_unsorted_row_order() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("bedcov.bam");
    run({
        let mut command = binary();
        command
            .args(["view", "--no-pg", "-b", "-o"])
            .arg(&bam)
            .arg(fixture("bedcov.sam"));
        command
    });
    run({
        let mut command = binary();
        command.arg("index").arg(&bam);
        command
    });
    let bed = directory.path().join("dense.bed");
    let mut rows = String::new();
    let mut expected = String::new();
    for index in 0..300 {
        if index % 2 == 0 {
            rows.push_str(&format!("chr1\t0\t13\tr{index}\n"));
            expected.push_str(&format!("chr1\t0\t13\tr{index}\t31\n"));
        } else {
            rows.push_str(&format!("chr1\t0\t5\tr{index}\n"));
            expected.push_str(&format!("chr1\t0\t5\tr{index}\t15\n"));
        }
    }
    fs::write(&bed, rows).unwrap();

    let output = run({
        let mut command = binary();
        command.arg("bedcov").arg("-@").arg("1").arg(&bed).arg(&bam);
        command
    });
    assert_eq!(output.stdout, expected.as_bytes());
}

#[test]
fn bedcov_supports_headers_gzip_and_explicit_indices() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("bedcov.bam");
    run({
        let mut command = binary();
        command
            .args(["view", "--no-pg", "-b", "-o"])
            .arg(&bam)
            .arg(fixture("bedcov.sam"));
        command
    });
    run({
        let mut command = binary();
        command.arg("index").arg(&bam);
        command
    });
    let default_index = PathBuf::from(format!("{}.bai", bam.display()));
    let custom_index = directory.path().join("custom.index");
    fs::rename(default_index, &custom_index).unwrap();
    let gzip_bed = directory.path().join("regions.bed.gz");
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder
        .write_all(&fs::read(fixture("bedcov.bed")).unwrap())
        .unwrap();
    fs::write(&gzip_bed, encoder.finish().unwrap()).unwrap();

    let output = run({
        let mut command = binary();
        command
            .arg("bedcov")
            .arg("-H")
            .arg("-X")
            .arg(&gzip_bed)
            .arg(&bam)
            .arg(&custom_index);
        command
    });
    let expected_header = format!(
        "#chrom\tchromStart\tchromEnd\tname\t{}_cov\n",
        bam.display()
    );
    assert!(output.stdout.starts_with(expected_header.as_bytes()));
    assert!(
        output
            .stdout
            .ends_with(b"chr1\t0\t13\tall\t31\nchr1\t0\t5\tleft\t15\n")
    );
}

#[test]
fn coverage_summaries_use_shared_json_and_preserve_transactions() {
    let directory = tempfile::tempdir().unwrap();
    let bam = indexed_bam(directory.path());
    let coverage_output = directory.path().join("coverage.tsv");
    let coverage = run({
        let mut command = binary();
        command
            .arg("--json")
            .arg("coverage")
            .arg("-o")
            .arg(&coverage_output)
            .arg(&bam);
        command
    });
    let value: serde_json::Value = serde_json::from_slice(&coverage.stdout).unwrap();
    assert_eq!(value["result"]["command"], "coverage");
    assert_eq!(value["result"]["report"]["references"][0]["reads"], 2);
    assert!(coverage_output.is_file());

    let bed_output = directory.path().join("bedcov.tsv");
    let bedcov = run({
        let mut command = binary();
        command
            .arg("--json")
            .arg("bedcov")
            .arg("-o")
            .arg(&bed_output)
            .arg(fixture("bedcov.bed"))
            .arg(&bam);
        command
    });
    let value: serde_json::Value = serde_json::from_slice(&bedcov.stdout).unwrap();
    assert_eq!(value["result"]["command"], "bedcov");
    assert_eq!(value["result"]["summary"]["regions"], 2);
    assert!(bed_output.is_file());

    let malformed = directory.path().join("malformed.bed");
    fs::write(&malformed, b"chr1\tbad\t10\n").unwrap();
    fs::write(&bed_output, b"keep\n").unwrap();
    run_failure({
        let mut command = binary();
        command
            .arg("bedcov")
            .arg("-o")
            .arg(&bed_output)
            .arg(&malformed)
            .arg(&bam);
        command
    });
    assert_eq!(fs::read(&bed_output).unwrap(), b"keep\n");
}

#[test]
fn coverage_emits_the_complete_reference_table() {
    let output = run({
        let mut command = binary();
        command.arg("coverage").arg(fixture("records.sam"));
        command
    });
    assert_eq!(
        output.stdout,
        b"#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\nchr1\t1\t40\t2\t16\t40\t0.4\t40\t60\n"
    );
}

#[test]
fn coverage_applies_record_base_and_depth_filters() {
    let cases: &[(&[&str], &[u8])] = &[
        (&["-l", "9"], b"chr1\t1\t40\t0\t0\t0\t0\t0\t0\n"),
        (&["-q", "61"], b"chr1\t1\t40\t0\t0\t0\t0\t0\t0\n"),
        (&["-Q", "41"], b"chr1\t1\t40\t2\t0\t0\t0\t0\t60\n"),
        (
            &["--rf", "PAIRED"],
            b"chr1\t1\t40\t2\t16\t40\t0.4\t40\t60\n",
        ),
        (&["--ff", "READ2"], b"chr1\t1\t40\t1\t8\t20\t0.2\t40\t60\n"),
        (&["--min-depth", "2"], b"chr1\t1\t40\t2\t0\t0\t0\t0\t60\n"),
        (&["-d", "1"], b"chr1\t1\t40\t2\t16\t40\t0.4\t40\t60\n"),
    ];
    for (options, expected) in cases {
        let output = run({
            let mut command = binary();
            command
                .arg("coverage")
                .arg("-H")
                .args(*options)
                .arg(fixture("records.sam"));
            command
        });
        assert_eq!(&output.stdout, expected, "{options:?}");
    }
}

#[test]
fn coverage_merges_inputs_and_supports_lists_stdin_and_regions() {
    let directory = tempfile::tempdir().unwrap();
    let bam = indexed_bam(directory.path());
    let list = directory.path().join("inputs.txt");
    fs::write(&list, format!("{}\n{}\n", bam.display(), bam.display())).unwrap();

    let expected_multi = b"chr1\t1\t40\t4\t16\t40\t0.8\t40\t60\n";
    let positional = run({
        let mut command = binary();
        command.arg("coverage").arg("-H").arg(&bam).arg(&bam);
        command
    });
    let listed = run({
        let mut command = binary();
        command.arg("coverage").arg("-H").arg("-b").arg(&list);
        command
    });
    assert_eq!(positional.stdout, expected_multi);
    assert_eq!(listed.stdout, expected_multi);

    let stdin = run_with_stdin(
        {
            let mut command = binary();
            command.arg("coverage").arg("-H").arg("-");
            command
        },
        &fs::read(fixture("records.sam")).unwrap(),
    );
    assert_eq!(stdin.stdout, b"chr1\t1\t40\t2\t16\t40\t0.4\t40\t60\n");

    let region = run({
        let mut command = binary();
        command
            .arg("coverage")
            .arg("-H")
            .args(["-r", "chr1:5-20"])
            .arg(&bam);
        command
    });
    assert_eq!(region.stdout, b"chr1\t5\t20\t2\t8\t50\t0.5\t40\t60\n");
}

#[test]
fn idxstats_reads_index_metadata_and_unplaced_counts() {
    let directory = tempfile::tempdir().unwrap();
    let bam = indexed_bam(directory.path());
    let output = run({
        let mut command = binary();
        command.arg("idxstats").arg(&bam);
        command
    });
    assert_eq!(output.stdout, b"chr1\t40\t2\t0\n*\t0\t0\t1\n");
}

#[test]
fn idxstats_scans_a_coordinate_sorted_file_without_an_index() {
    let output = run({
        let mut command = binary();
        command.arg("idxstats").arg(fixture("records.sam"));
        command
    });
    assert_eq!(output.stdout, b"chr1\t40\t2\t0\n*\t0\t0\t1\n");
}

#[test]
fn idxstats_accepts_an_explicit_index() {
    let directory = tempfile::tempdir().unwrap();
    let bam = indexed_bam(directory.path());
    let default_index = PathBuf::from(format!("{}.bai", bam.display()));
    let custom_index = directory.path().join("custom.index");
    fs::rename(default_index, &custom_index).unwrap();

    let output = run({
        let mut command = binary();
        command
            .arg("idxstats")
            .arg("-X")
            .arg(&bam)
            .arg(&custom_index);
        command
    });
    assert_eq!(output.stdout, b"chr1\t40\t2\t0\n*\t0\t0\t1\n");
}

#[test]
fn idxstats_rejects_a_broken_explicit_index_instead_of_scanning() {
    let directory = tempfile::tempdir().unwrap();
    let bam = indexed_bam(directory.path());
    let broken = directory.path().join("broken.bai");
    fs::write(&broken, b"not an index").unwrap();

    let output = run_failure({
        let mut command = binary();
        command.arg("idxstats").arg("-X").arg(&bam).arg(&broken);
        command
    });
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("index"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn idxstats_rejects_an_unsorted_scan() {
    let directory = tempfile::tempdir().unwrap();
    let sam = directory.path().join("unsorted.sam");
    fs::write(
        &sam,
        b"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:40\n@SQ\tSN:chr2\tLN:40\na\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\nb\t0\tchr2\t1\t60\t1M\t*\t0\t0\tA\tI\nc\t0\tchr1\t2\t60\t1M\t*\t0\t0\tA\tI\n",
    )
    .unwrap();

    let output = run_failure({
        let mut command = binary();
        command.arg("idxstats").arg(&sam);
        command
    });
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("coordinate sorted"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn idxstats_uses_the_shared_json_envelope() {
    let output = run({
        let mut command = binary();
        command
            .arg("--json")
            .arg("idxstats")
            .arg(fixture("records.sam"));
        command
    });
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["result"]["command"], "idxstats");
    assert_eq!(value["result"]["report"]["references"][0]["mapped"], 2);
    assert_eq!(value["result"]["report"]["unplaced_unmapped"], 1);
}

#[test]
fn idxstats_preserves_named_output_on_failure_and_rejects_input_aliases() {
    let directory = tempfile::tempdir().unwrap();
    let bam = indexed_bam(directory.path());
    let output_path = directory.path().join("stats.tsv");
    let broken = directory.path().join("broken.bai");
    fs::write(&output_path, b"keep\n").unwrap();
    fs::write(&broken, b"not an index").unwrap();

    run_failure({
        let mut command = binary();
        command
            .arg("idxstats")
            .arg("-X")
            .arg("-o")
            .arg(&output_path)
            .arg(&bam)
            .arg(&broken);
        command
    });
    assert_eq!(fs::read(&output_path).unwrap(), b"keep\n");

    let before = fs::read(&bam).unwrap();
    run_failure({
        let mut command = binary();
        command.arg("idxstats").arg("-o").arg(&bam).arg(&bam);
        command
    });
    assert_eq!(fs::read(&bam).unwrap(), before);
}
