use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn appended_extension(path: &Path, extension: &str) -> PathBuf {
    let mut path = path.as_os_str().to_owned();
    path.push(".");
    path.push(extension);
    PathBuf::from(path)
}

fn run(mut command: Command) -> Output {
    let output = command.output().expect("spawn command");
    assert!(
        output.status.success(),
        "status={:?}\nstdout={}\nstderr={}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn run_ours(arguments: &[&str]) -> Output {
    let mut command = binary();
    command.args(arguments);
    run(command)
}

fn run_samtools(arguments: &[&str]) -> Output {
    let mut command = Command::new("samtools");
    command.args(arguments);
    run(command)
}

fn run_with_stdin(mut command: Command, input: &[u8]) -> Output {
    command
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let mut child = command.spawn().expect("spawn command");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input)
        .expect("write stdin");
    let output = child.wait_with_output().expect("wait for command");
    assert!(
        output.status.success(),
        "status={:?}\nstdout={}\nstderr={}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn assert_samtools_1_24() {
    let output = run_samtools(&["--version"]);
    let version = String::from_utf8(output.stdout).unwrap();
    assert!(version.starts_with("samtools 1.24\n"), "{version}");
}

struct AlignmentSet {
    sam: PathBuf,
    bam: PathBuf,
    cram: PathBuf,
    reference: PathBuf,
}

fn build_alignment_set(directory: &Path) -> AlignmentSet {
    let sam = golden("records.sam");
    let reference = directory.join("reference.fa");
    let bam = directory.join("records.bam");
    let cram = directory.join("records.cram");
    fs::copy(golden("reference.fa"), &reference).unwrap();

    run_samtools(&["faidx", reference.to_str().unwrap()]);
    run_samtools(&[
        "view",
        "-b",
        "-o",
        bam.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]);
    run_samtools(&[
        "view",
        "-C",
        "-T",
        reference.to_str().unwrap(),
        "-o",
        cram.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]);

    AlignmentSet {
        sam,
        bam,
        cram,
        reference,
    }
}

#[test]
fn flagstat_matches_committed_samtools_text() {
    let input = golden("flagstat-small.bam");
    let output = run_ours(&["flagstat", input.to_str().unwrap()]);
    assert_eq!(
        output.stdout,
        fs::read(golden("flagstat-small.txt")).unwrap()
    );
}

#[test]
fn flagstat_accepts_sam() {
    let input = golden("records.sam");
    let output = run_ours(&["flagstat", input.to_str().unwrap()]);
    let text = String::from_utf8(output.stdout).unwrap();
    assert!(text.starts_with("3 + 0 in total"), "{text}");
    assert!(text.contains("2 + 0 mapped (66.67% : N/A)"), "{text}");
}

#[test]
fn truncated_alignment_fails_loudly() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("truncated.bam");
    let mut bytes = fs::read(golden("flagstat-small.bam")).unwrap();
    bytes.truncate(bytes.len() / 2);
    fs::write(&input, bytes).unwrap();

    let output = binary().args(["flagstat"]).arg(input).output().unwrap();
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("reading alignment record"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn view_rejects_malformed_raw_bam_records() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("malformed.bam");
    let output = directory.path().join("output.bam");

    let mut writer = noodles::bam::io::Writer::new(File::create(&input).unwrap());
    writer
        .write_header(&noodles::sam::Header::default())
        .unwrap();
    writer.get_mut().write_all(&31u32.to_le_bytes()).unwrap();
    writer.get_mut().write_all(&[0; 31]).unwrap();
    writer.try_finish().unwrap();

    for arguments in [
        vec!["view", "-c", input.to_str().unwrap()],
        vec![
            "view",
            "-b",
            "-o",
            output.to_str().unwrap(),
            input.to_str().unwrap(),
        ],
    ] {
        let result = binary().args(arguments).output().unwrap();
        assert!(!result.status.success());
        assert!(
            String::from_utf8_lossy(&result.stderr).contains("invalid BAM record"),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
    assert!(!output.exists());
}

#[test]
fn machine_output_uses_the_shared_envelope() {
    let input = golden("records.sam");
    let output = run_ours(&["--json", "flagstat", input.to_str().unwrap()]);
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["status"], "ok");
    assert_eq!(value["tool"], "rsomics-bam");
    assert_eq!(value["result"]["command"], "flagstat");
    assert_eq!(
        value["result"]["counts"]["total"],
        serde_json::json!([3, 0])
    );
}

#[test]
fn domain_and_envelope_json_are_not_mixed() {
    let input = golden("records.sam");
    let output = binary()
        .args(["--json", "flagstat", "--output-fmt", "json"])
        .arg(input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stderr).unwrap();
    assert_eq!(value["status"], "error");
}

#[test]
fn head_defaults_to_header_only() {
    let input = golden("records.sam");
    let output = run_ours(&["head", input.to_str().unwrap()]);
    let text = String::from_utf8(output.stdout).unwrap();
    assert_eq!(
        text,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:40\n@RG\tID:rg1\tSM:sample-a\n"
    );
}

#[test]
fn head_limits_headers_and_records() {
    let input = golden("records.sam");
    let output = run_ours(&["head", "-H", "1", "-n", "2", input.to_str().unwrap()]);
    let text = String::from_utf8(output.stdout).unwrap();
    assert_eq!(text.lines().count(), 3, "{text}");
    assert!(text.starts_with("@HD\tVN:1.6\tSO:coordinate\n"), "{text}");
    assert!(text.contains("read1\t99\tchr1\t1\t60\t8M"), "{text}");
    assert!(text.contains("read1\t147\tchr1\t17\t60\t8M"), "{text}");
}

#[test]
fn head_reads_alignment_from_stdin() {
    let mut command = binary();
    command.args(["head", "-n", "1", "-"]);
    let output = run_with_stdin(command, &fs::read(golden("records.sam")).unwrap());
    let text = String::from_utf8(output.stdout).unwrap();
    assert_eq!(text.lines().count(), 4, "{text}");
    assert!(
        text.ends_with("read1\t99\tchr1\t1\t60\t8M\t=\t17\t24\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n")
    );
}

#[test]
fn head_rejects_json_stream_mixing() {
    let input = golden("records.sam");
    let output = binary()
        .args(["--json", "head"])
        .arg(input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stderr).unwrap();
    assert_eq!(value["status"], "error");
}

#[test]
fn view_streams_sam_body_and_counts() {
    let input = golden("records.sam");
    let output = run_ours(&["view", input.to_str().unwrap()]);
    let text = String::from_utf8(output.stdout).unwrap();
    assert_eq!(text.lines().count(), 3, "{text}");
    assert!(!text.starts_with('@'), "{text}");

    let output = run_ours(&["view", "-c", "-f", "paired", input.to_str().unwrap()]);
    assert_eq!(output.stdout, b"2\n");
}

#[test]
fn view_controls_program_header_provenance() {
    let input = golden("records.sam");

    let default = run_ours(&["view", "-h", input.to_str().unwrap()]);
    let default = String::from_utf8(default.stdout).unwrap();
    let program = default
        .lines()
        .find(|line| line.starts_with("@PG\tID:rsomics-bam\t"))
        .unwrap();
    assert!(program.contains("\tPN:rsomics-bam\t"));
    assert!(program.contains(concat!("\tVN:", env!("CARGO_PKG_VERSION"), "\t")));
    assert!(program.contains("\tCL:"));

    let suppressed = run_ours(&["view", "-h", "--no-PG", input.to_str().unwrap()]);
    let suppressed = String::from_utf8(suppressed.stdout).unwrap();
    assert!(!suppressed.lines().any(|line| line.starts_with("@PG")));
}

#[test]
fn view_reads_alignment_from_stdin() {
    let mut command = binary();
    command.args(["view", "-c", "-"]);
    let output = run_with_stdin(command, &fs::read(golden("records.sam")).unwrap());
    assert_eq!(output.stdout, b"3\n");
}

#[test]
fn view_count_uses_the_shared_json_envelope() {
    let input = golden("records.sam");
    let output = run_ours(&["--json", "view", "-c", input.to_str().unwrap()]);
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["status"], "ok");
    assert_eq!(value["result"]["command"], "view");
    assert_eq!(value["result"]["summary"]["selected"], 3);
    assert_eq!(value["result"]["summary"]["rejected"], 0);
}

#[test]
fn view_saves_filter_counts_transactionally() {
    let directory = tempfile::tempdir().unwrap();
    let counts = directory.path().join("counts.json");
    let input = golden("records.sam");
    let output = run_ours(&[
        "view",
        "-c",
        "-f",
        "paired",
        "--save-counts",
        counts.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(output.stdout, b"2\n");

    let value: serde_json::Value = serde_json::from_slice(&fs::read(&counts).unwrap()).unwrap();
    assert_eq!(value["records_processed"], 3);
    assert_eq!(value["records_filter_accepted"], 2);
    assert_eq!(value["records_filter_rejected"], 1);

    fs::write(&counts, b"existing\n").unwrap();
    let invalid = directory.path().join("invalid.sam");
    fs::write(
        &invalid,
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\nbroken\t0\tchr1\n",
    )
    .unwrap();
    let failed = binary()
        .args(["view", "--save-counts"])
        .arg(&counts)
        .arg(&invalid)
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert_eq!(fs::read(&counts).unwrap(), b"existing\n");
}

#[test]
fn view_rejects_ambiguous_count_outputs() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("output.sam");
    let input = golden("records.sam");

    for arguments in [
        vec![
            "view",
            "-o",
            output.to_str().unwrap(),
            "--save-counts",
            output.to_str().unwrap(),
            input.to_str().unwrap(),
        ],
        vec!["view", "--save-counts", "-", input.to_str().unwrap()],
    ] {
        let result = binary().args(arguments).output().unwrap();
        assert!(!result.status.success());
    }
    assert!(!output.exists());
}

#[test]
fn view_rejects_json_stream_mixing() {
    let input = golden("records.sam");
    let output = binary()
        .args(["--json", "view"])
        .arg(input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stderr).unwrap();
    assert_eq!(value["status"], "error");
}

#[test]
fn view_commits_file_output_only_after_success() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("output.sam");
    fs::write(&output, b"existing\n").unwrap();

    let invalid = directory.path().join("invalid.sam");
    fs::write(
        &invalid,
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\nbroken\t0\tchr1\n",
    )
    .unwrap();
    let failed = binary()
        .args(["view", "-o"])
        .arg(&output)
        .arg(&invalid)
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert_eq!(fs::read(&output).unwrap(), b"existing\n");

    let succeeded = binary()
        .args(["view", "-o"])
        .arg(&output)
        .arg(golden("records.sam"))
        .output()
        .unwrap();
    assert!(succeeded.status.success());
    assert_eq!(fs::read_to_string(&output).unwrap().lines().count(), 3);
}

#[test]
fn view_writes_finished_bam_by_flag_and_extension() {
    let directory = tempfile::tempdir().unwrap();
    let input = golden("records.sam");

    for (output, options) in [
        (directory.path().join("inferred.bam"), Vec::new()),
        (directory.path().join("explicit.data"), vec!["-O", "bam"]),
    ] {
        let mut command = binary();
        command
            .arg("view")
            .args(options)
            .args(["-o"])
            .arg(&output)
            .arg(&input);
        let created = command.output().unwrap();
        assert!(
            created.status.success(),
            "{}",
            String::from_utf8_lossy(&created.stderr)
        );

        let decoded = run_ours(&["view", output.to_str().unwrap()]);
        assert_eq!(
            decoded.stdout,
            fs::read_to_string(&input)
                .unwrap()
                .lines()
                .filter(|line| !line.starts_with('@'))
                .map(|line| format!("{line}\n"))
                .collect::<String>()
                .as_bytes()
        );
        assert!(
            run_ours(&["quickcheck", output.to_str().unwrap()])
                .stdout
                .is_empty()
        );
    }
}

#[test]
fn view_controls_bam_compression() {
    let directory = tempfile::tempdir().unwrap();
    let input = golden("records.sam");
    let mut outputs = Vec::new();

    for (name, option) in [
        ("default.bam", "-b"),
        ("fast.bam", "-1"),
        ("uncompressed.bam", "-u"),
    ] {
        let output = directory.path().join(name);
        run_ours(&[
            "view",
            option,
            "-o",
            output.to_str().unwrap(),
            input.to_str().unwrap(),
        ]);
        assert_eq!(
            run_ours(&["view", output.to_str().unwrap()]).stdout,
            run_ours(&["view", input.to_str().unwrap()]).stdout
        );
        outputs.push(output);
    }

    assert!(fs::metadata(&outputs[2]).unwrap().len() > fs::metadata(&outputs[0]).unwrap().len());

    let rejected = binary()
        .args(["view", "-u", "-O", "sam"])
        .arg(input)
        .output()
        .unwrap();
    assert!(!rejected.status.success());
    assert!(
        String::from_utf8_lossy(&rejected.stderr)
            .contains("BAM compression options cannot be combined with SAM output")
    );
}

#[test]
fn view_writes_parallel_bam_within_the_thread_budget() {
    let directory = tempfile::tempdir().unwrap();
    let input = golden("records.sam");
    let bam = directory.path().join("parallel.bam");
    run_ours(&[
        "view",
        "-@",
        "2",
        "-b",
        "-o",
        bam.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(
        run_ours(&["view", bam.to_str().unwrap()]).stdout,
        run_ours(&["view", input.to_str().unwrap()]).stdout
    );
    assert!(
        run_ours(&["quickcheck", bam.to_str().unwrap()])
            .stdout
            .is_empty()
    );

    let index = noodles::bam::fs::index(&bam).unwrap();
    noodles::bam::bai::fs::write(appended_extension(&bam, "bai"), &index).unwrap();
    let region = directory.path().join("region.bam");
    run_ours(&[
        "view",
        "-@",
        "2",
        "-b",
        "-o",
        region.to_str().unwrap(),
        bam.to_str().unwrap(),
        "chr1:17-24",
    ]);
    let records = String::from_utf8(run_ours(&["view", region.to_str().unwrap()]).stdout).unwrap();
    assert_eq!(records.lines().count(), 1);
    assert!(records.starts_with("read1\t147\tchr1\t17\t"));
}

#[test]
fn view_queries_indexed_regions_in_order() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("records.bam");
    run_ours(&[
        "view",
        "-b",
        "-o",
        bam.to_str().unwrap(),
        golden("records.sam").to_str().unwrap(),
    ]);
    let index = noodles::bam::fs::index(&bam).unwrap();
    noodles::bam::bai::fs::write(bam.with_extension("bai"), &index).unwrap();

    let output = run_ours(&[
        "view",
        bam.to_str().unwrap(),
        "chr1:1-20",
        "chr1:17-24",
        "*",
    ]);
    let rows = String::from_utf8(output.stdout).unwrap();
    let names_and_positions = rows
        .lines()
        .map(|line| {
            let mut fields = line.split('\t');
            (
                fields.next().unwrap().to_owned(),
                fields.nth(2).unwrap().to_owned(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        names_and_positions,
        [
            ("read1".to_owned(), "1".to_owned()),
            ("read1".to_owned(), "17".to_owned()),
            ("read1".to_owned(), "17".to_owned()),
            ("unmapped".to_owned(), "0".to_owned()),
        ]
    );
}

#[test]
fn view_header_only_does_not_require_a_region_index() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("records.bam");
    run_ours(&[
        "view",
        "-b",
        "-o",
        bam.to_str().unwrap(),
        golden("records.sam").to_str().unwrap(),
    ]);

    let output = run_ours(&["view", "-H", bam.to_str().unwrap(), "chr1:1-8"]);
    assert!(
        String::from_utf8(output.stdout)
            .unwrap()
            .starts_with("@HD\t")
    );
}

#[test]
fn view_region_requires_an_index() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("records.bam");
    run_ours(&[
        "view",
        "-b",
        "-o",
        bam.to_str().unwrap(),
        golden("records.sam").to_str().unwrap(),
    ]);

    let output = binary()
        .args(["view"])
        .arg(&bam)
        .arg("chr1:1-8")
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("opening indexed alignment input"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn quickcheck_accepts_sam_and_bam() {
    for input in [golden("records.sam"), golden("flagstat-small.bam")] {
        let output = run_ours(&["quickcheck", input.to_str().unwrap()]);
        assert!(output.stdout.is_empty());
        assert!(output.stderr.is_empty());
    }
}

#[test]
fn quickcheck_detects_missing_bam_eof() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("missing-eof.bam");
    let mut bytes = fs::read(golden("flagstat-small.bam")).unwrap();
    bytes.truncate(bytes.len() - 28);
    fs::write(&input, bytes).unwrap();

    let output = binary()
        .args(["quickcheck", "-v"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        format!("{}\n", input.display())
    );
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("was missing EOF block when one should be present"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn quickcheck_unmapped_allows_header_without_targets() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("unmapped.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\nread1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n",
    )
    .unwrap();

    let rejected = binary().args(["quickcheck"]).arg(&input).output().unwrap();
    assert!(!rejected.status.success());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("had no targets in header"));

    let accepted = run_ours(&["quickcheck", "-u", input.to_str().unwrap()]);
    assert!(accepted.stdout.is_empty());
    assert!(accepted.stderr.is_empty());
}

#[test]
fn quickcheck_quiet_suppresses_per_input_warning() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("empty");
    fs::write(&input, []).unwrap();

    let output = binary()
        .args(["quickcheck", "-q"])
        .arg(input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(
        !stderr.contains("could not be opened for reading"),
        "{stderr}"
    );
    assert!(
        stderr.contains("quickcheck failed for 1 of 1 inputs"),
        "{stderr}"
    );
}

#[test]
fn quickcheck_uses_the_shared_json_envelope() {
    let valid = run_ours(&[
        "--json",
        "quickcheck",
        golden("records.sam").to_str().unwrap(),
    ]);
    let value: serde_json::Value = serde_json::from_slice(&valid.stdout).unwrap();
    assert_eq!(value["status"], "ok");
    assert_eq!(value["result"]["command"], "quickcheck");
    assert_eq!(
        value["result"]["report"]["files"][0]["problems"],
        serde_json::json!([])
    );

    let directory = tempfile::tempdir().unwrap();
    let invalid = directory.path().join("invalid");
    fs::write(&invalid, b"not alignment data").unwrap();
    let failed = binary()
        .args(["--json", "quickcheck"])
        .arg(invalid)
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert!(failed.stdout.is_empty());
    let value: serde_json::Value = serde_json::from_slice(&failed.stderr).unwrap();
    assert_eq!(value["status"], "error");
}

#[test]
fn samples_reads_sam_headers_and_custom_tags() {
    let input = golden("records.sam");
    let output = run_ours(&["samples", input.to_str().unwrap()]);
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        format!("sample-a\t{}\n", input.display())
    );

    let output = run_ours(&["samples", "-T", "ID", input.to_str().unwrap()]);
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        format!("rg1\t{}\n", input.display())
    );
}

#[test]
fn samples_deduplicates_values_and_ignores_missing_tags() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("samples.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\n\
          @RG\tID:r1\tSM:zeta\n@RG\tID:r2\tSM:alpha\n\
          @RG\tID:r3\tSM:zeta\n@RG\tID:r4\n@RG\tID:r5\tSM:beta\n",
    )
    .unwrap();

    let output = run_ours(&["samples", input.to_str().unwrap()]);
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        format!("alpha\t{0}\nbeta\t{0}\nzeta\t{0}\n", input.display())
    );
}

#[test]
fn samples_uses_the_shared_json_envelope() {
    let input = golden("records.sam");
    let output = run_ours(&["--json", "samples", input.to_str().unwrap()]);
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["status"], "ok");
    assert_eq!(value["result"]["command"], "samples");
    assert_eq!(value["result"]["report"]["entries"][0]["value"], "sample-a");
}

#[test]
fn samples_reads_input_paths_from_stdin() {
    let input = golden("records.sam");
    let mut command = binary();
    command.arg("samples");
    let output = run_with_stdin(command, format!("{}\n", input.display()).as_bytes());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        format!("sample-a\t{}\n", input.display())
    );
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn flags_match_samtools_1_24() {
    assert_samtools_1_24();
    for value in [
        "0",
        "16",
        "020",
        "0x10",
        "paired",
        "paired,read1",
        "SECONDARY,SUPPLEMENTARY",
        "4096",
    ] {
        let ours = run_ours(&["flags", value]);
        let oracle = run_samtools(&["flags", value]);
        assert_eq!(ours.stdout, oracle.stdout, "{value}");
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn flagstat_matches_samtools_1_24_for_sam_bam_and_cram() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());

    for input in [&inputs.sam, &inputs.bam, &inputs.cram] {
        for format in ["text", "tsv", "json"] {
            let mut our_arguments =
                vec!["flagstat", "--output-fmt", format, input.to_str().unwrap()];
            if input == &inputs.cram {
                our_arguments.splice(1..1, ["--reference", inputs.reference.to_str().unwrap()]);
            }
            let ours = run_ours(&our_arguments);
            let oracle =
                run_samtools(&["flagstat", "--output-fmt", format, input.to_str().unwrap()]);
            if format == "json" {
                let ours: serde_json::Value = serde_json::from_slice(&ours.stdout).unwrap();
                let oracle: serde_json::Value = serde_json::from_slice(&oracle.stdout).unwrap();
                assert_eq!(ours, oracle, "{} {format}", input.display());
            } else {
                assert_eq!(ours.stdout, oracle.stdout, "{} {format}", input.display());
            }
        }
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn head_matches_samtools_1_24_for_sam_bam_and_cram() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());

    for input in [&inputs.sam, &inputs.bam, &inputs.cram] {
        for (our_options, oracle_options) in [
            (vec![], vec![]),
            (vec!["-H", "1"], vec!["-h", "1"]),
            (vec!["-H", "2", "-n", "2"], vec!["-h", "2", "-n", "2"]),
            (vec!["-H", "0", "-n", "3"], vec!["-h", "0", "-n", "3"]),
        ] {
            let mut our_arguments = vec!["head"];
            our_arguments.extend(our_options);
            if input == &inputs.cram {
                our_arguments.extend(["--reference", inputs.reference.to_str().unwrap()]);
            }
            our_arguments.push(input.to_str().unwrap());

            let mut oracle_arguments = vec!["head"];
            oracle_arguments.extend(oracle_options);
            if input == &inputs.cram {
                oracle_arguments.extend(["-T", inputs.reference.to_str().unwrap()]);
            }
            oracle_arguments.push(input.to_str().unwrap());

            let ours = run_ours(&our_arguments);
            let oracle = run_samtools(&oracle_arguments);
            assert_eq!(ours.stdout, oracle.stdout, "{}", input.display());
        }
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn head_stdin_matches_samtools_1_24_for_sam_bam_and_cram() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());

    for input in [&inputs.sam, &inputs.bam, &inputs.cram] {
        let bytes = fs::read(input).unwrap();

        let mut ours = binary();
        ours.args(["head", "-n", "2"]);
        if input == &inputs.cram {
            ours.args(["--reference", inputs.reference.to_str().unwrap()]);
        }
        ours.arg("-");

        let mut oracle = Command::new("samtools");
        oracle.args(["head", "-n", "2"]);
        if input == &inputs.cram {
            oracle.args(["-T", inputs.reference.to_str().unwrap()]);
        }
        oracle.arg("-");

        let ours = run_with_stdin(ours, &bytes);
        let oracle = run_with_stdin(oracle, &bytes);
        assert_eq!(ours.stdout, oracle.stdout, "{}", input.display());
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn head_restores_cram_md_and_nm_like_samtools_1_24() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    let sam = directory.path().join("cigar-cases.sam");
    let cram = directory.path().join("cigar-cases.cram");
    fs::copy(golden("reference.fa"), &reference).unwrap();
    fs::write(
        &sam,
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\n\
mismatch\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTTCGT\t*\n\
insertion\t0\tchr1\t1\t60\t4M2I4M\t*\t0\t0\tACGTTTACGT\t*\n\
deletion\t0\tchr1\t1\t60\t4M2D4M\t*\t0\t0\tACGTGTAC\t*\n\
skip\t0\tchr1\t1\t60\t4M2N4M\t*\t0\t0\tACGTGTAC\t*\n\
ambiguous\t0\tchr1\t1\t60\t4M\t*\t0\t0\tNNNN\t*\n\
equal\t0\tchr1\t1\t60\t4M\t*\t0\t0\t====\t*\n",
    )
    .unwrap();

    run_samtools(&["faidx", reference.to_str().unwrap()]);
    run_samtools(&[
        "view",
        "-C",
        "-T",
        reference.to_str().unwrap(),
        "-o",
        cram.to_str().unwrap(),
        sam.to_str().unwrap(),
    ]);

    let ours = run_ours(&[
        "head",
        "-H",
        "0",
        "-n",
        "6",
        "--reference",
        reference.to_str().unwrap(),
        cram.to_str().unwrap(),
    ]);
    let oracle = run_samtools(&[
        "head",
        "-h",
        "0",
        "-n",
        "6",
        "-T",
        reference.to_str().unwrap(),
        cram.to_str().unwrap(),
    ]);
    assert_eq!(ours.stdout, oracle.stdout);
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn view_matches_samtools_1_24_for_streaming_filters() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());

    for input in [&inputs.sam, &inputs.bam, &inputs.cram] {
        for (our_options, oracle_options) in [
            (vec!["--no-PG"], vec!["--no-PG"]),
            (vec!["--no-PG", "-h"], vec!["--no-PG", "-h"]),
            (vec!["--no-PG", "-H"], vec!["--no-PG", "-H"]),
            (vec!["--no-PG", "-c"], vec!["--no-PG", "-c"]),
            (
                vec!["--no-PG", "-c", "-f", "paired"],
                vec!["--no-PG", "-c", "-f", "paired"],
            ),
            (
                vec!["--no-PG", "-c", "-F", "unmap"],
                vec!["--no-PG", "-c", "-F", "unmap"],
            ),
            (
                vec!["--no-PG", "-c", "--rf", "read1,read2"],
                vec!["--no-PG", "-c", "--rf", "read1,read2"],
            ),
            (
                vec!["--no-PG", "-c", "-G", "paired,proper_pair"],
                vec!["--no-PG", "-c", "-G", "paired,proper_pair"],
            ),
            (
                vec!["--no-PG", "-c", "-q", "60"],
                vec!["--no-PG", "-c", "-q", "60"],
            ),
        ] {
            let mut our_arguments = vec!["view"];
            our_arguments.extend(our_options);
            if input == &inputs.cram {
                our_arguments.extend(["--reference", inputs.reference.to_str().unwrap()]);
            }
            our_arguments.push(input.to_str().unwrap());

            let mut oracle_arguments = vec!["view"];
            oracle_arguments.extend(oracle_options);
            if input == &inputs.cram {
                oracle_arguments.extend(["-T", inputs.reference.to_str().unwrap()]);
            }
            oracle_arguments.push(input.to_str().unwrap());

            let ours = run_ours(&our_arguments);
            let oracle = run_samtools(&oracle_arguments);
            assert_eq!(
                ours.stdout,
                oracle.stdout,
                "{} {our_arguments:?}",
                input.display()
            );
        }
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn view_saved_counts_match_samtools_1_24() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());

    for (position, input) in [&inputs.sam, &inputs.bam, &inputs.cram]
        .into_iter()
        .enumerate()
    {
        let ours_counts = directory.path().join(format!("ours-{position}.json"));
        let oracle_counts = directory.path().join(format!("oracle-{position}.json"));
        let mut ours = vec![
            "view",
            "--no-PG",
            "-c",
            "-f",
            "paired",
            "--save-counts",
            ours_counts.to_str().unwrap(),
        ];
        let mut oracle = vec![
            "view",
            "--no-PG",
            "-c",
            "-f",
            "paired",
            "--save-counts",
            oracle_counts.to_str().unwrap(),
        ];
        if input == &inputs.cram {
            ours.extend(["-T", inputs.reference.to_str().unwrap()]);
            oracle.extend(["-T", inputs.reference.to_str().unwrap()]);
        }
        ours.push(input.to_str().unwrap());
        oracle.push(input.to_str().unwrap());

        assert_eq!(run_ours(&ours).stdout, run_samtools(&oracle).stdout);
        let ours: serde_json::Value =
            serde_json::from_slice(&fs::read(ours_counts).unwrap()).unwrap();
        let oracle: serde_json::Value =
            serde_json::from_slice(&fs::read(oracle_counts).unwrap()).unwrap();
        assert_eq!(ours, oracle, "{}", input.display());
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn view_bam_output_matches_samtools_1_24_for_sam_bam_and_cram() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());

    for (position, input) in [&inputs.sam, &inputs.bam, &inputs.cram]
        .into_iter()
        .enumerate()
    {
        let ours_path = directory.path().join(format!("ours-{position}.bam"));
        let oracle_path = directory.path().join(format!("oracle-{position}.bam"));

        let mut ours = binary();
        ours.args(["view", "--no-PG", "-@", "2", "-b", "-f", "paired", "-o"])
            .arg(&ours_path);
        if input == &inputs.cram {
            ours.args(["-T", inputs.reference.to_str().unwrap()]);
        }
        ours.arg(input);
        run(ours);

        let mut oracle = Command::new("samtools");
        oracle
            .args(["view", "--no-PG", "-@", "2", "-b", "-f", "paired", "-o"])
            .arg(&oracle_path);
        if input == &inputs.cram {
            oracle.args(["-T", inputs.reference.to_str().unwrap()]);
        }
        oracle.arg(input);
        run(oracle);

        let ours = run_samtools(&["view", "--no-PG", "-h", ours_path.to_str().unwrap()]);
        let oracle = run_samtools(&["view", "--no-PG", "-h", oracle_path.to_str().unwrap()]);
        assert_eq!(ours.stdout, oracle.stdout, "{}", input.display());
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn view_bam_compression_modes_match_samtools_1_24_records() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let input = golden("records.sam");

    for (position, option) in ["-1", "-u"].into_iter().enumerate() {
        let ours_path = directory.path().join(format!("ours-{position}.bam"));
        let oracle_path = directory.path().join(format!("oracle-{position}.bam"));
        run_ours(&[
            "view",
            "--no-PG",
            option,
            "-o",
            ours_path.to_str().unwrap(),
            input.to_str().unwrap(),
        ]);
        run_samtools(&[
            "view",
            "--no-PG",
            option,
            "-o",
            oracle_path.to_str().unwrap(),
            input.to_str().unwrap(),
        ]);

        let ours = run_samtools(&["view", "--no-PG", "-h", ours_path.to_str().unwrap()]);
        let oracle = run_samtools(&["view", "--no-PG", "-h", oracle_path.to_str().unwrap()]);
        assert_eq!(ours.stdout, oracle.stdout, "{option}");
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn view_regions_match_samtools_1_24_for_bam_and_cram() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());
    run_samtools(&["index", inputs.bam.to_str().unwrap()]);
    let alternative_crai = inputs.cram.with_extension("crai");
    run_samtools(&[
        "index",
        "-o",
        alternative_crai.to_str().unwrap(),
        inputs.cram.to_str().unwrap(),
    ]);

    for input in [&inputs.bam, &inputs.cram] {
        let mut ours = vec!["view"];
        let mut oracle = vec!["view", "--no-PG"];
        if input == &inputs.cram {
            ours.extend(["-T", inputs.reference.to_str().unwrap()]);
            oracle.extend(["-T", inputs.reference.to_str().unwrap()]);
        }
        ours.push(input.to_str().unwrap());
        oracle.push(input.to_str().unwrap());
        ours.extend(["chr1:1-20", "chr1:17-24", "*"]);
        oracle.extend(["chr1:1-20", "chr1:17-24", "*"]);

        assert_eq!(
            run_ours(&ours).stdout,
            run_samtools(&oracle).stdout,
            "{}",
            input.display()
        );
    }

    fs::rename(
        appended_extension(&inputs.bam, "bai"),
        directory.path().join("saved.bai"),
    )
    .unwrap();
    run_samtools(&["index", "-c", inputs.bam.to_str().unwrap()]);
    let ours = run_ours(&["view", inputs.bam.to_str().unwrap(), "chr1:1-8"]);
    let oracle = run_samtools(&["view", "--no-PG", inputs.bam.to_str().unwrap(), "chr1:1-8"]);
    assert_eq!(ours.stdout, oracle.stdout);
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn quickcheck_matches_samtools_1_24_decisions() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());
    let no_targets = directory.path().join("unmapped.sam");
    fs::write(
        &no_targets,
        b"@HD\tVN:1.6\nread1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n",
    )
    .unwrap();
    let missing_eof = directory.path().join("missing-eof.bam");
    let mut bytes = fs::read(&inputs.bam).unwrap();
    bytes.truncate(bytes.len() - 28);
    fs::write(&missing_eof, bytes).unwrap();
    let invalid = directory.path().join("invalid");
    fs::write(&invalid, b"not alignment data").unwrap();

    for input in [
        &inputs.sam,
        &inputs.bam,
        &inputs.cram,
        &no_targets,
        &missing_eof,
        &invalid,
    ] {
        for options in [Vec::new(), vec!["-u"], vec!["-q"], vec!["-q", "-v"]] {
            let mut ours = binary();
            ours.arg("quickcheck").args(&options).arg(input);
            let ours = ours.output().unwrap();

            let mut oracle = Command::new("samtools");
            oracle.arg("quickcheck").args(&options).arg(input);
            let oracle = oracle.output().unwrap();

            assert_eq!(
                ours.status.success(),
                oracle.status.success(),
                "{} {options:?}\nours={}\noracle={}",
                input.display(),
                String::from_utf8_lossy(&ours.stderr),
                String::from_utf8_lossy(&oracle.stderr)
            );
            if options.contains(&"-v") {
                assert_eq!(ours.stdout, oracle.stdout, "{}", input.display());
            }
        }
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn samples_matches_samtools_1_24_for_sam_bam_cram_and_metadata() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());
    let bam_index = directory.path().join("records.bam.bai");
    run_samtools(&[
        "index",
        "-o",
        bam_index.to_str().unwrap(),
        inputs.bam.to_str().unwrap(),
    ]);

    for input in [&inputs.sam, &inputs.bam, &inputs.cram] {
        for (ours, oracle) in [
            (Vec::new(), Vec::new()),
            (vec!["-H"], vec!["-h"]),
            (vec!["-T", "ID"], vec!["-T", "ID"]),
        ] {
            let mut our_arguments = vec!["samples"];
            our_arguments.extend(ours);
            our_arguments.push(input.to_str().unwrap());
            let ours = run_ours(&our_arguments);

            let mut oracle_arguments = vec!["samples"];
            oracle_arguments.extend(oracle);
            oracle_arguments.push(input.to_str().unwrap());
            let oracle = run_samtools(&oracle_arguments);
            assert_eq!(ours.stdout, oracle.stdout, "{}", input.display());
        }
    }

    let ours = run_ours(&[
        "samples",
        "-H",
        "-i",
        "-f",
        inputs.reference.to_str().unwrap(),
        "-X",
        inputs.bam.to_str().unwrap(),
        bam_index.to_str().unwrap(),
    ]);
    let oracle = run_samtools(&[
        "samples",
        "-h",
        "-i",
        "-f",
        inputs.reference.to_str().unwrap(),
        "-X",
        inputs.bam.to_str().unwrap(),
        bam_index.to_str().unwrap(),
    ]);
    assert_eq!(ours.stdout, oracle.stdout);

    let reference_list = directory.path().join("references.txt");
    fs::write(&reference_list, format!("{}\n", inputs.reference.display())).unwrap();
    let ours = run_ours(&[
        "samples",
        "-H",
        "-F",
        reference_list.to_str().unwrap(),
        inputs.bam.to_str().unwrap(),
    ]);
    let oracle = run_samtools(&[
        "samples",
        "-h",
        "-F",
        reference_list.to_str().unwrap(),
        inputs.bam.to_str().unwrap(),
    ]);
    assert_eq!(ours.stdout, oracle.stdout);

    let mut ours = binary();
    ours.args(["samples", "-i", "-X"]);
    let ours = run_with_stdin(
        ours,
        format!("{}\t{}\n", inputs.bam.display(), bam_index.display()).as_bytes(),
    );
    let mut oracle = Command::new("samtools");
    oracle.args(["samples", "-i", "-X"]);
    let oracle = run_with_stdin(
        oracle,
        format!("{}\t{}\n", inputs.bam.display(), bam_index.display()).as_bytes(),
    );
    assert_eq!(ours.stdout, oracle.stdout);

    let many = directory.path().join("many.sam");
    let mut header = String::from("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\n");
    for position in 0..20 {
        header.push_str(&format!(
            "@RG\tID:r{position}\tSM:sample-{}\tLB:library-{position}\n",
            position % 13
        ));
    }
    fs::write(&many, header).unwrap();
    for tag in ["SM", "LB"] {
        let ours = run_ours(&["samples", "-T", tag, many.to_str().unwrap()]);
        let oracle = run_samtools(&["samples", "-T", tag, many.to_str().unwrap()]);
        assert_eq!(ours.stdout, oracle.stdout, "tag={tag}");
    }
}
