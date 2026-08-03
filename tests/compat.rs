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

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn mpileup_matches_samtools_1_24_for_sam_bam_and_cram() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());

    for input in [&inputs.sam, &inputs.bam] {
        for options in [
            Vec::new(),
            vec!["-Q", "20"],
            vec!["-A"],
            vec!["-x"],
            vec!["-a"],
            vec!["-aa"],
        ] {
            let mut ours = vec!["mpileup"];
            ours.extend(options.iter().copied());
            ours.push(input.to_str().unwrap());
            let mut oracle = vec!["mpileup"];
            oracle.extend(options.iter().copied());
            oracle.push(input.to_str().unwrap());
            assert_eq!(
                run_ours(&ours).stdout,
                run_samtools(&oracle).stdout,
                "{} {options:?}",
                input.display()
            );
        }
    }

    for input in [&inputs.sam, &inputs.bam, &inputs.cram] {
        for options in [Vec::new(), vec!["-B"], vec!["-E"]] {
            let reference = inputs.reference.to_str().unwrap();
            let mut ours = vec!["mpileup", "-f", reference];
            ours.extend(options.iter().copied());
            ours.push(input.to_str().unwrap());
            let mut oracle = vec!["mpileup", "-f", reference];
            oracle.extend(options.iter().copied());
            oracle.push(input.to_str().unwrap());
            assert_eq!(
                run_ours(&ours).stdout,
                run_samtools(&oracle).stdout,
                "{} {options:?}",
                input.display()
            );
        }
    }

    let indels = golden("mpileup-records.sam");
    let reference = golden("mpileup-reference.fa");
    for options in [
        Vec::new(),
        vec!["-Q", "0"],
        vec!["-x"],
        vec!["-f", reference.to_str().unwrap()],
        vec!["-f", reference.to_str().unwrap(), "-B"],
        vec!["-f", reference.to_str().unwrap(), "-E"],
    ] {
        let mut ours = vec!["mpileup"];
        ours.extend(options.iter().copied());
        ours.push(indels.to_str().unwrap());
        let mut oracle = vec!["mpileup"];
        oracle.extend(options.iter().copied());
        oracle.push(indels.to_str().unwrap());
        assert_eq!(
            run_ours(&ours).stdout,
            run_samtools(&oracle).stdout,
            "indels {options:?}"
        );
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn index_matches_samtools_1_24_for_bai_csi_crai_and_bgzf_sam() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());

    for (name, options) in [("bai", vec!["-@", "2"]), ("csi", vec!["-c", "-m", "1"])] {
        let ours = directory.path().join(format!("ours-{name}.bam"));
        let oracle = directory.path().join(format!("oracle-{name}.bam"));
        fs::copy(&inputs.bam, &ours).unwrap();
        fs::copy(&inputs.bam, &oracle).unwrap();

        let mut ours_arguments = vec!["index"];
        ours_arguments.extend(options.iter().copied());
        ours_arguments.push(ours.to_str().unwrap());
        run_ours(&ours_arguments);
        let mut oracle_arguments = vec!["index"];
        oracle_arguments.extend(options.iter().copied());
        oracle_arguments.push(oracle.to_str().unwrap());
        run_samtools(&oracle_arguments);

        assert_eq!(
            run_ours(&["view", ours.to_str().unwrap(), "chr1:1-8"]).stdout,
            run_samtools(&["view", oracle.to_str().unwrap(), "chr1:1-8"]).stdout,
            "{name} region"
        );
        assert_eq!(
            run_samtools(&["idxstats", ours.to_str().unwrap()]).stdout,
            run_samtools(&["idxstats", oracle.to_str().unwrap()]).stdout,
            "{name} idxstats"
        );
    }

    let ours_cram = directory.path().join("ours.cram");
    let oracle_cram = directory.path().join("oracle.cram");
    fs::copy(&inputs.cram, &ours_cram).unwrap();
    fs::copy(&inputs.cram, &oracle_cram).unwrap();
    run_ours(&["index", "-@", "2", ours_cram.to_str().unwrap()]);
    run_samtools(&["index", "-@", "2", oracle_cram.to_str().unwrap()]);
    assert_eq!(
        run_ours(&[
            "view",
            "-T",
            inputs.reference.to_str().unwrap(),
            ours_cram.to_str().unwrap(),
            "chr1:1-8",
        ])
        .stdout,
        run_samtools(&[
            "view",
            "-T",
            inputs.reference.to_str().unwrap(),
            oracle_cram.to_str().unwrap(),
            "chr1:1-8",
        ])
        .stdout
    );
    assert_eq!(
        run_samtools(&["idxstats", ours_cram.to_str().unwrap()]).stdout,
        run_samtools(&["idxstats", oracle_cram.to_str().unwrap()]).stdout
    );

    let source = directory.path().join("records.sam.gz");
    run_samtools(&[
        "view",
        "-h",
        "-O",
        "sam",
        "-o",
        source.to_str().unwrap(),
        inputs.sam.to_str().unwrap(),
    ]);
    let ours_sam = directory.path().join("ours.sam.gz");
    let oracle_sam = directory.path().join("oracle.sam.gz");
    fs::copy(&source, &ours_sam).unwrap();
    fs::copy(&source, &oracle_sam).unwrap();
    run_ours(&["index", ours_sam.to_str().unwrap()]);
    run_samtools(&["index", oracle_sam.to_str().unwrap()]);
    assert_eq!(
        run_ours(&["view", ours_sam.to_str().unwrap(), "chr1:1-8"]).stdout,
        run_samtools(&["view", oracle_sam.to_str().unwrap(), "chr1:1-8"]).stdout
    );
    assert_eq!(
        run_samtools(&["idxstats", ours_sam.to_str().unwrap()]).stdout,
        run_samtools(&["idxstats", oracle_sam.to_str().unwrap()]).stdout
    );
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
fn mpileup_named_output_uses_the_shared_json_envelope() {
    let directory = tempfile::tempdir().unwrap();
    let output_path = directory.path().join("pileup.txt");
    let output = run_ours(&[
        "--json",
        "mpileup",
        "-o",
        output_path.to_str().unwrap(),
        golden("records.sam").to_str().unwrap(),
    ]);
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["status"], "ok");
    assert_eq!(value["result"]["command"], "mpileup");
    assert_eq!(value["result"]["summary"]["positions"], 16);
    assert_eq!(
        fs::read(output_path).unwrap(),
        fs::read(golden("mpileup-default.txt")).unwrap()
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
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:40\n\
         @RG\tID:rg1\tSM:sample-a\tLB:lib-a\n"
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
fn view_program_path_remains_available() {
    let program = rsomics_bam::Program::new("tool", "1.0.0", "tool view input.bam").unwrap();
    let _: rsomics_bam::view::Program<'_> = program;
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
fn view_filters_by_query_length() {
    let input = golden("records.sam");

    assert_eq!(
        run_ours(&["view", "-c", "-m", "1", input.to_str().unwrap()]).stdout,
        b"2\n"
    );
    assert_eq!(
        run_ours(&["view", "-c", "-m", "9", input.to_str().unwrap()]).stdout,
        b"0\n"
    );
}

#[test]
fn view_filters_records_and_headers_by_read_group() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("read-groups.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\n@RG\tID:rg1\tSM:a\n@RG\tID:rg2\tSM:b\n\
          r1\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\tRG:Z:rg1\n\
          r2\t0\tchr1\t2\t60\t1M\t*\t0\t0\tC\tI\tRG:Z:rg2\n\
          r3\t0\tchr1\t3\t60\t1M\t*\t0\t0\tG\tI\n",
    )
    .unwrap();

    let selected = String::from_utf8(
        run_ours(&[
            "view",
            "-h",
            "--no-PG",
            "-r",
            "rg1",
            input.to_str().unwrap(),
        ])
        .stdout,
    )
    .unwrap();
    assert!(selected.contains("@RG\tID:rg1\t"));
    assert!(!selected.contains("@RG\tID:rg2\t"));
    assert!(selected.contains("\nr1\t"));
    assert!(!selected.contains("\nr2\t"));
    assert!(selected.contains("\nr3\t"));

    assert_eq!(
        run_ours(&[
            "view",
            "-c",
            "-r",
            "rg1",
            "-r",
            "rg2",
            input.to_str().unwrap(),
        ])
        .stdout,
        b"3\n"
    );

    let bam = directory.path().join("read-groups.bam");
    let filtered_bam = directory.path().join("filtered.bam");
    run_ours(&[
        "view",
        "--no-PG",
        "-b",
        "-o",
        bam.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    run_ours(&[
        "view",
        "--no-PG",
        "-b",
        "-r",
        "rg1",
        "-o",
        filtered_bam.to_str().unwrap(),
        bam.to_str().unwrap(),
    ]);
    assert_eq!(
        run_ours(&["view", "--no-PG", "-h", filtered_bam.to_str().unwrap()]).stdout,
        selected.as_bytes()
    );

    let invalid_sam = directory.path().join("invalid-read-group.sam");
    fs::write(
        &invalid_sam,
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\nr1\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\tRG:i:1\n",
    )
    .unwrap();

    let invalid_bam = directory.path().join("invalid-read-group.bam");
    let mut writer = noodles::bam::io::Writer::new(File::create(&invalid_bam).unwrap());
    writer
        .write_header(&noodles::sam::Header::default())
        .unwrap();
    let mut record = rsomics_bamio::raw::RawRecord::default();
    record
        .append_aux(*b"RG", b'i', &1i32.to_le_bytes())
        .unwrap();
    rsomics_bamio::raw::write_record(writer.get_mut(), &record).unwrap();
    writer.try_finish().unwrap();

    for invalid in [invalid_sam, invalid_bam] {
        let failed = binary()
            .args(["view", "-c", "-r", "rg1"])
            .arg(invalid)
            .output()
            .unwrap();
        assert!(!failed.status.success());
        assert!(
            String::from_utf8_lossy(&failed.stderr).contains("RG tag must be a string"),
            "{}",
            String::from_utf8_lossy(&failed.stderr)
        );
    }
}

#[test]
fn view_filters_by_read_names_from_files() {
    let directory = tempfile::tempdir().unwrap();
    let read1 = directory.path().join("read1.txt");
    let read2 = directory.path().join("read2.txt");
    fs::write(&read1, b"read1\tread1\n").unwrap();
    fs::write(&read2, b"unmapped\n").unwrap();
    let input = golden("records.sam");

    assert_eq!(
        run_ours(&[
            "view",
            "-c",
            "-N",
            read1.to_str().unwrap(),
            input.to_str().unwrap(),
        ])
        .stdout,
        b"2\n"
    );
    assert_eq!(
        run_ours(&[
            "view",
            "-c",
            "-N",
            read1.to_str().unwrap(),
            "-N",
            read2.to_str().unwrap(),
            input.to_str().unwrap(),
        ])
        .stdout,
        b"3\n"
    );

    let excluded = format!("^{}", read1.display());
    assert_eq!(
        run_ours(&["view", "-c", "-N", &excluded, input.to_str().unwrap(),]).stdout,
        b"1\n"
    );

    let mixed = binary()
        .args(["view", "-c", "-N"])
        .arg(&read1)
        .args(["-N", &excluded])
        .arg(&input)
        .output()
        .unwrap();
    assert!(!mixed.status.success());
    assert!(
        String::from_utf8_lossy(&mixed.stderr)
            .contains("cannot mix include and exclude read-name files")
    );

    let missing = binary()
        .args(["view", "-c", "-N"])
        .arg(directory.path().join("missing.txt"))
        .arg(&input)
        .output()
        .unwrap();
    assert!(!missing.status.success());
    assert!(String::from_utf8_lossy(&missing.stderr).contains("opening read-name file"));

    let unnamed = directory.path().join("unnamed.sam");
    let star = directory.path().join("star.txt");
    fs::write(&unnamed, b"@HD\tVN:1.6\n*\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n").unwrap();
    fs::write(&star, b"*\n").unwrap();
    assert_eq!(
        run_ours(&[
            "view",
            "-c",
            "-N",
            star.to_str().unwrap(),
            unnamed.to_str().unwrap(),
        ])
        .stdout,
        b"1\n"
    );
}

#[test]
fn view_filters_by_library_from_read_group_headers() {
    let directory = tempfile::tempdir().unwrap();
    let input = golden("records.sam");

    assert_eq!(
        run_ours(&["view", "-c", "-l", "lib-a", input.to_str().unwrap()]).stdout,
        b"2\n"
    );
    assert_eq!(
        run_ours(&["view", "-c", "-l", "missing", input.to_str().unwrap()]).stdout,
        b"0\n"
    );

    let selected = String::from_utf8(
        run_ours(&[
            "view",
            "--no-PG",
            "-h",
            "-l",
            "lib-a",
            input.to_str().unwrap(),
        ])
        .stdout,
    )
    .unwrap();
    assert!(selected.contains("@RG\tID:rg1\tSM:sample-a\tLB:lib-a\n"));
    assert_eq!(
        selected
            .lines()
            .filter(|line| !line.starts_with('@'))
            .count(),
        2
    );

    let bam = directory.path().join("records.bam");
    run_ours(&[
        "view",
        "--no-PG",
        "-b",
        "-o",
        bam.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(
        run_ours(&["view", "-c", "-l", "lib-a", bam.to_str().unwrap()]).stdout,
        b"2\n"
    );

    let invalid = directory.path().join("invalid-library.sam");
    fs::write(
        &invalid,
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\n@RG\tID:rg1\tLB:lib-a\n\
          r1\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\tRG:i:1\n",
    )
    .unwrap();
    let failed = binary()
        .args(["view", "-c", "-l", "lib-a"])
        .arg(invalid)
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert!(String::from_utf8_lossy(&failed.stderr).contains("RG tag must be a string"));
}

#[test]
fn view_changes_flags_after_filtering() {
    let input = golden("records.sam");
    let changed = String::from_utf8(
        run_ours(&[
            "view",
            "--no-PG",
            "--add-flags",
            "0x400",
            "--remove-flags",
            "0x10",
            input.to_str().unwrap(),
        ])
        .stdout,
    )
    .unwrap();
    let flags = changed
        .lines()
        .map(|line| line.split('\t').nth(1).unwrap().parse::<u16>().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(flags, [1123, 1155, 1028]);

    let filtered = run_ours(&[
        "view",
        "--no-PG",
        "-f",
        "0x04",
        "--remove-flags",
        "0x04",
        input.to_str().unwrap(),
    ]);
    assert_eq!(
        String::from_utf8(filtered.stdout)
            .unwrap()
            .split('\t')
            .nth(1),
        Some("0")
    );

    assert_eq!(
        run_ours(&[
            "view",
            "-c",
            "-f",
            "0x400",
            "--add-flags",
            "0x400",
            input.to_str().unwrap(),
        ])
        .stdout,
        b"0\n"
    );
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
fn view_keeps_or_removes_selected_auxiliary_tags() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("tags.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\n@RG\tID:rg1\n\
          r1\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\tRG:Z:rg1\tNM:i:1\tAS:i:8\tXX:Z:x\n",
    )
    .unwrap();

    let removed = run_ours(&[
        "view",
        "--no-PG",
        "-x",
        "NM",
        "--remove-tag",
        "AS",
        input.to_str().unwrap(),
    ]);
    let removed = String::from_utf8(removed.stdout).unwrap();
    assert!(removed.contains("\tRG:Z:rg1\tXX:Z:x\n"), "{removed}");
    assert!(!removed.contains("\tNM:"), "{removed}");
    assert!(!removed.contains("\tAS:"), "{removed}");

    for arguments in [
        vec!["--keep-tag", "RG,XX"],
        vec!["-x", "^RG,XX"],
        vec!["--keep-tag", "RG", "--keep-tag", "XX"],
    ] {
        let mut command = vec!["view", "--no-PG"];
        command.extend(arguments);
        command.push(input.to_str().unwrap());
        assert_eq!(run_ours(&command).stdout, removed.as_bytes());
    }

    assert_eq!(
        run_ours(&["view", "-c", "--keep-tag", "", input.to_str().unwrap(),]).stdout,
        b"1\n"
    );

    for arguments in [
        vec!["view", "-x", "NM", "--keep-tag", "RG"],
        vec!["view", "-x", "NM", "-x", "^RG"],
        vec!["view", "-x", "N"],
    ] {
        let failed = binary().args(arguments).arg(&input).output().unwrap();
        assert!(!failed.status.success());
    }

    let bam = directory.path().join("tags.bam");
    run_ours(&[
        "view",
        "--no-PG",
        "-b",
        "--keep-tag",
        "RG,XX",
        "-o",
        bam.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(
        run_ours(&["view", "--no-PG", bam.to_str().unwrap()]).stdout,
        removed.as_bytes()
    );
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
fn view_minimum_query_length_matches_samtools_1_24() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());

    for input in [&inputs.sam, &inputs.bam, &inputs.cram] {
        for minimum in ["0", "1", "8", "9"] {
            let mut ours = vec!["view", "--no-PG", "-c", "-m", minimum];
            let mut oracle = vec!["view", "--no-PG", "-c", "-m", minimum];
            if input == &inputs.cram {
                ours.extend(["-T", inputs.reference.to_str().unwrap()]);
                oracle.extend(["-T", inputs.reference.to_str().unwrap()]);
            }
            ours.push(input.to_str().unwrap());
            oracle.push(input.to_str().unwrap());

            assert_eq!(
                run_ours(&ours).stdout,
                run_samtools(&oracle).stdout,
                "{} minimum={minimum}",
                input.display()
            );
        }
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn view_read_group_filters_match_samtools_1_24() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());

    for input in [&inputs.sam, &inputs.bam, &inputs.cram] {
        for filter in [
            &["-r", "rg1"][..],
            &["-r", "missing"][..],
            &["-r", "rg1", "-r", "missing"][..],
        ] {
            let mut ours = vec!["view", "--no-PG", "-h"];
            let mut oracle = vec!["view", "--no-PG", "-h"];
            ours.extend(filter.iter().copied());
            oracle.extend(filter.iter().copied());
            if input == &inputs.cram {
                ours.extend(["-T", inputs.reference.to_str().unwrap()]);
                oracle.extend(["-T", inputs.reference.to_str().unwrap()]);
            }
            ours.push(input.to_str().unwrap());
            oracle.push(input.to_str().unwrap());

            assert_eq!(
                run_ours(&ours).stdout,
                run_samtools(&oracle).stdout,
                "{} filter={filter:?}",
                input.display()
            );

            ours.retain(|argument| *argument != "-h");
            oracle.retain(|argument| *argument != "-h");
            ours.insert(2, "-c");
            oracle.insert(2, "-c");
            assert_eq!(
                run_ours(&ours).stdout,
                run_samtools(&oracle).stdout,
                "{} count filter={filter:?}",
                input.display()
            );
        }
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn view_qname_files_match_samtools_1_24() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());
    let read1 = directory.path().join("read1.txt");
    let read2 = directory.path().join("read2.txt");
    fs::write(&read1, b"read1\n").unwrap();
    fs::write(&read2, b"unmapped\n").unwrap();
    let excluded = format!("^{}", read1.display());

    for (position, input) in [&inputs.sam, &inputs.bam, &inputs.cram]
        .into_iter()
        .enumerate()
    {
        for names in [
            vec![read1.to_str().unwrap()],
            vec![read1.to_str().unwrap(), read2.to_str().unwrap()],
            vec![excluded.as_str()],
        ] {
            let mut ours = vec!["view", "--no-PG"];
            let mut oracle = ours.clone();
            for names in &names {
                ours.extend(["-N", names]);
                oracle.extend(["-N", names]);
            }
            if input == &inputs.cram {
                ours.extend(["-T", inputs.reference.to_str().unwrap()]);
                oracle.extend(["-T", inputs.reference.to_str().unwrap()]);
            }
            ours.push(input.to_str().unwrap());
            oracle.push(input.to_str().unwrap());
            assert_eq!(
                run_ours(&ours).stdout,
                run_samtools(&oracle).stdout,
                "{} names={names:?}",
                input.display()
            );

            ours.insert(2, "-c");
            oracle.insert(2, "-c");
            assert_eq!(
                run_ours(&ours).stdout,
                run_samtools(&oracle).stdout,
                "{} count names={names:?}",
                input.display()
            );
        }

        let ours_bam = directory.path().join(format!("ours-qname-{position}.bam"));
        let oracle_bam = directory
            .path()
            .join(format!("oracle-qname-{position}.bam"));
        let mut ours = vec![
            "view",
            "--no-PG",
            "-b",
            "-N",
            read1.to_str().unwrap(),
            "-o",
            ours_bam.to_str().unwrap(),
        ];
        let mut oracle = vec![
            "view",
            "--no-PG",
            "-b",
            "-N",
            read1.to_str().unwrap(),
            "-o",
            oracle_bam.to_str().unwrap(),
        ];
        if input == &inputs.cram {
            ours.extend(["-T", inputs.reference.to_str().unwrap()]);
            oracle.extend(["-T", inputs.reference.to_str().unwrap()]);
        }
        ours.push(input.to_str().unwrap());
        oracle.push(input.to_str().unwrap());
        run_ours(&ours);
        run_samtools(&oracle);
        assert_eq!(
            run_samtools(&["view", "--no-PG", "-h", ours_bam.to_str().unwrap()]).stdout,
            run_samtools(&["view", "--no-PG", "-h", oracle_bam.to_str().unwrap()]).stdout,
            "{} BAM output",
            input.display()
        );
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn view_library_filter_matches_samtools_1_24() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());

    for (position, input) in [&inputs.sam, &inputs.bam, &inputs.cram]
        .into_iter()
        .enumerate()
    {
        let mut ours = vec!["view", "--no-PG", "-h", "-l", "lib-a"];
        let mut oracle = ours.clone();
        if input == &inputs.cram {
            ours.extend(["-T", inputs.reference.to_str().unwrap()]);
            oracle.extend(["-T", inputs.reference.to_str().unwrap()]);
        }
        ours.push(input.to_str().unwrap());
        oracle.push(input.to_str().unwrap());
        assert_eq!(
            run_ours(&ours).stdout,
            run_samtools(&oracle).stdout,
            "{}",
            input.display()
        );

        ours.retain(|argument| *argument != "-h");
        oracle.retain(|argument| *argument != "-h");
        ours.insert(2, "-c");
        oracle.insert(2, "-c");
        assert_eq!(
            run_ours(&ours).stdout,
            run_samtools(&oracle).stdout,
            "{} count",
            input.display()
        );

        let ours_bam = directory
            .path()
            .join(format!("ours-library-{position}.bam"));
        let oracle_bam = directory
            .path()
            .join(format!("oracle-library-{position}.bam"));
        ours.retain(|argument| *argument != "-c");
        oracle.retain(|argument| *argument != "-c");
        ours.splice(2..2, ["-b", "-o", ours_bam.to_str().unwrap()]);
        oracle.splice(2..2, ["-b", "-o", oracle_bam.to_str().unwrap()]);
        run_ours(&ours);
        run_samtools(&oracle);
        assert_eq!(
            run_samtools(&["view", "--no-PG", "-h", ours_bam.to_str().unwrap()]).stdout,
            run_samtools(&["view", "--no-PG", "-h", oracle_bam.to_str().unwrap()]).stdout,
            "{} BAM output",
            input.display()
        );
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn view_flag_changes_match_samtools_1_24() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());

    for (position, input) in [&inputs.sam, &inputs.bam, &inputs.cram]
        .into_iter()
        .enumerate()
    {
        let mut ours = vec![
            "view",
            "--no-PG",
            "-h",
            "--add-flags",
            "0x400",
            "--remove-flags",
            "0x10",
        ];
        let mut oracle = ours.clone();
        if input == &inputs.cram {
            ours.extend(["-T", inputs.reference.to_str().unwrap()]);
            oracle.extend(["-T", inputs.reference.to_str().unwrap()]);
        }
        ours.push(input.to_str().unwrap());
        oracle.push(input.to_str().unwrap());
        assert_eq!(
            run_ours(&ours).stdout,
            run_samtools(&oracle).stdout,
            "{}",
            input.display()
        );

        let ours_bam = directory.path().join(format!("ours-flags-{position}.bam"));
        let oracle_bam = directory
            .path()
            .join(format!("oracle-flags-{position}.bam"));
        ours.retain(|argument| *argument != "-h");
        oracle.retain(|argument| *argument != "-h");
        ours.splice(2..2, ["-b", "-o", ours_bam.to_str().unwrap()]);
        oracle.splice(2..2, ["-b", "-o", oracle_bam.to_str().unwrap()]);
        run_ours(&ours);
        run_samtools(&oracle);
        assert_eq!(
            run_samtools(&["view", "--no-PG", "-h", ours_bam.to_str().unwrap()]).stdout,
            run_samtools(&["view", "--no-PG", "-h", oracle_bam.to_str().unwrap()]).stdout,
            "{} BAM output",
            input.display()
        );
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn view_auxiliary_tag_changes_match_samtools_1_24() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let inputs = build_alignment_set(directory.path());

    for (position, input) in [&inputs.sam, &inputs.bam, &inputs.cram]
        .into_iter()
        .enumerate()
    {
        for options in [
            &["-x", "RG"][..],
            &["--keep-tag", "RG"][..],
            &["-x", "^RG"][..],
        ] {
            let mut ours = vec!["view", "--no-PG", "-h"];
            let mut oracle = ours.clone();
            ours.extend(options);
            oracle.extend(options);
            if input == &inputs.cram {
                ours.extend(["-T", inputs.reference.to_str().unwrap()]);
                oracle.extend(["-T", inputs.reference.to_str().unwrap()]);
            }
            ours.push(input.to_str().unwrap());
            oracle.push(input.to_str().unwrap());
            assert_eq!(
                run_ours(&ours).stdout,
                run_samtools(&oracle).stdout,
                "{} options={options:?}",
                input.display()
            );
        }

        let ours_bam = directory.path().join(format!("ours-tags-{position}.bam"));
        let oracle_bam = directory.path().join(format!("oracle-tags-{position}.bam"));
        let mut ours = vec![
            "view",
            "--no-PG",
            "-b",
            "-x",
            "RG",
            "-o",
            ours_bam.to_str().unwrap(),
        ];
        let mut oracle = vec![
            "view",
            "--no-PG",
            "-b",
            "-x",
            "RG",
            "-o",
            oracle_bam.to_str().unwrap(),
        ];
        if input == &inputs.cram {
            ours.extend(["-T", inputs.reference.to_str().unwrap()]);
            oracle.extend(["-T", inputs.reference.to_str().unwrap()]);
        }
        ours.push(input.to_str().unwrap());
        oracle.push(input.to_str().unwrap());
        run_ours(&ours);
        run_samtools(&oracle);
        assert_eq!(
            run_samtools(&["view", "--no-PG", "-h", ours_bam.to_str().unwrap()]).stdout,
            run_samtools(&["view", "--no-PG", "-h", oracle_bam.to_str().unwrap()]).stdout,
            "{} BAM output",
            input.display()
        );
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
