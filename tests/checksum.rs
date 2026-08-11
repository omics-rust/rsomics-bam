use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const EXPECTED_DEFAULT: &str = "# Checksum 1.0 for file: <input>\n\
# Aux tags:          BC,FI,QT,RT,TC\n\
# BAM flags:         PAIRED,READ1,READ2\n\
\n\
# Group    QC          count  flag+seq  +name     +qual     +aux      combined\n\
all        all             3  2edf95c7  15d129dc  6f9d052a  2edf95c7  0e252c96\n\
-          all             1  50bd7e68  0909945c  5dd2abc7  50bd7e68  1b46c159\n\
rg1        all             2  3355f51b  7202f1e0  42002f0d  3355f51b  75c869ca\n";

const EXPECTED_MERGE: &str = "# Checksum 1.0 for file: merge\n\
# Aux tags:          BC,FI,QT,RT,TC\n\
# BAM flags:         PAIRED,READ1,READ2\n\
\n\
# Group    QC          count  flag+seq  +name     +qual     +aux      combined\n\
all        all             6  5d7d7ae4  7bc82c96  22f1c078  5d7d7ae4  0dcab84b\n\
-          all             2  16e880ca  3d15daa4  2084f10b  16e880ca  3506a33e\n\
rg1        all             4  053169e4  1194c2e9  5eadf72e  053169e4  07b8ee95\n";

const SANITIZE_INPUT: &str = "@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:10\n\
r1\t0\tchr1\t1\t42\t1=1X2M\t*\t0\t0\tACGT\tIIII\tNM:i:1\tMD:Z:1A2\n\
r2\t4\tchr1\t2\t30\t4M\t*\t0\t0\tTGCA\tJJJJ\tNM:i:0\tMD:Z:4\n";

const EXPECTED_SANITIZED: &str = "# Checksum 1.0 for file: <input>\n\
# Aux tags:          *,cF,MD,NM\n\
# BAM flags:         PAIRED,PROPER_PAIR,UNMAP,MUNMAP,REVERSE,MREVERSE,READ1,READ2,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY\n\
\n\
# Group    QC          count  flag+seq  +name     +qual     +aux      +chr/pos  +cigar    +mate     combined\n\
all        all             2  11c180db  53aefa63  5e20c336  11c180db  50643658  63f7ab40  14e7cfdf  298f32c7\n\
-          all             2  11c180db  53aefa63  5e20c336  11c180db  50643658  63f7ab40  14e7cfdf  298f32c7\n";

const EXPECTED_FASTQ: &str = "# Checksum 1.0 for file: <input>\n\
# Aux tags:          BC,FI,QT,RT,TC\n\
# BAM flags:         PAIRED,READ1,READ2\n\
\n\
# Group    QC          count  flag+seq  +name     +qual     +aux      combined\n\
all        all             2  5f26eee0  24d401e5  15e3883b  5f26eee0  34fe13f0\n\
-          all             2  5f26eee0  24d401e5  15e3883b  5f26eee0  34fe13f0\n";

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn corpus(name: &str) -> PathBuf {
    fixture("checksum").join(name)
}

fn run(arguments: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(arguments)
        .output()
        .unwrap()
}

fn run_with_stdin(arguments: &[&str], input: &[u8]) -> std::process::Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(arguments)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(input).unwrap();
    child.wait_with_output().unwrap()
}

fn normalized_report(output: &[u8], input: &Path) -> String {
    String::from_utf8(output.to_vec())
        .unwrap()
        .replace(&input.display().to_string(), "<input>")
}

fn normalized_source_header(output: &[u8]) -> Vec<u8> {
    let Some(newline) = output.iter().position(|&byte| byte == b'\n') else {
        return output.to_vec();
    };
    if !output.starts_with(b"# Checksum 1.0 for file:") {
        return output.to_vec();
    }
    let mut normalized = b"# Checksum 1.0 for file:\n".to_vec();
    normalized.extend_from_slice(&output[newline + 1..]);
    normalized
}

#[test]
fn checksum_matches_the_default_samtools_report() {
    let input = fixture("records.sam");
    let output = run(&["checksum", input.to_str().unwrap()]);

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(normalized_report(&output.stdout, &input), EXPECTED_DEFAULT);
}

#[test]
fn checksum_help_exposes_the_stable_surface() {
    let output = run(&["checksum", "--help"]);
    assert!(output.status.success());
    let help = String::from_utf8(output.stdout).unwrap();
    for option in [
        "-F, --exclude-flags",
        "-f, --require-flags",
        "-b, --flag-mask",
        "-t, --tags",
        "-O, --in-order",
        "-P, --check-pos",
        "-C, --check-cigar",
        "-M, --check-mate",
        "-z, --sanitize",
        "-N, --count",
        "-o, --output",
        "-q, --show-qc",
        "-v, --verbose",
        "-T, --tabs",
        "-m, --merge",
        "-a, --all",
        "-B, --bamseqchksum",
        "-@, --threads",
    ] {
        assert!(help.contains(option), "missing {option} in {help}");
    }
}

#[test]
fn merge_combines_native_reports() {
    let directory = tempfile::tempdir().unwrap();
    let input = fixture("records.sam");
    let first = directory.path().join("first.chk");
    let second = directory.path().join("second.chk");
    for report in [&first, &second] {
        let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
            .args(["checksum", "-o"])
            .arg(report)
            .arg(&input)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["checksum", "-m"])
        .arg(first)
        .arg(second)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, EXPECTED_MERGE.as_bytes());
}

#[test]
fn all_applies_the_samtools_sanitization_contract() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("sanitize.sam");
    std::fs::write(&input, SANITIZE_INPUT).unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["checksum", "--all"])
        .arg(&input)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        normalized_report(&output.stdout, &input),
        EXPECTED_SANITIZED
    );
}

#[test]
fn fastq_pair_suffixes_match_htslib_record_semantics() {
    let input = fixture("import-r1.fastq");
    let output = run(&["checksum", input.to_str().unwrap()]);

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(normalized_report(&output.stdout, &input), EXPECTED_FASTQ);
}

#[test]
fn fastq_standard_input_is_detected_without_an_extension() {
    let input = std::fs::read(fixture("import-r1.fastq")).unwrap();
    let output = run_with_stdin(&["checksum", "-"], &input);

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        EXPECTED_FASTQ.replace("<input>", "-")
    );
}

#[test]
fn named_fastq_is_detected_by_content() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("reads");
    std::fs::copy(fixture("import-r1.fastq"), &input).unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["checksum"])
        .arg(&input)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(normalized_report(&output.stdout, &input), EXPECTED_FASTQ);
}

#[test]
fn compressed_fastq_is_detected_by_content_for_files_and_stdin() {
    let input = std::fs::read(fixture("import-r1.fastq")).unwrap();
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(&input).unwrap();
    let compressed = encoder.finish().unwrap();

    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("reads.data");
    std::fs::write(&path, &compressed).unwrap();
    let named = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["checksum"])
        .arg(&path)
        .output()
        .unwrap();
    assert!(
        named.status.success(),
        "{}",
        String::from_utf8_lossy(&named.stderr)
    );
    assert_eq!(normalized_report(&named.stdout, &path), EXPECTED_FASTQ);

    let streamed = run_with_stdin(&["checksum", "-"], &compressed);
    assert!(
        streamed.status.success(),
        "{}",
        String::from_utf8_lossy(&streamed.stderr)
    );
    assert_eq!(
        String::from_utf8(streamed.stdout).unwrap(),
        EXPECTED_FASTQ.replace("<input>", "-")
    );
}

#[test]
fn sam_and_bam_standard_input_remain_alignment_streams() {
    let sam = std::fs::read(fixture("records.sam")).unwrap();
    let output = run_with_stdin(&["checksum", "-"], &sam);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        EXPECTED_DEFAULT.replace("<input>", "-")
    );

    let bam = std::fs::read(corpus("chk1.bam")).unwrap();
    let output = run_with_stdin(&["checksum", "-"], &bam);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        normalized_source_header(&output.stdout),
        std::fs::read(corpus("chk1.1.expected")).unwrap()
    );
}

#[test]
fn complete_upstream_compute_corpus_matches_samtools_1_24() {
    let cases: &[(&[&str], &str, &str)] = &[
        (&[], "chk1.bam", "chk1.1.expected"),
        (&["-qv"], "chk1.bam", "chk1.3.expected"),
        (&["-B"], "chk1.bam", "chk1.4.expected"),
        (&[], "chk2.cram", "chk2.1.expected"),
        (&["-a"], "chk2.cram", "chk2.2.expected"),
        (&["-qv"], "chk2.cram", "chk2.3.expected"),
        (&["-qv", "-a"], "chk2.cram", "chk2.4.expected"),
    ];
    for &(arguments, input_name, expected_name) in cases {
        let input = corpus(input_name);
        let mut command = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"));
        command.arg("checksum").args(arguments).arg(&input);
        let output = command.output().unwrap();
        assert!(
            output.status.success(),
            "{input_name}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let expected = std::fs::read(corpus(expected_name)).unwrap();
        assert_eq!(
            normalized_source_header(&output.stdout),
            expected,
            "{input_name} {arguments:?}"
        );
    }
}

#[test]
fn upstream_native_and_bamseq_merge_cases_match() {
    let cases: &[(&[&str], &[&str], &str)] = &[
        (
            &[],
            &["chk1.1.expected", "chk1.4.expected"],
            "chk1.7.expected",
        ),
        (
            &["-B"],
            &["chk1.1.expected", "chk1.4.expected"],
            "chk1.8.expected",
        ),
        (&["-B"], &["chk1.4.expected"], "chk1.4.expected"),
        (
            &["-B"],
            &["chk1.4.expected", "chk1.5.expected"],
            "chk1.6.expected",
        ),
    ];
    for &(arguments, inputs, expected_name) in cases {
        let mut command = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"));
        command.arg("checksum").args(arguments).arg("-m");
        for input in inputs {
            command.arg(corpus(input));
        }
        let output = command.output().unwrap();
        assert!(
            output.status.success(),
            "{inputs:?}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let expected = std::fs::read(corpus(expected_name)).unwrap();
        assert_eq!(
            normalized_source_header(&output.stdout),
            normalized_source_header(&expected),
            "{inputs:?} {arguments:?}"
        );
    }
}

#[test]
fn failures_do_not_replace_named_reports() {
    let directory = tempfile::tempdir().unwrap();
    let output_path = directory.path().join("output.chk");
    std::fs::write(&output_path, b"sentinel\n").unwrap();
    let malformed = directory.path().join("malformed.sam");
    std::fs::write(
        &malformed,
        b"@HD\tVN:1.6\nbroken\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tbad\n",
    )
    .unwrap();
    let failed = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["checksum", "-o"])
        .arg(&output_path)
        .arg(&malformed)
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert_eq!(std::fs::read(&output_path).unwrap(), b"sentinel\n");

    let malformed_auxiliary = directory.path().join("malformed-auxiliary.sam");
    std::fs::write(
        &malformed_auxiliary,
        b"@HD\tVN:1.6\nread\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tXY:Q:bad\n",
    )
    .unwrap();
    let failed = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["checksum", "-o"])
        .arg(&output_path)
        .arg(&malformed_auxiliary)
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert_eq!(std::fs::read(&output_path).unwrap(), b"sentinel\n");

    let incompatible = directory.path().join("incompatible.chk");
    let report = std::fs::read_to_string(corpus("chk1.1.expected"))
        .unwrap()
        .replace("BC,FI,QT,RT,TC", "AM");
    std::fs::write(&incompatible, report).unwrap();
    let failed = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["checksum", "-m", "-o"])
        .arg(&output_path)
        .arg(corpus("chk1.1.expected"))
        .arg(incompatible)
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert_eq!(std::fs::read(&output_path).unwrap(), b"sentinel\n");
}

#[test]
fn json_keeps_the_compatibility_report_separate() {
    let directory = tempfile::tempdir().unwrap();
    let report_path = directory.path().join("report.chk");
    let input = fixture("records.sam");
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["--json", "checksum", "-o"])
        .arg(&report_path)
        .arg(&input)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["result"]["command"], "checksum");
    assert_eq!(value["result"]["reports"][0]["groups"][0]["name"], "all");
    assert_eq!(
        value["result"]["reports"][0]["groups"][0]["rows"][0]["count"],
        3
    );
    assert_eq!(
        normalized_report(&std::fs::read(report_path).unwrap(), &input),
        EXPECTED_DEFAULT
    );
}

#[test]
fn incompatible_output_modes_and_empty_merge_fail_loud() {
    let input = fixture("records.sam");
    for arguments in [
        vec!["checksum", "-m"],
        vec!["checksum", "-B", "-a", input.to_str().unwrap()],
        vec!["checksum", "-a", "-F", "0x100", input.to_str().unwrap()],
        vec!["checksum", "-@", "257", input.to_str().unwrap()],
        vec!["--json", "checksum", input.to_str().unwrap()],
    ] {
        let output = run(&arguments);
        assert!(
            !output.status.success(),
            "{arguments:?}: {}",
            String::from_utf8_lossy(&output.stdout)
        );
    }
}

#[test]
fn malformed_merge_contracts_are_rejected() {
    let original = std::fs::read_to_string(corpus("chk1.1.expected")).unwrap();
    let cases = [
        (
            "repeated-version",
            format!("# Checksum 1.0 for file: duplicate\n{original}"),
        ),
        (
            "unsupported-version",
            original.replacen("# Checksum 1.0", "# Checksum 2.0", 1),
        ),
        (
            "missing-tags",
            without_line_prefix(&original, "# Aux tags:"),
        ),
        (
            "missing-flags",
            without_line_prefix(&original, "# BAM flags:"),
        ),
        (
            "repeated-columns",
            original.replacen(
                "# Group    QC          count  flag+seq  +name     +qual     +aux      combined\n",
                "# Group    QC          count  flag+seq  +name     +qual     +aux      combined\n# Group    QC          count  flag+seq  +name     +qual     +aux      combined\n",
                1,
            ),
        ),
        (
            "unknown-column",
            original.replacen("+aux      combined", "+aux      +unknown  combined", 1),
        ),
        (
            "bad-combined",
            original.replacen("435f5683", "00000001", 1),
        ),
        (
            "missing-component",
            without_line_prefix(&original, "ERR013140"),
        ),
        (
            "duplicate-tags",
            original.replacen("BC,FI,QT,RT,TC", "BC,BC", 1),
        ),
        (
            "duplicate-row",
            format!(
                "{original}{}",
                original
                    .lines()
                    .find(|line| line.starts_with("ERR013140"))
                    .unwrap()
            ),
        ),
        (
            "hybrid-header",
            format!(
                "# Aux tags: BC,FI,QT,RT,TC\n{}",
                std::fs::read_to_string(corpus("chk1.4.expected")).unwrap()
            ),
        ),
    ];

    for (name, report) in cases {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join(format!("{name}.chk"));
        std::fs::write(&path, report).unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
            .args(["checksum", "-m"])
            .arg(path)
            .output()
            .unwrap();
        assert!(
            !output.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&output.stdout)
        );
    }
}

fn without_line_prefix(input: &str, prefix: &str) -> String {
    input
        .split_inclusive('\n')
        .filter(|line| !line.starts_with(prefix))
        .collect()
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn live_samtools_1_24_oracle_matches_the_stable_surface() {
    let version = Command::new("samtools").arg("--version").output().unwrap();
    assert_eq!(
        String::from_utf8_lossy(&version.stdout)
            .lines()
            .next()
            .unwrap_or_default(),
        "samtools 1.24"
    );

    let cases: &[(&[&str], &str)] = &[
        (&[], "chk1.bam"),
        (&["-qv"], "chk1.bam"),
        (&["-B"], "chk1.bam"),
        (&["-T"], "chk1.bam"),
        (&["-O"], "chk1.bam"),
        (&["-OO"], "chk1.bam"),
        (&["-P"], "chk1.bam"),
        (&["-C"], "chk1.bam"),
        (&["-M"], "chk1.bam"),
        (&["-t", "*,cF,MD,NM"], "chk1.bam"),
        (&["-z", "all,cigarx"], "chk1.bam"),
        (
            &["-F", "0x100", "-f", "0x1", "-b", "0xfff", "-c", "-N", "3"],
            "chk1.bam",
        ),
        (&["-a"], "chk1.bam"),
        (&["-@", "2"], "chk1.bam"),
        (&[], "chk2.cram"),
        (&["-a"], "chk2.cram"),
        (&["-qvT"], "chk2.cram"),
    ];

    for &(arguments, input_name) in cases {
        assert_live_oracle(arguments, &corpus(input_name));
    }
    for (arguments, input) in [
        (&[][..], fixture("import-r1.fastq")),
        (&["-a"][..], fixture("import-iupac.fastq")),
        (&["-qvT"][..], fixture("import-r1.fastq")),
        (&[][..], fixture("reference.fa")),
        (&["-a"][..], fixture("reference.fa")),
    ] {
        assert_live_oracle(arguments, &input);
    }

    let long_cigar_directory = tempfile::tempdir().unwrap();
    let long_cigar = long_cigar_fixture(long_cigar_directory.path());
    for arguments in [&["-C"][..], &["-t", "*"][..], &["-a"][..]] {
        assert_live_oracle(arguments, &long_cigar);
    }

    let merge_cases: &[(&[&str], &[&str])] = &[
        (&[], &["chk1.1.expected", "chk1.4.expected"]),
        (&["-B"], &["chk1.1.expected", "chk1.4.expected"]),
        (&["-B"], &["chk1.4.expected", "chk1.5.expected"]),
    ];
    for &(arguments, inputs) in merge_cases {
        let mut ours = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"));
        ours.arg("checksum").args(arguments).arg("-m");
        let mut oracle = Command::new("samtools");
        oracle.arg("checksum").args(arguments).arg("-m");
        for input in inputs {
            ours.arg(corpus(input));
            oracle.arg(corpus(input));
        }
        let ours = ours.output().unwrap();
        let oracle = oracle.output().unwrap();
        assert_success_and_equal(&ours, &oracle, &format!("merge {arguments:?} {inputs:?}"));
    }
}

fn long_cigar_fixture(directory: &Path) -> PathBuf {
    let sam = directory.join("long.sam");
    let bam = directory.join("long.bam");
    let cigar = "1M1I".repeat(32_768);
    let sequence = "A".repeat(65_536);
    std::fs::write(
        &sam,
        format!(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40000\nlong\t0\tchr1\t1\t60\t{cigar}\t*\t0\t0\t{sequence}\t*\n"
        ),
    )
    .unwrap();
    let output = Command::new("samtools")
        .args(["view", "-b", "-o"])
        .arg(&bam)
        .arg(&sam)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    bam
}

fn assert_live_oracle(arguments: &[&str], input: &Path) {
    let ours = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .arg("checksum")
        .args(arguments)
        .arg(input)
        .output()
        .unwrap();
    let oracle = Command::new("samtools")
        .arg("checksum")
        .args(arguments)
        .arg(input)
        .output()
        .unwrap();
    assert_success_and_equal(
        &ours,
        &oracle,
        &format!("{arguments:?} {}", input.display()),
    );
}

fn assert_success_and_equal(
    ours: &std::process::Output,
    oracle: &std::process::Output,
    label: &str,
) {
    assert!(
        ours.status.success(),
        "{label}: rsomics: {}",
        String::from_utf8_lossy(&ours.stderr)
    );
    assert!(
        oracle.status.success(),
        "{label}: samtools: {}",
        String::from_utf8_lossy(&oracle.stderr)
    );
    assert_eq!(ours.stdout, oracle.stdout, "{label}");
}
