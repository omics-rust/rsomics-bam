use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn make_bam(dir: &Path, name: &str, sam: &str) -> PathBuf {
    let source = dir.join(format!("{name}.sam"));
    let output = dir.join(format!("{name}.bam"));
    std::fs::write(&source, sam).unwrap();
    let result = binary()
        .args(["view", "--no-pg", "-b", "-o"])
        .arg(&output)
        .arg(&source)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    output
}

fn render(path: &Path) -> String {
    let output = binary()
        .args(["view", "--no-pg", "-h"])
        .arg(path)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

fn assert_one_eof(path: &Path) {
    let data = std::fs::read(path).unwrap();
    assert!(data.ends_with(&BGZF_EOF));
    assert_eq!(
        data.windows(BGZF_EOF.len())
            .filter(|window| *window == BGZF_EOF)
            .count(),
        1
    );
}

fn remove_eof(path: &Path) {
    let file = std::fs::OpenOptions::new().write(true).open(path).unwrap();
    file.set_len(file.metadata().unwrap().len() - BGZF_EOF.len() as u64)
        .unwrap();
}

#[test]
fn cat_merges_read_groups_and_preserves_record_order() {
    let dir = tempfile::tempdir().unwrap();
    let a = make_bam(
        dir.path(),
        "a",
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:a\tSM:left\na1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:a\n",
    );
    let b = make_bam(
        dir.path(),
        "b",
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:b\tSM:right\nb1\t0\tchr1\t20\t60\t4M\t*\t0\t0\tTGCA\tIIII\tRG:Z:b\n",
    );
    let output = dir.path().join("cat.bam");
    let result = binary()
        .args(["--json", "cat", "--no-pg"])
        .arg(&a)
        .arg(&b)
        .args(["-o"])
        .arg(&output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let summary: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(summary["result"]["summary"]["inputs"], 2);

    let sam = render(&output);
    assert!(sam.contains("@RG\tID:a\tSM:left\n"));
    assert!(sam.contains("@RG\tID:b\tSM:right\n"));
    assert!(sam.find("a1\t").unwrap() < sam.find("b1\t").unwrap());
    assert_one_eof(&output);
}

#[test]
fn reheader_renames_references_without_reencoding_records() {
    let dir = tempfile::tempdir().unwrap();
    let input = make_bam(
        dir.path(),
        "input",
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\nr1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\nr2\t0\tchr2\t20\t60\t4M\t*\t0\t0\tTGCA\tIIII\n",
    );
    let header = dir.path().join("replacement.sam");
    std::fs::write(
        &header,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:one\tLN:1000\n@SQ\tSN:two\tLN:2000\n@CO\treplaced\n",
    )
    .unwrap();
    let output = dir.path().join("reheader.bam");
    let result = binary()
        .args(["--json", "reheader", "--no-pg"])
        .arg(&header)
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
    let summary: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(summary["result"]["summary"]["reference_sequences"], 2);

    let sam = render(&output);
    assert!(sam.contains("@SQ\tSN:one\tLN:1000\n"));
    assert!(sam.contains("@SQ\tSN:two\tLN:2000\n"));
    assert!(sam.contains("r1\t0\tone\t10\t"));
    assert!(sam.contains("r2\t0\ttwo\t20\t"));
    assert_one_eof(&output);
}

#[test]
fn cat_preflight_failures_preserve_existing_output() {
    let dir = tempfile::tempdir().unwrap();
    let good = make_bam(
        dir.path(),
        "good",
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
    );
    let different = make_bam(
        dir.path(),
        "different",
        "@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:1000\nr2\t0\tchr2\t10\t60\t4M\t*\t0\t0\tTGCA\tIIII\n",
    );
    let truncated = dir.path().join("truncated.bam");
    std::fs::copy(&good, &truncated).unwrap();
    remove_eof(&truncated);
    let output = dir.path().join("output.bam");

    for second in [&different, &truncated] {
        std::fs::write(&output, b"sentinel").unwrap();
        let result = binary()
            .args(["cat", "--no-pg"])
            .arg(&good)
            .arg(second)
            .args(["-o"])
            .arg(&output)
            .output()
            .unwrap();
        assert!(!result.status.success());
        assert_eq!(std::fs::read(&output).unwrap(), b"sentinel");
    }

    let original = std::fs::read(&good).unwrap();
    let result = binary()
        .args(["cat", "--no-pg"])
        .arg(&good)
        .args(["-o"])
        .arg(&good)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(&good).unwrap(), original);

    let header = dir.path().join("header.sam");
    std::fs::write(&header, "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n").unwrap();
    let original = std::fs::read(&header).unwrap();
    let result = binary()
        .args(["cat", "--no-pg", "--header"])
        .arg(&header)
        .arg(&good)
        .args(["-o"])
        .arg(&header)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(&header).unwrap(), original);
}

#[test]
fn reheader_preflight_failures_preserve_existing_output() {
    let dir = tempfile::tempdir().unwrap();
    let input = make_bam(
        dir.path(),
        "input",
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
    );
    let wrong = dir.path().join("wrong.sam");
    std::fs::write(
        &wrong,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n",
    )
    .unwrap();
    let output = dir.path().join("output.bam");
    std::fs::write(&output, b"sentinel").unwrap();
    let result = binary()
        .args(["reheader", "--no-pg"])
        .arg(&wrong)
        .arg(&input)
        .args(["-o"])
        .arg(&output)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(&output).unwrap(), b"sentinel");

    let truncated = dir.path().join("truncated.bam");
    std::fs::copy(&input, &truncated).unwrap();
    remove_eof(&truncated);
    let header = dir.path().join("header.sam");
    std::fs::write(&header, "@HD\tVN:1.6\n@SQ\tSN:one\tLN:1000\n").unwrap();
    let result = binary()
        .args(["reheader", "--no-pg"])
        .arg(&header)
        .arg(&truncated)
        .args(["-o"])
        .arg(&output)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(&output).unwrap(), b"sentinel");

    let original = std::fs::read(&input).unwrap();
    let result = binary()
        .args(["reheader", "--no-pg"])
        .arg(&header)
        .arg(&input)
        .args(["-o"])
        .arg(&input)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(&input).unwrap(), original);

    let original = std::fs::read(&header).unwrap();
    let result = binary()
        .args(["reheader", "--no-pg"])
        .arg(&header)
        .arg(&input)
        .args(["-o"])
        .arg(&header)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(&header).unwrap(), original);
}

#[test]
fn cat_lists_external_headers_and_stdout_are_complete() {
    let dir = tempfile::tempdir().unwrap();
    let a = make_bam(
        dir.path(),
        "a",
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\na\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
    );
    let b = make_bam(
        dir.path(),
        "b",
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nb\t0\tchr1\t20\t60\t4M\t*\t0\t0\tTGCA\tIIII\n",
    );
    let list = dir.path().join("inputs.txt");
    std::fs::write(&list, format!("{}\n\n", a.display())).unwrap();
    let header = dir.path().join("header.sam");
    std::fs::write(
        &header,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@CO\texternal\n",
    )
    .unwrap();
    let result = binary()
        .args(["cat", "--no-pg", "--list"])
        .arg(&list)
        .args(["--header"])
        .arg(&header)
        .arg(&b)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let output = dir.path().join("stdout.bam");
    std::fs::write(&output, result.stdout).unwrap();
    let sam = render(&output);
    assert!(sam.contains("@CO\texternal\n"));
    assert!(sam.find("a\t").unwrap() < sam.find("b\t").unwrap());
    assert_one_eof(&output);
}

#[test]
fn internal_eof_and_output_failures_are_fatal() {
    let dir = tempfile::tempdir().unwrap();
    let input = make_bam(
        dir.path(),
        "input",
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
    );
    let corrupt = dir.path().join("corrupt.bam");
    let mut bytes = std::fs::read(&input).unwrap();
    bytes.extend_from_slice(b"unexpected");
    bytes.extend_from_slice(&BGZF_EOF);
    std::fs::write(&corrupt, bytes).unwrap();
    let output = dir.path().join("output.bam");
    std::fs::write(&output, b"sentinel").unwrap();
    let result = binary()
        .args(["cat", "--no-pg"])
        .arg(&corrupt)
        .args(["-o"])
        .arg(&output)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(&output).unwrap(), b"sentinel");

    let cat_error = rsomics_bam::cat::write(
        std::slice::from_ref(&input),
        rsomics_bam::cat::Options::default(),
        FailingWriter { remaining: 8 },
    )
    .unwrap_err();
    assert!(matches!(cat_error, rsomics_common::RsomicsError::Io(_)));

    let header = dir.path().join("header.sam");
    std::fs::write(&header, "@HD\tVN:1.6\n@SQ\tSN:one\tLN:1000\n").unwrap();
    let reheader_error = rsomics_bam::reheader::write(
        &header,
        &input,
        rsomics_bam::reheader::Options::default(),
        FlushFail(Vec::new()),
    )
    .unwrap_err();
    assert!(matches!(
        reheader_error,
        rsomics_common::RsomicsError::Io(_)
    ));
}

#[test]
fn headers_spanning_multiple_bgzf_frames_are_supported() {
    let dir = tempfile::tempdir().unwrap();
    let mut sam = String::from("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n");
    for index in 0..2000 {
        sam.push_str(&format!("@CO\t{index:04}-{}\n", "x".repeat(48)));
    }
    sam.push_str("r1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
    let input = make_bam(dir.path(), "large", &sam);
    let output = dir.path().join("cat.bam");
    let result = binary()
        .args(["cat", "--no-pg"])
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
    let rendered = render(&output);
    assert!(rendered.contains("@CO\t1999-"));
    assert!(rendered.contains("r1\t0\tchr1\t10\t"));
}

#[test]
fn program_records_are_added_and_suppressed() {
    let dir = tempfile::tempdir().unwrap();
    let input = make_bam(
        dir.path(),
        "input",
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@PG\tID:aligner\tPN:aligner\nr1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
    );
    let cat_output = dir.path().join("cat.bam");
    let result = binary()
        .arg("cat")
        .arg(&input)
        .args(["-o"])
        .arg(&cat_output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let text = render(&cat_output);
    let program = text
        .lines()
        .find(|line| line.starts_with("@PG\tID:rsomics-bam\t"))
        .unwrap();
    assert!(program.contains("\tPN:rsomics-bam\t"));
    assert!(program.contains(concat!("\tVN:", env!("CARGO_PKG_VERSION"), "\t")));
    assert!(program.contains("\tPP:aligner"));

    let header = dir.path().join("header.sam");
    std::fs::write(
        &header,
        "@HD\tVN:1.6\n@SQ\tSN:one\tLN:1000\n@PG\tID:replacement\tPN:replacement\n",
    )
    .unwrap();
    let reheader_output = dir.path().join("reheader.bam");
    let result = binary()
        .args(["reheader", "--no-PG"])
        .arg(&header)
        .arg(&input)
        .args(["-o"])
        .arg(&reheader_output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let text = render(&reheader_output);
    assert!(text.contains("@PG\tID:replacement\tPN:replacement\n"));
    assert!(!text.contains("@PG\tID:rsomics-bam\t"));
}

#[test]
fn unsupported_stream_contracts_fail_before_replacing_output() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("input.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
    )
    .unwrap();
    let bam = make_bam(
        dir.path(),
        "input",
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
    );
    let output = dir.path().join("output.bam");

    for result in [
        binary()
            .args(["cat", "--no-PG"])
            .arg(&sam)
            .args(["-o"])
            .arg(&output)
            .output()
            .unwrap(),
        binary()
            .args(["reheader", "--no-PG"])
            .arg(&sam)
            .arg(&sam)
            .args(["-o"])
            .arg(&output)
            .output()
            .unwrap(),
        binary()
            .args(["--json", "cat", "--no-PG"])
            .arg(&bam)
            .output()
            .unwrap(),
        binary()
            .args(["--json", "reheader", "--no-PG"])
            .arg(&sam)
            .arg(&bam)
            .output()
            .unwrap(),
    ] {
        assert!(!result.status.success());
        assert!(!output.exists());
    }
}

struct FailingWriter {
    remaining: usize,
}

impl Write for FailingWriter {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        if self.remaining == 0 {
            return Err(io::Error::other("injected write failure"));
        }
        let len = buffer.len().min(self.remaining);
        self.remaining -= len;
        Ok(len)
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

struct FlushFail(Vec<u8>);

impl Write for FlushFail {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        self.0.extend_from_slice(buffer);
        Ok(buffer.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Err(io::Error::other("injected flush failure"))
    }
}
