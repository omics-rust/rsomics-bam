use std::fs;
#[cfg(unix)]
use std::fs::{File, OpenOptions};
#[cfg(unix)]
use std::io::Read;
use std::io::{self, Write};
#[cfg(unix)]
use std::os::fd::{AsRawFd, FromRawFd};
use std::path::{Path, PathBuf};
#[cfg(unix)]
use std::process::Stdio;
use std::process::{Command, Output};
#[cfg(unix)]
use std::thread;
#[cfg(unix)]
use std::time::{Duration, Instant};

use rsomics_bam::tview::{Format, Options};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
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

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn indexed_fixture(root: &Path) -> (PathBuf, PathBuf) {
    let bam = root.join("input.bam");
    let reference = root.join("reference.fa");
    fs::copy(golden("mpileup-reference.fa"), &reference).unwrap();
    fs::copy(
        golden("mpileup-reference.fa.fai"),
        root.join("reference.fa.fai"),
    )
    .unwrap();
    run({
        let mut command = binary();
        command
            .args(["view", "--bam", "--no-pg", "--output"])
            .arg(&bam)
            .arg(golden("mpileup-records.sam"));
        command
    });
    run({
        let mut command = binary();
        command.arg("index").arg(&bam);
        command
    });
    (bam, reference)
}

fn sampled_fixture(root: &Path) -> (PathBuf, PathBuf) {
    let sam = root.join("sampled.sam");
    let bam = root.join("sampled.bam");
    let reference = root.join("reference.fa");
    fs::copy(golden("mpileup-reference.fa"), &reference).unwrap();
    fs::copy(
        golden("mpileup-reference.fa.fai"),
        root.join("reference.fa.fai"),
    )
    .unwrap();
    fs::write(
        &sam,
        b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:120\n@RG\tID:rg1\tSM:alpha\n@RG\tID:rg2\tSM:beta\na\t0\tchr1\t5\t60\t12M\t*\t0\t0\tACGTACGTACGT\tIIIIIIIIIIII\tRG:Z:rg1\nb\t0\tchr1\t5\t60\t12M\t*\t0\t0\tTTTTTTTTTTTT\tIIIIIIIIIIII\tRG:Z:rg2\n",
    )
    .unwrap();
    run({
        let mut command = binary();
        command
            .args(["view", "--bam", "--no-pg", "--output"])
            .arg(&bam)
            .arg(&sam);
        command
    });
    run({
        let mut command = binary();
        command.arg("index").arg(&bam);
        command
    });
    (bam, reference)
}

fn edge_fixture(root: &Path) -> (PathBuf, PathBuf) {
    let sam = root.join("edge.sam");
    let bam = root.join("edge.bam");
    let reference = root.join("reference.fa");
    fs::copy(golden("mpileup-reference.fa"), &reference).unwrap();
    fs::copy(
        golden("mpileup-reference.fa.fai"),
        root.join("reference.fa.fai"),
    )
    .unwrap();
    fs::write(
        &sam,
        b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:120\nsoft-ins-del\t0\tchr1\t2\t60\t2S4M1I3M1D2M2S\t*\t0\t0\tTTACGTAGTACGAA\tIIIIIIIIIIIIII\nskip\t16\tchr1\t3\t25\t2M3N3M\t*\t0\t0\tACGTA\tIIIII\nskip-forward\t0\tchr1\t3\t25\t2M3N3M\t*\t0\t0\tTGCAT\tIIIII\npad-eq-x\t256\tchr1\t5\t15\t2=1P2X\t*\t0\t0\t==TT\tIIII\nmissing\t1\tchr1\t6\t255\t5M\t*\t0\t0\tACGTA\t*\nhard\t0\tchr1\t8\t40\t2H5M2H\t*\t0\t0\tTGCAT\tIIIII\ncolor\t0\tchr1\t10\t35\t5M\t*\t0\t0\tCGTAC\tIIIII\tCS:Z:A01230\tCQ:Z:IIIII\n",
    )
    .unwrap();
    run({
        let mut command = binary();
        command
            .args(["view", "--bam", "--no-pg", "--output"])
            .arg(&bam)
            .arg(&sam);
        command
    });
    run({
        let mut command = binary();
        command.arg("index").arg(&bam);
        command
    });
    (bam, reference)
}

fn dense_fixture(root: &Path) -> (PathBuf, PathBuf) {
    let sam = root.join("dense.sam");
    let bam = root.join("dense.bam");
    let reference = root.join("reference.fa");
    fs::copy(golden("mpileup-reference.fa"), &reference).unwrap();
    fs::copy(
        golden("mpileup-reference.fa.fai"),
        root.join("reference.fa.fai"),
    )
    .unwrap();
    let mut writer = io::BufWriter::new(File::create(&sam).unwrap());
    writer
        .write_all(b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:120\n")
        .unwrap();
    for index in 0..8001 {
        writeln!(
            writer,
            "dense{index:05}\t0\tchr1\t5\t60\t12M\t*\t0\t0\tACGTACGTACGT\tIIIIIIIIIIII"
        )
        .unwrap();
    }
    writer.flush().unwrap();
    run({
        let mut command = binary();
        command
            .args(["view", "--bam", "--no-pg", "--output"])
            .arg(&bam)
            .arg(&sam);
        command
    });
    run({
        let mut command = binary();
        command.arg("index").arg(&bam);
        command
    });
    (bam, reference)
}

fn padded(lines: &[&str], width: usize) -> Vec<u8> {
    let mut output = Vec::new();
    for line in lines {
        output.extend_from_slice(format!("{line:<width$}\n").as_bytes());
    }
    output
}

fn expected_with_insertions() -> Vec<u8> {
    padded(
        &[
            "1         11          21",
            "ACGTACGTACGTACGTA**CGTACGTACGTACGTACGTAC",
            "    .....MSKWMSKW  MSKWMSKYG..CGTACGTACG",
            "    .............**.......     gtacgtacg",
            "         ACGTACGTTTACGTACGT",
            "                   gtacgtacgtacgtacgtac",
        ],
        40,
    )
}

fn expected_without_insertions() -> Vec<u8> {
    padded(
        &[
            "1         11        21        31",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            "    .....MSKWMSKWMSKWMSKYG..CGTACGTACGT.",
            "    ....................     gtacgtacgt",
            "         ACGTACGTACGTACGT              A",
            "                 gtacgtacgtacgtacgtac",
        ],
        40,
    )
}

#[test]
fn public_text_renderer_matches_the_samtools_grid() {
    let directory = tempfile::tempdir().unwrap();
    let (bam, reference) = indexed_fixture(directory.path());
    let mut output = Vec::new();
    let summary = rsomics_bam::tview::write(
        &bam,
        Options {
            reference: Some(&reference),
            position: Some("chr1:1"),
            width: 40,
            ..Options::default()
        },
        Format::Text,
        &mut output,
    )
    .unwrap();

    assert_eq!(output, expected_with_insertions());
    assert_eq!(summary.reference, "chr1");
    assert_eq!(summary.start, 1);
    assert_eq!(summary.width, 40);
    assert_eq!(summary.alignment_rows, 3);
}

#[test]
fn command_text_modes_match_the_stable_grids() {
    let directory = tempfile::tempdir().unwrap();
    let (bam, reference) = indexed_fixture(directory.path());
    let default = run({
        let mut command = binary();
        command
            .args([
                "tview",
                "--display",
                "text",
                "--width",
                "40",
                "--position",
                "chr1:1",
                "--reference",
            ])
            .arg(&reference)
            .arg(&bam);
        command
    });
    assert_eq!(default.stdout, expected_with_insertions());

    let hidden = run({
        let mut command = binary();
        command
            .args(["tview", "-d", "T", "-w", "40", "-p", "chr1:1", "-i", "-T"])
            .arg(&reference)
            .arg(&bam);
        command
    });
    assert_eq!(hidden.stdout, expected_without_insertions());

    let threaded = run({
        let mut command = binary();
        command
            .args([
                "tview", "-d", "text", "-w", "40", "-p", "chr1:1", "-@", "2", "-T",
            ])
            .arg(&reference)
            .arg(&bam);
        command
    });
    assert_eq!(threaded.stdout, expected_with_insertions());
}

#[test]
fn columns_environment_controls_only_the_noninteractive_default() {
    let directory = tempfile::tempdir().unwrap();
    let (bam, _) = indexed_fixture(directory.path());
    let output = run({
        let mut command = binary();
        command
            .args(["tview", "-d", "T"])
            .arg(&bam)
            .env("COLUMNS", "23");
        command
    });
    assert!(
        output
            .stdout
            .split(|byte| *byte == b'\n')
            .all(|line| { line.is_empty() || line.len() == 23 })
    );

    let output = run({
        let mut command = binary();
        command
            .args(["tview", "-d", "T"])
            .arg(bam)
            .env("COLUMNS", "invalid");
        command
    });
    assert!(
        output
            .stdout
            .split(|byte| *byte == b'\n')
            .all(|line| { line.is_empty() || line.len() == 80 })
    );
}

#[test]
fn named_output_json_and_failures_are_transactional() {
    let directory = tempfile::tempdir().unwrap();
    let (bam, reference) = indexed_fixture(directory.path());
    let target = directory.path().join("view.txt");
    let output = run({
        let mut command = binary();
        command
            .args([
                "--json", "tview", "-d", "text", "-w", "40", "-p", "chr1:1", "-T",
            ])
            .arg(&reference)
            .args(["--output"])
            .arg(&target)
            .arg(&bam);
        command
    });
    assert_eq!(fs::read(&target).unwrap(), expected_with_insertions());
    let json: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(json["result"]["command"], "tview");
    assert_eq!(json["result"]["summary"]["reference"], "chr1");

    fs::write(&target, b"sentinel").unwrap();
    let missing = directory.path().join("missing.bam");
    let failed = binary()
        .args(["tview", "-d", "text", "--output"])
        .arg(&target)
        .arg(missing)
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert_eq!(fs::read(target).unwrap(), b"sentinel");
}

#[test]
fn public_renderer_propagates_write_failures() {
    struct FailedWriter;

    impl Write for FailedWriter {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::new(io::ErrorKind::BrokenPipe, "closed"))
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    let directory = tempfile::tempdir().unwrap();
    let (bam, reference) = indexed_fixture(directory.path());
    let error = rsomics_bam::tview::write(
        &bam,
        Options {
            reference: Some(&reference),
            ..Options::default()
        },
        Format::Text,
        FailedWriter,
    )
    .unwrap_err();
    assert!(error.to_string().contains("closed"));
}

#[test]
fn custom_index_and_html_escape_the_location() {
    let directory = tempfile::tempdir().unwrap();
    let (bam, reference) = indexed_fixture(directory.path());
    let custom_index = directory.path().join("custom.bai");
    fs::rename(bam.with_extension("bam.bai"), &custom_index).unwrap();
    let output = run({
        let mut command = binary();
        command
            .args([
                "tview",
                "--display",
                "html",
                "--width",
                "20",
                "--position",
                "chr1:2",
                "--index",
            ])
            .arg(custom_index)
            .args(["--reference"])
            .arg(reference)
            .arg(bam);
        command
    });
    let html = String::from_utf8(output.stdout).unwrap();
    assert!(html.starts_with("<!doctype html>"));
    assert!(html.contains("<pre"));
    assert!(html.contains("chr1:2"));
    assert!(!html.contains("\u{1b}["));
}

#[test]
fn html_preserves_cell_styles_and_escapes_alignment_symbols() {
    let directory = tempfile::tempdir().unwrap();
    let (bam, reference) = edge_fixture(directory.path());
    let output = run({
        let mut command = binary();
        command
            .args(["tview", "-d", "H", "-w", "40", "-p", "chr1:1", "-T"])
            .arg(reference)
            .arg(bam);
        command
    });
    let html = String::from_utf8(output.stdout).unwrap();
    assert!(html.contains("data-location=\"chr1:1\""));
    assert_eq!(html.matches("<span class=\"location\">").count(), 1);
    assert!(html.contains("class=\"green underline\""));
    assert!(html.contains("&lt;"));
    assert!(html.contains("&gt;"));
    assert!(!html.contains("\u{1b}["));
}

#[test]
fn malformed_alignment_and_index_inputs_fail_nonzero() {
    let directory = tempfile::tempdir().unwrap();
    let (bam, _) = indexed_fixture(directory.path());
    let truncated = directory.path().join("truncated.bam");
    let mut bytes = fs::read(&bam).unwrap();
    bytes.truncate(bytes.len() / 2);
    fs::write(&truncated, bytes).unwrap();
    fs::copy(
        bam.with_extension("bam.bai"),
        truncated.with_extension("bam.bai"),
    )
    .unwrap();
    let output = binary()
        .args(["tview", "-d", "T"])
        .arg(truncated)
        .output()
        .unwrap();
    assert!(!output.status.success());

    let wrong_index = directory.path().join("wrong.bai");
    fs::write(&wrong_index, b"not an index").unwrap();
    let output = binary()
        .args(["tview", "-d", "T", "-X"])
        .arg(wrong_index)
        .arg(bam)
        .output()
        .unwrap();
    assert!(!output.status.success());
}

#[test]
fn threaded_cram_matches_the_single_thread_viewport() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.cram");
    fs::copy(golden("cram-size/version-3.1.cram"), &input).unwrap();
    run({
        let mut command = binary();
        command.arg("index").arg(&input);
        command
    });
    let single = run({
        let mut command = binary();
        command.args(["tview", "-d", "T", "-w", "40"]).arg(&input);
        command
    });
    let threaded = run({
        let mut command = binary();
        command
            .args(["tview", "-d", "T", "-w", "40", "-@", "2"])
            .arg(input);
        command
    });
    assert_eq!(threaded.stdout, single.stdout);
}

#[test]
fn sample_and_read_group_selectors_form_a_checked_union() {
    let directory = tempfile::tempdir().unwrap();
    let (bam, reference) = sampled_fixture(directory.path());
    for (selector, expected_rows) in [("alpha", 1), ("rg2", 1), ("alpha,rg2", 2)] {
        let mut output = Vec::new();
        let summary = rsomics_bam::tview::write(
            &bam,
            Options {
                reference: Some(&reference),
                position: Some("chr1:1"),
                sample: Some(selector),
                width: 30,
                ..Options::default()
            },
            Format::Text,
            &mut output,
        )
        .unwrap();
        assert_eq!(summary.alignment_rows, expected_rows, "{selector}");
        assert_eq!(output.len(), (3 + expected_rows) * 31, "{selector}");
    }
    let error = rsomics_bam::tview::write(
        &bam,
        Options {
            sample: Some("missing"),
            ..Options::default()
        },
        Format::Text,
        Vec::new(),
    )
    .unwrap_err();
    assert!(error.to_string().contains("absent"));
}

#[test]
fn complete_cigar_and_missing_quality_cells_render_without_loss() {
    let directory = tempfile::tempdir().unwrap();
    let (bam, reference) = edge_fixture(directory.path());
    let output = run({
        let mut command = binary();
        command
            .args(["tview", "-d", "T", "-w", "40", "-p", "chr1:1", "-T"])
            .arg(reference)
            .arg(bam);
        command
    });
    assert!(output.stdout.contains(&b'>'));
    assert!(output.stdout.contains(&b'<'));
    assert!(output.stdout.contains(&b'*'));
    assert!(output.stdout.contains(&b'='));
}

#[test]
fn interactive_mode_requires_a_terminal_and_json_never_mixes_streams() {
    let directory = tempfile::tempdir().unwrap();
    let (bam, _) = indexed_fixture(directory.path());
    let interactive = binary().arg("tview").arg(&bam).output().unwrap();
    assert!(!interactive.status.success());
    assert!(String::from_utf8_lossy(&interactive.stderr).contains("terminal"));

    let width = binary()
        .args(["tview", "--width", "20"])
        .arg(&bam)
        .output()
        .unwrap();
    assert!(!width.status.success());
    assert!(String::from_utf8_lossy(&width.stderr).contains("text and HTML"));

    let json = binary()
        .args(["--json", "tview", "--display", "text"])
        .arg(bam)
        .output()
        .unwrap();
    assert!(!json.status.success());
    assert!(String::from_utf8_lossy(&json.stderr).contains("--output"));
}

#[cfg(unix)]
#[test]
fn interactive_terminal_draws_handles_events_and_restores_the_pty() {
    let directory = tempfile::tempdir().unwrap();
    let (bam, reference) = indexed_fixture(directory.path());
    let (mut master, slave, original) = pty(40, 10);
    let mut command = binary();
    command
        .args(["tview", "-d", "C", "-p", "chr1:1", "-T"])
        .arg(reference)
        .arg(bam)
        .env("TERM", "xterm-256color")
        .stdin(Stdio::from(slave.try_clone().unwrap()))
        .stdout(Stdio::from(slave.try_clone().unwrap()))
        .stderr(Stdio::from(slave.try_clone().unwrap()));
    let mut child = command.spawn().unwrap();
    thread::sleep(Duration::from_millis(300));
    master.write_all(b"?xgchr1:2\r").unwrap();
    thread::sleep(Duration::from_millis(150));
    resize(&master, 52, 12);
    master.write_all(b".imnbczNCsrvq").unwrap();

    let (status, transcript) = finish_pty(&mut child, &mut master);
    assert!(status.success(), "{}", String::from_utf8_lossy(&transcript));
    let transcript = String::from_utf8_lossy(&transcript);
    assert!(transcript.contains("\u{1b}[?1049h"), "{transcript:?}");
    assert!(transcript.contains("\u{1b}[?1049l"), "{transcript:?}");
    assert!(transcript.contains("\u{1b}[?25l"), "{transcript:?}");
    assert!(transcript.contains("\u{1b}[?25h"), "{transcript:?}");
    assert!(transcript.contains("rsomics-bam tview"), "{transcript:?}");
    assert!(
        transcript.contains("goto&gt;") || transcript.contains("goto>"),
        "{transcript:?}"
    );
    assert_eq!(terminal_flags(&slave), original);
}

#[cfg(unix)]
#[test]
fn interactive_application_failure_restores_the_pty() {
    let directory = tempfile::tempdir().unwrap();
    let (bam, _) = indexed_fixture(directory.path());
    let (mut master, slave, original) = pty(40, 10);
    let mut command = binary();
    command
        .arg("tview")
        .arg(&bam)
        .env("TERM", "xterm-256color")
        .stdin(Stdio::from(slave.try_clone().unwrap()))
        .stdout(Stdio::from(slave.try_clone().unwrap()))
        .stderr(Stdio::from(slave.try_clone().unwrap()));
    let mut child = command.spawn().unwrap();
    thread::sleep(Duration::from_millis(300));
    fs::remove_file(bam.with_extension("bam.bai")).unwrap();
    master.write_all(b"l").unwrap();

    let (status, transcript) = finish_pty(&mut child, &mut master);
    assert!(
        !status.success(),
        "{}",
        String::from_utf8_lossy(&transcript)
    );
    let transcript = String::from_utf8_lossy(&transcript);
    assert!(transcript.contains("\u{1b}[?1049l"), "{transcript:?}");
    assert!(transcript.contains("\u{1b}[?25h"), "{transcript:?}");
    assert_eq!(terminal_flags(&slave), original);
}

#[cfg(unix)]
fn finish_pty(
    child: &mut std::process::Child,
    master: &mut File,
) -> (std::process::ExitStatus, Vec<u8>) {
    let deadline = Instant::now() + Duration::from_secs(15);
    let mut transcript = Vec::new();
    loop {
        let mut buffer = [0; 16 * 1024];
        match master.read(&mut buffer) {
            Ok(0) => {}
            Ok(length) => transcript.extend_from_slice(&buffer[..length]),
            Err(error) if error.kind() == std::io::ErrorKind::WouldBlock => {}
            Err(error) => panic!("reading PTY: {error}"),
        }
        if let Some(status) = child.try_wait().unwrap() {
            return (status, transcript);
        }
        if Instant::now() >= deadline {
            child.kill().unwrap();
            panic!(
                "terminal tview did not exit: {}",
                String::from_utf8_lossy(&transcript)
            );
        }
        thread::sleep(Duration::from_millis(20));
    }
}

#[cfg(unix)]
fn pty(width: u16, height: u16) -> (File, File, libc::tcflag_t) {
    unsafe {
        let master_fd = libc::posix_openpt(libc::O_RDWR | libc::O_NOCTTY);
        assert!(master_fd >= 0);
        assert_eq!(libc::grantpt(master_fd), 0);
        assert_eq!(libc::unlockpt(master_fd), 0);
        let name = libc::ptsname(master_fd);
        assert!(!name.is_null());
        let path = std::ffi::CStr::from_ptr(name).to_str().unwrap();
        let slave = OpenOptions::new()
            .read(true)
            .write(true)
            .open(path)
            .unwrap();
        let master = File::from_raw_fd(master_fd);
        resize(&master, width, height);
        let flags = libc::fcntl(master.as_raw_fd(), libc::F_GETFL);
        assert!(flags >= 0);
        assert_eq!(
            libc::fcntl(master.as_raw_fd(), libc::F_SETFL, flags | libc::O_NONBLOCK),
            0
        );
        let original = terminal_flags(&slave);
        (master, slave, original)
    }
}

#[cfg(unix)]
fn resize(master: &File, width: u16, height: u16) {
    let size = libc::winsize {
        ws_row: height,
        ws_col: width,
        ws_xpixel: 0,
        ws_ypixel: 0,
    };
    unsafe {
        assert_eq!(libc::ioctl(master.as_raw_fd(), libc::TIOCSWINSZ, &size), 0);
    }
}

#[cfg(unix)]
fn terminal_flags(slave: &File) -> libc::tcflag_t {
    unsafe {
        let mut attributes = std::mem::zeroed();
        assert_eq!(libc::tcgetattr(slave.as_raw_fd(), &mut attributes), 0);
        attributes.c_lflag
    }
}

#[test]
fn invalid_display_width_position_and_index_fail_loudly() {
    let directory = tempfile::tempdir().unwrap();
    let (bam, _) = indexed_fixture(directory.path());
    for arguments in [
        vec!["tview", "--display", "unknown"],
        vec!["tview", "--display", "text", "--width", "0"],
        vec!["tview", "--display", "text", "--width", "1000001"],
        vec!["tview", "--display", "text", "--position", "chr1:0"],
    ] {
        let output = binary().args(arguments).arg(&bam).output().unwrap();
        assert!(!output.status.success());
    }

    let output = binary()
        .args(["tview", "--display", "text", "--index", "missing.bai"])
        .arg(bam)
        .output()
        .unwrap();
    assert!(!output.status.success());
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn text_oracle_matches_samtools_1_24() {
    let version = run({
        let mut command = Command::new("samtools");
        command.arg("--version");
        command
    });
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("samtools 1.24\n"));
    let directory = tempfile::tempdir().unwrap();
    let (bam, reference) = indexed_fixture(directory.path());

    for hide_insertions in [false, true] {
        let mut ours = binary();
        ours.args(["tview", "-d", "T", "-w", "40", "-p", "chr1:1", "-T"])
            .arg(&reference);
        let mut oracle = Command::new("samtools");
        oracle.args(["tview", "-d", "T", "-w", "40", "-p", "chr1:1"]);
        if hide_insertions {
            ours.arg("-i");
            oracle.arg("-i");
        }
        let ours = run({
            ours.arg(&bam);
            ours
        });
        let oracle = run({
            oracle.arg(&bam).arg(&reference);
            oracle
        });
        assert_eq!(ours.stdout, oracle.stdout);
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn sample_and_read_group_views_match_samtools_1_24() {
    let directory = tempfile::tempdir().unwrap();
    let (bam, reference) = sampled_fixture(directory.path());
    for selector in ["alpha", "rg2"] {
        let ours = run({
            let mut command = binary();
            command
                .args([
                    "tview", "-d", "T", "-w", "40", "-p", "chr1:1", "-s", selector, "-T",
                ])
                .arg(&reference)
                .arg(&bam);
            command
        });
        let oracle = run({
            let mut command = Command::new("samtools");
            command
                .args([
                    "tview", "-d", "T", "-w", "40", "-p", "chr1:1", "-s", selector,
                ])
                .arg(&bam)
                .arg(&reference);
            command
        });
        assert_eq!(ours.stdout, oracle.stdout, "{selector}");
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn depth_limit_matches_samtools_1_24() {
    let directory = tempfile::tempdir().unwrap();
    let (bam, reference) = dense_fixture(directory.path());
    let ours = run({
        let mut command = binary();
        command
            .args([
                "tview", "-@", "2", "-d", "T", "-w", "40", "-p", "chr1:1", "-T",
            ])
            .arg(&reference)
            .arg(&bam);
        command
    });
    let oracle = run({
        let mut command = Command::new("samtools");
        command
            .args(["tview", "-d", "T", "-w", "40", "-p", "chr1:1"])
            .arg(&bam)
            .arg(&reference);
        command
    });
    assert_eq!(ours.stdout, oracle.stdout);
    assert_eq!(
        ours.stdout.iter().filter(|byte| **byte == b'\n').count(),
        8003
    );
}

#[test]
#[ignore = "release oracle: requires the samtools 1.24 source tree"]
fn large_coordinate_fixture_matches_samtools_1_24() {
    let samtools = executable("samtools");
    let root = samtools.parent().unwrap();
    let fixture = root.join("test/large_pos/longref.sam");
    let expected = root.join("test/large_pos/tview.expected.out");
    assert!(fixture.is_file(), "missing {}", fixture.display());
    assert!(expected.is_file(), "missing {}", expected.display());
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("longref.sam.gz");
    run({
        let mut command = Command::new(&samtools);
        command
            .args(["view", "-h", "--no-PG", "-O", "sam.gz", "-o"])
            .arg(&input)
            .arg(&fixture);
        command
    });
    run({
        let mut command = Command::new(&samtools);
        command.args(["index", "-c"]).arg(&input);
        command
    });
    let ours = run({
        let mut command = binary();
        command
            .args([
                "tview",
                "-d",
                "T",
                "-w",
                "80",
                "-@",
                "2",
                "-p",
                "CHROMOSOME_I:10000000000",
            ])
            .arg(&input);
        command
    });
    let oracle = run({
        let mut command = Command::new(samtools);
        command
            .args([
                "tview",
                "-d",
                "T",
                "-w",
                "80",
                "-p",
                "CHROMOSOME_I:10000000000",
            ])
            .arg(input);
        command
    });
    let expected = fs::read(expected).unwrap();
    assert_eq!(ours.stdout, expected);
    assert_eq!(ours.stdout, oracle.stdout);
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn complete_cigar_text_matrix_matches_samtools_1_24() {
    let version = run({
        let mut command = Command::new("samtools");
        command.arg("--version");
        command
    });
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("samtools 1.24\n"));
    let directory = tempfile::tempdir().unwrap();
    let (bam, reference) = edge_fixture(directory.path());
    for (position, reference_path, hide_insertions) in [
        ("chr1:1", Some(reference.as_path()), false),
        ("chr1:1", Some(reference.as_path()), true),
        ("chr1:4", None, false),
        ("chr1:118", Some(reference.as_path()), false),
    ] {
        let mut ours = binary();
        ours.args(["tview", "-d", "T", "-w", "40", "-p", position]);
        let mut oracle = Command::new("samtools");
        oracle.args(["tview", "-d", "T", "-w", "40", "-p", position]);
        if hide_insertions {
            ours.arg("-i");
            oracle.arg("-i");
        }
        if let Some(reference_path) = reference_path {
            ours.arg("-T").arg(reference_path);
        }
        let ours = run({
            ours.arg(&bam);
            ours
        });
        let oracle = run({
            oracle.arg(&bam);
            if let Some(reference_path) = reference_path {
                oracle.arg(reference_path);
            }
            oracle
        });
        assert_eq!(ours.stdout, oracle.stdout, "{position}");
    }
}

fn executable(name: &str) -> PathBuf {
    std::env::split_paths(&std::env::var_os("PATH").unwrap())
        .map(|path| path.join(name))
        .find(|path| path.is_file())
        .unwrap_or_else(|| panic!("{name} is absent from PATH"))
}
