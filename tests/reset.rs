use std::io::Write;
use std::process::{Command, Stdio};

const INPUT: &str = "@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:20\n\
@RG\tID:rg1\tSM:sample1\n\
@PG\tID:aligner\tPN:aligner\n\
@CO\tdiscard me\n\
forward\t99\tchr1\t2\t60\t4M\t=\t8\t10\tACGT\tABCD\tRG:Z:rg1\tNM:i:0\tBC:Z:AA\n\
reverse\t147\tchr1\t8\t60\t4M\t=\t2\t-10\tAGTC\tABCD\tRG:Z:rg1\tNM:i:1\tBC:Z:BB\n\
secondary\t355\tchr1\t8\t60\t4M\t=\t2\t-10\tAGTC\tABCD\tRG:Z:rg1\tNM:i:1\n";

const EXPECTED: &str = "@HD\tVN:1.6\n\
@RG\tID:rg1\tSM:sample1\n\
@PG\tID:aligner\tPN:aligner\n\
forward\t77\t*\t0\t0\t*\t*\t0\t0\tACGT\tABCD\tRG:Z:rg1\tBC:Z:AA\n\
reverse\t141\t*\t0\t0\t*\t*\t0\t0\tGACT\tDCBA\tRG:Z:rg1\tBC:Z:BB\n";

const EXPECTED_KEEP_BC: &str = "@HD\tVN:1.6\n\
@RG\tID:rg1\tSM:sample1\n\
@PG\tID:aligner\tPN:aligner\n\
forward\t77\t*\t0\t0\t*\t*\t0\t0\tACGT\tABCD\tBC:Z:AA\n\
reverse\t141\t*\t0\t0\t*\t*\t0\t0\tGACT\tDCBA\tBC:Z:BB\n";

const EXPECTED_NO_RG: &str = "@HD\tVN:1.6\n\
@PG\tID:aligner\tPN:aligner\n\
forward\t77\t*\t0\t0\t*\t*\t0\t0\tACGT\tABCD\tBC:Z:AA\n\
reverse\t141\t*\t0\t0\t*\t*\t0\t0\tGACT\tDCBA\tBC:Z:BB\n";

const EXPECTED_REJECT_PG: &str = "@HD\tVN:1.6\n\
@RG\tID:rg1\tSM:sample1\n\
forward\t77\t*\t0\t0\t*\t*\t0\t0\tACGT\tABCD\tRG:Z:rg1\tBC:Z:AA\n\
reverse\t141\t*\t0\t0\t*\t*\t0\t0\tGACT\tDCBA\tRG:Z:rg1\tBC:Z:BB\n";

const DUPLICATE_INPUT: &str = "@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:20\n\
duplicate\t1123\tchr1\t2\t60\t4M\t=\t8\t10\tACGT\tABCD\n";

const AMBIGUOUS_INPUT: &str = "@HD\tVN:1.6\n\
@SQ\tSN:chr1\tLN:20\n\
ambiguous\t16\tchr1\t2\t60\t15M\t*\t0\t0\tACMGRSVTWYHKDBN\tABCDEFGHIJKLMNO\tNM:i:2\tBC:Z:odd\n";

const EXPECTED_DUPLICATE: &str = "@HD\tVN:1.6\n\
duplicate\t1101\t*\t0\t0\t*\t*\t0\t0\tACGT\tABCD\n";

#[test]
fn reset_restores_primary_reads_to_unaligned_sam() {
    let output = run_reset(&["reset", "--no-PG", "-O", "sam", "-"], INPUT);

    assert_eq!(output.stdout, EXPECTED.as_bytes());
}

#[test]
fn output_format_is_case_insensitive() {
    let output = run_reset(&["reset", "--no-PG", "-O", "SAM", "-"], INPUT);

    assert_eq!(output.stdout, EXPECTED.as_bytes());
}

#[test]
fn keep_tag_takes_precedence_over_remove_tag() {
    let output = run_reset(
        &[
            "reset",
            "--no-PG",
            "--keep-tag",
            "BC",
            "-x",
            "RG",
            "-O",
            "sam",
            "-",
        ],
        INPUT,
    );
    assert_eq!(output.stdout, EXPECTED_KEEP_BC.as_bytes());
}

#[test]
fn no_rg_drops_read_group_header_and_record_tags() {
    let output = run_reset(&["reset", "--no-PG", "--no-RG", "-O", "sam", "-"], INPUT);

    assert_eq!(output.stdout, EXPECTED_NO_RG.as_bytes());
}

#[test]
fn reject_pg_drops_the_matching_program_and_its_successors() {
    let output = run_reset(
        &[
            "reset",
            "--no-PG",
            "--reject-PG",
            "aligner",
            "-O",
            "sam",
            "-",
        ],
        INPUT,
    );

    assert_eq!(output.stdout, EXPECTED_REJECT_PG.as_bytes());
}

#[test]
fn dupflag_preserves_the_duplicate_bit() {
    let output = run_reset(
        &["reset", "--no-PG", "--dupflag", "-O", "sam", "-"],
        DUPLICATE_INPUT,
    );

    assert_eq!(output.stdout, EXPECTED_DUPLICATE.as_bytes());
}

#[test]
fn remove_tag_caret_switches_to_keep_mode() {
    let output = run_reset(&["reset", "--no-PG", "-x", "^BC", "-O", "sam", "-"], INPUT);

    assert_eq!(output.stdout, EXPECTED_KEEP_BC.as_bytes());
}

#[test]
fn named_sam_output_is_inferred_from_the_extension() {
    let directory = tempfile::tempdir().unwrap();
    let output_path = directory.path().join("reset.sam");
    let mut child = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["reset", "--no-PG", "-o"])
        .arg(&output_path)
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(INPUT.as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(std::fs::read(output_path).unwrap(), EXPECTED.as_bytes());
}

#[test]
fn named_bam_extension_is_case_insensitive() {
    let directory = tempfile::tempdir().unwrap();
    let output_path = directory.path().join("reset.BAM");
    let mut child = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["reset", "--no-PG", "-o"])
        .arg(&output_path)
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(INPUT.as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let bytes = std::fs::read(output_path).unwrap();
    assert_eq!(&bytes[..4], &[0x1f, 0x8b, 0x08, 0x04]);
}

#[test]
fn input_output_alias_is_rejected_without_replacing_the_input() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("input.sam");
    std::fs::write(&path, INPUT).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["reset", "--no-PG", "-o"])
        .arg(&path)
        .arg(&path)
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert_eq!(std::fs::read(path).unwrap(), INPUT.as_bytes());
}

#[test]
fn malformed_auxiliary_tag_list_is_rejected() {
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["reset", "--no-PG", "-x", "BC,,RG", "-"])
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("exactly two bytes"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn malformed_input_does_not_replace_an_existing_output() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("broken.sam");
    let output = directory.path().join("output.sam");
    std::fs::write(
        &input,
        "@HD\tVN:1.6\nread\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tbad\n",
    )
    .unwrap();
    std::fs::write(&output, b"previous output\n").unwrap();

    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["reset", "--no-PG", "-o"])
        .arg(&output)
        .arg(&input)
        .output()
        .unwrap();

    assert!(!result.status.success());
    assert_eq!(std::fs::read(&output).unwrap(), b"previous output\n");
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn reset_matches_samtools_1_24_for_sam_bam_and_cram() {
    let version = Command::new("samtools").arg("--version").output().unwrap();
    assert!(version.status.success());
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("samtools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    let sam = directory.path().join("input.sam");
    let bam = directory.path().join("input.bam");
    let cram = directory.path().join("input.cram");
    std::fs::write(&reference, ">chr1\nACGTACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&sam, INPUT).unwrap();
    require_success(Command::new("samtools").arg("faidx").arg(&reference));
    require_success(
        Command::new("samtools")
            .args(["view", "--no-PG", "-b", "-o"])
            .arg(&bam)
            .arg(&sam),
    );
    require_success(
        Command::new("samtools")
            .args(["view", "--no-PG", "-C", "-T"])
            .arg(&reference)
            .args(["-o"])
            .arg(&cram)
            .arg(&sam),
    );

    for input in [&sam, &bam, &cram] {
        for options in [
            vec!["--no-PG"],
            vec!["--no-PG", "--no-RG"],
            vec!["--no-PG", "--keep-tag", "BC", "-x", "RG"],
            vec!["--no-PG", "-x", "^BC"],
            vec!["--no-PG", "--reject-PG", "aligner"],
            vec!["--no-PG", "--dupflag"],
        ] {
            let mut ours = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"));
            ours.arg("reset")
                .args(&options)
                .args(["-O", "sam", "-@", "2", "-T"]);
            ours.arg(&reference).arg(input);
            let ours = require_output(&mut ours);

            let mut oracle = Command::new("samtools");
            oracle
                .arg("reset")
                .args(&options)
                .args(["-O", "sam", "-@", "2", "-T"]);
            oracle.arg(&reference).arg(input);
            let oracle = require_output(&mut oracle);

            assert_eq!(
                ours.stdout,
                oracle.stdout,
                "{} {options:?}",
                input.display()
            );
        }
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn cram_output_decodes_to_the_samtools_1_24_result() {
    let directory = tempfile::tempdir().unwrap();
    let input = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden/stats/1_map_cigar.cram");
    let ours = directory.path().join("ours.cram");
    let oracle = directory.path().join("oracle.cram");

    require_success(
        Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
            .args(["reset", "--no-PG", "-O", "cram", "-@", "2", "-o"])
            .arg(&ours)
            .arg(&input),
    );
    require_success(
        Command::new("samtools")
            .args(["reset", "--no-PG", "-O", "cram", "-@", "2", "-o"])
            .arg(&oracle)
            .arg(&input),
    );

    let ours = require_output(
        Command::new("samtools")
            .args(["view", "--no-PG", "-h"])
            .arg(&ours),
    );
    let oracle = require_output(
        Command::new("samtools")
            .args(["view", "--no-PG", "-h"])
            .arg(&oracle),
    );
    assert_eq!(ours.stdout, oracle.stdout);
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn bam_fast_path_matches_samtools_for_odd_ambiguous_reverse_reads() {
    let directory = tempfile::tempdir().unwrap();
    let sam = directory.path().join("input.sam");
    let bam = directory.path().join("input.bam");
    let ours = directory.path().join("ours.bam");
    let oracle = directory.path().join("oracle.bam");
    std::fs::write(&sam, AMBIGUOUS_INPUT).unwrap();
    require_success(
        Command::new("samtools")
            .args(["view", "--no-PG", "-b", "-o"])
            .arg(&bam)
            .arg(&sam),
    );

    require_success(
        Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
            .args(["reset", "--no-PG", "-@", "2", "-o"])
            .arg(&ours)
            .arg(&bam),
    );
    require_success(
        Command::new("samtools")
            .args(["reset", "--no-PG", "-@", "2", "-o"])
            .arg(&oracle)
            .arg(&bam),
    );

    let ours = require_output(
        Command::new("samtools")
            .args(["view", "--no-PG", "-h"])
            .arg(&ours),
    );
    let oracle = require_output(
        Command::new("samtools")
            .args(["view", "--no-PG", "-h"])
            .arg(&oracle),
    );
    assert_eq!(ours.stdout, oracle.stdout);
}

fn require_success(command: &mut Command) {
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn require_output(command: &mut Command) -> std::process::Output {
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn run_reset(arguments: &[&str], input: &str) -> std::process::Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(arguments)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}
