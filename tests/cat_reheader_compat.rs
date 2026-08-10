use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn samtools_available() -> bool {
    Command::new("samtools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|status| status.success())
}

fn make_bam(directory: &Path, name: &str, sam: &str) -> PathBuf {
    let source = directory.join(format!("{name}.sam"));
    let output = directory.join(format!("{name}.bam"));
    std::fs::write(&source, sam).unwrap();
    let status = Command::new("samtools")
        .args(["view", "-b", "--no-PG", "-o"])
        .arg(&output)
        .arg(&source)
        .status()
        .unwrap();
    assert!(status.success());
    output
}

fn sam_text(path: &Path) -> Vec<u8> {
    let output = Command::new("samtools")
        .args(["view", "-h", "--no-PG"])
        .arg(path)
        .output()
        .unwrap();
    assert!(output.status.success());
    output.stdout
}

fn require_success(command: &mut Command) {
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn inputs(directory: &Path) -> (PathBuf, PathBuf) {
    let a = make_bam(
        directory,
        "a",
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:a\tSM:left\na1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:a\n",
    );
    let b = make_bam(
        directory,
        "b",
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:b\tSM:right\nb1\t0\tchr1\t20\t60\t4M\t*\t0\t0\tTGCA\tIIII\tRG:Z:b\n",
    );
    (a, b)
}

#[test]
fn cat_default_and_read_group_merge_match_samtools_1_24() {
    if !samtools_available() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();
    let (a, b) = inputs(directory.path());
    let ours = directory.path().join("ours.bam");
    let oracle = directory.path().join("oracle.bam");
    require_success(
        binary()
            .args(["cat", "--no-PG"])
            .arg(&a)
            .arg(&b)
            .args(["-o"])
            .arg(&ours),
    );
    require_success(
        Command::new("samtools")
            .args(["cat", "--no-PG", "-o"])
            .arg(&oracle)
            .arg(&a)
            .arg(&b),
    );
    assert_eq!(sam_text(&ours), sam_text(&oracle));
}

#[test]
fn cat_list_and_external_header_match_samtools_1_24() {
    if !samtools_available() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();
    let (a, b) = inputs(directory.path());
    let list = directory.path().join("inputs.txt");
    std::fs::write(&list, format!("{}\n", a.display())).unwrap();
    let header = directory.path().join("header.sam");
    std::fs::write(
        &header,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@CO\texternal\n",
    )
    .unwrap();
    let ours = directory.path().join("ours.bam");
    let oracle = directory.path().join("oracle.bam");
    require_success(
        binary()
            .args(["cat", "--no-PG", "--list"])
            .arg(&list)
            .args(["--header"])
            .arg(&header)
            .arg(&b)
            .args(["-o"])
            .arg(&ours),
    );
    require_success(
        Command::new("samtools")
            .args(["cat", "--no-PG", "-b"])
            .arg(&list)
            .args(["-h"])
            .arg(&header)
            .args(["-o"])
            .arg(&oracle)
            .arg(&b),
    );
    assert_eq!(sam_text(&ours), sam_text(&oracle));
}

#[test]
fn reheader_reference_rename_matches_samtools_1_24() {
    if !samtools_available() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();
    let input = make_bam(
        directory.path(),
        "input",
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\nr1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\nr2\t0\tchr2\t20\t60\t4M\t*\t0\t0\tTGCA\tIIII\n",
    );
    let header = directory.path().join("header.sam");
    std::fs::write(
        &header,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:one\tLN:1000\n@SQ\tSN:two\tLN:2000\n@CO\treplaced\n",
    )
    .unwrap();
    let ours = directory.path().join("ours.bam");
    let oracle = directory.path().join("oracle.bam");
    require_success(
        binary()
            .args(["reheader", "--no-PG"])
            .arg(&header)
            .arg(&input)
            .args(["-o"])
            .arg(&ours),
    );
    let output = Command::new("samtools")
        .args(["reheader", "--no-PG"])
        .arg(&header)
        .arg(&input)
        .output()
        .unwrap();
    assert!(output.status.success());
    std::fs::write(&oracle, output.stdout).unwrap();
    assert_eq!(sam_text(&ours), sam_text(&oracle));
}

#[test]
fn bam_and_cram_header_sources_match_samtools_1_24() {
    if !samtools_available() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();
    let (a, b) = inputs(directory.path());
    let reference = directory.path().join("reference.fa");
    std::fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(1000))).unwrap();
    require_success(Command::new("samtools").arg("faidx").arg(&reference));

    let cat_header = directory.path().join("cat-header.sam");
    std::fs::write(
        &cat_header,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@CO\texternal-format\n",
    )
    .unwrap();
    let reheader = directory.path().join("replacement.sam");
    std::fs::write(
        &reheader,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@CO\treplacement-format\n",
    )
    .unwrap();

    for (label, format) in [("bam", "-b"), ("cram", "-C")] {
        let cat_source = directory.path().join(format!("cat-header.{label}"));
        let replacement_source = directory.path().join(format!("replacement.{label}"));
        for (source, output) in [(&cat_header, &cat_source), (&reheader, &replacement_source)] {
            let mut command = Command::new("samtools");
            command.args(["view", "--no-PG", format, "-o"]);
            command.arg(output);
            if format == "-C" {
                command.args(["-T"]).arg(&reference);
            }
            require_success(command.arg(source));
        }

        let ours_cat = directory.path().join(format!("ours-cat-{label}.bam"));
        let oracle_cat = directory.path().join(format!("oracle-cat-{label}.bam"));
        require_success(
            binary()
                .args(["cat", "--no-PG", "--header"])
                .arg(&cat_source)
                .arg(&a)
                .arg(&b)
                .args(["-o"])
                .arg(&ours_cat),
        );
        require_success(
            Command::new("samtools")
                .args(["cat", "--no-PG", "-h"])
                .arg(&cat_source)
                .args(["-o"])
                .arg(&oracle_cat)
                .arg(&a)
                .arg(&b),
        );
        assert_eq!(sam_text(&ours_cat), sam_text(&oracle_cat));

        let ours_reheader = directory.path().join(format!("ours-reheader-{label}.bam"));
        let oracle_reheader = directory
            .path()
            .join(format!("oracle-reheader-{label}.bam"));
        require_success(
            binary()
                .args(["reheader", "--no-PG"])
                .arg(&replacement_source)
                .arg(&a)
                .args(["-o"])
                .arg(&ours_reheader),
        );
        let output = Command::new("samtools")
            .args(["reheader", "--no-PG"])
            .arg(&replacement_source)
            .arg(&a)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        std::fs::write(&oracle_reheader, output.stdout).unwrap();
        assert_eq!(sam_text(&ours_reheader), sam_text(&oracle_reheader));
    }
}
