use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn fixture(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden/split")
        .join(path)
}

fn samtools_available() -> bool {
    Command::new("samtools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|status| status.success())
}

fn require_success(command: &mut Command) {
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn sam(path: &Path, header: bool) -> Vec<u8> {
    let mut command = Command::new("samtools");
    command.args(["view", "--no-PG"]);
    if header {
        command.arg("-h");
    }
    let output = command.arg(path).output().unwrap();
    assert!(output.status.success());
    output.stdout
}

#[test]
fn read_group_and_integer_tag_outputs_match_samtools_1_24() {
    if !samtools_available() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();
    let input = fixture("read-group/tworg.bam");
    let ours = directory.path().join("ours");
    let oracle_pattern = directory.path().join("oracle.%!.bam");
    require_success(
        binary()
            .args(["split", "--no-PG", "--output-prefix"])
            .arg(&ours)
            .arg(&input),
    );
    require_success(
        Command::new("samtools")
            .args(["split", "--no-PG", "-f"])
            .arg(&oracle_pattern)
            .arg(&input),
    );
    for label in ["rg1", "rg2"] {
        assert_eq!(
            sam(&directory.path().join(format!("ours.{label}.bam")), true),
            sam(&directory.path().join(format!("oracle.{label}.bam")), true)
        );
    }

    let source = directory.path().join("tags.sam");
    std::fs::write(
        &source,
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\na\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tF\tNM:i:0\nb\t0\tchr1\t2\t60\t1M\t*\t0\t0\tC\tF\tNM:i:6\nc\t0\tchr1\t3\t60\t1M\t*\t0\t0\tG\tF\tNM:i:4\nd\t0\tchr1\t4\t60\t1M\t*\t0\t0\tT\tF\tNM:i:3\n",
    )
    .unwrap();
    let tagged = directory.path().join("tags.bam");
    require_success(
        Command::new("samtools")
            .args(["view", "-b", "--no-PG", "-o"])
            .arg(&tagged)
            .arg(&source),
    );
    let ours = directory.path().join("tag-ours");
    let oracle_pattern = directory.path().join("tag-oracle.%!.bam");
    require_success(
        binary()
            .args(["split", "--no-PG", "--tag", "NM", "--output-prefix"])
            .arg(&ours)
            .arg(&tagged),
    );
    require_success(
        Command::new("samtools")
            .args(["split", "--no-PG", "-d", "NM", "-f"])
            .arg(&oracle_pattern)
            .arg(&tagged),
    );
    for label in ["0", "6", "4", "3"] {
        assert_eq!(
            sam(
                &directory.path().join(format!("tag-ours.{label}.bam")),
                true
            ),
            sam(
                &directory.path().join(format!("tag-oracle.{label}.bam")),
                true
            )
        );
    }
}

#[test]
#[ignore = "requires RSeQC 5.0.4 and samtools 1.24"]
fn gene_and_mate_outputs_match_rseqc_5_0_4() {
    assert!(samtools_available());
    let directory = tempfile::tempdir().unwrap();
    let rseqc = PathBuf::from(std::env::var_os("RSOMICS_RSEQC_BIN").unwrap());

    for fixture_name in ["paired", "flags"] {
        let input = fixture(&format!("mates/{fixture_name}.bam"));
        let ours = directory.path().join(format!("{fixture_name}-ours"));
        let oracle = directory.path().join(format!("{fixture_name}-oracle"));
        require_success(
            binary()
                .args(["split", "--no-PG", "--mates", "--output-prefix"])
                .arg(&ours)
                .arg(&input),
        );
        require_success(
            Command::new(rseqc.join("split_paired_bam.py"))
                .args(["-i"])
                .arg(&input)
                .args(["-o"])
                .arg(&oracle),
        );
        for label in ["R1", "R2", "unmap"] {
            assert_eq!(
                sam(
                    &directory
                        .path()
                        .join(format!("{fixture_name}-ours.{label}.bam")),
                    false
                ),
                sam(
                    &directory
                        .path()
                        .join(format!("{fixture_name}-oracle.{label}.bam")),
                    false
                )
            );
        }
    }

    let input = fixture("genes/reads.bam");
    let bed = fixture("genes/genes.strict.bed12");
    let ours = directory.path().join("genes-ours");
    let oracle = directory.path().join("genes-oracle");
    require_success(
        binary()
            .args(["split", "--no-PG", "--genes"])
            .arg(&bed)
            .args(["--output-prefix"])
            .arg(&ours)
            .arg(&input),
    );
    require_success(
        Command::new(rseqc.join("split_bam.py"))
            .args(["-i"])
            .arg(&input)
            .args(["-r"])
            .arg(&bed)
            .args(["-o"])
            .arg(&oracle),
    );
    for label in ["in", "ex", "junk"] {
        assert_eq!(
            sam(
                &directory.path().join(format!("genes-ours.{label}.bam")),
                false
            ),
            sam(
                &directory.path().join(format!("genes-oracle.{label}.bam")),
                false
            )
        );
    }
}
