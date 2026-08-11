use std::io::Write;
use std::process::{Command, Stdio};

fn root() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/upstream/samtools-consensus")
}

fn command(arguments: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(arguments)
        .output()
        .unwrap()
}

#[test]
fn consensus_reads_sam_from_stdin() {
    let mut child = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args([
            "consensus",
            "--mode",
            "simple",
            "--call-fract",
            "0.6",
            "--format",
            "fastq",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&std::fs::read(root().join("consen1.sam")).unwrap())
        .unwrap();

    let output = child.wait_with_output().unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        std::fs::read(root().join("expected/1q.out")).unwrap()
    );
}

#[test]
fn consensus_failure_preserves_existing_output() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("broken.sam");
    let output = directory.path().join("consensus.fa");
    std::fs::write(&input, b"@HD\tVN:1.6\nnot-a-record\n").unwrap();
    std::fs::write(&output, b"preserve me\n").unwrap();

    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["consensus", "--output"])
        .arg(&output)
        .arg(&input)
        .output()
        .unwrap();

    assert!(!result.status.success());
    assert_eq!(std::fs::read(output).unwrap(), b"preserve me\n");
}

#[test]
fn simple_default_uses_observable_heterozygous_fraction() {
    fn run(extra: &[&str]) -> Vec<u8> {
        let mut command = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"));
        command.args([
            "consensus",
            "--mode",
            "simple",
            "--ambig",
            "--format",
            "pileup",
        ]);
        command.args(extra);
        command.arg(root().join("consen1.sam"));
        let output = command.output().unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        output.stdout
    }

    let default = run(&[]);

    assert_eq!(default, run(&["--het-fract", "0.5"]));
    assert_ne!(default, run(&["--het-fract", "0.15"]));
}

#[test]
fn uncovered_indexed_region_outputs_match_samtools_1_24() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("consen2.bam");
    let bam_text = bam.to_str().unwrap();
    let sam = root().join("consen2.sam");

    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["view", "--bam", "--no-pg", "--output"])
        .arg(&bam)
        .arg(&sam)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let output = command(&["index", "--threads", "0", bam_text]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    for (expected, format, region, all_positions) in [
        ("empty.out", "fastq", "c2:1-2", false),
        ("empty.out", "pileup", "c2:13-14", false),
        ("14q1.out", "fastq", "c2:1-2", true),
        ("14q2.out", "fastq", "c2:14-15", true),
        ("14p.out", "pileup", "c2:1-2", true),
        ("15p.out", "pileup", "c2:14-15", true),
    ] {
        let mut consensus = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"));
        consensus.args([
            "consensus",
            "--mode",
            "simple",
            "--format",
            format,
            "--region",
            region,
        ]);
        if all_positions {
            consensus.arg("-a");
        }
        let output = consensus.arg(&bam).output().unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(
            output.stdout,
            std::fs::read(root().join("expected").join(expected)).unwrap(),
            "{expected}"
        );
    }
}
