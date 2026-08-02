use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use rsomics_bam::depth::{self, Options, PositionMode};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn write(options: Options<'_>) -> Vec<u8> {
    let mut output = Vec::new();
    let summary = depth::write(&[fixture("depth-records.sam")], options, &mut output).unwrap();
    assert_eq!(summary.inputs, 1);
    output
}

#[test]
fn default_depth_matches_committed_oracle() {
    assert_eq!(
        write(Options::default()),
        include_bytes!("golden/depth-default.tsv")
    );
}

#[test]
fn deletion_and_overlap_modes_match_committed_oracles() {
    assert_eq!(
        write(Options {
            include_deletions: true,
            ..Options::default()
        }),
        include_bytes!("golden/depth-deletions.tsv")
    );
    assert_eq!(
        write(Options {
            remove_overlaps: true,
            ..Options::default()
        }),
        include_bytes!("golden/depth-overlaps.tsv")
    );
}

#[test]
fn all_references_and_bed_selection_are_exact() {
    assert_eq!(
        write(Options {
            positions: PositionMode::AllReferences,
            ..Options::default()
        }),
        include_bytes!("golden/depth-all-references.tsv")
    );
    let bed = fixture("depth-regions.bed");
    let output = write(Options {
        bed: Some(&bed),
        ..Options::default()
    });
    assert_eq!(
        output,
        b"chr1\t2\t1\nchr1\t3\t2\nchr1\t4\t1\nchr1\t5\t0\nchr1\t10\t2\nchr1\t11\t2\nchr1\t12\t2\n"
    );
}

#[test]
fn multiple_inputs_keep_separate_depth_columns() {
    let input = fixture("depth-records.sam");
    let mut output = Vec::new();
    let summary = depth::write(&[input.clone(), input], Options::default(), &mut output).unwrap();
    assert_eq!(summary.inputs, 2);
    assert_eq!(summary.positions, 15);
    assert!(output.starts_with(b"chr1\t2\t1\t1\nchr1\t3\t2\t2\n"));
    assert!(output.ends_with(b"chr1\t16\t1\t1\n"));
}

#[test]
fn shorter_record_does_not_truncate_active_depth() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("nested.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:20\nlong\t0\tchr1\t2\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\nshort\t0\tchr1\t2\t60\t2M\t*\t0\t0\tCC\tII\n",
    )
    .unwrap();
    let mut output = Vec::new();
    depth::write(&[input], Options::default(), &mut output).unwrap();
    assert!(output.ends_with(b"chr1\t11\t1\n"));
}

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

#[test]
fn named_output_is_transactional_and_json_stays_separate() {
    let directory = tempfile::tempdir().unwrap();
    let target = directory.path().join("depth.tsv");
    fs::write(&target, b"original\n").unwrap();

    let output = binary()
        .args(["depth", "-o"])
        .arg(&target)
        .arg(fixture("depth-unsorted.sam"))
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_eq!(fs::read(&target).unwrap(), b"original\n");

    let output = run({
        let mut command = binary();
        command
            .args(["--json", "depth", "-o"])
            .arg(&target)
            .arg(fixture("depth-records.sam"));
        command
    });
    assert_eq!(
        fs::read(&target).unwrap(),
        include_bytes!("golden/depth-default.tsv")
    );
    let envelope: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(envelope["result"]["summary"]["inputs"], 1);
    assert_eq!(envelope["result"]["summary"]["positions"], 15);
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn depth_matches_samtools_1_24_for_sam_bam_and_cram() {
    let version = run({
        let mut command = Command::new("samtools");
        command.arg("--version");
        command
    });
    assert!(
        String::from_utf8(version.stdout)
            .unwrap()
            .starts_with("samtools 1.24\n")
    );

    let directory = tempfile::tempdir().unwrap();
    let sam = fixture("depth-records.sam");
    let reference = fixture("depth-reference.fa");
    let bam = directory.path().join("records.bam");
    let cram = directory.path().join("records.cram");
    run({
        let mut command = Command::new("samtools");
        command.args(["view", "-b", "-o"]).arg(&bam).arg(&sam);
        command
    });
    run({
        let mut command = Command::new("samtools");
        command
            .args(["view", "-C", "-T"])
            .arg(&reference)
            .args(["-o"])
            .arg(&cram)
            .arg(&sam);
        command
    });

    for input in [&sam, &bam, &cram] {
        for options in [
            Vec::new(),
            vec!["-J"],
            vec!["-s"],
            vec!["-q", "1"],
            vec!["-Q", "30"],
            vec!["-l", "3"],
            vec!["-g", "DUP"],
            vec!["--incl-flags", "PAIRED"],
            vec!["--require-flags", "READ1"],
            vec!["-a"],
            vec!["-a", "-a"],
            vec!["-b", fixture("depth-regions.bed").to_str().unwrap()],
            vec!["-H"],
        ] {
            let mut ours = binary();
            ours.arg("depth");
            let mut upstream = Command::new("samtools");
            upstream.arg("depth");
            if input == &cram {
                ours.args(["--reference"]).arg(&reference);
                upstream.args(["--reference"]).arg(&reference);
            }
            ours.args(&options).arg(input);
            upstream.args(&options).arg(input);
            assert_eq!(
                run(ours).stdout,
                run(upstream).stdout,
                "{input:?} {options:?}"
            );
        }
    }

    let mut ours = binary();
    ours.arg("depth").arg(&sam).arg(&sam);
    let mut upstream = Command::new("samtools");
    upstream.arg("depth").arg(&sam).arg(&sam);
    assert_eq!(run(ours).stdout, run(upstream).stdout);
}
