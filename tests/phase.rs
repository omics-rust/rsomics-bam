use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use noodles::bam;
use rsomics_bamio::raw::RecordReader;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn run(extra: &[&str], input: &str) -> std::process::Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"));
    command.arg("phase").args(extra).arg(fixture(input));
    command.output().unwrap()
}

fn boundary_fixture(directory: &Path) -> PathBuf {
    let sequences = [
        "ACGTACGTAA",
        "ACGTACGTAA",
        "TCGTACGTAA",
        "TCGTACGTAA",
        "ACGTTCGTAA",
        "ACGTTCGTAA",
        "TCGTTCGTAA",
        "TCGTTCGTAA",
        "ACGTACGTAA",
        "TCGTTCGTAA",
        "ACGTACGTAA",
        "TCGTTCGTAA",
    ];
    let mut source =
        String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n");
    for (group, reference, position) in [("a", "chr1", 100), ("b", "chr1", 300), ("c", "chr2", 50)]
    {
        for (index, sequence) in sequences.iter().enumerate() {
            let sequence = if group == "b" {
                if sequence.starts_with('T') {
                    "TCGTACGTAA"
                } else {
                    "ACGTACGTAA"
                }
            } else {
                sequence
            };
            source.push_str(&format!(
                "{group}{}\t0\t{reference}\t{position}\t60\t10M\t*\t0\t0\t{sequence}\tIIIIIIIIII\n",
                index + 1
            ));
        }
    }
    let path = directory.join("boundaries.sam");
    fs::write(&path, source).unwrap();
    path
}

fn chimera_fixture(directory: &Path) -> PathBuf {
    let mut source = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n");
    for index in 1..=6 {
        source.push_str(&format!(
            "a{index}\t0\tchr1\t100\t60\t12M\t*\t0\t0\tAAAAAAAAAAAA\tIIIIIIIIIIII\n"
        ));
        source.push_str(&format!(
            "t{index}\t0\tchr1\t100\t60\t12M\t*\t0\t0\tTTTTTTTTTTTT\tIIIIIIIIIIII\n"
        ));
    }
    source.push_str("head\t0\tchr1\t100\t60\t12M\t*\t0\t0\tAAAAAAATTTTT\tIIIIIIIIIIII\n");
    source.push_str("tail\t0\tchr1\t100\t60\t12M\t*\t0\t0\tTTTTTAAAAAAA\tIIIIIIIIIIII\n");
    let path = directory.join("chimeras.sam");
    fs::write(&path, source).unwrap();
    path
}

fn normalize_evidence(output: &[u8]) -> Vec<u8> {
    let text = String::from_utf8(output.to_vec()).unwrap();
    let mut normalized = String::new();
    let mut evidence = Vec::new();
    let flush = |normalized: &mut String, evidence: &mut Vec<&str>| {
        evidence.sort_unstable();
        for line in evidence.drain(..) {
            normalized.push_str(line);
            normalized.push('\n');
        }
    };
    for line in text.lines() {
        if line.starts_with("EV\t") {
            evidence.push(line);
        } else {
            flush(&mut normalized, &mut evidence);
            normalized.push_str(line);
            normalized.push('\n');
        }
    }
    flush(&mut normalized, &mut evidence);
    normalized.into_bytes()
}

#[test]
fn default_phase_records_match_samtools_1_24() {
    let output = run(&[], "phase.sam");

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        normalize_evidence(&output.stdout),
        normalize_evidence(include_bytes!("golden/phase-default.txt"))
    );
}

#[test]
fn default_het_lod_is_the_observable_37() {
    let implicit = run(&[], "phase-lod.sam");
    let q37 = run(&["-q", "37"], "phase-lod.sam");
    let q40 = run(&["-q", "40"], "phase-lod.sam");

    for output in [&implicit, &q37, &q40] {
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    assert_eq!(implicit.stdout, q37.stdout);
    assert!(implicit.stdout.windows(3).any(|window| window == b"M0\t"));
    assert!(!q40.stdout.windows(3).any(|window| window == b"M0\t"));
}

#[test]
fn one_marker_window_preserves_samtools_haplotype_orientation() {
    let output = run(&["-k", "1"], "phase.sam");

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report = String::from_utf8(output.stdout).unwrap();
    assert!(
        report.contains("M2\tchr1\t100\t100\tA\tT\t1\t2\t0\t6\t4\n"),
        "{report}"
    );
    assert!(
        report.contains("M1\tchr1\t100\t104\tT\tA\t2\t2\t0\t6\t4\n"),
        "{report}"
    );
}

#[test]
fn window_is_rejected_before_an_unbounded_workspace_is_possible() {
    let accepted = run(&["-k", "23", "-q", "1000"], "phase.sam");
    let rejected = run(&["-k", "24"], "phase.sam");
    let excessive_depth = run(&["-D", "65536"], "phase.sam");

    assert!(
        accepted.status.success(),
        "{}",
        String::from_utf8_lossy(&accepted.stderr)
    );
    assert!(!rejected.status.success());
    assert!(
        String::from_utf8_lossy(&rejected.stderr).contains("phase window must be between 1 and 23")
    );
    assert!(!excessive_depth.status.success());
    assert!(
        String::from_utf8_lossy(&excessive_depth.stderr)
            .contains("maximum phase depth must be between 1 and 65535")
    );
}

#[test]
fn marker_indexes_are_reference_local_across_phase_sets() {
    let directory = tempfile::tempdir().unwrap();
    let input = boundary_fixture(directory.path());
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .arg("phase")
        .arg(input)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report = String::from_utf8(output.stdout).unwrap();
    assert!(report.contains("M0\tchr1\t300\t300\tT\tA\t3\t0\t0\t0\t0\n"));
    assert!(report.contains("M1\tchr2\t50\t50\tT\tA\t1\t4\t0\t6\t2\n"));
    assert!(report.contains("M1\tchr2\t50\t54\tT\tA\t2\t4\t0\t6\t2\n"));
    assert!(report.contains("EV\t0\tchr2\t1\t40\t2M\t"));
}

fn bam_records(path: &Path) -> (usize, Vec<(String, bool)>) {
    let mut reader = bam::io::Reader::new(std::fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    let mut records = RecordReader::new(reader.get_mut());
    let mut values = Vec::new();
    while let Some(record) = records.next().unwrap() {
        values.push((
            String::from_utf8(record.name().to_vec()).unwrap(),
            record.aux_type(*b"ZP") == Some(b'A') && record.aux_value(*b"ZP") == Some(b"Y"),
        ));
    }
    (header.programs().as_ref().len(), values)
}

#[test]
fn bam_partitions_preserve_records_tags_and_program_suppression() {
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("haplotypes");
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["phase", "--no-pg", "-b"])
        .arg(&prefix)
        .arg(fixture("phase.sam"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let mut all = Vec::new();
    for (suffix, expected) in [
        ("0.bam", &["r1", "r2", "r3", "r9", "r11"][..]),
        ("1.bam", &["r4", "r5", "r6", "r7", "r8", "r10", "r12"][..]),
        ("chimera.bam", &[][..]),
    ] {
        let (programs, records) =
            bam_records(&PathBuf::from(format!("{}.{}", prefix.display(), suffix)));
        assert_eq!(programs, 0);
        assert_eq!(
            records
                .iter()
                .map(|(name, _)| name.as_str())
                .collect::<Vec<_>>(),
            expected
        );
        all.extend(records);
    }
    all.sort_unstable();
    let mut expected: Vec<_> = (1..=12).map(|index| format!("r{index}")).collect();
    expected.sort_unstable();
    assert_eq!(all.len(), 12);
    assert_eq!(all.iter().filter(|(_, tagged)| *tagged).count(), 8);
    assert_eq!(
        all.iter().map(|(name, _)| name).collect::<Vec<_>>(),
        expected.iter().collect::<Vec<_>>()
    );
}

#[test]
fn ambiguous_reads_are_routed_to_the_chimera_partition() {
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("ambiguous");
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["phase", "--no-pg", "-A", "-b"])
        .arg(&prefix)
        .arg(fixture("phase.sam"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let (_, records) = bam_records(&PathBuf::from(format!("{}.chimera.bam", prefix.display())));
    assert_eq!(
        records
            .iter()
            .map(|(name, tagged)| (name.as_str(), *tagged))
            .collect::<Vec<_>>(),
        [("r3", false), ("r4", false), ("r5", false), ("r6", false)]
    );
}

#[test]
fn head_and_tail_chimeras_are_repaired_and_partitioned() {
    let directory = tempfile::tempdir().unwrap();
    let input = chimera_fixture(directory.path());
    let repaired = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .arg("phase")
        .arg(&input)
        .output()
        .unwrap();
    let unrepaired = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["phase", "-F"])
        .arg(&input)
        .output()
        .unwrap();
    for output in [&repaired, &unrepaired] {
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    assert_eq!(
        repaired
            .stdout
            .windows(b"YF:i:1".len())
            .filter(|window| *window == b"YF:i:1")
            .count(),
        2
    );
    assert!(
        !unrepaired
            .stdout
            .windows(b"YF:i:1".len())
            .any(|window| window == b"YF:i:1")
    );

    let prefix = directory.path().join("repaired");
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["phase", "--no-pg", "-b"])
        .arg(&prefix)
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let partitioned: Vec<_> = ["0.bam", "1.bam", "chimera.bam"]
        .iter()
        .flat_map(|suffix| {
            bam_records(&PathBuf::from(format!("{}.{}", prefix.display(), suffix)))
                .1
                .into_iter()
                .map(move |record| (*suffix, record))
        })
        .collect();
    assert_eq!(
        partitioned
            .iter()
            .filter(|(_, (name, _))| name == "head" || name == "tail")
            .map(|(suffix, (name, tagged))| (*suffix, name.as_str(), *tagged))
            .collect::<Vec<_>>(),
        [
            ("chimera.bam", "head", false),
            ("chimera.bam", "tail", false)
        ]
    );
}

#[test]
fn no_variant_partition_does_not_discard_input_records() {
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("no-variants");
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["phase", "--no-pg", "-q", "1000", "-b"])
        .arg(&prefix)
        .arg(fixture("phase.sam"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let total: usize = ["0.bam", "1.bam", "chimera.bam"]
        .iter()
        .map(|suffix| {
            bam_records(&PathBuf::from(format!("{}.{}", prefix.display(), suffix)))
                .1
                .len()
        })
        .sum();
    assert_eq!(total, 12);
}

#[test]
fn phase_identity_uses_the_complete_read_name() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("colliding-names.sam");
    let source = fs::read_to_string(fixture("phase.sam"))
        .unwrap()
        .replacen("\nr1\t", "\nAa\t", 1)
        .replacen("\nr2\t", "\nBB\t", 1);
    fs::write(&input, source).unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .arg("phase")
        .arg(input)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output
            .stdout
            .split(|&byte| byte == b'\n')
            .filter(|line| line.starts_with(b"EV\t"))
            .count(),
        12
    );
}

#[test]
fn long_uninformative_reads_do_not_delay_phase_assignments() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("long-read.sam");
    let mut source = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n");
    source.push_str(&format!(
        "long\t0\tchr1\t90\t0\t100M\t*\t0\t0\t{}\t{}\n",
        "A".repeat(100),
        "!".repeat(100)
    ));
    let sequences = [
        "ACGTACGTAA",
        "ACGTACGTAA",
        "TCGTACGTAA",
        "TCGTACGTAA",
        "ACGTTCGTAA",
        "ACGTTCGTAA",
        "TCGTTCGTAA",
        "TCGTTCGTAA",
        "ACGTACGTAA",
        "TCGTTCGTAA",
        "ACGTACGTAA",
        "TCGTTCGTAA",
    ];
    for (prefix, position) in [("r", 100), ("s", 150)] {
        for (index, sequence) in sequences.iter().enumerate() {
            source.push_str(&format!(
                "{prefix}{}\t0\tchr1\t{position}\t60\t10M\t*\t0\t0\t{sequence}\tIIIIIIIIII\n",
                index + 1
            ));
        }
    }
    fs::write(&input, source).unwrap();

    let prefix = directory.path().join("partition");
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["phase", "--no-pg", "-b"])
        .arg(&prefix)
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let records: Vec<_> = ["0.bam", "1.bam", "chimera.bam"]
        .iter()
        .flat_map(|suffix| {
            bam_records(&PathBuf::from(format!("{}.{}", prefix.display(), suffix))).1
        })
        .collect();
    assert_eq!(records.len(), 25);
    assert_eq!(
        records
            .iter()
            .filter(|(name, tagged)| name.starts_with('r') && *tagged)
            .count(),
        8
    );
}

#[test]
fn cram_partitions_require_and_use_an_indexed_reference() {
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("haplotypes");
    let without_reference = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["phase", "--no-pg", "-O", "cram", "-b"])
        .arg(&prefix)
        .arg(fixture("phase.sam"))
        .output()
        .unwrap();
    assert!(!without_reference.status.success());
    assert!(
        String::from_utf8_lossy(&without_reference.stderr)
            .contains("CRAM partition output requires --reference")
    );
    for suffix in ["0.cram", "1.cram", "chimera.cram"] {
        assert!(!PathBuf::from(format!("{}.{}", prefix.display(), suffix)).exists());
    }

    let reference = directory.path().join("reference.fa");
    fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(1000))).unwrap();
    fs::write(
        PathBuf::from(format!("{}.fai", reference.display())),
        b"chr1\t1000\t6\t1000\t1001\n",
    )
    .unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["phase", "--no-pg", "-O", "cram", "-b"])
        .arg(&prefix)
        .arg("--reference")
        .arg(&reference)
        .arg(fixture("phase.sam"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let total: usize = ["0.cram", "1.cram", "chimera.cram"]
        .iter()
        .map(|suffix| {
            let path = PathBuf::from(format!("{}.{}", prefix.display(), suffix));
            let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
                .args(["view", "--no-pg", "--reference"])
                .arg(&reference)
                .arg(path)
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            output
                .stdout
                .split(|&byte| byte == b'\n')
                .filter(|line| !line.is_empty())
                .count()
        })
        .sum();
    assert_eq!(total, 12);
}

#[test]
fn stdin_and_json_report_output_use_the_same_phase_contract() {
    let mut child = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["phase", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&fs::read(fixture("phase.sam")).unwrap())
        .unwrap();
    let streamed = child.wait_with_output().unwrap();
    assert!(
        streamed.status.success(),
        "{}",
        String::from_utf8_lossy(&streamed.stderr)
    );
    assert_eq!(
        normalize_evidence(&streamed.stdout),
        normalize_evidence(include_bytes!("golden/phase-default.txt"))
    );

    let directory = tempfile::tempdir().unwrap();
    let report = directory.path().join("phase.txt");
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["--json", "phase", "-o"])
        .arg(&report)
        .arg(fixture("phase.sam"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let json: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(json["result"]["command"], "phase");
    assert_eq!(json["result"]["summary"]["phase_sets"], 1);
    assert_eq!(json["result"]["summary"]["heterozygous_sites"], 2);
    assert_eq!(
        normalize_evidence(&fs::read(report).unwrap()),
        normalize_evidence(include_bytes!("golden/phase-default.txt"))
    );
}

#[test]
fn report_and_partitions_commit_as_one_group() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("malformed.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nbroken\n",
    )
    .unwrap();
    let report = directory.path().join("phase.txt");
    let prefix = directory.path().join("haplotypes");
    fs::write(&report, b"report sentinel\n").unwrap();
    let partitions = ["0.sam", "1.sam", "chimera.sam"].map(|suffix| {
        let path = PathBuf::from(format!("{}.{}", prefix.display(), suffix));
        fs::write(&path, format!("{suffix} sentinel\n")).unwrap();
        path
    });

    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["phase", "-O", "sam", "-o"])
        .arg(&report)
        .arg("-b")
        .arg(&prefix)
        .arg(&input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_eq!(fs::read(&report).unwrap(), b"report sentinel\n");
    for (path, suffix) in partitions.iter().zip(["0.sam", "1.sam", "chimera.sam"]) {
        assert_eq!(
            fs::read(path).unwrap(),
            format!("{suffix} sentinel\n").as_bytes()
        );
    }
}

#[test]
fn phase_propagates_report_write_failures() {
    struct Closed;

    impl Write for Closed {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::new(
                std::io::ErrorKind::BrokenPipe,
                "closed",
            ))
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Err(std::io::Error::new(
                std::io::ErrorKind::BrokenPipe,
                "closed",
            ))
        }
    }

    let error = rsomics_bam::phase::write(
        &fixture("phase.sam"),
        rsomics_bam::phase::Options::default(),
        Closed,
    )
    .unwrap_err();
    assert!(error.to_string().contains("closed"));
}

#[test]
fn colliding_report_and_partition_targets_are_rejected_without_changes() {
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("partition");
    let report = PathBuf::from(format!("{}.0.bam", prefix.display()));
    fs::write(&report, b"sentinel\n").unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["phase", "-o"])
        .arg(&report)
        .arg("-b")
        .arg(&prefix)
        .arg(fixture("phase.sam"))
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("phase report and partition outputs require different files")
    );
    assert_eq!(fs::read(report).unwrap(), b"sentinel\n");
    for suffix in ["1.bam", "chimera.bam"] {
        assert!(!PathBuf::from(format!("{}.{}", prefix.display(), suffix)).exists());
    }
}

#[test]
fn unsorted_input_fails_instead_of_emitting_a_partial_report() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("unsorted.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
late\t0\tchr1\t200\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n\
early\t0\tchr1\t100\t60\t10M\t*\t0\t0\tTTTTTTTTTT\tIIIIIIIIII\n",
    )
    .unwrap();
    let report = directory.path().join("phase.txt");
    fs::write(&report, b"sentinel\n").unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["phase", "-o"])
        .arg(&report)
        .arg(input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_eq!(fs::read(report).unwrap(), b"sentinel\n");
}
