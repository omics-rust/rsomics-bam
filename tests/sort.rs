use std::fs::{self, File};
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use noodles::bam;
use noodles::sam::header::record::value::map::{header::tag, tag::Other};
use rsomics_bamio::raw::RecordReader;

type HeaderTag = Other<tag::Standard>;
type SortCase<'a> = (&'a [&'a str], &'a [&'a str], &'a [(HeaderTag, &'a [u8])]);

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/sort-records.sam")
}

fn run(arguments: &[&str]) -> Output {
    let output = bin().args(arguments).output().unwrap();
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn fail(arguments: &[&str]) -> Output {
    let output = bin().args(arguments).output().unwrap();
    assert!(!output.status.success());
    output
}

fn snapshot(path: &Path) -> (Vec<String>, noodles::sam::Header) {
    let mut reader = bam::io::reader::Builder.build_from_path(path).unwrap();
    let header = reader.read_header().unwrap();
    let mut records = RecordReader::new(reader.get_mut());
    let mut values = Vec::new();
    while let Some(record) = records.next().unwrap() {
        values.push(format!(
            "{}\t{}\t{}\t{}",
            String::from_utf8_lossy(record.name()),
            record.flags(),
            record.reference_sequence_id(),
            record.alignment_start()
        ));
    }
    (values, header)
}

fn header_value(header: &noodles::sam::Header, key: HeaderTag) -> Option<&[u8]> {
    header
        .header()?
        .other_fields()
        .get(&key)
        .map(|value| value.as_slice())
}

#[test]
fn coordinate_and_queryname_orders_match_committed_oracles() {
    let directory = tempfile::tempdir().unwrap();
    let cases: &[SortCase<'_>] = &[
        (
            &[],
            &[
                "pair2\t99\t0\t49",
                "single1\t0\t0\t99",
                "single01\t0\t0\t99",
                "tie-forward\t0\t0\t99",
                "tie-reverse\t16\t0\t99",
                "pair2\t147\t0\t99",
                "pair12\t99\t1\t199",
                "pair12\t147\t1\t249",
                "unmapped10\t4\t-1\t-1",
            ],
            &[(tag::SORT_ORDER, b"coordinate")],
        ),
        (
            &["-n"],
            &[
                "pair2\t99\t0\t49",
                "pair2\t147\t0\t99",
                "pair12\t99\t1\t199",
                "pair12\t147\t1\t249",
                "single1\t0\t0\t99",
                "single01\t0\t0\t99",
                "tie-forward\t0\t0\t99",
                "tie-reverse\t16\t0\t99",
                "unmapped10\t4\t-1\t-1",
            ],
            &[
                (tag::SORT_ORDER, b"queryname"),
                (tag::SUBSORT_ORDER, b"queryname:natural"),
            ],
        ),
        (
            &["-N"],
            &[
                "pair12\t99\t1\t199",
                "pair12\t147\t1\t249",
                "pair2\t99\t0\t49",
                "pair2\t147\t0\t99",
                "single01\t0\t0\t99",
                "single1\t0\t0\t99",
                "tie-forward\t0\t0\t99",
                "tie-reverse\t16\t0\t99",
                "unmapped10\t4\t-1\t-1",
            ],
            &[
                (tag::SORT_ORDER, b"queryname"),
                (tag::SUBSORT_ORDER, b"queryname:lexicographical"),
            ],
        ),
        (
            &["--template-coordinate"],
            &[
                "pair2\t99\t0\t49",
                "pair2\t147\t0\t99",
                "single01\t0\t0\t99",
                "single1\t0\t0\t99",
                "tie-forward\t0\t0\t99",
                "tie-reverse\t16\t0\t99",
                "pair12\t99\t1\t199",
                "pair12\t147\t1\t249",
                "unmapped10\t4\t-1\t-1",
            ],
            &[
                (tag::SORT_ORDER, b"unsorted"),
                (tag::GROUP_ORDER, b"query"),
                (tag::SUBSORT_ORDER, b"unsorted:template-coordinate"),
            ],
        ),
    ];

    let input = fixture();
    for (index, (flags, expected, headers)) in cases.iter().enumerate() {
        let output = directory.path().join(format!("case-{index}.bam"));
        let mut arguments = vec!["sort"];
        arguments.extend_from_slice(flags);
        arguments.push(input.to_str().unwrap());
        arguments.extend(["-o", output.to_str().unwrap(), "-@", "0"]);
        run(&arguments);
        let (actual, header) = snapshot(&output);
        assert_eq!(actual, *expected);
        for (key, value) in *headers {
            assert_eq!(header_value(&header, *key), Some(*value));
        }
        if !headers.iter().any(|(key, _)| *key == tag::GROUP_ORDER) {
            assert_eq!(header_value(&header, tag::GROUP_ORDER), None);
        }
        if !headers.iter().any(|(key, _)| *key == tag::SUBSORT_ORDER) {
            assert_eq!(header_value(&header, tag::SUBSORT_ORDER), None);
        }
    }
}

#[test]
fn external_runs_are_merged_and_reported() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("large.sam");
    write_tied_sam(&input, 2_048, 12_000);
    let output = directory.path().join("sorted.bam");
    let result = run(&[
        "--json",
        "sort",
        "-m",
        "1M",
        "-@",
        "0",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    let value: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(value["result"]["summary"]["records"], 2_048);
    assert!(
        value["result"]["summary"]["temporary_runs"]
            .as_u64()
            .unwrap()
            > 32
    );
    assert!(value["result"]["summary"]["merge_passes"].as_u64().unwrap() >= 2);

    let (records, _) = snapshot(&output);
    assert_eq!(records.len(), 2_048);
    assert!(records.windows(2).all(|pair| pair[0] > pair[1]));
    let prefix = output.file_name().unwrap().to_string_lossy();
    assert!(!fs::read_dir(directory.path()).unwrap().any(|entry| {
        entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .starts_with(&format!("{prefix}.tmp."))
    }));
}

#[test]
fn failures_preserve_outputs_and_aliases_are_rejected() {
    let directory = tempfile::tempdir().unwrap();
    let invalid = directory.path().join("missing-mc.sam");
    fs::write(
        &invalid,
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\npair\t1\tchr1\t1\t60\t10M\t=\t20\t0\tAAAAAAAAAA\tFFFFFFFFFF\n",
    )
    .unwrap();
    let output = directory.path().join("output.bam");
    fs::write(&output, b"sentinel").unwrap();
    let failed = fail(&[
        "sort",
        "--template-coordinate",
        invalid.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    assert!(String::from_utf8_lossy(&failed.stderr).contains("no MC tag"));
    assert_eq!(fs::read(&output).unwrap(), b"sentinel");

    let alias = directory.path().join("alias.sam");
    fs::hard_link(&invalid, &alias).unwrap();
    let failed = fail(&[
        "sort",
        invalid.to_str().unwrap(),
        "-o",
        alias.to_str().unwrap(),
    ]);
    assert!(String::from_utf8_lossy(&failed.stderr).contains("different files"));
}

#[test]
fn binary_stdout_is_a_complete_bam() {
    let output = run(&["sort", "-@", "0", fixture().to_str().unwrap()]);
    let mut reader = bam::io::Reader::new(output.stdout.as_slice());
    let header = reader.read_header().unwrap();
    assert_eq!(
        header_value(&header, tag::SORT_ORDER),
        Some(b"coordinate".as_slice())
    );
    let mut records = RecordReader::new(reader.get_mut());
    let mut count = 0;
    while records.next().unwrap().is_some() {
        count += 1;
    }
    assert_eq!(count, 9);
}

#[test]
fn uncompressed_bam_is_supported_and_missing_bgzf_eof_fails() {
    let directory = tempfile::tempdir().unwrap();
    let raw = directory.path().join("raw.bam");
    run(&[
        "view",
        "-u",
        fixture().to_str().unwrap(),
        "-o",
        raw.to_str().unwrap(),
    ]);
    let sorted = directory.path().join("raw-sorted.bam");
    run(&[
        "sort",
        "-@",
        "0",
        raw.to_str().unwrap(),
        "-o",
        sorted.to_str().unwrap(),
    ]);
    assert_eq!(snapshot(&sorted).0.len(), 9);

    let compressed = directory.path().join("compressed.bam");
    run(&[
        "view",
        "-b",
        fixture().to_str().unwrap(),
        "-o",
        compressed.to_str().unwrap(),
    ]);
    let mut bytes = fs::read(&compressed).unwrap();
    bytes.truncate(bytes.len() - 56);
    fs::write(&compressed, bytes).unwrap();
    fs::write(&sorted, b"sentinel").unwrap();
    let failed = fail(&[
        "sort",
        compressed.to_str().unwrap(),
        "-o",
        sorted.to_str().unwrap(),
    ]);
    assert!(String::from_utf8_lossy(&failed.stderr).contains("end-of-file marker is missing"));
    assert_eq!(fs::read(&sorted).unwrap(), b"sentinel");
}

#[test]
#[ignore = "requires samtools 1.24"]
fn all_orders_match_samtools_1_24() {
    assert_eq!(samtools_version(), "samtools 1.24");
    let directory = tempfile::tempdir().unwrap();
    let input = fixture();
    for (name, flags) in [
        ("coordinate", Vec::<&str>::new()),
        ("natural", vec!["-n"]),
        ("ascii", vec!["-N"]),
        ("template", vec!["--template-coordinate"]),
    ] {
        let ours = directory.path().join(format!("ours-{name}.bam"));
        let oracle = directory.path().join(format!("oracle-{name}.bam"));
        let mut ours_args = vec!["sort"];
        ours_args.extend_from_slice(&flags);
        ours_args.push(input.to_str().unwrap());
        ours_args.extend(["-o", ours.to_str().unwrap(), "-@", "0"]);
        run(&ours_args);

        let mut command = Command::new("samtools");
        command.arg("sort").args(&flags).args(["-@", "0", "-o"]);
        command.arg(&oracle).arg(&input);
        assert!(command.status().unwrap().success());
        assert_eq!(samtools_records(&ours), samtools_records(&oracle));
        assert_eq!(sort_header(&ours), sort_header(&oracle));
    }
}

#[test]
#[ignore = "requires samtools 1.24"]
fn bam_cram_and_external_runs_match_samtools_1_24() {
    assert_eq!(samtools_version(), "samtools 1.24");
    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    let mut fasta = File::create(&reference).unwrap();
    writeln!(
        fasta,
        ">chr1\n{}\n>chr2\n{}",
        "A".repeat(1000),
        "C".repeat(1000)
    )
    .unwrap();
    assert!(
        Command::new("samtools")
            .args(["faidx"])
            .arg(&reference)
            .status()
            .unwrap()
            .success()
    );

    let bam = directory.path().join("input.bam");
    assert!(
        Command::new("samtools")
            .args(["view", "-b", "-o"])
            .arg(&bam)
            .arg(fixture())
            .status()
            .unwrap()
            .success()
    );
    let cram = directory.path().join("input.cram");
    assert!(
        Command::new("samtools")
            .args(["view", "-C", "-T"])
            .arg(&reference)
            .args(["-o"])
            .arg(&cram)
            .arg(fixture())
            .status()
            .unwrap()
            .success()
    );

    for input in [&bam, &cram] {
        let ours = directory.path().join(format!(
            "ours-{}.bam",
            input.extension().unwrap().to_string_lossy()
        ));
        let oracle = directory.path().join(format!(
            "oracle-{}.bam",
            input.extension().unwrap().to_string_lossy()
        ));
        let mut arguments = vec!["sort", "-@", "0"];
        if input == &cram {
            arguments.extend(["--reference", reference.to_str().unwrap()]);
        }
        arguments.push(input.to_str().unwrap());
        arguments.extend(["-o", ours.to_str().unwrap()]);
        run(&arguments);
        assert!(
            Command::new("samtools")
                .args(["sort", "-@", "0", "-o"])
                .arg(&oracle)
                .arg(input)
                .status()
                .unwrap()
                .success()
        );
        assert_eq!(
            canonical_records(&samtools_records(&ours)),
            canonical_records(&samtools_records(&oracle))
        );
    }

    let large = directory.path().join("large.sam");
    write_large_sam(&large, 20_000, 100);
    let ours = directory.path().join("ours-large.bam");
    let oracle = directory.path().join("oracle-large.bam");
    run(&[
        "sort",
        "-N",
        "-m",
        "1M",
        "-@",
        "0",
        large.to_str().unwrap(),
        "-o",
        ours.to_str().unwrap(),
    ]);
    assert!(
        Command::new("samtools")
            .args(["sort", "-N", "-m", "1M", "-@", "0", "-o"])
            .arg(&oracle)
            .arg(&large)
            .status()
            .unwrap()
            .success()
    );
    assert_eq!(samtools_records(&ours), samtools_records(&oracle));
}

fn write_large_sam(path: &Path, records: usize, read_length: usize) {
    let mut file = File::create(path).unwrap();
    let sequence = "A".repeat(read_length);
    let quality = "F".repeat(read_length);
    writeln!(file, "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000000").unwrap();
    for i in (0..records).rev() {
        writeln!(
            file,
            "read{i:08}\t0\tchr1\t{}\t60\t{read_length}M\t*\t0\t0\t{}\t{}",
            i + 1,
            sequence,
            quality
        )
        .unwrap();
    }
}

fn write_tied_sam(path: &Path, records: usize, read_length: usize) {
    let mut file = File::create(path).unwrap();
    let sequence = "A".repeat(read_length);
    let quality = "F".repeat(read_length);
    writeln!(file, "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000000").unwrap();
    for i in (0..records).rev() {
        writeln!(
            file,
            "read{i:08}\t0\tchr1\t1\t60\t{read_length}M\t*\t0\t0\t{}\t{}",
            sequence, quality
        )
        .unwrap();
    }
}

fn samtools_version() -> String {
    let output = Command::new("samtools").arg("--version").output().unwrap();
    String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .to_owned()
}

fn samtools_records(path: &Path) -> Vec<u8> {
    let output = Command::new("samtools")
        .arg("view")
        .arg(path)
        .output()
        .unwrap();
    assert!(output.status.success());
    output.stdout
}

fn canonical_records(records: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(records)
        .lines()
        .map(|line| {
            let mut fields = line.split('\t').collect::<Vec<_>>();
            fields[11..].sort_unstable();
            fields.join("\t")
        })
        .collect()
}

fn sort_header(path: &Path) -> String {
    let output = Command::new("samtools")
        .args(["view", "-H"])
        .arg(path)
        .output()
        .unwrap();
    assert!(output.status.success());
    String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .find(|line| line.starts_with("@HD"))
        .unwrap()
        .to_owned()
}
