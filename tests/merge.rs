use std::fs;
use std::path::Path;
use std::process::{Command, Output};

use noodles::bam;
use noodles::sam::header::record::value::map::read_group::tag as read_group_tag;
use rsomics_bamio::raw::RecordReader;

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
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
    assert!(
        !output.status.success(),
        "stdout: {}",
        String::from_utf8_lossy(&output.stdout)
    );
    output
}

fn write_sam(path: &Path, header: &str, records: &[&str]) {
    let mut source = header.as_bytes().to_vec();
    if !source.ends_with(b"\n") {
        source.push(b'\n');
    }
    for record in records {
        source.extend_from_slice(record.as_bytes());
        source.push(b'\n');
    }
    fs::write(path, source).unwrap();
}

fn snapshot(path: &Path) -> (noodles::sam::Header, Vec<String>) {
    let mut reader = bam::io::reader::Builder.build_from_path(path).unwrap();
    let header = reader.read_header().unwrap();
    let mut records = RecordReader::new(reader.get_mut());
    let mut values = Vec::new();
    while let Some(record) = records.next().unwrap() {
        let rg = record
            .aux_value(*b"RG")
            .and_then(|value| value.strip_suffix(&[0]))
            .unwrap_or_default();
        let pg = record
            .aux_value(*b"PG")
            .and_then(|value| value.strip_suffix(&[0]))
            .unwrap_or_default();
        values.push(format!(
            "{}\t{}\t{}\t{}\t{}",
            String::from_utf8_lossy(record.name()),
            record.reference_sequence_id(),
            record.alignment_start(),
            String::from_utf8_lossy(rg),
            String::from_utf8_lossy(pg)
        ));
    }
    (header, values)
}

fn samtools(arguments: &[&str]) -> Output {
    let output = Command::new("samtools").args(arguments).output().unwrap();
    assert!(
        output.status.success(),
        "samtools stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn samtools_records(path: &Path) -> Vec<Vec<u8>> {
    samtools(&["view", path.to_str().unwrap()])
        .stdout
        .split(|byte| *byte == b'\n')
        .filter(|line| !line.is_empty())
        .map(|line| {
            let mut fields = line
                .split(|byte| *byte == b'\t')
                .map(<[u8]>::to_vec)
                .collect::<Vec<_>>();
            let mut auxiliary = fields.split_off(11);
            auxiliary.retain(|field| !field.starts_with(b"MD:Z:") && !field.starts_with(b"NM:i:"));
            auxiliary.sort();
            fields.extend(auxiliary);
            fields.join(&b'\t')
        })
        .collect()
}

fn samtools_core_records(path: &Path) -> Vec<Vec<u8>> {
    samtools(&["view", path.to_str().unwrap()])
        .stdout
        .split(|byte| *byte == b'\n')
        .filter(|line| !line.is_empty())
        .map(|line| {
            line.split(|byte| *byte == b'\t')
                .take(11)
                .collect::<Vec<_>>()
                .join(&b'\t')
        })
        .collect()
}

#[test]
fn coordinate_merge_reconciles_headers_and_translates_records() {
    let directory = tempfile::tempdir().unwrap();
    let a = directory.path().join("a.sam");
    let b = directory.path().join("b.sam");
    let output = directory.path().join("merged.bam");
    write_sam(
        &a,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n@PG\tID:align\tPN:first\n@RG\tID:lane\tSM:first\tPG:align\n@CO\tfirst",
        &[
            "a1\t0\tchr1\t10\t60\t1M\t*\t0\t0\tA\tF\tRG:Z:lane\tPG:Z:align",
            "a2\t0\tchr2\t20\t60\t1M\t*\t0\t0\tC\tF\tRG:Z:lane\tPG:Z:align",
        ],
    );
    write_sam(
        &b,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n@SQ\tSN:chr3\tLN:1000\n@PG\tID:align\tPN:second\n@RG\tID:lane\tSM:second\tPG:align\n@CO\tsecond",
        &[
            "b1\t0\tchr1\t5\t60\t1M\t*\t0\t0\tG\tF\tRG:Z:lane\tPG:Z:align",
            "b2\t0\tchr3\t30\t60\t1M\t*\t0\t0\tT\tF\tRG:Z:lane\tPG:Z:align",
        ],
    );

    run(&[
        "merge",
        "--no-PG",
        "-@",
        "0",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    let (header, records) = snapshot(&output);
    assert_eq!(
        header
            .reference_sequences()
            .keys()
            .map(|name| name.as_slice())
            .collect::<Vec<_>>(),
        [b"chr1".as_slice(), b"chr2".as_slice(), b"chr3".as_slice()]
    );
    assert_eq!(
        records,
        [
            "b1\t0\t4\tlane.1\talign.1",
            "a1\t0\t9\tlane\talign",
            "a2\t1\t19\tlane\talign",
            "b2\t2\t29\tlane.1\talign.1",
        ]
    );
    assert!(header.read_groups().contains_key(b"lane".as_slice()));
    assert!(header.read_groups().contains_key(b"lane.1".as_slice()));
    assert!(header.programs().as_ref().contains_key(b"align".as_slice()));
    assert!(
        header
            .programs()
            .as_ref()
            .contains_key(b"align.1".as_slice())
    );
    assert_eq!(
        header.read_groups()[b"lane.1".as_slice()]
            .other_fields()
            .get(&read_group_tag::PROGRAM)
            .map(|value| value.as_slice()),
        Some(b"align.1".as_slice())
    );
    assert_eq!(header.comments().len(), 2);
}

#[test]
fn combine_modes_keep_first_colliding_header_records() {
    let directory = tempfile::tempdir().unwrap();
    let a = directory.path().join("a.sam");
    let b = directory.path().join("b.sam");
    let output = directory.path().join("merged.bam");
    let header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@PG\tID:align\tPN:one\n@RG\tID:lane\tSM:one\tPG:align";
    write_sam(
        &a,
        header,
        &["a\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tF\tRG:Z:lane\tPG:Z:align"],
    );
    write_sam(
        &b,
        header,
        &["b\t0\tchr1\t2\t60\t1M\t*\t0\t0\tC\tF\tRG:Z:lane\tPG:Z:align"],
    );
    run(&[
        "merge",
        "-c",
        "-p",
        "--no-PG",
        "-@",
        "0",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    let (header, records) = snapshot(&output);
    assert_eq!(header.read_groups().len(), 1);
    assert_eq!(header.programs().as_ref().len(), 1);
    assert!(records.iter().all(|record| record.ends_with("lane\talign")));
}

#[test]
fn invalid_inputs_preserve_existing_output() {
    let directory = tempfile::tempdir().unwrap();
    let valid = directory.path().join("valid.sam");
    let conflict = directory.path().join("conflict.sam");
    let unknown = directory.path().join("unknown.sam");
    let output = directory.path().join("merged.bam");
    write_sam(
        &valid,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:lane",
        &["v\t0\tchr1\t2\t60\t1M\t*\t0\t0\tA\tF\tRG:Z:lane"],
    );
    write_sam(
        &conflict,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:999",
        &["c\t0\tchr1\t3\t60\t1M\t*\t0\t0\tC\tF"],
    );
    write_sam(
        &unknown,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000",
        &["u\t0\tchr1\t1\t60\t1M\t*\t0\t0\tG\tF\tRG:Z:missing"],
    );
    fs::write(&output, b"sentinel").unwrap();

    let result = fail(&[
        "merge",
        valid.to_str().unwrap(),
        conflict.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    assert!(String::from_utf8_lossy(&result.stderr).contains("conflicting header"));
    assert_eq!(fs::read(&output).unwrap(), b"sentinel");

    let result = fail(&[
        "merge",
        unknown.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    assert!(String::from_utf8_lossy(&result.stderr).contains("unknown read group"));
    assert_eq!(fs::read(&output).unwrap(), b"sentinel");
}

#[test]
fn actual_order_is_checked_and_stdout_is_complete_bam() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("unordered.sam");
    let output = directory.path().join("merged.bam");
    write_sam(
        &input,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000",
        &[
            "late\t0\tchr1\t20\t60\t1M\t*\t0\t0\tA\tF",
            "early\t0\tchr1\t10\t60\t1M\t*\t0\t0\tC\tF",
        ],
    );
    fs::write(&output, b"sentinel").unwrap();
    let result = fail(&[
        "merge",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    assert!(String::from_utf8_lossy(&result.stderr).contains("not ordered"));
    assert_eq!(fs::read(&output).unwrap(), b"sentinel");

    write_sam(
        &input,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000",
        &["only\t0\tchr1\t10\t60\t1M\t*\t0\t0\tA\tF"],
    );
    let result = run(&["merge", "--no-PG", "-@", "0", input.to_str().unwrap()]);
    let mut reader = bam::io::Reader::new(result.stdout.as_slice());
    reader.read_header().unwrap();
    let mut records = RecordReader::new(reader.get_mut());
    assert!(records.next().unwrap().is_some());
    assert!(records.next().unwrap().is_none());
}

#[test]
fn incompatible_dictionaries_headers_and_aliases_fail() {
    let directory = tempfile::tempdir().unwrap();
    let first = directory.path().join("first.sam");
    let reordered = directory.path().join("reordered.sam");
    let bad_program = directory.path().join("bad-program.sam");
    let wrong_order = directory.path().join("wrong-order.sam");
    let output = directory.path().join("output.bam");
    write_sam(
        &first,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:100",
        &[],
    );
    write_sam(
        &reordered,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr2\tLN:100\n@SQ\tSN:chr1\tLN:100",
        &[],
    );
    write_sam(
        &bad_program,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n@PG\tID:child\tPP:missing",
        &[],
    );
    write_sam(
        &wrong_order,
        "@HD\tVN:1.6\tSO:queryname\tSS:queryname:natural\n@SQ\tSN:chr1\tLN:100",
        &[],
    );

    let result = fail(&[
        "merge",
        first.to_str().unwrap(),
        reordered.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    assert!(String::from_utf8_lossy(&result.stderr).contains("destroy input ordering"));
    let result = fail(&[
        "merge",
        bad_program.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    assert!(String::from_utf8_lossy(&result.stderr).contains("unknown previous program"));
    let result = fail(&[
        "merge",
        wrong_order.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    assert!(String::from_utf8_lossy(&result.stderr).contains("does not declare coordinate"));

    let alias = directory.path().join("alias.sam");
    fs::hard_link(&first, &alias).unwrap();
    let result = fail(&[
        "merge",
        first.to_str().unwrap(),
        "-o",
        alias.to_str().unwrap(),
    ]);
    assert!(String::from_utf8_lossy(&result.stderr).contains("different files"));
}

#[test]
fn input_count_is_bounded() {
    let mut arguments = vec!["merge"];
    arguments.extend(std::iter::repeat_n("missing.sam", 33));
    let result = fail(&arguments);
    assert!(String::from_utf8_lossy(&result.stderr).contains("at most 32 inputs"));
}

#[test]
#[ignore = "requires samtools 1.24"]
fn all_orders_match_samtools_1_24() {
    assert!(String::from_utf8_lossy(&samtools(&["--version"]).stdout).starts_with("samtools 1.24"));
    let directory = tempfile::tempdir().unwrap();
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/sort-records.sam");
    let cases: &[(&[&str], &[&str])] = &[
        (&[], &[]),
        (&["-n"], &["-n"]),
        (&["-N"], &["-N"]),
        (&["--template-coordinate"], &["--template-coordinate"]),
    ];

    for (index, (sort_flags, merge_flags)) in cases.iter().enumerate() {
        let input = directory.path().join(format!("input-{index}.bam"));
        let ours = directory.path().join(format!("ours-{index}.bam"));
        let upstream = directory.path().join(format!("upstream-{index}.bam"));
        let mut sort_args = vec!["sort"];
        sort_args.extend_from_slice(sort_flags);
        sort_args.extend(["-@", "0", "-o", input.to_str().unwrap()]);
        sort_args.push(fixture.to_str().unwrap());
        samtools(&sort_args);

        let mut ours_args = vec!["merge"];
        ours_args.extend_from_slice(merge_flags);
        ours_args.extend([
            "--no-PG",
            "-c",
            "-p",
            "-@",
            "0",
            input.to_str().unwrap(),
            input.to_str().unwrap(),
            "-o",
            ours.to_str().unwrap(),
        ]);
        run(&ours_args);

        let mut upstream_args = vec!["merge"];
        upstream_args.extend_from_slice(merge_flags);
        upstream_args.extend([
            "--no-PG",
            "-c",
            "-p",
            "-@",
            "0",
            "-f",
            "-o",
            upstream.to_str().unwrap(),
            input.to_str().unwrap(),
            input.to_str().unwrap(),
        ]);
        samtools(&upstream_args);
        assert_eq!(samtools_records(&ours), samtools_records(&upstream));
    }
}

#[test]
#[ignore = "requires samtools 1.24"]
fn sam_bam_and_cram_inputs_match_samtools_1_24() {
    assert!(String::from_utf8_lossy(&samtools(&["--version"]).stdout).starts_with("samtools 1.24"));
    let directory = tempfile::tempdir().unwrap();
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/sort-records.sam");
    let bam = directory.path().join("input.bam");
    let sam = directory.path().join("input.sam");
    let cram = directory.path().join("input.cram");
    let reference = directory.path().join("reference.fa");
    let ours = directory.path().join("ours.bam");
    let upstream = directory.path().join("upstream.bam");

    samtools(&[
        "sort",
        "-@",
        "0",
        "-o",
        bam.to_str().unwrap(),
        fixture.to_str().unwrap(),
    ]);
    fs::write(
        &reference,
        format!(">chr1\n{}\n>chr2\n{}\n", "A".repeat(1000), "C".repeat(1000)),
    )
    .unwrap();
    samtools(&["faidx", reference.to_str().unwrap()]);
    fs::write(
        &sam,
        samtools(&["view", "-h", bam.to_str().unwrap()]).stdout,
    )
    .unwrap();
    samtools(&[
        "view",
        "-C",
        "-T",
        reference.to_str().unwrap(),
        "-o",
        cram.to_str().unwrap(),
        bam.to_str().unwrap(),
    ]);

    run(&[
        "merge",
        "--no-PG",
        "-c",
        "-p",
        "-@",
        "0",
        "--reference",
        reference.to_str().unwrap(),
        sam.to_str().unwrap(),
        bam.to_str().unwrap(),
        cram.to_str().unwrap(),
        "-o",
        ours.to_str().unwrap(),
    ]);
    samtools(&[
        "merge",
        "--no-PG",
        "-c",
        "-p",
        "-@",
        "0",
        "--reference",
        reference.to_str().unwrap(),
        "-f",
        "-o",
        upstream.to_str().unwrap(),
        sam.to_str().unwrap(),
        bam.to_str().unwrap(),
        cram.to_str().unwrap(),
    ]);
    let mut ours_records = samtools_core_records(&ours);
    let mut upstream_records = samtools_core_records(&upstream);
    ours_records.sort();
    upstream_records.sort();
    assert_eq!(ours_records, upstream_records);
}
