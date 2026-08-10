use std::fs;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use noodles::bam;
use rsomics_bamio::raw::RecordReader;

#[derive(Debug)]
struct Record {
    name: Vec<u8>,
    flags: u16,
    reference: i32,
    position: i32,
    mate_reference: i32,
    mate_position: i32,
    template_length: i32,
    mc: Option<Vec<u8>>,
    mq: Option<u32>,
    ms: Option<u32>,
}

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/fixmate.sam")
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

fn records(path: &Path) -> (noodles::sam::Header, Vec<Record>) {
    let mut reader = bam::io::reader::Builder.build_from_path(path).unwrap();
    let header = reader.read_header().unwrap();
    let mut raw = RecordReader::new(reader.get_mut());
    let mut records = Vec::new();
    while let Some(record) = raw.next().unwrap() {
        records.push(Record {
            name: record.name().to_vec(),
            flags: record.flags(),
            reference: record.reference_sequence_id(),
            position: record.alignment_start(),
            mate_reference: record.mate_reference_sequence_id(),
            mate_position: record.mate_alignment_start(),
            template_length: record.template_length(),
            mc: record.aux_value(*b"MC").map(ToOwned::to_owned),
            mq: integer_tag(&record, *b"MQ"),
            ms: integer_tag(&record, *b"ms"),
        });
    }
    (header, records)
}

fn integer_tag(record: &rsomics_bamio::raw::RecordRef<'_>, tag: [u8; 2]) -> Option<u32> {
    record
        .aux_value(tag)
        .and_then(|value| value.try_into().ok())
        .map(u32::from_le_bytes)
}

fn named<'a>(records: &'a [Record], name: &[u8]) -> Vec<&'a Record> {
    records
        .iter()
        .filter(|record| record.name == name)
        .collect()
}

#[test]
fn repairs_primary_supplementary_orphan_and_multi_primary_records() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("fixed.bam");
    let result = run(&[
        "--json",
        "fixmate",
        "-m",
        "--no-PG",
        "-@",
        "0",
        fixture().to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    let value: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(value["result"]["summary"]["records"], 21);
    assert_eq!(value["result"]["summary"]["written_records"], 21);
    assert_eq!(value["result"]["summary"]["templates"], 9);
    assert_eq!(value["result"]["summary"]["paired_templates"], 8);

    let (header, values) = records(&output);
    assert!(header.programs().as_ref().is_empty());

    let pair = named(&values, b"pair");
    assert_eq!(pair[0].mate_position, 999);
    assert_eq!(pair[0].template_length, 108);
    assert_eq!(pair[0].mc.as_deref(), Some(b"2S8M\0".as_slice()));
    assert_eq!(pair[0].mq, Some(45));
    assert_eq!(pair[0].ms, Some(400));
    assert_eq!(pair[1].template_length, -108);
    assert_eq!(pair[1].mc.as_deref(), Some(b"8M2S\0".as_slice()));
    assert_eq!(pair[1].mq, Some(60));
    assert_eq!(pair[1].ms, Some(315));

    let missing = named(&values, b"missing");
    assert_eq!(missing[0].ms, Some(2_550));
    assert_eq!(missing[1].ms, Some(2_550));

    let multi = named(&values, b"multi");
    assert_eq!(multi[0].mate_position, 699);
    assert_eq!(multi[0].mq, Some(42));
    assert_eq!(multi[3].mate_position, 699);
    assert_eq!(multi[3].mq, Some(42));
    assert_eq!(multi[3].mc.as_deref(), Some(b"10M\0".as_slice()));

    let orphan = named(&values, b"orphan");
    assert_eq!(orphan[0].flags, 65);
    assert_eq!(orphan[0].mate_reference, -1);
    assert_eq!(orphan[0].mate_position, -1);
    assert_eq!(orphan[0].template_length, 0);

    let supplementary = named(&values, b"supp");
    assert_eq!(supplementary[2].mate_reference, 0);
    assert_eq!(supplementary[2].mate_position, 1399);
    assert_ne!(supplementary[2].flags & 0x20, 0);
    assert_eq!(supplementary[2].mq, Some(55));
    assert_eq!(supplementary[2].mc.as_deref(), Some(b"8M2S\0".as_slice()));

    let cross = named(&values, b"cross");
    assert_eq!(cross[0].flags & 0x2, 0);
    assert_eq!(cross[1].flags & 0x2, 0);
    assert_eq!(cross[0].template_length, 0);
    assert_eq!(cross[1].template_length, 0);

    let unmapped = named(&values, b"unmapped");
    assert_eq!(unmapped[1].reference, 0);
    assert_eq!(unmapped[1].position, 1499);
    assert_eq!(unmapped[1].mate_reference, 0);
    assert_eq!(unmapped[1].mate_position, 1499);
    assert_eq!(unmapped[1].mq, Some(60));
    assert_eq!(unmapped[1].mc.as_deref(), Some(b"10M\0".as_slice()));

    let wrong = named(&values, b"wrong");
    assert_eq!(wrong[0].flags & 0x2, 0);
    assert_eq!(wrong[1].flags & 0x2, 0);
}

#[test]
fn remove_and_no_proper_pair_match_their_declared_scope() {
    let directory = tempfile::tempdir().unwrap();
    let removed = directory.path().join("removed.bam");
    run(&[
        "fixmate",
        "-r",
        "--no-pg",
        "-@",
        "0",
        fixture().to_str().unwrap(),
        "-o",
        removed.to_str().unwrap(),
    ]);
    let (_, removed_records) = records(&removed);
    assert_eq!(removed_records.len(), 19);
    assert!(
        removed_records
            .iter()
            .all(|record| record.flags & 0x104 == 0)
    );

    let unchecked = directory.path().join("unchecked.bam");
    run(&[
        "fixmate",
        "-p",
        "--no-pg",
        "-@",
        "0",
        fixture().to_str().unwrap(),
        "-o",
        unchecked.to_str().unwrap(),
    ]);
    let (_, unchecked_records) = records(&unchecked);
    let wrong = named(&unchecked_records, b"wrong");
    assert_ne!(wrong[0].flags & 0x2, 0);
    assert_ne!(wrong[1].flags & 0x2, 0);
}

#[test]
fn long_cigar_populates_the_complete_mate_cigar() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("long.sam");
    let output = directory.path().join("long.bam");
    let operation_count = usize::from(u16::MAX) + 1;
    let cigar = "1M".repeat(operation_count);
    let sequence = "A".repeat(operation_count);
    let quality = "I".repeat(operation_count);
    fs::write(
        &input,
        format!(
            "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:70000\nlong\t65\tchr1\t1\t60\t{cigar}\t*\t0\t0\t{sequence}\t{quality}\nlong\t145\tchr1\t65537\t50\t1M\t*\t0\t0\tA\tI\n"
        ),
    )
    .unwrap();
    run(&[
        "fixmate",
        "--no-pg",
        "-@",
        "0",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    let (_, records) = records(&output);
    assert_eq!(records[0].mc.as_deref(), Some(b"1M\0".as_slice()));
    assert_eq!(records[1].mc.as_ref().unwrap().len(), cigar.len() + 1);
    assert_eq!(records[1].mc.as_ref().unwrap().last(), Some(&0));
    assert_eq!(records[0].template_length, 65_537);
    assert_eq!(records[1].template_length, -65_537);
}

#[test]
fn coordinate_input_and_path_aliases_fail_transactionally() {
    let directory = tempfile::tempdir().unwrap();
    let coordinate = directory.path().join("coordinate.sam");
    let source = fs::read_to_string(fixture()).unwrap();
    fs::write(&coordinate, source.replace("SO:queryname", "SO:coordinate")).unwrap();
    let output = directory.path().join("output.bam");
    fs::write(&output, b"keep").unwrap();
    let result = fail(&[
        "fixmate",
        coordinate.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ]);
    assert!(String::from_utf8_lossy(&result.stderr).contains("coordinate sorted"));
    assert_eq!(fs::read(&output).unwrap(), b"keep");

    let original = fs::read(fixture()).unwrap();
    let result = fail(&[
        "fixmate",
        fixture().to_str().unwrap(),
        "-o",
        fixture().to_str().unwrap(),
    ]);
    assert!(String::from_utf8_lossy(&result.stderr).contains("different files"));
    assert_eq!(fs::read(fixture()).unwrap(), original);
}

#[test]
fn stdout_is_a_complete_bam_and_help_excludes_unimplemented_modes() {
    let directory = tempfile::tempdir().unwrap();
    let result = run(&["fixmate", "--no-pg", "-@", "0", fixture().to_str().unwrap()]);
    let output = directory.path().join("stdout.bam");
    fs::write(&output, result.stdout).unwrap();
    assert_eq!(records(&output).1.len(), 21);

    let help = run(&["fixmate", "--help"]);
    let text = String::from_utf8(help.stdout).unwrap();
    assert!(text.contains("-m, --mate-score"));
    assert!(text.contains("--no-PG"));
    assert!(!text.contains("template-cigar"));
    assert!(!text.contains("base-mod"));

    let piped = directory.path().join("piped.bam");
    let mut child = bin()
        .args([
            "fixmate",
            "--no-pg",
            "-@",
            "0",
            "-",
            "-o",
            piped.to_str().unwrap(),
        ])
        .stdin(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&fs::read(fixture()).unwrap())
        .unwrap();
    assert!(child.wait().unwrap().success());
    assert_eq!(records(&piped).1.len(), 21);
}

#[test]
#[ignore = "requires samtools 1.24"]
fn samtools_1_24_matrix_matches_all_supported_modes() {
    assert_eq!(samtools_version(), "samtools 1.24");
    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    write_reference(&reference);
    samtools(&["faidx", reference.to_str().unwrap()]);

    let bam = directory.path().join("input.bam");
    samtools(&[
        "view",
        "--no-PG",
        "-b",
        "-o",
        bam.to_str().unwrap(),
        fixture().to_str().unwrap(),
    ]);
    let cram = directory.path().join("input.cram");
    samtools(&[
        "view",
        "--no-PG",
        "-C",
        "-T",
        reference.to_str().unwrap(),
        "-o",
        cram.to_str().unwrap(),
        fixture().to_str().unwrap(),
    ]);

    for input in [fixture(), bam, cram] {
        for options in [
            Vec::<&str>::new(),
            vec!["-m"],
            vec!["-r"],
            vec!["-p"],
            vec!["-m", "-r", "-p"],
        ] {
            let stem = format!(
                "{}-{}",
                input.extension().unwrap().to_string_lossy(),
                options.join("").replace('-', "")
            );
            let expected = directory.path().join(format!("expected-{stem}.bam"));
            let actual = directory.path().join(format!("actual-{stem}.bam"));

            let mut upstream = vec!["fixmate", "-z", "off", "--no-PG"];
            upstream.extend(options.iter().copied());
            upstream.push(input.to_str().unwrap());
            upstream.push(expected.to_str().unwrap());
            samtools(&upstream);

            let mut ours = vec!["fixmate", "--no-pg", "-@", "0"];
            ours.extend(options.iter().copied());
            if input
                .extension()
                .is_some_and(|extension| extension == "cram")
            {
                ours.extend(["--reference", reference.to_str().unwrap()]);
            }
            ours.push(input.to_str().unwrap());
            ours.extend(["-o", actual.to_str().unwrap()]);
            run(&ours);

            assert_eq!(canonical_sam(&actual), canonical_sam(&expected), "{stem}");
        }
    }
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

fn samtools_version() -> String {
    let output = samtools(&["--version"]);
    String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .to_owned()
}

fn canonical_sam(path: &Path) -> Vec<String> {
    let output = samtools(&["view", "--no-PG", "-h", path.to_str().unwrap()]);
    String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|line| {
            if line.starts_with('@') {
                return line.to_owned();
            }
            let mut fields = line.split('\t').map(str::to_owned).collect::<Vec<_>>();
            fields[11..].sort();
            fields.join("\t")
        })
        .collect()
}

fn write_reference(path: &Path) {
    let sequence = "A".repeat(2_000);
    fs::write(path, format!(">chr1\n{sequence}\n>chr2\n{sequence}\n")).unwrap();
}
