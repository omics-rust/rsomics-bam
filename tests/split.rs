use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles::bam;
use rsomics_bam::split::{Format, Mode, Options, run};
use rsomics_bamio::raw::RecordReader;

type Records = Vec<(Vec<u8>, Vec<u8>)>;

fn fixture(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden/split")
        .join(path)
}

fn snapshot(path: &Path) -> (noodles::sam::Header, Records) {
    let mut reader = bam::io::reader::Builder.build_from_path(path).unwrap();
    let header = reader.read_header().unwrap();
    let mut raw = RecordReader::new(reader.get_mut());
    let mut records = Vec::new();
    while let Some(record) = raw.next().unwrap() {
        let read_group = record
            .aux_value(*b"RG")
            .and_then(|value| value.strip_suffix(&[0]))
            .unwrap_or_default();
        records.push((record.name().to_vec(), read_group.to_vec()));
    }
    (header, records)
}

fn names(path: &Path) -> Vec<Vec<u8>> {
    snapshot(path).1.into_iter().map(|(name, _)| name).collect()
}

fn view_body(path: &Path) -> Vec<u8> {
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
        .args(["view", "--no-pg"])
        .arg(path)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

fn body_names(body: &[u8]) -> Vec<Vec<u8>> {
    body.split(|byte| *byte == b'\n')
        .filter(|line| !line.is_empty())
        .map(|line| line.split(|byte| *byte == b'\t').next().unwrap().to_vec())
        .collect()
}

fn write_sam(path: &Path, header: &str, records: &[&str]) {
    let mut source = header.as_bytes().to_vec();
    source.push(b'\n');
    for record in records {
        source.extend_from_slice(record.as_bytes());
        source.push(b'\n');
    }
    fs::write(path, source).unwrap();
}

#[test]
fn read_group_split_creates_complete_filtered_outputs() {
    let directory = tempfile::tempdir().unwrap();
    let prefix = directory.path().join("sample");
    let summary = run(
        &fixture("read-group/tworg.bam"),
        Options {
            mode: Mode::ReadGroup,
            output_prefix: &prefix,
            unaccounted: None,
            unaccounted_header: None,
            format: Format::Bam,
            maximum_outputs: 100,
            zero_pad: 0,
            reference: None,
            additional_threads: 0,
            program: None,
        },
    )
    .unwrap();

    assert_eq!(summary.records, 9);
    assert_eq!(summary.skipped, 0);
    assert_eq!(summary.outputs.len(), 2);
    assert_eq!(summary.outputs[0].label, "rg1");
    assert_eq!(summary.outputs[0].records, 8);
    assert_eq!(summary.outputs[1].label, "rg2");
    assert_eq!(summary.outputs[1].records, 1);

    let (rg1_header, rg1_records) = snapshot(&directory.path().join("sample.rg1.bam"));
    assert_eq!(
        rg1_header
            .read_groups()
            .keys()
            .map(|id| id.as_slice())
            .collect::<Vec<_>>(),
        [b"rg1".as_slice()]
    );
    assert_eq!(
        rg1_records,
        [
            "read1", "read1", "read2", "read2", "read4", "read5", "read6", "read7"
        ]
        .map(|name| (name.as_bytes().to_vec(), b"rg1".to_vec()))
    );

    let (rg2_header, rg2_records) = snapshot(&directory.path().join("sample.rg2.bam"));
    assert_eq!(
        rg2_header
            .read_groups()
            .keys()
            .map(|id| id.as_slice())
            .collect::<Vec<_>>(),
        [b"rg2".as_slice()]
    );
    assert_eq!(rg2_records, [(b"read4".to_vec(), b"rg2".to_vec())]);
}

#[test]
fn missing_and_unknown_read_groups_are_explicit_and_transactional() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.sam");
    let prefix = directory.path().join("sample");
    let unaccounted = directory.path().join("other.bam");
    let rg1 = directory.path().join("sample.rg1.bam");
    write_sam(
        &input,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@RG\tID:rg1",
        &[
            "known\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tF\tRG:Z:rg1",
            "missing\t0\tchr1\t2\t60\t1M\t*\t0\t0\tC\tF",
            "unknown\t0\tchr1\t3\t60\t1M\t*\t0\t0\tG\tF\tRG:Z:ghost",
        ],
    );
    fs::write(&rg1, b"existing").unwrap();

    let base = Options {
        mode: Mode::ReadGroup,
        output_prefix: &prefix,
        unaccounted: None,
        unaccounted_header: None,
        format: Format::Bam,
        maximum_outputs: 100,
        zero_pad: 0,
        reference: None,
        additional_threads: 0,
        program: None,
    };
    assert!(run(&input, base).is_err());
    assert_eq!(fs::read(&rg1).unwrap(), b"existing");
    assert!(!directory.path().join("sample.ghost.bam").exists());

    let summary = run(
        &input,
        Options {
            unaccounted: Some(&unaccounted),
            ..base
        },
    )
    .unwrap();
    assert_eq!(summary.records, 3);
    assert_eq!(summary.outputs.len(), 2);
    assert_eq!(summary.outputs[0].records, 1);
    assert_eq!(summary.outputs[1].label, "unaccounted");
    assert_eq!(summary.outputs[1].records, 2);
    assert_eq!(
        snapshot(&unaccounted).1,
        [
            (b"missing".to_vec(), Vec::new()),
            (b"unknown".to_vec(), b"ghost".to_vec()),
        ]
    );
    assert!(!directory.path().join("sample.ghost.bam").exists());
}

#[test]
fn sam_and_cram_outputs_keep_the_same_partition_contract() {
    let directory = tempfile::tempdir().unwrap();
    let input = fixture("read-group/tworg.bam");
    let sam_prefix = directory.path().join("sam");
    let base = Options {
        mode: Mode::ReadGroup,
        output_prefix: &sam_prefix,
        unaccounted: None,
        unaccounted_header: None,
        format: Format::Sam,
        maximum_outputs: 100,
        zero_pad: 0,
        reference: None,
        additional_threads: 0,
        program: None,
    };
    let sam_summary = run(&input, base).unwrap();
    assert_eq!(
        sam_summary
            .outputs
            .iter()
            .map(|output| output.records)
            .sum::<u64>(),
        9
    );
    for (label, expected) in [("rg1", 8), ("rg2", 1)] {
        let source = fs::read_to_string(directory.path().join(format!("sam.{label}.sam"))).unwrap();
        assert_eq!(
            source.lines().filter(|line| !line.starts_with('@')).count(),
            expected
        );
        assert_eq!(
            source
                .lines()
                .filter(|line| line.starts_with("@RG"))
                .count(),
            1
        );
        assert!(source.contains(&format!("@RG\tID:{label}")));
    }

    let cram_prefix = directory.path().join("cram");
    let without_reference = run(
        &input,
        Options {
            output_prefix: &cram_prefix,
            format: Format::Cram,
            ..base
        },
    );
    assert!(without_reference.is_err());
    assert!(!directory.path().join("cram.rg1.cram").exists());

    let reference = directory.path().join("reference.fa");
    fs::write(
        &reference,
        format!(">chr1\n{}\n>chr2\n{}\n", "A".repeat(1000), "A".repeat(1000)),
    )
    .unwrap();
    fs::write(
        directory.path().join("reference.fa.fai"),
        b"chr1\t1000\t6\t1000\t1001\nchr2\t1000\t1013\t1000\t1001\n",
    )
    .unwrap();
    let cram_summary = run(
        &input,
        Options {
            output_prefix: &cram_prefix,
            format: Format::Cram,
            reference: Some(&reference),
            ..base
        },
    )
    .unwrap();
    assert_eq!(
        cram_summary
            .outputs
            .iter()
            .map(|output| output.records)
            .sum::<u64>(),
        9
    );
    for (label, expected) in [("rg1", 8), ("rg2", 1)] {
        let path = directory.path().join(format!("cram.{label}.cram"));
        let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
            .args(["view", "--no-pg", "--reference"])
            .arg(&reference)
            .arg(&path)
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
                .split(|byte| *byte == b'\n')
                .filter(|line| !line.is_empty())
                .count(),
            expected
        );
    }
}

#[test]
fn unaccounted_header_is_isolated_and_dictionary_checked() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.sam");
    let replacement = directory.path().join("replacement.sam");
    let incompatible = directory.path().join("incompatible.sam");
    let prefix = directory.path().join("sample");
    let unaccounted = directory.path().join("other.bam");
    write_sam(
        &input,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@RG\tID:rg1\n@CO\toriginal",
        &["missing\t0\tchr1\t2\t60\t1M\t*\t0\t0\tC\tF"],
    );
    fs::write(
        &replacement,
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@CO\treplacement\n",
    )
    .unwrap();
    fs::write(
        &incompatible,
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:99\n@CO\twrong\n",
    )
    .unwrap();

    let base = Options {
        mode: Mode::ReadGroup,
        output_prefix: &prefix,
        unaccounted: Some(&unaccounted),
        unaccounted_header: Some(&incompatible),
        format: Format::Bam,
        maximum_outputs: 100,
        zero_pad: 0,
        reference: None,
        additional_threads: 0,
        program: None,
    };
    assert!(run(&input, base).is_err());
    assert!(!unaccounted.exists());

    let summary = run(
        &input,
        Options {
            unaccounted_header: Some(&replacement),
            ..base
        },
    )
    .unwrap();
    assert_eq!(
        summary
            .outputs
            .iter()
            .map(|output| output.records)
            .sum::<u64>(),
        1
    );
    let (header, records) = snapshot(&unaccounted);
    assert!(header.read_groups().is_empty());
    assert_eq!(header.comments(), &[b"replacement".as_slice()]);
    assert_eq!(records, [(b"missing".to_vec(), Vec::new())]);
    let (primary_header, primary_records) = snapshot(&directory.path().join("sample.rg1.bam"));
    assert_eq!(primary_header.comments(), &[b"original".as_slice()]);
    assert!(primary_records.is_empty());
}

#[test]
fn explicit_tags_encode_paths_pad_integers_and_enforce_the_output_cap() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("tags.sam");
    let prefix = directory.path().join("tags");
    let unaccounted = directory.path().join("overflow.bam");
    write_sam(
        &input,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@RG\tID:rg1",
        &[
            "slash\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tF\tCB:Z:a/b",
            "percent\t0\tchr1\t2\t60\t1M\t*\t0\t0\tC\tF\tCB:Z:a%b",
            "missing\t0\tchr1\t3\t60\t1M\t*\t0\t0\tG\tF",
        ],
    );
    let summary = run(
        &input,
        Options {
            mode: Mode::Tag(*b"CB"),
            output_prefix: &prefix,
            unaccounted: Some(&unaccounted),
            unaccounted_header: None,
            format: Format::Bam,
            maximum_outputs: 1,
            zero_pad: 0,
            reference: None,
            additional_threads: 0,
            program: None,
        },
    )
    .unwrap();
    assert_eq!(summary.records, 3);
    assert_eq!(
        summary
            .outputs
            .iter()
            .map(|output| output.records)
            .sum::<u64>(),
        3
    );
    let encoded = directory.path().join("tags.a%2Fb.bam");
    assert!(encoded.exists());
    assert!(!directory.path().join("tags.a").exists());
    assert_eq!(snapshot(&encoded).1, [(b"slash".to_vec(), Vec::new())]);
    assert_eq!(
        snapshot(&unaccounted)
            .1
            .into_iter()
            .map(|(name, _)| name)
            .collect::<Vec<_>>(),
        [b"percent".to_vec(), b"missing".to_vec()]
    );
    assert!(!directory.path().join("tags.a%25b.bam").exists());
    assert_eq!(snapshot(&encoded).0.read_groups().len(), 1);

    let integer = directory.path().join("integer.sam");
    let integer_prefix = directory.path().join("integer");
    write_sam(
        &integer,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100",
        &["seven\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tF\tNM:i:7"],
    );
    run(
        &integer,
        Options {
            mode: Mode::Tag(*b"NM"),
            output_prefix: &integer_prefix,
            unaccounted: None,
            unaccounted_header: None,
            format: Format::Bam,
            maximum_outputs: 1,
            zero_pad: 3,
            reference: None,
            additional_threads: 0,
            program: None,
        },
    )
    .unwrap();
    assert!(directory.path().join("integer.007.bam").exists());
}

#[test]
fn outputs_cannot_alias_primary_or_auxiliary_inputs() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("sample.rg1.bam");
    fs::copy(fixture("read-group/tworg.bam"), &input).unwrap();
    let original = fs::read(&input).unwrap();
    assert!(
        run(
            &input,
            Options {
                mode: Mode::ReadGroup,
                output_prefix: &directory.path().join("sample"),
                unaccounted: None,
                unaccounted_header: None,
                format: Format::Bam,
                maximum_outputs: 100,
                zero_pad: 0,
                reference: None,
                additional_threads: 0,
                program: None,
            },
        )
        .is_err()
    );
    assert_eq!(fs::read(&input).unwrap(), original);

    let sam = directory.path().join("missing.sam");
    let header = directory.path().join("header.sam");
    write_sam(
        &sam,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@RG\tID:rg1",
        &["missing\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tF"],
    );
    fs::write(&header, b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n").unwrap();
    let original_header = fs::read(&header).unwrap();
    assert!(
        run(
            &sam,
            Options {
                mode: Mode::ReadGroup,
                output_prefix: &directory.path().join("unused"),
                unaccounted: Some(&header),
                unaccounted_header: Some(&header),
                format: Format::Bam,
                maximum_outputs: 100,
                zero_pad: 0,
                reference: None,
                additional_threads: 0,
                program: None,
            },
        )
        .is_err()
    );
    assert_eq!(fs::read(&header).unwrap(), original_header);
}

#[test]
fn explicit_read_group_tag_synthesizes_undeclared_headers() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.sam");
    let prefix = directory.path().join("sample");
    write_sam(
        &input,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@RG\tID:rg1\tSM:declared",
        &[
            "known\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tF\tRG:Z:rg1",
            "unknown\t0\tchr1\t2\t60\t1M\t*\t0\t0\tC\tF\tRG:Z:ghost",
        ],
    );
    let summary = run(
        &input,
        Options {
            mode: Mode::Tag(*b"RG"),
            output_prefix: &prefix,
            unaccounted: None,
            unaccounted_header: None,
            format: Format::Bam,
            maximum_outputs: 100,
            zero_pad: 0,
            reference: None,
            additional_threads: 0,
            program: None,
        },
    )
    .unwrap();
    assert_eq!(
        summary
            .outputs
            .iter()
            .map(|output| output.records)
            .sum::<u64>(),
        2
    );
    let (header, records) = snapshot(&directory.path().join("sample.ghost.bam"));
    assert_eq!(
        header
            .read_groups()
            .keys()
            .map(|id| id.as_slice())
            .collect::<Vec<_>>(),
        [b"ghost".as_slice()]
    );
    assert!(
        header.read_groups()[b"ghost".as_slice()]
            .other_fields()
            .is_empty()
    );
    assert_eq!(records, [(b"unknown".to_vec(), b"ghost".to_vec())]);
}

#[test]
fn seeded_parts_are_a_reproducible_exact_cover() {
    let directory = tempfile::tempdir().unwrap();
    let input = fixture("parts/diverse.bam");
    let first = directory.path().join("first");
    let second = directory.path().join("second");
    let base = Options {
        mode: Mode::Parts {
            count: 3,
            seed: 7,
            skip_unmapped: false,
        },
        output_prefix: &first,
        unaccounted: None,
        unaccounted_header: None,
        format: Format::Bam,
        maximum_outputs: 100,
        zero_pad: 2,
        reference: None,
        additional_threads: 0,
        program: None,
    };
    let first_summary = run(&input, base).unwrap();
    let second_summary = run(
        &input,
        Options {
            output_prefix: &second,
            ..base
        },
    )
    .unwrap();
    assert_eq!(first_summary.records, 9);
    assert_eq!(first_summary.outputs.len(), 3);
    assert_eq!(
        first_summary
            .outputs
            .iter()
            .map(|output| output.records)
            .sum::<u64>(),
        9
    );
    assert_eq!(
        second_summary
            .outputs
            .iter()
            .map(|output| output.records)
            .sum::<u64>(),
        9
    );

    let input_names = names(&input).into_iter().collect::<HashSet<_>>();
    let mut seen = HashSet::new();
    for part in 0..3 {
        let label = format!("{part:02}");
        let first_names = names(&directory.path().join(format!("first.{label}.bam")));
        let second_names = names(&directory.path().join(format!("second.{label}.bam")));
        assert_eq!(first_names, second_names);
        for name in first_names {
            assert!(seen.insert(name));
        }
    }
    assert_eq!(seen, input_names);
}

#[test]
fn parts_validate_limits_and_skip_only_unmapped_records() {
    let directory = tempfile::tempdir().unwrap();
    let input = fixture("parts/diverse.bam");
    let prefix = directory.path().join("parts");
    let base = Options {
        mode: Mode::Parts {
            count: 0,
            seed: 0,
            skip_unmapped: false,
        },
        output_prefix: &prefix,
        unaccounted: None,
        unaccounted_header: None,
        format: Format::Bam,
        maximum_outputs: 2,
        zero_pad: 0,
        reference: None,
        additional_threads: 0,
        program: None,
    };
    assert!(run(&input, base).is_err());
    assert!(!directory.path().join("parts.0.bam").exists());
    assert!(
        run(
            &input,
            Options {
                mode: Mode::Parts {
                    count: 3,
                    seed: 0,
                    skip_unmapped: false,
                },
                ..base
            },
        )
        .is_err()
    );
    assert!(!directory.path().join("parts.1.bam").exists());

    let summary = run(
        &input,
        Options {
            mode: Mode::Parts {
                count: 2,
                seed: 0,
                skip_unmapped: true,
            },
            ..base
        },
    )
    .unwrap();
    assert_eq!(summary.records, 9);
    assert_eq!(summary.skipped, 1);
    assert_eq!(
        summary
            .outputs
            .iter()
            .map(|output| output.records)
            .sum::<u64>(),
        8
    );
    let retained = [
        names(&directory.path().join("parts.0.bam")),
        names(&directory.path().join("parts.1.bam")),
    ]
    .concat()
    .into_iter()
    .collect::<HashSet<_>>();
    assert!(!retained.contains(b"r5".as_slice()));
    assert_eq!(retained.len(), 8);
}

#[test]
fn gene_mode_matches_the_retained_rseqc_record_bodies() {
    let directory = tempfile::tempdir().unwrap();
    let input = fixture("genes/reads.bam");
    let bed = fixture("genes/genes.strict.bed12");
    let prefix = directory.path().join("genes");
    let summary = run(
        &input,
        Options {
            mode: Mode::Genes(&bed),
            output_prefix: &prefix,
            unaccounted: None,
            unaccounted_header: None,
            format: Format::Bam,
            maximum_outputs: 100,
            zero_pad: 0,
            reference: None,
            additional_threads: 0,
            program: None,
        },
    )
    .unwrap();
    assert_eq!(summary.records, 9);
    assert_eq!(summary.skipped, 0);
    assert_eq!(
        summary
            .outputs
            .iter()
            .map(|output| (output.label.as_str(), output.records))
            .collect::<Vec<_>>(),
        [("in", 5), ("ex", 2), ("junk", 2)]
    );
    for label in ["in", "ex", "junk"] {
        assert_eq!(
            view_body(&directory.path().join(format!("genes.{label}.bam"))),
            fs::read(fixture(&format!("genes/{label}.sam"))).unwrap()
        );
    }
}

#[test]
fn gene_mode_uses_leftmost_start_and_preserves_targets_on_bed_failure() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("reads.sam");
    let bed = directory.path().join("genes.bed12");
    let prefix = directory.path().join("genes");
    write_sam(
        &input,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:500",
        &[
            "overlap\t0\tchr1\t51\t60\t10M50N\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF",
            "inside\t0\tchr1\t101\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF",
            "qcfail\t512\tchr1\t101\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF",
            "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*",
        ],
    );
    fs::write(
        &bed,
        b"chr1\t100\t200\tgene\t0\t+\t100\t200\t0\t1\t100,\t0,\n",
    )
    .unwrap();
    run(
        &input,
        Options {
            mode: Mode::Genes(&bed),
            output_prefix: &prefix,
            unaccounted: None,
            unaccounted_header: None,
            format: Format::Sam,
            maximum_outputs: 100,
            zero_pad: 0,
            reference: None,
            additional_threads: 0,
            program: None,
        },
    )
    .unwrap();
    assert_eq!(
        body_names(&view_body(&directory.path().join("genes.in.sam"))),
        [b"inside".to_vec()]
    );
    assert_eq!(
        body_names(&view_body(&directory.path().join("genes.ex.sam"))),
        [b"overlap".to_vec()]
    );
    assert_eq!(
        body_names(&view_body(&directory.path().join("genes.junk.sam"))),
        [b"qcfail".to_vec(), b"unmapped".to_vec()]
    );

    for label in ["in", "ex", "junk"] {
        fs::write(directory.path().join(format!("prior.{label}.bam")), label).unwrap();
    }
    let prior = directory.path().join("prior");
    assert!(
        run(
            &input,
            Options {
                mode: Mode::Genes(&fixture("genes/genes.bed12")),
                output_prefix: &prior,
                unaccounted: None,
                unaccounted_header: None,
                format: Format::Bam,
                maximum_outputs: 100,
                zero_pad: 0,
                reference: None,
                additional_threads: 0,
                program: None,
            },
        )
        .is_err()
    );
    for label in ["in", "ex", "junk"] {
        assert_eq!(
            fs::read(directory.path().join(format!("prior.{label}.bam"))).unwrap(),
            label.as_bytes()
        );
    }
}
