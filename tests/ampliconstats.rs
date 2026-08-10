use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden/ampliconstats")
        .join(name)
}

fn make_multireference_input(directory: &Path) -> (PathBuf, PathBuf) {
    let bed = directory.join("multi.bed");
    let sam = directory.join("multi.sam");
    let bam = directory.join("multi.bam");
    fs::write(
        &bed,
        b"chr1\t10\t20\ta\t0\t+\nchr1\t12\t22\ta-alt\t0\t+\nchr1\t80\t90\tb\t0\t-\nchr2\t10\t20\tc\t0\t+\nchr2\t70\t80\td\t0\t-\n",
    )
    .unwrap();
    fs::write(
        &sam,
        b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:200\n@SQ\tSN:chr2\tLN:150\n@RG\tID:rg1\tSM:multi\nq1\t99\tchr1\t21\t60\t20M\t=\t62\t61\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1\nq1\t147\tchr1\t62\t60\t20M\t=\t21\t-61\tTTTTTTTTTTTTTTTTTTTT\tIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1\nq2\t99\tchr2\t21\t60\t20M\t=\t52\t51\tCCCCCCCCCCCCCCCCCCCC\tIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1\nq2\t147\tchr2\t52\t60\t20M\t=\t21\t-51\tGGGGGGGGGGGGGGGGGGGG\tIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1\n",
    )
    .unwrap();
    run({
        let mut command = binary();
        command
            .args(["view", "--no-pg", "-b", "-o"])
            .arg(&bam)
            .arg(&sam);
        command
    });
    (bed, bam)
}

fn make_unsorted_input(directory: &Path) -> PathBuf {
    let decoded = run({
        let mut command = binary();
        command
            .args(["view", "--no-pg", "-h"])
            .arg(fixture("amplicons.bam"));
        command
    });
    let text = String::from_utf8(decoded.stdout).unwrap();
    let mut headers = Vec::new();
    let mut records = Vec::new();
    for line in text.lines() {
        if line.starts_with('@') {
            headers.push(line);
        } else {
            records.push(line);
        }
    }
    let last = records.len() - 1;
    records.swap(0, last);
    let sam = directory.join("unsorted.sam");
    let bam = directory.join("unsorted.bam");
    fs::write(
        &sam,
        headers
            .into_iter()
            .chain(records)
            .collect::<Vec<_>>()
            .join("\n")
            + "\n",
    )
    .unwrap();
    run({
        let mut command = binary();
        command
            .args(["view", "--no-pg", "-b", "-o"])
            .arg(&bam)
            .arg(&sam);
        command
    });
    bam
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

fn assert_samtools_1_24() {
    let output = run({
        let mut command = Command::new("samtools");
        command.arg("--version");
        command
    });
    assert_eq!(
        String::from_utf8_lossy(&output.stdout).lines().next(),
        Some("samtools 1.24")
    );
}

fn normalized(text: &str) -> Vec<&str> {
    text.lines()
        .filter(|line| {
            !line.starts_with("SS\tSamtools version:") && !line.starts_with("SS\tCommand line:")
        })
        .collect()
}

#[test]
fn default_output_matches_the_committed_samtools_oracle() {
    let result = run({
        let mut command = binary();
        command
            .arg("ampliconstats")
            .arg(fixture("primers.bed"))
            .arg(fixture("amplicons.bam"));
        command
    });
    let actual = String::from_utf8(result.stdout).unwrap();
    let expected = fs::read_to_string(fixture("expected_ampliconstats.txt")).unwrap();
    assert_eq!(normalized(&actual), normalized(&expected));
}

#[test]
fn sample_names_named_output_and_json_use_the_product_contract() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("amplicons.txt");
    let result = run({
        let mut command = binary();
        command
            .arg("--json")
            .args(["ampliconstats", "-s", "-o"])
            .arg(&output)
            .arg(fixture("primers.bed"))
            .arg(fixture("amplicons.bam"));
        command
    });
    let json: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(json["result"]["command"], "ampliconstats");
    assert_eq!(json["result"]["summary"]["files"], 1);
    assert_eq!(json["result"]["summary"]["amplicons"], 2);
    assert_eq!(json["result"]["summary"]["records"], 50);
    let text = fs::read_to_string(output).unwrap();
    assert!(text.contains("FSS\tsample1\tchr1\traw total sequences:\t50"));
}

#[test]
fn failures_preserve_named_output_and_reject_aliases() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("stats.txt");
    fs::write(&output, b"sentinel").unwrap();
    let failed = binary()
        .args(["ampliconstats", "-o"])
        .arg(&output)
        .arg(fixture("primers.bed"))
        .arg(fixture("expected_ampliconstats.txt"))
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert_eq!(fs::read(&output).unwrap(), b"sentinel");

    let alias = binary()
        .args(["ampliconstats", "-o"])
        .arg(fixture("amplicons.bam"))
        .arg(fixture("primers.bed"))
        .arg(fixture("amplicons.bam"))
        .output()
        .unwrap();
    assert!(!alias.status.success());

    fs::write(&output, b"sentinel").unwrap();
    let oversized = binary()
        .args(["ampliconstats", "--max-amplicon-length", "20", "-o"])
        .arg(&output)
        .arg(fixture("primers.bed"))
        .arg(fixture("amplicons.bam"))
        .output()
        .unwrap();
    assert!(!oversized.status.success());
    assert_eq!(fs::read(&output).unwrap(), b"sentinel");

    let unsorted = make_unsorted_input(directory.path());
    let bad_order = binary()
        .args(["ampliconstats", "-o"])
        .arg(&output)
        .arg(fixture("primers.bed"))
        .arg(unsorted)
        .output()
        .unwrap();
    assert!(!bad_order.status.success());
    assert_eq!(fs::read(&output).unwrap(), b"sentinel");
}

#[test]
fn multiple_references_and_alternative_primers_use_one_model() {
    let directory = tempfile::tempdir().unwrap();
    let (bed, bam) = make_multireference_input(directory.path());
    let result = run({
        let mut command = binary();
        command.arg("ampliconstats").arg(&bed).arg(&bam);
        command
    });
    let text = String::from_utf8(result.stdout).unwrap();
    assert!(text.contains("AMPLICON\tchr1\t1\t11-20,13-22\t81-90"));
    assert!(text.contains("AMPLICON\tchr2\t2\t11-20\t71-80"));
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn live_option_matrix_matches_samtools_1_24() {
    assert_samtools_1_24();
    let cases: &[(&str, &[&str])] = &[
        ("default", &[]),
        ("required", &["-f", "PAIRED"]),
        ("filter", &["-F", "SECONDARY,QCFAIL"]),
        ("limits", &["-a", "20", "-l", "200"]),
        ("depths", &["-d", "1,2,3"]),
        ("margin", &["-m", "5"]),
        ("sample", &["-s"]),
        ("template", &["-t", "10", "-b", "5", "-c", "1"]),
        ("depth_bin", &["-D", "0"]),
        ("single_ref", &["-S"]),
        ("threads", &["-@", "2"]),
    ];
    for (label, options) in cases {
        let ours = run({
            let mut command = binary();
            command
                .arg("ampliconstats")
                .args(*options)
                .arg(fixture("primers.bed"))
                .arg(fixture("amplicons.bam"));
            command
        });
        let oracle = run({
            let mut command = Command::new("samtools");
            command
                .arg("ampliconstats")
                .args(*options)
                .arg(fixture("primers.bed"))
                .arg(fixture("amplicons.bam"));
            command
        });
        let ours = String::from_utf8(ours.stdout).unwrap();
        let oracle = String::from_utf8(oracle.stdout).unwrap();
        assert_eq!(normalized(&ours), normalized(&oracle), "case={label}");
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn live_multiple_input_output_matches_samtools_1_24() {
    assert_samtools_1_24();
    let input = fixture("amplicons.bam");
    let ours = run({
        let mut command = binary();
        command
            .args(["ampliconstats", "-s"])
            .arg(fixture("primers.bed"))
            .arg(&input)
            .arg(&input);
        command
    });
    let oracle = run({
        let mut command = Command::new("samtools");
        command
            .args(["ampliconstats", "-s"])
            .arg(fixture("primers.bed"))
            .arg(&input)
            .arg(&input);
        command
    });
    let ours = String::from_utf8(ours.stdout).unwrap();
    let oracle = String::from_utf8(oracle.stdout).unwrap();
    assert_eq!(normalized(&ours), normalized(&oracle));
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn live_multiple_references_and_alternative_primers_match_samtools_1_24() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let (bed, bam) = make_multireference_input(directory.path());
    let ours = run({
        let mut command = binary();
        command
            .args(["ampliconstats", "-c", "1"])
            .arg(&bed)
            .arg(&bam);
        command
    });
    let oracle = run({
        let mut command = Command::new("samtools");
        command
            .args(["ampliconstats", "-c", "1"])
            .arg(&bed)
            .arg(&bam);
        command
    });
    let ours = String::from_utf8(ours.stdout).unwrap();
    let oracle = String::from_utf8(oracle.stdout).unwrap();
    assert_eq!(normalized(&ours), normalized(&oracle));
}
