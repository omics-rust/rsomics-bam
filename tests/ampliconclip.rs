use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

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

fn decoded_with_samtools(path: &Path) -> Vec<u8> {
    run({
        let mut command = Command::new("samtools");
        command.args(["view", "--no-PG"]).arg(path);
        command
    })
    .stdout
}

fn header_with_samtools(path: &Path) -> Vec<u8> {
    run({
        let mut command = Command::new("samtools");
        command.args(["view", "--no-PG", "-H"]).arg(path);
        command
    })
    .stdout
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden/ampliconclip")
        .join(name)
}

fn make_input(directory: &Path) -> PathBuf {
    let input = directory.join("input.bam");
    run({
        let mut command = binary();
        command
            .args(["sort", "--no-pg", "-o"])
            .arg(&input)
            .arg(fixture("reads.sam"));
        command
    });
    input
}

const MODES: &[(&str, &[&str])] = &[
    ("soft_default", &[]),
    ("hard", &["--hard-clip"]),
    ("both_soft", &["--both-ends"]),
    ("both_hard", &["--both-ends", "--hard-clip"]),
    ("strand_soft", &["--strand"]),
    ("strand_hard", &["--strand", "--hard-clip"]),
    ("strand_both", &["--strand", "--both-ends"]),
    ("keep_tag", &["--keep-tag"]),
    ("fail", &["--fail"]),
    ("clipped", &["--clipped"]),
    ("no_excluded", &["--no-excluded"]),
    ("tol0", &["--tolerance", "0"]),
    ("tol20", &["--tolerance", "20"]),
    ("filter_len", &["--filter-len", "50"]),
    ("fail_len", &["--fail-len", "50"]),
    ("unmap_len", &["--unmap-len", "90"]),
    (
        "combo",
        &[
            "--both-ends",
            "--hard-clip",
            "--strand",
            "--fail-len",
            "30",
            "--keep-tag",
            "--unmap-len",
            "40",
        ],
    ),
];

#[test]
fn committed_alignment_matrix_matches_samtools() {
    let directory = tempfile::tempdir().unwrap();
    let input = make_input(directory.path());
    for (label, options) in MODES {
        let output = directory.path().join(format!("{label}.bam"));
        run({
            let mut command = binary();
            command
                .args(["ampliconclip", "--no-pg", "-b"])
                .arg(fixture("primers.bed"))
                .args(*options)
                .args(["-o"])
                .arg(&output)
                .arg(&input);
            command
        });
        let decoded = run({
            let mut command = binary();
            command.args(["view", "--no-pg"]).arg(&output);
            command
        });
        let expected = fs::read(fixture(&format!("expected/{label}.sam"))).unwrap();
        assert_eq!(decoded.stdout, expected, "mode={label}");
    }
}

#[test]
fn json_counts_and_transactional_outputs_share_the_product_contract() {
    let directory = tempfile::tempdir().unwrap();
    let input = make_input(directory.path());
    let output = directory.path().join("clipped.bam");
    let stats = directory.path().join("stats.txt");
    let counts = directory.path().join("counts.bedgraph");
    let result = run({
        let mut command = binary();
        command
            .arg("--json")
            .args(["ampliconclip", "--no-pg", "-b"])
            .arg(fixture("primers.bed"))
            .args(["-o"])
            .arg(&output)
            .args(["-f"])
            .arg(&stats)
            .arg("--primer-counts")
            .arg(&counts)
            .arg(&input);
        command
    });
    let json: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(json["result"]["command"], "ampliconclip");
    assert_eq!(json["result"]["run"]["summary"]["total"], 9);
    assert_eq!(json["result"]["run"]["summary"]["written"], 9);
    assert!(
        fs::read_to_string(stats)
            .unwrap()
            .contains("TOTAL CLIPPED: 7")
    );
    assert_eq!(fs::read_to_string(counts).unwrap().lines().count(), 5);
}

#[test]
fn failures_preserve_named_outputs_and_reject_aliases() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("output.bam");
    fs::write(&output, b"sentinel").unwrap();
    let failed = binary()
        .args(["ampliconclip", "-b"])
        .arg(fixture("primers.bed"))
        .args(["-o"])
        .arg(&output)
        .arg(fixture("reads.sam"))
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert_eq!(fs::read(&output).unwrap(), b"sentinel");

    let input = make_input(directory.path());
    let alias = binary()
        .args(["ampliconclip", "-b"])
        .arg(fixture("primers.bed"))
        .args(["-o"])
        .arg(&input)
        .arg(&input)
        .output()
        .unwrap();
    assert!(!alias.status.success());
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn live_alignment_matrix_matches_samtools_1_24() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let input = make_input(directory.path());
    for (label, options) in MODES {
        let ours = directory.path().join(format!("ours-{label}.bam"));
        let oracle = directory.path().join(format!("oracle-{label}.bam"));
        run({
            let mut command = binary();
            command
                .args(["ampliconclip", "--no-pg", "-b"])
                .arg(fixture("primers.bed"))
                .args(*options)
                .arg("-o")
                .arg(&ours)
                .arg(&input);
            command
        });
        run({
            let mut command = Command::new("samtools");
            command
                .args(["ampliconclip", "--no-PG", "-b"])
                .arg(fixture("primers.bed"))
                .args(*options)
                .arg("-o")
                .arg(&oracle)
                .arg(&input);
            command
        });
        assert_eq!(
            decoded_with_samtools(&ours),
            decoded_with_samtools(&oracle),
            "mode={label}"
        );
        assert_eq!(
            header_with_samtools(&ours),
            header_with_samtools(&oracle),
            "header mode={label}"
        );
    }
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn live_auxiliary_outputs_match_samtools_1_24() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let input = make_input(directory.path());
    let ours = directory.path().join("ours.bam");
    let oracle = directory.path().join("oracle.bam");
    let ours_rejects = directory.path().join("ours-rejects.bam");
    let oracle_rejects = directory.path().join("oracle-rejects.bam");
    let ours_stats = directory.path().join("ours-stats.txt");
    let oracle_stats = directory.path().join("oracle-stats.txt");
    let ours_counts = directory.path().join("ours-counts.txt");
    let oracle_counts = directory.path().join("oracle-counts.txt");
    let options = [
        "-u",
        "--original",
        "--filter-len",
        "50",
        "--strand",
        "-@",
        "2",
    ];

    run({
        let mut command = binary();
        command
            .args(["ampliconclip", "--no-pg", "-b"])
            .arg(fixture("primers.bed"))
            .args(options)
            .arg("--rejects-file")
            .arg(&ours_rejects)
            .arg("--primer-counts")
            .arg(&ours_counts)
            .arg("-f")
            .arg(&ours_stats)
            .arg("-o")
            .arg(&ours)
            .arg(&input);
        command
    });
    run({
        let mut command = Command::new("samtools");
        command
            .args(["ampliconclip", "--no-PG", "-b"])
            .arg(fixture("primers.bed"))
            .args(options)
            .arg("--rejects-file")
            .arg(&oracle_rejects)
            .arg("--primer-counts")
            .arg(&oracle_counts)
            .arg("-f")
            .arg(&oracle_stats)
            .arg("-o")
            .arg(&oracle)
            .arg(&input);
        command
    });

    assert_eq!(decoded_with_samtools(&ours), decoded_with_samtools(&oracle));
    assert_eq!(
        decoded_with_samtools(&ours_rejects),
        decoded_with_samtools(&oracle_rejects)
    );
    assert_eq!(
        fs::read(&ours_counts).unwrap(),
        fs::read(&oracle_counts).unwrap()
    );
    let normalize = |path: &Path| {
        fs::read_to_string(path)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with("COMMAND:"))
            .map(str::to_owned)
            .collect::<Vec<_>>()
    };
    assert_eq!(normalize(&ours_stats), normalize(&oracle_stats));
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn live_cigar_and_auxiliary_types_match_samtools_1_24() {
    assert_samtools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let sam = directory.path().join("types.sam");
    let input = directory.path().join("types.bam");
    let ours = directory.path().join("ours-types.bam");
    let oracle = directory.path().join("oracle-types.bam");
    fs::write(
        &sam,
        b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n@SQ\tSN:chr2\tLN:2000\ntypes\t0\tchr1\t101\t60\t5S10=2X3I4D5N1P20M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\tXA:A:C\tXZ:Z:text\tXH:H:0A0B\tXc:B:c,-1,2\tXC:B:C,1,2\tXs:B:s,-3,4\tXS:B:S,3,4\tXi:B:i,-5,6\tXI:B:I,5,6\tXf:B:f,1.5,2.5\tNM:i:6\tMD:Z:10AC0^AAAA20\n",
    )
    .unwrap();
    run({
        let mut command = Command::new("samtools");
        command.args(["view", "-b", "-o"]).arg(&input).arg(&sam);
        command
    });
    for (mut command, output) in [(binary(), &ours), (Command::new("samtools"), &oracle)] {
        run({
            command
                .args(["ampliconclip", "--no-PG", "--hard-clip", "--keep-tag", "-b"])
                .arg(fixture("primers.bed"))
                .arg("-o")
                .arg(output)
                .arg(&input);
            command
        });
    }
    assert_eq!(decoded_with_samtools(&ours), decoded_with_samtools(&oracle));
}
