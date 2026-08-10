use std::fs;
use std::path::Path;
use std::process::{Command, Output};

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

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn samtools() -> Command {
    Command::new("samtools")
}

fn stable_sam(bytes: &[u8]) -> Vec<&[u8]> {
    bytes
        .split(|byte| *byte == b'\n')
        .filter(|line| !line.is_empty() && !line.starts_with(b"@PG\t"))
        .collect()
}

fn decode(path: &Path) -> Output {
    run({
        let mut command = samtools();
        command.args(["view", "-h"]).arg(path);
        command
    })
}

fn compare(input: &Path, options: &[&str]) {
    let ours = run({
        let mut command = binary();
        command
            .arg("addreplacerg")
            .args(options)
            .arg("--no-PG")
            .arg(input);
        command
    });
    let oracle = run({
        let mut command = samtools();
        command
            .arg("addreplacerg")
            .args(options)
            .arg("--no-PG")
            .arg(input);
        command
    });
    assert_eq!(stable_sam(&ours.stdout), stable_sam(&oracle.stdout));
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn addreplacerg_matches_samtools_1_24_source_and_mode_matrix() {
    let version = run({
        let mut command = samtools();
        command.arg("--version");
        command
    });
    assert!(version.stdout.starts_with(b"samtools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.sam");
    fs::write(
        &input,
        b"@HD\tVN:1.6\tSO:unknown\n\
@SQ\tSN:chr1\tLN:20\n\
@RG\tID:old\tSM:before\n\
one\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:old\tNM:i:0\n\
two\t4\t*\t0\t0\t*\t*\t0\t0\tN\t!\n",
    )
    .unwrap();

    for options in [
        vec!["-R", "old"],
        vec![],
        vec!["-r", "ID:new\\tSM:after"],
        vec!["-r", "ID:new", "-r", "SM:after"],
        vec!["-r", "ID:new\\tSM:after", "-m", "orphan_only"],
    ] {
        compare(&input, &options);
    }

    let replaced = directory.path().join("replaced.sam");
    run({
        let mut command = samtools();
        command
            .args([
                "addreplacerg",
                "-w",
                "-r",
                "ID:old\\tSM:replacement",
                "--no-PG",
                "-o",
            ])
            .arg(&replaced)
            .arg(&input);
        command
    });
    compare(&input, &["-w", "-r", "ID:old\\tSM:replacement"]);

    let numeric = directory.path().join("numeric.sam");
    fs::write(
        &numeric,
        b"@HD\tVN:1.6\n@RG\tID:old\nread\t4\t*\t0\t0\t*\t*\t0\t0\tN\t!\tRG:i:7\n",
    )
    .unwrap();
    compare(&numeric, &["-R", "old"]);
    compare(&numeric, &["-R", "old", "-m", "orphan_only"]);
}

#[test]
#[ignore = "release oracle: requires samtools 1.24"]
fn addreplacerg_bam_and_cram_paths_match_samtools_1_24() {
    let version = run({
        let mut command = samtools();
        command.arg("--version");
        command
    });
    assert!(version.stdout.starts_with(b"samtools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    fs::write(&reference, b">chr1\nACGTACGTACGTACGTACGT\n").unwrap();
    run({
        let mut command = samtools();
        command.arg("faidx").arg(&reference);
        command
    });
    let sam = directory.path().join("input.sam");
    fs::write(
        &sam,
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:20\n@RG\tID:old\tSM:before\nread\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:old\n",
    )
    .unwrap();
    let bam = directory.path().join("input.bam");
    let cram = directory.path().join("input.cram");
    run({
        let mut command = samtools();
        command.args(["view", "-b", "-o"]).arg(&bam).arg(&sam);
        command
    });
    run({
        let mut command = samtools();
        command
            .args(["view", "-C", "-T"])
            .arg(&reference)
            .args(["-o"])
            .arg(&cram)
            .arg(&sam);
        command
    });

    for input in [&bam, &cram] {
        let ours = directory.path().join(format!(
            "{}.ours.bam",
            input.file_stem().unwrap().to_string_lossy()
        ));
        let oracle = directory.path().join(format!(
            "{}.oracle.bam",
            input.file_stem().unwrap().to_string_lossy()
        ));
        run({
            let mut command = binary();
            command.args(["addreplacerg", "-r", "ID:new\\tSM:after", "--no-PG"]);
            if input == &cram {
                command.arg("--reference").arg(&reference);
            }
            command.args(["-O", "bam", "-o"]).arg(&ours).arg(input);
            command
        });
        run({
            let mut command = samtools();
            command.args(["addreplacerg", "-r", "ID:new\\tSM:after", "--no-PG"]);
            if input == &cram {
                command.arg("--reference").arg(&reference);
            }
            command.args(["-O", "bam", "-o"]).arg(&oracle).arg(input);
            command
        });
        assert_eq!(
            stable_sam(&decode(&ours).stdout),
            stable_sam(&decode(&oracle).stdout),
            "{}",
            input.display()
        );
    }
}
