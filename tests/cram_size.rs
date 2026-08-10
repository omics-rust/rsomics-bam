use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam"))
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden/cram-size")
        .join(name)
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
fn all_text_modes_match_samtools_1_24() {
    for (arguments, expected) in [
        (&[][..], "default.txt"),
        (&["-v"][..], "verbose.txt"),
        (&["-e"][..], "encodings.txt"),
        (&["-vv"][..], "verbose.txt"),
        (&["-ee"][..], "encodings.txt"),
    ] {
        let actual = run({
            let mut command = binary();
            command
                .arg("cram-size")
                .args(arguments)
                .arg(fixture("input.cram"));
            command
        });
        assert_eq!(actual.stdout, fs::read(fixture(expected)).unwrap());
    }
}

#[test]
fn supported_versions_match_samtools_1_24() {
    for version in ["2.1", "3.0", "3.1"] {
        for (arguments, suffix) in [(&[][..], ""), (&["-e"][..], "-encodings")] {
            let actual = run({
                let mut command = binary();
                command
                    .arg("cram-size")
                    .args(arguments)
                    .arg(fixture(&format!("version-{version}.cram")));
                command
            });
            assert_eq!(
                actual.stdout,
                fs::read(fixture(&format!("version-{version}{suffix}.txt"))).unwrap(),
                "CRAM {version}"
            );
        }
    }
}

#[test]
fn standard_input_matches_file_input() {
    let mut child = binary()
        .args(["cram-size", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&fs::read(fixture("input.cram")).unwrap())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, fs::read(fixture("default.txt")).unwrap());
}

#[test]
fn named_output_and_json_use_product_contract() {
    let directory = tempfile::tempdir().unwrap();
    let output_path = directory.path().join("sizes.txt");
    let output = run({
        let mut command = binary();
        command
            .arg("--json")
            .args(["cram-size", "-o"])
            .arg(&output_path)
            .arg(fixture("input.cram"));
        command
    });
    assert_eq!(
        fs::read(output_path).unwrap(),
        fs::read(fixture("default.txt")).unwrap()
    );
    let json: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(json["result"]["command"], "cram-size");
    let report = &json["result"]["report"];
    assert_eq!(report["containers"], 2);
    assert_eq!(report["slices"], 6);
    assert_eq!(report["sequences"], 569);
    assert_eq!(report["bases"], 57_572);
    assert_eq!(report["file_size"], 60_965);
    assert_eq!(report["format_overhead_size"], 15_133);
    assert_eq!(report["blocks"].as_array().unwrap().len(), 30);
}

#[test]
fn invalid_inputs_fail_without_replacing_output() {
    let directory = tempfile::tempdir().unwrap();
    let output_path = directory.path().join("sizes.txt");
    fs::write(&output_path, b"sentinel").unwrap();

    let bam = binary()
        .args(["cram-size", "-o"])
        .arg(&output_path)
        .arg(fixture("../flagstat-small.bam"))
        .output()
        .unwrap();
    assert!(!bam.status.success());
    assert_eq!(fs::read(&output_path).unwrap(), b"sentinel");

    let mut corrupted = fs::read(fixture("input.cram")).unwrap();
    let middle = corrupted.len() / 2;
    corrupted[middle] ^= 0x80;
    let corrupt_path = directory.path().join("corrupt.cram");
    fs::write(&corrupt_path, corrupted).unwrap();
    let corrupt = binary()
        .args(["cram-size", "-o"])
        .arg(&output_path)
        .arg(&corrupt_path)
        .output()
        .unwrap();
    assert!(!corrupt.status.success());
    assert_eq!(fs::read(&output_path).unwrap(), b"sentinel");

    let mut unsupported = fs::read(fixture("input.cram")).unwrap();
    unsupported[4] = 4;
    let unsupported_path = directory.path().join("unsupported.cram");
    fs::write(&unsupported_path, unsupported).unwrap();
    let unsupported = binary()
        .args(["cram-size", "-o"])
        .arg(&output_path)
        .arg(&unsupported_path)
        .output()
        .unwrap();
    assert!(!unsupported.status.success());
    assert_eq!(fs::read(&output_path).unwrap(), b"sentinel");

    let mut truncated = fs::read(fixture("input.cram")).unwrap();
    truncated.pop();
    let truncated_path = directory.path().join("truncated.cram");
    fs::write(&truncated_path, truncated).unwrap();
    let truncated = binary()
        .args(["cram-size", "-o"])
        .arg(&output_path)
        .arg(&truncated_path)
        .output()
        .unwrap();
    assert!(!truncated.status.success());
    assert_eq!(fs::read(&output_path).unwrap(), b"sentinel");

    let alias = binary()
        .args(["cram-size", "-o"])
        .arg(fixture("input.cram"))
        .arg(fixture("input.cram"))
        .output()
        .unwrap();
    assert!(!alias.status.success());
}

#[test]
fn json_requires_named_compatibility_output() {
    let output = binary()
        .args(["--json", "cram-size"])
        .arg(fixture("input.cram"))
        .output()
        .unwrap();
    assert!(!output.status.success());
}
