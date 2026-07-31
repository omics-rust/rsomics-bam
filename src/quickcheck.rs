use std::io::Write;
use std::path::{Path, PathBuf};

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::hts_quickcheck;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Options {
    pub allow_no_targets: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Problem {
    Open,
    NotSequence,
    Header,
    NoTargets,
    EofCheck,
    MissingEof,
    Close,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct FileReport {
    pub path: PathBuf,
    pub problems: Vec<Problem>,
}

impl FileReport {
    #[must_use]
    pub fn is_ok(&self) -> bool {
        self.problems.is_empty()
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Report {
    pub files: Vec<FileReport>,
}

impl Report {
    #[must_use]
    pub fn is_ok(&self) -> bool {
        self.files.iter().all(FileReport::is_ok)
    }

    #[must_use]
    pub fn failed(&self) -> usize {
        self.files.iter().filter(|file| !file.is_ok()).count()
    }

    pub fn write_diagnostics(
        &self,
        verbose: u8,
        quiet: bool,
        mut stdout: impl Write,
        mut stderr: impl Write,
    ) -> Result<()> {
        for file in self.files.iter().filter(|file| !file.is_ok()) {
            if verbose > 0 {
                writeln!(stdout, "{}", file.path.display()).map_err(RsomicsError::Io)?;
            }
            if !quiet {
                for problem in &file.problems {
                    writeln!(stderr, "{}", diagnostic(&file.path, *problem))
                        .map_err(RsomicsError::Io)?;
                }
            }
        }
        stdout.flush().map_err(RsomicsError::Io)?;
        stderr.flush().map_err(RsomicsError::Io)
    }
}

#[must_use]
pub fn check(path: &Path, options: Options) -> FileReport {
    hts_quickcheck::check(path, options.allow_no_targets)
}

#[must_use]
pub fn check_all(paths: &[PathBuf], options: Options) -> Report {
    Report {
        files: paths.iter().map(|path| check(path, options)).collect(),
    }
}

fn diagnostic(path: &Path, problem: Problem) -> String {
    let path = path.display();
    match problem {
        Problem::Open => format!("{path} could not be opened for reading."),
        Problem::NotSequence => format!("{path} was not identified as sequence data."),
        Problem::Header => format!("{path} caused an error whilst reading its header."),
        Problem::NoTargets => format!("{path} had no targets in header."),
        Problem::EofCheck => format!("{path} caused an error whilst checking for EOF block."),
        Problem::MissingEof => {
            format!("{path} was missing EOF block when one should be present.")
        }
        Problem::Close => format!("{path} did not close cleanly."),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagnostics_list_only_failed_inputs() {
        let report = Report {
            files: vec![
                FileReport {
                    path: PathBuf::from("ok.bam"),
                    problems: Vec::new(),
                },
                FileReport {
                    path: PathBuf::from("bad.bam"),
                    problems: vec![Problem::Header, Problem::MissingEof],
                },
            ],
        };
        let mut stdout = Vec::new();
        let mut stderr = Vec::new();
        report
            .write_diagnostics(1, false, &mut stdout, &mut stderr)
            .unwrap();

        assert_eq!(stdout, b"bad.bam\n");
        let stderr = String::from_utf8(stderr).unwrap();
        assert!(stderr.contains("bad.bam caused an error whilst reading its header."));
        assert!(stderr.contains("bad.bam was missing EOF block"));
        assert!(!stderr.contains("ok.bam"));
    }
}
