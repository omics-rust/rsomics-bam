use std::io::Write;
use std::path::{Path, PathBuf};

use rsomics_common::Result;
use serde::Serialize;

use crate::{Program, sort};

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub memory_limit: u64,
    pub additional_threads: Option<usize>,
    pub temporary_prefix: Option<&'a Path>,
    pub reference: Option<&'a Path>,
    pub destination: Option<&'a Path>,
    pub program: Option<Program<'a>>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub input: PathBuf,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output: Option<PathBuf>,
    pub records: u64,
    pub memory_limit: u64,
    pub additional_threads: usize,
    pub temporary_runs: u64,
    pub merge_passes: u32,
}

pub fn write<W>(input_path: &Path, options: Options<'_>, output: W) -> Result<Summary>
where
    W: Write + Send + 'static,
{
    let stats = sort::write_collated(
        input_path,
        sort::CollateOptions {
            memory_limit: options.memory_limit,
            additional_threads: options.additional_threads,
            temporary_prefix: options.temporary_prefix,
            reference: options.reference,
            program: options.program,
        },
        output,
    )?;
    Ok(Summary {
        input: input_path.to_path_buf(),
        output: options.destination.map(Path::to_path_buf),
        records: stats.records,
        memory_limit: stats.memory_limit,
        additional_threads: stats.additional_threads,
        temporary_runs: stats.temporary_runs,
        merge_passes: stats.merge_passes,
    })
}
