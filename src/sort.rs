use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread;

use noodles::{bam, bgzf, sam};
use rayon::prelude::*;
use rsomics_bamio::raw::{self, RawRecord, RawRecordEncoder, RecordRef};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

pub use crate::alignment_order::Order;
use crate::alignment_order::{
    OrderedRecord, compare_ordered_records, library_lookup, ordered_record, set_sort_order,
};
use crate::hts_quickcheck::{require_bgzf_eof, require_cram_eof};
use crate::{Program, input, md, output};

const MIN_MEMORY: u64 = 1 << 20;
const MAX_FAN_IN: usize = 32;

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub order: Order,
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
    pub order: Order,
    pub records: u64,
    pub memory_limit: u64,
    pub additional_threads: usize,
    pub temporary_runs: u64,
    pub merge_passes: u32,
}

struct RunFile {
    file: tempfile::NamedTempFile,
}

struct TempLayout {
    directory: PathBuf,
    prefix: OsString,
}

type BamReader = bam::io::Reader<bgzf::io::Reader<BufReader<File>>>;

struct RunReader {
    reader: BamReader,
}

struct HeapEntry {
    entry: OrderedRecord,
    order: Order,
    run: usize,
    sequence: u64,
}

impl PartialEq for HeapEntry {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for HeapEntry {}

impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        compare_ordered_records(self.order, &self.entry, &other.entry)
            .then_with(|| self.run.cmp(&other.run))
            .then_with(|| self.sequence.cmp(&other.sequence))
            .reverse()
    }
}

impl TempLayout {
    fn new(prefix: Option<&Path>) -> Result<Self> {
        match prefix {
            Some(prefix) => {
                let name = prefix.file_name().ok_or_else(|| {
                    RsomicsError::ConfigError(format!(
                        "temporary prefix has no file name: {}",
                        prefix.display()
                    ))
                })?;
                let directory = prefix
                    .parent()
                    .filter(|parent| !parent.as_os_str().is_empty())
                    .unwrap_or_else(|| Path::new("."));
                if !directory.is_dir() {
                    return Err(RsomicsError::ConfigError(format!(
                        "temporary directory does not exist: {}",
                        directory.display()
                    )));
                }
                let mut name = name.to_os_string();
                name.push(".");
                Ok(Self {
                    directory: directory.to_path_buf(),
                    prefix: name,
                })
            }
            None => Ok(Self {
                directory: std::env::temp_dir(),
                prefix: OsString::from("rsomics-bam-sort."),
            }),
        }
    }

    fn create(&self) -> Result<RunFile> {
        tempfile::Builder::new()
            .prefix(&self.prefix)
            .suffix(".bam")
            .tempfile_in(&self.directory)
            .map(|file| RunFile { file })
            .map_err(RsomicsError::Io)
    }
}

impl RunReader {
    fn open(run: &RunFile, expected_header: &sam::Header) -> Result<Self> {
        let file = run.file.reopen().map_err(RsomicsError::Io)?;
        let mut reader = bam::io::Reader::from(bgzf::io::Reader::new(BufReader::new(file)));
        let header = reader.read_header().map_err(RsomicsError::Io)?;
        if &header != expected_header {
            return Err(RsomicsError::InvalidInput(
                "temporary sort run header changed during merge".to_owned(),
            ));
        }
        Ok(Self { reader })
    }

    fn next(&mut self) -> Result<Option<RawRecord>> {
        let mut record = RawRecord::default();
        match raw::read_record(self.reader.get_mut(), &mut record)? {
            0 => Ok(None),
            _ => Ok(Some(record)),
        }
    }
}

pub fn write<W>(input_path: &Path, options: Options<'_>, output: W) -> Result<Summary>
where
    W: Write + Send + 'static,
{
    if options.memory_limit < MIN_MEMORY {
        return Err(RsomicsError::ConfigError(
            "sort memory must be at least 1 MiB".to_owned(),
        ));
    }
    let additional_threads = options
        .additional_threads
        .unwrap_or_else(default_additional_threads);
    if additional_threads > 256 {
        return Err(RsomicsError::ConfigError(
            "sort additional thread count cannot exceed 256".to_owned(),
        ));
    }
    let total_threads = additional_threads
        .checked_add(1)
        .ok_or_else(|| RsomicsError::ConfigError("sort thread count overflows".to_owned()))?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(total_threads)
        .build()
        .map_err(|error| RsomicsError::ConfigError(format!("creating sort workers: {error}")))?;

    let named_format = if input_path == Path::new("-") {
        None
    } else {
        let format = input::detect_format(input_path)?;
        match format {
            input::Format::Bam | input::Format::Sam if input::is_bgzf(input_path)? => {
                require_bgzf_eof(input_path)?;
            }
            input::Format::Cram => require_cram_eof(input_path)?,
            input::Format::Bam | input::Format::Sam => {}
        }
        Some(format)
    };
    let input_threads = match named_format {
        Some(input::Format::Bam) => additional_threads,
        _ => 0,
    };
    let mut reader = input::open(input_path, options.reference, input_threads)?;
    let input_format = reader.format();
    let mut header = reader.read_header(input_path)?;
    set_sort_order(&mut header, options.order);
    if let Some(program) = options.program {
        program.add_to(&mut header)?;
    }
    let libraries = library_lookup(&header);
    let layout = TempLayout::new(options.temporary_prefix)?;
    let mut entries = Vec::new();
    let mut runs = Vec::new();
    let mut memory = 0u64;
    let mut records = 0u64;

    let mut ingest = |record| {
        let ordinal = records;
        records = records.checked_add(1).ok_or_else(count_overflow)?;
        let entry = ordered_record(record, options.order, &header, &libraries, ordinal)?;
        let entry_memory = entry.memory()?;
        if !entries.is_empty()
            && memory
                .checked_add(entry_memory)
                .is_none_or(|value| value > options.memory_limit)
        {
            runs.push(write_run(
                &mut entries,
                &header,
                options.order,
                &pool,
                &layout,
                additional_threads,
            )?);
            memory = 0;
        }
        memory = memory.saturating_add(entry_memory);
        entries.push(entry);
        Ok(true)
    };
    if input_format == input::Format::Cram {
        let mut reference = options
            .reference
            .map(md::ReferenceCache::open)
            .transpose()?;
        let mut encoder = RawRecordEncoder::new();
        reader.visit_records(&header, input_path, |record| {
            let record = md::complete(&header, record, reference.as_mut())?;
            ingest(encoder.encode(&header, &record)?)
        })?;
    } else {
        reader.visit_owned_raw_records(&header, input_path, ingest)?;
    }

    let initial_runs;
    let merge_passes;
    if runs.is_empty() {
        sort_entries(&mut entries, options.order, &pool);
        write_entries(output, &header, &entries, additional_threads)?;
        initial_runs = 0;
        merge_passes = 0;
    } else {
        if !entries.is_empty() {
            runs.push(write_run(
                &mut entries,
                &header,
                options.order,
                &pool,
                &layout,
                additional_threads,
            )?);
        }
        initial_runs = u64::try_from(runs.len()).map_err(|_| {
            RsomicsError::InvalidInput("temporary run count exceeds u64".to_owned())
        })?;
        let (runs, consolidation_passes) = consolidate_runs(
            runs,
            &header,
            options.order,
            &libraries,
            &layout,
            additional_threads,
        )?;
        merge_runs_to_writer(
            &runs,
            &header,
            options.order,
            &libraries,
            output,
            additional_threads,
        )?;
        merge_passes = consolidation_passes.checked_add(1).ok_or_else(|| {
            RsomicsError::InvalidInput("sort merge pass count exceeds u32".to_owned())
        })?;
    }

    Ok(Summary {
        input: input_path.to_path_buf(),
        output: options.destination.map(Path::to_path_buf),
        order: options.order,
        records,
        memory_limit: options.memory_limit,
        additional_threads,
        temporary_runs: initial_runs,
        merge_passes,
    })
}

fn write_entries<W>(
    output: W,
    header: &sam::Header,
    entries: &[OrderedRecord],
    additional_threads: usize,
) -> Result<()>
where
    W: Write + Send + 'static,
{
    let mut writer = output::Writer::new(
        output::Format::Bam,
        output::Compression::Default,
        additional_threads,
        output,
    );
    writer.write_header(header)?;
    for entry in entries {
        write_record(&mut writer, &entry.record)?;
    }
    writer.finish(header)
}

fn write_run(
    entries: &mut Vec<OrderedRecord>,
    header: &sam::Header,
    order: Order,
    pool: &rayon::ThreadPool,
    layout: &TempLayout,
    additional_threads: usize,
) -> Result<RunFile> {
    sort_entries(entries, order, pool);
    let run = layout.create()?;
    let file = run.file.reopen().map_err(RsomicsError::Io)?;
    let mut writer = output::Writer::new(
        output::Format::Bam,
        output::Compression::Fast,
        additional_threads,
        file,
    );
    writer.write_header(header)?;
    for entry in entries.iter() {
        write_record(&mut writer, &entry.record)?;
    }
    writer.finish(header)?;
    entries.clear();
    Ok(run)
}

fn write_record<W>(writer: &mut output::Writer<W>, record: &RawRecord) -> Result<()>
where
    W: Write + Send + 'static,
{
    let record = RecordRef::from_bytes(record.as_bytes())?;
    writer.write_raw_record(&record)
}

fn sort_entries(entries: &mut [OrderedRecord], order: Order, pool: &rayon::ThreadPool) {
    pool.install(|| {
        entries.par_sort_unstable_by(|a, b| {
            compare_ordered_records(order, a, b).then_with(|| a.ordinal.cmp(&b.ordinal))
        });
    });
}

fn consolidate_runs(
    mut runs: Vec<RunFile>,
    header: &sam::Header,
    order: Order,
    libraries: &HashMap<Vec<u8>, Arc<[u8]>>,
    layout: &TempLayout,
    additional_threads: usize,
) -> Result<(Vec<RunFile>, u32)> {
    let mut passes = 0u32;
    while runs.len() > MAX_FAN_IN {
        let mut source = runs.into_iter();
        let mut next = Vec::new();
        loop {
            let group = source.by_ref().take(MAX_FAN_IN).collect::<Vec<_>>();
            if group.is_empty() {
                break;
            }
            if group.len() == 1 {
                next.extend(group);
                continue;
            }
            let merged = layout.create()?;
            let file = merged.file.reopen().map_err(RsomicsError::Io)?;
            merge_runs_to_writer(&group, header, order, libraries, file, additional_threads)?;
            next.push(merged);
        }
        runs = next;
        passes = passes.checked_add(1).ok_or_else(|| {
            RsomicsError::InvalidInput("sort merge pass count exceeds u32".to_owned())
        })?;
    }
    Ok((runs, passes))
}

fn merge_runs_to_writer<W>(
    runs: &[RunFile],
    header: &sam::Header,
    order: Order,
    libraries: &HashMap<Vec<u8>, Arc<[u8]>>,
    output: W,
    additional_threads: usize,
) -> Result<()>
where
    W: Write + Send + 'static,
{
    let mut readers = runs
        .iter()
        .map(|run| RunReader::open(run, header))
        .collect::<Result<Vec<_>>>()?;
    let mut heap = BinaryHeap::new();
    for (run, reader) in readers.iter_mut().enumerate() {
        if let Some(record) = reader.next()? {
            heap.push(HeapEntry {
                entry: ordered_record(record, order, header, libraries, 0)?,
                order,
                run,
                sequence: 0,
            });
        }
    }

    let mut writer = output::Writer::new(
        output::Format::Bam,
        output::Compression::Default,
        additional_threads,
        output,
    );
    writer.write_header(header)?;
    while let Some(item) = heap.pop() {
        write_record(&mut writer, &item.entry.record)?;
        let sequence = item.sequence.checked_add(1).ok_or_else(count_overflow)?;
        if let Some(record) = readers[item.run].next()? {
            heap.push(HeapEntry {
                entry: ordered_record(record, order, header, libraries, sequence)?,
                order,
                run: item.run,
                sequence,
            });
        }
    }
    writer.finish(header)
}

pub(crate) fn default_additional_threads() -> usize {
    thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1)
        .saturating_sub(1)
        .min(4)
}

fn count_overflow() -> RsomicsError {
    RsomicsError::InvalidInput("alignment record count exceeds u64".to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn automatic_threads_are_bounded() {
        assert!(default_additional_threads() <= 4);
    }
}
