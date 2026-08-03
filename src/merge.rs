use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, VecDeque};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, mpsc};
use std::thread;

use noodles::sam;
use noodles::sam::header::record::value::map::header::tag as header_tag;
use rsomics_bamio::raw::{RawRecord, RawRecordEncoder, RecordRef};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

pub use crate::alignment_order::Order;
use crate::alignment_order::{
    OrderedRecord, compare_ordered_records, library_lookup, ordered_record, set_sort_order,
};
use crate::hts_quickcheck::{require_bgzf_eof, require_cram_eof};
use crate::{Program, header_merge, input, md, output};

const MAX_INPUTS: usize = 32;
const INPUT_BATCH_BYTES: usize = 1 << 20;

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub order: Order,
    pub additional_threads: Option<usize>,
    pub reference: Option<&'a Path>,
    pub destination: Option<&'a Path>,
    pub combine_read_groups: bool,
    pub combine_programs: bool,
    pub program: Option<Program<'a>>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub inputs: Vec<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output: Option<PathBuf>,
    pub order: Order,
    pub records: u64,
    pub additional_threads: usize,
}

struct InputState {
    path: PathBuf,
    translation: header_merge::Translation,
    receiver: mpsc::Receiver<Vec<RawRecord>>,
    pending: VecDeque<RawRecord>,
    worker: Option<thread::JoinHandle<Result<()>>>,
    previous: Option<OrderedRecord>,
    sequence: u64,
}

struct HeapEntry {
    record: OrderedRecord,
    order: Order,
    input: usize,
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
        compare_ordered_records(self.order.into(), &self.record, &other.record)
            .then_with(|| self.input.cmp(&other.input))
            .then_with(|| self.sequence.cmp(&other.sequence))
            .reverse()
    }
}

impl InputState {
    fn next(
        &mut self,
        order: Order,
        output_header: &sam::Header,
        libraries: &HashMap<Vec<u8>, Arc<[u8]>>,
    ) -> Result<Option<OrderedRecord>> {
        let mut record = loop {
            if let Some(record) = self.pending.pop_front() {
                break record;
            }
            match self.receiver.recv() {
                Ok(batch) => self.pending = batch.into(),
                Err(_) => {
                    self.join_worker()?;
                    return Ok(None);
                }
            }
        };
        header_merge::translate(&mut record, &self.translation)?;
        let current = ordered_record(
            record,
            order.into(),
            output_header,
            libraries,
            self.sequence,
        )?;
        if self.previous.as_ref().is_some_and(|previous| {
            compare_ordered_records(order.into(), previous, &current) == Ordering::Greater
        }) {
            return Err(RsomicsError::InvalidInput(format!(
                "input {} is not ordered as {}",
                self.path.display(),
                order_name(order)
            )));
        }
        self.previous = Some(current.clone());
        self.sequence = self
            .sequence
            .checked_add(1)
            .ok_or_else(record_count_overflow)?;
        Ok(Some(current))
    }

    fn join_worker(&mut self) -> Result<()> {
        let Some(worker) = self.worker.take() else {
            return Ok(());
        };
        worker.join().map_err(|_| {
            RsomicsError::InvalidInput(format!(
                "alignment reader for {} panicked",
                self.path.display()
            ))
        })?
    }
}

impl Drop for InputState {
    fn drop(&mut self) {
        let (_, disconnected) = mpsc::sync_channel(1);
        self.receiver = disconnected;
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}

pub fn write<W>(inputs: &[PathBuf], options: Options<'_>, output: W) -> Result<Summary>
where
    W: Write + Send + 'static,
{
    validate_inputs(inputs)?;
    validate_eof_markers(inputs)?;

    let mut headers = Vec::with_capacity(inputs.len());
    for path in inputs {
        let mut reader = input::open(path, options.reference, 0)?;
        let header = reader.read_header(path)?;
        validate_declared_order(&header, options.order, path)?;
        headers.push(header);
    }

    let (mut output_header, translations) = header_merge::reconcile(
        &headers,
        header_merge::Options {
            combine_read_groups: options.combine_read_groups,
            combine_programs: options.combine_programs,
        },
    )?;
    if matches!(options.order, Order::Coordinate | Order::TemplateCoordinate)
        && translations
            .iter()
            .any(|translation| !translation.preserves_reference_order())
    {
        return Err(RsomicsError::InvalidInput(
            "merged reference dictionary would destroy input ordering".to_owned(),
        ));
    }
    set_sort_order(&mut output_header, options.order);
    if let Some(program) = options.program {
        program.add_to(&mut output_header)?;
    }

    let libraries = library_lookup(&output_header);
    let mut states = headers
        .into_iter()
        .zip(translations)
        .zip(inputs)
        .map(|((header, translation), path)| {
            let (sender, receiver) = mpsc::sync_channel(1);
            let worker_path = path.clone();
            let reference = options.reference.map(Path::to_path_buf);
            let worker = thread::Builder::new()
                .name(format!("merge-{}", path.display()))
                .spawn(move || stream_records(worker_path, header, reference, sender))
                .map_err(|error| {
                    RsomicsError::ConfigError(format!(
                        "starting alignment reader for {}: {error}",
                        path.display()
                    ))
                })?;
            Ok(InputState {
                path: path.clone(),
                translation,
                receiver,
                pending: VecDeque::new(),
                worker: Some(worker),
                previous: None,
                sequence: 0,
            })
        })
        .collect::<Result<Vec<_>>>()?;

    let additional_threads = options
        .additional_threads
        .unwrap_or_else(crate::sort::default_additional_threads);
    let mut writer = output::Writer::new(
        output::Format::Bam,
        output::Compression::Default,
        additional_threads,
        output,
    );
    writer.write_header(&output_header)?;

    let mut heap = BinaryHeap::new();
    for (input, state) in states.iter_mut().enumerate() {
        if let Some(record) = state.next(options.order, &output_header, &libraries)? {
            heap.push(HeapEntry {
                record,
                order: options.order,
                input,
                sequence: 0,
            });
        }
    }

    let mut records = 0u64;
    while let Some(entry) = heap.pop() {
        let record = RecordRef::from_bytes(entry.record.record.as_bytes())?;
        writer.write_raw_record(&record)?;
        records = records.checked_add(1).ok_or_else(record_count_overflow)?;
        if let Some(record) = states[entry.input].next(options.order, &output_header, &libraries)? {
            heap.push(HeapEntry {
                record,
                order: options.order,
                input: entry.input,
                sequence: entry
                    .sequence
                    .checked_add(1)
                    .ok_or_else(record_count_overflow)?,
            });
        }
    }
    writer.finish(&output_header)?;

    Ok(Summary {
        inputs: inputs.to_vec(),
        output: options.destination.map(Path::to_path_buf),
        order: options.order,
        records,
        additional_threads,
    })
}

fn stream_records(
    path: PathBuf,
    expected_header: sam::Header,
    reference: Option<PathBuf>,
    sender: mpsc::SyncSender<Vec<RawRecord>>,
) -> Result<()> {
    let mut reader = input::open(&path, reference.as_deref(), 0)?;
    let header = reader.read_header(&path)?;
    if header != expected_header {
        return Err(RsomicsError::InvalidInput(format!(
            "alignment header changed while opening {}",
            path.display()
        )));
    }
    let mut batches = BatchSender::new(sender);
    if reader.format() == input::Format::Cram {
        let mut cache = reference
            .as_deref()
            .map(md::ReferenceCache::open)
            .transpose()?;
        let mut encoder = RawRecordEncoder::new();
        reader.visit_records(&header, &path, |record| {
            let record = md::complete(&header, record, cache.as_mut())?;
            Ok(batches.push(encoder.encode(&header, &record)?))
        })?;
    } else {
        reader.visit_owned_raw_records(&header, &path, |record| Ok(batches.push(record)))?;
    }
    batches.finish();
    Ok(())
}

struct BatchSender {
    sender: mpsc::SyncSender<Vec<RawRecord>>,
    records: Vec<RawRecord>,
    bytes: usize,
}

impl BatchSender {
    fn new(sender: mpsc::SyncSender<Vec<RawRecord>>) -> Self {
        Self {
            sender,
            records: Vec::new(),
            bytes: 0,
        }
    }

    fn push(&mut self, record: RawRecord) -> bool {
        self.bytes = self.bytes.saturating_add(record.as_bytes().len());
        self.records.push(record);
        self.bytes < INPUT_BATCH_BYTES || self.flush()
    }

    fn finish(mut self) {
        self.flush();
    }

    fn flush(&mut self) -> bool {
        if self.records.is_empty() {
            return true;
        }
        self.bytes = 0;
        self.sender.send(std::mem::take(&mut self.records)).is_ok()
    }
}

fn validate_inputs(inputs: &[PathBuf]) -> Result<()> {
    if inputs.is_empty() {
        return Err(RsomicsError::ConfigError(
            "merge requires at least one input".to_owned(),
        ));
    }
    if inputs.len() > MAX_INPUTS {
        return Err(RsomicsError::ConfigError(format!(
            "merge accepts at most {MAX_INPUTS} inputs"
        )));
    }
    if inputs.iter().any(|path| path == Path::new("-")) {
        return Err(RsomicsError::ConfigError(
            "merge inputs must be named files".to_owned(),
        ));
    }
    Ok(())
}

fn validate_eof_markers(inputs: &[PathBuf]) -> Result<()> {
    for path in inputs {
        let format = input::detect_format(path)?;
        match format {
            input::Format::Bam | input::Format::Sam if input::is_bgzf(path)? => {
                require_bgzf_eof(path)?;
            }
            input::Format::Cram => require_cram_eof(path)?,
            input::Format::Bam | input::Format::Sam => {}
        }
    }
    Ok(())
}

fn validate_declared_order(header: &sam::Header, order: Order, path: &Path) -> Result<()> {
    let fields = header.header().map(|header| header.other_fields());
    let sort_order = fields
        .and_then(|fields| fields.get(&header_tag::SORT_ORDER))
        .map(|value| value.as_slice());
    let group_order = fields
        .and_then(|fields| fields.get(&header_tag::GROUP_ORDER))
        .map(|value| value.as_slice());
    let subsort = fields
        .and_then(|fields| fields.get(&header_tag::SUBSORT_ORDER))
        .map(|value| value.as_slice());
    let valid = match order {
        Order::Coordinate => sort_order == Some(b"coordinate"),
        Order::QueryNameNatural => {
            sort_order == Some(b"queryname") && matches!(subsort, None | Some(b"queryname:natural"))
        }
        Order::QueryNameLexicographical => {
            sort_order == Some(b"queryname") && subsort == Some(b"queryname:lexicographical")
        }
        Order::TemplateCoordinate => {
            sort_order == Some(b"unsorted")
                && group_order == Some(b"query")
                && subsort == Some(b"unsorted:template-coordinate")
        }
    };
    if valid {
        Ok(())
    } else {
        Err(RsomicsError::InvalidInput(format!(
            "input {} does not declare {} order",
            path.display(),
            order_name(order)
        )))
    }
}

fn order_name(order: Order) -> &'static str {
    match order {
        Order::Coordinate => "coordinate",
        Order::QueryNameNatural => "natural query-name",
        Order::QueryNameLexicographical => "lexicographical query-name",
        Order::TemplateCoordinate => "template-coordinate",
    }
}

fn record_count_overflow() -> RsomicsError {
    RsomicsError::InvalidInput("alignment record count exceeds u64".to_owned())
}
