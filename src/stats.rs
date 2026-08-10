use std::collections::{BTreeMap, HashSet};
use std::path::Path;
use std::str::FromStr;

use noodles::sam::header::record::value::map::read_group::tag as read_group_tag;
use rsomics_common::{Result, RsomicsError};
use rust_htslib::bam::Read as _;
use serde::ser::{SerializeMap, SerializeStruct};
use serde::{Serialize, Serializer};

use crate::input;

mod barcode;
mod checksum;
mod coverage;
mod record;
mod record_data;
mod ref_stats;
mod reference;
mod regions;
mod render;

const MAX_COVERAGE_INTERVALS: usize = 1_000_000;
const MAX_SPLIT_REPORTS: usize = 4096;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub struct CoverageBins {
    pub minimum: usize,
    pub maximum: usize,
    pub step: usize,
}

impl Default for CoverageBins {
    fn default() -> Self {
        Self {
            minimum: 1,
            maximum: 1000,
            step: 1,
        }
    }
}

impl CoverageBins {
    fn normalized(self) -> std::result::Result<Self, String> {
        if self.minimum > self.maximum || self.step == 0 {
            return Err(format!(
                "invalid coverage bins: {},{},{}",
                self.minimum, self.maximum, self.step
            ));
        }
        let difference = self.maximum - self.minimum;
        let step = if difference == usize::MAX {
            self.step
        } else {
            self.step.min(difference + 1)
        };
        let quotient = difference / step;
        if quotient >= MAX_COVERAGE_INTERVALS {
            return Err(format!(
                "coverage bins exceed {MAX_COVERAGE_INTERVALS} intervals"
            ));
        }
        let maximum = self
            .minimum
            .checked_add(quotient * step)
            .and_then(|value| value.checked_add(step - 1))
            .ok_or_else(|| "coverage bin boundary overflows".to_owned())?;
        Ok(Self {
            minimum: self.minimum,
            maximum,
            step,
        })
    }

    pub(crate) fn count(self) -> usize {
        3 + (self.maximum - self.minimum) / self.step
    }

    pub(crate) fn index(self, depth: usize) -> usize {
        if depth < self.minimum {
            0
        } else if depth > self.maximum {
            self.count() - 1
        } else {
            1 + (depth - self.minimum) / self.step
        }
    }
}

impl FromStr for CoverageBins {
    type Err = String;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        let values = value
            .split(',')
            .map(str::parse::<usize>)
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(|_| format!("invalid coverage bins: {value}"))?;
        if values.len() != 3 || values[0] > values[1] || values[2] == 0 {
            return Err(format!("invalid coverage bins: {value}"));
        }
        Self {
            minimum: values[0],
            maximum: values[1],
            step: values[2],
        }
        .normalized()
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub reference: Option<&'a Path>,
    pub additional_threads: usize,
    pub required_flags: u16,
    pub filtered_flags: u16,
    pub read_length: Option<usize>,
    pub coverage: CoverageBins,
    pub maximum_insert_size: usize,
    pub insert_bulk: f64,
    pub trim_quality: u8,
    pub gc_depth: f64,
    pub sparse: bool,
    pub coverage_threshold: usize,
    pub targets: Option<&'a Path>,
    pub regions: &'a [String],
    pub index: Option<&'a Path>,
    pub id: Option<&'a str>,
    pub split_tag: Option<[u8; 2]>,
    pub remove_overlaps: bool,
    pub reference_stats: bool,
    pub reference_stats_chunk_mib: usize,
}

impl Default for Options<'_> {
    fn default() -> Self {
        Self {
            reference: None,
            additional_threads: 0,
            required_flags: 0,
            filtered_flags: 0,
            read_length: None,
            coverage: CoverageBins::default(),
            maximum_insert_size: 8000,
            insert_bulk: 0.99,
            trim_quality: 0,
            gc_depth: 20_000.0,
            sparse: false,
            coverage_threshold: 0,
            targets: None,
            regions: &[],
            index: None,
            id: None,
            split_tag: None,
            remove_overlaps: false,
            reference_stats: false,
            reference_stats_chunk_mib: 1,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Report {
    stats: record::Accumulator,
    splits: BTreeMap<Vec<u8>, record::Accumulator>,
    options: OwnedOptions,
}

impl Serialize for Report {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let options = self.options.borrowed();
        let mut report = serializer.serialize_struct("Report", 2)?;
        report.serialize_field("stats", &StatsView::new(&self.stats, options))?;
        report.serialize_field(
            "splits",
            &SplitReports {
                splits: &self.splits,
                options,
            },
        )?;
        report.end()
    }
}

#[derive(Serialize)]
struct StatsView<'a> {
    #[serde(flatten)]
    stats: &'a record::Accumulator,
    coverage: CoverageView,
}

impl<'a> StatsView<'a> {
    fn new(stats: &'a record::Accumulator, options: Options<'_>) -> Self {
        let histogram = stats
            .sorted
            .then(|| stats.coverage_histogram(options.coverage));
        let bases_above = stats
            .target_bases
            .map(|_| stats.coverage.bases_above(options.coverage_threshold));
        Self {
            stats,
            coverage: CoverageView {
                minimum: options.coverage.minimum,
                maximum: options.coverage.maximum,
                step: options.coverage.step,
                histogram,
                threshold: options.coverage_threshold,
                bases_above,
            },
        }
    }
}

#[derive(Serialize)]
struct CoverageView {
    minimum: usize,
    maximum: usize,
    step: usize,
    histogram: Option<Vec<u64>>,
    threshold: usize,
    bases_above: Option<u64>,
}

struct SplitReports<'a> {
    splits: &'a BTreeMap<Vec<u8>, record::Accumulator>,
    options: Options<'a>,
}

impl Serialize for SplitReports<'_> {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut map = serializer.serialize_map(Some(self.splits.len()))?;
        for (value, stats) in self.splits {
            let value = std::str::from_utf8(value).map_err(serde::ser::Error::custom)?;
            map.serialize_entry(value, &StatsView::new(stats, self.options))?;
        }
        map.end()
    }
}

#[derive(Clone, Debug)]
struct OwnedOptions {
    reference: Option<std::path::PathBuf>,
    additional_threads: usize,
    required_flags: u16,
    filtered_flags: u16,
    read_length: Option<usize>,
    coverage: CoverageBins,
    maximum_insert_size: usize,
    insert_bulk: f64,
    trim_quality: u8,
    gc_depth: f64,
    sparse: bool,
    coverage_threshold: usize,
    targets: Option<std::path::PathBuf>,
    regions: Vec<String>,
    index: Option<std::path::PathBuf>,
    id: Option<String>,
    split_tag: Option<[u8; 2]>,
    remove_overlaps: bool,
    reference_stats: bool,
    reference_stats_chunk_mib: usize,
}

impl OwnedOptions {
    fn borrowed(&self) -> Options<'_> {
        Options {
            reference: self.reference.as_deref(),
            additional_threads: self.additional_threads,
            required_flags: self.required_flags,
            filtered_flags: self.filtered_flags,
            read_length: self.read_length,
            coverage: self.coverage,
            maximum_insert_size: self.maximum_insert_size,
            insert_bulk: self.insert_bulk,
            trim_quality: self.trim_quality,
            gc_depth: self.gc_depth,
            sparse: self.sparse,
            coverage_threshold: self.coverage_threshold,
            targets: self.targets.as_deref(),
            regions: &self.regions,
            index: self.index.as_deref(),
            id: self.id.as_deref(),
            split_tag: self.split_tag,
            remove_overlaps: self.remove_overlaps,
            reference_stats: self.reference_stats,
            reference_stats_chunk_mib: self.reference_stats_chunk_mib,
        }
    }
}

pub fn collect(path: &Path, mut options: Options<'_>) -> Result<Report> {
    options.coverage = options
        .coverage
        .normalized()
        .map_err(RsomicsError::ConfigError)?;
    if !options.insert_bulk.is_finite() || !(0.0..=1.0).contains(&options.insert_bulk) {
        return Err(RsomicsError::ConfigError(
            "insert bulk must be between zero and one".to_owned(),
        ));
    }
    if !options.gc_depth.is_finite() || options.gc_depth <= 0.0 {
        return Err(RsomicsError::ConfigError(
            "GC-depth must be greater than zero".to_owned(),
        ));
    }
    let reference_chunk = options
        .reference_stats_chunk_mib
        .checked_mul(1024 * 1024)
        .ok_or_else(|| RsomicsError::ConfigError("reference chunk size overflows".to_owned()))?;
    let mut reference = options
        .reference
        .map(|path| reference::Reference::open(path, reference_chunk))
        .transpose()?;
    let has_reference = reference.is_some();
    let format = (path != Path::new("-"))
        .then(|| input::detect_format(path))
        .transpose()?;
    let alignment_reference = if format.is_none_or(|format| format == input::Format::Cram) {
        options.reference
    } else {
        None
    };
    if !options.regions.is_empty() || options.index.is_some() {
        let mut indexed = if let Some(index) = options.index {
            input::open_indexed_with_index(path, index, alignment_reference)?
        } else {
            input::open_indexed(path, alignment_reference)?
        };
        indexed.read_header(path)?;
    }
    let reader_threads = if format == Some(input::Format::Cram) {
        0
    } else {
        options.additional_threads
    };
    let mut reader = input::open(path, alignment_reference, reader_threads)?;
    let header = reader.read_header(path)?;
    let references = header
        .reference_sequences()
        .iter()
        .map(|(name, map)| {
            let bytes: &[u8] = name.as_ref();
            (
                bytes.to_vec(),
                u64::try_from(usize::from(map.length())).unwrap(),
            )
        })
        .collect::<Vec<_>>();
    let read_groups = options.id.map(|id| {
        header
            .read_groups()
            .iter()
            .filter_map(|(name, read_group)| {
                let name_bytes: &[u8] = name.as_ref();
                let id_matches = name_bytes == id.as_bytes();
                let sample_matches = read_group
                    .other_fields()
                    .get(&read_group_tag::SAMPLE)
                    .is_some_and(|value| {
                        let value: &[u8] = value.as_ref();
                        value == id.as_bytes()
                    });
                (id_matches || sample_matches).then(|| name.to_vec())
            })
            .collect::<HashSet<_>>()
    });
    let target_regions = options
        .targets
        .map(|path| regions::Regions::from_targets(path, &references))
        .transpose()?;
    let cli_regions = (!options.regions.is_empty())
        .then(|| regions::Regions::from_cli(options.regions, &references))
        .transpose()?;
    let selected = match (target_regions, cli_regions) {
        (Some(targets), Some(cli)) => Some(targets.intersect(cli)?),
        (Some(targets), None) => Some(targets),
        (None, Some(cli)) => Some(cli),
        (None, None) => None,
    };
    if options.coverage_threshold > 0 && selected.is_none() {
        return Err(RsomicsError::ConfigError(
            "coverage threshold requires target regions or indexed regions".to_owned(),
        ));
    }
    let mut stats = record::Accumulator::new(
        has_reference,
        selected.as_ref().map(regions::Regions::bases),
    );
    let mut splits = BTreeMap::<Vec<u8>, record::Accumulator>::new();
    let mut process = |record: &record_data::RecordData| {
        let split_value = record.split_value.as_ref();
        if let Some(value) = split_value {
            validate_split_value(value)?;
            if splits.len() == MAX_SPLIT_REPORTS && !splits.contains_key(value) {
                return Err(RsomicsError::InvalidInput(format!(
                    "split tag has more than {MAX_SPLIT_REPORTS} distinct values"
                )));
            }
            splits.entry(value.clone()).or_insert_with(|| {
                record::Accumulator::new(
                    has_reference,
                    selected.as_ref().map(regions::Regions::bases),
                )
            });
        }
        let selected_read_group = read_groups.as_ref().is_none_or(|read_groups| {
            record
                .read_group
                .as_ref()
                .is_some_and(|value| read_groups.contains(value))
        });
        if !selected_read_group {
            return Ok(true);
        }
        if let Some(selected) = &selected {
            let end = record.reference_end()?;
            if !selected.overlaps(record.reference, record.position, end) {
                return Ok(true);
            }
        }
        let sequence = if has_reference && record.reference >= 0 {
            let name = usize::try_from(record.reference)
                .ok()
                .and_then(|id| references.get(id))
                .map(|(name, _)| name)
                .ok_or_else(|| {
                    RsomicsError::InvalidInput(format!(
                        "mapped record has unknown reference ID {}",
                        record.reference
                    ))
                })?;
            let reference = reference.as_mut().unwrap();
            let length = reference.length(name).ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "reference FASTA is missing mapped sequence {}",
                    String::from_utf8_lossy(name)
                ))
            })?;
            let start = usize::try_from(record.position).map_err(|_| {
                RsomicsError::InvalidInput("mapped read has a negative position".to_owned())
            })?;
            let alignment_end = usize::try_from(record.reference_end()?)
                .map_err(|_| RsomicsError::InvalidInput("alignment end is negative".to_owned()))?;
            if alignment_end > length {
                return Err(RsomicsError::InvalidInput(
                    "alignment extends beyond its reference".to_owned(),
                ));
            }
            let gc_end = start
                .saturating_add(options.gc_depth.ceil() as usize)
                .min(length);
            Some(reference.get(name, start, alignment_end.max(gc_end))?)
        } else {
            None
        };
        stats.collect(record, sequence, selected.as_ref(), options)?;
        if let Some(value) = split_value {
            splits
                .get_mut(value)
                .expect("split accumulator was inserted")
                .collect(record, sequence, selected.as_ref(), options)?;
        }
        Ok(true)
    };
    if !options.regions.is_empty() && format == Some(input::Format::Bam) {
        visit_indexed_hts_records(path, options, selected.as_ref().unwrap(), &mut process)?;
    } else if format == Some(input::Format::Cram) {
        visit_hts_records(path, options, &mut process)?;
    } else if format == Some(input::Format::Bam) {
        let mut record = record_data::RecordData::default();
        reader.visit_raw_bam_records(path, |source| {
            record.decode_raw(&source, options.split_tag)?;
            process(&record)
        })?;
    } else {
        reader.visit_records(&header, path, |source| {
            let record = record_data::RecordData::decode(&header, source, options.split_tag)?;
            process(&record)
        })?;
    }
    if options.reference_stats {
        stats.reference_stats = Some(ref_stats::ReferenceStats::collect(
            &references,
            selected.as_ref(),
            reference.as_mut(),
        )?);
    }
    Ok(Report {
        stats,
        splits,
        options: OwnedOptions {
            reference: options.reference.map(Path::to_path_buf),
            additional_threads: options.additional_threads,
            required_flags: options.required_flags,
            filtered_flags: options.filtered_flags,
            read_length: options.read_length,
            coverage: options.coverage,
            maximum_insert_size: options.maximum_insert_size,
            insert_bulk: options.insert_bulk,
            trim_quality: options.trim_quality,
            gc_depth: options.gc_depth,
            sparse: options.sparse,
            coverage_threshold: options.coverage_threshold,
            targets: options.targets.map(Path::to_path_buf),
            regions: options.regions.to_vec(),
            index: options.index.map(Path::to_path_buf),
            id: options.id.map(str::to_owned),
            split_tag: options.split_tag,
            remove_overlaps: options.remove_overlaps,
            reference_stats: options.reference_stats,
            reference_stats_chunk_mib: options.reference_stats_chunk_mib,
        },
    })
}

fn visit_indexed_hts_records(
    path: &Path,
    options: Options<'_>,
    regions: &regions::Regions,
    mut process: impl FnMut(&record_data::RecordData) -> Result<bool>,
) -> Result<()> {
    let mut reader = if let Some(index) = options.index {
        rust_htslib::bam::IndexedReader::from_path_and_index(path, index)
    } else {
        rust_htslib::bam::IndexedReader::from_path(path)
    }
    .map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "opening indexed alignment {}: {error}",
            path.display()
        ))
    })?;
    if let Some(reference) = options.reference {
        reader.set_reference(reference).map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "setting alignment reference {}: {error}",
                reference.display()
            ))
        })?;
    }
    if options.additional_threads > 0 {
        reader
            .set_threads(options.additional_threads)
            .map_err(|error| {
                RsomicsError::ConfigError(format!("setting alignment threads: {error}"))
            })?;
    }

    let intervals = regions.iter().collect::<Vec<_>>();
    let mut previous_boundary = HashSet::new();
    let mut source = rust_htslib::bam::Record::new();
    let mut record = record_data::RecordData::default();
    for (index, (reference, interval)) in intervals.iter().enumerate() {
        reader
            .fetch((*reference, interval.start, interval.end))
            .map_err(|error| {
                RsomicsError::InvalidInput(format!(
                    "querying {}:{}-{}: {error}",
                    path.display(),
                    interval.start + 1,
                    interval.end
                ))
            })?;
        let next = intervals.get(index + 1);
        let mut current_boundary = HashSet::new();
        while let Some(result) = reader.read(&mut source) {
            result.map_err(|error| {
                RsomicsError::InvalidInput(format!(
                    "reading indexed alignment record from {}: {error}",
                    path.display()
                ))
            })?;
            let offset = reader.tell();
            record.decode_hts(&source, options.split_tag)?;
            let duplicate = previous_boundary.remove(&offset);
            if !duplicate && !process(&record)? {
                return Ok(());
            }
            if let Some((next_reference, next_interval)) = next
                && next_reference == reference
                && record.reference_end()? > next_interval.start
            {
                current_boundary.insert(offset);
            }
        }
        previous_boundary = current_boundary;
    }
    Ok(())
}

fn visit_hts_records(
    path: &Path,
    options: Options<'_>,
    mut process: impl FnMut(&record_data::RecordData) -> Result<bool>,
) -> Result<()> {
    let mut reader = rust_htslib::bam::Reader::from_path(path).map_err(|error| {
        RsomicsError::InvalidInput(format!("opening alignment {}: {error}", path.display()))
    })?;
    if let Some(reference) = options.reference {
        reader.set_reference(reference).map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "setting CRAM reference {}: {error}",
                reference.display()
            ))
        })?;
    }
    if options.additional_threads > 0 {
        reader
            .set_threads(options.additional_threads)
            .map_err(|error| {
                RsomicsError::ConfigError(format!("setting alignment threads: {error}"))
            })?;
    }
    let mut source = rust_htslib::bam::Record::new();
    let mut record = record_data::RecordData::default();
    while let Some(result) = reader.read(&mut source) {
        result.map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "reading alignment record from {}: {error}",
                path.display()
            ))
        })?;
        record.decode_hts(&source, options.split_tag)?;
        if !process(&record)? {
            break;
        }
    }
    Ok(())
}

impl Report {
    pub fn write(&self, output: impl std::io::Write) -> Result<()> {
        render::write(output, &self.stats, self.options.borrowed(), None)
    }

    pub(crate) fn split_values(&self) -> impl Iterator<Item = &[u8]> {
        self.splits.keys().map(Vec::as_slice)
    }

    pub(crate) fn write_split(&self, value: &[u8], output: impl std::io::Write) -> Result<()> {
        let stats = self
            .splits
            .get(value)
            .ok_or_else(|| RsomicsError::ConfigError("unknown split report value".to_owned()))?;
        render::write(
            output,
            stats,
            self.options.borrowed(),
            self.options.split_tag.map(|tag| (tag, value)),
        )
    }
}

fn validate_split_value(value: &[u8]) -> Result<()> {
    if value.is_empty()
        || value == b"."
        || value == b".."
        || value
            .iter()
            .any(|byte| !byte.is_ascii_alphanumeric() && !matches!(byte, b'.' | b'_' | b'-'))
    {
        return Err(RsomicsError::InvalidInput(format!(
            "split tag value is not a safe file-name component: {:?}",
            String::from_utf8_lossy(value)
        )));
    }
    Ok(())
}
