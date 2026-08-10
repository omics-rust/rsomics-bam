use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::{collections::HashMap, ops::Range};

use flate2::read::MultiGzDecoder;
use noodles::core::{Position, Region as QueryRegion};
use noodles::sam;
use rayon::prelude::*;
use rsomics_bamio::raw::RawRecordEncoder;
use rsomics_common::{Result, RsomicsError};
use rsomics_pileup::{FlagFilter, PileupEngine, PileupError, PileupOptions, RecordFilter};
use rust_htslib::bam::Read as _;
use serde::Serialize;

use crate::input;

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub reference: Option<&'a Path>,
    pub indices: Option<&'a [PathBuf]>,
    pub additional_threads: usize,
    pub minimum_mapping_quality: u8,
    pub excluded_flags: u16,
    pub skip_deletions_and_skips: bool,
    pub depth_threshold: Option<usize>,
    pub maximum_depth: usize,
    pub read_count: bool,
    pub header: bool,
}

impl Default for Options<'_> {
    fn default() -> Self {
        Self {
            reference: None,
            indices: None,
            additional_threads: 0,
            minimum_mapping_quality: 0,
            excluded_flags: 0x704,
            skip_deletions_and_skips: false,
            depth_threshold: None,
            maximum_depth: i32::MAX as usize,
            read_count: false,
            header: false,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub regions: usize,
    pub inputs: usize,
}

struct Bed {
    header: Option<Vec<u8>>,
    regions: Vec<Region>,
}

struct Region {
    line: Vec<u8>,
    name: Vec<u8>,
    start: u64,
    end: u64,
    fields: usize,
}

struct Source {
    path: PathBuf,
    reader: input::Reader,
    header: sam::Header,
}

#[derive(Default)]
struct RegionStats {
    coverage: u64,
    threshold_bases: u64,
    reads: u64,
}

pub fn write(
    bed_path: &Path,
    inputs: &[PathBuf],
    options: Options<'_>,
    output: impl Write,
) -> Result<Summary> {
    if inputs.is_empty() {
        return Err(RsomicsError::ConfigError(
            "bedcov requires at least one alignment input".to_owned(),
        ));
    }
    if let Some(indices) = options.indices
        && indices.len() != inputs.len()
    {
        return Err(RsomicsError::ConfigError(
            "bedcov requires one custom index per alignment input".to_owned(),
        ));
    }

    let bed = Bed::read(bed_path)?;
    let mut sources = open_sources(inputs, options)?;
    let dictionary = reference_dictionary(&sources[0].header);
    for source in sources.iter().skip(1) {
        if reference_dictionary(&source.header) != dictionary {
            return Err(RsomicsError::InvalidInput(format!(
                "reference dictionary in {} differs from {}",
                source.path.display(),
                sources[0].path.display()
            )));
        }
    }
    let sweep = bed.regions.len() >= 256
        && options.depth_threshold.is_none()
        && !options.read_count
        && options.maximum_depth == i32::MAX as usize
        && sources
            .iter()
            .all(|source| source.reader.format() == input::Format::Bam);
    let rows = if sweep {
        collect_sweep(inputs, &bed.regions, options)?
    } else if span_mode(options) {
        collect_indexed_spans(inputs, &bed.regions, options)?
    } else if options.additional_threads == 0 {
        bed.regions
            .iter()
            .map(|region| collect_row(&mut sources, region, options))
            .collect::<Result<Vec<_>>>()?
    } else {
        let workers = options
            .additional_threads
            .checked_add(1)
            .ok_or_else(|| RsomicsError::ConfigError("bedcov thread count overflows".to_owned()))?;
        rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .map_err(|error| {
                RsomicsError::ConfigError(format!("building bedcov worker pool: {error}"))
            })?
            .install(|| {
                bed.regions
                    .par_iter()
                    .map_init(
                        || open_sources(inputs, options),
                        |sources, region| match sources {
                            Ok(sources) => collect_row(sources, region, options),
                            Err(error) => Err(RsomicsError::InvalidInput(format!(
                                "opening bedcov worker inputs: {error}"
                            ))),
                        },
                    )
                    .collect::<Result<Vec<_>>>()
            })?
    };

    let mut output = BufWriter::new(output);
    let mut emitted_header = false;
    if options.header
        && let Some(header) = bed.header.as_deref()
    {
        write_header(&mut output, header, inputs, options)?;
        emitted_header = true;
    }
    for (region, stats) in bed.regions.iter().zip(rows) {
        if options.header && !emitted_header {
            let header = generated_header(region.fields);
            write_header(&mut output, &header, inputs, options)?;
            emitted_header = true;
        }

        output.write_all(&region.line).map_err(RsomicsError::Io)?;
        for stat in &stats {
            write!(output, "\t{}", stat.coverage).map_err(RsomicsError::Io)?;
        }
        if options.depth_threshold.is_some() {
            for stat in &stats {
                write!(output, "\t{}", stat.threshold_bases).map_err(RsomicsError::Io)?;
            }
        }
        if options.read_count {
            for stat in &stats {
                write!(output, "\t{}", stat.reads).map_err(RsomicsError::Io)?;
            }
        }
        writeln!(output).map_err(RsomicsError::Io)?;
    }
    output.flush().map_err(RsomicsError::Io)?;
    Ok(Summary {
        regions: bed.regions.len(),
        inputs: inputs.len(),
    })
}

fn collect_sweep(
    inputs: &[PathBuf],
    regions: &[Region],
    options: Options<'_>,
) -> Result<Vec<Vec<RegionStats>>> {
    let mut rows = (0..regions.len())
        .map(|_| {
            (0..inputs.len())
                .map(|_| RegionStats::default())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    for (source, path) in inputs.iter().enumerate() {
        let mut reader = input::open(path, options.reference, options.additional_threads)?;
        let header = reader.read_header(path)?;
        let names = header
            .reference_sequences()
            .keys()
            .enumerate()
            .map(|(reference_id, name)| (name.as_ref(), reference_id))
            .collect::<HashMap<_, _>>();
        let mut by_reference: HashMap<usize, Vec<(u64, u64, usize)>> = HashMap::new();
        for (region_index, region) in regions.iter().enumerate() {
            let reference_id = names.get(region.name.as_slice()).copied().ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "BED reference is absent from the alignment header: {}",
                    String::from_utf8_lossy(&region.name)
                ))
            })?;
            by_reference.entry(reference_id).or_default().push((
                region.start,
                region.end,
                region_index,
            ));
        }
        for regions in by_reference.values_mut() {
            regions.sort_unstable_by_key(|&(start, end, index)| (start, end, index));
        }

        let mut cigar = Vec::with_capacity(16);
        let mut previous = None;
        let mut current_reference = None;
        let mut cursor = 0usize;
        reader.visit_owned_raw_records(&header, path, |record| {
            let reference_id = record.reference_sequence_id();
            let start = record.alignment_start();
            if reference_id < 0 || start < 0 {
                return Ok(true);
            }
            let coordinate = (reference_id, start);
            if previous.is_some_and(|previous| coordinate < previous) {
                return Err(RsomicsError::InvalidInput(format!(
                    "alignment input is not coordinate sorted: {}",
                    path.display()
                )));
            }
            previous = Some(coordinate);
            let reference_id = usize::try_from(reference_id).unwrap();
            if current_reference != Some(reference_id) {
                current_reference = Some(reference_id);
                cursor = 0;
            }
            if record.flags() & options.excluded_flags != 0
                || record.mapping_quality() < options.minimum_mapping_quality
            {
                return Ok(true);
            }
            let Some(reference_regions) = by_reference.get(&reference_id) else {
                return Ok(true);
            };
            let start = u64::try_from(start).unwrap();
            while cursor < reference_regions.len() && reference_regions[cursor].1 <= start {
                cursor += 1;
            }
            record.decode_cigar_into(&mut cigar)?;
            let mut reference_position = start;
            if options.skip_deletions_and_skips {
                for &(kind, length) in &cigar {
                    let length = u64::from(length);
                    if matches!(kind, 0 | 7 | 8) {
                        let end = reference_position.checked_add(length).ok_or_else(|| {
                            RsomicsError::InvalidInput(
                                "alignment reference span overflows".to_owned(),
                            )
                        })?;
                        add_span(
                            reference_position..end,
                            reference_regions,
                            cursor,
                            &mut rows,
                            source,
                        );
                        reference_position = end;
                    } else if matches!(kind, 2 | 3) {
                        reference_position =
                            reference_position.checked_add(length).ok_or_else(|| {
                                RsomicsError::InvalidInput(
                                    "alignment reference span overflows".to_owned(),
                                )
                            })?;
                    }
                }
            } else {
                let span = cigar.iter().try_fold(0u64, |span, &(kind, length)| {
                    if matches!(kind, 0 | 2 | 3 | 7 | 8) {
                        span.checked_add(u64::from(length))
                    } else {
                        Some(span)
                    }
                });
                let end = start
                    .checked_add(span.ok_or_else(|| {
                        RsomicsError::InvalidInput("alignment reference span overflows".to_owned())
                    })?)
                    .ok_or_else(|| {
                        RsomicsError::InvalidInput("alignment reference span overflows".to_owned())
                    })?;
                add_span(start..end, reference_regions, cursor, &mut rows, source);
            }
            Ok(true)
        })?;
    }
    Ok(rows)
}

fn add_span(
    span: Range<u64>,
    regions: &[(u64, u64, usize)],
    cursor: usize,
    rows: &mut [Vec<RegionStats>],
    source: usize,
) {
    for &(start, end, region_index) in &regions[cursor..] {
        if start >= span.end {
            break;
        }
        rows[region_index][source].coverage +=
            end.min(span.end).saturating_sub(start.max(span.start));
    }
}

fn collect_row(
    sources: &mut [Source],
    region: &Region,
    options: Options<'_>,
) -> Result<Vec<RegionStats>> {
    let reference_id = sources[0]
        .header
        .reference_sequences()
        .get_index_of(region.name.as_slice())
        .ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "BED reference is absent from the alignment header: {}",
                String::from_utf8_lossy(&region.name)
            ))
        })?;
    sources
        .iter_mut()
        .map(|source| collect_region(source, reference_id, region, options))
        .collect()
}

fn span_mode(options: Options<'_>) -> bool {
    options.depth_threshold.is_none()
        && !options.read_count
        && options.maximum_depth == i32::MAX as usize
}

fn collect_indexed_spans(
    inputs: &[PathBuf],
    regions: &[Region],
    options: Options<'_>,
) -> Result<Vec<Vec<RegionStats>>> {
    let mut rows = (0..regions.len())
        .map(|_| {
            (0..inputs.len())
                .map(|_| RegionStats::default())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    for (source, path) in inputs.iter().enumerate() {
        let mut reader = match options.indices {
            Some(indices) => {
                rust_htslib::bam::IndexedReader::from_path_and_index(path, &indices[source])
            }
            None => rust_htslib::bam::IndexedReader::from_path(path),
        }
        .map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "opening indexed alignment {}: {error}",
                path.display()
            ))
        })?;
        if let Some(reference) = options.reference {
            reader.set_reference(reference).map_err(|error| {
                RsomicsError::ConfigError(format!(
                    "attaching reference {} to {}: {error}",
                    reference.display(),
                    path.display()
                ))
            })?;
        }
        if options.additional_threads > 0 {
            reader
                .set_threads(options.additional_threads)
                .map_err(|error| {
                    RsomicsError::ConfigError(format!(
                        "configuring {} bedcov threads for {}: {error}",
                        options.additional_threads,
                        path.display()
                    ))
                })?;
        }
        for (region_index, region) in regions.iter().enumerate() {
            if region.start == region.end {
                continue;
            }
            let reference_id = reader.header().tid(&region.name).ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "BED reference is absent from the alignment header: {}",
                    String::from_utf8_lossy(&region.name)
                ))
            })?;
            let start = i64::try_from(region.start).map_err(|_| {
                RsomicsError::InvalidInput("BED start exceeds the supported range".to_owned())
            })?;
            let end = i64::try_from(region.end).map_err(|_| {
                RsomicsError::InvalidInput("BED end exceeds the supported range".to_owned())
            })?;
            reader.fetch((reference_id, start, end)).map_err(|error| {
                RsomicsError::InvalidInput(format!(
                    "querying {}:{}-{} from {}: {error}",
                    String::from_utf8_lossy(&region.name),
                    region.start,
                    region.end,
                    path.display()
                ))
            })?;
            for result in reader.records() {
                let record = result.map_err(|error| {
                    RsomicsError::InvalidInput(format!(
                        "reading alignment record from {}: {error}",
                        path.display()
                    ))
                })?;
                if record.flags() & options.excluded_flags != 0
                    || record.mapq() < options.minimum_mapping_quality
                {
                    continue;
                }
                let mut reference_position = u64::try_from(record.pos()).map_err(|_| {
                    RsomicsError::InvalidInput("mapped alignment has no position".to_owned())
                })?;
                if options.skip_deletions_and_skips {
                    for operation in record.cigar().iter() {
                        let length = u64::from(operation.len());
                        match operation {
                            rust_htslib::bam::record::Cigar::Match(_)
                            | rust_htslib::bam::record::Cigar::Equal(_)
                            | rust_htslib::bam::record::Cigar::Diff(_) => {
                                let end =
                                    reference_position.checked_add(length).ok_or_else(|| {
                                        RsomicsError::InvalidInput(
                                            "alignment reference span overflows".to_owned(),
                                        )
                                    })?;
                                rows[region_index][source].coverage += end
                                    .min(region.end)
                                    .saturating_sub(reference_position.max(region.start));
                                reference_position = end;
                            }
                            rust_htslib::bam::record::Cigar::Del(_)
                            | rust_htslib::bam::record::Cigar::RefSkip(_) => {
                                reference_position =
                                    reference_position.checked_add(length).ok_or_else(|| {
                                        RsomicsError::InvalidInput(
                                            "alignment reference span overflows".to_owned(),
                                        )
                                    })?;
                            }
                            _ => {}
                        }
                    }
                } else {
                    let span = record.cigar().iter().try_fold(0u64, |span, operation| {
                        if matches!(
                            operation,
                            rust_htslib::bam::record::Cigar::Match(_)
                                | rust_htslib::bam::record::Cigar::Del(_)
                                | rust_htslib::bam::record::Cigar::RefSkip(_)
                                | rust_htslib::bam::record::Cigar::Equal(_)
                                | rust_htslib::bam::record::Cigar::Diff(_)
                        ) {
                            span.checked_add(u64::from(operation.len())).ok_or_else(|| {
                                RsomicsError::InvalidInput(
                                    "alignment reference span overflows".to_owned(),
                                )
                            })
                        } else {
                            Ok(span)
                        }
                    })?;
                    let end = reference_position.checked_add(span).ok_or_else(|| {
                        RsomicsError::InvalidInput("alignment reference span overflows".to_owned())
                    })?;
                    rows[region_index][source].coverage += end
                        .min(region.end)
                        .saturating_sub(reference_position.max(region.start));
                }
            }
        }
    }
    Ok(rows)
}

fn open_sources(inputs: &[PathBuf], options: Options<'_>) -> Result<Vec<Source>> {
    inputs
        .iter()
        .enumerate()
        .map(|(index, path)| {
            if path == Path::new("-") {
                return Err(RsomicsError::ConfigError(
                    "bedcov inputs must be indexed named files".to_owned(),
                ));
            }
            let mut reader = match options.indices {
                Some(indices) => {
                    input::open_indexed_with_index(path, &indices[index], options.reference)?
                }
                None => input::open_indexed(path, options.reference)?,
            };
            let header = reader.read_header(path)?;
            Ok(Source {
                path: path.clone(),
                reader,
                header,
            })
        })
        .collect()
}

fn collect_region(
    source: &mut Source,
    reference_id: usize,
    region: &Region,
    options: Options<'_>,
) -> Result<RegionStats> {
    if region.start == region.end {
        return Ok(RegionStats::default());
    }
    let start =
        Position::try_from(usize::try_from(region.start).unwrap() + 1).map_err(|error| {
            RsomicsError::InvalidInput(format!("BED start exceeds the supported range: {error}"))
        })?;
    let end = Position::try_from(usize::try_from(region.end).unwrap()).map_err(|error| {
        RsomicsError::InvalidInput(format!("BED end exceeds the supported range: {error}"))
    })?;
    let query = QueryRegion::new(region.name.clone(), start..=end);
    let reference_lengths = source
        .header
        .reference_sequences()
        .values()
        .map(|reference| u64::try_from(usize::from(reference.length())).unwrap());
    let maximum_depth = options
        .depth_threshold
        .map_or(options.maximum_depth, |threshold| {
            options.maximum_depth.max(threshold)
        });
    let pileup_options = PileupOptions {
        filter: RecordFilter {
            flags: FlagFilter {
                skip_any_set: options.excluded_flags,
                ..FlagFilter::default()
            },
            minimum_mapping_quality: options.minimum_mapping_quality,
            include_anomalous_pairs: true,
        },
        adjust_overlaps: false,
        maximum_depth_per_source: (maximum_depth != 0).then_some(maximum_depth),
    };
    let mut pileup = PileupEngine::new(reference_lengths, pileup_options);
    let mut stats = RegionStats::default();
    let mut encoder = RawRecordEncoder::new();
    source
        .reader
        .visit_region(&source.header, &source.path, Some(&query), |record| {
            pileup
                .push(encoder.encode(&source.header, record)?)
                .map_err(pileup_error)?;
            drain_region(&mut pileup, reference_id, region, options, &mut stats)?;
            Ok(true)
        })?;
    pileup.finish().map_err(pileup_error)?;
    drain_region(&mut pileup, reference_id, region, options, &mut stats)?;
    Ok(stats)
}

fn drain_region(
    pileup: &mut PileupEngine,
    reference_id: usize,
    region: &Region,
    options: Options<'_>,
    stats: &mut RegionStats,
) -> Result<()> {
    pileup.drain(|column| {
        for entry in column.entries() {
            if entry.projection().is_head {
                stats.reads += 1;
            }
        }
        let position = u64::try_from(column.position()).unwrap();
        if usize::try_from(column.reference_id()).ok() != Some(reference_id)
            || position < region.start
            || position >= region.end
        {
            return Ok(());
        }
        let exclude_gaps = options.skip_deletions_and_skips || options.depth_threshold.is_some();
        let depth = if exclude_gaps {
            column
                .entries()
                .filter(|entry| {
                    let projection = entry.projection();
                    !projection.is_deletion && !projection.is_reference_skip
                })
                .count()
        } else {
            column.len()
        };
        stats.coverage += u64::try_from(depth).unwrap();
        if options
            .depth_threshold
            .is_some_and(|threshold| depth >= threshold)
        {
            stats.threshold_bases += 1;
        }
        Ok::<_, RsomicsError>(())
    })?;
    Ok(())
}

fn write_header(
    output: &mut impl Write,
    bed_header: &[u8],
    inputs: &[PathBuf],
    options: Options<'_>,
) -> Result<()> {
    output.write_all(bed_header).map_err(RsomicsError::Io)?;
    for input in inputs {
        write!(output, "\t{}_cov", input.display()).map_err(RsomicsError::Io)?;
    }
    if options.depth_threshold.is_some() {
        for input in inputs {
            write!(output, "\t{}_depth", input.display()).map_err(RsomicsError::Io)?;
        }
    }
    if options.read_count {
        for input in inputs {
            write!(output, "\t{}_count", input.display()).map_err(RsomicsError::Io)?;
        }
    }
    writeln!(output).map_err(RsomicsError::Io)
}

fn generated_header(fields: usize) -> Vec<u8> {
    const NAMES: [&str; 12] = [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRgb",
        "blockCount",
        "blockSizes",
        "blockStarts",
    ];
    let mut header = Vec::new();
    for index in 0..fields {
        header.extend_from_slice(if index == 0 { b"#" } else { b"\t" });
        header.extend_from_slice(NAMES.get(index).copied().unwrap_or(".").as_bytes());
    }
    header
}

impl Bed {
    fn read(path: &Path) -> Result<Self> {
        let mut source = File::open(path).map_err(RsomicsError::Io)?;
        let mut magic = [0; 2];
        let bytes = source.read(&mut magic).map_err(RsomicsError::Io)?;
        drop(source);
        let source = File::open(path).map_err(RsomicsError::Io)?;
        let reader: Box<dyn BufRead> = if bytes == 2 && magic == [0x1f, 0x8b] {
            Box::new(BufReader::new(MultiGzDecoder::new(source)))
        } else {
            Box::new(BufReader::new(source))
        };
        let mut header = None;
        let mut regions = Vec::new();
        for (index, result) in reader.split(b'\n').enumerate() {
            let mut line = result.map_err(RsomicsError::Io)?;
            if line.last() == Some(&b'\r') {
                line.pop();
            }
            if line.is_empty() {
                continue;
            }
            if line.starts_with(b"#") {
                if line.starts_with(b"#chrom\t") && header.is_none() {
                    header = Some(line);
                }
                continue;
            }
            if line.starts_with(b"track ") || line.starts_with(b"browser ") {
                continue;
            }
            let fields = line.iter().filter(|&&byte| byte == b'\t').count() + 1;
            let mut values = line.split(|byte| byte.is_ascii_whitespace());
            let name = values.next().unwrap().to_vec();
            let start = parse_coordinate(path, index + 1, values.next(), "start")?;
            let end = parse_coordinate(path, index + 1, values.next(), "end")?;
            if end < start {
                return Err(invalid_bed(path, index + 1, "end precedes start"));
            }
            regions.push(Region {
                line,
                name,
                start,
                end,
                fields,
            });
        }
        Ok(Self { header, regions })
    }
}

fn reference_dictionary(header: &sam::Header) -> Vec<(Vec<u8>, u64)> {
    header
        .reference_sequences()
        .iter()
        .map(|(name, reference)| {
            (
                name.to_vec(),
                u64::try_from(usize::from(reference.length())).unwrap(),
            )
        })
        .collect()
}

fn parse_coordinate(path: &Path, line: usize, value: Option<&[u8]>, field: &str) -> Result<u64> {
    let value = value.ok_or_else(|| invalid_bed(path, line, &format!("missing {field}")))?;
    std::str::from_utf8(value)
        .ok()
        .and_then(|value| value.parse().ok())
        .ok_or_else(|| invalid_bed(path, line, &format!("invalid {field}")))
}

fn invalid_bed(path: &Path, line: usize, reason: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{}:{line}: {reason}", path.display()))
}

fn pileup_error(error: PileupError) -> RsomicsError {
    RsomicsError::InvalidInput(format!("building BED coverage pileup: {error}"))
}
