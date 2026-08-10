use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use noodles::core::Region;
use noodles::sam;
use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use rsomics_pileup::{PileupEngine, PileupError, PileupOptions};
use serde::Serialize;

use crate::{alignment_stream, coverage_hts, input};

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ReferenceCoverage {
    pub name: String,
    pub start: u64,
    pub end: u64,
    pub reads: u64,
    pub covered_bases: u64,
    pub coverage: f64,
    pub mean_depth: f64,
    pub mean_base_quality: f64,
    pub mean_mapping_quality: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Report {
    pub inputs: usize,
    pub references: Vec<ReferenceCoverage>,
}

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub reference: Option<&'a Path>,
    pub additional_threads: usize,
    pub minimum_read_length: usize,
    pub minimum_mapping_quality: u8,
    pub minimum_base_quality: u8,
    pub required_flags: u16,
    pub excluded_flags: u16,
    pub maximum_depth: usize,
    pub minimum_depth: usize,
    pub region: Option<&'a str>,
}

impl Default for Options<'_> {
    fn default() -> Self {
        Self {
            reference: None,
            additional_threads: 0,
            minimum_read_length: 0,
            minimum_mapping_quality: 0,
            minimum_base_quality: 0,
            required_flags: 0,
            excluded_flags: 0x704,
            maximum_depth: 1_000_000,
            minimum_depth: 1,
            region: None,
        }
    }
}

#[derive(Clone, Eq, PartialEq)]
struct Reference {
    name: Vec<u8>,
    length: u64,
}

#[derive(Default)]
struct Stats {
    reads: u64,
    covered_bases: u64,
    depth_sum: u64,
    base_quality_sum: u64,
    quality_bases: u64,
    mapping_quality_sum: u64,
}

#[derive(Clone, Copy)]
struct Selection {
    reference_id: Option<usize>,
    start: u64,
    end: u64,
}

struct State<'a> {
    options: Options<'a>,
    references: Vec<Reference>,
    selection: Selection,
    stats: Vec<Stats>,
    pileup: PileupEngine,
    cigar: Vec<(u8, u32)>,
    sources: usize,
    quality_available: Vec<bool>,
}

pub fn collect(inputs: &[PathBuf], options: Options<'_>) -> Result<Report> {
    if inputs.len() == 1 && inputs[0] != Path::new("-") && options.region.is_none() {
        let scan = coverage_hts::collect(&inputs[0], options)?;
        let references = scan
            .references
            .into_iter()
            .map(|reference| Reference {
                name: reference.name,
                length: reference.length,
            })
            .collect::<Vec<_>>();
        let stats = scan
            .stats
            .into_iter()
            .map(|stats| Stats {
                reads: stats.reads,
                covered_bases: stats.covered_bases,
                depth_sum: stats.depth_sum,
                base_quality_sum: stats.base_quality_sum,
                quality_bases: stats.quality_bases,
                mapping_quality_sum: stats.mapping_quality_sum,
            })
            .collect::<Vec<_>>();
        return Ok(build_report(
            &references,
            &stats,
            Selection {
                reference_id: None,
                start: 0,
                end: 0,
            },
            1,
        ));
    }
    let region = options
        .region
        .map(str::parse::<Region>)
        .transpose()
        .map_err(|error| RsomicsError::ConfigError(format!("invalid region: {error}")))?;
    let mut state = alignment_stream::merge(
        inputs,
        options.reference,
        options.additional_threads,
        region.as_ref(),
        |headers| State::new(headers, options, region.as_ref(), inputs.len()),
        |state, source, record| state.push(source, record),
    )?;
    state.finish()
}

impl<'a> State<'a> {
    fn new(
        headers: &[alignment_stream::StreamHeader],
        options: Options<'a>,
        region: Option<&Region>,
        sources: usize,
    ) -> Result<Self> {
        let references = reference_dictionary(&headers[0].header);
        for (index, header) in headers.iter().enumerate().skip(1) {
            if reference_dictionary(&header.header) != references {
                return Err(RsomicsError::InvalidInput(format!(
                    "reference dictionary in input {} differs from the first input",
                    index + 1
                )));
            }
        }
        let selection = Selection::resolve(region, &headers[0].header, &references)?;
        let pileup = PileupEngine::new(
            references.iter().map(|reference| reference.length),
            PileupOptions::default(),
        );
        let stats = (0..references.len()).map(|_| Stats::default()).collect();
        Ok(Self {
            options,
            references,
            selection,
            stats,
            pileup,
            cigar: Vec::with_capacity(16),
            sources,
            quality_available: headers
                .iter()
                .map(|header| {
                    header.format != input::Format::Cram || options.minimum_base_quality > 0
                })
                .collect(),
        })
    }

    fn push(&mut self, source: usize, record: RawRecord) -> Result<()> {
        let flags = record.flags();
        if flags & self.options.excluded_flags != 0
            || (self.options.required_flags != 0 && flags & self.options.required_flags == 0)
            || record.mapping_quality() < self.options.minimum_mapping_quality
            || record.reference_sequence_id() < 0
        {
            return Ok(());
        }
        record.decode_cigar_into(&mut self.cigar)?;
        let query_length = self
            .cigar
            .iter()
            .try_fold(0usize, |length, &(kind, count)| {
                if matches!(kind, 0 | 1 | 4 | 7 | 8) {
                    length.checked_add(usize::try_from(count).ok()?)
                } else {
                    Some(length)
                }
            });
        let query_length = query_length.ok_or_else(|| {
            RsomicsError::InvalidInput("alignment query length overflows".to_owned())
        })?;
        if query_length < self.options.minimum_read_length {
            return Ok(());
        }
        let reference_id = usize::try_from(record.reference_sequence_id()).unwrap();
        let stats = self.stats.get_mut(reference_id).ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "alignment reference ID {reference_id} is absent from the header"
            ))
        })?;
        stats.reads += 1;
        stats.mapping_quality_sum += u64::from(record.mapping_quality());
        self.pileup
            .push_with_source(u32::try_from(source).unwrap(), record)
            .map_err(pileup_error)?;
        self.drain()
    }

    fn drain(&mut self) -> Result<()> {
        let selection = self.selection;
        let maximum_depth = self.options.maximum_depth;
        let minimum_depth = self.options.minimum_depth;
        let minimum_base_quality = self.options.minimum_base_quality;
        let sources = self.sources;
        let quality_available = &self.quality_available;
        let stats = &mut self.stats;
        self.pileup.drain(|column| {
            let reference_id = usize::try_from(column.reference_id()).unwrap();
            let position = u64::try_from(column.position()).unwrap();
            if !selection.contains(reference_id, position) {
                return Ok(());
            }
            let mut source_depths = vec![0usize; sources];
            let mut depth = 0u64;
            let mut base_quality_sum = 0u64;
            let mut quality_bases = 0u64;
            for entry in column.entries() {
                let source = usize::try_from(entry.source_id()).unwrap();
                if maximum_depth != 0 && source_depths[source] >= maximum_depth {
                    continue;
                }
                source_depths[source] += 1;
                let projection = entry.projection();
                if projection.is_deletion || projection.is_reference_skip {
                    continue;
                }
                let quality = if quality_available[source] {
                    entry
                        .record()
                        .quality_scores()
                        .get(projection.qpos)
                        .copied()
                        .unwrap_or(u8::MAX)
                } else {
                    u8::MAX
                };
                if quality < minimum_base_quality {
                    continue;
                }
                depth += 1;
                base_quality_sum += u64::from(quality);
                quality_bases += 1;
            }
            if depth >= u64::try_from(minimum_depth).unwrap() {
                let stats = &mut stats[reference_id];
                stats.covered_bases += 1;
                stats.depth_sum += depth;
                stats.base_quality_sum += base_quality_sum;
                stats.quality_bases += quality_bases;
            }
            Ok::<_, RsomicsError>(())
        })?;
        Ok(())
    }

    fn finish(&mut self) -> Result<Report> {
        self.pileup.finish().map_err(pileup_error)?;
        self.drain()?;
        Ok(build_report(
            &self.references,
            &self.stats,
            self.selection,
            self.sources,
        ))
    }
}

fn build_report(
    references: &[Reference],
    stats: &[Stats],
    selection: Selection,
    sources: usize,
) -> Report {
    let references = references
        .iter()
        .enumerate()
        .filter_map(|(reference_id, reference)| {
            selection
                .bounds(reference_id, reference.length)
                .map(|(start, end)| (reference_id, reference, start, end))
        })
        .map(|(reference_id, reference, start, end)| {
            let stats = &stats[reference_id];
            let length = end - start;
            ReferenceCoverage {
                name: String::from_utf8_lossy(&reference.name).into_owned(),
                start: start + 1,
                end,
                reads: stats.reads,
                covered_bases: stats.covered_bases,
                coverage: ratio(stats.covered_bases, length) * 100.0,
                mean_depth: ratio(stats.depth_sum, length),
                mean_base_quality: ratio(stats.base_quality_sum, stats.quality_bases),
                mean_mapping_quality: ratio(stats.mapping_quality_sum, stats.reads),
            }
        })
        .collect();
    Report {
        inputs: sources,
        references,
    }
}

impl Selection {
    fn resolve(
        region: Option<&Region>,
        header: &sam::Header,
        references: &[Reference],
    ) -> Result<Self> {
        let Some(region) = region else {
            return Ok(Self {
                reference_id: None,
                start: 0,
                end: 0,
            });
        };
        let reference_id = header
            .reference_sequences()
            .get_index_of(region.name())
            .ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "region reference is absent from the alignment header: {}",
                    String::from_utf8_lossy(region.name())
                ))
            })?;
        let reference_length = references[reference_id].length;
        let start = region
            .interval()
            .start()
            .map(|position| u64::try_from(usize::from(position) - 1).unwrap())
            .unwrap_or(0);
        let end = region
            .interval()
            .end()
            .map(|position| u64::try_from(usize::from(position)).unwrap())
            .unwrap_or(reference_length)
            .min(reference_length);
        if start >= end {
            return Err(RsomicsError::InvalidInput(format!(
                "region starts outside reference length {reference_length}"
            )));
        }
        Ok(Self {
            reference_id: Some(reference_id),
            start,
            end,
        })
    }

    fn contains(self, reference_id: usize, position: u64) -> bool {
        self.bounds(reference_id, u64::MAX)
            .is_some_and(|(start, end)| start <= position && position < end)
    }

    fn bounds(self, reference_id: usize, length: u64) -> Option<(u64, u64)> {
        match self.reference_id {
            Some(selected) if selected == reference_id => Some((self.start, self.end)),
            Some(_) => None,
            None => Some((0, length)),
        }
    }
}

fn reference_dictionary(header: &sam::Header) -> Vec<Reference> {
    header
        .reference_sequences()
        .iter()
        .map(|(name, reference)| Reference {
            name: name.to_vec(),
            length: u64::try_from(usize::from(reference.length())).unwrap(),
        })
        .collect()
}

fn ratio(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

fn pileup_error(error: PileupError) -> RsomicsError {
    RsomicsError::InvalidInput(format!("building coverage pileup: {error}"))
}

impl Report {
    pub fn write(&self, header: bool, output: impl Write) -> Result<()> {
        let mut output = BufWriter::new(output);
        if header {
            writeln!(
                output,
                "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"
            )
            .map_err(RsomicsError::Io)?;
        }
        for reference in &self.references {
            writeln!(
                output,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                reference.name,
                reference.start,
                reference.end,
                reference.reads,
                reference.covered_bases,
                format_general(reference.coverage, 6),
                format_general(reference.mean_depth, 6),
                format_general(reference.mean_base_quality, 3),
                format_general(reference.mean_mapping_quality, 3),
            )
            .map_err(RsomicsError::Io)?;
        }
        output.flush().map_err(RsomicsError::Io)
    }
}

fn format_general(value: f64, precision: usize) -> String {
    if value == 0.0 {
        return "0".to_owned();
    }
    let exponent = value.abs().log10().floor() as i32;
    if exponent < -4 || exponent >= precision as i32 {
        let mut value = format!("{:.*e}", precision - 1, value);
        if let Some((mantissa, exponent)) = value.split_once('e') {
            let mantissa = mantissa.trim_end_matches('0').trim_end_matches('.');
            let exponent = exponent.parse::<i32>().unwrap();
            value = format!("{mantissa}e{exponent:+03}");
        }
        value
    } else {
        let decimals = (precision as i32 - exponent - 1).max(0) as usize;
        format!("{value:.decimals$}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_owned()
    }
}
