use std::collections::{BTreeMap, HashMap};
use std::ops::Range;

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;
use smallvec::SmallVec;

use super::barcode::BarcodeStats;
use super::checksum::Checksums;
use super::coverage::Coverage;
use super::record_data::RecordData;
use super::ref_stats::ReferenceStats;
use super::reference::Slice as ReferenceSlice;
use super::regions::Regions;
use super::{CoverageBins, Options};

const PAIRED: u16 = 0x001;
const PROPER_PAIR: u16 = 0x002;
const UNMAPPED: u16 = 0x004;
const MATE_UNMAPPED: u16 = 0x008;
const REVERSE: u16 = 0x010;
const MATE_REVERSE: u16 = 0x020;
const READ1: u16 = 0x040;
const READ2: u16 = 0x080;
const SECONDARY: u16 = 0x100;
const QC_FAIL: u16 = 0x200;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize)]
pub(crate) struct BaseCounts {
    pub(crate) a: u64,
    pub(crate) c: u64,
    pub(crate) g: u64,
    pub(crate) t: u64,
    pub(crate) n: u64,
    pub(crate) other: u64,
}

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub(crate) struct GcDepth {
    pub(crate) gc: f32,
    pub(crate) depth: u32,
}

#[derive(Clone, Debug, Default, Serialize)]
#[serde(transparent)]
pub(crate) struct QualityCycles(Vec<SmallVec<[(u16, u64); 1]>>);

impl QualityCycles {
    pub(crate) fn ensure_length(&mut self, length: usize) -> Result<()> {
        if self.0.len() < length {
            self.0.try_reserve(length - self.0.len()).map_err(|_| {
                RsomicsError::InvalidInput(
                    "quality-cycle count exceeds available memory".to_owned(),
                )
            })?;
            self.0.resize_with(length, SmallVec::new);
        }
        Ok(())
    }

    pub(crate) fn increment(&mut self, cycle: usize, quality: usize) {
        let row = &mut self.0[cycle];
        if let Some((_, count)) = row
            .iter_mut()
            .find(|(value, _)| usize::from(*value) == quality)
        {
            *count += 1;
        } else {
            row.push((quality as u16, 1));
        }
    }

    pub(crate) fn get(&self, cycle: usize, quality: usize) -> u64 {
        self.0
            .get(cycle)
            .and_then(|row| row.iter().find(|(value, _)| usize::from(*value) == quality))
            .map_or(0, |(_, count)| *count)
    }

    pub(crate) fn len(&self) -> usize {
        self.0.len()
    }
}

impl BaseCounts {
    pub(crate) fn acgt(self) -> u64 {
        self.a + self.c + self.g + self.t
    }

    pub(crate) fn add(self, other: Self) -> Self {
        Self {
            a: self.a + other.a,
            c: self.c + other.c,
            g: self.g + other.g,
            t: self.t + other.t,
            n: self.n + other.n,
            other: self.other + other.other,
        }
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub(crate) struct Summary {
    pub(crate) filtered: u64,
    pub(crate) first: u64,
    pub(crate) last: u64,
    pub(crate) other: u64,
    pub(crate) mapped_paired: u64,
    pub(crate) mapped_single: u64,
    pub(crate) unmapped: u64,
    pub(crate) properly_paired: u64,
    pub(crate) paired: u64,
    pub(crate) duplicated: u64,
    pub(crate) mq0: u64,
    pub(crate) qc_failed: u64,
    pub(crate) secondary: u64,
    pub(crate) supplementary: u64,
    pub(crate) total_length: u64,
    pub(crate) first_length: u64,
    pub(crate) last_length: u64,
    pub(crate) mapped_bases: u64,
    pub(crate) cigar_bases: u64,
    pub(crate) trimmed_bases: u64,
    pub(crate) duplicated_bases: u64,
    pub(crate) mismatches: u64,
    pub(crate) quality_sum: f64,
    pub(crate) anomalous: u64,
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct Accumulator {
    pub(crate) checksum: Checksums,
    pub(crate) summary: Summary,
    pub(crate) sorted: bool,
    pub(crate) max_length: usize,
    pub(crate) max_first_length: usize,
    pub(crate) max_last_length: usize,
    pub(crate) max_quality: usize,
    pub(crate) first_qualities: QualityCycles,
    pub(crate) last_qualities: QualityCycles,
    pub(crate) mismatch_cycles: Option<QualityCycles>,
    pub(crate) first_bases: Vec<BaseCounts>,
    pub(crate) last_bases: Vec<BaseCounts>,
    pub(crate) oriented_bases: Vec<BaseCounts>,
    pub(crate) first_gc: Vec<u64>,
    pub(crate) last_gc: Vec<u64>,
    pub(crate) insert_sizes: BTreeMap<usize, [u64; 3]>,
    pub(crate) read_lengths: BTreeMap<usize, u64>,
    pub(crate) first_lengths: BTreeMap<usize, u64>,
    pub(crate) last_lengths: BTreeMap<usize, u64>,
    pub(crate) mapping_qualities: Vec<u64>,
    pub(crate) insertions: BTreeMap<usize, u64>,
    pub(crate) deletions: BTreeMap<usize, u64>,
    pub(crate) insertion_cycles: Vec<[u64; 2]>,
    pub(crate) deletion_cycles: Vec<[u64; 2]>,
    #[serde(skip)]
    pub(crate) coverage: Coverage,
    pub(crate) barcodes: Vec<BarcodeStats>,
    pub(crate) target_bases: Option<u64>,
    pub(crate) reference_stats: Option<ReferenceStats>,
    pub(crate) gc_depth: Vec<GcDepth>,
    #[serde(skip)]
    pairs: HashMap<Vec<u8>, PairCoverage>,
    #[serde(skip)]
    pair_expiry: BTreeMap<i64, Vec<Vec<u8>>>,
    #[serde(skip)]
    previous_coordinate: Option<(i32, i64)>,
    #[serde(skip)]
    gc_reference: i32,
    #[serde(skip)]
    gc_position: i64,
}

impl Accumulator {
    pub(crate) fn new(reference: bool, target_bases: Option<u64>) -> Self {
        Self {
            checksum: Checksums::default(),
            summary: Summary::default(),
            sorted: true,
            max_length: 0,
            max_first_length: 0,
            max_last_length: 0,
            max_quality: 0,
            first_qualities: QualityCycles::default(),
            last_qualities: QualityCycles::default(),
            mismatch_cycles: reference.then(QualityCycles::default),
            first_bases: Vec::new(),
            last_bases: Vec::new(),
            oriented_bases: Vec::new(),
            first_gc: vec![0; 200],
            last_gc: vec![0; 200],
            insert_sizes: BTreeMap::new(),
            read_lengths: BTreeMap::new(),
            first_lengths: BTreeMap::new(),
            last_lengths: BTreeMap::new(),
            mapping_qualities: vec![0; 256],
            insertions: BTreeMap::new(),
            deletions: BTreeMap::new(),
            insertion_cycles: Vec::new(),
            deletion_cycles: Vec::new(),
            coverage: Coverage::default(),
            barcodes: [
                (*b"BC", *b"QT"),
                (*b"CR", *b"CY"),
                (*b"OX", *b"BZ"),
                (*b"RX", *b"QX"),
            ]
            .into_iter()
            .map(|(sequence, quality)| BarcodeStats::new(sequence, quality))
            .collect(),
            target_bases,
            reference_stats: None,
            gc_depth: vec![GcDepth::default()],
            pairs: HashMap::new(),
            pair_expiry: BTreeMap::new(),
            previous_coordinate: None,
            gc_reference: -1,
            gc_position: -1,
        }
    }

    pub(crate) fn collect(
        &mut self,
        record: &RecordData,
        reference: Option<ReferenceSlice<'_>>,
        regions: Option<&Regions>,
        options: Options<'_>,
    ) -> Result<()> {
        let flags = record.flags;
        if flags & options.required_flags != options.required_flags
            || flags & options.filtered_flags != 0
        {
            self.summary.filtered += 1;
            return Ok(());
        }
        if options
            .read_length
            .is_some_and(|length| length != record.sequence.len())
        {
            return Ok(());
        }
        self.checksum.update(record);
        if flags & SECONDARY != 0 {
            self.summary.secondary += 1;
            return Ok(());
        }
        if flags & SUPPLEMENTARY != 0 {
            self.summary.supplementary += 1;
        }
        let sequence_length = record.sequence.len();
        if sequence_length == 0 {
            return Ok(());
        }
        if flags & DUPLICATE != 0 {
            self.summary.duplicated += 1;
            self.summary.duplicated_bases += sequence_length as u64;
        }
        let order = read_order(flags);
        let cigar = &record.cigar;
        let hard_clipped = cigar.iter().filter(|&&(kind, _)| kind == 5).try_fold(
            sequence_length,
            |length, &(_, count)| {
                length
                    .checked_add(count as usize)
                    .ok_or_else(|| RsomicsError::InvalidInput("read length overflows".to_owned()))
            },
        )?;
        self.ensure_cycles(hard_clipped, order)?;
        self.max_length = self.max_length.max(hard_clipped);
        if order == 1 {
            self.max_first_length = self.max_first_length.max(hard_clipped);
        } else if order == 2 {
            self.max_last_length = self.max_last_length.max(hard_clipped);
        }
        if flags & (UNMAPPED | SECONDARY | SUPPLEMENTARY | QC_FAIL | DUPLICATE) == 0 {
            self.mapping_qualities[usize::from(record.mapping_quality)] += 1;
        }
        if flags & SUPPLEMENTARY == 0 {
            *self.read_lengths.entry(hard_clipped).or_default() += 1;
            if order == 1 {
                *self.first_lengths.entry(hard_clipped).or_default() += 1;
            } else if order == 2 {
                *self.last_lengths.entry(hard_clipped).or_default() += 1;
            }
            self.collect_original(record, flags, order, options)?;
        }
        if flags & UNMAPPED != 0 {
            return Ok(());
        }
        self.collect_indels(record, flags, order, hard_clipped, cigar)?;
        self.collect_insert(record, flags, options.maximum_insert_size);
        self.summary.mismatches += record.edit_distance.unwrap_or(0);
        self.collect_alignment(record, cigar, regions, options.remove_overlaps)?;
        self.collect_mismatches(record, flags, hard_clipped, cigar, reference)?;
        self.collect_gc_depth(record, reference, options.gc_depth)?;
        Ok(())
    }

    fn collect_original(
        &mut self,
        record: &RecordData,
        flags: u16,
        order: usize,
        options: Options<'_>,
    ) -> Result<()> {
        let length = record.sequence.len();
        self.summary.total_length += length as u64;
        self.summary.qc_failed += u64::from(flags & QC_FAIL != 0);
        self.summary.paired += u64::from(flags & PAIRED != 0);
        let reverse = flags & REVERSE != 0;
        let (gc, maximum_quality, quality_sum) = match order {
            1 => collect_cycle_stats(
                &record.sequence,
                &record.qualities,
                reverse,
                &mut self.first_bases,
                &mut self.oriented_bases,
                &mut self.first_qualities,
            ),
            2 => collect_cycle_stats(
                &record.sequence,
                &record.qualities,
                reverse,
                &mut self.last_bases,
                &mut self.oriented_bases,
                &mut self.last_qualities,
            ),
            _ => (0, 0, 0.0),
        };
        self.max_quality = self.max_quality.max(maximum_quality);
        self.summary.quality_sum += quality_sum;
        if order == 1 || order == 2 {
            let start = gc * 199 / length;
            let end = ((gc + 1) * 199 / length).min(199);
            let histogram = if order == 1 {
                &mut self.first_gc
            } else {
                &mut self.last_gc
            };
            for value in &mut histogram[start..end] {
                *value += 1;
            }
        }
        if order == 1 {
            for index in 0..self.barcodes.len() {
                self.barcodes[index].collect(
                    &record.name,
                    record.barcodes[index].as_deref(),
                    record.barcode_qualities[index].as_deref(),
                )?;
            }
        }
        match order {
            1 => {
                self.summary.first += 1;
                self.summary.first_length += length as u64;
            }
            2 => {
                self.summary.last += 1;
                self.summary.last_length += length as u64;
            }
            _ => self.summary.other += 1,
        }
        self.summary.trimmed_bases +=
            bwa_trimmed(&record.qualities, length, reverse, options.trim_quality) as u64;
        if flags & UNMAPPED != 0 {
            self.summary.unmapped += 1;
        } else {
            self.summary.mapped_bases += length as u64;
            self.summary.mq0 += u64::from(record.mapping_quality == 0);
            if flags & PAIRED != 0 && flags & MATE_UNMAPPED == 0 {
                self.summary.mapped_paired += 1;
                self.summary.properly_paired += u64::from(flags & PROPER_PAIR != 0);
                self.summary.anomalous += u64::from(record.reference != record.mate_reference);
            } else {
                self.summary.mapped_single += 1;
            }
        }
        Ok(())
    }

    fn collect_indels(
        &mut self,
        record: &RecordData,
        flags: u16,
        order: usize,
        read_length: usize,
        cigar: &[(u8, u32)],
    ) -> Result<()> {
        let forward = flags & REVERSE == 0;
        let mut cycle = 0usize;
        for &(kind, raw_count) in cigar {
            let count = raw_count as usize;
            if count == 0 {
                continue;
            }
            match kind {
                1 => {
                    let index = if forward {
                        cycle
                    } else {
                        read_length.checked_sub(cycle + count).ok_or_else(|| {
                            RsomicsError::InvalidInput("insertion cycle exceeds read".to_owned())
                        })?
                    };
                    if order == 1 || order == 2 {
                        self.insertion_cycles[index][order - 1] += 1;
                    }
                    cycle += count;
                    if count <= 300 {
                        *self.insertions.entry(count).or_default() += 1;
                    }
                }
                2 => {
                    let index = if forward {
                        cycle.checked_sub(1)
                    } else {
                        read_length.checked_sub(cycle + 1)
                    };
                    if let Some(index) = index
                        && (order == 1 || order == 2)
                    {
                        self.deletion_cycles[index][order - 1] += 1;
                    }
                    if count <= 300 {
                        *self.deletions.entry(count).or_default() += 1;
                    }
                }
                3 | 5 | 6 => {}
                _ => cycle += count,
            }
        }
        let _ = record;
        Ok(())
    }

    fn collect_insert(&mut self, record: &RecordData, flags: u16, maximum: usize) {
        if flags & (PAIRED | MATE_UNMAPPED | SECONDARY | SUPPLEMENTARY) != PAIRED
            || flags & UNMAPPED != 0
        {
            return;
        }
        let size = usize::try_from(record.template_length.unsigned_abs()).unwrap_or(usize::MAX);
        let size = if maximum == 0 {
            size
        } else {
            size.min(maximum)
        };
        if size == 0 && record.reference != record.mate_reference {
            return;
        }
        let first = if flags & READ1 != 0 { 1 } else { -1 };
        let forward = if flags & REVERSE == 0 { 1 } else { -1 };
        let mate_forward = if flags & MATE_REVERSE == 0 { 1 } else { -1 };
        let delta = record.mate_position - record.position;
        let orientation = if forward * mate_forward > 0 {
            2
        } else if first * delta > 0 {
            usize::from(first * forward <= 0)
        } else if first * delta < 0 {
            usize::from(first * forward > 0)
        } else {
            0
        };
        self.insert_sizes.entry(size).or_default()[orientation] += 1;
    }

    fn collect_alignment(
        &mut self,
        record: &RecordData,
        cigar: &[(u8, u32)],
        regions: Option<&Regions>,
        remove_overlaps: bool,
    ) -> Result<()> {
        let active_region = if let Some(regions) = regions {
            regions.first_overlap(record.reference, record.position, record.reference_end()?)
        } else {
            None
        };
        let coordinate = (record.reference, record.position);
        if self
            .previous_coordinate
            .is_some_and(|previous| previous.0 == coordinate.0 && coordinate.1 < previous.1)
        {
            if regions.is_some() {
                return Err(RsomicsError::InvalidInput(
                    "target-region statistics require coordinate-sorted input".to_owned(),
                ));
            }
            self.sorted = false;
        }
        self.previous_coordinate = Some(coordinate);
        if self.sorted {
            self.coverage.advance(record.reference, record.position);
        }
        let mut position = record.position;
        let mut coverage_chunks = Vec::new();
        for &(kind, count) in cigar {
            let end = position + i64::from(count);
            if matches!(kind, 0 | 7 | 8) {
                if let Some(regions) = regions {
                    for range in regions.intersections(record.reference, position, end) {
                        coverage_chunks.push(range);
                    }
                    if let Some(active) = &active_region {
                        self.summary.cigar_bases += u64::try_from(
                            active
                                .end
                                .min(end)
                                .saturating_sub(active.start.max(position)),
                        )
                        .unwrap();
                    }
                } else {
                    self.summary.cigar_bases += u64::from(count);
                    coverage_chunks.push(position..end);
                }
            } else if kind == 1
                && (regions.is_none()
                    || active_region
                        .as_ref()
                        .is_some_and(|range| range.contains(&position)))
            {
                self.summary.cigar_bases += u64::from(count);
            }
            if matches!(kind, 0 | 2 | 3 | 7 | 8) {
                position = end;
            }
        }
        if !self.sorted {
            return Ok(());
        }
        if remove_overlaps {
            let overlap = self.insert_without_mate_overlap(record, coverage_chunks)?;
            self.summary.cigar_bases =
                self.summary
                    .cigar_bases
                    .checked_sub(overlap)
                    .ok_or_else(|| {
                        RsomicsError::InvalidInput(
                            "overlapping mate bases exceed mapped CIGAR bases".to_owned(),
                        )
                    })?;
        } else {
            self.insert_coverage_chunks(&coverage_chunks);
        }
        Ok(())
    }

    fn insert_without_mate_overlap(
        &mut self,
        record: &RecordData,
        chunks: Vec<Range<i64>>,
    ) -> Result<u64> {
        self.expire_pairs(record.position);
        let order = read_order(record.flags);
        let eligible = record.flags & PAIRED != 0
            && record.flags & MATE_UNMAPPED == 0
            && record.template_length.unsigned_abs() < 2 * record.sequence.len() as u64
            && matches!(order, 1 | 2);
        if !eligible {
            self.insert_coverage_chunks(&chunks);
            return Ok(0);
        }

        if let Some(pair) = self.pairs.get(&record.name)
            && pair.order != order
        {
            let stored = pair.chunks.clone();
            let mut overlap = 0u64;
            for chunk in chunks {
                let (kept, removed) = subtract_ranges(chunk, &stored);
                overlap += removed;
                self.insert_coverage_chunks(&kept);
            }
            self.pairs.remove(&record.name);
            return Ok(overlap);
        }

        self.insert_coverage_chunks(&chunks);
        let pair = self
            .pairs
            .entry(record.name.clone())
            .or_insert_with(|| PairCoverage {
                order,
                chunks: Vec::new(),
                maximum_end: i64::MIN,
            });
        pair.maximum_end = pair.maximum_end.max(
            chunks
                .iter()
                .map(|range| range.end)
                .max()
                .unwrap_or(i64::MIN),
        );
        pair.chunks.extend(chunks);
        pair.chunks
            .sort_unstable_by_key(|range| (range.start, range.end));
        self.pair_expiry
            .entry(pair.maximum_end)
            .or_default()
            .push(record.name.clone());
        Ok(0)
    }

    fn insert_coverage_chunks(&mut self, chunks: &[Range<i64>]) {
        for range in chunks {
            self.coverage.insert(range.start, range.end);
        }
    }

    fn expire_pairs(&mut self, position: i64) {
        while self
            .pair_expiry
            .first_key_value()
            .is_some_and(|(&end, _)| end < position)
        {
            let (end, names) = self.pair_expiry.pop_first().unwrap();
            for name in names {
                if self
                    .pairs
                    .get(&name)
                    .is_some_and(|pair| pair.maximum_end == end)
                {
                    self.pairs.remove(&name);
                }
            }
        }
    }

    fn collect_mismatches(
        &mut self,
        record: &RecordData,
        flags: u16,
        read_length: usize,
        cigar: &[(u8, u32)],
        reference: Option<ReferenceSlice<'_>>,
    ) -> Result<()> {
        let Some(reference) = reference else {
            return Ok(());
        };
        let reverse = flags & REVERSE != 0;
        let mut query = 0usize;
        let mut cycle = 0usize;
        let mut reference_index = usize::try_from(record.position).map_err(|_| {
            RsomicsError::InvalidInput("mapped read has a negative position".to_owned())
        })?;
        for &(kind, raw_count) in cigar {
            let count = raw_count as usize;
            match kind {
                1 => {
                    query += count;
                    cycle += count;
                }
                2 => reference_index += count,
                4 => {
                    query += count;
                    cycle += count;
                }
                5 => cycle += count,
                3 | 6 => {}
                0 | 7 | 8 => {
                    let local_reference = reference_index.checked_sub(reference.start);
                    if local_reference
                        .and_then(|start| start.checked_add(count))
                        .is_none_or(|end| end > reference.sequence.len())
                        || query + count > record.sequence.len()
                    {
                        return Err(RsomicsError::InvalidInput(
                            "alignment extends beyond its sequence or reference".to_owned(),
                        ));
                    }
                    let local_reference = local_reference.unwrap();
                    for offset in 0..count {
                        let base = record.sequence[query + offset];
                        let reference_base =
                            reference_code(reference.sequence[local_reference + offset]);
                        let output_cycle = if reverse {
                            read_length - cycle - offset - 1
                        } else {
                            cycle + offset
                        };
                        if base == 15 {
                            self.mismatch_cycles
                                .as_mut()
                                .unwrap()
                                .increment(output_cycle, 0);
                        } else if base != 0 && reference_base != 0 && base != reference_base {
                            let quality = record
                                .qualities
                                .get(query + offset)
                                .copied()
                                .unwrap_or(255)
                                .wrapping_add(1) as usize;
                            self.mismatch_cycles
                                .as_mut()
                                .unwrap()
                                .increment(output_cycle, quality);
                        }
                    }
                    query += count;
                    cycle += count;
                    reference_index += count;
                }
                _ => {
                    return Err(RsomicsError::InvalidInput(format!(
                        "unsupported CIGAR operation {kind}"
                    )));
                }
            }
        }
        Ok(())
    }

    fn collect_gc_depth(
        &mut self,
        record: &RecordData,
        reference: Option<ReferenceSlice<'_>>,
        raw_bin_size: f64,
    ) -> Result<()> {
        if !self.sorted {
            return Ok(());
        }
        let bin_size = raw_bin_size as i64;
        if bin_size <= 0 {
            return Err(RsomicsError::ConfigError(
                "GC-depth bin size must be at least one".to_owned(),
            ));
        }
        let read_span = record
            .cigar
            .iter()
            .filter(|&&(kind, _)| matches!(kind, 0 | 2 | 3 | 7 | 8))
            .try_fold(0i64, |total, &(_, count)| {
                total.checked_add(i64::from(count))
            })
            .ok_or_else(|| {
                RsomicsError::InvalidInput("alignment reference span overflows".to_owned())
            })?;
        let end = record.position.checked_add(read_span).ok_or_else(|| {
            RsomicsError::InvalidInput("alignment end position overflows".to_owned())
        })?;
        let new_bin = self.gc_position < 0
            || self.gc_reference != record.reference
            || if reference.is_some() {
                self.gc_position + bin_size < end
            } else {
                record.position - self.gc_position > bin_size
            };
        if new_bin {
            self.gc_depth.push(GcDepth::default());
            self.gc_reference = record.reference;
            self.gc_position = record.position;
            if let Some(reference) = reference {
                let start = usize::try_from(record.position).map_err(|_| {
                    RsomicsError::InvalidInput("mapped read has a negative position".to_owned())
                })?;
                let local_start = start.checked_sub(reference.start).ok_or_else(|| {
                    RsomicsError::InvalidInput("reference window starts after the read".to_owned())
                })?;
                let local_end = local_start
                    .saturating_add(bin_size as usize)
                    .min(reference.sequence.len());
                let mut gc = 0u64;
                let mut total = 0u64;
                for base in &reference.sequence[local_start..local_end] {
                    match base.to_ascii_uppercase() {
                        b'G' | b'C' => {
                            gc += 1;
                            total += 1;
                        }
                        b'A' | b'T' => total += 1,
                        _ => {}
                    }
                }
                self.gc_depth.last_mut().unwrap().gc = if total == 0 {
                    0.0
                } else {
                    gc as f32 / total as f32
                };
            }
        }
        let bin = self.gc_depth.last_mut().unwrap();
        bin.depth = bin.depth.checked_add(1).ok_or_else(|| {
            RsomicsError::InvalidInput("GC-depth read count overflows".to_owned())
        })?;
        if reference.is_none() {
            let gc = record
                .sequence
                .iter()
                .filter(|&&base| matches!(base, 2 | 4))
                .count();
            bin.gc += gc as f32 / record.sequence.len() as f32;
        }
        Ok(())
    }

    fn ensure_cycles(&mut self, length: usize, order: usize) -> Result<()> {
        if order == 1 {
            self.first_qualities.ensure_length(length)?;
            ensure_vec_length(&mut self.first_bases, length)?;
        } else if order == 2 {
            self.last_qualities.ensure_length(length)?;
            ensure_vec_length(&mut self.last_bases, length)?;
        }
        if matches!(order, 1 | 2) {
            ensure_vec_length(&mut self.oriented_bases, length)?;
        }
        if let Some(cycles) = &mut self.mismatch_cycles {
            cycles.ensure_length(length)?;
        }
        let indel_length = length.saturating_add(1);
        ensure_vec_length(&mut self.insertion_cycles, indel_length)?;
        ensure_vec_length(&mut self.deletion_cycles, indel_length)?;
        Ok(())
    }

    pub(crate) fn coverage_histogram(&self, bins: CoverageBins) -> Vec<u64> {
        self.coverage.histogram(bins)
    }
}

fn ensure_vec_length<T: Default>(values: &mut Vec<T>, length: usize) -> Result<()> {
    if values.len() < length {
        values.try_reserve(length - values.len()).map_err(|_| {
            RsomicsError::InvalidInput("read cycle count exceeds available memory".to_owned())
        })?;
        values.resize_with(length, T::default);
    }
    Ok(())
}

fn collect_cycle_stats(
    sequence: &[u8],
    quality_scores: &[u8],
    reverse: bool,
    bases: &mut [BaseCounts],
    oriented: &mut [BaseCounts],
    qualities: &mut QualityCycles,
) -> (usize, usize, f64) {
    let mut gc = 0;
    let mut maximum_quality = 0;
    let mut quality_sum = 0.0;
    for (index, &base) in sequence.iter().enumerate() {
        let cycle = if reverse {
            sequence.len() - index - 1
        } else {
            index
        };
        increment_base(&mut bases[cycle], base);
        increment_oriented(&mut oriented[cycle], base, reverse);
        gc += usize::from(matches!(base, 2 | 4));
        let quality = usize::from(quality_scores.get(index).copied().unwrap_or(255));
        qualities.increment(cycle, quality);
        maximum_quality = maximum_quality.max(quality);
        quality_sum += quality as f64;
    }
    (gc, maximum_quality, quality_sum)
}

#[derive(Clone, Debug)]
struct PairCoverage {
    order: usize,
    chunks: Vec<Range<i64>>,
    maximum_end: i64,
}

fn subtract_ranges(chunk: Range<i64>, stored: &[Range<i64>]) -> (Vec<Range<i64>>, u64) {
    let mut kept = vec![chunk.clone()];
    for mask in stored {
        let mut next = Vec::new();
        for range in kept {
            if mask.end <= range.start || mask.start >= range.end {
                next.push(range);
                continue;
            }
            if range.start < mask.start {
                next.push(range.start..mask.start);
            }
            if mask.end < range.end {
                next.push(mask.end..range.end);
            }
        }
        kept = next;
    }
    let kept_length = kept
        .iter()
        .map(|range| u64::try_from(range.end - range.start).unwrap())
        .sum::<u64>();
    (
        kept,
        u64::try_from(chunk.end - chunk.start).unwrap() - kept_length,
    )
}

fn read_order(flags: u16) -> usize {
    if flags & PAIRED == 0 {
        1
    } else {
        usize::from(flags & READ1 != 0) + 2 * usize::from(flags & READ2 != 0)
    }
}

fn increment_base(counts: &mut BaseCounts, code: u8) {
    match code {
        1 => counts.a += 1,
        2 => counts.c += 1,
        4 => counts.g += 1,
        8 => counts.t += 1,
        15 => counts.n += 1,
        _ => counts.other += 1,
    }
}

fn increment_oriented(counts: &mut BaseCounts, code: u8, reverse: bool) {
    increment_base(
        counts,
        if reverse {
            match code {
                1 => 8,
                2 => 4,
                4 => 2,
                8 => 1,
                value => value,
            }
        } else {
            code
        },
    );
}

fn reference_code(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        _ => 0,
    }
}

fn bwa_trimmed(qualities: &[u8], length: usize, reverse: bool, threshold: u8) -> usize {
    if threshold == 0 || length < 35 || qualities.is_empty() {
        return 0;
    }
    let mut sum = 0i32;
    let mut maximum = 0i32;
    let mut trimmed = 0;
    for offset in 0..=length - 35 {
        let quality = qualities[if reverse { offset } else { length - offset - 1 }];
        sum += i32::from(threshold) - i32::from(quality);
        if sum < 0 {
            break;
        }
        if sum > maximum {
            maximum = sum;
            trimmed = offset;
        }
    }
    trimmed
}
