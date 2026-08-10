mod model;
mod output;

use std::collections::{BTreeMap, HashMap};
use std::hash::{BuildHasherDefault, Hasher};
use std::io::Write;
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::alignment_order::validate_coordinate_fields;
use crate::amplicon::PrimerBed;
use crate::input;
use model::Amplicon;

type FastMap<K, V> = HashMap<K, V, BuildHasherDefault<FnvHasher>>;

struct FnvHasher(u64);

impl Default for FnvHasher {
    fn default() -> Self {
        Self(0xcbf2_9ce4_8422_2325)
    }
}

impl Hasher for FnvHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        let mut hash = self.0;
        for byte in bytes {
            hash ^= u64::from(*byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
        self.0 = hash;
    }
}

pub const MAX_DEPTH_LEVELS: usize = 5;

pub const DEFAULT_FLAG_FILTER: u16 = 0x0B04;
pub const DEFAULT_MAX_DELTA: i64 = 30;
pub const DEFAULT_MIN_DEPTH: [u32; MAX_DEPTH_LEVELS] = [1, 0, 0, 0, 0];
pub const DEFAULT_TCOORD_MIN_COUNT: u32 = 10;
pub const DEFAULT_DEPTH_BIN: f64 = 0.01;
pub const MAX_AMP_DEFAULT: usize = 1000;
pub const MAX_AMP_LEN_DEFAULT: usize = 1000;

#[derive(Clone, Debug)]
pub struct Options {
    pub flag_require: u16,
    pub flag_filter: u16,
    pub max_delta: i64,
    pub min_depth: [u32; MAX_DEPTH_LEVELS],
    pub max_amp: usize,
    pub max_amp_len: usize,
    pub tlen_adj: i64,
    pub depth_bin: f64,
    pub tcoord_min_count: u32,
    pub tcoord_bin: i64,
    pub additional_threads: usize,
    pub single_ref: bool,
    pub use_sample_name: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            flag_require: 0,
            flag_filter: DEFAULT_FLAG_FILTER,
            max_delta: DEFAULT_MAX_DELTA,
            min_depth: DEFAULT_MIN_DEPTH,
            max_amp: MAX_AMP_DEFAULT,
            max_amp_len: MAX_AMP_LEN_DEFAULT,
            tlen_adj: 0,
            depth_bin: DEFAULT_DEPTH_BIN,
            tcoord_min_count: DEFAULT_TCOORD_MIN_COUNT,
            tcoord_bin: 1,
            additional_threads: 0,
            single_ref: false,
            use_sample_name: false,
        }
    }
}

struct AmpStats {
    nseq: i64,
    nfiltered: i64,
    nfailprimer: i64,

    max_amp_len: usize,
    ref_len: i64,

    nreads: Vec<i64>,
    nreads2: Vec<i64>,
    nfull_reads: Vec<f64>,
    nrperc: Vec<f64>,
    nrperc2: Vec<f64>,
    nbases: Vec<i64>,
    nbases2: Vec<i64>,
    coverage: Vec<i64>,
    covered_perc: Vec<[f64; MAX_DEPTH_LEVELS]>,
    covered_perc2: Vec<[f64; MAX_DEPTH_LEVELS]>,
    tcoord: Vec<FastMap<u64, u64>>,
    amp_dist: Vec<[i64; 3]>,
    depth_all: Vec<i64>,
    depth_valid: Vec<i64>,
}

impl AmpStats {
    fn new(ref_len: i64, max_amp: usize, max_amp_len: usize) -> Result<Self> {
        let ref_len = usize::try_from(ref_len)
            .map_err(|_| RsomicsError::InvalidInput("reference length exceeds usize".to_owned()))?;
        let coverage_len = max_amp.checked_mul(max_amp_len).ok_or_else(|| {
            RsomicsError::ConfigError("amplicon coverage allocation overflows".to_owned())
        })?;
        Ok(AmpStats {
            nseq: 0,
            nfiltered: 0,
            nfailprimer: 0,
            max_amp_len,
            ref_len: i64::try_from(ref_len).map_err(|_| {
                RsomicsError::InvalidInput("reference length exceeds i64".to_owned())
            })?,
            nreads: vec![0; max_amp],
            nreads2: vec![0; max_amp],
            nfull_reads: vec![0.0; max_amp],
            nrperc: vec![0.0; max_amp],
            nrperc2: vec![0.0; max_amp],
            nbases: vec![0; max_amp],
            nbases2: vec![0; max_amp],
            coverage: vec![0; coverage_len],
            covered_perc: vec![[0.0; MAX_DEPTH_LEVELS]; max_amp],
            covered_perc2: vec![[0.0; MAX_DEPTH_LEVELS]; max_amp],
            tcoord: (0..=max_amp).map(|_| FastMap::default()).collect(),
            amp_dist: vec![[0; 3]; max_amp],
            depth_all: vec![0; ref_len],
            depth_valid: vec![0; ref_len],
        })
    }

    fn reset(&mut self) {
        self.nseq = 0;
        self.nfiltered = 0;
        self.nfailprimer = 0;
        self.nreads.fill(0);
        self.nreads2.fill(0);
        self.nfull_reads.fill(0.0);
        self.nrperc.fill(0.0);
        self.nrperc2.fill(0.0);
        self.nbases.fill(0);
        self.nbases2.fill(0);
        self.coverage.fill(0);
        for x in &mut self.covered_perc {
            *x = [0.0; MAX_DEPTH_LEVELS];
        }
        for x in &mut self.covered_perc2 {
            *x = [0.0; MAX_DEPTH_LEVELS];
        }
        for h in &mut self.tcoord {
            h.clear();
        }
        for x in &mut self.amp_dist {
            *x = [0; 3];
        }
        self.depth_all.fill(0);
        self.depth_valid.fill(0);
    }
}

struct RefData {
    ref_name: String,
    ref_len: i64,
    amps: Vec<Amplicon>,
    lookup: PosLookup,
    local: AmpStats,
    global: AmpStats,
    first_amp_idx: usize,
}

struct PosLookup {
    pos2start: Vec<i32>,
    pos2end: Vec<i32>,
}

impl PosLookup {
    fn build(amps: &[Amplicon], ref_len: i64, max_delta: i64) -> Self {
        let len = ref_len as usize;
        let mut pos2start = vec![-1i32; len];
        let mut pos2end = vec![-1i32; len];

        for (i, amp) in amps.iter().enumerate() {
            for &lpos in &amp.lefts {
                let lo = ((lpos - max_delta).max(1)) as usize;
                let hi = ((lpos + max_delta) as usize).min(len);
                for p in lo..=hi {
                    if p <= len {
                        pos2start[p - 1] = i as i32;
                    }
                }
            }
            for &rpos in &amp.rights {
                let lo = ((rpos - max_delta).max(1)) as usize;
                let hi = ((rpos + max_delta) as usize).min(len);
                for p in lo..=hi {
                    if p <= len {
                        pos2end[p - 1] = i as i32;
                    }
                }
            }
        }
        PosLookup { pos2start, pos2end }
    }

    #[inline]
    fn get_start(&self, pos: i64) -> i32 {
        if pos >= 0 && (pos as usize) < self.pos2start.len() {
            self.pos2start[pos as usize]
        } else {
            -1
        }
    }

    #[inline]
    fn get_end(&self, pos: i64) -> i32 {
        if pos >= 0 && (pos as usize) < self.pos2end.len() {
            self.pos2end[pos as usize]
        } else {
            -1
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct PendingPair {
    start: i64,
    end: i64,
}

#[derive(Default)]
struct PairTracker {
    pending: BTreeMap<i64, FastMap<Vec<u8>, PendingPair>>,
}

impl PairTracker {
    fn observe(
        &mut self,
        name: &[u8],
        start: i64,
        end: i64,
        same_reference: bool,
        mate_start: i64,
    ) -> Option<(i64, i64)> {
        self.evict_before(start);
        let mut empty = false;
        let pair = same_reference
            .then(|| {
                self.pending.get_mut(&start).and_then(|pairs| {
                    let pair = pairs.remove(name);
                    empty = pairs.is_empty();
                    pair
                })
            })
            .flatten();
        if empty {
            self.pending.remove(&start);
        }
        if let Some(pair) = pair {
            return Some((pair.start, pair.end));
        }
        if same_reference && mate_start >= start {
            self.pending
                .entry(mate_start)
                .or_default()
                .insert(name.to_vec(), PendingPair { start, end });
        }
        None
    }

    fn evict_before(&mut self, position: i64) {
        while self
            .pending
            .first_key_value()
            .is_some_and(|(&expiry, _)| expiry < position)
        {
            self.pending.pop_first();
        }
    }
}

#[inline]
fn bam_endpos(start: i64, cigar_ops: impl Iterator<Item = (u8, u32)>) -> Result<i64> {
    let mut pos = start;
    for (op, len) in cigar_ops {
        if matches!(op, 0 | 2 | 3 | 7 | 8) {
            pos = pos.checked_add(i64::from(len)).ok_or_else(|| {
                RsomicsError::InvalidInput("CIGAR reference span overflows i64".to_owned())
            })?;
        }
    }
    Ok(pos)
}

fn validate_record_span(start: i64, end: i64, reference_length: i64, name: &[u8]) -> Result<()> {
    if start < 0 && end > start {
        return Err(RsomicsError::InvalidInput(format!(
            "read {} has a reference-consuming CIGAR without an alignment position",
            String::from_utf8_lossy(name)
        )));
    }
    if end > reference_length {
        return Err(RsomicsError::InvalidInput(format!(
            "read {} overhangs the end of its reference",
            String::from_utf8_lossy(name)
        )));
    }
    Ok(())
}

fn require_coordinate_order(
    previous: &mut Option<(i32, i64)>,
    coordinate: (i32, i64),
) -> Result<()> {
    if previous.is_some_and(|previous| coordinate < previous) {
        return Err(RsomicsError::InvalidInput(
            "ampliconstats input is not coordinate ordered".to_owned(),
        ));
    }
    *previous = Some(coordinate);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn accumulate_record(
    flags: u16,
    start: i64,
    end: i64,
    tlen: i32,
    qname: &[u8],
    args: &Options,
    amps: &[Amplicon],
    lookup: &PosLookup,
    stats: &mut AmpStats,
    ref_len: i64,
    pair_tracker: &mut PairTracker,
    same_mate_reference: bool,
    mate_start: i64,
) -> Result<()> {
    const BAM_FUNMAP: u16 = 0x4;
    const BAM_FPAIRED: u16 = 0x1;
    const BAM_FREVERSE: u16 = 0x10;
    const BAM_FMUNMAP: u16 = 0x8;
    const BAM_FSUPPLEMENTARY: u16 = 0x800;
    const BAM_FSECONDARY: u16 = 0x100;

    stats.nseq += 1;

    if (u32::from(flags) & u32::from(args.flag_require)) != u32::from(args.flag_require)
        || (u32::from(flags) & u32::from(args.flag_filter)) != 0
    {
        stats.nfiltered += 1;
        return Ok(());
    }

    if end == start && (args.flag_filter & BAM_FUNMAP) != 0 {
        stats.nfiltered += 1;
        return Ok(());
    }

    let is_paired = (flags & BAM_FPAIRED) != 0;
    let is_reverse = (flags & BAM_FREVERSE) != 0;
    let is_secondary = (flags & BAM_FSECONDARY) != 0;
    let is_supplementary = (flags & BAM_FSUPPLEMENTARY) != 0;
    let mate_unmapped = (flags & BAM_FMUNMAP) != 0;
    let is_unmapped = (flags & BAM_FUNMAP) != 0;

    let mut mstart = start;
    let mut prev_start = 0i64;
    let mut prev_end = 0i64;
    if is_paired
        && !is_supplementary
        && !is_secondary
        && let Some((ps, pe)) = pair_tracker.observe(
            qname,
            start,
            end,
            same_mate_reference && !mate_unmapped,
            mate_start,
        )
    {
        prev_start = ps;
        prev_end = pe;
        mstart = mstart.max(prev_end);
    }

    let depth_end = end.min(ref_len);
    for i in mstart..depth_end {
        stats.depth_all[i as usize] += 1;
    }

    let anum = if is_reverse || !is_paired {
        lookup.get_end(end - 1)
    } else {
        lookup.get_start(start)
    };

    if anum == -1 {
        stats.nfailprimer += 1;
    }

    if anum >= 0 {
        let a = anum as usize;
        let amp = &amps[a];
        let c_start = start.max(amp.max_left);
        let c_end = end.min(amp.min_right + 1);
        if c_end > c_start {
            stats.nreads[a] += 1;
            stats.nbases[a] += c_end - c_start;

            let ostart = start.max(amp.min_left - 1);
            let oend = end.min(amp.max_right);
            let offset = amp.min_left - 1;
            for i in ostart..oend {
                let apos = a * stats.max_amp_len + (i - offset) as usize;
                if apos < stats.coverage.len() {
                    stats.coverage[apos] += 1;
                }
            }
        } else {
            stats.nfailprimer += 1;
        }
    }

    let mut oth_anum: i32 = -1;
    let t_end: i64;

    if is_paired {
        let raw_t_end = if is_reverse { end } else { start } + tlen as i64;
        t_end = raw_t_end
            + if tlen > 0 {
                -args.tlen_adj
            } else {
                args.tlen_adj
            };
        if t_end > 0 && t_end < ref_len && tlen != 0 {
            oth_anum = if is_reverse {
                lookup.get_start(t_end)
            } else {
                lookup.get_end(t_end)
            };
        }
    } else {
        oth_anum = lookup.get_start(start);
        t_end = end;
    }

    let astatus: usize;
    if anum != -1 && oth_anum != -1 {
        astatus = if oth_anum == anum { 0 } else { 1 };
        if start <= t_end {
            stats.amp_dist[anum as usize][astatus] += 1;
        }
    } else if anum >= 0 {
        astatus = 2;
        stats.amp_dist[anum as usize][2] += 1;
    } else {
        astatus = 2;
    }

    if astatus == 0 && !is_unmapped && !mate_unmapped {
        if prev_end != 0 && mstart > prev_end {
            for i in prev_start..prev_end {
                if (i as usize) < stats.depth_valid.len() {
                    stats.depth_valid[i as usize] -= 1;
                }
            }
            if anum >= 0 {
                stats.nfull_reads[anum as usize] -= if is_paired { 0.5 } else { 1.0 };
            }
        } else {
            for i in mstart..depth_end {
                stats.depth_valid[i as usize] += 1;
            }
            if anum >= 0 {
                stats.nfull_reads[anum as usize] += if is_paired { 0.5 } else { 1.0 };
            }
        }
    }

    if is_paired && tlen <= 0 {
        return Ok(());
    }

    let tc_start = start;
    let tc_end = if is_paired {
        start + tlen as i64 - 1
    } else {
        end
    };
    let reported_start = tc_start.saturating_add(1).clamp(0, i64::from(u32::MAX)) as u64;
    let reported_end = tc_end.saturating_add(1).clamp(0, i64::from(u32::MAX)) as u64;
    let tcoord_key = reported_start | (reported_end << 32);

    let bucket = if anum >= 0 { anum as usize + 1 } else { 0 };
    if bucket < stats.tcoord.len() {
        let entry = stats.tcoord[bucket].entry(tcoord_key).or_insert(0);
        let count = (*entry & 0xFFFF_FFFF).wrapping_add(1);
        *entry = count | ((astatus as u64) << 32);
    }
    Ok(())
}

fn compute_covered_perc(stats: &mut AmpStats, amps: &[Amplicon], args: &Options) {
    for (a, amp) in amps.iter().enumerate() {
        let alen = (amp.min_right - amp.max_left + 1).max(1) as f64;
        let offset = amp.min_left - 1;
        for d in 0..MAX_DEPTH_LEVELS {
            if d > 0 && args.min_depth[d] == 0 {
                break;
            }
            let mut covered = 0i64;
            for j in (amp.max_left - 1)..amp.min_right {
                let apos = a * stats.max_amp_len + (j - offset) as usize;
                if apos < stats.coverage.len() && stats.coverage[apos] >= args.min_depth[d] as i64 {
                    covered += 1;
                }
            }
            stats.covered_perc[a][d] = 100.0 * covered as f64 / alen;
        }
    }
}

fn append_to_global(local: &AmpStats, global: &mut AmpStats, namp: usize, all_nseq: i64) {
    global.nseq += local.nseq;
    global.nfiltered += local.nfiltered;
    global.nfailprimer += local.nfailprimer;

    for a in 0..=namp {
        if a >= local.tcoord.len() || a >= global.tcoord.len() {
            break;
        }
        for (&key, &lval) in &local.tcoord[a] {
            let lcount = lval & 0xFFFF_FFFF;
            if lcount == 0 {
                continue;
            }
            let gentry = global.tcoord[a].entry(key).or_insert(0);
            *gentry =
                ((*gentry & 0xFFFF_FFFF).wrapping_add(lcount)) | (lval & 0xFFFF_FFFF_0000_0000);
        }
    }

    for a in 0..namp {
        global.nreads[a] += local.nreads[a];
        global.nreads2[a] += local.nreads[a] * local.nreads[a];
        global.nfull_reads[a] += local.nfull_reads[a];

        let rperc = if all_nseq > 0 {
            100.0 * local.nreads[a] as f64 / all_nseq as f64
        } else {
            0.0
        };
        global.nrperc[a] += rperc;
        global.nrperc2[a] += rperc * rperc;

        global.nbases[a] += local.nbases[a];
        global.nbases2[a] += local.nbases[a] * local.nbases[a];

        for d in 0..MAX_DEPTH_LEVELS {
            global.covered_perc[a][d] += local.covered_perc[a][d];
            global.covered_perc2[a][d] += local.covered_perc[a][d] * local.covered_perc[a][d];
        }

        for d in 0..3 {
            global.amp_dist[a][d] += local.amp_dist[a][d];
        }
    }

    for i in 0..local.ref_len as usize {
        if i < global.depth_all.len() {
            global.depth_all[i] += local.depth_all[i];
        }
        if i < global.depth_valid.len() {
            global.depth_valid[i] += local.depth_valid[i];
        }
    }
}

fn sample_name_from_path(path: &Path) -> String {
    path.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string()
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub files: usize,
    pub references: usize,
    pub amplicons: usize,
    pub records: u64,
}

pub fn write(
    args: &Options,
    bed_path: &Path,
    bam_paths: &[&Path],
    argv_str: &str,
    out: &mut impl Write,
) -> Result<Summary> {
    if bam_paths.is_empty() {
        return Err(RsomicsError::ConfigError(
            "ampliconstats requires at least one BAM input".to_owned(),
        ));
    }
    validate_options(args)?;
    let primer_bed = PrimerBed::read(bed_path)?;
    if primer_bed
        .references()
        .iter()
        .flat_map(|reference| &reference.primers)
        .any(|primer| primer.strand.is_none())
    {
        return Err(RsomicsError::InvalidInput(
            "ampliconstats requires strand in column six of every primer row".to_owned(),
        ));
    }
    if args.single_ref && primer_bed.references().len() != 1 {
        return Err(RsomicsError::ConfigError(
            "--single-ref requires a BED containing exactly one reference".to_owned(),
        ));
    }

    let mut first_reader = input::open(bam_paths[0], None, args.additional_threads)?;
    require_bam(first_reader.format())?;
    let header = first_reader.read_header(bam_paths[0])?;

    let nref_bam = header.reference_sequences().len();
    let multi_ref = !args.single_ref;

    writeln!(out, "# Summary statistics, used for scaling the plots.")?;
    writeln!(
        out,
        "SS\tSamtools version: rsomics-bam {}",
        env!("CARGO_PKG_VERSION")
    )?;
    writeln!(out, "SS\tCommand line: {argv_str}")?;
    writeln!(out, "SS\tNumber of files:\t{}", bam_paths.len())?;

    let mut refs: Vec<Option<RefData>> = Vec::with_capacity(nref_bam);
    let mut amp_offset = 0usize;
    let allocation_amplicons = args
        .max_amp
        .checked_add(1)
        .ok_or_else(|| RsomicsError::ConfigError("--max-amplicons is too large".to_owned()))?;
    let allocation_length = args.max_amp_len.checked_add(1).ok_or_else(|| {
        RsomicsError::ConfigError("--max-amplicon-length is too large".to_owned())
    })?;

    for (name, seq) in header.reference_sequences().iter() {
        let ref_name = name.to_string();
        let ref_len = usize::from(seq.length()) as i64;

        if let Some(primer_list) = primer_bed.get(&ref_name) {
            let namp = model::count(&primer_list.primers)?;
            if namp > args.max_amp {
                return Err(RsomicsError::ConfigError(format!(
                    "reference {ref_name} has {namp} amplicons, exceeding --max-amplicons {}",
                    args.max_amp
                )));
            }

            if multi_ref {
                writeln!(out, "SS\tNumber of amplicons:\t{}\t{}", ref_name, namp)?;
                writeln!(out, "SS\tReference length:\t{}\t{}", ref_name, ref_len)?;
            } else {
                writeln!(out, "SS\tNumber of amplicons:\t{}", namp)?;
                writeln!(out, "SS\tReference length:\t{}", ref_len)?;
            }

            let amps_placeholder: Vec<Amplicon> = Vec::new();
            let lookup = PosLookup {
                pos2start: vec![],
                pos2end: vec![],
            };
            let local = AmpStats::new(ref_len, allocation_amplicons, allocation_length)?;
            let global = AmpStats::new(ref_len, allocation_amplicons, allocation_length)?;

            refs.push(Some(RefData {
                ref_name,
                ref_len,
                amps: amps_placeholder,
                lookup,
                local,
                global,
                first_amp_idx: amp_offset,
            }));
            amp_offset += namp;
        } else {
            refs.push(None);
        }
    }
    for reference in primer_bed.references() {
        if !header
            .reference_sequences()
            .contains_key(reference.name.as_bytes())
        {
            return Err(RsomicsError::InvalidInput(format!(
                "primer reference is absent from the alignment dictionary: {}",
                reference.name
            )));
        }
    }

    writeln!(out, "SS\tEnd of summary")?;

    let mut first_with_primers = true;
    for slot in refs.iter_mut().flatten() {
        let ref_arg = if multi_ref {
            Some(slot.ref_name.as_str())
        } else {
            None
        };
        let primer_list = primer_bed.get(&slot.ref_name).ok_or_else(|| {
            RsomicsError::InvalidInput(format!("missing primers for {}", slot.ref_name))
        })?;

        let do_title = first_with_primers;
        first_with_primers = false;

        let amps = model::build(
            &primer_list.primers,
            out,
            ref_arg,
            slot.first_amp_idx,
            do_title,
            args.max_amp_len,
        )?;

        let lookup = PosLookup::build(&amps, slot.ref_len, args.max_delta);
        slot.amps = amps;
        slot.lookup = lookup;
    }

    let dictionary = reference_dictionary(&header);
    let mut record_count = 0u64;
    for bam_path in bam_paths {
        for slot in refs.iter_mut().flatten() {
            slot.local.reset();
        }

        let mut reader = input::open(bam_path, None, args.additional_threads)?;
        require_bam(reader.format())?;
        let input_header = reader.read_header(bam_path)?;
        if reference_dictionary(&input_header) != dictionary {
            return Err(RsomicsError::InvalidInput(format!(
                "alignment dictionary differs from the first input: {}",
                bam_path.display()
            )));
        }

        let sample_name = if args.use_sample_name {
            first_sample_name(&input_header).ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "--use-sample-name requires an SM field in the first read group: {}",
                    bam_path.display()
                ))
            })?
        } else {
            sample_name_from_path(bam_path)
        };

        let mut pair_trackers: Vec<PairTracker> =
            (0..nref_bam).map(|_| PairTracker::default()).collect();
        let mut previous_coordinate = None;
        let mut cigar = Vec::new();
        reader.visit_raw_bam_records(bam_path, |rec| {
            record_count = record_count.checked_add(1).ok_or_else(|| {
                RsomicsError::InvalidInput("alignment record count exceeds u64".to_owned())
            })?;
            validate_coordinate_fields(
                rec.reference_sequence_id(),
                rec.alignment_start(),
                rec.mate_reference_sequence_id(),
                rec.mate_alignment_start(),
                nref_bam,
            )?;
            let ref_id = rec.reference_sequence_id();
            if ref_id < 0 {
                return Ok(true);
            }
            let ref_idx = ref_id as usize;
            let start = i64::from(rec.alignment_start());
            let coordinate = (ref_id, start);
            require_coordinate_order(&mut previous_coordinate, coordinate)?;
            if refs[ref_idx].is_none() {
                return Ok(true);
            }

            let flags = rec.flags();
            let tlen = rec.template_length();
            let qname = rec.name();
            rec.decode_cigar_into(&mut cigar)?;
            let end = bam_endpos(start, cigar.iter().copied())?;
            let reference_length = refs[ref_idx].as_ref().unwrap().ref_len;
            validate_record_span(start, end, reference_length, qname)?;
            let mate_reference = rec.mate_reference_sequence_id();
            let mate_start = i64::from(rec.mate_alignment_start());

            let slot = refs[ref_idx].as_mut().ok_or_else(|| {
                RsomicsError::InvalidInput("reference statistics slot is absent".to_owned())
            })?;
            accumulate_record(
                flags,
                start,
                end,
                tlen,
                qname,
                args,
                &slot.amps,
                &slot.lookup,
                &mut slot.local,
                slot.ref_len,
                &mut pair_trackers[ref_idx],
                mate_reference == ref_id,
                mate_start,
            )?;
            Ok(true)
        })?;

        for slot in refs.iter_mut().flatten() {
            compute_covered_perc(&mut slot.local, &slot.amps, args);
        }

        let all_nseq: i64 = refs
            .iter()
            .flatten()
            .map(|s| s.local.nseq - s.local.nfiltered - s.local.nfailprimer)
            .sum();

        output::dump_stats(
            'F',
            &sample_name,
            bam_paths.len(),
            &refs,
            args,
            true,
            multi_ref,
            out,
        )?;

        for slot in refs.iter_mut().flatten() {
            let namp = slot.amps.len();
            append_to_global(&slot.local, &mut slot.global, namp, all_nseq);
        }
    }

    output::dump_stats(
        'C',
        "COMBINED",
        bam_paths.len(),
        &refs,
        args,
        false,
        multi_ref,
        out,
    )?;

    out.flush().map_err(RsomicsError::Io)?;
    Ok(Summary {
        files: bam_paths.len(),
        references: refs.iter().flatten().count(),
        amplicons: refs.iter().flatten().map(|slot| slot.amps.len()).sum(),
        records: record_count,
    })
}

fn validate_options(args: &Options) -> Result<()> {
    if args.max_delta < 0 {
        return Err(RsomicsError::ConfigError(
            "--pos-margin must not be negative".to_owned(),
        ));
    }
    if args.max_amp == 0 || args.max_amp_len == 0 {
        return Err(RsomicsError::ConfigError(
            "amplicon limits must be greater than zero".to_owned(),
        ));
    }
    if !(0.0..=1.0).contains(&args.depth_bin) {
        return Err(RsomicsError::ConfigError(
            "--depth-bin must be between zero and one".to_owned(),
        ));
    }
    if args.tcoord_bin <= 0 {
        return Err(RsomicsError::ConfigError(
            "--tcoord-bin must be greater than zero".to_owned(),
        ));
    }
    Ok(())
}

fn require_bam(format: input::Format) -> Result<()> {
    if format == input::Format::Bam {
        Ok(())
    } else {
        Err(RsomicsError::ConfigError(
            "ampliconstats 0.19 requires BAM input".to_owned(),
        ))
    }
}

fn reference_dictionary(header: &noodles::sam::Header) -> Vec<(String, usize)> {
    header
        .reference_sequences()
        .iter()
        .map(|(name, sequence)| (name.to_string(), usize::from(sequence.length())))
        .collect()
}

fn first_sample_name(header: &noodles::sam::Header) -> Option<String> {
    use noodles::sam::header::record::value::map::read_group::tag;

    header.read_groups().values().next().and_then(|read_group| {
        read_group
            .other_fields()
            .get(&tag::SAMPLE)
            .map(|value| String::from_utf8_lossy(AsRef::<[u8]>::as_ref(value)).into_owned())
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_and_cross_reference_mates_do_not_accumulate() {
        let mut tracker = PairTracker::default();
        assert_eq!(tracker.observe(b"missing", 10, 20, true, 50), None);
        assert_eq!(tracker.pending.len(), 1);

        assert_eq!(tracker.observe(b"cross", 51, 60, false, 80), None);
        assert!(tracker.pending.is_empty());
    }

    #[test]
    fn matching_mates_return_the_first_alignment() {
        let mut tracker = PairTracker::default();
        tracker.observe(b"pair", 10, 30, true, 20);
        assert_eq!(tracker.observe(b"pair", 20, 40, true, 10), Some((10, 30)));
        assert!(tracker.pending.is_empty());
    }

    #[test]
    fn mates_sharing_an_expected_coordinate_are_tracked_independently() {
        let mut tracker = PairTracker::default();
        tracker.observe(b"first", 10, 30, true, 20);
        tracker.observe(b"second", 11, 31, true, 20);

        assert_eq!(tracker.observe(b"first", 20, 40, true, 10), Some((10, 30)));
        assert_eq!(tracker.observe(b"second", 20, 41, true, 11), Some((11, 31)));
        assert!(tracker.pending.is_empty());
    }

    #[test]
    fn inconsistent_mate_metadata_does_not_form_a_pair() {
        let mut tracker = PairTracker::default();
        tracker.observe(b"pair", 10, 30, true, 20);

        assert_eq!(tracker.observe(b"pair", 20, 40, false, 10), None);
    }

    #[test]
    fn end_position_overflow_is_an_error() {
        let error = bam_endpos(i64::MAX, [(0, 1)].into_iter()).unwrap_err();
        assert!(error.to_string().contains("CIGAR reference span overflows"));
    }

    #[test]
    fn invalid_coordinates_and_order_are_errors() {
        use crate::alignment_order::validate_coordinate_fields;

        assert!(validate_coordinate_fields(-2, -1, -1, -1, 1).is_err());
        assert!(validate_coordinate_fields(1, 0, -1, -1, 1).is_err());
        assert!(validate_coordinate_fields(0, -2, 0, 0, 1).is_err());
        assert!(validate_record_span(-1, 9, 100, b"bad").is_err());
        assert!(validate_record_span(90, 101, 100, b"overhang").is_err());

        let mut previous = None;
        require_coordinate_order(&mut previous, (0, 20)).unwrap();
        let error = require_coordinate_order(&mut previous, (0, 19)).unwrap_err();
        assert!(error.to_string().contains("not coordinate ordered"));
    }
}
