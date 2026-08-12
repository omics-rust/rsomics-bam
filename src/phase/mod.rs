mod errmod;
mod model;
mod output;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use noodles::sam;
use rsomics_bamio::raw::{RawRecord, RawRecordEncoder};
use rsomics_common::{Result, RsomicsError};
use rsomics_pileup::{FlagFilter, PileupEngine, PileupOptions, RecordFilter};
use serde::Serialize;

use self::errmod::Errmod;
use self::model::{Fragment, Site};
use self::output::Router;
pub(crate) use self::output::{Format as PartitionFormat, paths as partition_paths};
use crate::Program;

const EXCLUDED_FLAGS: u16 = 0x4 | 0x100 | 0x200 | 0x400;
pub(crate) const MAX_DEPTH: usize = u16::MAX as usize;
pub(crate) const MAX_WINDOW: usize = 23;

#[derive(Clone)]
pub struct Options {
    pub window: usize,
    pub minimum_lod: i32,
    pub minimum_base_quality: u8,
    pub maximum_depth: usize,
    pub fix_chimeras: bool,
    pub reference: Option<PathBuf>,
    pub additional_threads: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            window: 13,
            minimum_lod: 37,
            minimum_base_quality: 13,
            maximum_depth: 256,
            fix_chimeras: true,
            reference: None,
            additional_threads: 0,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct Summary {
    pub phase_sets: u64,
    pub heterozygous_sites: u64,
}

pub fn write(input_path: &Path, options: Options, mut output: impl Write) -> Result<Summary> {
    run(input_path, options, &mut output, None)
}

pub(crate) fn write_partitioned(
    input_path: &Path,
    options: Options,
    mut output: impl Write,
    files: [File; 3],
    format: PartitionFormat,
    drop_ambiguous: bool,
    program: Option<Program<'_>>,
) -> Result<Summary> {
    run(
        input_path,
        options,
        &mut output,
        Some(Partition {
            files,
            format,
            drop_ambiguous,
            program,
        }),
    )
}

struct Partition<'a> {
    files: [File; 3],
    format: PartitionFormat,
    drop_ambiguous: bool,
    program: Option<Program<'a>>,
}

fn run(
    input_path: &Path,
    options: Options,
    output: &mut impl Write,
    partition: Option<Partition<'_>>,
) -> Result<Summary> {
    validate(&options)?;
    let mut reader = crate::input::open(
        input_path,
        options.reference.as_deref(),
        options.additional_threads,
    )?;
    let header = reader.read_header(input_path)?;
    let references: Vec<_> = header
        .reference_sequences()
        .iter()
        .map(|(name, reference)| (name.to_string(), usize::from(reference.length()) as u64))
        .collect();
    let pileup_options = PileupOptions {
        filter: RecordFilter {
            flags: FlagFilter {
                skip_any_set: EXCLUDED_FLAGS,
                ..FlagFilter::default()
            },
            ..RecordFilter::default()
        },
        adjust_overlaps: false,
        maximum_depth_per_source: None,
    };
    let mut pileup =
        PileupEngine::new(references.iter().map(|(_, length)| *length), pileup_options);
    let router = partition
        .map(|partition| {
            let repository = options
                .reference
                .as_deref()
                .map(crate::input::reference_repository)
                .transpose()?
                .unwrap_or_default();
            let mut output_header = header.clone();
            if let Some(program) = partition.program {
                program.add_to(&mut output_header)?;
            }
            Router::new(
                partition.files,
                partition.format,
                repository,
                output_header,
                partition.drop_ambiguous,
            )
        })
        .transpose()?;
    let model = Errmod::new()?;
    let mut output = BufWriter::with_capacity(256 * 1024, output);
    write_header(&mut output)?;
    let mut state = State::new(options, &references, &model, &mut output, router);

    if state.router.is_some() {
        let mut encoder = RawRecordEncoder::new();
        reader.visit_records(&header, input_path, |source| {
            let record = encoder.encode(&header, source)?;
            if record.flags() & EXCLUDED_FLAGS == 0 {
                let record_buf =
                    sam::alignment::RecordBuf::try_from_alignment_record(&header, source)
                        .map_err(RsomicsError::Io)?;
                state.queue(&record, record_buf)?;
            }
            pileup.push(record).map_err(pileup_error)?;
            pileup.drain(|column| state.column(column))?;
            Ok(true)
        })?;
    } else {
        reader.visit_owned_raw_records(&header, input_path, |record| {
            pileup.push(record).map_err(pileup_error)?;
            pileup.drain(|column| state.column(column))?;
            Ok(true)
        })?;
    }
    pileup.finish().map_err(pileup_error)?;
    pileup.drain(|column| state.column(column))?;
    state.finish()?;
    let summary = state.summary;
    drop(state);
    output.flush().map_err(RsomicsError::Io)?;
    Ok(summary)
}

fn validate(options: &Options) -> Result<()> {
    if options.window == 0 || options.window > MAX_WINDOW {
        return Err(RsomicsError::ConfigError(format!(
            "phase window must be between 1 and {MAX_WINDOW}"
        )));
    }
    if options.maximum_depth == 0 || options.maximum_depth > MAX_DEPTH {
        return Err(RsomicsError::ConfigError(format!(
            "maximum phase depth must be between 1 and {}",
            MAX_DEPTH
        )));
    }
    Ok(())
}

struct State<'a, W> {
    options: Options,
    references: &'a [(String, u64)],
    model: &'a Errmod,
    output: &'a mut W,
    reference_id: Option<i32>,
    bases: Vec<u16>,
    sites: Vec<Site>,
    fragments: Vec<Fragment>,
    fragment_indices: HashMap<Arc<[u8]>, usize>,
    marker_offset: usize,
    summary: Summary,
    router: Option<Router>,
}

impl<'a, W: Write> State<'a, W> {
    fn new(
        options: Options,
        references: &'a [(String, u64)],
        model: &'a Errmod,
        output: &'a mut W,
        router: Option<Router>,
    ) -> Self {
        Self {
            options,
            references,
            model,
            output,
            reference_id: None,
            bases: Vec::new(),
            sites: Vec::new(),
            fragments: Vec::new(),
            fragment_indices: HashMap::new(),
            marker_offset: 0,
            summary: Summary::default(),
            router,
        }
    }

    fn queue(&mut self, record: &RawRecord, record_buf: sam::alignment::RecordBuf) -> Result<()> {
        let end = alignment_end(record)?;
        if let Some(router) = self.router.as_mut() {
            router.push(
                record.reference_sequence_id(),
                end,
                record.name().into(),
                record_buf,
            );
        }
        Ok(())
    }

    fn column(&mut self, column: &rsomics_pileup::Column<'_>) -> Result<()> {
        if self.reference_id != Some(column.reference_id()) {
            self.flush_block(i64::MAX)?;
            self.reference_id = Some(column.reference_id());
            self.marker_offset = 0;
        }
        if column.len() > self.options.maximum_depth {
            return Ok(());
        }
        self.bases.clear();
        self.bases.reserve(column.len());
        let mut observed = 0_u8;
        for entry in column.entries() {
            let projection = entry.projection();
            if projection.is_deletion || projection.is_reference_skip {
                continue;
            }
            let record = entry.record();
            let Some(&base_quality) = record.quality_scores().get(projection.qpos) else {
                continue;
            };
            if base_quality < self.options.minimum_base_quality {
                continue;
            }
            let Some(base) = base_index(record.seq_nibble(projection.qpos)) else {
                continue;
            };
            observed |= 1 << base;
            let quality = base_quality.min(record.mapping_quality()).clamp(4, 63);
            let reverse = u16::from(record.flags() & 0x10 != 0);
            self.bases
                .push(u16::from(quality) << 5 | reverse << 4 | u16::from(base));
        }
        if observed.count_ones() < 2 {
            return Ok(());
        }
        let Some(call) = self.model.call(&mut self.bases)? else {
            return Ok(());
        };
        if call.lod < self.options.minimum_lod {
            return Ok(());
        }
        let mut observations = Vec::with_capacity(column.len());
        let mut overlaps = false;
        for entry in column.entries() {
            let projection = entry.projection();
            if projection.is_deletion || projection.is_reference_skip {
                continue;
            }
            let record = entry.record();
            if record.mapping_quality() == 0 {
                continue;
            }
            let name = record.name();
            overlaps |= self.fragment_indices.contains_key(name);
            let base = base_index(record.seq_nibble(projection.qpos));
            let allele = match base {
                Some(base) if base == call.alleles[0] => 1,
                Some(base) if base == call.alleles[1] => 2,
                _ => 0,
            };
            observations.push((name, allele, record.alignment_start()));
        }
        if !overlaps && !self.sites.is_empty() {
            self.flush_block(column.position())?;
        } else if !self.sites.is_empty() {
            model::ensure_workspace(self.options.window, self.sites.len() + 1)?;
        }
        let site_index = self.sites.len();
        self.sites.push(Site {
            position: column.position(),
            alleles: call.alleles,
        });
        self.summary.heterozygous_sites += 1;
        for (name, allele, alignment_start) in observations {
            if let Some(&index) = self.fragment_indices.get(name) {
                self.fragments[index].push(site_index, allele);
            } else {
                let key: Arc<[u8]> = name.into();
                let index = self.fragments.len();
                self.fragments.push(Fragment::new(
                    key.clone(),
                    site_index,
                    allele,
                    alignment_start,
                ));
                self.fragment_indices.insert(key, index);
            }
        }
        Ok(())
    }

    fn flush_block(&mut self, cutoff: i64) -> Result<()> {
        let Some(reference_id) = self.reference_id else {
            self.clear();
            return Ok(());
        };
        if self.sites.is_empty() {
            if let Some(router) = self.router.as_mut() {
                router.route(reference_id, cutoff, &[])?;
            }
            self.clear();
            return Ok(());
        }
        let reference = self
            .references
            .get(reference_id as usize)
            .ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "pileup reference ID {reference_id} is absent from the header"
                ))
            })?
            .0
            .as_str();
        let count = model::write_block(
            self.output,
            reference,
            &self.sites,
            &mut self.fragments,
            self.options.window,
            self.options.fix_chimeras,
            self.marker_offset,
        )?;
        self.marker_offset += count;
        self.summary.phase_sets += 1;
        if let Some(router) = self.router.as_mut() {
            router.route(reference_id, cutoff, &self.fragments)?;
        }
        self.clear();
        Ok(())
    }

    fn clear(&mut self) {
        self.sites.clear();
        self.fragments.clear();
        self.fragment_indices.clear();
    }

    fn finish(&mut self) -> Result<()> {
        self.flush_block(i64::MAX)?;
        if let Some(router) = self.router.take() {
            router.finish()?;
        }
        Ok(())
    }
}

fn write_header(output: &mut impl Write) -> Result<()> {
    output
        .write_all(
            b"CC\n\
CC\tDescriptions:\n\
CC\n\
CC\t  CC      comments\n\
CC\t  PS      start of a phase set\n\
CC\t  FL      filtered region\n\
CC\t  M[012]  markers; 0 for singletons, 1 for phased and 2 for filtered\n\
CC\t  EV      supporting reads; SAM format\n\
CC\t  //      end of a phase set\n\
CC\n\
CC\tFormats of PS, FL and M[012] lines (1-based coordinates):\n\
CC\n\
CC\t  PS  chr  phaseSetStart  phaseSetEnd\n\
CC\t  FL  chr  filterStart    filterEnd\n\
CC\t  M?  chr  PS  pos  allele0  allele1  hetIndex  #supports0  #errors0  #supp1  #err1\n\
CC\n\
CC\n",
        )
        .map_err(RsomicsError::Io)
}

fn base_index(nibble: u8) -> Option<u8> {
    match nibble {
        1 => Some(0),
        2 => Some(1),
        4 => Some(2),
        8 => Some(3),
        _ => None,
    }
}

fn alignment_end(record: &RawRecord) -> Result<i64> {
    record
        .decoded_cigar()?
        .into_iter()
        .filter(|(kind, _)| matches!(kind, 0 | 2 | 3 | 7 | 8))
        .try_fold(i64::from(record.alignment_start()), |end, (_, length)| {
            end.checked_add(i64::from(length))
                .ok_or_else(|| RsomicsError::InvalidInput("alignment end exceeds i64".to_owned()))
        })
}

fn pileup_error(error: rsomics_pileup::PileupError) -> RsomicsError {
    RsomicsError::InvalidInput(format!("building phase pileup: {error}"))
}
