#![allow(unsafe_code)]

use std::ffi::c_void;
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rust_htslib::bam::{self, Read as _};
use rust_htslib::htslib;

use crate::coverage;
use crate::input;

pub(crate) struct Reference {
    pub(crate) name: Vec<u8>,
    pub(crate) length: u64,
}

#[derive(Default)]
pub(crate) struct Stats {
    pub(crate) reads: u64,
    pub(crate) covered_bases: u64,
    pub(crate) depth_sum: u64,
    pub(crate) base_quality_sum: u64,
    pub(crate) quality_bases: u64,
    pub(crate) mapping_quality_sum: u64,
}

pub(crate) struct Scan {
    pub(crate) references: Vec<Reference>,
    pub(crate) stats: Vec<Stats>,
}

struct ReadState {
    file: *mut htslib::samFile,
    header: *mut htslib::sam_hdr_t,
    stats: *mut Stats,
    stats_len: usize,
    minimum_mapping_quality: u8,
    minimum_read_length: usize,
    required_flags: u16,
    excluded_flags: u16,
}

struct Pileup(htslib::bam_plp_t);

impl Drop for Pileup {
    fn drop(&mut self) {
        unsafe { htslib::bam_plp_destroy(self.0) };
    }
}

pub(crate) fn collect(input_path: &Path, options: coverage::Options<'_>) -> Result<Scan> {
    let format = input::detect_format(input_path)?;
    let mut reader = bam::Reader::from_path(input_path)
        .map_err(|error| hts_error("opening alignment", input_path, error))?;
    if let Some(reference) = options.reference {
        reader
            .set_reference(reference)
            .map_err(|error| hts_error("attaching CRAM reference", input_path, error))?;
    }
    if options.additional_threads > 0 {
        reader
            .set_threads(options.additional_threads)
            .map_err(|error| hts_error("configuring alignment threads", input_path, error))?;
    }
    if format == input::Format::Cram {
        configure_cram(&mut reader, options, input_path)?;
    }

    let references = reference_dictionary(reader.header())?;
    let mut stats: Vec<_> = (0..references.len()).map(|_| Stats::default()).collect();
    let mut state = ReadState {
        file: reader.htsfile(),
        header: reader.header().inner_ptr() as *mut htslib::sam_hdr_t,
        stats: stats.as_mut_ptr(),
        stats_len: stats.len(),
        minimum_mapping_quality: options.minimum_mapping_quality,
        minimum_read_length: options.minimum_read_length,
        required_flags: options.required_flags,
        excluded_flags: options.excluded_flags,
    };
    // The reader and fixed-size stats buffer outlive every synchronous pileup callback.
    let raw = unsafe {
        htslib::bam_plp_init(
            Some(read_filtered),
            (&mut state as *mut ReadState).cast::<c_void>(),
        )
    };
    if raw.is_null() {
        return Err(RsomicsError::InvalidInput(
            "initializing coverage pileup failed".to_owned(),
        ));
    }
    let pileup = Pileup(raw);
    let maximum_depth = if options.maximum_depth == 0 {
        i32::MAX
    } else {
        i32::try_from(options.maximum_depth).map_err(|_| {
            RsomicsError::ConfigError(format!("maximum coverage depth exceeds {}", i32::MAX))
        })?
    };
    unsafe { htslib::bam_plp_set_maxcnt(pileup.0, maximum_depth) };

    loop {
        let mut reference_id = 0;
        let mut _position = 0;
        let mut count = 0;
        let alignments = unsafe {
            htslib::bam_plp_auto(pileup.0, &mut reference_id, &mut _position, &mut count)
        };
        if alignments.is_null() {
            if count == -1 {
                return Err(RsomicsError::InvalidInput(format!(
                    "building coverage pileup from {} failed",
                    input_path.display()
                )));
            }
            break;
        }
        let reference_id = usize::try_from(reference_id).map_err(|_| {
            RsomicsError::InvalidInput("coverage pileup has a negative reference ID".to_owned())
        })?;
        let stats = stats.get_mut(reference_id).ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "alignment reference ID {reference_id} is absent from the header"
            ))
        })?;
        let alignments =
            unsafe { std::slice::from_raw_parts(alignments, usize::try_from(count).unwrap()) };
        let mut depth = u64::try_from(count).unwrap();
        let mut base_quality_sum = 0u64;
        let mut quality_bases = 0u64;
        for alignment in alignments {
            if alignment.is_del() != 0 || alignment.is_refskip() != 0 {
                depth -= 1;
                continue;
            }
            let record = unsafe { &*alignment.b };
            if alignment.qpos < 0 || alignment.qpos >= record.core.l_qseq {
                continue;
            }
            let quality = unsafe { quality_at(record, usize::try_from(alignment.qpos).unwrap()) };
            if quality < options.minimum_base_quality {
                depth -= 1;
            } else {
                base_quality_sum += u64::from(quality);
                quality_bases += 1;
            }
        }
        if depth >= u64::try_from(options.minimum_depth).unwrap() {
            stats.covered_bases += 1;
            stats.depth_sum += depth;
            stats.base_quality_sum += base_quality_sum;
            stats.quality_bases += quality_bases;
        }
    }

    Ok(Scan { references, stats })
}

unsafe extern "C" fn read_filtered(data: *mut c_void, record: *mut htslib::bam1_t) -> i32 {
    let state = unsafe { &mut *data.cast::<ReadState>() };
    loop {
        let result = unsafe { htslib::sam_read1(state.file, state.header, record) };
        if result < 0 {
            return result;
        }
        let record = unsafe { &*record };
        let core = &record.core;
        if core.tid < 0
            || core.flag & state.excluded_flags != 0
            || (state.required_flags != 0 && core.flag & state.required_flags == 0)
            || core.qual < state.minimum_mapping_quality
            || (state.minimum_read_length != 0
                && unsafe { query_length(record) } < state.minimum_read_length)
        {
            continue;
        }
        let reference_id = match usize::try_from(core.tid) {
            Ok(reference_id) if reference_id < state.stats_len => reference_id,
            _ => return -2,
        };
        let stats = unsafe { &mut *state.stats.add(reference_id) };
        stats.reads += 1;
        stats.mapping_quality_sum += u64::from(core.qual);
        return result;
    }
}

unsafe fn query_length(record: &htslib::bam1_t) -> usize {
    let cigar = unsafe {
        record
            .data
            .add(usize::from(record.core.l_qname))
            .cast::<u32>()
    };
    let mut length = 0usize;
    for index in 0..record.core.n_cigar as usize {
        let operation = unsafe { *cigar.add(index) };
        if matches!(operation & 0x0f, 0 | 1 | 4 | 7 | 8) {
            length = length.saturating_add((operation >> 4) as usize);
        }
    }
    length
}

unsafe fn quality_at(record: &htslib::bam1_t, query_position: usize) -> u8 {
    let core = &record.core;
    let cigar_bytes = core.n_cigar as usize * std::mem::size_of::<u32>();
    let sequence_bytes = usize::try_from(core.l_qseq).unwrap().div_ceil(2);
    // BAM stores CIGAR, packed sequence, and qualities consecutively after QNAME.
    let offset = usize::from(core.l_qname) + cigar_bytes + sequence_bytes + query_position;
    unsafe { *record.data.add(offset) }
}

fn configure_cram(
    reader: &mut bam::Reader,
    options: coverage::Options<'_>,
    input_path: &Path,
) -> Result<()> {
    let mut fields = htslib::sam_fields_SAM_FLAG
        | htslib::sam_fields_SAM_RNAME
        | htslib::sam_fields_SAM_POS
        | htslib::sam_fields_SAM_MAPQ
        | htslib::sam_fields_SAM_CIGAR
        | htslib::sam_fields_SAM_SEQ;
    if options.minimum_base_quality > 0 {
        fields |= htslib::sam_fields_SAM_QUAL;
    }
    reader
        .set_cram_options(htslib::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS, fields)
        .map_err(|error| hts_error("configuring CRAM fields", input_path, error))?;
    reader
        .set_cram_options(htslib::hts_fmt_option_CRAM_OPT_DECODE_MD, 0)
        .map_err(|error| hts_error("disabling CRAM MD decoding", input_path, error))?;
    Ok(())
}

fn reference_dictionary(header: &bam::HeaderView) -> Result<Vec<Reference>> {
    header
        .target_names()
        .iter()
        .enumerate()
        .map(|(reference_id, name)| {
            let length = header
                .target_len(u32::try_from(reference_id).unwrap())
                .ok_or_else(|| {
                    RsomicsError::InvalidInput(format!(
                        "reference {reference_id} has no declared length"
                    ))
                })?;
            Ok(Reference {
                name: name.to_vec(),
                length,
            })
        })
        .collect()
}

fn hts_error(action: &str, input_path: &Path, error: rust_htslib::errors::Error) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{action} from {}: {error}", input_path.display()))
}
