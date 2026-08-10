mod record;

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use noodles::sam;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::{self, header::tag as header_tag};
use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::amplicon::{PrimerBed, Strand};
use crate::output::{Compression, Format, Writer};
use crate::{Program, input};
use record::{Clipping, active_query_len, clip_left, clip_right, end_position, unmap};

const FLAG_UNMAPPED: u16 = 0x04;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_QC_FAIL: u16 = 0x200;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ClipMode {
    Soft,
    Hard,
}

#[derive(Clone, Debug)]
pub struct Options<'a> {
    pub mode: ClipMode,
    pub both_ends: bool,
    pub use_strand: bool,
    pub tolerance: i64,
    pub mark_fail: bool,
    pub clipped_only: bool,
    pub exclude_flagged: bool,
    pub filter_length: Option<i64>,
    pub fail_length: Option<i64>,
    pub unmap_length: Option<i64>,
    pub keep_tags: bool,
    pub original: bool,
    pub uncompressed: bool,
    pub additional_threads: usize,
    pub program: Option<Program<'a>>,
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub total: u64,
    pub clipped: u64,
    pub forward_clipped: u64,
    pub reverse_clipped: u64,
    pub both_clipped: u64,
    pub not_clipped: u64,
    pub excluded: u64,
    pub filtered: u64,
    pub failed: u64,
    pub written: u64,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct PrimerCount {
    pub reference: String,
    pub start: i64,
    pub end: i64,
    pub name: String,
    pub score: String,
    pub strand: Option<String>,
    pub count: u64,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Run {
    pub summary: Summary,
    pub primer_counts: Vec<PrimerCount>,
}

struct Site {
    start: i64,
    end: i64,
    strand: Option<Strand>,
    original_index: usize,
}

struct SiteIndex {
    sites: Vec<Site>,
    longest: i64,
}

pub fn write(
    input_path: &Path,
    bed_path: &Path,
    options: Options<'_>,
    output: Box<dyn Write + Send>,
    rejects: Option<Box<dyn Write + Send>>,
) -> Result<Run> {
    if options.tolerance < 0 {
        return Err(RsomicsError::ConfigError(
            "--tolerance must not be negative".to_owned(),
        ));
    }
    let bed = PrimerBed::read(bed_path)?;
    if options.use_strand
        && bed
            .references()
            .iter()
            .flat_map(|reference| &reference.primers)
            .any(|primer| primer.strand.is_none())
    {
        return Err(RsomicsError::InvalidInput(
            "--strand requires column six for every primer BED row".to_owned(),
        ));
    }

    let mut reader = input::open(input_path, None, options.additional_threads)?;
    if reader.format() != input::Format::Bam {
        return Err(RsomicsError::ConfigError(
            "ampliconclip 0.19 requires BAM input".to_owned(),
        ));
    }
    let mut header = reader.read_header(input_path)?;
    let reference_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|name| String::from_utf8_lossy(name.as_ref()).into_owned())
        .collect();
    for reference in bed.references() {
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
    set_unknown_order(&mut header);
    if let Some(program) = options.program {
        program.add_to(&mut header)?;
    }

    let compression = if options.uncompressed {
        Compression::Uncompressed
    } else {
        Compression::Default
    };
    let mut writer = Writer::new(Format::Bam, compression, options.additional_threads, output);
    writer.write_header(&header)?;
    let mut reject_writer = rejects
        .map(|output| {
            let mut writer =
                Writer::new(Format::Bam, compression, options.additional_threads, output);
            writer.write_header(&header).map(|()| writer)
        })
        .transpose()?;

    let (mut indices, reference_lookup, mut counts) = build_indices(&bed);
    let clipping = match options.mode {
        ClipMode::Soft => Clipping::Soft,
        ClipMode::Hard => Clipping::Hard,
    };
    let mut summary = Summary::default();
    let mut previous_coordinate = None;

    reader.visit_owned_raw_records(&header, input_path, |mut record| {
        summary.total += 1;
        let coordinate = (record.reference_sequence_id(), record.alignment_start());
        if coordinate.0 >= 0 {
            if previous_coordinate.is_some_and(|previous| coordinate < previous) {
                return Err(RsomicsError::InvalidInput(
                    "ampliconclip input is not coordinate ordered".to_owned(),
                ));
            }
            previous_coordinate = Some(coordinate);
        }

        let excluded = record.flags() & (FLAG_UNMAPPED | FLAG_QC_FAIL) != 0;
        let reference_id = record.reference_sequence_id();
        let site_index = usize::try_from(reference_id)
            .ok()
            .and_then(|index| reference_names.get(index))
            .and_then(|name| reference_lookup.get(name).copied());
        let mut filtered = false;
        let mut was_clipped = false;

        if !excluded {
            if let Some(site_index) = site_index {
                let sites = &mut indices[site_index];
                if options.both_ends {
                    let left = matching_site(
                        sites,
                        i64::from(record.alignment_start()),
                        false,
                        options.use_strand,
                        options.tolerance,
                    );
                    if let Some((bases, primer)) = left {
                        if options.original {
                            add_original_tag(&mut record, &reference_names)?;
                        }
                        record = clip_left(&record, bases, clipping)?;
                        counts[site_index][primer] += 1;
                        summary.forward_clipped += 1;
                        was_clipped = true;
                    }
                    let right = matching_site(
                        sites,
                        end_position(&record)?,
                        true,
                        options.use_strand,
                        options.tolerance,
                    );
                    if let Some((bases, primer)) = right {
                        if options.original && !was_clipped {
                            add_original_tag(&mut record, &reference_names)?;
                        }
                        record = clip_right(&record, bases, clipping)?;
                        counts[site_index][primer] += 1;
                        summary.reverse_clipped += 1;
                        if was_clipped {
                            summary.both_clipped += 1;
                        }
                        was_clipped = true;
                    }
                } else {
                    let reverse = record.flags() & FLAG_REVERSE != 0;
                    let position = if reverse {
                        end_position(&record)?
                    } else {
                        i64::from(record.alignment_start())
                    };
                    if let Some((bases, primer)) = matching_site(
                        sites,
                        position,
                        reverse,
                        options.use_strand,
                        options.tolerance,
                    ) {
                        if options.original {
                            add_original_tag(&mut record, &reference_names)?;
                        }
                        record = if reverse {
                            summary.reverse_clipped += 1;
                            clip_right(&record, bases, clipping)?
                        } else {
                            summary.forward_clipped += 1;
                            clip_left(&record, bases, clipping)?
                        };
                        counts[site_index][primer] += 1;
                        was_clipped = true;
                    }
                }
            }

            if !was_clipped {
                summary.not_clipped += 1;
                if options.mark_fail {
                    record.set_flag_bits(FLAG_QC_FAIL);
                }
            } else if !options.keep_tags {
                record.remove_aux(*b"NM");
                record.remove_aux(*b"MD");
            }

            let length = active_query_len(&record)?;
            if options.fail_length.is_some_and(|limit| length <= limit) {
                record.set_flag_bits(FLAG_QC_FAIL);
            }
            if options.filter_length.is_some_and(|limit| length <= limit) {
                filtered = true;
            }
            if options.unmap_length.is_some_and(|limit| length <= limit) {
                record = unmap(&record)?;
            }
            if record.flags() & FLAG_QC_FAIL != 0 {
                summary.failed += 1;
            }
            if options.clipped_only && !was_clipped {
                filtered = true;
            }
        } else {
            summary.excluded += 1;
            if options.exclude_flagged {
                filtered = true;
            }
        }

        if filtered {
            summary.filtered += 1;
            if let Some(writer) = &mut reject_writer {
                writer.write_owned_raw_record(&record)?;
            }
        } else {
            summary.written += 1;
            writer.write_owned_raw_record(&record)?;
        }
        Ok(true)
    })?;

    writer.finish(&header)?;
    if let Some(writer) = reject_writer {
        writer.finish(&header)?;
    }
    summary.clipped = summary.forward_clipped + summary.reverse_clipped;
    Ok(Run {
        summary,
        primer_counts: flatten_counts(&bed, &counts, options.use_strand),
    })
}

fn build_indices(bed: &PrimerBed) -> (Vec<SiteIndex>, HashMap<String, usize>, Vec<Vec<u64>>) {
    let mut indices = Vec::with_capacity(bed.references().len());
    let mut lookup = HashMap::new();
    let mut counts = Vec::with_capacity(bed.references().len());
    for (reference_index, reference) in bed.references().iter().enumerate() {
        lookup.insert(reference.name.clone(), reference_index);
        counts.push(vec![0; reference.primers.len()]);
        let longest = reference
            .primers
            .iter()
            .map(|primer| primer.end - primer.start)
            .max()
            .unwrap_or(0);
        let mut sites: Vec<_> = reference
            .primers
            .iter()
            .enumerate()
            .map(|(original_index, primer)| Site {
                start: primer.start,
                end: primer.end,
                strand: primer.strand,
                original_index,
            })
            .collect();
        sites.sort_by_key(|site| site.end);
        indices.push(SiteIndex { sites, longest });
    }
    (indices, lookup, counts)
}

fn matching_site(
    index: &SiteIndex,
    position: i64,
    reverse: bool,
    use_strand: bool,
    tolerance: i64,
) -> Option<(u32, usize)> {
    let threshold = if reverse {
        position.saturating_sub(tolerance).max(0)
    } else {
        position
    };
    let mut best = None;
    let first = index.sites.partition_point(|site| site.end <= threshold);
    for site in &index.sites[first.saturating_sub(1)..] {
        if use_strand
            && site.strand
                != Some(if reverse {
                    Strand::Reverse
                } else {
                    Strand::Forward
                })
        {
            continue;
        }
        let left = if reverse {
            site.start
        } else {
            site.start.saturating_sub(tolerance).max(0)
        };
        let right = if reverse {
            site.end.saturating_add(tolerance)
        } else {
            site.end
        };
        if position
            .saturating_add(index.longest)
            .saturating_add(tolerance)
            < right
        {
            break;
        }
        if (left..=right).contains(&position) {
            let size = if reverse {
                position - site.start
            } else {
                site.end - position
            };
            if size > 0 && best.is_none_or(|(current, _)| size > current) {
                best = Some((size, site.original_index));
            }
        }
    }
    best.and_then(|(size, index)| u32::try_from(size).ok().map(|size| (size, index)))
}

fn set_unknown_order(header: &mut sam::Header) {
    let fields = header
        .header_mut()
        .get_or_insert_with(Map::<map::Header>::default)
        .other_fields_mut();
    if fields
        .get(&header_tag::SORT_ORDER)
        .is_some_and(|value| AsRef::<[u8]>::as_ref(value) == b"coordinate")
    {
        fields.insert(header_tag::SORT_ORDER, "unknown".into());
        fields.shift_remove(&header_tag::SUBSORT_ORDER);
    }
}

fn add_original_tag(record: &mut RawRecord, _references: &[String]) -> Result<()> {
    if record.aux_value(*b"OA").is_some() {
        return Ok(());
    }
    let reference = String::from_utf8_lossy(record.name());
    let strand = if record.flags() & FLAG_REVERSE != 0 {
        '-'
    } else {
        '+'
    };
    let cigar = cigar_text(record)?;
    let nm = integer_aux(record, *b"NM").map_or(String::new(), |value| value.to_string());
    let value = format!(
        "{reference},{},{strand},{cigar},{},{nm};\0",
        i64::from(record.alignment_start()) + 1,
        record.mapping_quality()
    );
    record.set_aux(*b"OA", b'Z', value.as_bytes())
}

fn cigar_text(record: &RawRecord) -> Result<String> {
    let mut text = String::new();
    for (kind, length) in record.decoded_cigar()? {
        text.push_str(&length.to_string());
        text.push(match kind {
            0 => 'M',
            1 => 'I',
            2 => 'D',
            3 => 'N',
            4 => 'S',
            5 => 'H',
            6 => 'P',
            7 => '=',
            8 => 'X',
            _ => {
                return Err(RsomicsError::InvalidInput(
                    "invalid CIGAR operation".to_owned(),
                ));
            }
        });
    }
    if text.is_empty() {
        text.push('*');
    }
    Ok(text)
}

fn integer_aux(record: &RawRecord, tag: [u8; 2]) -> Option<i64> {
    let value = record.aux_value(tag)?;
    match record.aux_type(tag)? {
        b'c' => Some(i64::from(i8::from_le_bytes(value.try_into().ok()?))),
        b'C' => Some(i64::from(u8::from_le_bytes(value.try_into().ok()?))),
        b's' => Some(i64::from(i16::from_le_bytes(value.try_into().ok()?))),
        b'S' => Some(i64::from(u16::from_le_bytes(value.try_into().ok()?))),
        b'i' => Some(i64::from(i32::from_le_bytes(value.try_into().ok()?))),
        b'I' => Some(i64::from(u32::from_le_bytes(value.try_into().ok()?))),
        _ => None,
    }
}

fn flatten_counts(bed: &PrimerBed, counts: &[Vec<u64>], use_strand: bool) -> Vec<PrimerCount> {
    bed.references()
        .iter()
        .zip(counts)
        .flat_map(|(reference, counts)| {
            reference
                .primers
                .iter()
                .zip(counts)
                .map(|(primer, &count)| PrimerCount {
                    reference: reference.name.clone(),
                    start: primer.start,
                    end: primer.end,
                    name: primer.name.clone(),
                    score: primer.score.clone(),
                    strand: use_strand.then(|| match primer.strand {
                        Some(Strand::Forward) => "+".to_owned(),
                        Some(Strand::Reverse) => "-".to_owned(),
                        None => ".".to_owned(),
                    }),
                    count,
                })
        })
        .collect()
}

impl Summary {
    pub fn write(&self, mut writer: impl Write, command_line: Option<&str>) -> Result<()> {
        if let Some(command_line) = command_line {
            writeln!(writer, "COMMAND: {command_line}")?;
        }
        writeln!(writer, "TOTAL READS: {}", self.total)?;
        writeln!(writer, "TOTAL CLIPPED: {}", self.clipped)?;
        writeln!(writer, "FORWARD CLIPPED: {}", self.forward_clipped)?;
        writeln!(writer, "REVERSE CLIPPED: {}", self.reverse_clipped)?;
        writeln!(writer, "BOTH CLIPPED: {}", self.both_clipped)?;
        writeln!(writer, "NOT CLIPPED: {}", self.not_clipped)?;
        writeln!(writer, "EXCLUDED: {}", self.excluded)?;
        writeln!(writer, "FILTERED: {}", self.filtered)?;
        writeln!(writer, "FAILED: {}", self.failed)?;
        writeln!(writer, "WRITTEN: {}", self.written)?;
        writer.flush().map_err(RsomicsError::Io)
    }
}

pub fn write_primer_counts(counts: &[PrimerCount], mut writer: impl Write) -> Result<()> {
    writeln!(
        writer,
        "#CHR\tLEFT\tRIGHT\tNAME\tSCORE\tSTRAND\tNUM_CLIPPED"
    )?;
    for count in counts {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            count.reference,
            count.start,
            count.end,
            count.name,
            count.score,
            count.strand.as_deref().unwrap_or("."),
            count.count
        )?;
    }
    writer.flush().map_err(RsomicsError::Io)
}
