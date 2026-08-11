use std::io::Write;
use std::path::Path;

use noodles::sam;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::header::record::value::{Map, map};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::{Program, input, output};

mod cram;
mod raw;

const FLAG_PAIRED: u16 = 0x01;
const FLAG_PROPER_PAIR: u16 = 0x02;
const FLAG_UNMAPPED: u16 = 0x04;
const FLAG_MATE_UNMAPPED: u16 = 0x08;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_MATE_REVERSE: u16 = 0x20;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_DUPLICATE: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

const DEFAULT_REMOVE_TAGS: [[u8; 2]; 15] = [
    *b"AS", *b"CC", *b"CG", *b"CP", *b"H1", *b"H2", *b"HI", *b"H0", *b"IH", *b"MC", *b"MD", *b"MQ",
    *b"NM", *b"SA", *b"TS",
];

#[derive(Clone, Copy)]
struct TagFilter<'a> {
    remove: &'a [[u8; 2]],
    keep: Option<&'a [[u8; 2]]>,
    no_rg: bool,
}

impl<'a> TagFilter<'a> {
    fn new(remove: &'a [[u8; 2]], keep: Option<&'a [[u8; 2]]>, no_rg: bool) -> Self {
        Self {
            remove,
            keep,
            no_rg,
        }
    }

    fn remove(self, tag: [u8; 2]) -> bool {
        self.no_rg && tag == *b"RG"
            || self.keep.map_or_else(
                || DEFAULT_REMOVE_TAGS.contains(&tag) || self.remove.contains(&tag),
                |keep| !keep.contains(&tag),
            )
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    Sam,
    Bam,
    Cram,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Options<'a> {
    pub format: Format,
    pub remove_tags: &'a [[u8; 2]],
    pub keep_tags: Option<&'a [[u8; 2]]>,
    pub no_rg: bool,
    pub reject_pg: Option<&'a str>,
    pub keep_duplicate: bool,
    pub reference: Option<&'a Path>,
    pub additional_threads: usize,
    pub program: Option<Program<'a>>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub written: u64,
}

pub fn write<W>(input_path: &Path, options: Options<'_>, output: W) -> Result<Summary>
where
    W: Write + Send + 'static,
{
    let input_threads = if input_path == Path::new("-")
        || options.format == Format::Bam
        || input::detect_format(input_path)? == input::Format::Cram
    {
        0
    } else {
        options.additional_threads
    };
    let mut reader = input::open(input_path, options.reference, input_threads)?;
    let input_format = reader.format();
    let input_header = reader.read_header(input_path)?;
    let mut output_header = reset_header(&input_header, options.no_rg, options.reject_pg);
    if let Some(program) = options.program {
        program.add_to(&mut output_header)?;
    }

    let output_format = match options.format {
        Format::Sam => output::Format::Sam,
        Format::Bam => output::Format::Bam,
        Format::Cram => {
            return Err(RsomicsError::ConfigError(
                "CRAM output requires a file-backed destination".to_owned(),
            ));
        }
    };
    let mut writer = output::Writer::new(
        output_format,
        output::Compression::Default,
        options.additional_threads,
        output,
    );
    writer.write_header(&output_header)?;

    let tags = TagFilter::new(options.remove_tags, options.keep_tags, options.no_rg);
    let mut written = 0u64;
    if input_format == input::Format::Bam
        && options.format == Format::Bam
        && input_path != Path::new("-")
    {
        let mut transformed = Vec::new();
        reader.visit_raw_bam_records(input_path, |source| {
            if raw::reset(&source, tags, options.keep_duplicate, &mut transformed)? {
                let record = rsomics_bamio::raw::RecordRef::from_bytes(&transformed)?;
                writer.write_raw_record(&record)?;
                written = increment_record_count(written)?;
            }
            Ok(true)
        })?;
    } else {
        reader.visit_records(&input_header, input_path, |source| {
            let mut record =
                sam::alignment::RecordBuf::try_from_alignment_record(&input_header, source)
                    .map_err(RsomicsError::Io)?;
            if reset_record(
                &mut record,
                tags,
                options.keep_duplicate,
                input_format == input::Format::Cram,
            ) {
                writer.write_record(&output_header, &record)?;
                written = increment_record_count(written)?;
            }
            Ok(true)
        })?;
    }
    writer.finish(&output_header)?;

    Ok(Summary { written })
}

pub fn write_cram_path(
    input_path: &Path,
    options: Options<'_>,
    output_path: &Path,
) -> Result<Summary> {
    cram::write(input_path, options, output_path)
}

fn reset_header(input: &sam::Header, no_rg: bool, reject_pg: Option<&str>) -> sam::Header {
    let mut output = input.clone();
    *output.header_mut() = Some(Map::<map::Header>::default());
    output.reference_sequences_mut().clear();
    output.comments_mut().clear();
    if no_rg {
        output.read_groups_mut().clear();
    }
    if let Some(reject_pg) = reject_pg {
        let programs = output.programs_mut().as_mut();
        let rejected = programs
            .keys()
            .skip_while(|id| AsRef::<[u8]>::as_ref(*id) != reject_pg.as_bytes())
            .cloned()
            .collect::<Vec<_>>();
        for id in rejected {
            programs.shift_remove(AsRef::<[u8]>::as_ref(&id));
        }
    }
    output
}

fn reset_record(
    record: &mut sam::alignment::RecordBuf,
    tags: TagFilter<'_>,
    keep_duplicate: bool,
    cram_input: bool,
) -> bool {
    let flags = u16::from(record.flags());
    if flags & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
        return false;
    }

    let (reset_flags, reverse) = transformed_flags(flags, keep_duplicate);
    if reverse {
        record.sequence_mut().as_mut().reverse();
        for base in record.sequence_mut().as_mut() {
            *base = complement(*base);
        }
        record.quality_scores_mut().as_mut().reverse();
    }
    *record.flags_mut() = sam::alignment::record::Flags::from(reset_flags);
    *record.reference_sequence_id_mut() = None;
    *record.alignment_start_mut() = None;
    *record.mapping_quality_mut() = Some(sam::alignment::record::MappingQuality::MIN);
    record.cigar_mut().as_mut().clear();
    *record.mate_reference_sequence_id_mut() = None;
    *record.mate_alignment_start_mut() = None;
    *record.template_length_mut() = 0;

    let mut retained = record
        .data()
        .iter()
        .filter(|(tag, _)| !tags.remove((*tag).into()))
        .map(|(tag, value)| (tag, value.clone()))
        .collect::<Vec<_>>();
    if cram_input && let Some(index) = retained.iter().position(|(tag, _)| *tag == Tag::READ_GROUP)
    {
        let read_group = retained.remove(index);
        retained.push(read_group);
    }
    *record.data_mut() = retained.into_iter().collect();
    true
}

fn transformed_flags(flags: u16, keep_duplicate: bool) -> (u16, bool) {
    let mut output = (flags & !(FLAG_PROPER_PAIR | FLAG_MATE_REVERSE)) | FLAG_UNMAPPED;
    if !keep_duplicate {
        output &= !FLAG_DUPLICATE;
    }
    if flags & FLAG_PAIRED != 0 {
        output |= FLAG_MATE_UNMAPPED;
    }
    let reverse = flags & FLAG_REVERSE != 0;
    if reverse {
        output &= !FLAG_REVERSE;
    }
    (output, reverse)
}

fn increment_record_count(count: u64) -> Result<u64> {
    count
        .checked_add(1)
        .ok_or_else(|| RsomicsError::InvalidInput("reset record count exceeds u64".to_owned()))
}

fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'M' => b'K',
        b'R' => b'Y',
        b'S' => b'S',
        b'V' => b'B',
        b'W' => b'W',
        b'Y' => b'R',
        b'H' => b'D',
        b'K' => b'M',
        b'D' => b'H',
        b'B' => b'V',
        other => other,
    }
}
