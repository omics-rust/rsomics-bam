mod merge;
mod record;
mod report;

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use rsomics_bamio::raw::{RawRecord, RawRecordEncoder, RecordRef};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::input;

pub const DEFAULT_FLAG_MASK: u16 = 0x0c1;
pub const DEFAULT_EXCLUDED_FLAGS: u16 = 0x900;
pub const DEFAULT_TAGS: [[u8; 2]; 5] = [*b"BC", *b"FI", *b"QT", *b"RT", *b"TC"];

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Sanitize(u8);

impl Sanitize {
    const POSITION: u8 = 2;
    const MAPPING_QUALITY: u8 = 4;
    const UNMAPPED: u8 = 8;
    const CIGAR: u8 = 16;
    const AUXILIARY: u8 = 32;
    const CIGAR_DUPLICATES: u8 = 64;
    const CIGAR_EQX: u8 = 128;

    pub fn parse(value: &str) -> std::result::Result<Self, String> {
        let mut bits = 0;
        for token in value.split(',') {
            match token {
                "all" | "*" => bits = 127,
                "none" | "off" => bits = 0,
                "on" => {
                    bits = Self::MAPPING_QUALITY
                        | Self::UNMAPPED
                        | Self::CIGAR
                        | Self::AUXILIARY
                        | Self::CIGAR_DUPLICATES
                }
                "pos" => bits |= Self::POSITION,
                "mqual" => bits |= Self::MAPPING_QUALITY,
                "unmap" => bits |= Self::UNMAPPED,
                "cigar" => bits |= Self::CIGAR,
                "aux" => bits |= Self::AUXILIARY,
                "cigdup" => bits |= Self::CIGAR_DUPLICATES,
                "cigarx" => bits |= Self::CIGAR_EQX | Self::CIGAR_DUPLICATES,
                _ => return Err(format!("unrecognized sanitization keyword {token:?}")),
            }
        }
        Ok(Self(bits))
    }

    pub fn all() -> Self {
        Self(255)
    }

    fn contains(self, bit: u8) -> bool {
        self.0 & bit != 0
    }

    fn is_empty(self) -> bool {
        self.0 == 0
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TagSelection {
    Listed(Vec<[u8; 2]>),
    AllExcept(Vec<[u8; 2]>),
}

impl Default for TagSelection {
    fn default() -> Self {
        Self::Listed(DEFAULT_TAGS.to_vec())
    }
}

impl TagSelection {
    pub fn display(&self) -> String {
        let (wildcard, tags) = match self {
            Self::Listed(tags) => (false, tags),
            Self::AllExcept(tags) => (true, tags),
        };
        std::iter::once(wildcard.then_some("*".to_owned()))
            .flatten()
            .chain(
                tags.iter()
                    .map(|tag| String::from_utf8_lossy(tag).into_owned()),
            )
            .collect::<Vec<_>>()
            .join(",")
    }
}

#[derive(Clone, Debug)]
pub struct Options {
    pub required_flags: u16,
    pub excluded_flags: u16,
    pub flag_mask: u16,
    pub reverse_complement: bool,
    pub tags: TagSelection,
    pub order: u8,
    pub check_position: bool,
    pub check_cigar: bool,
    pub check_mate: bool,
    pub maximum_records: Option<u64>,
    pub show_qc: bool,
    pub verbose: bool,
    pub tabs: bool,
    pub bamseqchksum: bool,
    pub sanitize: Sanitize,
    pub additional_threads: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            required_flags: 0,
            excluded_flags: DEFAULT_EXCLUDED_FLAGS,
            flag_mask: DEFAULT_FLAG_MASK,
            reverse_complement: true,
            tags: TagSelection::default(),
            order: 0,
            check_position: false,
            check_cigar: false,
            check_mate: false,
            maximum_records: None,
            show_qc: false,
            verbose: false,
            tabs: false,
            bamseqchksum: false,
            sanitize: Sanitize::default(),
            additional_threads: 0,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct ChecksumValues {
    pub sequence: u64,
    pub name: u64,
    pub quality: u64,
    pub auxiliary: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub position: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cigar: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mate: Option<u64>,
    pub combined: u64,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Row {
    pub qc: Qc,
    pub count: u64,
    pub checksums: ChecksumValues,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Qc {
    All,
    Pass,
    Fail,
}

impl Qc {
    fn label(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::Pass => "pass",
            Self::Fail => "fail",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Group {
    pub name: String,
    pub rows: Vec<Row>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Report {
    pub source: String,
    pub auxiliary_tags: String,
    pub flag_mask: u16,
    pub groups: Vec<Group>,
    #[serde(skip)]
    layout: Layout,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct Layout {
    position: bool,
    cigar: bool,
    mate: bool,
    tabs: bool,
    bamseqchksum: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct RecordChecksums {
    sequence: u32,
    name: u32,
    quality: u32,
    auxiliary: u32,
    position: u32,
    cigar: u32,
    mate: u32,
}

#[derive(Clone, Copy, Debug)]
struct Sums {
    sequence: [u64; 3],
    name: [u64; 3],
    quality: [u64; 3],
    auxiliary: [u64; 3],
    position: [u64; 3],
    cigar: [u64; 3],
    mate: [u64; 3],
    count: [u64; 3],
}

impl Default for Sums {
    fn default() -> Self {
        Self {
            sequence: [1; 3],
            name: [1; 3],
            quality: [1; 3],
            auxiliary: [1; 3],
            position: [1; 3],
            cigar: [1; 3],
            mate: [1; 3],
            count: [0; 3],
        }
    }
}

impl Sums {
    fn update(&mut self, checksums: RecordChecksums, qc_fail: bool, options: &Options, count: u64) {
        let order_crc = if options.order == 0 {
            0
        } else {
            let index = if options.order == 1 {
                count
            } else {
                self.count[0]
            };
            crc32fast::hash(&index.to_le_bytes())
        };
        self.update_row(0, checksums, order_crc, 1);
        if options.show_qc && !qc_fail {
            self.update_row(1, checksums, order_crc, 1);
        }
        if options.show_qc && qc_fail {
            self.update_row(2, checksums, order_crc, 1);
        }
    }

    fn update_row(&mut self, row: usize, checksums: RecordChecksums, order_crc: u32, count: u64) {
        self.sequence[row] = fold(self.sequence[row], order_crc ^ checksums.sequence);
        self.name[row] = fold(self.name[row], order_crc ^ checksums.name);
        self.quality[row] = fold(self.quality[row], order_crc ^ checksums.quality);
        self.auxiliary[row] = fold(self.auxiliary[row], order_crc ^ checksums.auxiliary);
        self.position[row] = fold(self.position[row], order_crc ^ checksums.position);
        self.cigar[row] = fold(self.cigar[row], order_crc ^ checksums.cigar);
        self.mate[row] = fold(self.mate[row], order_crc ^ checksums.mate);
        self.count[row] += count;
    }
}

#[derive(Default)]
struct Accumulator {
    all: Sums,
    without_read_group: Sums,
    read_groups: BTreeMap<Vec<u8>, Sums>,
}

impl Accumulator {
    fn update(
        &mut self,
        checksums: RecordChecksums,
        read_group: Option<&[u8]>,
        qc_fail: bool,
        options: &Options,
    ) {
        let group_count = if let Some(read_group) = read_group {
            if let Some(group) = self.read_groups.get_mut(read_group) {
                let count = group.count[0];
                group.update(checksums, qc_fail, options, count);
                count
            } else {
                let mut group = Sums::default();
                group.update(checksums, qc_fail, options, 0);
                self.read_groups.insert(read_group.to_vec(), group);
                0
            }
        } else {
            let count = self.without_read_group.count[0];
            self.without_read_group
                .update(checksums, qc_fail, options, count);
            count
        };
        self.all.update(checksums, qc_fail, options, group_count);
    }

    fn report(self, source: String, options: &Options) -> Report {
        let layout = Layout {
            position: options.check_position,
            cigar: options.check_cigar,
            mate: options.check_mate,
            tabs: options.tabs || options.bamseqchksum,
            bamseqchksum: options.bamseqchksum,
        };
        let mut groups = vec![group("all", self.all, options, layout)];
        if options.verbose || self.without_read_group.count[0] > 0 || options.bamseqchksum {
            groups.push(group("-", self.without_read_group, options, layout));
        }
        groups.extend(
            self.read_groups
                .into_iter()
                .map(|(name, sums)| group(&String::from_utf8_lossy(&name), sums, options, layout)),
        );
        Report {
            source,
            auxiliary_tags: options.tags.display(),
            flag_mask: options.flag_mask,
            groups,
            layout,
        }
    }

    fn merge_row(&mut self, group_name: &[u8], qc: Qc, checksums: RecordChecksums, count: u64) {
        let row = match qc {
            Qc::All => 0,
            Qc::Pass => 1,
            Qc::Fail => 2,
        };
        let group = if group_name == b"-" {
            &mut self.without_read_group
        } else {
            self.read_groups.entry(group_name.to_vec()).or_default()
        };
        group.update_row(row, checksums, 0, count);
        self.all.update_row(row, checksums, 0, count);
    }
}

pub use merge::merge;

pub fn collect(path: &Path, options: &Options) -> Result<Report> {
    if path == Path::new("-") {
        return collect_standard_input(options);
    }
    if is_sequence_file(path)? {
        return collect_sequence(path, options);
    }
    let mut reader = input::open(path, None, options.additional_threads)?;
    let header = reader.read_header(path)?;
    let mut accumulator = Accumulator::default();
    let mut processed = 0u64;
    let mut scratch = record::Scratch::default();
    let reference_lengths = header
        .reference_sequences()
        .values()
        .map(|reference| i64::try_from(usize::from(reference.length())).unwrap())
        .collect::<Vec<_>>();
    let mut visit = |record: RecordRef<'_>| {
        if let Some(alignment) =
            record::alignment_with_scratch(record, &reference_lengths, options, &mut scratch)?
        {
            accumulator.update(
                alignment.checksums,
                alignment.read_group,
                alignment.qc_fail,
                options,
            );
            processed += 1;
        }
        Ok(options
            .maximum_records
            .is_none_or(|limit| processed < limit))
    };
    if reader.format() == input::Format::Bam {
        reader.visit_raw_bam_records(path, &mut visit)?;
    } else {
        reader.visit_owned_raw_records(&header, path, |record: RawRecord| {
            visit(RecordRef::from_bytes(record.as_bytes())?)
        })?;
    }
    Ok(accumulator.report(path.display().to_string(), options))
}

fn collect_standard_input(options: &Options) -> Result<Report> {
    if options.additional_threads > 0 {
        return Err(RsomicsError::ConfigError(
            "additional decoding threads require a file-backed BAM input".to_owned(),
        ));
    }
    let mut input = BufReader::with_capacity(256 * 1024, io::stdin());
    if sequence_stream_prefix(input.fill_buf().map_err(RsomicsError::Io)?) {
        collect_sequence_stream(input, options)
    } else {
        collect_alignment_stream(input, options)
    }
}

fn collect_sequence_stream(input: impl Read, options: &Options) -> Result<Report> {
    let mut reader = rsomics_seqio::open_reader(input)?;
    let mut accumulator = Accumulator::default();
    let mut processed = 0u64;
    let mut scratch = record::Scratch::default();
    while let Some(sequence) = reader.read_record()? {
        if let Some(checksums) = record::sequence(sequence, options, &mut scratch) {
            accumulator.update(checksums, None, false, options);
            processed += 1;
        }
        if options
            .maximum_records
            .is_some_and(|limit| processed >= limit)
        {
            break;
        }
    }
    Ok(accumulator.report("-".to_owned(), options))
}

fn collect_alignment_stream(input: impl Read + 'static, options: &Options) -> Result<Report> {
    let mut reader = noodles_util::alignment::io::reader::Builder::default()
        .build_from_reader(Box::new(input) as Box<dyn Read>)
        .map_err(|error| {
            RsomicsError::InvalidInput(format!("opening alignment input -: {error}"))
        })?;
    let header = reader.read_header().map_err(|error| {
        RsomicsError::InvalidInput(format!("reading alignment header from -: {error}"))
    })?;
    let reference_lengths = header
        .reference_sequences()
        .values()
        .map(|reference| i64::try_from(usize::from(reference.length())).unwrap())
        .collect::<Vec<_>>();
    let mut accumulator = Accumulator::default();
    let mut encoder = RawRecordEncoder::new();
    let mut processed = 0u64;
    let mut scratch = record::Scratch::default();
    for result in reader.records(&header) {
        let source = result.map_err(|error| {
            RsomicsError::InvalidInput(format!("reading alignment record from -: {error}"))
        })?;
        let raw = encoder.encode(&header, source.as_ref())?;
        let record = RecordRef::from_bytes(raw.as_bytes())?;
        if let Some(alignment) =
            record::alignment_with_scratch(record, &reference_lengths, options, &mut scratch)?
        {
            accumulator.update(
                alignment.checksums,
                alignment.read_group,
                alignment.qc_fail,
                options,
            );
            processed += 1;
        }
        if options
            .maximum_records
            .is_some_and(|limit| processed >= limit)
        {
            break;
        }
    }
    Ok(accumulator.report("-".to_owned(), options))
}

fn sequence_prefix(bytes: &[u8]) -> bool {
    if bytes.first() == Some(&b'>') {
        return true;
    }
    if bytes.first() != Some(&b'@') {
        return false;
    }
    let first_line = bytes.split(|&byte| byte == b'\n').next().unwrap_or(bytes);
    ![b"@HD\t", b"@SQ\t", b"@RG\t", b"@PG\t", b"@CO\t"]
        .iter()
        .any(|prefix| first_line.starts_with(*prefix))
        && first_line.iter().filter(|&&byte| byte == b'\t').count() < 10
}

fn collect_sequence(path: &Path, options: &Options) -> Result<Report> {
    if options.additional_threads > 0 {
        return Err(RsomicsError::ConfigError(
            "additional decoding threads require BAM input".to_owned(),
        ));
    }
    let mut reader = rsomics_seqio::open_path(path)?;
    let mut accumulator = Accumulator::default();
    let mut processed = 0u64;
    let mut scratch = record::Scratch::default();
    while let Some(sequence) = reader.read_record()? {
        if let Some(checksums) = record::sequence(sequence, options, &mut scratch) {
            accumulator.update(checksums, None, false, options);
            processed += 1;
        }
        if options
            .maximum_records
            .is_some_and(|limit| processed >= limit)
        {
            break;
        }
    }
    Ok(accumulator.report(path.display().to_string(), options))
}

fn is_sequence_file(path: &Path) -> Result<bool> {
    let file = File::open(path).map_err(|error| {
        RsomicsError::Io(io::Error::new(
            error.kind(),
            format!("opening {}: {error}", path.display()),
        ))
    })?;
    let mut input = BufReader::with_capacity(256 * 1024, file);
    Ok(sequence_stream_prefix(
        input.fill_buf().map_err(RsomicsError::Io)?,
    ))
}

fn sequence_stream_prefix(bytes: &[u8]) -> bool {
    if bytes.starts_with(&[0x1f, 0x8b]) {
        let decoder = flate2::read::GzDecoder::new(bytes);
        let mut decoded = Vec::with_capacity(64 * 1024);
        return decoder
            .take(64 * 1024)
            .read_to_end(&mut decoded)
            .is_ok_and(|_| sequence_prefix(&decoded));
    }
    sequence_prefix(bytes)
}

fn group(name: &str, sums: Sums, options: &Options, layout: Layout) -> Group {
    let mut rows = Vec::new();
    for (index, qc) in [(0, Qc::All), (1, Qc::Pass), (2, Qc::Fail)] {
        if index > 0 && !options.show_qc && !(options.bamseqchksum && index == 1) {
            continue;
        }
        if !options.verbose && sums.count[index] == 0 && !options.bamseqchksum {
            continue;
        }
        let mut combined = 1;
        for value in [
            sums.count[index] >> 32,
            sums.count[index] & u64::from(u32::MAX),
            sums.sequence[index],
            sums.name[index],
            sums.sequence[index],
            sums.auxiliary[index],
        ] {
            combined = fold(combined, value as u32);
        }
        for value in [
            layout.position.then_some(sums.position[index]),
            layout.cigar.then_some(sums.cigar[index]),
            layout.mate.then_some(sums.mate[index]),
        ]
        .into_iter()
        .flatten()
        {
            combined = fold(combined, value as u32);
        }
        rows.push(Row {
            qc,
            count: sums.count[index],
            checksums: ChecksumValues {
                sequence: sums.sequence[index],
                name: sums.name[index],
                quality: sums.quality[index],
                auxiliary: sums.auxiliary[index],
                position: layout.position.then_some(sums.position[index]),
                cigar: layout.cigar.then_some(sums.cigar[index]),
                mate: layout.mate.then_some(sums.mate[index]),
                combined,
            },
        });
    }
    Group {
        name: name.to_owned(),
        rows,
    }
}

fn fold(hash: u64, crc: u32) -> u64 {
    const PRIME: u64 = (1u64 << 31) - 1;
    let mut value = u64::from(crc) & PRIME;
    if value == 0 || value == PRIME {
        value = 1;
    }
    hash * value % PRIME
}
