use std::cmp::Ordering;
use std::collections::HashMap;
use std::mem;
use std::sync::Arc;

use noodles::sam;
use noodles::sam::header::record::value::{
    Map,
    map::{self, header::tag as header_tag, read_group::tag as read_group_tag},
};
use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

const FLAG_PAIRED: u16 = 0x01;
const FLAG_UNMAPPED: u16 = 0x04;
const FLAG_MATE_UNMAPPED: u16 = 0x08;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_MATE_REVERSE: u16 = 0x20;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Order {
    #[default]
    Coordinate,
    QueryNameNatural,
    QueryNameLexicographical,
    TemplateCoordinate,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum InternalOrder {
    Coordinate,
    QueryNameNatural,
    QueryNameLexicographical,
    TemplateCoordinate,
    Collate,
}

impl From<Order> for InternalOrder {
    fn from(order: Order) -> Self {
        match order {
            Order::Coordinate => Self::Coordinate,
            Order::QueryNameNatural => Self::QueryNameNatural,
            Order::QueryNameLexicographical => Self::QueryNameLexicographical,
            Order::TemplateCoordinate => Self::TemplateCoordinate,
        }
    }
}

#[derive(Clone)]
pub(crate) struct OrderedRecord {
    pub(crate) record: RawRecord,
    key: EntryKey,
    pub(crate) ordinal: u64,
}

#[derive(Clone)]
enum EntryKey {
    Coordinate((u32, u32, bool)),
    QueryName,
    Template(TemplateKey),
    Collate(u32),
}

#[derive(Clone)]
struct TemplateKey {
    tid1: i32,
    tid2: i32,
    pos1: i64,
    pos2: i64,
    neg1: bool,
    neg2: bool,
    library: Arc<[u8]>,
    upper: bool,
}

impl OrderedRecord {
    pub(crate) fn memory(&self) -> Result<u64> {
        let payload = u64::try_from(self.record.as_bytes().len())
            .map_err(|_| RsomicsError::InvalidInput("record size exceeds u64".to_owned()))?;
        let key = match &self.key {
            EntryKey::Template(key) => u64::try_from(key.library.len()).unwrap_or(u64::MAX),
            EntryKey::Coordinate(_) | EntryKey::QueryName | EntryKey::Collate(_) => 0,
        };
        payload
            .checked_add(mem::size_of::<Self>() as u64)
            .and_then(|value| value.checked_add(key))
            .ok_or_else(|| {
                RsomicsError::InvalidInput("sort memory accounting overflows".to_owned())
            })
    }
}

pub(crate) fn ordered_record(
    record: RawRecord,
    order: InternalOrder,
    header: &sam::Header,
    libraries: &HashMap<Vec<u8>, Arc<[u8]>>,
    ordinal: u64,
) -> Result<OrderedRecord> {
    validate_record_coordinates(&record, header.reference_sequences().len())?;
    let key = match order {
        InternalOrder::Coordinate => EntryKey::Coordinate(coordinate_key(&record)),
        InternalOrder::QueryNameNatural | InternalOrder::QueryNameLexicographical => {
            EntryKey::QueryName
        }
        InternalOrder::TemplateCoordinate => EntryKey::Template(template_key(&record, libraries)?),
        InternalOrder::Collate => EntryKey::Collate(collate_hash(record.name())),
    };
    Ok(OrderedRecord {
        record,
        key,
        ordinal,
    })
}

pub(crate) fn validate_record_coordinates(
    record: &RawRecord,
    reference_count: usize,
) -> Result<()> {
    validate_coordinate_fields(
        record.reference_sequence_id(),
        record.alignment_start(),
        record.mate_reference_sequence_id(),
        record.mate_alignment_start(),
        reference_count,
    )
}

pub(crate) fn validate_coordinate_fields(
    reference: i32,
    position: i32,
    mate_reference: i32,
    mate_position: i32,
    reference_count: usize,
) -> Result<()> {
    validate_reference_id(reference, reference_count, "reference")?;
    validate_reference_id(mate_reference, reference_count, "mate reference")?;
    if position < -1 || mate_position < -1 {
        return Err(RsomicsError::InvalidInput(
            "alignment position is below -1".to_owned(),
        ));
    }
    Ok(())
}

fn validate_reference_id(value: i32, count: usize, field: &str) -> Result<()> {
    if value < -1 || usize::try_from(value).is_ok_and(|value| value >= count) {
        return Err(RsomicsError::InvalidInput(format!(
            "{field} sequence ID {value} is outside the header dictionary"
        )));
    }
    Ok(())
}

pub(crate) fn compare_ordered_records(
    order: InternalOrder,
    a: &OrderedRecord,
    b: &OrderedRecord,
) -> Ordering {
    match (order, &a.key, &b.key) {
        (InternalOrder::Coordinate, EntryKey::Coordinate(a), EntryKey::Coordinate(b)) => a.cmp(b),
        (InternalOrder::QueryNameNatural, EntryKey::QueryName, EntryKey::QueryName) => {
            natural_cmp(a.record.name(), b.record.name())
                .then_with(|| name_flag_key(a.record.flags()).cmp(&name_flag_key(b.record.flags())))
        }
        (InternalOrder::QueryNameLexicographical, EntryKey::QueryName, EntryKey::QueryName) => a
            .record
            .name()
            .cmp(b.record.name())
            .then_with(|| name_flag_key(a.record.flags()).cmp(&name_flag_key(b.record.flags()))),
        (
            InternalOrder::TemplateCoordinate,
            EntryKey::Template(a_key),
            EntryKey::Template(b_key),
        ) => compare_template(a_key, &a.record, b_key, &b.record),
        (InternalOrder::Collate, EntryKey::Collate(a_key), EntryKey::Collate(b_key)) => a_key
            .cmp(b_key)
            .then_with(|| a.record.name().cmp(b.record.name()))
            .then_with(|| ((a.record.flags() >> 6) & 3).cmp(&((b.record.flags() >> 6) & 3))),
        _ => unreachable!("entry key matches the selected order"),
    }
}

fn collate_hash(name: &[u8]) -> u32 {
    let Some((&first, rest)) = name.split_first() else {
        return 0;
    };
    let mut hash = u32::from(first);
    for &byte in rest {
        hash = hash.wrapping_mul(31).wrapping_add(u32::from(byte));
    }
    hash = hash.wrapping_add(!(hash << 15));
    hash ^= hash >> 10;
    hash = hash.wrapping_add(hash << 3);
    hash ^= hash >> 6;
    hash = hash.wrapping_add(!(hash << 11));
    hash ^ (hash >> 16)
}

fn coordinate_key(record: &RawRecord) -> (u32, u32, bool) {
    let tid = u32::try_from(record.reference_sequence_id()).unwrap_or(u32::MAX);
    let pos = u32::try_from(i64::from(record.alignment_start()) + 1).unwrap_or(0);
    (tid, pos, record.flags() & FLAG_REVERSE != 0)
}

fn name_flag_key(flags: u16) -> u32 {
    (u32::from(flags & 0x00c0) << 8)
        | (u32::from(flags & 0x0100) << 3)
        | (u32::from(flags & 0x0800) >> 3)
}

fn natural_cmp(a: &[u8], b: &[u8]) -> Ordering {
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if !a[i].is_ascii_digit() || !b[j].is_ascii_digit() {
            match a[i].cmp(&b[j]) {
                Ordering::Equal => {
                    i += 1;
                    j += 1;
                }
                order => return order,
            }
            continue;
        }
        while i < a.len() && a[i] == b'0' {
            i += 1;
        }
        while j < b.len() && b[j] == b'0' {
            j += 1;
        }
        while i < a.len() && j < b.len() && a[i].is_ascii_digit() && a[i] == b[j] {
            i += 1;
            j += 1;
        }
        let difference = match (a.get(i), b.get(j)) {
            (Some(left), Some(right)) => left.cmp(right),
            (Some(_), None) => Ordering::Greater,
            (None, Some(_)) => Ordering::Less,
            (None, None) => Ordering::Equal,
        };
        while i < a.len() && j < b.len() && a[i].is_ascii_digit() && b[j].is_ascii_digit() {
            i += 1;
            j += 1;
        }
        if i < a.len() && a[i].is_ascii_digit() {
            return Ordering::Greater;
        }
        if j < b.len() && b[j].is_ascii_digit() {
            return Ordering::Less;
        }
        if difference != Ordering::Equal {
            return difference;
        }
    }
    (a.len() - i).cmp(&(b.len() - j))
}

pub(crate) fn library_lookup(header: &sam::Header) -> HashMap<Vec<u8>, Arc<[u8]>> {
    header
        .read_groups()
        .iter()
        .filter_map(|(id, read_group)| {
            read_group
                .other_fields()
                .get(&read_group_tag::LIBRARY)
                .map(|library| {
                    (
                        id.as_slice().to_vec(),
                        Arc::<[u8]>::from(library.as_slice()),
                    )
                })
        })
        .collect()
}

fn template_key(
    record: &RawRecord,
    libraries: &HashMap<Vec<u8>, Arc<[u8]>>,
) -> Result<TemplateKey> {
    let mut key = TemplateKey {
        tid1: i32::MAX,
        tid2: i32::MAX,
        pos1: i64::MAX,
        pos2: i64::MAX,
        neg1: false,
        neg2: false,
        library: z_aux(record, *b"RG", "RG")?
            .and_then(|id| libraries.get(id))
            .cloned()
            .unwrap_or_else(|| Arc::from([])),
        upper: false,
    };

    if record.flags() & FLAG_UNMAPPED == 0 {
        key.tid1 = record.reference_sequence_id();
        key.neg1 = record.flags() & FLAG_REVERSE != 0;
        key.pos1 = if key.neg1 {
            unclipped_end(record)?
        } else {
            unclipped_start(record)?
        };
    }
    if record.flags() & FLAG_PAIRED != 0 && record.flags() & FLAG_MATE_UNMAPPED == 0 {
        let cigar = z_aux(record, *b"MC", "MC")?.ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "read {} has a mapped mate but no MC tag",
                String::from_utf8_lossy(record.name())
            ))
        })?;
        let cigar = parse_cigar(cigar)?;
        key.tid2 = record.mate_reference_sequence_id();
        key.neg2 = record.flags() & FLAG_MATE_REVERSE != 0;
        key.pos2 = if key.neg2 {
            unclipped_other_end(record.mate_alignment_start(), &cigar)?
        } else {
            unclipped_other_start(record.mate_alignment_start(), &cigar)?
        };
    }
    z_aux(record, *b"CB", "CB")?;
    z_aux(record, *b"MI", "MI")?;

    let current_first = key.tid1 < key.tid2
        || (key.tid1 == key.tid2 && key.pos1 < key.pos2)
        || (key.tid1 == key.tid2 && key.pos1 == key.pos2 && !key.neg1);
    if !current_first {
        key.upper = true;
        mem::swap(&mut key.tid1, &mut key.tid2);
        mem::swap(&mut key.pos1, &mut key.pos2);
        mem::swap(&mut key.neg1, &mut key.neg2);
    }
    Ok(key)
}

fn compare_template(
    a: &TemplateKey,
    a_record: &RawRecord,
    b: &TemplateKey,
    b_record: &RawRecord,
) -> Ordering {
    a.tid1
        .cmp(&b.tid1)
        .then_with(|| a.tid2.cmp(&b.tid2))
        .then_with(|| a.pos1.cmp(&b.pos1))
        .then_with(|| a.pos2.cmp(&b.pos2))
        .then_with(|| b.neg1.cmp(&a.neg1))
        .then_with(|| b.neg2.cmp(&a.neg2))
        .then_with(|| a.library.cmp(&b.library))
        .then_with(|| aux_or_empty(a_record, *b"CB").cmp(aux_or_empty(b_record, *b"CB")))
        .then_with(|| {
            molecular_cmp(
                aux_or_empty(a_record, *b"MI"),
                aux_or_empty(b_record, *b"MI"),
            )
        })
        .then_with(|| a_record.name().cmp(b_record.name()))
        .then_with(|| a.upper.cmp(&b.upper))
}

fn aux_or_empty(record: &RawRecord, tag: [u8; 2]) -> &[u8] {
    record
        .aux_value(tag)
        .and_then(|value| value.strip_suffix(&[0]))
        .unwrap_or_default()
}

fn z_aux<'a>(record: &'a RawRecord, tag: [u8; 2], label: &str) -> Result<Option<&'a [u8]>> {
    match (record.aux_type(tag), record.aux_value(tag)) {
        (None, None) => Ok(None),
        (Some(b'Z'), Some(value)) => value.strip_suffix(&[0]).map(Some).ok_or_else(|| {
            RsomicsError::InvalidInput(format!("read has unterminated {label}:Z tag"))
        }),
        (Some(_), Some(_)) => Err(RsomicsError::InvalidInput(format!(
            "read {} has a non-string {label} tag",
            String::from_utf8_lossy(record.name())
        ))),
        _ => Err(RsomicsError::InvalidInput(format!(
            "read {} has a malformed {label} tag",
            String::from_utf8_lossy(record.name())
        ))),
    }
}

fn molecular_cmp(a: &[u8], b: &[u8]) -> Ordering {
    fn trim(value: &[u8]) -> &[u8] {
        if value.len() >= 2 && value[value.len() - 2] == b'/' {
            &value[..value.len() - 2]
        } else {
            value
        }
    }
    trim(a).cmp(trim(b))
}

fn unclipped_start(record: &RawRecord) -> Result<i64> {
    let cigar = record.decoded_cigar()?;
    let clipped = leading_soft(&cigar)?;
    i64::from(record.alignment_start())
        .checked_add(1)
        .and_then(|value| value.checked_sub(clipped))
        .ok_or_else(|| RsomicsError::InvalidInput("unclipped start overflows".to_owned()))
}

fn unclipped_end(record: &RawRecord) -> Result<i64> {
    let cigar = record.decoded_cigar()?;
    let span = reference_span(&cigar)?.max(1);
    let clipped = trailing_soft(&cigar)?;
    i64::from(record.alignment_start())
        .checked_add(span)
        .and_then(|value| value.checked_add(clipped))
        .ok_or_else(|| RsomicsError::InvalidInput("unclipped end overflows".to_owned()))
}

fn leading_soft(cigar: &[(u8, u32)]) -> Result<i64> {
    let mut clipped = 0i64;
    for &(kind, length) in cigar {
        match kind {
            4 => {
                clipped = clipped
                    .checked_add(i64::from(length))
                    .ok_or_else(cigar_overflow)?
            }
            5 => {}
            _ => break,
        }
    }
    Ok(clipped)
}

fn trailing_soft(cigar: &[(u8, u32)]) -> Result<i64> {
    let mut clipped = 0i64;
    for &(kind, length) in cigar.iter().rev() {
        match kind {
            4 => {
                clipped = clipped
                    .checked_add(i64::from(length))
                    .ok_or_else(cigar_overflow)?
            }
            5 => {}
            _ => break,
        }
    }
    Ok(clipped)
}

fn reference_span(cigar: &[(u8, u32)]) -> Result<i64> {
    cigar
        .iter()
        .filter(|(kind, _)| matches!(kind, 0 | 2 | 3 | 7 | 8))
        .try_fold(0i64, |span, (_, length)| {
            span.checked_add(i64::from(*length))
                .ok_or_else(cigar_overflow)
        })
}

fn parse_cigar(value: &[u8]) -> Result<Vec<(u8, u64)>> {
    if value == b"*" {
        return Ok(Vec::new());
    }
    let mut cigar = Vec::new();
    let mut length = None::<u64>;
    for &byte in value {
        if byte.is_ascii_digit() {
            let digit = u64::from(byte - b'0');
            length = Some(
                length
                    .unwrap_or_default()
                    .checked_mul(10)
                    .and_then(|value| value.checked_add(digit))
                    .ok_or_else(cigar_overflow)?,
            );
            continue;
        }
        let length = length.take().filter(|length| *length > 0).ok_or_else(|| {
            RsomicsError::InvalidInput("MC tag contains an invalid CIGAR length".to_owned())
        })?;
        let kind = match byte {
            b'M' => 0,
            b'I' => 1,
            b'D' => 2,
            b'N' => 3,
            b'S' => 4,
            b'H' => 5,
            b'P' => 6,
            b'=' => 7,
            b'X' => 8,
            _ => {
                return Err(RsomicsError::InvalidInput(format!(
                    "MC tag contains unsupported CIGAR operation {}",
                    char::from(byte)
                )));
            }
        };
        cigar.push((kind, length));
    }
    if length.is_some() || cigar.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "MC tag contains an incomplete CIGAR".to_owned(),
        ));
    }
    Ok(cigar)
}

fn unclipped_other_start(position: i32, cigar: &[(u8, u64)]) -> Result<i64> {
    let mut clipped = 0i64;
    for &(kind, length) in cigar {
        match kind {
            4 => {
                clipped = clipped
                    .checked_add(i64::try_from(length).map_err(|_| cigar_overflow())?)
                    .ok_or_else(cigar_overflow)?;
            }
            5 => {}
            _ => break,
        }
    }
    i64::from(position)
        .checked_add(1)
        .and_then(|value| value.checked_sub(clipped))
        .ok_or_else(cigar_overflow)
}

fn unclipped_other_end(position: i32, cigar: &[(u8, u64)]) -> Result<i64> {
    let mut span = 0i64;
    let mut initial_clips = true;
    for &(kind, length) in cigar {
        let length = i64::try_from(length).map_err(|_| cigar_overflow())?;
        match kind {
            0 | 2 | 3 | 7 | 8 => {
                span = span.checked_add(length).ok_or_else(cigar_overflow)?;
                initial_clips = false;
            }
            4 if !initial_clips => {
                span = span.checked_add(length).ok_or_else(cigar_overflow)?;
            }
            _ => {}
        }
    }
    i64::from(position)
        .checked_add(span)
        .ok_or_else(cigar_overflow)
}

pub(crate) fn set_sort_order(header: &mut sam::Header, order: Order) {
    let hd = header
        .header_mut()
        .get_or_insert_with(Map::<map::Header>::default);
    let fields = hd.other_fields_mut();
    let (sort, group, subsort) = match order {
        Order::Coordinate => ("coordinate", None, None),
        Order::QueryNameNatural => ("queryname", None, Some("queryname:natural")),
        Order::QueryNameLexicographical => ("queryname", None, Some("queryname:lexicographical")),
        Order::TemplateCoordinate => (
            "unsorted",
            Some("query"),
            Some("unsorted:template-coordinate"),
        ),
    };
    fields.insert(header_tag::SORT_ORDER, sort.into());
    match group {
        Some(group) => {
            fields.insert(header_tag::GROUP_ORDER, group.into());
        }
        None => {
            fields.shift_remove(&header_tag::GROUP_ORDER);
        }
    }
    match subsort {
        Some(subsort) => {
            fields.insert(header_tag::SUBSORT_ORDER, subsort.into());
        }
        None => {
            fields.shift_remove(&header_tag::SUBSORT_ORDER);
        }
    }
}

pub(crate) fn set_collate_order(header: &mut sam::Header) {
    let hd = header
        .header_mut()
        .get_or_insert_with(Map::<map::Header>::default);
    let fields = hd.other_fields_mut();
    fields.insert(header_tag::SORT_ORDER, "unsorted".into());
    fields.insert(header_tag::GROUP_ORDER, "query".into());
}
fn cigar_overflow() -> RsomicsError {
    RsomicsError::InvalidInput("CIGAR coordinate calculation overflows".to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn natural_names_match_samtools_numeric_order() {
        assert_eq!(natural_cmp(b"read7", b"read12"), Ordering::Less);
        assert_eq!(natural_cmp(b"read01", b"read1"), Ordering::Equal);
        assert_eq!(natural_cmp(b"read2a", b"read2b"), Ordering::Less);
        assert_eq!(natural_cmp(b"read10", b"read2"), Ordering::Greater);
    }

    #[test]
    fn queryname_flags_order_read_ends_and_alignment_classes() {
        assert!(name_flag_key(0x40) < name_flag_key(0x80));
        assert!(name_flag_key(0) < name_flag_key(0x800));
        assert!(name_flag_key(0x800) < name_flag_key(0x100));
    }

    #[test]
    fn collate_hash_matches_samtools_1_24() {
        assert_eq!(collate_hash(b""), 0);
        assert_eq!(collate_hash(b"read1"), 4_022_420_600);
        assert_eq!(collate_hash(b"pair2"), 2_054_658_801);
        assert_eq!(collate_hash(b"pair12"), 2_390_552_450);
        assert_eq!(collate_hash(b"unmapped10"), 1_792_842_080);
    }

    #[test]
    fn parses_template_cigars() {
        assert_eq!(
            parse_cigar(b"5S10M2D3M").unwrap(),
            vec![(4, 5), (0, 10), (2, 2), (0, 3)]
        );
        assert!(parse_cigar(b"10M0S").is_err());
        assert!(parse_cigar(b"10").is_err());
    }

    #[test]
    fn template_molecular_suffix_is_ignored() {
        assert_eq!(molecular_cmp(b"molecule/1", b"molecule/2"), Ordering::Equal);
    }
}
