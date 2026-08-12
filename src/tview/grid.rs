use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use rsomics_pileup::{Column, ColumnEntry};

use super::call::Caller;
use super::rows::ReadState;
use super::{Cell, CellColor, Settings, Viewport, line_count};

const REVERSE: u16 = 0x10;
const PAIRED: u16 = 0x01;
const PROPER_PAIR: u16 = 0x02;
const SECONDARY: u16 = 0x100;

#[derive(Clone, Debug)]
pub(super) struct Entry {
    row: usize,
    flags: u16,
    name_symbol: Option<u8>,
    color_symbol: Option<u8>,
    color_matches: bool,
    color_quality: Option<u8>,
    pub(super) base: u8,
    pub(super) quality: u8,
    pub(super) mapping_quality: u8,
    pub(super) reverse: bool,
    pub(super) deletion: bool,
    pub(super) reference_skip: bool,
    insertion: Vec<Inserted>,
}

#[derive(Clone, Debug)]
struct Inserted {
    base: Option<u8>,
    name_symbol: Option<u8>,
    color_symbol: Option<u8>,
    color_matches: bool,
}

#[derive(Clone, Debug)]
pub(super) struct PositionData {
    position: i64,
    entries: Vec<Entry>,
    insertion_width: usize,
}

pub(super) fn own_column(
    column: &Column<'_, ReadState>,
    position: i64,
    rows: &[usize],
) -> PositionData {
    let entries = column
        .entries()
        .zip(rows)
        .map(|(entry, row)| own_entry(entry, column.position(), *row))
        .collect::<Vec<_>>();
    let insertion_width = entries
        .iter()
        .map(|entry| entry.insertion.len())
        .max()
        .unwrap_or(0);
    PositionData {
        position,
        entries,
        insertion_width,
    }
}

fn own_entry(entry: ColumnEntry<'_, ReadState>, position: i64, row: usize) -> Entry {
    let projection = entry.projection();
    let qualities = entry.record().quality_scores();
    let reverse = entry.record().flags() & REVERSE != 0;
    let hard_clip = entry
        .cigar()
        .next()
        .filter(|(kind, _)| *kind == 5)
        .map(|(_, length)| length as usize)
        .unwrap_or(0);
    let color = color_data(entry.record(), projection.qpos, reverse, hard_clip);
    let quality = if projection.is_deletion {
        let previous = projection
            .qpos
            .checked_sub(1)
            .and_then(|position| qualities.get(position))
            .copied()
            .unwrap_or(255);
        let next = qualities.get(projection.qpos).copied().unwrap_or(255);
        previous.min(next)
    } else {
        qualities.get(projection.qpos).copied().unwrap_or(255)
    };
    Entry {
        row,
        flags: entry.record().flags(),
        name_symbol: name_symbol(entry.record(), projection.qpos),
        color_symbol: color.symbol,
        color_matches: color.matches,
        color_quality: color.quality,
        base: if projection.is_deletion {
            16
        } else {
            entry.record().seq_nibble(projection.qpos)
        },
        quality,
        mapping_quality: entry.record().mapping_quality(),
        reverse,
        deletion: projection.is_deletion,
        reference_skip: projection.is_reference_skip,
        insertion: insertion(entry, position, hard_clip),
    }
}

fn insertion(entry: ColumnEntry<'_, ReadState>, position: i64, hard_clip: usize) -> Vec<Inserted> {
    let projection = entry.projection();
    if !at_cigar_end(entry, position) {
        return Vec::new();
    }
    let mut query = projection.qpos + usize::from(!projection.is_deletion);
    let mut bases = Vec::new();
    let mut has_insertion = false;
    for (kind, length) in entry.cigar().skip(projection.cigar_index + 1) {
        if !matches!(kind, 1 | 6) {
            break;
        }
        for _ in 0..length {
            let display_query = projection
                .qpos
                .checked_add(bases.len())
                .and_then(|query| query.checked_add(1));
            let name_symbol = display_query.and_then(|query| name_symbol(entry.record(), query));
            let color = display_query
                .map(|query| {
                    color_data(
                        entry.record(),
                        query,
                        entry.record().flags() & REVERSE != 0,
                        hard_clip,
                    )
                })
                .unwrap_or(ColorData {
                    symbol: None,
                    matches: false,
                    quality: None,
                });
            if kind == 1 {
                has_insertion = true;
                bases.push(Inserted {
                    base: Some(entry.record().seq_nibble(query)),
                    name_symbol,
                    color_symbol: color.symbol,
                    color_matches: color.matches,
                });
                query += 1;
            } else {
                bases.push(Inserted {
                    base: None,
                    name_symbol,
                    color_symbol: color.symbol,
                    color_matches: color.matches,
                });
            }
        }
    }
    if has_insertion { bases } else { Vec::new() }
}

fn name_symbol(record: &RawRecord, query: usize) -> Option<u8> {
    query
        .checked_add(1)
        .filter(|end| *end < record.name().len())
        .and_then(|_| record.name().get(query).copied())
}

struct ColorData {
    symbol: Option<u8>,
    matches: bool,
    quality: Option<u8>,
}

fn color_data(record: &RawRecord, query: usize, reverse: bool, hard_clip: usize) -> ColorData {
    let color_space = aux_string(record, *b"CS");
    let symbol = color_space.and_then(|values| {
        color_index(values.len(), query, reverse, hard_clip)
            .and_then(|index| values.get(index).copied())
    });
    let matches = color_space.zip(symbol).is_some_and(|(values, symbol)| {
        corrected_color(record, values, query, reverse, hard_clip) == Some(symbol)
    });
    let quality = aux_string(record, *b"CQ").and_then(|values| {
        let index = if reverse {
            values.len().checked_sub(1 + query + hard_clip)
        } else {
            Some(query)
        }?;
        values.get(index).copied()
    });
    ColorData {
        symbol,
        matches,
        quality,
    }
}

fn aux_string(record: &RawRecord, tag: [u8; 2]) -> Option<&[u8]> {
    (record.aux_type(tag) == Some(b'Z'))
        .then(|| record.aux_value(tag))
        .flatten()
        .map(|value| value.strip_suffix(&[0]).unwrap_or(value))
}

fn color_index(length: usize, query: usize, reverse: bool, hard_clip: usize) -> Option<usize> {
    if reverse {
        length.checked_sub(1 + query + hard_clip)
    } else {
        query.checked_add(1).filter(|index| *index < length)
    }
}

fn corrected_color(
    record: &RawRecord,
    colors: &[u8],
    query: usize,
    reverse: bool,
    hard_clip: usize,
) -> Option<u8> {
    let index = color_index(colors.len(), query, reverse, hard_clip)?;
    let current = sequence_symbol(record.seq_nibble(query), false);
    let previous = if reverse {
        if index == 1 {
            complement(*colors.first()?)
        } else {
            let next = query.checked_add(1)?;
            (next < record.sequence_len())
                .then(|| sequence_symbol(record.seq_nibble(next), false))?
        }
    } else if query == 0 {
        *colors.first()?
    } else {
        sequence_symbol(record.seq_nibble(query - 1), false)
    };
    Some((nt_to_int(previous)? ^ nt_to_int(current)?) + b'0')
}

fn nt_to_int(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

fn at_cigar_end(entry: ColumnEntry<'_, ReadState>, position: i64) -> bool {
    let reference_length = entry
        .cigar()
        .take(entry.projection().cigar_index + 1)
        .filter(|(kind, _)| matches!(kind, 0 | 2 | 3 | 7 | 8))
        .map(|(_, length)| i64::from(length))
        .sum::<i64>();
    i64::from(entry.record().alignment_start()) + reference_length - 1 == position
}

pub(super) struct GridBuilder {
    reference_name: String,
    reference_length: u64,
    references: Vec<String>,
    start: u64,
    width: usize,
    settings: Settings,
    reference: Option<Vec<u8>>,
    lines: Vec<Vec<Cell>>,
    caller: Caller,
    x: usize,
    offset: usize,
    alignment_rows: usize,
}

impl GridBuilder {
    pub(super) fn new(
        reference_name: String,
        reference_length: u64,
        references: Vec<String>,
        start: u64,
        width: usize,
        settings: Settings,
        reference: Option<Vec<u8>>,
    ) -> Result<Self> {
        line_count(width, 0)?;
        Ok(Self {
            reference_name,
            reference_length,
            references,
            start,
            width,
            settings,
            reference,
            lines: vec![vec![Cell::blank(); width]; 3],
            caller: Caller::new()?,
            x: 0,
            offset: 0,
            alignment_rows: 0,
        })
    }

    pub(super) fn column(&mut self, column: PositionData) -> Result<()> {
        while self.x < self.width && self.position()? < column.position {
            self.draw(None)?;
        }
        if self.x < self.width && self.position()? == column.position {
            self.ensure_rows(&column)?;
            self.draw(Some(&column))?;
        }
        Ok(())
    }

    pub(super) fn finish(mut self) -> Result<Viewport> {
        while self.x < self.width {
            self.draw(None)?;
        }
        if self.alignment_rows == 0 {
            self.lines.truncate(2);
        }
        Ok(Viewport {
            reference: self.reference_name,
            reference_length: self.reference_length,
            references: self.references,
            start: self.start,
            width: self.width,
            alignment_rows: self.alignment_rows,
            lines: self.lines,
        })
    }

    fn position(&self) -> Result<i64> {
        self.start
            .checked_sub(1)
            .and_then(|position| position.checked_add(u64::try_from(self.offset).unwrap()))
            .and_then(|position| i64::try_from(position).ok())
            .ok_or_else(|| {
                RsomicsError::InvalidInput(
                    "tview coordinate exceeds signed 64-bit range".to_owned(),
                )
            })
    }

    fn ensure_rows(&mut self, column: &PositionData) -> Result<()> {
        let required = column
            .entries
            .iter()
            .map(|entry| entry.row + 1)
            .max()
            .unwrap_or(0);
        if required > self.alignment_rows {
            let rows = line_count(self.width, required)?;
            self.lines.resize(rows, vec![Cell::blank(); self.width]);
            self.alignment_rows = required;
        }
        Ok(())
    }

    fn draw(&mut self, column: Option<&PositionData>) -> Result<()> {
        let position = self.position()?;
        let reference_base = self
            .reference
            .as_deref()
            .and_then(|sequence| sequence.get(self.offset))
            .copied()
            .unwrap_or(b'N')
            .to_ascii_uppercase();
        let interval = if position < 1_000_000_000 { 10 } else { 20 };
        if position % interval == 0 && self.width - self.x >= 10 {
            let label = (position + 1).to_string();
            for (cell, symbol) in self.lines[0][self.x..self.x + label.len()]
                .iter_mut()
                .zip(label.bytes())
            {
                cell.symbol = symbol;
            }
        }
        self.lines[1][self.x].symbol = reference_base;
        if let Some(column) = column {
            for entry in &column.entries {
                self.lines[3 + entry.row][self.x] = read_cell(entry, reference_base, self.settings);
            }
            let call = self.caller.call(&column.entries, reference_base)?;
            self.lines[2][self.x] = Cell {
                symbol: if call.base.to_ascii_uppercase() == reference_base {
                    b'.'
                } else {
                    call.base
                },
                color: quality_color(call.confidence),
                underline: true,
            };
            if !self.settings.hide_insertions {
                for insertion_offset in 0..column.insertion_width {
                    if self.x + 1 >= self.width {
                        break;
                    }
                    self.x += 1;
                    self.lines[1][self.x] = Cell {
                        symbol: b'*',
                        color: CellColor::Red,
                        underline: false,
                    };
                    for entry in &column.entries {
                        let symbol = entry
                            .insertion
                            .get(insertion_offset)
                            .map_or(b'*', |inserted| {
                                insertion_symbol(inserted, entry, self.settings)
                            });
                        self.lines[3 + entry.row][self.x] = Cell {
                            symbol,
                            color: entry_color(entry, self.settings),
                            underline: underlined(entry),
                        };
                    }
                }
            }
        }
        self.x += 1;
        self.offset += 1;
        Ok(())
    }
}

fn read_cell(entry: &Entry, reference: u8, settings: Settings) -> Cell {
    Cell {
        symbol: read_symbol(entry, reference, settings),
        color: entry_color(entry, settings),
        underline: underlined(entry),
    }
}

fn read_symbol(entry: &Entry, reference: u8, settings: Settings) -> u8 {
    if entry.reference_skip {
        if settings.skips_as_deletions {
            return b'*';
        }
        return if entry.reverse { b'<' } else { b'>' };
    }
    if entry.deletion {
        return b'*';
    }
    if settings.base == super::BaseMode::ColorSpace
        && let Some(symbol) = entry.color_symbol
    {
        if settings.dots && entry.color_matches {
            return if entry.reverse { b',' } else { b'.' };
        }
        return strand_symbol(symbol, entry.reverse);
    }
    if settings.show_names {
        return entry.name_symbol.unwrap_or(b' ');
    }
    let base = sequence_symbol(entry.base, entry.reverse);
    if settings.dots && base.to_ascii_uppercase() == reference {
        if entry.reverse { b',' } else { b'.' }
    } else {
        base
    }
}

fn insertion_symbol(inserted: &Inserted, entry: &Entry, settings: Settings) -> u8 {
    if settings.base == super::BaseMode::ColorSpace
        && let Some(symbol) = inserted.color_symbol
    {
        if settings.dots && inserted.color_matches {
            return if entry.reverse { b',' } else { b'.' };
        }
        return strand_symbol(symbol, entry.reverse);
    }
    if settings.show_names {
        return inserted.name_symbol.unwrap_or(b' ');
    }
    inserted
        .base
        .map(|base| sequence_symbol(base, entry.reverse))
        .unwrap_or(b'*')
}

fn entry_color(entry: &Entry, settings: Settings) -> CellColor {
    match settings.color {
        super::ColorMode::MappingQuality => quality_color(entry.mapping_quality / 10 + 1),
        super::ColorMode::BaseQuality => quality_color(entry.quality / 10 + 1),
        super::ColorMode::Nucleotide => nucleotide_color(entry.base),
        super::ColorMode::ColorSpace => entry
            .color_symbol
            .and_then(color_space_color)
            .unwrap_or_else(|| nucleotide_color(entry.base)),
        super::ColorMode::ColorQuality => {
            quality_color(entry.color_quality.unwrap_or(entry.quality) / 10 + 1)
        }
    }
}

fn quality_color(bucket: u8) -> CellColor {
    match bucket.min(4) {
        0 | 1 => CellColor::Blue,
        2 => CellColor::Green,
        3 => CellColor::Yellow,
        _ => CellColor::White,
    }
}

fn nucleotide_color(base: u8) -> CellColor {
    match base {
        1 => CellColor::Green,
        2 => CellColor::Cyan,
        4 => CellColor::Magenta,
        8 => CellColor::Red,
        _ => CellColor::Blue,
    }
}

fn underlined(entry: &Entry) -> bool {
    entry.flags & SECONDARY != 0 || (entry.flags & PAIRED != 0 && entry.flags & PROPER_PAIR == 0)
}

fn color_space_color(symbol: u8) -> Option<CellColor> {
    match symbol {
        b'0' => Some(CellColor::Green),
        b'1' => Some(CellColor::Cyan),
        b'2' => Some(CellColor::Magenta),
        b'3' => Some(CellColor::Red),
        b'4' => Some(CellColor::Blue),
        _ => None,
    }
}

fn sequence_symbol(base: u8, reverse: bool) -> u8 {
    let symbol = b"=ACMGRSVTWYHKDBN"
        .get(usize::from(base))
        .copied()
        .unwrap_or(b'N');
    strand_symbol(symbol, reverse)
}

fn strand_symbol(symbol: u8, reverse: bool) -> u8 {
    if reverse {
        symbol.to_ascii_lowercase()
    } else {
        symbol.to_ascii_uppercase()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tview::{BaseMode, ColorMode};

    fn entry() -> Entry {
        Entry {
            row: 0,
            flags: 0,
            name_symbol: Some(b'Q'),
            color_symbol: Some(b'2'),
            color_matches: true,
            color_quality: Some(9),
            base: 1,
            quality: 19,
            mapping_quality: 29,
            reverse: false,
            deletion: false,
            reference_skip: false,
            insertion: Vec::new(),
        }
    }

    #[test]
    fn symbols_cover_strand_names_skips_deletions_and_color_space() {
        let mut entry = entry();
        let mut settings = Settings::default();
        assert_eq!(read_symbol(&entry, b'A', settings), b'.');
        assert_eq!(read_symbol(&entry, b'C', settings), b'A');
        entry.reverse = true;
        assert_eq!(read_symbol(&entry, b'A', settings), b',');
        assert_eq!(read_symbol(&entry, b'C', settings), b'a');

        settings.show_names = true;
        assert_eq!(read_symbol(&entry, b'A', settings), b'Q');
        settings.show_names = false;
        entry.reference_skip = true;
        assert_eq!(read_symbol(&entry, b'A', settings), b'<');
        settings.skips_as_deletions = true;
        assert_eq!(read_symbol(&entry, b'A', settings), b'*');
        entry.reference_skip = false;
        entry.deletion = true;
        assert_eq!(read_symbol(&entry, b'A', settings), b'*');

        entry.deletion = false;
        entry.reverse = false;
        settings.base = BaseMode::ColorSpace;
        settings.skips_as_deletions = false;
        assert_eq!(read_symbol(&entry, b'A', settings), b'.');
        settings.dots = false;
        assert_eq!(read_symbol(&entry, b'A', settings), b'2');
    }

    #[test]
    fn insertion_symbols_and_style_modes_are_explicit() {
        let mut entry = entry();
        let inserted = Inserted {
            base: Some(4),
            name_symbol: Some(b'R'),
            color_symbol: Some(b'3'),
            color_matches: false,
        };
        let pad = Inserted {
            base: None,
            name_symbol: Some(b'P'),
            color_symbol: Some(b'1'),
            color_matches: false,
        };
        let mut settings = Settings::default();
        assert_eq!(insertion_symbol(&inserted, &entry, settings), b'G');
        assert_eq!(insertion_symbol(&pad, &entry, settings), b'*');
        settings.show_names = true;
        assert_eq!(insertion_symbol(&pad, &entry, settings), b'P');
        settings.show_names = false;
        settings.base = BaseMode::ColorSpace;
        assert_eq!(insertion_symbol(&pad, &entry, settings), b'1');

        settings.color = ColorMode::MappingQuality;
        assert_eq!(entry_color(&entry, settings), CellColor::Yellow);
        settings.color = ColorMode::BaseQuality;
        assert_eq!(entry_color(&entry, settings), CellColor::Green);
        settings.color = ColorMode::Nucleotide;
        assert_eq!(entry_color(&entry, settings), CellColor::Green);
        settings.color = ColorMode::ColorSpace;
        assert_eq!(entry_color(&entry, settings), CellColor::Magenta);
        settings.color = ColorMode::ColorQuality;
        assert_eq!(entry_color(&entry, settings), CellColor::Blue);

        entry.flags = SECONDARY;
        assert!(underlined(&entry));
        entry.flags = PAIRED;
        assert!(underlined(&entry));
        entry.flags = PAIRED | PROPER_PAIR;
        assert!(!underlined(&entry));
    }
}
