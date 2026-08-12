use std::io::Write;

use rsomics_common::{Result, RsomicsError};

use super::Viewport;

pub(super) fn text(viewport: &Viewport, mut output: impl Write) -> Result<()> {
    let mut symbols = Vec::with_capacity(viewport.width + 1);
    for line in &viewport.lines {
        symbols.clear();
        symbols.extend(line.iter().map(|cell| cell.symbol));
        symbols.push(b'\n');
        output.write_all(&symbols).map_err(RsomicsError::Io)?;
    }
    output.flush().map_err(RsomicsError::Io)
}

pub(super) fn html(viewport: &Viewport, mut output: impl Write) -> Result<()> {
    write!(
        output,
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><title>{}:{}</title><style>body{{background:#000;color:#fff}}.blue{{color:#44f}}.green{{color:#0c0}}.yellow{{color:#cc0}}.white{{color:#fff}}.cyan{{color:#0cc}}.magenta{{color:#c0c}}.red{{color:#f44}}.underline{{text-decoration:underline}}.location{{position:absolute;left:-10000px}}</style></head><body><pre data-location=\"{}:{}\">",
        escape(&viewport.reference),
        viewport.start,
        escape(&viewport.reference),
        viewport.start
    )
    .map_err(RsomicsError::Io)?;
    writeln!(
        output,
        "<span class=\"location\">{}:{}</span>",
        escape(&viewport.reference),
        viewport.start
    )
    .map_err(RsomicsError::Io)?;
    for line in &viewport.lines {
        let mut active = "";
        for cell in line {
            let next = class(*cell);
            if next != active {
                if !active.is_empty() {
                    output.write_all(b"</span>").map_err(RsomicsError::Io)?;
                }
                if !next.is_empty() {
                    write!(output, "<span class=\"{next}\">").map_err(RsomicsError::Io)?;
                }
                active = next;
            }
            write_byte(&mut output, cell.symbol)?;
        }
        if !active.is_empty() {
            output.write_all(b"</span>").map_err(RsomicsError::Io)?;
        }
        output.write_all(b"\n").map_err(RsomicsError::Io)?;
    }
    output
        .write_all(b"</pre></body></html>\n")
        .map_err(RsomicsError::Io)?;
    output.flush().map_err(RsomicsError::Io)
}

fn class(cell: super::Cell) -> &'static str {
    let color = match cell.color {
        super::CellColor::Default => "",
        super::CellColor::Blue => "blue",
        super::CellColor::Green => "green",
        super::CellColor::Yellow => "yellow",
        super::CellColor::White => "white",
        super::CellColor::Cyan => "cyan",
        super::CellColor::Magenta => "magenta",
        super::CellColor::Red => "red",
    };
    match (color.is_empty(), cell.underline) {
        (true, false) => "",
        (true, true) => "underline",
        (false, false) => color,
        (false, true) => match cell.color {
            super::CellColor::Default => unreachable!(),
            super::CellColor::Blue => "blue underline",
            super::CellColor::Green => "green underline",
            super::CellColor::Yellow => "yellow underline",
            super::CellColor::White => "white underline",
            super::CellColor::Cyan => "cyan underline",
            super::CellColor::Magenta => "magenta underline",
            super::CellColor::Red => "red underline",
        },
    }
}

fn write_byte(output: &mut impl Write, byte: u8) -> Result<()> {
    match byte {
        b'&' => output.write_all(b"&amp;"),
        b'<' => output.write_all(b"&lt;"),
        b'>' => output.write_all(b"&gt;"),
        b'\"' => output.write_all(b"&quot;"),
        byte if byte.is_ascii() => output.write_all(&[byte]),
        _ => output.write_all("�".as_bytes()),
    }
    .map_err(RsomicsError::Io)
}

fn escape(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '&' => escaped.push_str("&amp;"),
            '<' => escaped.push_str("&lt;"),
            '>' => escaped.push_str("&gt;"),
            '"' => escaped.push_str("&quot;"),
            _ => escaped.push(character),
        }
    }
    escaped
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tview::{Cell, CellColor};

    #[test]
    fn html_groups_adjacent_cells_with_the_same_style() {
        let viewport = Viewport {
            reference: "chr1".to_owned(),
            reference_length: 10,
            references: vec!["chr1".to_owned()],
            start: 1,
            width: 5,
            alignment_rows: 0,
            lines: vec![vec![
                Cell {
                    symbol: b'A',
                    color: CellColor::Blue,
                    underline: false,
                },
                Cell {
                    symbol: b'C',
                    color: CellColor::Blue,
                    underline: false,
                },
                Cell {
                    symbol: b'G',
                    color: CellColor::Default,
                    underline: false,
                },
                Cell {
                    symbol: b'T',
                    color: CellColor::Default,
                    underline: true,
                },
                Cell {
                    symbol: b'&',
                    color: CellColor::Default,
                    underline: true,
                },
            ]],
        };
        let mut output = Vec::new();
        html(&viewport, &mut output).unwrap();
        let output = String::from_utf8(output).unwrap();
        assert!(
            output
                .contains("<span class=\"blue\">AC</span>G<span class=\"underline\">T&amp;</span>")
        );
    }
}
