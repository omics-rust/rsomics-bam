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
        for cell in line {
            let class = class(*cell);
            if class.is_empty() {
                write_byte(&mut output, cell.symbol)?;
            } else {
                write!(output, "<span class=\"{class}\">").map_err(RsomicsError::Io)?;
                write_byte(&mut output, cell.symbol)?;
                output.write_all(b"</span>").map_err(RsomicsError::Io)?;
            }
        }
        output.write_all(b"\n").map_err(RsomicsError::Io)?;
    }
    output
        .write_all(b"</pre></body></html>\n")
        .map_err(RsomicsError::Io)?;
    output.flush().map_err(RsomicsError::Io)
}

fn class(cell: super::Cell) -> String {
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
        (true, false) => String::new(),
        (true, true) => "underline".to_owned(),
        (false, false) => color.to_owned(),
        (false, true) => format!("{color} underline"),
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
