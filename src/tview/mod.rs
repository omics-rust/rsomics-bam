mod call;
mod grid;
mod model;
mod render;
mod rows;
mod terminal;

use std::io::Write;
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

const MAX_WIDTH: usize = 1_000_000;
const MAX_CELLS: usize = 16 * 1024 * 1024;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum ColorMode {
    #[default]
    MappingQuality,
    BaseQuality,
    Nucleotide,
    ColorSpace,
    ColorQuality,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum BaseMode {
    #[default]
    Nucleotide,
    ColorSpace,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct Settings {
    pub(crate) color: ColorMode,
    pub(crate) base: BaseMode,
    pub(crate) dots: bool,
    pub(crate) hide_insertions: bool,
    pub(crate) skips_as_deletions: bool,
    pub(crate) show_names: bool,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            color: ColorMode::MappingQuality,
            base: BaseMode::Nucleotide,
            dots: true,
            hide_insertions: false,
            skips_as_deletions: false,
            show_names: false,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum CellColor {
    #[default]
    Default,
    Blue,
    Green,
    Yellow,
    White,
    Cyan,
    Magenta,
    Red,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct Cell {
    pub(crate) symbol: u8,
    pub(crate) color: CellColor,
    pub(crate) underline: bool,
}

impl Cell {
    fn blank() -> Self {
        Self {
            symbol: b' ',
            ..Self::default()
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Format {
    #[default]
    Text,
    Html,
}

#[derive(Clone, Copy, Debug)]
pub struct Options<'a> {
    pub reference: Option<&'a Path>,
    pub index: Option<&'a Path>,
    pub position: Option<&'a str>,
    pub sample: Option<&'a str>,
    pub width: usize,
    pub hide_insertions: bool,
    pub additional_threads: usize,
}

impl Default for Options<'_> {
    fn default() -> Self {
        Self {
            reference: None,
            index: None,
            position: None,
            sample: None,
            width: 80,
            hide_insertions: false,
            additional_threads: 0,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub reference: String,
    pub start: u64,
    pub width: usize,
    pub alignment_rows: usize,
}

pub(crate) struct Viewport {
    reference: String,
    reference_length: u64,
    references: Vec<String>,
    start: u64,
    width: usize,
    alignment_rows: usize,
    lines: Vec<Vec<Cell>>,
}

pub fn write(
    input: &Path,
    options: Options<'_>,
    format: Format,
    output: impl Write,
) -> Result<Summary> {
    validate_width(options.width)?;
    let viewport = model::load(
        input,
        options,
        Settings {
            hide_insertions: options.hide_insertions,
            ..Settings::default()
        },
    )?;
    match format {
        Format::Text => render::text(&viewport, output)?,
        Format::Html => render::html(&viewport, output)?,
    }
    Ok(Summary {
        reference: viewport.reference,
        start: viewport.start,
        width: viewport.width,
        alignment_rows: viewport.alignment_rows,
    })
}

pub(crate) fn interactive(input: &Path, options: Options<'_>) -> Result<Summary> {
    terminal::run(input, options)
}

fn validate_width(width: usize) -> Result<()> {
    if width == 0 || width > MAX_WIDTH {
        return Err(RsomicsError::ConfigError(format!(
            "tview width must be between 1 and {MAX_WIDTH}"
        )));
    }
    Ok(())
}

fn line_count(width: usize, alignment_rows: usize) -> Result<usize> {
    let header_rows = if alignment_rows == 0 { 2 } else { 3 };
    let rows = alignment_rows.checked_add(header_rows).ok_or_else(|| {
        RsomicsError::ConfigError("tview row count exceeds this platform".to_owned())
    })?;
    let cells = rows.checked_mul(width).ok_or_else(|| {
        RsomicsError::ConfigError("tview viewport size exceeds this platform".to_owned())
    })?;
    if cells > MAX_CELLS {
        return Err(RsomicsError::ConfigError(format!(
            "tview viewport requires {cells} cells; the limit is {MAX_CELLS}"
        )));
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn viewport_limits_fail_before_allocation() {
        assert!(validate_width(MAX_WIDTH).is_ok());
        assert!(validate_width(MAX_WIDTH + 1).is_err());
        assert_eq!(line_count(80, 0).unwrap(), 2);
        assert_eq!(line_count(80, 10).unwrap(), 13);
        assert!(line_count(MAX_WIDTH, 20).is_err());
        assert!(line_count(usize::MAX, usize::MAX).is_err());
    }
}
