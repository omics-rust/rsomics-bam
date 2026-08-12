use std::io::{self, Write};
use std::path::{Path, PathBuf};

use crossterm::cursor::{Hide, MoveTo, Show};
use crossterm::event::{self, Event, KeyCode, KeyEvent, KeyEventKind, KeyModifiers};
use crossterm::style::{
    Attribute, Color, Print, ResetColor, SetAttribute, SetBackgroundColor, SetForegroundColor,
};
use crossterm::terminal::{
    self, Clear, ClearType, EnterAlternateScreen, LeaveAlternateScreen, disable_raw_mode,
    enable_raw_mode,
};
use crossterm::{execute, queue};
use rsomics_common::{Result, RsomicsError};

use super::{
    BaseMode, Cell, CellColor, ColorMode, Options, Settings, Summary, Viewport, model,
    validate_width,
};

struct Session {
    input: PathBuf,
    reference: Option<PathBuf>,
    index: Option<PathBuf>,
    sample: Option<String>,
    threads: usize,
}

#[derive(Clone, Debug)]
struct State {
    reference: String,
    start: u64,
    row_shift: usize,
    settings: Settings,
    inverse: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Action {
    None,
    Redraw,
    Help,
    Goto,
    Quit,
}

impl State {
    fn key(&mut self, key: KeyEvent, viewport: &Viewport, visible_rows: usize) -> Action {
        if key.kind == KeyEventKind::Release {
            return Action::None;
        }
        if key.modifiers.contains(KeyModifiers::CONTROL) {
            return match key.code {
                KeyCode::Char('c') => Action::Quit,
                KeyCode::Char('h') => {
                    self.move_left(1000);
                    Action::Redraw
                }
                KeyCode::Char('l') => {
                    self.move_right(1000, viewport.reference_length);
                    Action::Redraw
                }
                _ => Action::None,
            };
        }
        match key.code {
            KeyCode::Char('q') | KeyCode::Esc => Action::Quit,
            KeyCode::Char('?') => Action::Help,
            KeyCode::Char('g' | '/') => Action::Goto,
            KeyCode::Left | KeyCode::Char('h') => {
                self.move_left(1);
                Action::Redraw
            }
            KeyCode::Right | KeyCode::Char('l') => {
                self.move_right(1, viewport.reference_length);
                Action::Redraw
            }
            KeyCode::Char('H') => {
                self.move_left(20);
                Action::Redraw
            }
            KeyCode::Char('L') => {
                self.move_right(20, viewport.reference_length);
                Action::Redraw
            }
            KeyCode::Char(' ') => {
                self.move_right(viewport.width as u64, viewport.reference_length);
                Action::Redraw
            }
            KeyCode::Backspace => {
                self.move_left(viewport.width as u64);
                Action::Redraw
            }
            KeyCode::Up | KeyCode::Char('j') => {
                self.row_shift = self.row_shift.saturating_sub(1);
                Action::Redraw
            }
            KeyCode::Down | KeyCode::Char('k') => {
                self.row_shift = (self.row_shift + 1).min(max_row_shift(viewport, visible_rows));
                Action::Redraw
            }
            KeyCode::Char('J') => {
                self.row_shift = self.row_shift.saturating_sub(20);
                Action::Redraw
            }
            KeyCode::Char('K') => {
                self.row_shift = (self.row_shift + 20).min(max_row_shift(viewport, visible_rows));
                Action::Redraw
            }
            KeyCode::Char('m') => {
                self.settings.color = ColorMode::MappingQuality;
                Action::Redraw
            }
            KeyCode::Char('b') => {
                self.settings.color = ColorMode::BaseQuality;
                Action::Redraw
            }
            KeyCode::Char('n') => {
                self.settings.color = ColorMode::Nucleotide;
                Action::Redraw
            }
            KeyCode::Char('c') => {
                self.settings.color = ColorMode::ColorSpace;
                Action::Redraw
            }
            KeyCode::Char('z') => {
                self.settings.color = ColorMode::ColorQuality;
                Action::Redraw
            }
            KeyCode::Char('.') => {
                self.settings.dots = !self.settings.dots;
                Action::Redraw
            }
            KeyCode::Char('s') => {
                self.settings.skips_as_deletions = !self.settings.skips_as_deletions;
                Action::Redraw
            }
            KeyCode::Char('r') => {
                self.settings.show_names = !self.settings.show_names;
                Action::Redraw
            }
            KeyCode::Char('N') => {
                self.settings.base = BaseMode::Nucleotide;
                Action::Redraw
            }
            KeyCode::Char('C') => {
                self.settings.base = BaseMode::ColorSpace;
                Action::Redraw
            }
            KeyCode::Char('i') => {
                self.settings.hide_insertions = !self.settings.hide_insertions;
                Action::Redraw
            }
            KeyCode::Char('v') => {
                self.inverse = !self.inverse;
                Action::Redraw
            }
            _ => Action::None,
        }
    }

    fn move_left(&mut self, distance: u64) {
        self.start = self.start.saturating_sub(distance).max(1);
    }

    fn move_right(&mut self, distance: u64, reference_length: u64) {
        self.start = self.start.saturating_add(distance).min(reference_length);
    }
}

pub(super) fn run(input: &Path, options: Options<'_>) -> Result<Summary> {
    let session = Session {
        input: input.to_owned(),
        reference: options.reference.map(Path::to_owned),
        index: options.index.map(Path::to_owned),
        sample: options.sample.map(str::to_owned),
        threads: options.additional_threads,
    };
    let (width, height) = dimensions()?;
    let settings = Settings {
        hide_insertions: options.hide_insertions,
        ..Settings::default()
    };
    let mut viewport = load(&session, options.position, width, settings)?;
    let mut state = State {
        reference: viewport.reference.clone(),
        start: viewport.start,
        row_shift: 0,
        settings,
        inverse: false,
    };
    let mut screen = Screen::enter()?;
    let result = (|| {
        let mut dimensions = (width, height);
        loop {
            state.row_shift = state
                .row_shift
                .min(max_row_shift(&viewport, visible_rows(dimensions.1)));
            screen.draw(&viewport, &state, dimensions.1)?;
            match event::read().map_err(RsomicsError::Io)? {
                Event::Resize(width, height) => {
                    dimensions = checked_dimensions(width, height)?;
                    viewport = load(
                        &session,
                        Some(&format!("{}:{}", state.reference, state.start)),
                        dimensions.0,
                        state.settings,
                    )?;
                }
                Event::Key(key) => match state.key(key, &viewport, visible_rows(dimensions.1)) {
                    Action::None => {}
                    Action::Quit => break,
                    Action::Help => screen.help(dimensions)?,
                    Action::Goto => {
                        if let Some(value) = screen.goto(&state, &viewport, dimensions)? {
                            match load(&session, Some(&value), dimensions.0, state.settings) {
                                Ok(next) => {
                                    state.reference = next.reference.clone();
                                    state.start = next.start;
                                    state.row_shift = 0;
                                    viewport = next;
                                }
                                Err(error) => screen.message(&error.to_string(), dimensions)?,
                            }
                        }
                    }
                    Action::Redraw => {
                        viewport = load(
                            &session,
                            Some(&format!("{}:{}", state.reference, state.start)),
                            dimensions.0,
                            state.settings,
                        )?;
                    }
                },
                _ => {}
            }
        }
        Ok(Summary {
            reference: viewport.reference,
            start: viewport.start,
            width: viewport.width,
            alignment_rows: viewport.alignment_rows,
        })
    })();
    let restore = screen.restore();
    match (result, restore) {
        (Err(error), _) => Err(error),
        (Ok(_), Err(error)) => Err(error),
        (Ok(summary), Ok(())) => Ok(summary),
    }
}

fn load(
    session: &Session,
    position: Option<&str>,
    width: usize,
    settings: Settings,
) -> Result<Viewport> {
    validate_width(width)?;
    model::load(
        &session.input,
        Options {
            reference: session.reference.as_deref(),
            index: session.index.as_deref(),
            position,
            sample: session.sample.as_deref(),
            width,
            hide_insertions: settings.hide_insertions,
            additional_threads: session.threads,
        },
        settings,
    )
}

fn dimensions() -> Result<(usize, usize)> {
    let (width, height) = terminal::size().map_err(RsomicsError::Io)?;
    checked_dimensions(width, height)
}

fn checked_dimensions(width: u16, height: u16) -> Result<(usize, usize)> {
    if width == 0 || height < 4 {
        return Err(RsomicsError::ConfigError(
            "terminal tview requires at least 1 column and 4 rows".to_owned(),
        ));
    }
    Ok((usize::from(width), usize::from(height)))
}

fn visible_rows(height: usize) -> usize {
    height.saturating_sub(3)
}

fn max_row_shift(viewport: &Viewport, visible_rows: usize) -> usize {
    viewport.alignment_rows.saturating_sub(visible_rows)
}

struct Screen {
    output: io::Stdout,
    active: bool,
}

impl Screen {
    fn enter() -> Result<Self> {
        enable_raw_mode().map_err(RsomicsError::Io)?;
        let mut screen = Self {
            output: io::stdout(),
            active: true,
        };
        if let Err(error) = execute!(screen.output, EnterAlternateScreen, Hide) {
            let _ = screen.restore();
            return Err(RsomicsError::Io(error));
        }
        Ok(screen)
    }

    fn draw(&mut self, viewport: &Viewport, state: &State, height: usize) -> Result<()> {
        queue!(self.output, MoveTo(0, 0), Clear(ClearType::All)).map_err(RsomicsError::Io)?;
        for (screen_row, line) in viewport.lines.iter().take(3).enumerate() {
            self.line(screen_row, line, state.inverse)?;
        }
        for (screen_row, line) in viewport
            .lines
            .iter()
            .skip(3 + state.row_shift)
            .take(visible_rows(height))
            .enumerate()
        {
            self.line(3 + screen_row, line, state.inverse)?;
        }
        queue!(
            self.output,
            ResetColor,
            SetAttribute(Attribute::NoUnderline)
        )
        .map_err(RsomicsError::Io)?;
        self.output.flush().map_err(RsomicsError::Io)
    }

    fn line(&mut self, row: usize, line: &[Cell], inverse: bool) -> Result<()> {
        queue!(self.output, MoveTo(0, u16::try_from(row).unwrap())).map_err(RsomicsError::Io)?;
        for cell in line {
            self.style(*cell, inverse)?;
            queue!(self.output, Print(char::from(cell.symbol))).map_err(RsomicsError::Io)?;
        }
        Ok(())
    }

    fn style(&mut self, cell: Cell, inverse: bool) -> Result<()> {
        queue!(self.output, ResetColor).map_err(RsomicsError::Io)?;
        if cell.color != CellColor::Default {
            let color = terminal_color(cell.color);
            if inverse {
                let foreground = match cell.color {
                    CellColor::Blue | CellColor::Magenta | CellColor::Red => Color::White,
                    _ => Color::Black,
                };
                queue!(
                    self.output,
                    SetForegroundColor(foreground),
                    SetBackgroundColor(color)
                )
                .map_err(RsomicsError::Io)?;
            } else {
                queue!(self.output, SetForegroundColor(color)).map_err(RsomicsError::Io)?;
            }
        }
        queue!(
            self.output,
            SetAttribute(if cell.underline {
                Attribute::Underlined
            } else {
                Attribute::NoUnderline
            })
        )
        .map_err(RsomicsError::Io)
    }

    fn help(&mut self, dimensions: (usize, usize)) -> Result<()> {
        const HELP: &[&str] = &[
            "rsomics-bam tview",
            "? help   q/Esc/Ctrl-C exit   g or / goto",
            "arrows h/j/k/l move   H/J/K/L move 20",
            "Ctrl-H/Ctrl-L move 1 kb   space/backspace move a viewport",
            "m mapping quality   b base quality   n nucleotide color",
            "c color-space color   z color-space quality   v inverse",
            ". dots   s reference skips   r read names   i insertions",
            "N nucleotide bases   C color-space bases",
            "underline: secondary or improperly paired",
        ];
        queue!(self.output, MoveTo(0, 0), Clear(ClearType::All), ResetColor)
            .map_err(RsomicsError::Io)?;
        for (row, line) in HELP.iter().take(dimensions.1).enumerate() {
            let end = line.len().min(dimensions.0);
            queue!(
                self.output,
                MoveTo(0, u16::try_from(row).unwrap()),
                Print(&line[..end])
            )
            .map_err(RsomicsError::Io)?;
        }
        self.output.flush().map_err(RsomicsError::Io)?;
        event::read().map(drop).map_err(RsomicsError::Io)
    }

    fn goto(
        &mut self,
        state: &State,
        viewport: &Viewport,
        dimensions: (usize, usize),
    ) -> Result<Option<String>> {
        let mut value = String::new();
        loop {
            let prompt = format!("goto> {value}");
            queue!(
                self.output,
                MoveTo(0, u16::try_from(dimensions.1 - 1).unwrap()),
                Clear(ClearType::CurrentLine),
                ResetColor,
                Print(prefix(&prompt, dimensions.0))
            )
            .map_err(RsomicsError::Io)?;
            self.output.flush().map_err(RsomicsError::Io)?;
            let Event::Key(key) = event::read().map_err(RsomicsError::Io)? else {
                continue;
            };
            if key.kind == KeyEventKind::Release {
                continue;
            }
            match key.code {
                KeyCode::Esc => return Ok(None),
                KeyCode::Enter => {
                    if value.is_empty() {
                        return Ok(None);
                    }
                    return Ok(Some(normalize_goto(&value, state)));
                }
                KeyCode::Backspace => {
                    value.pop();
                }
                KeyCode::Tab => complete(&mut value, &viewport.references),
                KeyCode::Char(character)
                    if !key.modifiers.contains(KeyModifiers::CONTROL) && value.len() < 4096 =>
                {
                    value.push(character);
                }
                _ => {}
            }
        }
    }

    fn message(&mut self, message: &str, dimensions: (usize, usize)) -> Result<()> {
        queue!(
            self.output,
            MoveTo(0, u16::try_from(dimensions.1 - 1).unwrap()),
            Clear(ClearType::CurrentLine),
            ResetColor,
            Print(prefix(message, dimensions.0))
        )
        .map_err(RsomicsError::Io)?;
        self.output.flush().map_err(RsomicsError::Io)?;
        event::read().map(drop).map_err(RsomicsError::Io)
    }

    fn restore(&mut self) -> Result<()> {
        if !self.active {
            return Ok(());
        }
        let screen = execute!(
            self.output,
            ResetColor,
            SetAttribute(Attribute::Reset),
            Show,
            LeaveAlternateScreen
        )
        .map_err(RsomicsError::Io);
        let raw = disable_raw_mode().map_err(RsomicsError::Io);
        self.active = false;
        screen.and(raw)
    }
}

fn prefix(value: &str, characters: usize) -> &str {
    value
        .char_indices()
        .nth(characters)
        .map_or(value, |(end, _)| &value[..end])
}

impl Drop for Screen {
    fn drop(&mut self) {
        let _ = self.restore();
    }
}

fn normalize_goto(value: &str, state: &State) -> String {
    value.strip_prefix('=').map_or_else(
        || value.to_owned(),
        |position| format!("{}:{position}", state.reference),
    )
}

fn complete(value: &mut String, references: &[String]) {
    let prefix = value
        .split_once(':')
        .map_or(value.as_str(), |(name, _)| name);
    if let Some(reference) = references.iter().find(|name| name.starts_with(prefix)) {
        *value = format!("{reference}:");
    }
}

fn terminal_color(color: CellColor) -> Color {
    match color {
        CellColor::Default | CellColor::White => Color::White,
        CellColor::Blue => Color::Blue,
        CellColor::Green => Color::Green,
        CellColor::Yellow => Color::Yellow,
        CellColor::Cyan => Color::Cyan,
        CellColor::Magenta => Color::Magenta,
        CellColor::Red => Color::Red,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn viewport() -> Viewport {
        Viewport {
            reference: "chr1".to_owned(),
            reference_length: 10_000,
            references: vec!["chr1".to_owned(), "chr2".to_owned()],
            start: 100,
            width: 80,
            alignment_rows: 30,
            lines: vec![vec![Cell::blank(); 80]; 33],
        }
    }

    fn state() -> State {
        State {
            reference: "chr1".to_owned(),
            start: 100,
            row_shift: 0,
            settings: Settings::default(),
            inverse: false,
        }
    }

    fn key(code: KeyCode) -> KeyEvent {
        KeyEvent::new(code, KeyModifiers::NONE)
    }

    #[test]
    fn movement_saturates_and_scrolls_by_the_viewport() {
        let viewport = viewport();
        let mut state = state();
        assert_eq!(
            state.key(key(KeyCode::Char('h')), &viewport, 20),
            Action::Redraw
        );
        assert_eq!(state.start, 99);
        state.start = 1;
        state.key(key(KeyCode::Char('H')), &viewport, 20);
        assert_eq!(state.start, 1);
        state.key(key(KeyCode::Char(' ')), &viewport, 20);
        assert_eq!(state.start, 81);
        state.key(key(KeyCode::Char('K')), &viewport, 20);
        assert_eq!(state.row_shift, 10);
    }

    #[test]
    fn every_display_mode_and_toggle_has_a_deterministic_key() {
        let viewport = viewport();
        let mut state = state();
        for (code, mode) in [
            ('m', ColorMode::MappingQuality),
            ('b', ColorMode::BaseQuality),
            ('n', ColorMode::Nucleotide),
            ('c', ColorMode::ColorSpace),
            ('z', ColorMode::ColorQuality),
        ] {
            state.key(key(KeyCode::Char(code)), &viewport, 20);
            assert_eq!(state.settings.color, mode);
        }
        for code in ['.', 's', 'r', 'i', 'v'] {
            assert_eq!(
                state.key(key(KeyCode::Char(code)), &viewport, 20),
                Action::Redraw
            );
        }
        state.key(key(KeyCode::Char('C')), &viewport, 20);
        assert_eq!(state.settings.base, BaseMode::ColorSpace);
        state.key(key(KeyCode::Char('N')), &viewport, 20);
        assert_eq!(state.settings.base, BaseMode::Nucleotide);
    }

    #[test]
    fn goto_normalization_and_completion_are_bounded() {
        let state = state();
        assert_eq!(normalize_goto("=42", &state), "chr1:42");
        let mut value = "chr2".to_owned();
        complete(&mut value, &viewport().references);
        assert_eq!(value, "chr2:");
        assert_eq!(prefix("αβγ", 2), "αβ");
    }
}
