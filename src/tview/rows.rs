use rsomics_common::{Result, RsomicsError};
use rsomics_pileup::Column;

#[derive(Clone, Debug)]
pub(super) struct ReadState {
    pub(super) start: i64,
    pub(super) end: i64,
}

#[derive(Clone, Copy, Debug)]
struct FreeRow {
    row: usize,
    delay: u8,
}

#[derive(Default)]
pub(super) struct RowPacker {
    free: Vec<FreeRow>,
    previous: Vec<usize>,
    recycled: usize,
    tail_delay: u8,
    level_count: usize,
}

impl RowPacker {
    pub(super) fn pack(&mut self, column: &Column<'_, ReadState>) -> Result<Vec<usize>> {
        for row in &mut self.free {
            row.delay = row.delay.saturating_sub(1);
        }
        let position = column.position();
        let mut previous_index = 0;
        let mut rows = Vec::with_capacity(column.entries().len());
        let mut next_previous = Vec::with_capacity(self.previous.len());
        for entry in column.entries() {
            let state = entry.state();
            let row = if position == state.start {
                self.allocate()
            } else {
                let row = self.previous.get(previous_index).copied().ok_or_else(|| {
                    RsomicsError::InvalidInput(
                        "tview row state does not match the alignment pileup".to_owned(),
                    )
                })?;
                previous_index += 1;
                row
            };
            rows.push(row);
            if position == state.end {
                self.release(row);
            } else {
                next_previous.push(row);
            }
        }
        if previous_index != self.previous.len() {
            return Err(RsomicsError::InvalidInput(
                "tview row state retained an absent alignment".to_owned(),
            ));
        }
        let current_levels = rows.iter().map(|row| row + 1).max().unwrap_or(0);
        let before = self.free.len();
        self.free.retain(|row| row.row < current_levels);
        self.recycled += before - self.free.len();
        self.free.sort_by_key(|row| (row.delay, row.row));
        self.previous = next_previous;
        self.level_count = current_levels;
        Ok(rows)
    }

    fn allocate(&mut self) -> usize {
        if self.free.first().is_some_and(|row| row.delay == 0) {
            self.recycled += 1;
            return self.free.remove(0).row;
        }
        let row = self.level_count;
        self.level_count += 1;
        row
    }

    fn release(&mut self, row: usize) {
        self.free.push(FreeRow {
            row,
            delay: self.tail_delay,
        });
        if self.recycled == 0 {
            self.tail_delay = 0;
        } else {
            self.recycled -= 1;
            self.tail_delay = 2;
        }
    }
}
