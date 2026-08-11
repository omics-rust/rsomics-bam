use std::io::Write;

use rsomics_common::{Result, RsomicsError};

use super::record::{READ1, READ2, REVERSE};

#[derive(Clone, Copy)]
pub(super) struct Bed<'a> {
    pub reference: &'a str,
    pub start: u64,
    pub end: u64,
    pub name: &'a [u8],
    pub flags: u16,
    pub score: i64,
}

pub(super) fn bed6(output: &mut impl Write, bed: Bed<'_>, cigar: Option<&str>) -> Result<()> {
    output
        .write_all(bed.reference.as_bytes())
        .map_err(RsomicsError::Io)?;
    write!(output, "\t{}\t{}\t", bed.start, bed.end).map_err(RsomicsError::Io)?;
    output.write_all(bed.name).map_err(RsomicsError::Io)?;
    if bed.flags & READ1 != 0 {
        output.write_all(b"/1").map_err(RsomicsError::Io)?;
    }
    if bed.flags & READ2 != 0 {
        output.write_all(b"/2").map_err(RsomicsError::Io)?;
    }
    write!(
        output,
        "\t{}\t{}",
        bed.score,
        if bed.flags & REVERSE != 0 { '-' } else { '+' }
    )
    .map_err(RsomicsError::Io)?;
    if let Some(cigar) = cigar {
        write!(output, "\t{cigar}").map_err(RsomicsError::Io)?;
    }
    output.write_all(b"\n").map_err(RsomicsError::Io)
}

pub(super) fn bed12(
    output: &mut impl Write,
    bed: Bed<'_>,
    color: &str,
    blocks: &[(u64, u64)],
) -> Result<()> {
    output
        .write_all(bed.reference.as_bytes())
        .map_err(RsomicsError::Io)?;
    write!(output, "\t{}\t{}\t", bed.start, bed.end).map_err(RsomicsError::Io)?;
    output.write_all(bed.name).map_err(RsomicsError::Io)?;
    if bed.flags & READ1 != 0 {
        output.write_all(b"/1").map_err(RsomicsError::Io)?;
    }
    if bed.flags & READ2 != 0 {
        output.write_all(b"/2").map_err(RsomicsError::Io)?;
    }
    write!(
        output,
        "\t{}\t{}\t{}\t{}\t{color}\t{}\t",
        bed.score,
        if bed.flags & REVERSE != 0 { '-' } else { '+' },
        bed.start,
        bed.end,
        blocks.len()
    )
    .map_err(RsomicsError::Io)?;
    for (index, &(block_start, block_end)) in blocks.iter().enumerate() {
        if index > 0 {
            output.write_all(b",").map_err(RsomicsError::Io)?;
        }
        write!(output, "{}", block_end - block_start).map_err(RsomicsError::Io)?;
    }
    output.write_all(b"\t").map_err(RsomicsError::Io)?;
    for (index, &(block_start, _)) in blocks.iter().enumerate() {
        if index > 0 {
            output.write_all(b",").map_err(RsomicsError::Io)?;
        }
        write!(output, "{}", block_start - bed.start).map_err(RsomicsError::Io)?;
    }
    output.write_all(b"\n").map_err(RsomicsError::Io)
}
