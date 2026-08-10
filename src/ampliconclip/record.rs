use std::borrow::Cow;

use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use smallvec::SmallVec;

const POS: usize = 4;
const L_READ_NAME: usize = 8;
const MAPQ: usize = 9;
const N_CIGAR: usize = 12;
const BIN: usize = 10;
const FLAG: usize = 14;
const L_SEQ: usize = 16;
const FIXED_HEAD: usize = 32;

const FLAG_UNMAPPED: u16 = 0x4;

const CIGAR_SOFT_CLIP: u32 = 4;
const CIGAR_HARD_CLIP: u32 = 5;

fn cigar_type(op: u32) -> u32 {
    (0x0003_C1A7u32 >> (op << 1)) & 3
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Clipping {
    Soft,
    Hard,
}

struct RawBam<'a> {
    bytes: Cow<'a, [u8]>,
}

impl RawBam<'_> {
    pub fn flags(&self) -> u16 {
        u16::from_le_bytes([self.bytes[FLAG], self.bytes[FLAG + 1]])
    }

    pub fn set_flag_bits(&mut self, bits: u16) {
        let new = self.flags() | bits;
        self.bytes.to_mut()[FLAG..FLAG + 2].copy_from_slice(&new.to_le_bytes());
    }

    pub fn pos(&self) -> i64 {
        i64::from(i32::from_le_bytes(
            self.bytes[POS..POS + 4].try_into().unwrap(),
        ))
    }

    fn name_len(&self) -> usize {
        usize::from(self.bytes[L_READ_NAME])
    }

    fn n_cigar(&self) -> usize {
        usize::from(u16::from_le_bytes([
            self.bytes[N_CIGAR],
            self.bytes[N_CIGAR + 1],
        ]))
    }

    pub fn l_qseq(&self) -> usize {
        u32::from_le_bytes(self.bytes[L_SEQ..L_SEQ + 4].try_into().unwrap()) as usize
    }

    fn cigar_start(&self) -> usize {
        FIXED_HEAD + self.name_len()
    }

    fn seq_start(&self) -> usize {
        self.cigar_start() + self.n_cigar() * 4
    }

    fn qual_start(&self) -> usize {
        self.seq_start() + self.l_qseq().div_ceil(2)
    }

    fn aux_start(&self) -> usize {
        self.qual_start() + self.l_qseq()
    }

    fn cigar(&self) -> SmallVec<[(u32, u32); 8]> {
        let start = self.cigar_start();
        (0..self.n_cigar())
            .map(|i| {
                let off = start + i * 4;
                let raw = u32::from_le_bytes(self.bytes[off..off + 4].try_into().unwrap());
                (raw & 0xf, raw >> 4)
            })
            .collect()
    }

    pub fn endpos(&self) -> i64 {
        let span: i64 = self
            .cigar()
            .iter()
            .filter(|(op, _)| cigar_type(*op) & 2 != 0)
            .map(|(_, len)| i64::from(*len))
            .sum();
        self.pos() + if span > 0 { span } else { 1 }
    }

    pub fn unmap(&self) -> RawBam<'static> {
        let name_len = self.name_len();
        let l_qseq = self.l_qseq();
        let seq_start = self.seq_start();
        let aux_start = self.aux_start();

        let mut out = Vec::with_capacity(
            FIXED_HEAD + name_len + l_qseq.div_ceil(2) + l_qseq + self.bytes.len() - aux_start,
        );
        out.extend_from_slice(&self.bytes[..FIXED_HEAD + name_len]);
        out.extend_from_slice(&self.bytes[seq_start..aux_start]);
        out.extend_from_slice(&self.bytes[aux_start..]);

        let mut rec = RawBam {
            bytes: Cow::Owned(out),
        };
        rec.bytes.to_mut()[N_CIGAR..N_CIGAR + 2].copy_from_slice(&0u16.to_le_bytes());
        rec.bytes.to_mut()[MAPQ] = 0;
        rec.set_flag_bits(FLAG_UNMAPPED);
        rec
    }
}

fn cigar_gen(len: u32, op: u32) -> u32 {
    (len << 4) | op
}

fn rebuild(
    src: &RawBam<'_>,
    new_cigar: &[u32],
    qry_removed: usize,
    new_pos: i64,
    from_left: bool,
) -> Result<RawBam<'static>> {
    let name_len = src.name_len();
    let l_qseq = src.l_qseq();
    let new_l_qseq = l_qseq.checked_sub(qry_removed).ok_or_else(|| {
        RsomicsError::InvalidInput("primer clip exceeds the query length".to_owned())
    })?;
    let aux_start = src.aux_start();
    let aux = &src.bytes[aux_start..];

    let mut out = Vec::with_capacity(
        FIXED_HEAD
            + name_len
            + new_cigar.len() * 4
            + new_l_qseq.div_ceil(2)
            + new_l_qseq
            + aux.len(),
    );

    out.extend_from_slice(&src.bytes[..FIXED_HEAD + name_len]);

    for &op in new_cigar {
        out.extend_from_slice(&op.to_le_bytes());
    }

    let src_seq_start = src.seq_start();
    let src_seq = &src.bytes[src_seq_start..src_seq_start + l_qseq.div_ceil(2)];
    if from_left {
        append_seq_drop_head(src_seq, &mut out, l_qseq, qry_removed);
    } else {
        out.extend_from_slice(&src_seq[..new_l_qseq.div_ceil(2)]);
        if !new_l_qseq.is_multiple_of(2) {
            let last = out.len() - 1;
            out[last] &= 0xf0;
        }
    }

    let src_qual_start = src.qual_start();
    let src_qual = &src.bytes[src_qual_start..src_qual_start + l_qseq];
    if from_left {
        out.extend_from_slice(&src_qual[qry_removed..]);
    } else {
        out.extend_from_slice(&src_qual[..new_l_qseq]);
    }

    out.extend_from_slice(aux);

    let mut rec = RawBam {
        bytes: Cow::Owned(out),
    };
    let n = u16::try_from(new_cigar.len()).map_err(|_| {
        RsomicsError::InvalidInput("clipped CIGAR exceeds the BAM operation limit".to_owned())
    })?;
    rec.bytes.to_mut()[N_CIGAR..N_CIGAR + 2].copy_from_slice(&n.to_le_bytes());
    let lq = u32::try_from(new_l_qseq)
        .map_err(|_| RsomicsError::InvalidInput("query length exceeds u32".to_owned()))?;
    rec.bytes.to_mut()[L_SEQ..L_SEQ + 4].copy_from_slice(&lq.to_le_bytes());
    let p = i32::try_from(new_pos)
        .map_err(|_| RsomicsError::InvalidInput("clipped position exceeds i32".to_owned()))?;
    rec.bytes.to_mut()[POS..POS + 4].copy_from_slice(&p.to_le_bytes());
    Ok(rec)
}

fn append_seq_drop_head(src_seq: &[u8], out: &mut Vec<u8>, l_qseq: usize, drop: usize) {
    let new_len = l_qseq - drop;
    if drop.is_multiple_of(2) {
        out.extend_from_slice(&src_seq[drop / 2..drop / 2 + new_len.div_ceil(2)]);
    } else {
        let mut in_idx = drop / 2;
        let mut i = drop;
        while i < l_qseq - 1 {
            out.push(((src_seq[in_idx] & 0x0f) << 4) | ((src_seq[in_idx + 1] & 0xf0) >> 4));
            in_idx += 1;
            i += 2;
        }
        if i < l_qseq {
            out.push((src_seq[in_idx] & 0x0f) << 4);
        }
    }
}

fn trim_left(src: &RawBam<'_>, bases: u32, clipping: Clipping) -> Result<RawBam<'static>> {
    let cigar = src.cigar();
    let mut ref_remove = bases;
    let mut qry_removed: u32 = 0;
    let mut hardclip: u32 = 0;
    let mut new_pos = src.pos();
    let n = cigar.len();

    let mut i = 0;
    while i < n {
        let (op, oplen) = cigar[i];
        let ctype = cigar_type(op);
        if op == CIGAR_HARD_CLIP {
            hardclip += oplen;
        } else {
            if ctype & 2 != 0 {
                if oplen <= ref_remove {
                    ref_remove -= oplen;
                } else {
                    break;
                }
                new_pos += i64::from(oplen);
            }
            if ctype & 1 != 0 {
                qry_removed += oplen;
            }
        }
        i += 1;
    }

    if i < n {
        let (op, _) = cigar[i];
        let ctype = cigar_type(op);
        if ctype & 2 != 0 {
            new_pos += i64::from(ref_remove);
        }
        if ctype & 1 != 0 {
            qry_removed += ref_remove;
        }
    } else {
        if clipping == Clipping::Hard {
            return empty_read(src);
        }
        qry_removed = src.l_qseq() as u32;
    }

    let mut new_cigar = SmallVec::<[u32; 8]>::with_capacity(n + 2);
    match clipping {
        Clipping::Hard => {
            if hardclip + qry_removed > 0 {
                new_cigar.push(cigar_gen(hardclip + qry_removed, CIGAR_HARD_CLIP));
            }
        }
        Clipping::Soft => {
            if hardclip > 0 {
                new_cigar.push(cigar_gen(hardclip, CIGAR_HARD_CLIP));
            }
            if qry_removed > 0 {
                new_cigar.push(cigar_gen(qry_removed, CIGAR_SOFT_CLIP));
            }
        }
    }

    if i < n {
        let (op, oplen) = cigar[i];
        if oplen > ref_remove {
            new_cigar.push(cigar_gen(oplen - ref_remove, op));
            for &(o, l) in &cigar[i + 1..] {
                new_cigar.push(cigar_gen(l, o));
            }
        }
    }

    let phys_removed = if clipping == Clipping::Soft {
        0
    } else {
        qry_removed as usize
    };

    rebuild(src, &new_cigar, phys_removed, new_pos, true)
}

fn trim_right(src: &RawBam<'_>, bases: u32, clipping: Clipping) -> Result<RawBam<'static>> {
    let cigar = src.cigar();
    let mut ref_remove = bases;
    let mut qry_removed: u32 = 0;
    let mut hardclip: u32 = 0;
    let n = cigar.len() as i64;

    let mut i: i64 = n - 1;
    while i >= 0 {
        let (op, oplen) = cigar[i as usize];
        let ctype = cigar_type(op);
        if op == CIGAR_HARD_CLIP {
            hardclip += oplen;
        } else {
            if ctype & 2 != 0 {
                if oplen <= ref_remove {
                    ref_remove -= oplen;
                } else {
                    break;
                }
            }
            if ctype & 1 != 0 {
                qry_removed += oplen;
            }
        }
        i -= 1;
    }

    let cap = cigar.len() + 2;
    let mut slots = SmallVec::<[u32; 16]>::from_elem(0, cap);
    let mut new_n_cigar: usize = 0;
    let mut j: i64;

    if i >= 0 {
        let ctype = cigar_type(cigar[i as usize].0);
        if ctype & 1 != 0 {
            qry_removed += ref_remove;
        }
        j = i;
        if qry_removed > 0 {
            j += 1;
        }
        if hardclip > 0 && (clipping == Clipping::Soft || qry_removed == 0) {
            j += 1;
        }
    } else {
        if clipping == Clipping::Hard {
            return empty_read(src);
        }
        qry_removed = src.l_qseq() as u32;
        j = 0;
        if hardclip > 0 && clipping == Clipping::Soft {
            j += 1;
        }
    }

    if clipping == Clipping::Hard && hardclip + qry_removed > 0 {
        slots[j as usize] = cigar_gen(hardclip + qry_removed, CIGAR_HARD_CLIP);
        new_n_cigar += 1;
    }
    if clipping == Clipping::Soft {
        if hardclip > 0 {
            slots[j as usize] = cigar_gen(hardclip, CIGAR_HARD_CLIP);
            new_n_cigar += 1;
            if qry_removed > 0 {
                j -= 1;
            }
        }
        if qry_removed > 0 {
            slots[j as usize] = cigar_gen(qry_removed, CIGAR_SOFT_CLIP);
            new_n_cigar += 1;
        }
    }

    if j > 0 {
        j -= 1;
        let (op, oplen) = cigar[i as usize];
        slots[j as usize] = cigar_gen(oplen - ref_remove, op);
        new_n_cigar += 1;
    }

    while j > 0 {
        j -= 1;
        i -= 1;
        let (op, oplen) = cigar[i as usize];
        slots[j as usize] = cigar_gen(oplen, op);
        new_n_cigar += 1;
    }

    let phys_removed = if clipping == Clipping::Soft {
        0
    } else {
        qry_removed as usize
    };

    rebuild(src, &slots[..new_n_cigar], phys_removed, src.pos(), false)
}

fn empty_read(src: &RawBam<'_>) -> Result<RawBam<'static>> {
    let name_len = src.name_len();
    let aux_start = src.aux_start();
    let aux = &src.bytes[aux_start..];
    let mut out = Vec::with_capacity(FIXED_HEAD + name_len + aux.len());
    out.extend_from_slice(&src.bytes[..FIXED_HEAD + name_len]);
    out.extend_from_slice(aux);
    let mut rec = RawBam {
        bytes: Cow::Owned(out),
    };
    rec.bytes.to_mut()[N_CIGAR..N_CIGAR + 2].copy_from_slice(&0u16.to_le_bytes());
    rec.bytes.to_mut()[L_SEQ..L_SEQ + 4].copy_from_slice(&0u32.to_le_bytes());
    Ok(rec)
}

pub(crate) fn end_position(record: &RawRecord) -> Result<i64> {
    let span = record
        .cigar_ops()
        .filter(|(op, _)| cigar_type(u32::from(*op)) & 2 != 0)
        .map(|(_, len)| i64::from(len))
        .sum::<i64>();
    Ok(i64::from(record.alignment_start()) + span.max(1))
}

pub(crate) fn active_query_len(record: &RawRecord) -> Result<i64> {
    Ok(record
        .cigar_ops()
        .filter(|(op, _)| cigar_type(u32::from(*op)) & 1 != 0 && *op != CIGAR_SOFT_CLIP as u8)
        .map(|(_, len)| i64::from(len))
        .sum())
}

pub(crate) fn clip_left(record: &RawRecord, bases: u32, clipping: Clipping) -> Result<RawRecord> {
    finalize(trim_left(&from_record(record)?, bases, clipping)?)
}

pub(crate) fn clip_right(record: &RawRecord, bases: u32, clipping: Clipping) -> Result<RawRecord> {
    finalize(trim_right(&from_record(record)?, bases, clipping)?)
}

pub(crate) fn unmap(record: &RawRecord) -> Result<RawRecord> {
    finalize(from_record(record)?.unmap())
}

fn from_record(record: &RawRecord) -> Result<RawBam<'_>> {
    Ok(RawBam {
        bytes: Cow::Borrowed(record.as_bytes()),
    })
}

pub(crate) fn validate(record: &RawRecord) -> Result<()> {
    let mut cigar = record.cigar_ops();
    let first = cigar.next();
    let second = cigar.next();
    let third = cigar.next();
    if first == Some((CIGAR_SOFT_CLIP as u8, record.sequence_len() as u32))
        && second.is_some_and(|(op, _)| op == 3)
        && third.is_none()
        && record.aux_type(*b"CG") == Some(b'B')
    {
        return Err(RsomicsError::InvalidInput(format!(
            "read {}: long CIGAR clipping is not supported",
            String::from_utf8_lossy(record.name())
        )));
    }
    if record.cigar_ops().any(|(op, len)| op > 8 || len == 0) {
        return Err(RsomicsError::InvalidInput(format!(
            "read {}: invalid CIGAR operation",
            String::from_utf8_lossy(record.name())
        )));
    }
    Ok(())
}

fn finalize(mut record: RawBam<'_>) -> Result<RawRecord> {
    let start = record.pos();
    let end = record.endpos();
    let bin = if start < 0 {
        4680
    } else {
        reg2bin(start, end)?
    };
    record.bytes.to_mut()[BIN..BIN + 2].copy_from_slice(&bin.to_le_bytes());
    RawRecord::try_from(record.bytes.into_owned())
}

fn reg2bin(start: i64, end: i64) -> Result<u16> {
    if end <= start {
        return Err(RsomicsError::InvalidInput(
            "cannot calculate BAM bin for invalid coordinates".to_owned(),
        ));
    }
    let end = end - 1;
    let bin = if start >> 14 == end >> 14 {
        ((1 << 15) - 1) / 7 + (start >> 14)
    } else if start >> 17 == end >> 17 {
        ((1 << 12) - 1) / 7 + (start >> 17)
    } else if start >> 20 == end >> 20 {
        ((1 << 9) - 1) / 7 + (start >> 20)
    } else if start >> 23 == end >> 23 {
        ((1 << 6) - 1) / 7 + (start >> 23)
    } else if start >> 26 == end >> 26 {
        1 + (start >> 26)
    } else {
        0
    };
    u16::try_from(bin).map_err(|_| RsomicsError::InvalidInput("BAM bin exceeds u16".to_owned()))
}
