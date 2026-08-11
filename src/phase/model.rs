use std::io::Write;
use std::sync::Arc;

use rsomics_common::{Result, RsomicsError};

const MAX_VARS: usize = 256;
const MAX_COUNT_CELLS: usize = 16 * 1024 * 1024;
const FLIP_PENALTY: i32 = 2;
const FLIP_THRESHOLD: i32 = 4;
const MASK_THRESHOLD: i32 = 3;

#[derive(Clone, Copy)]
pub(super) struct Site {
    pub(super) position: i64,
    pub(super) alleles: [u8; 2],
}

pub(super) struct Fragment {
    pub(super) key: Arc<[u8]>,
    pub(super) first_site: usize,
    pub(super) sequence: Vec<u8>,
    pub(super) alignment_start: i32,
    pub(super) phase: u8,
    pub(super) phased: bool,
    pub(super) ambiguous: bool,
    pub(super) flipped: bool,
    pub(super) in_phase: u16,
    pub(super) out_phase: u16,
}

impl Fragment {
    pub(super) fn new(key: Arc<[u8]>, site: usize, allele: u8, alignment_start: i32) -> Self {
        Self {
            key,
            first_site: site,
            sequence: vec![allele],
            alignment_start,
            phase: 0,
            phased: false,
            ambiguous: false,
            flipped: false,
            in_phase: 0,
            out_phase: 0,
        }
    }

    pub(super) fn push(&mut self, site: usize, allele: u8) {
        let length = site - self.first_site + 1;
        if length < MAX_VARS {
            self.sequence.resize(length, 0);
            self.sequence[length - 1] = allele;
        }
    }
}

pub(super) fn write_block(
    output: &mut impl Write,
    reference: &str,
    sites: &[Site],
    fragments: &mut Vec<Fragment>,
    window: usize,
    fix_chimeras: bool,
    marker_offset: usize,
) -> Result<usize> {
    if sites.is_empty() {
        return Ok(0);
    }
    clean_fragments(fragments);
    let start = sites[0].position + 1;
    let end = sites.last().unwrap().position + 1;
    writeln!(output, "PS\t{reference}\t{start}\t{end}").map_err(RsomicsError::Io)?;

    if sites.len() == 1 {
        let site = sites[0];
        writeln!(
            output,
            "M0\t{reference}\t{start}\t{start}\t{}\t{}\t{}\t0\t0\t0\t0",
            base(site.alleles[0]),
            base(site.alleles[1]),
            marker_offset + 1
        )
        .map_err(RsomicsError::Io)?;
        for fragment in fragments.iter_mut() {
            if fragment.first_site != 0 {
                continue;
            }
            fragment.flipped = false;
            fragment.phased = fragment.sequence[0] != 0;
            if fragment.phased {
                fragment.phase = fragment.sequence[0] - 1;
            }
        }
        writeln!(output, "//").map_err(RsomicsError::Io)?;
        return Ok(1);
    }

    let weights = count_all(window, sites.len(), fragments)?;
    let path = dynamic_program(window, sites.len(), &weights)?;
    let counts = phase_fragments(sites.len(), &path, fragments, false);
    let masks = masks(sites.len(), &counts);
    let mut filtered = vec![false; sites.len()];
    for &(first, last) in &masks {
        for value in &mut filtered[first..=last] {
            *value = true;
        }
        writeln!(
            output,
            "FL\t{reference}\t{}\t{}",
            sites[first].position + 1,
            sites[last].position + 1
        )
        .map_err(RsomicsError::Io)?;
    }
    let counts = if fix_chimeras {
        phase_fragments(sites.len(), &path, fragments, true)
    } else {
        counts
    };

    for (index, site) in sites.iter().enumerate() {
        let packed = counts[index];
        let ordered = [
            site.alleles[path[index] as usize],
            site.alleles[1 - path[index] as usize],
        ];
        writeln!(
            output,
            "M{}\t{reference}\t{start}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            usize::from(filtered[index]) + 1,
            site.position + 1,
            base(ordered[0]),
            base(ordered[1]),
            marker_offset + index + 1,
            packed & 0xffff,
            packed >> 16 & 0xffff,
            packed >> 32 & 0xffff,
            packed >> 48 & 0xffff,
        )
        .map_err(RsomicsError::Io)?;
    }

    fragments.sort_by(|left, right| {
        left.first_site
            .cmp(&right.first_site)
            .then_with(|| left.key.cmp(&right.key))
    });
    for fragment in fragments
        .iter()
        .filter(|fragment| fragment.sequence.len() > 1)
    {
        let sequence: String = fragment
            .sequence
            .iter()
            .enumerate()
            .map(|(offset, &allele)| {
                if allele == 0 {
                    'N'
                } else {
                    base(sites[fragment.first_site + offset].alleles[(allele - 1) as usize])
                }
            })
            .collect();
        writeln!(
            output,
            "EV\t0\t{reference}\t{}\t40\t{}M\t*\t0\t0\t{sequence}\t*\tYP:i:{}\tYF:i:{}\tYI:i:{}\tYO:i:{}\tYS:i:{}",
            fragment.first_site + marker_offset + 1,
            fragment.sequence.len(),
            fragment.phase,
            u8::from(fragment.flipped),
            fragment.in_phase,
            fragment.out_phase,
            fragment.alignment_start + 1,
        )
        .map_err(RsomicsError::Io)?;
    }
    writeln!(output, "//").map_err(RsomicsError::Io)?;
    Ok(sites.len())
}

fn clean_fragments(fragments: &mut Vec<Fragment>) {
    fragments.retain_mut(|fragment| {
        let Some(first) = fragment.sequence.iter().position(|&allele| allele != 0) else {
            return false;
        };
        let last = fragment
            .sequence
            .iter()
            .rposition(|&allele| allele != 0)
            .unwrap();
        fragment.first_site += first;
        fragment.sequence = fragment.sequence[first..=last].to_vec();
        true
    });
}

fn count_all(window: usize, site_count: usize, fragments: &[Fragment]) -> Result<Vec<Vec<i32>>> {
    let (states, cells) = workspace(window, site_count)?;
    let mut result = Vec::new();
    result.try_reserve_exact(site_count).map_err(|_| {
        RsomicsError::ConfigError(format!("phase window requires {cells} count cells"))
    })?;
    for _ in 0..site_count {
        let mut row = Vec::new();
        row.try_reserve_exact(states).map_err(|_| {
            RsomicsError::ConfigError(format!("phase window requires {cells} count cells"))
        })?;
        row.resize(states, 0_i32);
        result.push(row);
    }
    let mut local = vec![0_u8; window];

    for fragment in fragments {
        if fragment.sequence.len() <= 1 || fragment.first_site >= site_count {
            continue;
        }
        for offset in 1..fragment.sequence.len() {
            for (index, value) in local.iter_mut().enumerate() {
                *value = if offset < window - 1 - index {
                    0
                } else {
                    fragment.sequence[offset - (window - 1 - index)]
                };
            }
            count_pattern(&local, &mut result[fragment.first_site + offset]);
        }
    }
    Ok(result)
}

pub(super) fn ensure_workspace(window: usize, site_count: usize) -> Result<()> {
    if site_count > 1 {
        workspace(window, site_count)?;
    }
    Ok(())
}

fn workspace(window: usize, site_count: usize) -> Result<(usize, usize)> {
    let states = 1usize
        .checked_shl(window as u32)
        .ok_or_else(|| RsomicsError::ConfigError("phase window is too large".to_owned()))?;
    let cells = states
        .checked_mul(site_count)
        .ok_or_else(|| RsomicsError::ConfigError("phase workspace size overflows".to_owned()))?;
    if cells > MAX_COUNT_CELLS {
        return Err(RsomicsError::ConfigError(format!(
            "phase window requires {cells} count cells; the limit is {MAX_COUNT_CELLS}"
        )));
    }
    Ok((states, cells))
}

fn count_pattern(sequence: &[u8], counts: &mut [i32]) {
    if sequence.last() == Some(&0) {
        return;
    }
    let mut fixed = 0;
    let mut variable = 0;
    let mut observed = 0;
    for (index, &allele) in sequence.iter().enumerate() {
        let bit = sequence.len() - index - 1;
        if allele == 0 {
            variable |= 1 << bit;
        } else {
            fixed |= usize::from(allele - 1) << bit;
            observed += 1;
        }
    }
    if observed <= 1 {
        return;
    }
    let mut subset = variable;
    loop {
        counts[fixed | subset] += 1;
        if subset == 0 {
            break;
        }
        subset = (subset - 1) & variable;
    }
}

fn dynamic_program(window: usize, sites: usize, weights: &[Vec<i32>]) -> Result<Vec<u8>> {
    let state_count = 1usize
        .checked_shl((window - 1) as u32)
        .ok_or_else(|| RsomicsError::ConfigError("phase window is too large".to_owned()))?;
    let mask = (1usize << window) - 1;
    let mut previous = vec![0_i64; state_count];
    let mut current = vec![0_i64; state_count];
    let mut backtrack = vec![vec![false; state_count]; sites];
    for site in 0..sites {
        for state in 0..state_count {
            let complement = !state & mask;
            let score = i64::from(weights[site][state]) + i64::from(weights[site][complement]);
            let left = previous[state >> 1] + score;
            let right = previous[complement >> 1] + score;
            if left > right {
                current[state] = left;
            } else {
                current[state] = right;
                backtrack[site][state] = true;
            }
        }
        std::mem::swap(&mut previous, &mut current);
    }
    let mut state = 0;
    let mut best = 0;
    for (index, &score) in previous.iter().enumerate() {
        if score > best {
            best = score;
            state = index;
        }
    }
    let mut which = false;
    let mut path = vec![0_u8; sites];
    for site in (0..sites).rev() {
        path[site] = if which {
            (!state & 1) as u8
        } else {
            (state & 1) as u8
        };
        if backtrack[site][state] {
            which = !which;
            state = (!state & mask) >> 1;
        } else {
            state >>= 1;
        }
    }
    Ok(path)
}

fn phase_fragments(
    sites: usize,
    path: &[u8],
    fragments: &mut [Fragment],
    repair: bool,
) -> Vec<u64> {
    let mut counts = vec![0_u64; sites];
    for fragment in fragments {
        if fragment.first_site >= sites {
            continue;
        }
        let mut matches = [0_u16; 2];
        for (offset, &allele) in fragment.sequence.iter().enumerate() {
            if allele != 0 {
                matches[usize::from(allele != path[fragment.first_site + offset] + 1)] += 1;
            }
        }
        fragment.phase = u8::from(matches[0] <= matches[1]);
        fragment.in_phase = matches[fragment.phase as usize];
        fragment.out_phase = matches[1 - fragment.phase as usize];
        fragment.phased = fragment.in_phase != fragment.out_phase;
        fragment.ambiguous = fragment.in_phase != 0
            && fragment.out_phase != 0
            && fragment.out_phase < 3
            && fragment.in_phase <= fragment.out_phase + 1;
        fragment.flipped = false;
        if repair && matches[0] >= 3 && matches[1] >= 3 {
            repair_chimera(fragment, path, matches);
        }
        if fragment.sequence.len() == 1 {
            continue;
        }
        for (offset, &allele) in fragment.sequence.iter().enumerate() {
            if allele == 0 {
                continue;
            }
            let site = fragment.first_site + offset;
            let observed = if fragment.phase == 0 {
                allele - 1
            } else {
                2 - allele
            };
            let shift = match (fragment.phase, observed == path[site]) {
                (0, true) => 0,
                (0, false) => 16,
                (1, true) => 32,
                (1, false) => 48,
                _ => unreachable!(),
            };
            counts[site] += 1_u64 << shift;
        }
    }
    counts
}

fn repair_chimera(fragment: &mut Fragment, path: &[u8], matches: [u16; 2]) {
    let length = fragment.sequence.len();
    let mut left = vec![[0_i32; 2]; length];
    let mut right = vec![[0_i32; 2]; length];
    let mut sum = [0_i32; 2];
    for (offset, &allele) in fragment.sequence.iter().enumerate() {
        if allele != 0 {
            let observed = if fragment.phase == 0 {
                allele - 1
            } else {
                2 - allele
            };
            sum[usize::from(observed != path[fragment.first_site + offset])] += 1;
        }
        left[offset] = sum;
    }
    sum = [0, 0];
    for offset in (0..length).rev() {
        let allele = fragment.sequence[offset];
        if allele != 0 {
            let observed = if fragment.phase == 0 {
                allele - 1
            } else {
                2 - allele
            };
            sum[usize::from(observed != path[fragment.first_site + offset])] += 1;
        }
        right[offset] = sum;
    }
    let mut best = 0;
    let mut split = None;
    let mut tail = false;
    for offset in 0..length - 1 {
        let scores = [
            left[offset][0] + right[offset + 1][1] - right[offset + 1][0] * FLIP_PENALTY,
            left[offset][1] + right[offset + 1][0] - right[offset + 1][1] * FLIP_PENALTY,
        ];
        let (score, flip_tail) = if scores[0] > scores[1] {
            (scores[0], true)
        } else {
            (scores[1], false)
        };
        if score > best {
            best = score;
            split = Some(offset);
            tail = flip_tail;
        }
    }
    if best - i32::from(matches[0]) >= FLIP_THRESHOLD
        && best - i32::from(matches[1]) >= FLIP_THRESHOLD
    {
        fragment.flipped = true;
        let split = split.unwrap();
        let range = if tail {
            split + 1..length
        } else {
            0..split + 1
        };
        for allele in &mut fragment.sequence[range] {
            if *allele != 0 {
                *allele = 3 - *allele;
            }
        }
    }
}

fn masks(sites: usize, counts: &[u64]) -> Vec<(usize, usize)> {
    let mut result = Vec::new();
    let mut score = 0;
    let mut best = 0;
    let mut best_index = 0;
    let mut start = 0;
    let mut index = 0;
    while index < sites {
        let value = counts[index];
        let c = [
            (value & 0xffff) as i32,
            (value >> 16 & 0xffff) as i32,
            (value >> 32 & 0xffff) as i32,
            (value >> 48 & 0xffff) as i32,
        ];
        let previous = score;
        let mut delta = if c[1] + c[3] == 0 {
            -(c[0] + c[2])
        } else {
            c[1] + c[3] - 1
        };
        delta += (c[3] - c[2]).max(0) + (c[1] - c[0]).max(0);
        score = (score + delta).max(0);
        if previous == 0 && score > 0 {
            start = index;
        }
        if (index + 1 == sites || score == 0) && best >= MASK_THRESHOLD {
            result.push((start, best_index));
            index = best_index;
            score = 0;
        } else if score > best {
            best = score;
            best_index = index;
        }
        if score == 0 {
            best = 0;
        }
        index += 1;
    }
    result
}

fn base(value: u8) -> char {
    b"ACGTX"[value as usize] as char
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_count_pattern(sequence: &[u8], counts: &mut [i32]) {
        if sequence.last() == Some(&0) {
            return;
        }
        let ambiguous = sequence.iter().filter(|&&allele| allele == 0).count();
        if sequence.len() - ambiguous <= 1 {
            return;
        }
        for choice in 0..1usize << ambiguous {
            let mut value = 0;
            let mut bit = 0;
            for &allele in sequence {
                let allele = if allele == 0 {
                    let selected = (choice >> bit & 1) as u8 + 1;
                    bit += 1;
                    selected
                } else {
                    allele
                };
                value = value << 1 | usize::from(allele - 1);
            }
            counts[value] += 1;
        }
    }

    #[test]
    fn count_workspace_is_bounded_before_allocation() {
        let error = count_all(13, MAX_COUNT_CELLS / (1 << 13) + 1, &[]).unwrap_err();
        assert!(error.to_string().contains("count cells"));
    }

    #[test]
    fn pattern_enumeration_matches_the_naive_definition() {
        for length in 1..=7 {
            for encoded in 0..3usize.pow(length as u32) {
                let mut value = encoded;
                let mut sequence = vec![0; length];
                for allele in sequence.iter_mut().rev() {
                    *allele = (value % 3) as u8;
                    value /= 3;
                }
                let mut actual = vec![0; 1 << length];
                let mut expected = vec![0; 1 << length];
                count_pattern(&sequence, &mut actual);
                naive_count_pattern(&sequence, &mut expected);
                assert_eq!(actual, expected, "{sequence:?}");
            }
        }
    }

    #[test]
    fn long_phase_scores_do_not_overflow() {
        let sites = 16_385;
        let weights = vec![vec![65_535; 2]; sites];
        assert_eq!(dynamic_program(1, sites, &weights).unwrap().len(), sites);
    }
}
