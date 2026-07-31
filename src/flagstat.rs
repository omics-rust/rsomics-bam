use std::fmt;
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rust_htslib::bam::{Read, Record};
use serde::Serialize;

use crate::input;

const PAIRED: u16 = 0x1;
const PROPER_PAIR: u16 = 0x2;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const READ1: u16 = 0x40;
const READ2: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const QCFAIL: u16 = 0x200;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize)]
pub struct Counts {
    pub total: [u64; 2],
    pub primary: [u64; 2],
    pub secondary: [u64; 2],
    pub supplementary: [u64; 2],
    pub duplicates: [u64; 2],
    pub primary_duplicates: [u64; 2],
    pub mapped: [u64; 2],
    pub primary_mapped: [u64; 2],
    pub paired: [u64; 2],
    pub read1: [u64; 2],
    pub read2: [u64; 2],
    pub properly_paired: [u64; 2],
    pub both_mapped: [u64; 2],
    pub singletons: [u64; 2],
    pub mate_different_reference: [u64; 2],
    pub mate_different_reference_mapq5: [u64; 2],
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Options<'a> {
    pub reference: Option<&'a Path>,
    pub additional_threads: usize,
}

pub fn count(input_path: &Path, options: Options<'_>) -> Result<Counts> {
    let mut reader = input::open(input_path, options.reference, options.additional_threads)?;
    let mut counts = Counts::default();
    let mut record = Record::new();

    while let Some(result) = reader.read(&mut record) {
        result.map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "reading alignment record from {}: {error}",
                input_path.display()
            ))
        })?;
        counts.tally(&record);
    }

    Ok(counts)
}

impl Counts {
    pub fn tally(&mut self, record: &Record) {
        let flags = record.flags();
        let category = usize::from(flags & QCFAIL != 0);
        let secondary = flags & SECONDARY != 0;
        let supplementary = flags & SUPPLEMENTARY != 0;
        let primary = !secondary && !supplementary;
        let mapped = flags & UNMAPPED == 0;

        self.total[category] += 1;

        if secondary {
            self.secondary[category] += 1;
        } else if supplementary {
            self.supplementary[category] += 1;
        } else {
            self.primary[category] += 1;
            if flags & PAIRED != 0 {
                self.paired[category] += 1;
                self.read1[category] += u64::from(flags & READ1 != 0);
                self.read2[category] += u64::from(flags & READ2 != 0);

                if flags & PROPER_PAIR != 0 && mapped {
                    self.properly_paired[category] += 1;
                }

                let mate_mapped = flags & MATE_UNMAPPED == 0;
                if mapped && mate_mapped {
                    self.both_mapped[category] += 1;
                    if record.tid() != record.mtid() {
                        self.mate_different_reference[category] += 1;
                        if record.mapq() >= 5 {
                            self.mate_different_reference_mapq5[category] += 1;
                        }
                    }
                } else if mapped {
                    self.singletons[category] += 1;
                }
            }
        }

        if mapped {
            self.mapped[category] += 1;
            if primary {
                self.primary_mapped[category] += 1;
            }
        }

        if flags & DUPLICATE != 0 {
            self.duplicates[category] += 1;
            if primary {
                self.primary_duplicates[category] += 1;
            }
        }
    }

    pub fn to_tsv(&self) -> String {
        let mut output = String::new();
        for (values, label) in [
            (
                pair(self.total),
                "total (QC-passed reads + QC-failed reads)",
            ),
            (pair(self.primary), "primary"),
            (pair(self.secondary), "secondary"),
            (pair(self.supplementary), "supplementary"),
            (pair(self.duplicates), "duplicates"),
            (pair(self.primary_duplicates), "primary duplicates"),
            (pair(self.mapped), "mapped"),
            (percent_pair(self.mapped, self.total), "mapped %"),
            (pair(self.primary_mapped), "primary mapped"),
            (
                percent_pair(self.primary_mapped, self.primary),
                "primary mapped %",
            ),
            (pair(self.paired), "paired in sequencing"),
            (pair(self.read1), "read1"),
            (pair(self.read2), "read2"),
            (pair(self.properly_paired), "properly paired"),
            (
                percent_pair(self.properly_paired, self.paired),
                "properly paired %",
            ),
            (pair(self.both_mapped), "with itself and mate mapped"),
            (pair(self.singletons), "singletons"),
            (percent_pair(self.singletons, self.paired), "singletons %"),
            (
                pair(self.mate_different_reference),
                "with mate mapped to a different chr",
            ),
            (
                pair(self.mate_different_reference_mapq5),
                "with mate mapped to a different chr (mapQ>=5)",
            ),
        ] {
            output.push_str(&values[0]);
            output.push('\t');
            output.push_str(&values[1]);
            output.push('\t');
            output.push_str(label);
            output.push('\n');
        }
        output
    }

    pub fn to_json(&self) -> serde_json::Value {
        serde_json::json!({
            "QC-passed reads": category_json(self, 0),
            "QC-failed reads": category_json(self, 1),
        })
    }
}

impl fmt::Display for Counts {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            formatter,
            "{} + {} in total (QC-passed reads + QC-failed reads)",
            self.total[0], self.total[1]
        )?;
        writeln!(
            formatter,
            "{} + {} primary",
            self.primary[0], self.primary[1]
        )?;
        writeln!(
            formatter,
            "{} + {} secondary",
            self.secondary[0], self.secondary[1]
        )?;
        writeln!(
            formatter,
            "{} + {} supplementary",
            self.supplementary[0], self.supplementary[1]
        )?;
        writeln!(
            formatter,
            "{} + {} duplicates",
            self.duplicates[0], self.duplicates[1]
        )?;
        writeln!(
            formatter,
            "{} + {} primary duplicates",
            self.primary_duplicates[0], self.primary_duplicates[1]
        )?;
        metric_line(formatter, "mapped", self.mapped, Some(self.total))?;
        metric_line(
            formatter,
            "primary mapped",
            self.primary_mapped,
            Some(self.primary),
        )?;
        metric_line(formatter, "paired in sequencing", self.paired, None)?;
        metric_line(formatter, "read1", self.read1, None)?;
        metric_line(formatter, "read2", self.read2, None)?;
        metric_line(
            formatter,
            "properly paired",
            self.properly_paired,
            Some(self.paired),
        )?;
        metric_line(
            formatter,
            "with itself and mate mapped",
            self.both_mapped,
            None,
        )?;
        metric_line(formatter, "singletons", self.singletons, Some(self.paired))?;
        metric_line(
            formatter,
            "with mate mapped to a different chr",
            self.mate_different_reference,
            None,
        )?;
        metric_line(
            formatter,
            "with mate mapped to a different chr (mapQ>=5)",
            self.mate_different_reference_mapq5,
            None,
        )
    }
}

fn metric_line(
    formatter: &mut fmt::Formatter<'_>,
    label: &str,
    values: [u64; 2],
    denominator: Option<[u64; 2]>,
) -> fmt::Result {
    if let Some(denominator) = denominator {
        writeln!(
            formatter,
            "{} + {} {} ({} : {})",
            values[0],
            values[1],
            label,
            percent(values[0], denominator[0]),
            percent(values[1], denominator[1])
        )
    } else {
        writeln!(formatter, "{} + {} {label}", values[0], values[1])
    }
}

fn category_json(counts: &Counts, index: usize) -> serde_json::Value {
    serde_json::json!({
        "total": counts.total[index],
        "primary": counts.primary[index],
        "secondary": counts.secondary[index],
        "supplementary": counts.supplementary[index],
        "duplicates": counts.duplicates[index],
        "primary duplicates": counts.primary_duplicates[index],
        "mapped": counts.mapped[index],
        "mapped %": percentage(counts.mapped[index], counts.total[index]),
        "primary mapped": counts.primary_mapped[index],
        "primary mapped %": percentage(counts.primary_mapped[index], counts.primary[index]),
        "paired in sequencing": counts.paired[index],
        "read1": counts.read1[index],
        "read2": counts.read2[index],
        "properly paired": counts.properly_paired[index],
        "properly paired %": percentage(counts.properly_paired[index], counts.paired[index]),
        "with itself and mate mapped": counts.both_mapped[index],
        "singletons": counts.singletons[index],
        "singletons %": percentage(counts.singletons[index], counts.paired[index]),
        "with mate mapped to a different chr": counts.mate_different_reference[index],
        "with mate mapped to a different chr (mapQ >= 5)": counts.mate_different_reference_mapq5[index],
    })
}

fn pair(values: [u64; 2]) -> [String; 2] {
    [values[0].to_string(), values[1].to_string()]
}

fn percent_pair(numerator: [u64; 2], denominator: [u64; 2]) -> [String; 2] {
    [
        percent(numerator[0], denominator[0]),
        percent(numerator[1], denominator[1]),
    ]
}

fn percent(numerator: u64, denominator: u64) -> String {
    percentage(numerator, denominator)
        .map_or_else(|| "N/A".to_owned(), |value| format!("{value:.2}%"))
}

#[allow(clippy::cast_precision_loss)]
fn percentage(numerator: u64, denominator: u64) -> Option<f64> {
    (denominator != 0).then(|| {
        let value = numerator as f64 / denominator as f64 * 100.0;
        (value * 100.0).round() / 100.0
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(flags: u16, tid: i32, mtid: i32, mapq: u8) -> Record {
        let mut record = Record::new();
        record.set_flags(flags);
        record.set_tid(tid);
        record.set_mtid(mtid);
        record.set_mapq(mapq);
        record
    }

    #[test]
    fn secondary_wins_over_supplementary() {
        let mut counts = Counts::default();
        counts.tally(&record(SECONDARY | SUPPLEMENTARY, -1, -1, 0));
        assert_eq!(counts.secondary, [1, 0]);
        assert_eq!(counts.supplementary, [0, 0]);
        assert_eq!(counts.primary, [0, 0]);
    }

    #[test]
    fn qc_categories_and_primary_counts_are_independent() {
        let mut counts = Counts::default();
        counts.tally(&record(0, 0, -1, 60));
        counts.tally(&record(QCFAIL | DUPLICATE, 0, -1, 60));
        assert_eq!(counts.total, [1, 1]);
        assert_eq!(counts.primary, [1, 1]);
        assert_eq!(counts.primary_duplicates, [0, 1]);
    }

    #[test]
    fn paired_metrics_follow_samtools_branching() {
        let mut counts = Counts::default();
        counts.tally(&record(PAIRED | PROPER_PAIR | READ1, 0, 0, 60));
        counts.tally(&record(PAIRED | READ2 | MATE_UNMAPPED, 1, -1, 4));
        assert_eq!(counts.paired, [2, 0]);
        assert_eq!(counts.properly_paired, [1, 0]);
        assert_eq!(counts.both_mapped, [1, 0]);
        assert_eq!(counts.singletons, [1, 0]);
    }

    #[test]
    fn zero_denominator_uses_samtools_na_marker() {
        assert_eq!(percent(0, 0), "N/A");
        assert_eq!(percentage(0, 0), None);
    }

    #[test]
    fn percentages_are_rounded_like_samtools_json() {
        assert_eq!(percentage(2, 3), Some(66.67));
    }
}
