use std::fmt;
use std::io::{self, Write};

use super::{Qc, Report};

impl Report {
    pub fn write(&self, mut output: impl Write) -> io::Result<()> {
        if self.layout.bamseqchksum {
            self.write_bamseqchksum(&mut output)
        } else if self.layout.tabs {
            self.write_tabs(&mut output)
        } else {
            self.write_aligned(&mut output)
        }
    }

    fn write_aligned(&self, output: &mut impl Write) -> io::Result<()> {
        writeln!(output, "# Checksum 1.0 for file: {}", self.source)?;
        writeln!(output, "# Aux tags:          {}", self.auxiliary_tags)?;
        writeln!(
            output,
            "# BAM flags:         {}",
            flag_names(self.flag_mask)
        )?;
        write!(
            output,
            "\n# Group    QC          count  flag+seq  +name     +qual     +aux    "
        )?;
        if self.layout.position {
            write!(output, "  +chr/pos")?;
        }
        if self.layout.cigar {
            write!(output, "  +cigar  ")?;
        }
        if self.layout.mate {
            write!(output, "  +mate   ")?;
        }
        writeln!(output, "  combined")?;
        for group in &self.groups {
            for row in &group.rows {
                write!(
                    output,
                    "{:<10} {:<4} {:>12}  {:08x}  {:08x}  {:08x}  {:08x}",
                    group.name,
                    row.qc.label(),
                    row.count,
                    row.checksums.sequence,
                    row.checksums.name,
                    row.checksums.quality,
                    row.checksums.auxiliary
                )?;
                for value in [
                    row.checksums.position,
                    row.checksums.cigar,
                    row.checksums.mate,
                ]
                .into_iter()
                .flatten()
                {
                    write!(output, "  {value:08x}")?;
                }
                writeln!(output, "  {:08x}", row.checksums.combined)?;
            }
        }
        Ok(())
    }

    fn write_tabs(&self, output: &mut impl Write) -> io::Result<()> {
        writeln!(output, "# Checksum 1.0 for file:\t{}", self.source)?;
        writeln!(output, "# Aux tags:\t{}", self.auxiliary_tags)?;
        writeln!(output, "# BAM flags:\t{}", flag_names(self.flag_mask))?;
        write!(output, "\n# Group\tQC\tcount\tflag+seq\t+name\t+qual\t+aux")?;
        if self.layout.position {
            write!(output, "\t+chr/pos")?;
        }
        if self.layout.cigar {
            write!(output, "\t+cigar")?;
        }
        if self.layout.mate {
            write!(output, "\t+mate")?;
        }
        writeln!(output, "\tcombined")?;
        for group in &self.groups {
            for row in &group.rows {
                write!(
                    output,
                    "{}\t{}\t{}\t{:x}\t{:x}\t{:x}\t{:x}",
                    group.name,
                    row.qc.label(),
                    row.count,
                    row.checksums.sequence,
                    row.checksums.name,
                    row.checksums.quality,
                    row.checksums.auxiliary
                )?;
                for value in [
                    row.checksums.position,
                    row.checksums.cigar,
                    row.checksums.mate,
                ]
                .into_iter()
                .flatten()
                {
                    write!(output, "\t{value:x}")?;
                }
                writeln!(output, "\t{:x}", row.checksums.combined)?;
            }
        }
        Ok(())
    }

    fn write_bamseqchksum(&self, output: &mut impl Write) -> io::Result<()> {
        writeln!(
            output,
            "###\tset\tcount\t\tb_seq\tname_b_seq\tb_seq_qual\tb_seq_tags(BC,FI,QT,RT,TC)"
        )?;
        for group in &self.groups {
            let name = if group.name == "-" { "" } else { &group.name };
            for row in &group.rows {
                if row.qc == Qc::Fail {
                    continue;
                }
                writeln!(
                    output,
                    "{}\t{}\t{}\t\t{:x}\t{:x}\t{:x}\t{:x}",
                    name,
                    row.qc.label(),
                    row.count,
                    row.checksums.sequence,
                    row.checksums.name,
                    row.checksums.quality,
                    row.checksums.auxiliary
                )?;
            }
        }
        Ok(())
    }
}

fn flag_names(mask: u16) -> impl fmt::Display {
    struct Names(u16);
    impl fmt::Display for Names {
        fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
            let names = crate::flags::describe(self.0).names;
            if names.is_empty() {
                formatter.write_str("0")
            } else {
                formatter.write_str(&names.join(","))
            }
        }
    }
    Names(mask)
}
