use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

use super::{Accumulator, Options, Qc, RecordChecksums, Report, Sums, TagSelection, fold};

#[derive(Clone, Debug, Eq, PartialEq)]
struct Columns {
    position: bool,
    cigar: bool,
    mate: bool,
}

#[derive(Default)]
struct Contract {
    tags: Option<TagSelection>,
    flag_mask: Option<u16>,
    columns: Option<Columns>,
}

pub fn merge(paths: &[impl AsRef<Path>], options: &Options) -> Result<Report> {
    let mut accumulator = Accumulator::default();
    let mut contract = Contract::default();
    for path in paths {
        parse(path.as_ref(), &mut contract, &mut accumulator)?;
    }
    let mut merged_options = options.clone();
    if let Some(tags) = contract.tags {
        merged_options.tags = tags;
    }
    if let Some(flag_mask) = contract.flag_mask {
        merged_options.flag_mask = flag_mask;
    }
    if let Some(columns) = contract.columns {
        merged_options.check_position = columns.position;
        merged_options.check_cigar = columns.cigar;
        merged_options.check_mate = columns.mate;
    }
    Ok(accumulator.report("merge".to_owned(), &merged_options))
}

fn parse(path: &Path, contract: &mut Contract, accumulator: &mut Accumulator) -> Result<()> {
    let input = File::open(path).map_err(|error| contextual_io(path, error))?;
    let mut format = None;
    let mut columns = None;
    let mut saw_version = false;
    let mut saw_tags = false;
    let mut saw_flags = false;
    let mut saw_rows = false;
    let mut rows = BTreeSet::new();
    let mut totals = BTreeMap::new();
    let mut components = Accumulator::default();
    for result in BufReader::new(input).lines() {
        let line = result.map_err(|error| contextual_io(path, error))?;
        if let Some(version) = line.strip_prefix("# Checksum ") {
            if saw_version {
                return invalid(path, "repeated checksum version header");
            }
            let version = version.split_whitespace().next().unwrap_or_default();
            if version != "1.0" {
                return invalid(
                    path,
                    format!("unsupported checksum output version {version:?}"),
                );
            }
            saw_version = true;
        } else if let Some(tags) = line.strip_prefix("# Aux tags:") {
            if saw_tags {
                return invalid(path, "repeated auxiliary-tag header");
            }
            saw_tags = true;
            let tags =
                parse_tag_contract(tags.trim()).map_err(|error| invalid_error(path, error))?;
            merge_value(&mut contract.tags, tags, path, "auxiliary tags")?;
        } else if let Some(flags) = line.strip_prefix("# BAM flags:") {
            if saw_flags {
                return invalid(path, "repeated flag-mask header");
            }
            saw_flags = true;
            let mask =
                crate::flags::parse(flags.trim()).map_err(|error| invalid_error(path, error))?;
            merge_value(&mut contract.flag_mask, mask, path, "flag mask")?;
        } else if line.starts_with("# Group") {
            if format.is_some() {
                return invalid(path, "repeated checksum column header");
            }
            if !saw_version || !saw_tags || !saw_flags {
                return invalid(path, "native checksum metadata is incomplete");
            }
            let parsed = native_columns(&line, path)?;
            merge_value(&mut contract.columns, parsed.clone(), path, "columns")?;
            columns = Some(parsed);
            format = Some(false);
        } else if line.starts_with("###\t") {
            if format.is_some() {
                return invalid(path, "repeated checksum column header");
            }
            if saw_version || saw_tags || saw_flags {
                return invalid(path, "native and bamseqchksum headers cannot be mixed");
            }
            let tags = parse_tag_contract(&bamseq_header(&line, path)?)
                .map_err(|error| invalid_error(path, error))?;
            saw_tags = true;
            merge_value(&mut contract.tags, tags, path, "auxiliary tags")?;
            let parsed = Columns {
                position: false,
                cigar: false,
                mate: false,
            };
            merge_value(&mut contract.columns, parsed.clone(), path, "columns")?;
            columns = Some(parsed);
            format = Some(true);
        } else if line.is_empty() || line.starts_with('#') {
            continue;
        } else {
            let is_bamseq = format.ok_or_else(|| {
                invalid_error(path, "data row appears before a recognized header")
            })?;
            let parsed_columns = columns
                .as_ref()
                .expect("recognized report format has columns");
            let row = if is_bamseq {
                parse_bamseq_row(&line, path)?
            } else {
                parse_native_row(&line, parsed_columns, path)?
            };
            saw_rows = true;
            let row_index = qc_index(row.qc);
            if !rows.insert((row.group.clone(), row_index)) {
                return invalid(path, format!("duplicate checksum row in {line:?}"));
            }
            if row.group == b"all" {
                totals.insert(row_index, (row.checksums, row.count));
            } else {
                components.merge_row(&row.group, row.qc, row.checksums, row.count);
                accumulator.merge_row(&row.group, row.qc, row.checksums, row.count);
            }
        }
    }
    if format.is_none() || !saw_rows {
        return invalid(path, "checksum report has no recognized data rows");
    }
    if !format.unwrap() && !saw_version {
        return invalid(path, "native checksum report has no version header");
    }
    if !saw_tags {
        return invalid(path, "checksum report has no auxiliary-tag header");
    }
    if format == Some(false) && !saw_flags {
        return invalid(path, "native checksum report has no flag-mask header");
    }
    if !totals.contains_key(&0) {
        return invalid(path, "checksum report has no all/all total row");
    }
    for (row, (checksums, count)) in totals {
        validate_total(path, row, checksums, count, &components.all)?;
    }
    Ok(())
}

fn qc_index(qc: Qc) -> u8 {
    match qc {
        Qc::All => 0,
        Qc::Pass => 1,
        Qc::Fail => 2,
    }
}

fn validate_total(
    path: &Path,
    row: u8,
    checksums: RecordChecksums,
    count: u64,
    sums: &Sums,
) -> Result<()> {
    let index = usize::from(row);
    let actual = RecordChecksums {
        sequence: sums.sequence[index] as u32,
        name: sums.name[index] as u32,
        quality: sums.quality[index] as u32,
        auxiliary: sums.auxiliary[index] as u32,
        position: sums.position[index] as u32,
        cigar: sums.cigar[index] as u32,
        mate: sums.mate[index] as u32,
    };
    if sums.count[index] != count || actual != checksums {
        return invalid(path, "all-group totals differ from component groups");
    }
    Ok(())
}

struct ParsedRow {
    group: Vec<u8>,
    qc: Qc,
    count: u64,
    checksums: RecordChecksums,
}

fn parse_native_row(line: &str, columns: &Columns, path: &Path) -> Result<ParsedRow> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    let expected =
        8 + usize::from(columns.position) + usize::from(columns.cigar) + usize::from(columns.mate);
    if fields.len() != expected {
        return invalid(path, format!("incorrect number of columns in {line:?}"));
    }
    let mut index = 0;
    let group = fields[index].as_bytes().to_vec();
    index += 1;
    let qc = parse_qc(fields[index], path)?;
    index += 1;
    let count = parse_count(fields[index], path)?;
    index += 1;
    let sequence = parse_hex(fields[index], path)?;
    index += 1;
    let name = parse_hex(fields[index], path)?;
    index += 1;
    let quality = parse_hex(fields[index], path)?;
    index += 1;
    let auxiliary = parse_hex(fields[index], path)?;
    index += 1;
    let position = if columns.position {
        let value = parse_hex(fields[index], path)?;
        index += 1;
        value
    } else {
        1
    };
    let cigar = if columns.cigar {
        let value = parse_hex(fields[index], path)?;
        index += 1;
        value
    } else {
        1
    };
    let mate = if columns.mate {
        let value = parse_hex(fields[index], path)?;
        index += 1;
        value
    } else {
        1
    };
    let combined = parse_hex(fields[index], path)?;
    let mut expected_combined = 1;
    for value in [
        count >> 32,
        count & u64::from(u32::MAX),
        u64::from(sequence),
        u64::from(name),
        u64::from(sequence),
        u64::from(auxiliary),
    ] {
        expected_combined = fold(expected_combined, value as u32);
    }
    for value in [
        columns.position.then_some(position),
        columns.cigar.then_some(cigar),
        columns.mate.then_some(mate),
    ]
    .into_iter()
    .flatten()
    {
        expected_combined = fold(expected_combined, value);
    }
    if u64::from(combined) != expected_combined {
        return invalid(path, format!("combined checksum is invalid in {line:?}"));
    }
    Ok(ParsedRow {
        group,
        qc,
        count,
        checksums: RecordChecksums {
            sequence,
            name,
            quality,
            auxiliary,
            position,
            cigar,
            mate,
        },
    })
}

fn parse_bamseq_row(line: &str, path: &Path) -> Result<ParsedRow> {
    let mut fields = line.split('\t');
    let group = fields.next().unwrap_or_default();
    let qc = fields
        .next()
        .ok_or_else(|| invalid_error(path, "missing QC column"))?;
    let count = fields
        .next()
        .ok_or_else(|| invalid_error(path, "missing count column"))?;
    if fields.next() != Some("") {
        return invalid(
            path,
            "bamseqchksum count must be followed by an empty column",
        );
    }
    let checksums = fields.collect::<Vec<_>>();
    if checksums.len() != 4 {
        return invalid(path, "bamseqchksum row must have four checksum columns");
    }
    Ok(ParsedRow {
        group: if group.is_empty() {
            b"-".to_vec()
        } else {
            group.as_bytes().to_vec()
        },
        qc: parse_qc(qc, path)?,
        count: parse_count(count, path)?,
        checksums: RecordChecksums {
            sequence: parse_hex(checksums[0], path)?,
            name: parse_hex(checksums[1], path)?,
            quality: parse_hex(checksums[2], path)?,
            auxiliary: parse_hex(checksums[3], path)?,
            position: 1,
            cigar: 1,
            mate: 1,
        },
    })
}

fn native_columns(line: &str, path: &Path) -> Result<Columns> {
    let fields = line
        .trim_start_matches('#')
        .split_whitespace()
        .collect::<Vec<_>>();
    let required = ["Group", "QC", "count", "flag+seq", "+name", "+qual", "+aux"];
    if fields.get(..required.len()) != Some(required.as_slice())
        || fields.last() != Some(&"combined")
    {
        return invalid(path, format!("unrecognized checksum header {line:?}"));
    }
    let optional = &fields[required.len()..fields.len() - 1];
    for field in optional {
        if !matches!(*field, "+chr/pos" | "+cigar" | "+mate") {
            return invalid(path, format!("unrecognized checksum column {field:?}"));
        }
    }
    let expected = ["+chr/pos", "+cigar", "+mate"]
        .into_iter()
        .filter(|field| optional.contains(field))
        .collect::<Vec<_>>();
    if optional != expected {
        return invalid(path, "optional checksum columns are out of order");
    }
    Ok(Columns {
        position: optional.contains(&"+chr/pos"),
        cigar: optional.contains(&"+cigar"),
        mate: optional.contains(&"+mate"),
    })
}

fn bamseq_header(line: &str, path: &Path) -> Result<String> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    if fields.len() != 7
        || fields[..6] != ["###", "set", "count", "b_seq", "name_b_seq", "b_seq_qual"]
    {
        return invalid(path, format!("unrecognized bamseqchksum header {line:?}"));
    }
    fields[6]
        .strip_prefix("b_seq_tags(")
        .and_then(|value| value.strip_suffix(')'))
        .map(str::to_owned)
        .ok_or_else(|| invalid_error(path, "invalid bamseqchksum auxiliary-tag header"))
}

fn parse_tag_contract(value: &str) -> std::result::Result<TagSelection, String> {
    let wildcard = value == "*" || value.starts_with("*,");
    let values = value
        .split(',')
        .skip(usize::from(wildcard))
        .map(|tag| {
            let bytes: [u8; 2] = tag
                .as_bytes()
                .try_into()
                .map_err(|_| format!("invalid auxiliary tag in report: {tag:?}"))?;
            if !bytes.iter().all(|byte| (b'0'..=b'z').contains(byte)) {
                return Err(format!("invalid auxiliary tag in report: {tag:?}"));
            }
            Ok(bytes)
        })
        .collect::<std::result::Result<Vec<[u8; 2]>, String>>()?;
    if !wildcard && values.is_empty() {
        return Err("checksum report has no auxiliary tags".to_owned());
    }
    let mut unique = values.clone();
    unique.sort_unstable();
    unique.dedup();
    if unique.len() != values.len() {
        return Err("checksum report contains duplicate auxiliary tags".to_owned());
    }
    Ok(if wildcard {
        TagSelection::AllExcept(values)
    } else {
        TagSelection::Listed(values)
    })
}

fn parse_qc(value: &str, path: &Path) -> Result<Qc> {
    match value {
        "all" => Ok(Qc::All),
        "pass" => Ok(Qc::Pass),
        "fail" => Ok(Qc::Fail),
        _ => invalid(path, format!("invalid QC value {value:?}")),
    }
}

fn parse_count(value: &str, path: &Path) -> Result<u64> {
    value
        .parse()
        .map_err(|_| invalid_error(path, format!("invalid record count {value:?}")))
}

fn parse_hex(value: &str, path: &Path) -> Result<u32> {
    u32::from_str_radix(value, 16)
        .map_err(|_| invalid_error(path, format!("invalid checksum value {value:?}")))
}

fn merge_value<T: Eq>(current: &mut Option<T>, value: T, path: &Path, label: &str) -> Result<()> {
    if current.as_ref().is_some_and(|current| current != &value) {
        return invalid(path, format!("{label} differ from earlier reports"));
    }
    current.get_or_insert(value);
    Ok(())
}

fn contextual_io(path: &Path, error: std::io::Error) -> RsomicsError {
    RsomicsError::Io(std::io::Error::new(
        error.kind(),
        format!("reading checksum report {}: {error}", path.display()),
    ))
}

fn invalid<T>(path: &Path, message: impl std::fmt::Display) -> Result<T> {
    Err(invalid_error(path, message))
}

fn invalid_error(path: &Path, message: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!("checksum report {}: {message}", path.display()))
}
