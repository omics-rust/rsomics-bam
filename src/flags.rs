use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

const DEFINITIONS: [FlagDefinition; 12] = [
    FlagDefinition::new(
        0x1,
        "PAIRED",
        "paired-end / multiple-segment sequencing technology",
    ),
    FlagDefinition::new(
        0x2,
        "PROPER_PAIR",
        "each segment properly aligned according to aligner",
    ),
    FlagDefinition::new(0x4, "UNMAP", "segment unmapped"),
    FlagDefinition::new(0x8, "MUNMAP", "next segment in the template unmapped"),
    FlagDefinition::new(0x10, "REVERSE", "SEQ is reverse complemented"),
    FlagDefinition::new(
        0x20,
        "MREVERSE",
        "SEQ of next segment in template is rev.complemented",
    ),
    FlagDefinition::new(0x40, "READ1", "the first segment in the template"),
    FlagDefinition::new(0x80, "READ2", "the last segment in the template"),
    FlagDefinition::new(0x100, "SECONDARY", "secondary alignment"),
    FlagDefinition::new(
        0x200,
        "QCFAIL",
        "not passing quality controls or other filters",
    ),
    FlagDefinition::new(0x400, "DUP", "PCR or optical duplicate"),
    FlagDefinition::new(0x800, "SUPPLEMENTARY", "supplementary alignment"),
];

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub struct FlagDefinition {
    pub bit: u16,
    pub name: &'static str,
    pub description: &'static str,
}

impl FlagDefinition {
    const fn new(bit: u16, name: &'static str, description: &'static str) -> Self {
        Self {
            bit,
            name,
            description,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct FlagValue {
    pub hexadecimal: String,
    pub decimal: u16,
    pub names: Vec<&'static str>,
}

pub fn definitions() -> &'static [FlagDefinition] {
    &DEFINITIONS
}

pub fn parse(token: &str) -> Result<u16> {
    if let Some(hexadecimal) = token
        .strip_prefix("0x")
        .or_else(|| token.strip_prefix("0X"))
    {
        return u16::from_str_radix(hexadecimal, 16).map_err(|_| invalid_flag(token));
    }

    if token.len() > 1 && token.starts_with('0') && token.bytes().all(|byte| byte.is_ascii_digit())
    {
        return u16::from_str_radix(&token[1..], 8).map_err(|_| invalid_flag(token));
    }

    if !token.is_empty() && token.bytes().all(|byte| byte.is_ascii_digit()) {
        return token.parse::<u16>().map_err(|_| invalid_flag(token));
    }

    let mut value = 0;
    for name in token.split(',') {
        let definition = DEFINITIONS
            .iter()
            .find(|definition| definition.name.eq_ignore_ascii_case(name))
            .ok_or_else(|| invalid_flag(token))?;
        value |= definition.bit;
    }
    Ok(value)
}

pub fn describe(value: u16) -> FlagValue {
    let names = DEFINITIONS
        .iter()
        .filter(|definition| value & definition.bit != 0)
        .map(|definition| definition.name)
        .collect();

    FlagValue {
        hexadecimal: format!("0x{value:x}"),
        decimal: value,
        names,
    }
}

pub fn render(value: &FlagValue) -> String {
    format!(
        "{}\t{}\t{}",
        value.hexadecimal,
        value.decimal,
        value.names.join(",")
    )
}

fn invalid_flag(token: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!("could not parse FLAG {token:?}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numeric_forms_match_strtol_base_zero() {
        assert_eq!(parse("16").unwrap(), 16);
        assert_eq!(parse("020").unwrap(), 16);
        assert_eq!(parse("0x10").unwrap(), 16);
        assert_eq!(parse("0X10").unwrap(), 16);
    }

    #[test]
    fn names_are_case_insensitive() {
        assert_eq!(parse("paired,READ1").unwrap(), 65);
        assert_eq!(parse("proper_pair").unwrap(), 2);
    }

    #[test]
    fn malformed_values_fail() {
        for value in ["", "nonsense", "08", "0x10000", "paired,unknown"] {
            assert!(parse(value).is_err(), "{value}");
        }
    }

    #[test]
    fn rendering_matches_samtools_shape() {
        assert_eq!(render(&describe(16)), "0x10\t16\tREVERSE");
        assert_eq!(render(&describe(3)), "0x3\t3\tPAIRED,PROPER_PAIR");
        assert_eq!(render(&describe(0)), "0x0\t0\t");
    }
}
