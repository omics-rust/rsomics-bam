pub(crate) mod addreplacerg;
pub(crate) mod cat;
pub(crate) mod collate;
pub(crate) mod depth;
pub(crate) mod fastx;
pub(crate) mod fixmate;
pub(crate) mod flags;
pub(crate) mod flagstat;
pub(crate) mod head;
pub(crate) mod import;
pub(crate) mod index;
pub(crate) mod markdup;
pub(crate) mod merge;
pub(crate) mod mpileup;
pub(crate) mod quickcheck;
pub(crate) mod reheader;
pub(crate) mod samples;
pub(crate) mod sort;
pub(crate) mod view;

pub(crate) fn parse_memory(value: &str) -> std::result::Result<u64, String> {
    let (number, multiplier) = match value.as_bytes().last().copied() {
        Some(b'k' | b'K') => (&value[..value.len() - 1], 1u64 << 10),
        Some(b'm' | b'M') => (&value[..value.len() - 1], 1u64 << 20),
        Some(b'g' | b'G') => (&value[..value.len() - 1], 1u64 << 30),
        _ => (value, 1),
    };
    let bytes = number
        .parse::<u64>()
        .map_err(|_| format!("invalid memory size: {value}"))?
        .checked_mul(multiplier)
        .ok_or_else(|| format!("memory size overflows: {value}"))?;
    if bytes < 1 << 20 {
        return Err("memory must be at least 1M".to_owned());
    }
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_memory_units() {
        assert_eq!(parse_memory("1M"), Ok(1 << 20));
        assert_eq!(parse_memory("2G"), Ok(2 << 30));
        assert!(parse_memory("1023K").is_err());
        assert!(parse_memory("1.5G").is_err());
    }
}
