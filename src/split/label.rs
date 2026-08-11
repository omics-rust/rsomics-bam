use rsomics_common::{Result, RsomicsError};

pub(super) fn encode(value: &[u8]) -> Result<String> {
    if value.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "split label cannot be empty".to_owned(),
        ));
    }

    const HEX: &[u8; 16] = b"0123456789ABCDEF";
    let mut output = String::with_capacity(value.len());
    for &byte in value {
        if byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'-' | b'_') {
            output.push(char::from(byte));
        } else {
            output.push('%');
            output.push(char::from(HEX[usize::from(byte >> 4)]));
            output.push(char::from(HEX[usize::from(byte & 0x0f)]));
        }
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_components_are_reversible_and_path_safe() {
        for (value, expected) in [
            (b"rg1".as_slice(), "rg1"),
            (b"a/b".as_slice(), "a%2Fb"),
            (b"a%b".as_slice(), "a%25b"),
            (b"a\\b".as_slice(), "a%5Cb"),
            ([0xff].as_slice(), "%FF"),
        ] {
            assert_eq!(encode(value).unwrap(), expected);
        }
        assert!(encode(b"").is_err());
    }
}
