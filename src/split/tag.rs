use noodles::sam;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};

#[derive(Debug, Eq, PartialEq)]
pub(super) enum Outcome {
    Present(Vec<u8>),
    Missing,
    Invalid,
}

pub(super) fn read(record: &RawRecord, tag: [u8; 2], width: usize) -> Outcome {
    let Some(value) = record.aux_value(tag) else {
        return Outcome::Missing;
    };
    let type_code = record.aux_type(tag);
    match type_code {
        Some(b'Z' | b'H') => value
            .strip_suffix(&[0])
            .map_or(Outcome::Invalid, |value| Outcome::Present(value.to_vec())),
        Some(b'c') if value.len() == 1 => integer(i64::from(i8::from_le_bytes([value[0]])), width),
        Some(b'C') if value.len() == 1 => integer(i64::from(value[0]), width),
        Some(b's') if value.len() == 2 => integer(
            i64::from(i16::from_le_bytes(value.try_into().unwrap())),
            width,
        ),
        Some(b'S') if value.len() == 2 => integer(
            i64::from(u16::from_le_bytes(value.try_into().unwrap())),
            width,
        ),
        Some(b'i') if value.len() == 4 => integer(
            i64::from(i32::from_le_bytes(value.try_into().unwrap())),
            width,
        ),
        Some(b'I') if value.len() == 4 => integer(
            i64::from(u32::from_le_bytes(value.try_into().unwrap())),
            width,
        ),
        _ => Outcome::Invalid,
    }
}

pub(super) fn read_string(record: &RawRecord, tag: [u8; 2]) -> Outcome {
    let Some(value) = record.aux_value(tag) else {
        return Outcome::Missing;
    };
    match record.aux_type(tag) {
        Some(b'Z') => value
            .strip_suffix(&[0])
            .map_or(Outcome::Invalid, |value| Outcome::Present(value.to_vec())),
        _ => Outcome::Invalid,
    }
}

pub(super) fn read_record(
    record: &dyn sam::alignment::Record,
    tag: [u8; 2],
    width: usize,
) -> Result<Outcome> {
    let data = record.data();
    let Some(value) = data.get(&Tag::from(tag)) else {
        return Ok(Outcome::Missing);
    };
    let value = value.map_err(RsomicsError::Io)?;
    Ok(match value {
        Value::String(value) | Value::Hex(value) => Outcome::Present(value.to_vec()),
        value => value
            .as_int()
            .map_or(Outcome::Invalid, |value| integer(value, width)),
    })
}

pub(super) fn read_string_record(
    record: &dyn sam::alignment::Record,
    tag: [u8; 2],
) -> Result<Outcome> {
    let data = record.data();
    let Some(value) = data.get(&Tag::from(tag)) else {
        return Ok(Outcome::Missing);
    };
    let value = value.map_err(RsomicsError::Io)?;
    Ok(match value {
        Value::String(value) => Outcome::Present(value.to_vec()),
        _ => Outcome::Invalid,
    })
}

fn integer(value: i64, width: usize) -> Outcome {
    let value = if width == 0 {
        value.to_string()
    } else if value < 0 {
        format!("-{value:0width$}", value = value.unsigned_abs())
    } else {
        format!("{value:0width$}")
    };
    Outcome::Present(value.into_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(type_code: u8, value: &[u8]) -> RawRecord {
        let mut record = RawRecord::default();
        record.set_aux(*b"ZZ", type_code, value).unwrap();
        record
    }

    #[test]
    fn tag_values_preserve_strings_hexadecimal_and_integer_ranges() {
        for (type_code, value, width, expected) in [
            (b'Z', b"group\0".as_slice(), 0, b"group".as_slice()),
            (b'H', b"0A7f\0".as_slice(), 0, b"0A7f".as_slice()),
            (
                b'c',
                (-5_i8).to_le_bytes().as_slice(),
                3,
                b"-005".as_slice(),
            ),
            (b'C', 7_u8.to_le_bytes().as_slice(), 3, b"007".as_slice()),
            (
                b's',
                i16::MIN.to_le_bytes().as_slice(),
                0,
                b"-32768".as_slice(),
            ),
            (
                b'S',
                u16::MAX.to_le_bytes().as_slice(),
                0,
                b"65535".as_slice(),
            ),
            (
                b'i',
                i32::MIN.to_le_bytes().as_slice(),
                0,
                b"-2147483648".as_slice(),
            ),
            (
                b'I',
                u32::MAX.to_le_bytes().as_slice(),
                0,
                b"4294967295".as_slice(),
            ),
        ] {
            assert_eq!(
                read(&record(type_code, value), *b"ZZ", width),
                Outcome::Present(expected.to_vec())
            );
        }
    }

    #[test]
    fn missing_and_non_groupable_types_are_unaccounted() {
        assert_eq!(read(&RawRecord::default(), *b"ZZ", 0), Outcome::Missing);
        for (type_code, value) in [
            (b'A', b"x".as_slice()),
            (b'f', 1.5_f32.to_le_bytes().as_slice()),
            (b'B', [b'C', 1, 0, 0, 0, 9].as_slice()),
        ] {
            assert_eq!(read(&record(type_code, value), *b"ZZ", 0), Outcome::Invalid);
        }
    }
}
