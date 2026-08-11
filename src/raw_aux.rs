use rsomics_bamio::raw::{RawRecord, RecordRef};

pub(crate) enum Integer {
    Missing,
    Value(i64),
    Invalid,
}

pub(crate) trait Fields {
    fn value(&self, tag: [u8; 2]) -> Option<&[u8]>;
    fn type_code(&self, tag: [u8; 2]) -> Option<u8>;
}

impl Fields for RawRecord {
    fn value(&self, tag: [u8; 2]) -> Option<&[u8]> {
        self.aux_value(tag)
    }

    fn type_code(&self, tag: [u8; 2]) -> Option<u8> {
        self.aux_type(tag)
    }
}

impl Fields for RecordRef<'_> {
    fn value(&self, tag: [u8; 2]) -> Option<&[u8]> {
        self.aux_value(tag)
    }

    fn type_code(&self, tag: [u8; 2]) -> Option<u8> {
        self.aux_type(tag)
    }
}

pub(crate) fn integer(record: &impl Fields, tag: [u8; 2]) -> Integer {
    let Some(value) = record.value(tag) else {
        return Integer::Missing;
    };
    let value = match record.type_code(tag) {
        Some(b'c') if value.len() == 1 => i64::from(i8::from_le_bytes([value[0]])),
        Some(b'C') if value.len() == 1 => i64::from(value[0]),
        Some(b's') if value.len() == 2 => i64::from(i16::from_le_bytes(value.try_into().unwrap())),
        Some(b'S') if value.len() == 2 => i64::from(u16::from_le_bytes(value.try_into().unwrap())),
        Some(b'i') if value.len() == 4 => i64::from(i32::from_le_bytes(value.try_into().unwrap())),
        Some(b'I') if value.len() == 4 => i64::from(u32::from_le_bytes(value.try_into().unwrap())),
        _ => return Integer::Invalid,
    };
    Integer::Value(value)
}
