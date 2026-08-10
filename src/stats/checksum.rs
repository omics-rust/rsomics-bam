use serde::Serialize;

use super::record_data::RecordData;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize)]
pub(crate) struct Checksums {
    pub(crate) names: u32,
    pub(crate) sequences: u32,
    pub(crate) qualities: u32,
}

impl Checksums {
    pub(crate) fn update(&mut self, record: &RecordData) {
        self.names = self.names.wrapping_add(crc32(&record.name));
        self.sequences = self.sequences.wrapping_add(crc32(&record.packed_sequence));
        let qualities = &record.qualities;
        self.qualities = self.qualities.wrapping_add(if qualities.is_empty() {
            crc32_repeated(0xff, record.sequence.len())
        } else {
            crc32(qualities)
        });
    }
}

fn crc32(bytes: &[u8]) -> u32 {
    crc32fast::hash(bytes)
}

fn crc32_repeated(byte: u8, length: usize) -> u32 {
    let block = [byte; 256];
    let mut crc = crc32fast::Hasher::new();
    for _ in 0..length / block.len() {
        crc.update(&block);
    }
    crc.update(&block[..length % block.len()]);
    crc.finalize()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crc_matches_standard_vector() {
        assert_eq!(crc32(b"123456789"), 0xcbf4_3926);
    }
}
