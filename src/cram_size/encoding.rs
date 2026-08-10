use std::collections::BTreeMap;
use std::io::{self, Cursor, Read};

use super::varint::{read_itf8, read_nonnegative_itf8, read_u8};
use super::{ContainerEncodings, EncodingSummary, invalid};

const DATA_SERIES_ORDER: [&str; 32] = [
    "RN", "QS", "IN", "SC", "BF", "CF", "AP", "RG", "MQ", "NS", "MF", "TS", "NP", "NF", "RL", "FN",
    "FC", "FP", "DL", "BA", "BS", "TL", "RI", "RS", "PD", "HC", "BB", "QQ", "TN", "TC", "TM", "TV",
];

#[derive(Debug)]
struct Encoding {
    description: String,
    content_ids: Vec<i32>,
}

pub(super) fn parse(data: &[u8]) -> io::Result<(ContainerEncodings, Vec<(i32, String)>)> {
    let mut reader = Cursor::new(data);
    skip_map(&mut reader, "preservation map")?;

    let record_payload = read_map(&mut reader, "record encoding map")?;
    let mut records = BTreeMap::new();
    let mut record_reader = Cursor::new(record_payload.as_slice());
    let count = read_nonnegative_itf8(&mut record_reader, "record encoding count")?;
    check_entry_count(count)?;
    for _ in 0..count {
        let mut key = [0; 2];
        record_reader.read_exact(&mut key)?;
        let encoding = parse_framed_encoding(&mut record_reader)?;
        if encoding.description != "NULL" {
            records.insert(String::from_utf8_lossy(&key).into_owned(), encoding);
        }
    }
    require_end(&record_reader, "record encoding map")?;

    let tag_payload = read_map(&mut reader, "tag encoding map")?;
    let mut tags: [Vec<(String, Encoding)>; 32] = std::array::from_fn(|_| Vec::new());
    let mut tag_reader = Cursor::new(tag_payload.as_slice());
    let count = read_nonnegative_itf8(&mut tag_reader, "tag encoding count")?;
    check_entry_count(count)?;
    for _ in 0..count {
        let packed = read_itf8(&mut tag_reader)? as u32;
        let bytes = [(packed >> 16) as u8, (packed >> 8) as u8, packed as u8];
        let key = if bytes[0] == 0 {
            String::from_utf8_lossy(&bytes[1..]).into_owned()
        } else {
            String::from_utf8_lossy(&bytes).into_owned()
        };
        let encoding = parse_framed_encoding(&mut tag_reader)?;
        let bucket = (usize::from(bytes[0]) * 3 + usize::from(bytes[1])) & 31;
        tags[bucket].push((key, encoding));
    }
    require_end(&tag_reader, "tag encoding map")?;
    require_end(&reader, "compression header")?;

    let mut ordered = Vec::new();
    for key in DATA_SERIES_ORDER {
        if let Some(encoding) = records.remove(key) {
            ordered.push((key.to_owned(), encoding));
        }
    }
    if !records.is_empty() {
        return Err(invalid("compression header contains unknown data series"));
    }
    for bucket in tags {
        ordered.extend(bucket.into_iter().rev());
    }

    let mut mappings = Vec::new();
    let entries = ordered
        .into_iter()
        .map(|(data_series, encoding)| {
            for content_id in encoding.content_ids {
                mappings.push((content_id, data_series.clone()));
            }
            EncodingSummary {
                data_series,
                encoding: encoding.description,
            }
        })
        .collect();
    Ok((ContainerEncodings { entries }, mappings))
}

fn skip_map(reader: &mut Cursor<&[u8]>, field: &str) -> io::Result<()> {
    let _ = read_map(reader, field)?;
    Ok(())
}

fn read_map(reader: &mut Cursor<&[u8]>, field: &str) -> io::Result<Vec<u8>> {
    let size = read_nonnegative_itf8(reader, field)?;
    if size > 64 * 1024 * 1024 {
        return Err(invalid(format!("oversized {field}: {size}")));
    }
    let mut payload = vec![0; size];
    reader.read_exact(&mut payload)?;
    Ok(payload)
}

fn parse_framed_encoding(reader: &mut Cursor<&[u8]>) -> io::Result<Encoding> {
    let kind = read_itf8(reader)?;
    let size = read_nonnegative_itf8(reader, "encoding size")?;
    if size > 64 * 1024 * 1024 {
        return Err(invalid(format!("oversized encoding: {size}")));
    }
    let mut payload = vec![0; size];
    reader.read_exact(&mut payload)?;
    parse_encoding(kind, &payload)
}

fn parse_encoding(kind: i32, data: &[u8]) -> io::Result<Encoding> {
    let mut reader = Cursor::new(data);
    let (description, content_ids) = match kind {
        0 => ("NULL".to_owned(), Vec::new()),
        1 => {
            let id = read_itf8(&mut reader)?;
            (format!("EXTERNAL(id={id})"), vec![id])
        }
        2 | 8 => {
            let _offset = read_itf8(&mut reader)?;
            let _parameter = read_itf8(&mut reader)?;
            ("?".to_owned(), vec![-1])
        }
        3 => {
            let count = read_nonnegative_itf8(&mut reader, "Huffman symbol count")?;
            check_entry_count(count)?;
            let mut codes = Vec::with_capacity(count);
            for _ in 0..count {
                codes.push(read_itf8(&mut reader)?);
            }
            let length_count = read_nonnegative_itf8(&mut reader, "Huffman length count")?;
            if length_count != count {
                return Err(invalid("Huffman symbol and length counts differ"));
            }
            let mut lengths = Vec::with_capacity(count);
            for _ in 0..count {
                lengths.push(read_nonnegative_itf8(&mut reader, "Huffman code length")?);
            }
            let codes = codes
                .iter()
                .map(i32::to_string)
                .collect::<Vec<_>>()
                .join(",");
            let lengths = lengths
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",");
            (
                format!("HUFFMAN(codes={{{codes}}},lengths={{{lengths}}})"),
                if count == 1 { Vec::new() } else { vec![-1] },
            )
        }
        4 => {
            let len = parse_framed_encoding(&mut reader)?;
            let value = parse_framed_encoding(&mut reader)?;
            let mut ids = len.content_ids;
            for id in value.content_ids {
                if !ids.contains(&id) {
                    ids.push(id);
                }
            }
            (
                format!(
                    "BYTE_ARRAY_LEN(len_codec={{{}}},val_codec={{{}}}",
                    len.description, value.description
                ),
                ids,
            )
        }
        5 => {
            let stop = read_u8(&mut reader)?;
            let id = read_itf8(&mut reader)?;
            (format!("BYTE_ARRAY_STOP(stop={stop},id={id})"), vec![id])
        }
        6 => {
            let offset = read_itf8(&mut reader)?;
            let bits = read_nonnegative_itf8(&mut reader, "beta bit count")?;
            (format!("BETA(offset={offset}, nbits={bits})"), vec![-1])
        }
        7 => {
            let offset = read_itf8(&mut reader)?;
            let k = read_nonnegative_itf8(&mut reader, "subexponential parameter")?;
            (format!("SUBEXP(offset={offset},k={k})"), vec![-1])
        }
        9 => {
            let offset = read_itf8(&mut reader)?;
            (format!("GAMMA(offset={offset})"), vec![-1])
        }
        _ => return Err(invalid(format!("unknown CRAM encoding {kind}"))),
    };
    require_end(&reader, "encoding parameters")?;
    Ok(Encoding {
        description,
        content_ids,
    })
}

fn check_entry_count(count: usize) -> io::Result<()> {
    if count > 100_000 {
        Err(invalid(format!("oversized encoding map: {count} entries")))
    } else {
        Ok(())
    }
}

fn require_end(reader: &Cursor<impl AsRef<[u8]>>, field: &str) -> io::Result<()> {
    if reader.position() == reader.get_ref().as_ref().len() as u64 {
        Ok(())
    } else {
        Err(invalid(format!("trailing bytes in {field}")))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn describes_core_and_external_encoding_families() {
        for (kind, data, description, ids) in [
            (1, &[11][..], "EXTERNAL(id=11)", &[11][..]),
            (2, &[0, 3][..], "?", &[-1][..]),
            (6, &[2, 5][..], "BETA(offset=2, nbits=5)", &[-1][..]),
            (7, &[2, 5][..], "SUBEXP(offset=2,k=5)", &[-1][..]),
            (8, &[0, 3][..], "?", &[-1][..]),
            (9, &[2][..], "GAMMA(offset=2)", &[-1][..]),
        ] {
            let encoding = parse_encoding(kind, data).unwrap();
            assert_eq!(encoding.description, description);
            assert_eq!(encoding.content_ids, ids);
        }
    }

    #[test]
    fn rejects_trailing_encoding_parameters() {
        let error = parse_encoding(1, &[11, 0]).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    }
}
