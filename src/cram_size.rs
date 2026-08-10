mod encoding;
mod parser;
mod render;
mod varint;

use std::collections::BTreeMap;
use std::io::{self, BufReader, Read, Write};

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Method {
    Raw,
    Gzip,
    Bzip2,
    Lzma,
    Fqzcomp,
    GzipMin,
    GzipMax,
    Bzip2_1,
    Bzip2_2,
    Bzip2_3,
    Bzip2_4,
    Bzip2_5,
    Bzip2_6,
    Bzip2_7,
    Bzip2_8,
    Bzip2_9,
    Rans4x8O0,
    Rans4x8O1,
    Rans4x16O0,
    Rans4x16O1,
    Rans4x16O0R,
    Rans4x16O1R,
    Rans4x16O0P,
    Rans4x16O1P,
    Rans4x16O0PR,
    Rans4x16O1PR,
    Rans32x16O0,
    Rans32x16O1,
    Rans32x16O0R,
    Rans32x16O1R,
    Rans32x16O0P,
    Rans32x16O1P,
    Rans32x16O0PR,
    Rans32x16O1PR,
    RansNx16Stripe,
    RansNx16Cat,
    ArithmeticO0,
    ArithmeticO1,
    ArithmeticO0R,
    ArithmeticO1R,
    ArithmeticO0P,
    ArithmeticO1P,
    ArithmeticO0PR,
    ArithmeticO1PR,
    ArithmeticStripe,
    ArithmeticCat,
    ArithmeticExternal,
    TokenizerRans,
    TokenizerArithmetic,
}

impl Method {
    fn classify(code: u8, src: &[u8]) -> io::Result<Self> {
        let first = src.first().copied().unwrap_or_default();
        match code {
            0 => Ok(Self::Raw),
            1 => Ok(match src.get(8) {
                Some(4) => Self::GzipMin,
                Some(2) => Self::GzipMax,
                _ => Self::Gzip,
            }),
            2 => Ok(match src.get(3).copied() {
                Some(b'1') => Self::Bzip2_1,
                Some(b'2') => Self::Bzip2_2,
                Some(b'3') => Self::Bzip2_3,
                Some(b'4') => Self::Bzip2_4,
                Some(b'5') => Self::Bzip2_5,
                Some(b'6') => Self::Bzip2_6,
                Some(b'7') => Self::Bzip2_7,
                Some(b'8') => Self::Bzip2_8,
                Some(b'9') => Self::Bzip2_9,
                _ => Self::Bzip2,
            }),
            3 => Ok(Self::Lzma),
            4 => Ok(if first == 1 {
                Self::Rans4x8O1
            } else {
                Self::Rans4x8O0
            }),
            5 => Ok(classify_rans_nx16(first)),
            6 => Ok(classify_arithmetic(first)),
            7 => Ok(Self::Fqzcomp),
            8 => Ok(if src.get(8) == Some(&1) {
                Self::TokenizerArithmetic
            } else {
                Self::TokenizerRans
            }),
            _ => Err(invalid(format!(
                "unknown CRAM block compression method {code}"
            ))),
        }
    }

    fn short(self) -> &'static str {
        match self {
            Self::Raw => ".",
            Self::Gzip => "g",
            Self::Bzip2 => "b",
            Self::Lzma => "l",
            Self::Fqzcomp => "f",
            Self::GzipMin => "_",
            Self::GzipMax => "G",
            Self::Bzip2_1
            | Self::Bzip2_2
            | Self::Bzip2_3
            | Self::Bzip2_4
            | Self::Bzip2_5
            | Self::Bzip2_6
            | Self::Bzip2_7
            | Self::Bzip2_8 => "b",
            Self::Bzip2_9 => "B",
            Self::Rans4x8O0 => "r",
            Self::Rans4x8O1 => "R",
            Self::Rans4x16O0 | Self::Rans4x16O0R | Self::Rans4x16O0P | Self::Rans4x16O0PR => "0",
            Self::Rans4x16O1 | Self::Rans4x16O1R | Self::Rans4x16O1P | Self::Rans4x16O1PR => "1",
            Self::Rans32x16O0 | Self::Rans32x16O0R | Self::Rans32x16O0P | Self::Rans32x16O0PR => {
                "4"
            }
            Self::Rans32x16O1 | Self::Rans32x16O1R | Self::Rans32x16O1P | Self::Rans32x16O1PR => {
                "5"
            }
            Self::RansNx16Stripe => "8",
            Self::RansNx16Cat => "2",
            Self::ArithmeticO0
            | Self::ArithmeticO0R
            | Self::ArithmeticO0P
            | Self::ArithmeticO0PR
            | Self::ArithmeticStripe
            | Self::ArithmeticCat
            | Self::ArithmeticExternal => "a",
            Self::ArithmeticO1
            | Self::ArithmeticO1R
            | Self::ArithmeticO1P
            | Self::ArithmeticO1PR => "A",
            Self::TokenizerRans => "n",
            Self::TokenizerArithmetic => "N",
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Raw => "raw",
            Self::Gzip => "gzip",
            Self::Bzip2 => "bzip2",
            Self::Lzma => "lzma",
            Self::Fqzcomp => "fqzcomp",
            Self::GzipMin => "gzip-min",
            Self::GzipMax => "gzip-max",
            Self::Bzip2_1 => "bzip2-1",
            Self::Bzip2_2 => "bzip2-2",
            Self::Bzip2_3 => "bzip2-3",
            Self::Bzip2_4 => "bzip2-4",
            Self::Bzip2_5 => "bzip2-5",
            Self::Bzip2_6 => "bzip2-6",
            Self::Bzip2_7 => "bzip2-7",
            Self::Bzip2_8 => "bzip2-8",
            Self::Bzip2_9 => "bzip2-9",
            Self::Rans4x8O0 => "r4x8-o0",
            Self::Rans4x8O1 => "r4x8-o1",
            Self::Rans4x16O0 => "r4x16-o0",
            Self::Rans4x16O1 => "r4x16-o1",
            Self::Rans4x16O0R => "r4x16-o0R",
            Self::Rans4x16O1R => "r4x16-o1R",
            Self::Rans4x16O0P => "r4x16-o0P",
            Self::Rans4x16O1P => "r4x16-o1P",
            Self::Rans4x16O0PR => "r4x16-o0PR",
            Self::Rans4x16O1PR => "r4x16-o1PR",
            Self::Rans32x16O0 => "r32x16-o0",
            Self::Rans32x16O1 => "r32x16-o1",
            Self::Rans32x16O0R => "r32x16-o0R",
            Self::Rans32x16O1R => "r32x16-o1R",
            Self::Rans32x16O0P => "r32x16-o0P",
            Self::Rans32x16O1P => "r32x16-o1P",
            Self::Rans32x16O0PR => "r32x16-o0PR",
            Self::Rans32x16O1PR => "r32x16-o1PR",
            Self::RansNx16Stripe => "rNx16-xo0",
            Self::RansNx16Cat => "rNx16-cat",
            Self::ArithmeticO0 => "arith-o0",
            Self::ArithmeticO1 => "arith-o1",
            Self::ArithmeticO0R => "arith-o0R",
            Self::ArithmeticO1R => "arith-o1R",
            Self::ArithmeticO0P => "arith-o0P",
            Self::ArithmeticO1P => "arith-o1P",
            Self::ArithmeticO0PR => "arith-o0PR",
            Self::ArithmeticO1PR => "arith-o1PR",
            Self::ArithmeticStripe => "arith-stripe",
            Self::ArithmeticCat => "arith-cat",
            Self::ArithmeticExternal => "arith-ext",
            Self::TokenizerRans => "tok3-rans",
            Self::TokenizerArithmetic => "tok3-arith",
        }
    }

    fn rank(self) -> u8 {
        match self {
            Self::Raw => 0,
            Self::Gzip => 1,
            Self::Bzip2 => 2,
            Self::Lzma => 3,
            Self::Fqzcomp => 7,
            Self::GzipMin => 9,
            Self::GzipMax => 10,
            Self::Bzip2_1 => 11,
            Self::Bzip2_2 => 12,
            Self::Bzip2_3 => 13,
            Self::Bzip2_4 => 14,
            Self::Bzip2_5 => 15,
            Self::Bzip2_6 => 16,
            Self::Bzip2_7 => 17,
            Self::Bzip2_8 => 18,
            Self::Bzip2_9 => 19,
            Self::Rans4x8O0 => 20,
            Self::Rans4x8O1 => 21,
            Self::Rans4x16O0 => 22,
            Self::Rans4x16O1 => 23,
            Self::Rans4x16O0R => 24,
            Self::Rans4x16O1R => 25,
            Self::Rans4x16O0P => 26,
            Self::Rans4x16O1P => 27,
            Self::Rans4x16O0PR => 28,
            Self::Rans4x16O1PR => 29,
            Self::Rans32x16O0 => 30,
            Self::Rans32x16O1 => 31,
            Self::Rans32x16O0R => 32,
            Self::Rans32x16O1R => 33,
            Self::Rans32x16O0P => 34,
            Self::Rans32x16O1P => 35,
            Self::Rans32x16O0PR => 36,
            Self::Rans32x16O1PR => 37,
            Self::RansNx16Stripe => 38,
            Self::RansNx16Cat => 39,
            Self::ArithmeticO0 => 40,
            Self::ArithmeticO1 => 41,
            Self::ArithmeticO0R => 42,
            Self::ArithmeticO1R => 43,
            Self::ArithmeticO0P => 44,
            Self::ArithmeticO1P => 45,
            Self::ArithmeticO0PR => 46,
            Self::ArithmeticO1PR => 47,
            Self::ArithmeticStripe => 48,
            Self::ArithmeticCat => 49,
            Self::ArithmeticExternal => 50,
            Self::TokenizerRans => 51,
            Self::TokenizerArithmetic => 52,
        }
    }
}

fn classify_rans_nx16(flags: u8) -> Method {
    if flags & 0x08 != 0 {
        return Method::RansNx16Stripe;
    }
    if flags & 0x20 != 0 {
        return Method::RansNx16Cat;
    }
    let order = flags & 0x01 != 0;
    let x32 = flags & 0x04 != 0;
    let rle = flags & 0x40 != 0;
    let pack = flags & 0x80 != 0;
    match (x32, order, rle, pack) {
        (false, false, false, false) => Method::Rans4x16O0,
        (false, true, false, false) => Method::Rans4x16O1,
        (false, false, true, false) => Method::Rans4x16O0R,
        (false, true, true, false) => Method::Rans4x16O1R,
        (false, false, false, true) => Method::Rans4x16O0P,
        (false, true, false, true) => Method::Rans4x16O1P,
        (false, false, true, true) => Method::Rans4x16O0PR,
        (false, true, true, true) => Method::Rans4x16O1PR,
        (true, false, false, false) => Method::Rans32x16O0,
        (true, true, false, false) => Method::Rans32x16O1,
        (true, false, true, false) => Method::Rans32x16O0R,
        (true, true, true, false) => Method::Rans32x16O1R,
        (true, false, false, true) => Method::Rans32x16O0P,
        (true, true, false, true) => Method::Rans32x16O1P,
        (true, false, true, true) => Method::Rans32x16O0PR,
        (true, true, true, true) => Method::Rans32x16O1PR,
    }
}

fn classify_arithmetic(flags: u8) -> Method {
    if flags & 0x08 != 0 {
        return Method::ArithmeticStripe;
    }
    if flags & 0x20 != 0 {
        return Method::ArithmeticCat;
    }
    if flags & 0x04 != 0 {
        return Method::ArithmeticExternal;
    }
    let variant = usize::from(flags & 0x03)
        + usize::from(flags & 0x40 != 0) * 2
        + usize::from(flags & 0x80 != 0) * 4;
    match variant {
        0 => Method::ArithmeticO0,
        1 => Method::ArithmeticO1,
        2 => Method::ArithmeticO0R,
        3 => Method::ArithmeticO1R,
        4 => Method::ArithmeticO0P,
        5 => Method::ArithmeticO1P,
        6 => Method::ArithmeticO0PR,
        7 => Method::ArithmeticO1PR,
        _ => unreachable!(),
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct EncodingSummary {
    pub data_series: String,
    pub encoding: String,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct ContainerEncodings {
    pub entries: Vec<EncodingSummary>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct MethodSummary {
    pub method: String,
    pub short: String,
    pub uncompressed_size: u64,
    pub compressed_size: u64,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct BlockSummary {
    pub content_id: Option<i32>,
    pub methods: Vec<MethodSummary>,
    pub data_series: Vec<String>,
    pub embedded_reference: bool,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Report {
    pub version: String,
    pub blocks: Vec<BlockSummary>,
    pub encodings: Vec<ContainerEncodings>,
    pub containers: u64,
    pub slices: u64,
    pub sequences: u64,
    pub bases: u64,
    pub file_size: u64,
    pub format_overhead_size: u64,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Options {
    pub verbose: bool,
    pub encodings: bool,
}

#[derive(Default)]
struct Accumulator {
    methods: BTreeMap<i32, Vec<(Method, u64, u64)>>,
    data_series: BTreeMap<i32, Vec<String>>,
    encodings: Vec<ContainerEncodings>,
    embedded_reference: Option<i32>,
    containers: u64,
    slices: u64,
    sequences: u64,
    bases: u64,
}

pub fn analyze(reader: impl Read) -> Result<Report> {
    parser::parse(BufReader::new(reader)).map_err(|error| match error.kind() {
        io::ErrorKind::InvalidData | io::ErrorKind::UnexpectedEof => {
            RsomicsError::InvalidInput(format!("invalid CRAM input: {error}"))
        }
        _ => RsomicsError::Io(error),
    })
}

impl Report {
    pub fn write(&self, options: Options, writer: impl Write) -> Result<()> {
        render::write(self, options, writer).map_err(RsomicsError::Io)
    }
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classifies_every_cram_compression_family() {
        for (code, data, expected) in [
            (0, &[][..], "raw"),
            (1, &[0; 9][..], "gzip"),
            (1, &[0, 0, 0, 0, 0, 0, 0, 0, 4][..], "gzip-min"),
            (1, &[0, 0, 0, 0, 0, 0, 0, 0, 2][..], "gzip-max"),
            (2, b"BZh6", "bzip2-6"),
            (3, &[][..], "lzma"),
            (4, &[0][..], "r4x8-o0"),
            (4, &[1][..], "r4x8-o1"),
            (5, &[0][..], "r4x16-o0"),
            (5, &[0x45][..], "r32x16-o1R"),
            (5, &[0x08][..], "rNx16-xo0"),
            (5, &[0x20][..], "rNx16-cat"),
            (6, &[0][..], "arith-o0"),
            (6, &[2][..], "arith-o0R"),
            (6, &[0x81][..], "arith-o1P"),
            (6, &[0x08][..], "arith-stripe"),
            (6, &[0x20][..], "arith-cat"),
            (6, &[0x04][..], "arith-ext"),
            (7, &[][..], "fqzcomp"),
            (8, &[0; 9][..], "tok3-rans"),
            (8, &[0, 0, 0, 0, 0, 0, 0, 0, 1][..], "tok3-arith"),
        ] {
            assert_eq!(Method::classify(code, data).unwrap().name(), expected);
        }
    }
}
