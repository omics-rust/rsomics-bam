use std::io::{self, Read};

pub(super) fn read_itf8(reader: &mut impl Read) -> io::Result<i32> {
    let b0 = u32::from(read_u8(reader)?);
    let value = if b0 & 0x80 == 0 {
        b0
    } else if b0 & 0x40 == 0 {
        (b0 & 0x7f) << 8 | u32::from(read_u8(reader)?)
    } else if b0 & 0x20 == 0 {
        (b0 & 0x3f) << 16 | u32::from(read_u16_be(reader)?)
    } else if b0 & 0x10 == 0 {
        (b0 & 0x1f) << 24 | read_u24_be(reader)?
    } else {
        let tail = read_u32_be(reader)?;
        (b0 & 0x0f) << 28 | (tail & 0xffff_fff0) >> 4 | tail & 0x0f
    };
    Ok(value as i32)
}

pub(super) fn read_nonnegative_itf8(reader: &mut impl Read, field: &str) -> io::Result<usize> {
    let value = read_itf8(reader)?;
    usize::try_from(value)
        .map_err(|_| super::invalid(format!("negative or oversized {field}: {value}")))
}

pub(super) fn read_ltf8(reader: &mut impl Read) -> io::Result<i64> {
    let b0 = u64::from(read_u8(reader)?);
    let value = if b0 & 0x80 == 0 {
        b0
    } else if b0 & 0x40 == 0 {
        (b0 & 0x7f) << 8 | u64::from(read_u8(reader)?)
    } else if b0 & 0x20 == 0 {
        (b0 & 0x3f) << 16 | u64::from(read_u16_be(reader)?)
    } else if b0 & 0x10 == 0 {
        (b0 & 0x1f) << 24 | u64::from(read_u24_be(reader)?)
    } else if b0 & 0x08 == 0 {
        (b0 & 0x0f) << 32 | u64::from(read_u32_be(reader)?)
    } else if b0 & 0x04 == 0 {
        (b0 & 0x07) << 40 | read_uint_be(reader, 5)?
    } else if b0 & 0x02 == 0 {
        (b0 & 0x03) << 48 | read_uint_be(reader, 6)?
    } else if b0 & 0x01 == 0 {
        read_uint_be(reader, 7)?
    } else {
        read_uint_be(reader, 8)?
    };
    Ok(value as i64)
}

pub(super) fn read_nonnegative_ltf8(reader: &mut impl Read, field: &str) -> io::Result<u64> {
    let value = read_ltf8(reader)?;
    u64::try_from(value).map_err(|_| super::invalid(format!("negative {field}: {value}")))
}

pub(super) fn read_u8(reader: &mut impl Read) -> io::Result<u8> {
    let mut byte = [0];
    reader.read_exact(&mut byte)?;
    Ok(byte[0])
}

pub(super) fn read_u32_le(reader: &mut impl Read) -> io::Result<u32> {
    let mut bytes = [0; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_u16_be(reader: &mut impl Read) -> io::Result<u16> {
    let mut bytes = [0; 2];
    reader.read_exact(&mut bytes)?;
    Ok(u16::from_be_bytes(bytes))
}

fn read_u24_be(reader: &mut impl Read) -> io::Result<u32> {
    let mut bytes = [0; 4];
    reader.read_exact(&mut bytes[1..])?;
    Ok(u32::from_be_bytes(bytes))
}

fn read_u32_be(reader: &mut impl Read) -> io::Result<u32> {
    let mut bytes = [0; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_be_bytes(bytes))
}

fn read_uint_be(reader: &mut impl Read, len: usize) -> io::Result<u64> {
    let mut bytes = [0; 8];
    reader.read_exact(&mut bytes[8 - len..])?;
    Ok(u64::from_be_bytes(bytes))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_itf8_boundaries_and_signed_values() {
        for (bytes, expected) in [
            (&[0x00][..], 0),
            (&[0x87, 0x55], 1877),
            (&[0xc7, 0x55, 0x99], 480_665),
            (&[0xe7, 0x55, 0x99, 0x66], 123_050_342),
            (&[0xf7, 0x55, 0x99, 0x66, 0x02], 1_968_805_474),
            (&[0xff, 0xff, 0xff, 0xff, 0x0f], -1),
        ] {
            assert_eq!(read_itf8(&mut &*bytes).unwrap(), expected);
        }
    }

    #[test]
    fn reads_ltf8_boundaries_and_signed_values() {
        for (bytes, expected) in [
            (&[0x00][..], 0),
            (&[0x81, 0x00], 256),
            (&[0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff], -1),
        ] {
            assert_eq!(read_ltf8(&mut &*bytes).unwrap(), expected);
        }
    }
}
