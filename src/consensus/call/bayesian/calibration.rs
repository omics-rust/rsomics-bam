use super::QualityCalibration;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum CalibrationPreset {
    Flat,
    Hifi,
    Hiseq,
    R10_4Sup,
    R10_4Dup,
    Ultima,
}

pub(super) fn preset(preset: CalibrationPreset) -> QualityCalibration {
    let maps = match preset {
        CalibrationPreset::Flat => {
            let identity = std::array::from_fn(|quality| quality as u8);
            return QualityCalibration {
                substitution: identity,
                undercall: identity,
                overcall: identity,
            };
        }
        CalibrationPreset::Hifi | CalibrationPreset::R10_4Dup => HIFI,
        CalibrationPreset::Hiseq => HISEQ,
        CalibrationPreset::R10_4Sup => R10_4_SUP,
        CalibrationPreset::Ultima => ULTIMA,
    };
    QualityCalibration {
        substitution: maps[0],
        undercall: maps[1],
        overcall: maps[2],
    }
}

pub(super) fn read(path: &Path) -> Result<QualityCalibration> {
    let input = File::open(path).map_err(|error| io_error(path, "opening", error))?;
    let mut calibration = QualityCalibration::default();
    let mut maximum = 0usize;
    for (index, result) in BufReader::new(input).lines().enumerate() {
        let line = result.map_err(|error| io_error(path, "reading", error))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split_ascii_whitespace().collect::<Vec<_>>();
        if fields.len() != 5 || fields[0] != "QUAL" {
            return Err(invalid(
                path,
                index + 1,
                "expected QUAL q substitution undercall overcall",
            ));
        }
        let quality = parse::<usize>(path, index + 1, fields[1], "quality")?;
        if quality > 100 {
            return Err(invalid(path, index + 1, "quality exceeds 100"));
        }
        if quality < maximum {
            return Err(invalid(path, index + 1, "qualities are not ascending"));
        }
        let substitution = parse::<u8>(path, index + 1, fields[2], "substitution")?;
        let undercall = parse::<u8>(path, index + 1, fields[3], "undercall")?;
        let overcall = parse::<u8>(path, index + 1, fields[4], "overcall")?;
        for current in maximum + 1..=quality {
            calibration.substitution[current] = calibration.substitution[current - 1];
            calibration.undercall[current] = calibration.undercall[current - 1];
            calibration.overcall[current] = calibration.overcall[current - 1];
        }
        if quality < 100 {
            calibration.substitution[quality] = substitution;
            calibration.undercall[quality] = undercall;
            calibration.overcall[quality] = overcall;
        }
        maximum = quality;
    }
    for quality in maximum + 1..=100 {
        calibration.substitution[quality] = calibration.substitution[maximum];
        calibration.undercall[quality] = calibration.undercall[maximum];
        calibration.overcall[quality] = calibration.overcall[maximum];
    }
    Ok(calibration)
}

fn parse<T: std::str::FromStr>(path: &Path, line: usize, value: &str, field: &str) -> Result<T> {
    value
        .parse()
        .map_err(|_| invalid(path, line, &format!("invalid {field}")))
}

fn invalid(path: &Path, line: usize, reason: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{}:{line}: {reason}", path.display()))
}

fn io_error(path: &Path, action: &str, error: std::io::Error) -> RsomicsError {
    RsomicsError::Io(std::io::Error::new(
        error.kind(),
        format!("{action} quality calibration {}: {error}", path.display()),
    ))
}

const HIFI: [[u8; 101]; 3] = [
    [
        10, 11, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 27, 28, 29, 30, 31, 32, 33,
        33, 34, 35, 36, 36, 37, 38, 38, 39, 39, 40, 40, 41, 41, 41, 41, 42, 42, 42, 42, 43, 43, 43,
        43, 43, 43, 43, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
        44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
        44, 44, 44, 44, 44, 44, 44, 44, 0,
    ],
    [
        4, 4, 4, 4, 5, 6, 6, 7, 8, 9, 10, 11, 11, 12, 13, 14, 15, 15, 16, 17, 18, 19, 19, 20, 20,
        21, 22, 23, 23, 24, 25, 25, 25, 26, 26, 26, 27, 27, 28, 28, 28, 28, 27, 27, 27, 28, 28, 28,
        28, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 26, 26, 25, 26, 26, 27, 27, 27, 26,
        26, 26, 26, 26, 26, 26, 26, 27, 27, 28, 29, 28, 28, 28, 27, 27, 27, 27, 27, 27, 28, 28, 30,
        30, 30, 30, 30, 30, 30, 0,
    ],
    [
        8, 8, 8, 8, 9, 10, 11, 12, 13, 14, 15, 15, 16, 17, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23,
        23, 23, 24, 24, 24, 25, 25, 25, 25, 25, 25, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 28, 28,
        28, 28, 28, 29, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
        30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
        30, 30, 30, 30, 30, 30, 30, 0,
    ],
];

const HISEQ: [[u8; 101]; 3] = [
    [
        2, 2, 2, 3, 3, 4, 5, 5, 6, 7, 8, 9, 10, 11, 11, 12, 13, 14, 15, 16, 17, 17, 18, 19, 20, 21,
        22, 22, 23, 24, 25, 26, 27, 28, 28, 29, 30, 31, 32, 33, 34, 34, 35, 36, 37, 38, 39, 39, 40,
        41, 42, 43, 44, 45, 45, 46, 47, 48, 49, 50, 51, 51, 52, 53, 54, 55, 56, 56, 57, 58, 59, 60,
        61, 62, 62, 63, 64, 65, 66, 67, 68, 68, 69, 70, 71, 72, 73, 73, 74, 75, 76, 77, 78, 79, 79,
        80, 81, 82, 83, 84, 0,
    ],
    [
        1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 25, 26, 27, 28, 29,
        31, 32, 33, 34, 35, 37, 38, 39, 40, 41, 43, 44, 45, 46, 47, 49, 50, 51, 52, 53, 55, 56, 57,
        58, 59, 61, 62, 63, 64, 65, 67, 68, 69, 70, 71, 73, 74, 75, 76, 77, 79, 80, 81, 82, 83, 85,
        86, 87, 88, 89, 91, 92, 93, 94, 95, 97, 98, 99, 100, 101, 103, 104, 105, 106, 107, 109,
        110, 111, 112, 113, 115, 116, 117, 118, 119, 0,
    ],
    [
        1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 25, 26, 27, 28, 29,
        31, 32, 33, 34, 35, 37, 38, 39, 40, 41, 43, 44, 45, 46, 47, 49, 50, 51, 52, 53, 55, 56, 57,
        58, 59, 61, 62, 63, 64, 65, 67, 68, 69, 70, 71, 73, 74, 75, 76, 77, 79, 80, 81, 82, 83, 85,
        86, 87, 88, 89, 91, 92, 93, 94, 95, 97, 98, 99, 100, 101, 103, 104, 105, 106, 107, 109,
        110, 111, 112, 113, 115, 116, 117, 118, 119, 0,
    ],
];

const R10_4_SUP: [[u8; 101]; 3] = [
    [
        0, 2, 2, 2, 3, 4, 4, 5, 6, 7, 7, 8, 9, 12, 13, 14, 15, 15, 16, 17, 18, 19, 20, 22, 24, 25,
        26, 27, 28, 29, 30, 31, 33, 34, 36, 37, 38, 38, 39, 39, 40, 40, 40, 40, 40, 40, 40, 41, 40,
        40, 41, 41, 40, 40, 40, 40, 41, 40, 40, 40, 40, 41, 41, 40, 40, 41, 40, 40, 39, 41, 40, 41,
        40, 40, 41, 41, 41, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
        40, 40, 40, 40, 40, 0,
    ],
    [
        0, 2, 2, 2, 3, 4, 5, 6, 7, 8, 8, 9, 9, 10, 10, 10, 11, 12, 12, 13, 13, 13, 14, 14, 15, 16,
        16, 17, 18, 18, 19, 19, 20, 21, 22, 23, 24, 25, 25, 25, 25, 25, 25, 25, 25, 25, 26, 26, 26,
        26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28,
        28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
        28, 28, 28, 28, 28, 0,
    ],
    [
        0, 4, 6, 6, 6, 7, 7, 8, 9, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 15, 15, 15, 16, 16,
        17, 17, 18, 18, 19, 19, 20, 20, 21, 22, 22, 23, 23, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
        24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
        24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
        24, 24, 24, 24, 24, 24, 0,
    ],
];

const ULTIMA: [[u8; 101]; 3] = [
    [
        2, 2, 3, 4, 5, 6, 6, 7, 8, 9, 10, 10, 11, 12, 13, 14, 14, 15, 16, 17, 18, 18, 19, 21, 22,
        23, 23, 24, 25, 26, 27, 27, 28, 29, 30, 31, 31, 32, 33, 34, 35, 35, 36, 37, 38, 39, 39, 40,
        42, 43, 44, 44, 45, 46, 47, 48, 48, 49, 50, 51, 52, 52, 53, 54, 55, 56, 56, 57, 58, 59, 60,
        60, 61, 63, 64, 65, 65, 66, 67, 68, 69, 69, 70, 71, 72, 73, 73, 74, 75, 76, 77, 77, 78, 79,
        80, 81, 81, 82, 84, 85, 0,
    ],
    [
        1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 10, 10, 10, 11, 12, 13, 13, 13,
        14, 15, 16, 16, 16, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 22, 22, 23, 23, 24, 24, 25,
        25, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 28, 28,
        28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
        28, 28, 28, 28, 0,
    ],
    [
        1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 10, 10, 10, 11, 12, 13, 13, 13,
        14, 15, 16, 16, 16, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 22, 22, 23, 23, 24, 24, 25,
        25, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 28, 28,
        28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
        28, 28, 28, 28, 0,
    ],
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preset_sentinels_match_samtools_1_24() {
        let flat = preset(CalibrationPreset::Flat);
        assert_eq!(flat.substitution[100], 100);
        assert_eq!(QualityCalibration::default().substitution[100], 0);

        let hifi = preset(CalibrationPreset::Hifi);
        assert_eq!(hifi.substitution[0], 10);
        assert_eq!(hifi.substitution[50], 44);
        assert_eq!(hifi.undercall[42], 27);
        assert_eq!(hifi.overcall[56], 30);

        let hiseq = preset(CalibrationPreset::Hiseq);
        assert_eq!(hiseq.substitution[99], 84);
        assert_eq!(hiseq.undercall[99], 119);

        let sup = preset(CalibrationPreset::R10_4Sup);
        assert_eq!(sup.substitution[32], 33);
        assert_eq!(sup.undercall[67], 28);
        assert_eq!(sup.overcall[40], 24);

        let dup = preset(CalibrationPreset::R10_4Dup);
        assert_eq!(dup.substitution, hifi.substitution);
        assert_eq!(dup.undercall, hifi.undercall);
        assert_eq!(dup.overcall, hifi.overcall);

        let ultima = preset(CalibrationPreset::Ultima);
        assert_eq!(ultima.substitution[48], 42);
        assert_eq!(ultima.undercall[71], 28);
        assert_eq!(ultima.overcall[99], 28);

        for preset in [
            CalibrationPreset::Hifi,
            CalibrationPreset::Hiseq,
            CalibrationPreset::R10_4Sup,
            CalibrationPreset::R10_4Dup,
            CalibrationPreset::Ultima,
        ] {
            let calibration = super::preset(preset);
            assert_eq!(calibration.substitution[100], 0);
            assert_eq!(calibration.undercall[100], 0);
            assert_eq!(calibration.overcall[100], 0);
        }
    }

    #[test]
    fn custom_file_interpolates_and_extends_like_samtools_1_24() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("calibration.txt");
        std::fs::write(
            &path,
            "# calibration\nQUAL 0 3 4 5\nQUAL 3 13 14 15\nQUAL 5 23 24 25\n",
        )
        .unwrap();

        let calibration = read(&path).unwrap();

        assert_eq!(&calibration.substitution[..7], &[3, 3, 3, 13, 13, 23, 23]);
        assert_eq!(&calibration.undercall[..7], &[4, 4, 4, 14, 14, 24, 24]);
        assert_eq!(&calibration.overcall[..7], &[5, 5, 5, 15, 15, 25, 25]);
        assert_eq!(calibration.substitution[100], 23);
    }

    #[test]
    fn custom_file_rejects_ambiguous_or_unsafe_rows() {
        let directory = tempfile::tempdir().unwrap();
        for (name, contents) in [
            ("descending", "QUAL 5 1 1 1\nQUAL 4 2 2 2\n"),
            ("quality", "QUAL 101 1 1 1\n"),
            ("mapping", "QUAL 1 256 1 1\n"),
            ("shape", "QUAL 1 2 3\n"),
            ("tag", "OTHER 1 2 3 4\n"),
        ] {
            let path = directory.path().join(name);
            std::fs::write(&path, contents).unwrap();
            assert!(read(&path).is_err(), "{name}");
        }
    }
}
