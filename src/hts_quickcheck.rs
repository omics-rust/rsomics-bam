#![allow(unsafe_code)]

use std::ffi::CString;
use std::path::Path;
use std::sync::Mutex;

use rust_htslib::htslib;

use crate::quickcheck::{FileReport, Problem};

static HTSLIB_LOG_LOCK: Mutex<()> = Mutex::new(());

pub(crate) fn check(path: &Path, allow_no_targets: bool) -> FileReport {
    let Ok(path_string) = CString::new(path.as_os_str().as_encoded_bytes()) else {
        return FileReport {
            path: path.to_path_buf(),
            problems: vec![Problem::Open],
        };
    };

    let log_guard = match HTSLIB_LOG_LOCK.lock() {
        Ok(guard) => guard,
        Err(poisoned) => poisoned.into_inner(),
    };

    // HTSlib logging is process-global, so all temporary changes are serialized and restored.
    let previous_log_level = unsafe { htslib::hts_get_log_level() };
    unsafe { htslib::hts_set_log_level(htslib::htsLogLevel_HTS_LOG_OFF) };
    let problems = unsafe { inspect(path_string.as_ptr(), allow_no_targets) };
    unsafe { htslib::hts_set_log_level(previous_log_level) };
    drop(log_guard);

    FileReport {
        path: path.to_path_buf(),
        problems,
    }
}

unsafe fn inspect(path: *const libc::c_char, allow_no_targets: bool) -> Vec<Problem> {
    let mut problems = Vec::new();

    // All pointers originate from HTSlib and are destroyed before this function returns.
    unsafe {
        let file = htslib::hts_open(path, c"r".as_ptr());
        if file.is_null() {
            problems.push(Problem::Open);
            return problems;
        }

        let format = htslib::hts_get_format(file);
        if format.is_null() || (*format).category != htslib::htsFormatCategory_sequence_data {
            problems.push(Problem::NotSequence);
        } else {
            let header = htslib::sam_hdr_read(file);
            if header.is_null() {
                problems.push(Problem::Header);
            } else {
                if !allow_no_targets && htslib::sam_hdr_nref(header) <= 0 {
                    problems.push(Problem::NoTargets);
                }
                htslib::sam_hdr_destroy(header);
            }
        }

        match htslib::hts_check_EOF(file) {
            i32::MIN..=-1 => problems.push(Problem::EofCheck),
            0 => problems.push(Problem::MissingEof),
            1..=3 => {}
            _ => problems.push(Problem::EofCheck),
        }

        if htslib::hts_close(file) < 0 {
            problems.push(Problem::Close);
        }
    }

    problems
}
