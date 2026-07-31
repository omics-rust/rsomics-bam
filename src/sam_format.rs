#![allow(unsafe_code)]

use std::ptr;
use std::slice;

use rsomics_common::{Result, RsomicsError};
use rust_htslib::bam::{HeaderView, Record};
use rust_htslib::htslib;

struct Buffer(htslib::kstring_t);

impl Buffer {
    fn new() -> Self {
        Self(htslib::kstring_t {
            l: 0,
            m: 0,
            s: ptr::null_mut(),
        })
    }

    fn as_bytes(&self) -> &[u8] {
        if self.0.l == 0 {
            return &[];
        }

        // HTSlib keeps at least `l` initialized bytes in `s` until the buffer is freed.
        unsafe { slice::from_raw_parts(self.0.s.cast(), self.0.l) }
    }
}

impl Drop for Buffer {
    fn drop(&mut self) {
        // `sam_format1` allocates `s` with the C allocator.
        unsafe { libc::free(self.0.s.cast()) }
    }
}

pub(crate) fn record(header: &HeaderView, record: &Record) -> Result<Vec<u8>> {
    let mut output = Buffer::new();

    // Both pointers remain live for the call, and `output` owns the returned allocation.
    let status = unsafe {
        htslib::sam_format1(
            header.inner_ptr(),
            record.inner() as *const htslib::bam1_t,
            &mut output.0,
        )
    };

    if status < 0 {
        Err(RsomicsError::InvalidInput(
            "formatting alignment record as SAM".to_owned(),
        ))
    } else {
        Ok(output.as_bytes().to_vec())
    }
}
