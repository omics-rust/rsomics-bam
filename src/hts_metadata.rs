#![allow(unsafe_code)]

use std::ffi::{CStr, CString};
use std::path::{Path, PathBuf};
use std::ptr;

use rsomics_common::{Result, RsomicsError};
use rust_htslib::bam::{Read, Reader};
use rust_htslib::htslib;

pub(crate) struct ReferenceDictionary {
    pub path: PathBuf,
    pub targets: Vec<(Vec<u8>, u64)>,
}

pub(crate) fn load_reference(path: &Path) -> Result<ReferenceDictionary> {
    let path_string = c_path(path)?;
    let mut targets = Vec::new();

    // The faidx handle owns every returned name and remains live until all names are copied.
    unsafe {
        let index = htslib::fai_load(path_string.as_ptr());
        if index.is_null() {
            return Err(RsomicsError::InvalidInput(format!(
                "loading reference index for {}",
                path.display()
            )));
        }

        let count = htslib::faidx_nseq(index);
        if count < 0 {
            htslib::fai_destroy(index);
            return Err(RsomicsError::InvalidInput(format!(
                "reading reference index for {}",
                path.display()
            )));
        }

        for position in 0..count {
            let name = htslib::faidx_iseq(index, position);
            if name.is_null() {
                htslib::fai_destroy(index);
                return Err(RsomicsError::InvalidInput(format!(
                    "reading reference name {position} from {}",
                    path.display()
                )));
            }
            let length = htslib::faidx_seq_len64(index, name);
            if length < 0 {
                htslib::fai_destroy(index);
                return Err(RsomicsError::InvalidInput(format!(
                    "reading reference length {position} from {}",
                    path.display()
                )));
            }
            targets.push((CStr::from_ptr(name).to_bytes().to_vec(), length as u64));
        }

        htslib::fai_destroy(index);
    }

    Ok(ReferenceDictionary {
        path: path.to_path_buf(),
        targets,
    })
}

pub(crate) fn has_index(reader: &Reader, alignment: &Path, index: Option<&Path>) -> Result<bool> {
    let alignment = c_path(alignment)?;
    let index = index.map(c_path).transpose()?;

    // The loaded index is independent of the reader and is destroyed before returning.
    let loaded = unsafe {
        htslib::sam_index_load3(
            reader.htsfile(),
            alignment.as_ptr(),
            index.as_ref().map_or(ptr::null(), |path| path.as_ptr()),
            htslib::HTS_IDX_SILENT_FAIL as libc::c_int,
        )
    };
    if loaded.is_null() {
        Ok(false)
    } else {
        unsafe { htslib::hts_idx_destroy(loaded) };
        Ok(true)
    }
}

fn c_path(path: &Path) -> Result<CString> {
    CString::new(path.as_os_str().as_encoded_bytes()).map_err(|_| {
        RsomicsError::InvalidInput(format!("path contains a null byte: {}", path.display()))
    })
}
