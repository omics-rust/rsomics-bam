#![allow(unsafe_code)]

use std::io::{BufWriter, Write};
use std::path::Path;
use std::ptr::NonNull;
use std::str::FromStr;

use noodles::core::Region;
use rsomics_common::{Result, RsomicsError};
use rust_htslib::bam::{self, Read as _};
use rust_htslib::htslib;

use super::{Builder, Selection, Summary, hts_error, hts_references};

pub(super) fn write(
    input_path: &Path,
    region: Option<&str>,
    additional_threads: usize,
    output: impl Write,
) -> Result<Summary> {
    let region = region
        .map(Region::from_str)
        .transpose()
        .map_err(|error| RsomicsError::ConfigError(format!("invalid region: {error}")))?;
    if region.is_some() && input_path == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "--region requires a file-backed indexed input".to_owned(),
        ));
    }

    if let Some(region) = region.as_ref() {
        let mut reader = bam::IndexedReader::from_path(input_path)
            .map_err(|error| hts_error("opening indexed CRAM", input_path, error))?;
        configure(&mut reader, additional_threads, input_path)?;
        require_cram(&reader)?;
        let references = hts_references(reader.header())?;
        let selection = Selection::new(&references, region)?;
        reader
            .fetch((
                i32::try_from(selection.reference_id).unwrap(),
                i64::try_from(selection.start).unwrap(),
                i64::try_from(selection.end).unwrap(),
            ))
            .map_err(|error| hts_error("querying CRAM region", input_path, error))?;
        let builder = Builder::new(references, Some(selection), BufWriter::new(output));
        extract(&reader, builder)
    } else {
        let mut reader = if input_path == Path::new("-") {
            bam::Reader::from_stdin()
        } else {
            bam::Reader::from_path(input_path)
        }
        .map_err(|error| hts_error("opening CRAM", input_path, error))?;
        configure(&mut reader, additional_threads, input_path)?;
        require_cram(&reader)?;
        let references = hts_references(reader.header())?;
        let builder = Builder::new(references, None, BufWriter::new(output));
        extract(&reader, builder)
    }
}

fn configure(
    reader: &mut impl bam::Read,
    additional_threads: usize,
    input_path: &Path,
) -> Result<()> {
    if additional_threads > 0 {
        reader
            .set_threads(additional_threads)
            .map_err(|error| hts_error("configuring CRAM threads", input_path, error))?;
    }
    Ok(())
}

fn require_cram(reader: &impl bam::Read) -> Result<()> {
    let file = unsafe { reader.htsfile().as_ref() }.ok_or_else(|| {
        RsomicsError::InvalidInput("alignment reader has no underlying file".to_owned())
    })?;
    if file.format.format == htslib::htsExactFormat_cram {
        Ok(())
    } else {
        Err(RsomicsError::ConfigError(
            "--embedded requires CRAM input".to_owned(),
        ))
    }
}

fn extract<W: Write>(reader: &impl bam::Read, mut builder: Builder<W>) -> Result<Summary> {
    let file = unsafe { reader.htsfile().as_ref() }.unwrap();
    let cram = unsafe { file.fp.cram };
    if cram.is_null() {
        return Err(RsomicsError::InvalidInput(
            "CRAM reader has no decoder state".to_owned(),
        ));
    }

    loop {
        let Some(container) = Container::read(cram) else {
            break;
        };
        if unsafe { htslib::cram_container_is_empty(cram) } != 0 {
            Block::read(cram)?;
            continue;
        }
        Block::read(cram)?;
        let mut slice_count = 0;
        unsafe { htslib::cram_container_get_landmarks(container.0.as_ptr(), &mut slice_count) };
        let slice_count = usize::try_from(slice_count).map_err(|_| {
            RsomicsError::InvalidInput("CRAM container has a negative slice count".to_owned())
        })?;

        for _ in 0..slice_count {
            let block = Block::read(cram)?;
            let slice = SliceHeader::decode(cram, &block)?;
            let block_count =
                usize::try_from(unsafe { htslib::cram_slice_hdr_get_num_blocks(slice.0.as_ptr()) })
                    .map_err(|_| {
                        RsomicsError::InvalidInput(
                            "CRAM slice has a negative block count".to_owned(),
                        )
                    })?;
            let embedded_id = unsafe { htslib::cram_slice_hdr_get_embed_ref_id(slice.0.as_ptr()) };
            let mut reference_id = 0;
            let mut start = 0;
            let mut span = 0;
            unsafe {
                htslib::cram_slice_hdr_get_coords(
                    slice.0.as_ptr(),
                    &mut reference_id,
                    &mut start,
                    &mut span,
                )
            };

            if let Some(selection) = &builder.selection
                && (reference_id != i32::try_from(selection.reference_id).unwrap()
                    || start > i64::try_from(selection.end).unwrap())
            {
                return builder.finish();
            }
            if embedded_id < 0 && reference_id != -1 {
                return Err(RsomicsError::InvalidInput(
                    "CRAM slice has no embedded reference block".to_owned(),
                ));
            }

            let mut found = false;
            for _ in 0..block_count {
                let block = Block::read(cram)?;
                if block.content_id() != embedded_id {
                    continue;
                }
                block.uncompress()?;
                builder.add_embedded(reference_id, start, block.data()?)?;
                found = true;
            }
            if embedded_id >= 0 && !found {
                return Err(RsomicsError::InvalidInput(format!(
                    "CRAM slice is missing embedded reference block {embedded_id}"
                )));
            }
        }
    }
    if unsafe { htslib::cram_eof(cram) } == 0 {
        return Err(RsomicsError::InvalidInput(
            "reading CRAM container failed before EOF".to_owned(),
        ));
    }
    builder.finish()
}

struct Container(NonNull<htslib::cram_container>);

impl Container {
    fn read(cram: *mut htslib::cram_fd) -> Option<Self> {
        NonNull::new(unsafe { htslib::cram_read_container(cram) }).map(Self)
    }
}

impl Drop for Container {
    fn drop(&mut self) {
        unsafe { htslib::cram_free_container(self.0.as_ptr()) };
    }
}

struct Block(NonNull<htslib::cram_block>);

impl Block {
    fn read(cram: *mut htslib::cram_fd) -> Result<Self> {
        NonNull::new(unsafe { htslib::cram_read_block(cram) })
            .map(Self)
            .ok_or_else(|| RsomicsError::InvalidInput("reading CRAM block failed".to_owned()))
    }

    fn content_id(&self) -> i32 {
        unsafe { htslib::cram_block_get_content_id(self.0.as_ptr()) }
    }

    fn uncompress(&self) -> Result<()> {
        if unsafe { htslib::cram_uncompress_block(self.0.as_ptr()) } == 0 {
            Ok(())
        } else {
            Err(RsomicsError::InvalidInput(
                "decompressing embedded CRAM reference failed".to_owned(),
            ))
        }
    }

    fn data(&self) -> Result<&[u8]> {
        let length =
            usize::try_from(unsafe { htslib::cram_block_get_uncomp_size(self.0.as_ptr()) })
                .map_err(|_| {
                    RsomicsError::InvalidInput(
                        "embedded CRAM reference has a negative length".to_owned(),
                    )
                })?;
        let data = unsafe { htslib::cram_block_get_data(self.0.as_ptr()) }.cast::<u8>();
        if length == 0 {
            return Ok(&[]);
        }
        if data.is_null() {
            return Err(RsomicsError::InvalidInput(
                "embedded CRAM reference has no data".to_owned(),
            ));
        }
        Ok(unsafe { std::slice::from_raw_parts(data, length) })
    }
}

impl Drop for Block {
    fn drop(&mut self) {
        unsafe { htslib::cram_free_block(self.0.as_ptr()) };
    }
}

struct SliceHeader(NonNull<htslib::cram_block_slice_hdr>);

impl SliceHeader {
    fn decode(cram: *mut htslib::cram_fd, block: &Block) -> Result<Self> {
        NonNull::new(unsafe { htslib::cram_decode_slice_header(cram, block.0.as_ptr()) })
            .map(Self)
            .ok_or_else(|| {
                RsomicsError::InvalidInput("decoding CRAM slice header failed".to_owned())
            })
    }
}

impl Drop for SliceHeader {
    fn drop(&mut self) {
        unsafe { htslib::cram_free_slice_header(self.0.as_ptr()) };
    }
}
