use std::collections::HashSet;
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::{Path, PathBuf};

use noodles::sam::alignment::Record;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::sam::header::record::value::map::read_group::tag as read_group_tag;
use rsomics_bamio::raw::RecordRef;
use rsomics_common::{Context, Result, RsomicsError};

#[derive(Debug, Eq, PartialEq)]
pub(crate) struct ReadGroupFilter {
    ids: HashSet<Vec<u8>>,
}

impl ReadGroupFilter {
    pub(crate) fn new(ids: &[String]) -> Option<Self> {
        if ids.is_empty() {
            return None;
        }

        Some(Self {
            ids: ids.iter().map(|id| id.as_bytes().to_vec()).collect(),
        })
    }

    pub(crate) fn retain_header(&self, header: &mut noodles::sam::Header) {
        header
            .read_groups_mut()
            .retain(|id, _| self.ids.contains::<[u8]>(id.as_ref()));
    }

    fn accepts(&self, record: &dyn Record) -> Result<bool> {
        with_read_group(record, |read_group| self.accepts_value(read_group))
    }

    fn accepts_raw(&self, record: &RecordRef<'_>) -> Result<bool> {
        with_raw_read_group(record, |read_group| self.accepts_value(read_group))
    }

    fn accepts_value(&self, value: Option<&[u8]>) -> bool {
        match value {
            Some(value) => self.ids.contains(value),
            None => true,
        }
    }
}

#[derive(Debug, Eq, PartialEq)]
pub(crate) struct LibraryFilter {
    read_groups: HashSet<Vec<u8>>,
}

impl LibraryFilter {
    pub(crate) fn new(header: &noodles::sam::Header, library: Option<&str>) -> Option<Self> {
        let library = library?;
        let read_groups = header
            .read_groups()
            .iter()
            .filter_map(|(id, read_group)| {
                read_group
                    .other_fields()
                    .get(&read_group_tag::LIBRARY)
                    .filter(|value| {
                        let value: &[u8] = value.as_ref();
                        value == library.as_bytes()
                    })
                    .map(|_| id.to_vec())
            })
            .collect();
        Some(Self { read_groups })
    }

    fn accepts(&self, record: &dyn Record) -> Result<bool> {
        with_read_group(record, |read_group| {
            read_group.is_some_and(|id| self.read_groups.contains(id))
        })
    }

    fn accepts_raw(&self, record: &RecordRef<'_>) -> Result<bool> {
        with_raw_read_group(record, |read_group| {
            read_group.is_some_and(|id| self.read_groups.contains(id))
        })
    }
}

fn with_read_group(
    record: &dyn Record,
    accept: impl FnOnce(Option<&[u8]>) -> bool,
) -> Result<bool> {
    let data = record.data();
    let read_group = match data
        .get(&Tag::READ_GROUP)
        .transpose()
        .map_err(RsomicsError::Io)?
    {
        Some(Value::String(value)) => Some(value.as_ref()),
        Some(_) => {
            return Err(RsomicsError::InvalidInput(
                "alignment RG tag must be a string".to_owned(),
            ));
        }
        None => None,
    };
    Ok(accept(read_group))
}

fn with_raw_read_group(
    record: &RecordRef<'_>,
    accept: impl FnOnce(Option<&[u8]>) -> bool,
) -> Result<bool> {
    let Some(value) = record.aux_value(*b"RG") else {
        return Ok(accept(None));
    };
    if record.aux_type(*b"RG") != Some(b'Z') {
        return Err(RsomicsError::InvalidInput(
            "alignment RG tag must be a string".to_owned(),
        ));
    }
    Ok(accept(Some(value.strip_suffix(&[0]).unwrap())))
}

#[derive(Debug, Eq, PartialEq)]
pub(crate) struct QnameFilter {
    names: HashSet<Vec<u8>>,
    exclude: bool,
}

impl QnameFilter {
    pub(crate) fn from_files(files: &[PathBuf]) -> Result<Option<Self>> {
        let Some(first) = files.first() else {
            return Ok(None);
        };
        let (exclude, _) = qname_file(first);
        let mut names = HashSet::new();

        for file in files {
            let (file_excludes, path) = qname_file(file);
            if file_excludes != exclude {
                return Err(RsomicsError::ConfigError(
                    "cannot mix include and exclude read-name files".to_owned(),
                ));
            }
            read_names(&path, &mut names)?;
        }

        Ok(Some(Self { names, exclude }))
    }

    fn accepts(&self, record: &dyn Record) -> bool {
        self.accepts_name(record.name().map_or(b"*".as_slice(), AsRef::as_ref))
    }

    fn accepts_raw(&self, record: &RecordRef<'_>) -> bool {
        self.accepts_name(record.name())
    }

    fn accepts_name(&self, name: &[u8]) -> bool {
        self.names.contains(name) != self.exclude
    }
}

fn qname_file(file: &Path) -> (bool, PathBuf) {
    let bytes = file.as_os_str().as_bytes();
    match bytes.strip_prefix(b"^") {
        Some(path) => (true, OsString::from_vec(path.to_vec()).into()),
        None => (false, file.to_owned()),
    }
}

fn read_names(path: &Path, names: &mut HashSet<Vec<u8>>) -> Result<()> {
    let file = File::open(path)
        .rs_with_context(|| format!("opening read-name file {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut line = Vec::new();

    loop {
        line.clear();
        if reader
            .read_until(b'\n', &mut line)
            .rs_with_context(|| format!("reading read-name file {}", path.display()))?
            == 0
        {
            return Ok(());
        }
        names.extend(
            line.split(|byte| byte.is_ascii_whitespace())
                .filter(|name| !name.is_empty())
                .map(<[u8]>::to_vec),
        );
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct Filter<'a> {
    pub require_all: u16,
    pub exclude_any: u16,
    pub include_any: u16,
    pub exclude_all: u16,
    pub read_groups: Option<&'a ReadGroupFilter>,
    pub qnames: Option<&'a QnameFilter>,
    pub library: Option<&'a LibraryFilter>,
    pub minimum_mapping_quality: u8,
    pub minimum_query_length: u64,
}

impl Filter<'_> {
    pub(crate) fn accepts(self, record: &dyn Record) -> Result<bool> {
        let record_flags = record.flags().map_err(RsomicsError::Io)?;
        let flags = u16::from(record_flags);
        let mapping_quality = record
            .mapping_quality()
            .transpose()
            .map_err(RsomicsError::Io)?
            .map_or_else(
                || {
                    if record_flags.is_unmapped() {
                        0
                    } else {
                        u8::MAX
                    }
                },
                |quality| quality.get(),
            );

        if !self.accepts_fields(flags, mapping_quality) {
            return Ok(false);
        }
        if let Some(read_groups) = self.read_groups
            && !read_groups.accepts(record)?
        {
            return Ok(false);
        }
        if let Some(qnames) = self.qnames
            && !qnames.accepts(record)
        {
            return Ok(false);
        }
        if let Some(library) = self.library
            && !library.accepts(record)?
        {
            return Ok(false);
        }
        if self.minimum_query_length == 0 {
            return Ok(true);
        }

        let query_length = record.cigar().read_length().map_err(RsomicsError::Io)? as u64;
        Ok(query_length >= self.minimum_query_length)
    }

    pub(crate) fn accepts_raw(self, record: &RecordRef<'_>) -> Result<bool> {
        let flags = record.flags();
        if !self.accepts_raw_fields(flags, record.mapping_quality()) {
            return Ok(false);
        }
        if let Some(read_groups) = self.read_groups
            && !read_groups.accepts_raw(record)?
        {
            return Ok(false);
        }
        if let Some(qnames) = self.qnames
            && !qnames.accepts_raw(record)
        {
            return Ok(false);
        }
        if let Some(library) = self.library
            && !library.accepts_raw(record)?
        {
            return Ok(false);
        }
        if self.minimum_query_length == 0 {
            return Ok(true);
        }

        let query_length = record
            .cigar_ops()
            .filter(|(kind, _)| matches!(kind, 0 | 1 | 4 | 7 | 8))
            .map(|(_, length)| u64::from(length))
            .sum::<u64>();
        Ok(query_length >= self.minimum_query_length)
    }

    fn accepts_raw_fields(self, flags: u16, mapping_quality: u8) -> bool {
        let mapping_quality = if mapping_quality == u8::MAX && flags & 0x04 != 0 {
            0
        } else {
            mapping_quality
        };
        self.accepts_fields(flags, mapping_quality)
    }

    fn accepts_fields(self, flags: u16, mapping_quality: u8) -> bool {
        (self.require_all == 0 || flags & self.require_all == self.require_all)
            && flags & self.exclude_any == 0
            && (self.include_any == 0 || flags & self.include_any != 0)
            && (self.exclude_all == 0 || flags & self.exclude_all != self.exclude_all)
            && mapping_quality >= self.minimum_mapping_quality
    }
}

#[cfg(test)]
mod tests {
    use noodles::sam::alignment::{
        RecordBuf,
        record::{Flags, MappingQuality},
    };

    use super::*;

    fn record(flags: u16, mapping_quality: Option<u8>) -> RecordBuf {
        let mut builder = RecordBuf::builder().set_flags(Flags::from(flags));
        if let Some(mapping_quality) = mapping_quality {
            builder = builder.set_mapping_quality(MappingQuality::new(mapping_quality).unwrap());
        }
        builder.build()
    }

    #[test]
    fn flag_predicates_follow_samtools_combinations() {
        let filter = Filter {
            require_all: 0x03,
            exclude_any: 0x100,
            include_any: 0xc0,
            exclude_all: 0x30,
            read_groups: None,
            qnames: None,
            library: None,
            minimum_mapping_quality: 20,
            minimum_query_length: 0,
        };

        assert!(filter.accepts(&record(0x43, Some(20))).unwrap());
        assert!(!filter.accepts(&record(0x41, Some(20))).unwrap());
        assert!(!filter.accepts(&record(0x143, Some(20))).unwrap());
        assert!(!filter.accepts(&record(0x03, Some(20))).unwrap());
        assert!(!filter.accepts(&record(0x73, Some(20))).unwrap());
        assert!(!filter.accepts(&record(0x43, Some(19))).unwrap());
    }

    #[test]
    fn missing_mapping_quality_distinguishes_mapped_and_unmapped_records() {
        let mapped = Filter {
            minimum_mapping_quality: u8::MAX,
            ..Filter::default()
        };
        let unmapped = Filter {
            minimum_mapping_quality: 1,
            ..Filter::default()
        };

        assert!(mapped.accepts(&record(0, None)).unwrap());
        assert!(!unmapped.accepts(&record(0x04, None)).unwrap());
    }

    #[test]
    fn minimum_query_length_uses_read_consuming_cigar_operations() {
        use noodles::sam::alignment::{
            record::cigar::{Op, op::Kind},
            record_buf::Cigar,
        };

        let cigar: Cigar = [
            Op::new(Kind::SoftClip, 2),
            Op::new(Kind::Match, 4),
            Op::new(Kind::Insertion, 1),
            Op::new(Kind::Deletion, 3),
            Op::new(Kind::SequenceMatch, 2),
            Op::new(Kind::SequenceMismatch, 1),
            Op::new(Kind::HardClip, 5),
        ]
        .into_iter()
        .collect();
        let record = RecordBuf::builder().set_cigar(cigar).build();

        assert!(
            Filter {
                minimum_query_length: 10,
                ..Filter::default()
            }
            .accepts(&record)
            .unwrap()
        );
        assert!(
            !Filter {
                minimum_query_length: 11,
                ..Filter::default()
            }
            .accepts(&record)
            .unwrap()
        );
    }

    #[test]
    fn raw_and_decoded_flag_predicates_match() {
        let filter = Filter {
            require_all: 0x03,
            exclude_any: 0x100,
            include_any: 0xc0,
            exclude_all: 0x30,
            read_groups: None,
            qnames: None,
            library: None,
            minimum_mapping_quality: 20,
            minimum_query_length: 0,
        };

        for (flags, mapping_quality) in [
            (0x43, 20),
            (0x41, 20),
            (0x143, 20),
            (0x03, 20),
            (0x73, 20),
            (0x43, 19),
        ] {
            assert_eq!(
                filter
                    .accepts(&record(flags, Some(mapping_quality)))
                    .unwrap(),
                filter.accepts_raw_fields(flags, mapping_quality)
            );
        }
    }
}
