use noodles::sam::alignment::Record;
use rsomics_bamio::raw::RecordRef;
use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct Filter {
    pub require_all: u16,
    pub exclude_any: u16,
    pub include_any: u16,
    pub exclude_all: u16,
    pub minimum_mapping_quality: u8,
    pub minimum_query_length: u64,
}

impl Filter {
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
        if self.minimum_query_length == 0 {
            return Ok(true);
        }

        let query_length = record.cigar().read_length().map_err(RsomicsError::Io)? as u64;
        Ok(query_length >= self.minimum_query_length)
    }

    pub(crate) fn accepts_raw(self, record: &RecordRef<'_>) -> bool {
        let flags = record.flags();
        if !self.accepts_raw_fields(flags, record.mapping_quality()) {
            return false;
        }
        if self.minimum_query_length == 0 {
            return true;
        }

        let query_length = record
            .cigar_ops()
            .filter(|(kind, _)| matches!(kind, 0 | 1 | 4 | 7 | 8))
            .map(|(_, length)| u64::from(length))
            .sum::<u64>();
        query_length >= self.minimum_query_length
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
