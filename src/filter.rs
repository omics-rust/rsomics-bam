use noodles::sam::alignment::Record;
use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct Filter {
    pub require_all: u16,
    pub exclude_any: u16,
    pub include_any: u16,
    pub exclude_all: u16,
    pub minimum_mapping_quality: u8,
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

        Ok(
            (self.require_all == 0 || flags & self.require_all == self.require_all)
                && flags & self.exclude_any == 0
                && (self.include_any == 0 || flags & self.include_any != 0)
                && (self.exclude_all == 0 || flags & self.exclude_all != self.exclude_all)
                && mapping_quality >= self.minimum_mapping_quality,
        )
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
        assert!(
            Filter {
                minimum_mapping_quality: u8::MAX,
                ..Filter::default()
            }
            .accepts(&record(0, None))
            .unwrap()
        );
        assert!(
            !Filter {
                minimum_mapping_quality: 1,
                ..Filter::default()
            }
            .accepts(&record(0x04, None))
            .unwrap()
        );
    }
}
