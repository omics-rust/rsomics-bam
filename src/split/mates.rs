use noodles::sam;
use rsomics_bamio::raw::RawRecord;

const UNMAPPED: u16 = 0x04;
const READ1: u16 = 0x40;
const RETAINED: u16 = 0x10 | 0x100 | 0x200 | 0x400;

pub(super) fn destination(flags: u16) -> usize {
    if flags & UNMAPPED != 0 {
        2
    } else if flags & READ1 != 0 {
        0
    } else {
        1
    }
}

pub(super) fn project_raw(record: &mut RawRecord) {
    record.clear_flag_bits(!RETAINED);
    record.set_mate_reference_sequence_id(-1);
    record.set_mate_alignment_start(-1);
    record.set_template_length(0);
}

pub(super) fn project(record: &mut sam::alignment::RecordBuf) {
    let flags = u16::from(record.flags()) & RETAINED;
    *record.flags_mut() = sam::alignment::record::Flags::from_bits_retain(flags);
    *record.mate_reference_sequence_id_mut() = None;
    *record.mate_alignment_start_mut() = None;
    *record.template_length_mut() = 0;
}
