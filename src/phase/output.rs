use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::path::Path;

use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::{fasta, sam};
use noodles_util::alignment;
use rsomics_common::{Result, RsomicsError};

use super::model::Fragment;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Format {
    Sam,
    Bam,
    Cram,
}

impl Format {
    pub(crate) fn extension(self) -> &'static str {
        match self {
            Self::Sam => "sam",
            Self::Bam => "bam",
            Self::Cram => "cram",
        }
    }

    fn noodles(self) -> alignment::io::Format {
        match self {
            Self::Sam => alignment::io::Format::Sam,
            Self::Bam => alignment::io::Format::Bam,
            Self::Cram => alignment::io::Format::Cram,
        }
    }
}

pub(super) struct Router {
    header: sam::Header,
    writers: [alignment::io::Writer<File>; 3],
    pending: VecDeque<Pending>,
    random: Drand48,
    drop_ambiguous: bool,
}

struct Pending {
    reference_id: i32,
    end: i64,
    key: Box<[u8]>,
    record: RecordBuf,
    assignment: Option<Assignment>,
}

#[derive(Clone, Copy)]
struct Assignment {
    destination: Option<usize>,
    tagged: bool,
}

impl Router {
    pub(super) fn new(
        files: [File; 3],
        format: Format,
        repository: fasta::Repository,
        header: sam::Header,
        drop_ambiguous: bool,
    ) -> Result<Self> {
        let build = |file| {
            alignment::io::writer::Builder::default()
                .set_format(format.noodles())
                .set_reference_sequence_repository(repository.clone())
                .build_from_writer(file)
                .map_err(RsomicsError::Io)
        };
        let [first, second, chimera] = files;
        let mut writers = [build(first)?, build(second)?, build(chimera)?];
        for writer in &mut writers {
            writer.write_header(&header).map_err(RsomicsError::Io)?;
        }
        Ok(Self {
            header,
            writers,
            pending: VecDeque::new(),
            random: Drand48::default(),
            drop_ambiguous,
        })
    }

    pub(super) fn push(&mut self, reference_id: i32, end: i64, key: Box<[u8]>, record: RecordBuf) {
        self.pending.push_back(Pending {
            reference_id,
            end,
            key,
            record,
            assignment: None,
        });
    }

    pub(super) fn route(
        &mut self,
        reference_id: i32,
        cutoff: i64,
        fragments: &[Fragment],
    ) -> Result<()> {
        let flip = self.random.less_than_half();
        let assignments: HashMap<&[u8], _> = fragments
            .iter()
            .map(|fragment| {
                let assignment = if fragment.ambiguous {
                    Assignment {
                        destination: self.drop_ambiguous.then_some(2),
                        tagged: false,
                    }
                } else if fragment.phased && fragment.flipped {
                    Assignment {
                        destination: Some(2),
                        tagged: false,
                    }
                } else if fragment.phased {
                    Assignment {
                        destination: Some(usize::from(fragment.phase) ^ usize::from(flip)),
                        tagged: true,
                    }
                } else {
                    Assignment {
                        destination: None,
                        tagged: false,
                    }
                };
                (fragment.key.as_ref(), assignment)
            })
            .collect();

        for pending in &mut self.pending {
            if pending.assignment.is_some() {
                continue;
            }
            if let Some(&assignment) = assignments.get(pending.key.as_ref()) {
                pending.assignment = Some(assignment);
            } else if pending.reference_id < reference_id
                || (pending.reference_id == reference_id && pending.end <= cutoff)
            {
                pending.assignment = Some(Assignment {
                    destination: None,
                    tagged: false,
                });
            }
        }

        while self
            .pending
            .front()
            .is_some_and(|pending| pending.assignment.is_some())
        {
            let mut pending = self.pending.pop_front().unwrap();
            let assignment = pending.assignment.unwrap();
            let destination = assignment
                .destination
                .unwrap_or_else(|| usize::from(self.random.less_than_half()));
            if assignment.tagged {
                pending
                    .record
                    .data_mut()
                    .insert(Tag::from(*b"ZP"), Value::Character(b'Y'));
            }
            self.writers[destination]
                .write_record(&self.header, &pending.record)
                .map_err(RsomicsError::Io)?;
        }
        Ok(())
    }

    pub(super) fn finish(mut self) -> Result<()> {
        while let Some(mut pending) = self.pending.pop_front() {
            let assignment = pending.assignment.unwrap_or(Assignment {
                destination: None,
                tagged: false,
            });
            let destination = assignment
                .destination
                .unwrap_or_else(|| usize::from(self.random.less_than_half()));
            if assignment.tagged {
                pending
                    .record
                    .data_mut()
                    .insert(Tag::from(*b"ZP"), Value::Character(b'Y'));
            }
            self.writers[destination]
                .write_record(&self.header, &pending.record)
                .map_err(RsomicsError::Io)?;
        }
        for writer in &mut self.writers {
            writer.finish(&self.header).map_err(RsomicsError::Io)?;
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct Drand48 {
    state: u64,
}

impl Default for Drand48 {
    fn default() -> Self {
        Self {
            state: 0x1234_abcd_330e,
        }
    }
}

impl Drand48 {
    fn less_than_half(&mut self) -> bool {
        self.state = self.state.wrapping_mul(0x5deece66d).wrapping_add(0xb) & ((1_u64 << 48) - 1);
        self.state < 1_u64 << 47
    }
}

pub(crate) fn paths(prefix: &Path, format: Format) -> [std::path::PathBuf; 3] {
    let extension = format.extension();
    [
        format!("{}.0.{extension}", prefix.display()).into(),
        format!("{}.1.{extension}", prefix.display()).into(),
        format!("{}.chimera.{extension}", prefix.display()).into(),
    ]
}
