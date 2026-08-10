use noodles::sam;
use noodles::sam::header::record::value::{
    Map,
    map::{Program as SamProgram, program::tag},
};
use rsomics_common::{Result, RsomicsError};

pub(crate) fn command_line() -> String {
    std::env::args_os()
        .map(|argument| argument.to_string_lossy().replace(['\t', '\r', '\n'], " "))
        .collect::<Vec<_>>()
        .join(" ")
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Program<'a> {
    name: &'a str,
    version: &'a str,
    command_line: &'a str,
}

impl<'a> Program<'a> {
    pub fn new(name: &'a str, version: &'a str, command_line: &'a str) -> Result<Self> {
        if name.is_empty() || version.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "program name and version cannot be empty".to_owned(),
            ));
        }
        if [name, version, command_line]
            .iter()
            .any(|value| value.contains(['\t', '\r', '\n']))
        {
            return Err(RsomicsError::InvalidInput(
                "program header fields cannot contain tabs or line breaks".to_owned(),
            ));
        }
        Ok(Self {
            name,
            version,
            command_line,
        })
    }

    pub(crate) fn add_to(self, header: &mut sam::Header) -> Result<Vec<u8>> {
        let programs = header.programs_mut().as_mut();
        let previous = programs.last().map(|(id, _)| id.clone());
        let mut id = self.name.to_owned();
        let mut suffix = 0u64;

        while programs.contains_key(id.as_bytes()) {
            suffix = suffix.checked_add(1).ok_or_else(|| {
                RsomicsError::InvalidInput("program ID suffix exceeds u64".to_owned())
            })?;
            id = format!("{}.{}", self.name, suffix);
        }

        let mut builder = Map::<SamProgram>::builder()
            .insert(tag::NAME, self.name)
            .insert(tag::VERSION, self.version)
            .insert(tag::COMMAND_LINE, self.command_line);
        if let Some(previous) = previous {
            builder = builder.insert(tag::PREVIOUS_PROGRAM_ID, previous);
        }
        let map = builder.build().map_err(|error| {
            RsomicsError::InvalidInput(format!("building program record: {error}"))
        })?;
        let id = id.into_bytes();
        programs.insert(id.clone().into(), map);
        Ok(id)
    }
}

#[cfg(test)]
mod tests {
    use noodles::sam::header::record::value::map::program::tag;

    use super::*;

    #[test]
    fn adds_after_the_last_existing_record_with_a_unique_id() {
        let mut header = sam::Header::builder()
            .add_program("rsomics-bam", Map::default())
            .add_program("aligner", Map::default())
            .build();
        Program::new("rsomics-bam", "1.2.3", "rsomics-bam view input.bam")
            .unwrap()
            .add_to(&mut header)
            .unwrap();

        let program = &header.programs().as_ref()[b"rsomics-bam.1".as_slice()];
        assert_eq!(
            program
                .other_fields()
                .get(&tag::PREVIOUS_PROGRAM_ID)
                .map(|value| value.as_ref()),
            Some(b"aligner".as_slice())
        );
        assert_eq!(
            program
                .other_fields()
                .get(&tag::NAME)
                .map(|value| value.as_ref()),
            Some(b"rsomics-bam".as_slice())
        );
        assert_eq!(
            program
                .other_fields()
                .get(&tag::VERSION)
                .map(|value| value.as_ref()),
            Some(b"1.2.3".as_slice())
        );
    }

    #[test]
    fn fields_cannot_create_header_fields_or_lines() {
        assert!(Program::new("rsomics-bam", "1.2.3", "view\tinput.bam").is_err());
        assert!(Program::new("", "1.2.3", "view input.bam").is_err());
    }
}
