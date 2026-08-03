use std::collections::{HashMap, HashSet};

use noodles::sam;
use noodles::sam::header::record::value::map::{
    program::tag as program_tag, read_group::tag as read_group_tag,
};
use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Options {
    pub(crate) combine_read_groups: bool,
    pub(crate) combine_programs: bool,
}

#[derive(Clone, Debug)]
pub(crate) struct Translation {
    references: Vec<i32>,
    read_groups: HashMap<Vec<u8>, Vec<u8>>,
    programs: HashMap<Vec<u8>, Vec<u8>>,
}

impl Translation {
    pub(crate) fn preserves_reference_order(&self) -> bool {
        self.references.windows(2).all(|pair| pair[0] < pair[1])
    }
}

pub(crate) fn reconcile(
    headers: &[sam::Header],
    options: Options,
) -> Result<(sam::Header, Vec<Translation>)> {
    let mut output = sam::Header::default();
    if let Some(first) = headers.first() {
        *output.header_mut() = first.header().cloned();
    }

    let references = reconcile_references(headers, &mut output)?;
    let programs = reconcile_programs(headers, &mut output, options.combine_programs)?;
    let read_groups =
        reconcile_read_groups(headers, &mut output, &programs, options.combine_read_groups)?;

    for header in headers {
        output
            .comments_mut()
            .extend(header.comments().iter().cloned());
    }

    let translations = references
        .into_iter()
        .zip(read_groups)
        .zip(programs)
        .map(|((references, read_groups), programs)| Translation {
            references,
            read_groups,
            programs,
        })
        .collect();
    Ok((output, translations))
}

fn reconcile_references(
    headers: &[sam::Header],
    output: &mut sam::Header,
) -> Result<Vec<Vec<i32>>> {
    let mut translations = Vec::with_capacity(headers.len());
    for header in headers {
        let mut translation = Vec::with_capacity(header.reference_sequences().len());
        for (name, reference) in header.reference_sequences() {
            let id = match output.reference_sequences().get_index_of(name) {
                Some(id) => {
                    let existing = output
                        .reference_sequences_mut()
                        .get_mut(name)
                        .expect("resolved reference exists");
                    if existing.length() != reference.length() {
                        return Err(RsomicsError::InvalidInput(format!(
                            "reference {} has conflicting header records",
                            String::from_utf8_lossy(name)
                        )));
                    }
                    for (tag, value) in reference.other_fields() {
                        match existing.other_fields().get(tag) {
                            Some(current) if current != value => {
                                return Err(RsomicsError::InvalidInput(format!(
                                    "reference {} has conflicting header records",
                                    String::from_utf8_lossy(name)
                                )));
                            }
                            Some(_) => {}
                            None => {
                                existing.other_fields_mut().insert(*tag, value.clone());
                            }
                        }
                    }
                    id
                }
                None => {
                    let id = output.reference_sequences().len();
                    output
                        .reference_sequences_mut()
                        .insert(name.clone(), reference.clone());
                    id
                }
            };
            translation.push(i32::try_from(id).map_err(|_| {
                RsomicsError::InvalidInput("reference dictionary exceeds i32".to_owned())
            })?);
        }
        translations.push(translation);
    }
    Ok(translations)
}

fn reconcile_programs(
    headers: &[sam::Header],
    output: &mut sam::Header,
    combine: bool,
) -> Result<Vec<HashMap<Vec<u8>, Vec<u8>>>> {
    let mut used = HashSet::<Vec<u8>>::new();
    let mut translations = Vec::with_capacity(headers.len());

    for header in headers {
        let mut translation = HashMap::new();
        for (id, _) in header.programs().as_ref() {
            let original = id.as_slice().to_vec();
            let translated = if combine && used.contains(&original) {
                original.clone()
            } else {
                unique_id(&original, &used)?
            };
            used.insert(translated.clone());
            translation.insert(original, translated);
        }
        translations.push(translation);
    }

    for (header, translation) in headers.iter().zip(&translations) {
        for (id, program) in header.programs().as_ref() {
            let translated = &translation[id.as_slice()];
            if combine
                && output
                    .programs()
                    .as_ref()
                    .contains_key(translated.as_slice())
            {
                continue;
            }
            let mut program = program.clone();
            if let Some(previous) = program
                .other_fields_mut()
                .get_mut(&program_tag::PREVIOUS_PROGRAM_ID)
            {
                let mapped = translation.get(previous.as_slice()).ok_or_else(|| {
                    RsomicsError::InvalidInput(format!(
                        "program {} references unknown previous program {}",
                        String::from_utf8_lossy(id),
                        String::from_utf8_lossy(previous)
                    ))
                })?;
                *previous = mapped.clone().into();
            }
            output
                .programs_mut()
                .as_mut()
                .insert(translated.clone().into(), program);
        }
    }
    Ok(translations)
}

fn reconcile_read_groups(
    headers: &[sam::Header],
    output: &mut sam::Header,
    programs: &[HashMap<Vec<u8>, Vec<u8>>],
    combine: bool,
) -> Result<Vec<HashMap<Vec<u8>, Vec<u8>>>> {
    let mut used = HashSet::<Vec<u8>>::new();
    let mut translations = Vec::with_capacity(headers.len());

    for (header, program_translation) in headers.iter().zip(programs) {
        let mut translation = HashMap::new();
        for (id, read_group) in header.read_groups() {
            let original = id.as_slice().to_vec();
            let translated = if combine && used.contains(&original) {
                original.clone()
            } else {
                unique_id(&original, &used)?
            };
            used.insert(translated.clone());
            translation.insert(original, translated.clone());

            if combine && output.read_groups().contains_key(translated.as_slice()) {
                continue;
            }
            let mut read_group = read_group.clone();
            if let Some(program) = read_group
                .other_fields_mut()
                .get_mut(&read_group_tag::PROGRAM)
            {
                let mapped = program_translation.get(program.as_slice()).ok_or_else(|| {
                    RsomicsError::InvalidInput(format!(
                        "read group {} references unknown program {}",
                        String::from_utf8_lossy(id),
                        String::from_utf8_lossy(program)
                    ))
                })?;
                *program = mapped.clone().into();
            }
            output
                .read_groups_mut()
                .insert(translated.into(), read_group);
        }
        translations.push(translation);
    }
    Ok(translations)
}

fn unique_id(id: &[u8], used: &HashSet<Vec<u8>>) -> Result<Vec<u8>> {
    if !used.contains(id) {
        return Ok(id.to_vec());
    }
    for suffix in 1..=u64::MAX {
        let mut candidate = id.to_vec();
        candidate.push(b'.');
        candidate.extend_from_slice(suffix.to_string().as_bytes());
        if !used.contains(&candidate) {
            return Ok(candidate);
        }
    }
    Err(RsomicsError::InvalidInput(
        "header ID suffix exceeds u64".to_owned(),
    ))
}

pub(crate) fn translate(record: &mut RawRecord, translation: &Translation) -> Result<()> {
    record.set_reference_sequence_id(translate_reference(
        record.reference_sequence_id(),
        &translation.references,
        "reference",
    )?);
    record.set_mate_reference_sequence_id(translate_reference(
        record.mate_reference_sequence_id(),
        &translation.references,
        "mate reference",
    )?);
    translate_aux(record, *b"RG", &translation.read_groups, "read group")?;
    translate_aux(record, *b"PG", &translation.programs, "program")
}

fn translate_reference(value: i32, translation: &[i32], label: &str) -> Result<i32> {
    if value == -1 {
        return Ok(-1);
    }
    let id = usize::try_from(value)
        .ok()
        .and_then(|id| translation.get(id));
    id.copied().ok_or_else(|| {
        RsomicsError::InvalidInput(format!(
            "{label} sequence ID {value} is outside the input header dictionary"
        ))
    })
}

fn translate_aux(
    record: &mut RawRecord,
    tag: [u8; 2],
    translation: &HashMap<Vec<u8>, Vec<u8>>,
    label: &str,
) -> Result<()> {
    let Some(value) = record.aux_value(tag) else {
        return Ok(());
    };
    if record.aux_type(tag) != Some(b'Z') {
        return Err(RsomicsError::InvalidInput(format!(
            "read {} has a non-string {label} tag",
            String::from_utf8_lossy(record.name())
        )));
    }
    let value = value
        .strip_suffix(&[0])
        .ok_or_else(|| RsomicsError::InvalidInput(format!("read has unterminated {label} tag")))?;
    let mapped = translation.get(value).ok_or_else(|| {
        RsomicsError::InvalidInput(format!(
            "read {} references unknown {label} {}",
            String::from_utf8_lossy(record.name()),
            String::from_utf8_lossy(value)
        ))
    })?;
    if mapped.as_slice() != value {
        let mut encoded = mapped.clone();
        encoded.push(0);
        record.set_aux(tag, b'Z', &encoded)?;
    }
    Ok(())
}
