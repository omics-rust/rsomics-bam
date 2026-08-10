use std::path::{Path, PathBuf};
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};

use noodles::core::Region;
use noodles::sam;
use rsomics_bamio::raw::{RawRecord, RawRecordEncoder};
use rsomics_common::{Result, RsomicsError};

use crate::input;

enum Message {
    Header(Box<sam::Header>, input::Format),
    Record(RawRecord),
    Finished,
    Error(RsomicsError),
}

pub(crate) struct StreamHeader {
    pub(crate) header: sam::Header,
    pub(crate) format: input::Format,
}

struct Stream {
    input: PathBuf,
    receiver: Receiver<Message>,
    next: Option<RawRecord>,
    finished: bool,
    last_coordinate: Option<(i32, i32)>,
}

pub(crate) fn merge<T>(
    inputs: &[PathBuf],
    reference: Option<&Path>,
    additional_threads: usize,
    region: Option<&Region>,
    initialize: impl FnOnce(&[StreamHeader]) -> Result<T>,
    mut visit: impl FnMut(&mut T, usize, RawRecord) -> Result<()>,
) -> Result<T> {
    if inputs.is_empty() {
        return Err(RsomicsError::ConfigError(
            "at least one alignment input is required".to_owned(),
        ));
    }
    if inputs.len() > 1 && inputs.iter().any(|input| input == Path::new("-")) {
        return Err(RsomicsError::ConfigError(
            "standard input cannot be combined with other alignment inputs".to_owned(),
        ));
    }

    std::thread::scope(|scope| {
        let mut streams = Vec::with_capacity(inputs.len());
        for input in inputs {
            let (sender, receiver) = sync_channel(1);
            let input = input.clone();
            let worker_input = input.clone();
            let worker_region = region.cloned();
            scope.spawn(move || {
                if let Err(error) = read_input(
                    &worker_input,
                    reference,
                    additional_threads,
                    worker_region.as_ref(),
                    &sender,
                ) {
                    let _ = sender.send(Message::Error(error));
                }
            });
            streams.push(Stream {
                input,
                receiver,
                next: None,
                finished: false,
                last_coordinate: None,
            });
        }

        let mut headers = Vec::with_capacity(streams.len());
        for stream in &mut streams {
            headers.push(receive_header(stream)?);
        }
        let mut state = initialize(&headers)?;
        for stream in &mut streams {
            receive_next(stream)?;
        }

        while let Some(index) = streams
            .iter()
            .enumerate()
            .filter_map(|(index, stream)| {
                stream.next.as_ref().map(|record| {
                    (
                        index,
                        record.reference_sequence_id(),
                        record.alignment_start(),
                    )
                })
            })
            .min_by_key(|&(index, reference_id, position)| (reference_id, position, index))
            .map(|(index, _, _)| index)
        {
            let record = streams[index].next.take().unwrap();
            let coordinate = (record.reference_sequence_id(), record.alignment_start());
            if coordinate.0 >= 0 {
                if streams[index]
                    .last_coordinate
                    .is_some_and(|previous| coordinate < previous)
                {
                    return Err(RsomicsError::InvalidInput(format!(
                        "alignment input is not coordinate sorted: {}",
                        streams[index].input.display()
                    )));
                }
                streams[index].last_coordinate = Some(coordinate);
            }
            visit(&mut state, index, record)?;
            receive_next(&mut streams[index])?;
        }
        Ok(state)
    })
}

fn read_input(
    path: &Path,
    reference: Option<&Path>,
    additional_threads: usize,
    region: Option<&Region>,
    sender: &SyncSender<Message>,
) -> Result<()> {
    let mut reader = if region.is_some() {
        input::open_indexed(path, reference)?
    } else {
        input::open(path, reference, additional_threads)?
    };
    let header = reader.read_header(path)?;
    if sender
        .send(Message::Header(Box::new(header.clone()), reader.format()))
        .is_err()
    {
        return Ok(());
    }

    if let Some(region) = region {
        let mut encoder = RawRecordEncoder::new();
        reader.visit_region(&header, path, Some(region), |record| {
            let record = encoder.encode(&header, record)?;
            Ok(sender.send(Message::Record(record)).is_ok())
        })?;
    } else {
        reader.visit_owned_raw_records(&header, path, |record| {
            Ok(sender.send(Message::Record(record)).is_ok())
        })?;
    }
    let _ = sender.send(Message::Finished);
    Ok(())
}

fn receive_header(stream: &mut Stream) -> Result<StreamHeader> {
    match stream.receiver.recv() {
        Ok(Message::Header(header, format)) => Ok(StreamHeader {
            header: *header,
            format,
        }),
        Ok(Message::Error(error)) => Err(error),
        Ok(Message::Record(_) | Message::Finished) => Err(RsomicsError::InvalidInput(format!(
            "alignment stream ended before its header: {}",
            stream.input.display()
        ))),
        Err(_) => Err(RsomicsError::InvalidInput(format!(
            "alignment reader stopped unexpectedly: {}",
            stream.input.display()
        ))),
    }
}

fn receive_next(stream: &mut Stream) -> Result<()> {
    if stream.finished {
        stream.next = None;
        return Ok(());
    }
    match stream.receiver.recv() {
        Ok(Message::Record(record)) => stream.next = Some(record),
        Ok(Message::Finished) => {
            stream.next = None;
            stream.finished = true;
        }
        Ok(Message::Error(error)) => return Err(error),
        Ok(Message::Header(_, _)) => {
            return Err(RsomicsError::InvalidInput(format!(
                "alignment stream emitted more than one header: {}",
                stream.input.display()
            )));
        }
        Err(_) => {
            return Err(RsomicsError::InvalidInput(format!(
                "alignment reader stopped unexpectedly: {}",
                stream.input.display()
            )));
        }
    }
    Ok(())
}
