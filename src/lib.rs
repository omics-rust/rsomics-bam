#![deny(unsafe_code)]

pub mod flags;
pub mod flagstat;
pub mod head;
pub mod quickcheck;
pub mod samples;
pub mod view;

mod cli;
mod commands;
mod filter;
mod hts_metadata;
mod hts_quickcheck;
mod input;
mod md;

#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
