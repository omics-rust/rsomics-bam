#![deny(unsafe_code)]

pub mod flags;
pub mod flagstat;
pub mod head;

mod cli;
mod input;
mod sam_format;

#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
