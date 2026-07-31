#![forbid(unsafe_code)]

pub mod flags;
pub mod flagstat;

mod cli;
mod input;

#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
