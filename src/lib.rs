#![deny(unsafe_code)]

//! Programmatic interfaces for the rsomics SAM, BAM, and CRAM product.

/// Per-position alignment depth.
pub mod depth;
/// SAM flag parsing and rendering.
pub mod flags;
/// Alignment flag statistics.
pub mod flagstat;
/// Header and leading-record inspection.
pub mod head;
/// Alignment random-access index construction.
pub mod index;
/// Ordered alignment merging.
pub mod merge;
/// Per-position pileup generation.
pub mod mpileup;
/// Alignment-header program provenance.
mod program;
/// Alignment integrity checks.
pub mod quickcheck;
/// Read-group sample metadata.
pub mod samples;
/// Bounded-memory alignment sorting.
pub mod sort;
/// Alignment filtering and format conversion.
pub mod view;

pub use program::Program;

mod alignment_order;
mod cli;
mod commands;
mod filter;
mod header_merge;
mod hts_metadata;
mod hts_quickcheck;
mod input;
mod md;
mod output;

#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
