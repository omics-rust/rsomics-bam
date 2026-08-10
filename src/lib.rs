#![deny(unsafe_code)]

//! Programmatic interfaces for the rsomics SAM, BAM, and CRAM product.

/// Alignment read-group editing.
pub mod addreplacerg;
/// Amplicon primer clipping.
pub mod ampliconclip;
/// Amplicon sequencing statistics.
pub mod ampliconstats;
/// Coverage totals over BED regions.
pub mod bedcov;
/// Alignment MD and NM tag recalculation.
pub mod calmd;
/// Compressed-block BAM concatenation.
pub mod cat;
/// Bounded-memory read-name collation.
pub mod collate;
/// Per-reference alignment coverage summaries.
pub mod coverage;
mod coverage_hts;
/// Padded-reference alignment projection.
pub mod depad;
/// Per-position alignment depth.
pub mod depth;
/// Name-grouped FASTA and FASTQ extraction.
pub mod fastx;
/// Mate-field repair for name-grouped alignments.
pub mod fixmate;
/// SAM flag parsing and rendering.
pub mod flags;
/// Alignment flag statistics.
pub mod flagstat;
/// Header and leading-record inspection.
pub mod head;
/// Per-reference index or stream statistics.
pub mod idxstats;
/// FASTQ conversion to unmapped alignment records.
pub mod import;
/// Alignment random-access index construction.
pub mod index;
/// Coordinate-sorted duplicate marking.
pub mod markdup;
/// Ordered alignment merging.
pub mod merge;
/// Per-position pileup generation.
pub mod mpileup;
/// Alignment-header program provenance.
mod program;
/// Alignment integrity checks.
pub mod quickcheck;
/// Compressed-block BAM header replacement.
pub mod reheader;
/// Read-group sample metadata.
pub mod samples;
/// Bounded-memory alignment sorting.
pub mod sort;
/// Alignment filtering and format conversion.
pub mod view;

pub use program::Program;

mod alignment_order;
mod alignment_stream;
mod amplicon;
mod bgzf_rewrite;
mod cli;
mod commands;
mod filter;
mod header_merge;
mod header_source;
mod hts_metadata;
mod hts_quickcheck;
mod input;
mod md;
mod output;

#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
