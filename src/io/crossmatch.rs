//! Crossmatch-style pairwise output, re-exported from `aln_core::crossmatch`.
//!
//! Moved to dfam-lib because that is where the *writer* already lived
//! (`aln_core::fmt::to_crossmatch`); a format the library could emit but not
//! parse meant any other consumer had to depend on a curation tool to read it
//! back.

pub use aln_core::crossmatch::*;
