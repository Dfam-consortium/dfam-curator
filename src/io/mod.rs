pub mod clustal;
pub mod crossmatch;
pub mod fasta;
pub mod linup_fmt;
pub mod stockholm;
pub mod twobit;

use crate::alignment::MultiAlign;
use crate::build::Reference;
use std::io;
use std::path::Path;

/// Detected input format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Stockholm,
    Fasta,
    /// Clustal ALN interleaved format.
    Clustal,
    /// Crossmatch / RepeatMasker `.align` pairwise format.
    Crossmatch,
}

/// Auto-detect the format by scanning the first non-blank meaningful line of the file.
///
/// Crossmatch `.align` files begin with a multi-line preamble (cross_match version,
/// matrix, run parameters) before the first integer score line, so we scan up to
/// 200 lines before giving up.
pub fn detect_format(path: &Path) -> io::Result<Format> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    let f = BufReader::new(File::open(path)?);
    for line in f.lines().take(200) {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with("# STOCKHOLM") {
            return Ok(Format::Stockholm);
        }
        if trimmed.to_ascii_uppercase().starts_with("CLUSTAL") {
            return Ok(Format::Clustal);
        }
        if trimmed.starts_with('>') {
            return Ok(Format::Fasta);
        }
        // Crossmatch score lines start with an integer score.
        if trimmed.chars().next().map(|c| c.is_ascii_digit()).unwrap_or(false) {
            return Ok(Format::Crossmatch);
        }
        // Non-matching lines are preamble — keep scanning.
    }
    Err(io::Error::new(io::ErrorKind::InvalidData, "unrecognised alignment format"))
}

/// Read a pre-built MSA from `path`, auto-detecting the format.
///
/// For crossmatch `.align` files, this performs transitive alignment
/// treating the **subject** as the reference (the typical RepeatMasker
/// convention where the repeat consensus is the subject).
pub fn read_alignment(path: &Path) -> io::Result<MultiAlign> {
    match detect_format(path)? {
        Format::Stockholm => stockholm::read(path),
        Format::Fasta     => fasta::read(path),
        Format::Clustal   => clustal::read(path),
        Format::Crossmatch => read_crossmatch_as_multialign(path),
    }
}

fn read_crossmatch_as_multialign(path: &Path) -> io::Result<MultiAlign> {
    let hits = crossmatch::read(path)?;
    if hits.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "crossmatch file is empty"));
    }
    crate::build::build_from_pairwise(&hits, Reference::Subject, None)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))
}
