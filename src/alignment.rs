use std::fmt;
use crate::matrix;

/// Strand orientation of an aligned sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Orientation {
    #[default]
    Forward,
    Reverse,
}

impl fmt::Display for Orientation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Orientation::Forward => write!(f, "+"),
            Orientation::Reverse => write!(f, "-"),
        }
    }
}

/// One entry in the multiple alignment — either the reference (index 0) or an
/// aligned instance.
///
/// All byte slices are stored in *gapped* form (gaps = b'-'), padded so every
/// sequence has the same total length as the reference.  Sequences that begin
/// after the left edge or end before the right edge are left-padded / right-
/// padded with spaces (b' ') rather than gaps so they can be distinguished from
/// internal deletions.
#[derive(Debug, Clone)]
pub struct SequenceRow {
    /// Sequence identifier (e.g. "chr1:12345-12500(+)").
    pub name: String,

    /// Gapped sequence bytes.  Length == alignment width.
    /// Interior gaps = b'-', flanking padding = b' '.
    pub seq: Vec<u8>,

    /// 0-based start position in the gapped alignment (first non-space column).
    pub start: usize,

    /// 0-based end position in the gapped alignment (last non-space column, inclusive).
    pub end: usize,

    /// 1-based, fully-closed start in the original (ungapped) source sequence.
    pub seq_start: u64,

    /// 1-based, fully-closed end in the original (ungapped) source sequence.
    pub seq_end: u64,

    pub orient: Orientation,

    /// Optional left-flanking genomic sequence (not part of the alignment).
    pub lf_seq: Option<Vec<u8>>,

    /// Optional right-flanking genomic sequence.
    pub rf_seq: Option<Vec<u8>>,

    /// GC fraction of the source region (0.0–1.0).
    pub gc_background: Option<f64>,

    /// Raw divergence from reference.
    pub div: Option<f64>,

    /// Kimura divergence.
    pub kdiv: Option<f64>,

    /// Transition count vs. reference.
    pub trans_i: Option<u32>,

    /// Transversion count vs. reference.
    pub trans_v: Option<u32>,

    /// Source divergence from original pairwise alignment data.
    pub src_div: Option<f64>,
}

impl SequenceRow {
    /// Construct a minimal entry from a name and a gapped sequence.
    pub fn new(name: impl Into<String>, seq: Vec<u8>) -> Self {
        let start = seq.iter().position(|&b| b != b' ').unwrap_or(0);
        let end = seq.iter().rposition(|&b| b != b' ').unwrap_or(0);
        SequenceRow {
            name: name.into(),
            seq,
            start,
            end,
            seq_start: 0,
            seq_end: 0,
            orient: Orientation::Forward,
            lf_seq: None,
            rf_seq: None,
            gc_background: None,
            div: None,
            kdiv: None,
            trans_i: None,
            trans_v: None,
            src_div: None,
        }
    }

    /// Ungapped sequence length (excludes gap and space characters).
    pub fn ungapped_len(&self) -> usize {
        self.seq.iter().filter(|&&b| b != b'-' && b != b' ').count()
    }
}

/// The multiple sequence alignment.
///
/// `sequences[0]` is the **reference** sequence (consensus or representative).
/// `sequences[1..]` are the aligned instances.
///
/// All sequences share the same `width()`.
#[derive(Debug, Clone)]
pub struct MultiAlign {
    /// Reference followed by aligned instances.
    pub sequences: Vec<SequenceRow>,

    /// Total gapped width (cached; equals sequences[0].seq.len()).
    width: usize,

    /// Cached consensus from the last call to consensus() / set_consensus().
    cached_consensus: Option<Vec<u8>>,

    /// Low-quality alignment blocks: list of (start, end) in gapped coords.
    pub low_quality_blocks: Vec<(usize, usize)>,
}

impl MultiAlign {
    /// Create an empty alignment.
    pub fn new() -> Self {
        MultiAlign {
            sequences: Vec::new(),
            width: 0,
            cached_consensus: None,
            low_quality_blocks: Vec::new(),
        }
    }

    /// Build from a reference sequence and a set of aligned instance sequences.
    ///
    /// All sequences must already be gapped and share the same length.
    pub fn from_sequences(reference: SequenceRow, instances: Vec<SequenceRow>) -> Self {
        let width = reference.seq.len();
        let mut sequences = Vec::with_capacity(1 + instances.len());
        sequences.push(reference);
        sequences.extend(instances);
        MultiAlign {
            sequences,
            width,
            cached_consensus: None,
            low_quality_blocks: Vec::new(),
        }
    }

    // ── Dimensions ────────────────────────────────────────────────────────────

    /// Gapped width of the alignment (number of columns).
    pub fn width(&self) -> usize {
        self.width
    }

    /// Number of aligned instances (does **not** count the reference).
    pub fn num_instances(&self) -> usize {
        self.sequences.len().saturating_sub(1)
    }

    /// Reference sequence entry (index 0).
    pub fn reference(&self) -> Option<&SequenceRow> {
        self.sequences.first()
    }

    /// Reference sequence bytes.
    pub fn reference_seq(&self) -> Option<&[u8]> {
        self.sequences.first().map(|s| s.seq.as_slice())
    }

    /// Instance at 0-based index `i` (among instances only; 0 = first non-reference).
    pub fn instance(&self, i: usize) -> Option<&SequenceRow> {
        self.sequences.get(i + 1)
    }

    // ── Profile ───────────────────────────────────────────────────────────────

    /// Build a position-frequency profile over all alignment columns.
    ///
    /// `profile[col][base_idx]` = count of that IUPAC character at that column.
    /// If `include_reference` is true, the reference row is counted; otherwise
    /// only instances are counted.
    pub fn build_profile(&self, include_reference: bool) -> Vec<[u32; matrix::ALPHA_LEN]> {
        let mut profile = vec![[0u32; matrix::ALPHA_LEN]; self.width];
        let start_idx = if include_reference { 0 } else { 1 };
        for seq in &self.sequences[start_idx..] {
            for (col, &b) in seq.seq.iter().enumerate() {
                if b == b' ' {
                    continue; // flanking padding — not part of this sequence
                }
                if let Some(idx) = matrix::alpha_idx(b.to_ascii_uppercase()) {
                    profile[col][idx] += 1;
                }
            }
        }
        profile
    }

    /// Per-column coverage depth (number of sequences spanning each column).
    pub fn coverage(&self) -> Vec<u32> {
        self.build_profile(false)
            .iter()
            .map(|col| col.iter().sum())
            .collect()
    }

    // ── Consensus ─────────────────────────────────────────────────────────────

    /// Return the cached consensus, or compute and cache it.
    ///
    /// Equivalent to the Perl `consensus(inclRef => 0)` call.
    pub fn consensus(&mut self, params: &crate::consensus::ConsensusParams) -> &[u8] {
        if self.cached_consensus.is_none() {
            let start = if params.include_reference { 0 } else { 1 };
            let raw: Vec<&[u8]> = self.sequences[start..]
                .iter()
                .map(|s| s.seq.as_slice())
                .collect();
            let cons = crate::consensus::build_consensus_from_sequences(&raw, params);
            self.cached_consensus = Some(cons);
        }
        self.cached_consensus.as_deref().unwrap()
    }

    /// Invalidate the cached consensus (call after modifying the alignment).
    pub fn invalidate_consensus(&mut self) {
        self.cached_consensus = None;
    }

    // ── Editing ───────────────────────────────────────────────────────────────

    /// Trim `left_bp` ungapped consensus positions from the left and `right_bp`
    /// from the right, then drop any instance sequences that become empty.
    ///
    /// Passing 0 for either side skips that trim.
    pub fn trim(&mut self, left_bp: usize, right_bp: usize) {
        if left_bp == 0 && right_bp == 0 {
            return;
        }
        // Determine gapped column boundaries by walking the reference.
        let ref_seq = match self.sequences.first() {
            Some(s) => s.seq.clone(),
            None => return,
        };

        let left_col = ungapped_to_gapped_col(&ref_seq, left_bp);
        let right_col = if right_bp == 0 {
            self.width
        } else {
            let ungapped_len = ref_seq.iter().filter(|&&b| b != b'-' && b != b' ').count();
            ungapped_to_gapped_col(&ref_seq, ungapped_len.saturating_sub(right_bp))
        };

        self.slice_columns(left_col, right_col);
    }

    /// Extract a sub-alignment from gapped column `col_start` (inclusive) to
    /// `col_end` (exclusive), dropping sequences that are entirely outside that
    /// range.
    pub fn slice_columns(&mut self, col_start: usize, col_end: usize) {
        let col_end = col_end.min(self.width);
        for seq in &mut self.sequences {
            seq.seq = seq.seq[col_start..col_end].to_vec();
            // Recalculate start/end within the new slice.
            seq.start = seq.seq.iter().position(|&b| b != b' ').unwrap_or(0);
            seq.end = seq.seq.iter().rposition(|&b| b != b' ').unwrap_or(0);
        }
        self.width = col_end - col_start;
        // Drop instances that became all-spaces.
        self.sequences.retain(|s| {
            // Always keep the reference (index 0 by position; we keep it unless it is also empty).
            s.seq.iter().any(|&b| b != b' ')
        });
        self.invalidate_consensus();
    }

    /// Reverse-complement the entire alignment in place.
    pub fn reverse_complement(&mut self) {
        for seq in &mut self.sequences {
            rc_in_place(&mut seq.seq);
            seq.orient = match seq.orient {
                Orientation::Forward => Orientation::Reverse,
                Orientation::Reverse => Orientation::Forward,
            };
            std::mem::swap(&mut seq.lf_seq, &mut seq.rf_seq);
            if let Some(ref mut s) = seq.lf_seq {
                rc_in_place(s);
            }
            if let Some(ref mut s) = seq.rf_seq {
                rc_in_place(s);
            }
        }
        self.invalidate_consensus();
    }
}

impl Default for MultiAlign {
    fn default() -> Self {
        Self::new()
    }
}

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Convert an ungapped position count to the corresponding gapped column index
/// in `seq`.  Returns `seq.len()` if `bp` exceeds the ungapped length.
fn ungapped_to_gapped_col(seq: &[u8], bp: usize) -> usize {
    let mut ungapped = 0usize;
    for (col, &b) in seq.iter().enumerate() {
        if b != b'-' && b != b' ' {
            if ungapped == bp {
                return col;
            }
            ungapped += 1;
        }
    }
    seq.len()
}

/// IUPAC reverse complement lookup table.
static RC_TABLE: [u8; 256] = {
    let mut t = [b'N'; 256];
    // Complement map (upper-case only; lower-case handled by to_ascii_uppercase before lookup).
    const PAIRS: &[(u8, u8)] = &[
        (b'A', b'T'), (b'T', b'A'), (b'G', b'C'), (b'C', b'G'),
        (b'R', b'Y'), (b'Y', b'R'), (b'K', b'M'), (b'M', b'K'),
        (b'S', b'S'), (b'W', b'W'), (b'B', b'V'), (b'V', b'B'),
        (b'D', b'H'), (b'H', b'D'), (b'N', b'N'), (b'X', b'X'),
        (b'-', b'-'), (b' ', b' '),
    ];
    let mut i = 0;
    while i < PAIRS.len() {
        let (a, b) = PAIRS[i];
        t[a as usize] = b;
        t[a.to_ascii_lowercase() as usize] = b;
        i += 1;
    }
    t
};

fn complement(b: u8) -> u8 {
    RC_TABLE[b as usize]
}

/// Reverse-complement `seq` in-place.
fn rc_in_place(seq: &mut Vec<u8>) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = complement(*b);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rc_simple() {
        let mut s = b"ACGT".to_vec();
        rc_in_place(&mut s);
        assert_eq!(&s, b"ACGT");
    }

    #[test]
    fn rc_with_gaps() {
        let mut s = b"AC-GT".to_vec();
        rc_in_place(&mut s);
        assert_eq!(&s, b"AC-GT");
    }

    #[test]
    fn ungapped_to_gapped() {
        let seq = b"AC-GT";
        assert_eq!(ungapped_to_gapped_col(seq, 0), 0);
        assert_eq!(ungapped_to_gapped_col(seq, 2), 3); // skip the gap at col 2
        assert_eq!(ungapped_to_gapped_col(seq, 4), 5); // past end
    }

    #[test]
    fn coverage_basic() {
        // Reference: ACG, two instances covering different ranges
        let ref_seq = SequenceRow::new("ref", b"ACG".to_vec());
        let inst1 = SequenceRow::new("s1", b"ACG".to_vec());
        let inst2 = SequenceRow::new("s2", b"A  ".to_vec()); // only column 0
        let msa = MultiAlign::from_sequences(ref_seq, vec![inst1, inst2]);
        let cov = msa.coverage();
        assert_eq!(cov[0], 2); // both instances at col 0
        assert_eq!(cov[1], 1);
        assert_eq!(cov[2], 1);
    }
}
