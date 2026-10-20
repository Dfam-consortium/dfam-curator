//! Compare a called consensus with the library consensus.
//!
//! We align the two sequences globally, so bases the called consensus has
//! beyond the library's ends land in terminal gaps instead of disappearing.
//! Those ends are then reported on their own, because they mean different
//! things: bases past an end extend the family, while missing bases at an end
//! may only mean no copy reached that far.

use aln_core::fmt;
use aln_core::matrix::SubstMatrix;
use aln_core::result::SearchResult;
use aln_core::seq::Sequence;
use aln_engine::params::{AlignMode, AlignParams};
use aln_engine::PairwiseAligner;
use aln_parasail::ParasailAligner;
use anyhow::Result;

#[derive(Debug, Default, Clone)]
pub struct ConsCmp {
    /// Positions where both have a base and the bases differ.
    pub sub: usize,
    /// Positions where the called base is not A, C, G or T and a copy was
    /// there, so the column really is ambiguous.
    pub amb: usize,
    /// Positions where the called base is not A, C, G or T because no copy
    /// covered the column.  These are missing data, not disagreement.
    pub unoccupied: usize,
    /// Internal bases only the called consensus has.
    pub ins: usize,
    /// Internal bases only the library consensus has.
    pub del: usize,
    /// Of those, the ones no sampled copy covered, which are missing data
    /// rather than a call the copies support.
    pub del_uncovered: usize,
    /// Bases the called consensus adds beyond the library's 5' and 3' ends.
    pub ext5: usize,
    pub ext3: usize,
    /// Library bases the called consensus lacks at its 5' and 3' ends.
    pub trim5: usize,
    pub trim3: usize,
    /// Of those, the ones no sampled copy covered.
    pub trim5_uncovered: usize,
    pub trim3_uncovered: usize,
    pub called_len: usize,
    pub library_len: usize,
    /// Positions compared, i.e. where both have a base.
    pub aligned: usize,
    /// Longest run of consecutive unoccupied reference positions.
    pub longest_unoccupied: usize,
    /// Unoccupied stretches in reference coordinates, 1-based inclusive.
    pub unoccupied_ranges: Vec<(usize, usize)>,
    /// Crossmatch-style alignment string, as ConsCAF carries it.
    pub caf: String,
    /// Library positions (0-based) where the called consensus carries a
    /// different unambiguous base, as (position, library base, called base).
    pub disputed: Vec<(usize, u8, u8)>,
    /// Library positions where both consensi carry the same unambiguous base.
    pub agreed: Vec<u32>,
}

impl ConsCmp {
    pub fn identity(&self) -> f64 {
        if self.aligned == 0 {
            return f64::NAN;
        }
        100.0 * (self.aligned - self.sub - self.amb - self.unoccupied) as f64 / self.aligned as f64
    }
    /// Identity over the columns a copy actually covered, which is what the
    /// sample can speak to.  Unoccupied columns are left out of both terms.
    pub fn identity_occupied(&self) -> f64 {
        let occupied = self.aligned - self.unoccupied;
        if occupied == 0 {
            return f64::NAN;
        }
        100.0 * (occupied - self.sub - self.amb) as f64 / occupied as f64
    }

    /// Called length minus library length, which the fields must account for.
    pub fn change(&self) -> i64 {
        self.called_len as i64 - self.library_len as i64
    }
}

/// A percentage rounded to one decimal, which is what the Perl writes.
fn pct(n: usize, total: usize) -> f64 {
    if total == 0 {
        return 0.0;
    }
    (1000.0 * n as f64 / total as f64).round() / 10.0
}

fn is_acgt(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// Align `called` to `library` and describe the difference.
///
/// `covered` marks, per library position, whether any sampled copy reached it.
pub fn compare(
    aligner: &ParasailAligner,
    library: &[u8],
    called: &[u8],
    covered: &[bool],
) -> Result<ConsCmp> {
    let mut out = ConsCmp {
        called_len: called.len(),
        library_len: library.len(),
        ..Default::default()
    };
    if called.is_empty() || library.is_empty() {
        return Ok(out);
    }
    // generateSeedAlignments runs parasail with the library consensus as the
    // query and the called consensus as the database.  The lineup matrix is
    // asymmetric, so the two orders score differently and can place gaps
    // differently; this keeps the same order.
    // generateSeedAlignments names the rows OLD and NEW in its ConsCAF line.
    let query = Sequence::new("OLD", library.to_vec());
    let subject = Sequence::new("NEW", called.to_vec());
    let Some(aln) = aligner
        .align(&query, &subject)
        .map_err(|e| anyhow::anyhow!("parasail: {e}"))?
    else {
        return Ok(out);
    };
    let (gapped_lib, gapped_called) = aln
        .gapped(&query.seq, &subject.seq)
        .map_err(|e| anyhow::anyhow!("gapped alignment: {e}"))?;
    let (gapped_subj, gapped_query) = (gapped_lib, gapped_called);
    let caf_result = aln.clone();

    // Terminal runs first: a gap in either row at the very start or end.
    let n = gapped_subj.len();
    let lead = (0..n)
        .take_while(|&i| gapped_subj[i] == b'-' || gapped_query[i] == b'-')
        .count();
    let trail = (0..n)
        .rev()
        .take_while(|&i| gapped_subj[i] == b'-' || gapped_query[i] == b'-')
        .count()
        .min(n - lead);

    let mut lib_pos = 0usize; // 0-based index into library
    for i in 0..n {
        let s = gapped_subj[i];
        let q = gapped_query[i];
        let terminal = i < lead || i >= n - trail;
        match (s == b'-', q == b'-') {
            (true, false) => {
                // Called consensus has a base the library lacks.
                if terminal {
                    if i < lead {
                        out.ext5 += 1;
                    } else {
                        out.ext3 += 1;
                    }
                } else {
                    out.ins += 1;
                }
            }
            (false, true) => {
                // Library base missing from the called consensus.
                let uncovered = covered.get(lib_pos).is_some_and(|c| !c);
                if terminal {
                    if i < lead {
                        out.trim5 += 1;
                        out.trim5_uncovered += usize::from(uncovered);
                    } else {
                        out.trim3 += 1;
                        out.trim3_uncovered += usize::from(uncovered);
                    }
                } else {
                    out.del += 1;
                    out.del_uncovered += usize::from(uncovered);
                }
                lib_pos += 1;
            }
            (false, false) => {
                out.aligned += 1;
                if !is_acgt(q) {
                    if covered.get(lib_pos).is_some_and(|c| *c) {
                        out.amb += 1;
                    } else {
                        out.unoccupied += 1;
                    }
                } else if !is_acgt(s) {
                    // A library ambiguity code is not counted as a substitution.
                } else if s.to_ascii_uppercase() != q.to_ascii_uppercase() {
                    out.sub += 1;
                    out.disputed.push((lib_pos, s.to_ascii_uppercase(), q.to_ascii_uppercase()));
                } else {
                    out.agreed.push(lib_pos as u32);
                }
                if s.to_ascii_uppercase() == q.to_ascii_uppercase() {
                } else {
                }
                lib_pos += 1;
            }
            (true, true) => {}
        }
    }
    // Unoccupied stretches, for the record's occupancy note.
    let mut run: Option<(usize, usize)> = None;
    for (i, c) in covered.iter().enumerate() {
        match (c, run) {
            (false, None) => run = Some((i + 1, i + 1)),
            (false, Some((a, _))) => run = Some((a, i + 1)),
            (true, Some(r)) => {
                out.unoccupied_ranges.push(r);
                run = None;
            }
            (true, None) => {}
        }
    }
    if let Some(r) = run {
        out.unoccupied_ranges.push(r);
    }
    // aln-core writes the CAF string, including the Perl's habit of leaving a
    // gap run open when the alignment ends in one.  The three percentages are
    // scaled as the Perl scales them: substitutions and deletions over the
    // library length, insertions over the called length.
    let mut sr = SearchResult::new(caf_result);
    let del_all = out.del + out.trim5 + out.trim3;
    sr.pct_diverge = pct(out.sub, out.library_len);
    sr.pct_delete = pct(del_all, out.library_len);
    sr.pct_insert = pct(out.ins + out.ext5 + out.ext3, out.called_len);
    out.caf = fmt::to_caf(&sr, &query.seq, &subject.seq).map_err(|e| anyhow::anyhow!("caf: {e}"))?;
    out.longest_unoccupied = out
        .unoccupied_ranges
        .iter()
        .map(|(a, b)| b - a + 1)
        .max()
        .unwrap_or(0);
    Ok(out)
}

/// A global aligner over the substitution matrix RepeatMasker used.
pub fn global_aligner(matrix: SubstMatrix, gap_open: u32, gap_extend: u32) -> Result<ParasailAligner> {
    let params = AlignParams {
        mode: AlignMode::Global,
        gap_open,
        gap_extend,
        min_score: i32::MIN / 2,
        traceback: true,
        ..Default::default()
    };
    ParasailAligner::new(matrix, params).map_err(|e| anyhow::anyhow!("parasail: {e}"))
}
