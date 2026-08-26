//! Tile first, then vote: the unified MHSP selection for both phases.
//!
//! The shipping pipeline selects HSPs in three steps — cull at 80% overlap,
//! pick a strand by summed score, then trim the winner into a tiling — and the
//! order makes two things awkward.
//!
//! * The **cull is largely redundant with the trim**. When one HSP contains
//!   another, both drop it. They diverge only in a narrow band (the trim keeps
//!   a remainder the cull would have deleted whole), and where they diverge the
//!   trim is the better rule: it keeps the novel bases as evidence rather than
//!   discarding an HSP because most of it was already claimed, and it thresholds
//!   in absolute bp (`min_row_len`) rather than as a percentage.
//! * The **vote reads pre-trim scores**. Strand selection sums HSP scores as
//!   the aligner reported them, before any trimming, so a strand can win on
//!   evidence that never becomes a row.
//!
//! This module inverts the order: trim each strand into a tiling, then sum what
//! actually became rows, then keep the better tiling. The vote then reflects
//! real evidence, and the cull is unnecessary because the trim subsumes it.
//!
//! It also applies to **both phases**. The shipping pipeline trims in phase 2
//! only, so the bootstrap consensus is called from an alignment that can
//! double-count an instance's bases without limit — phase 1 measures overlap on
//! the *reference* axis, so nothing constrains instance-axis reuse at all. That
//! makes the bootstrap MSA qualitatively different from every refinement MSA
//! after it, for no stated reason.

use crate::AlignmentSource;
use aln_core::{Alignment, Sequence, Strand};
use aln_engine::Result;

/// Trim each strand's HSPs into a tiling, then keep the higher-scoring tiling.
///
/// Wraps the raw engine directly — there is no cull underneath, which is the
/// point. `min_row_len` is the shortest row worth keeping after trimming.
pub struct TileVoteFilter<S> {
    inner: S,
    min_row_len: usize,
    /// Choose a strand per instance, in phase 1. Off leaves both tilings.
    best_orientation_phase1: bool,
    /// The same, in phase 2. Separate so the phases can be measured apart.
    best_orientation_phase2: bool,
    /// Tile-then-vote in phase 1 as well as phase 2.
    ///
    /// Off leaves phase 1 to whatever is underneath, which is how the phase-2
    /// half is measured on its own. Phase 1 is the half most likely to move
    /// results, because trimming there changes which candidate scores highest
    /// and so which reference seeds the whole run.
    tile_phase1: bool,
}

impl<S> TileVoteFilter<S> {
    pub fn new(inner: S, min_row_len: usize, phase1: bool, phase2: bool) -> Self {
        TileVoteFilter {
            inner,
            min_row_len,
            best_orientation_phase1: phase1,
            best_orientation_phase2: phase2,
            tile_phase1: true,
        }
    }

    /// Apply tile-then-vote in phase 2 only.
    pub fn phase2_only(mut self) -> Self {
        self.tile_phase1 = false;
        self
    }
}

/// Per instance: tile each strand, sum the tiled scores, keep the better set.
fn tile_then_vote(
    hits: Vec<(usize, Alignment)>,
    min_row_len: usize,
    best_orientation: bool,
) -> Vec<(usize, Alignment)> {
    let mut by: std::collections::HashMap<usize, Vec<Alignment>> =
        std::collections::HashMap::new();
    let mut order: Vec<usize> = Vec::new();
    for (i, a) in hits {
        if !by.contains_key(&i) {
            order.push(i);
        }
        by.entry(i).or_default().push(a);
    }

    let mut out = Vec::new();
    for i in order {
        let group = by.remove(&i).unwrap_or_default();
        let (fwd, rev): (Vec<Alignment>, Vec<Alignment>) =
            group.into_iter().partition(|a| a.strand == Strand::Plus);

        // Tile each orientation independently. `trim_reused_instance` rescales
        // a trimmed alignment's score in proportion to what survived, so the
        // sums below are of evidence that actually becomes rows.
        let fwd = crate::trim::trim_reused_instance(fwd, min_row_len);
        let rev = crate::trim::trim_reused_instance(rev, min_row_len);

        let keep: Vec<Alignment> = if !best_orientation {
            let mut both = fwd;
            both.extend(rev);
            both
        } else {
            let fs: i64 = fwd.iter().map(|a| a.score as i64).sum();
            let rs: i64 = rev.iter().map(|a| a.score as i64).sum();
            // Forward wins ties, as the shipping orientation rule does.
            if rs > fs {
                rev
            } else {
                fwd
            }
        };
        for a in keep {
            out.push((i, a));
        }
    }
    out
}

impl<S: AlignmentSource> AlignmentSource for TileVoteFilter<S> {
    fn against_reference(
        &self,
        reference: &Sequence,
        seqs: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>> {
        let hits = self.inner.against_reference(reference, seqs, skip)?;
        Ok(tile_then_vote(hits, self.min_row_len, self.best_orientation_phase2))
    }

    fn against_every_reference(
        &self,
        seqs: &[Sequence],
    ) -> Result<Vec<Vec<(usize, Alignment)>>> {
        let per = self.inner.against_every_reference(seqs)?;
        if !self.tile_phase1 {
            return Ok(per);
        }
        Ok(per
            .into_iter()
            .map(|h| tile_then_vote(h, self.min_row_len, self.best_orientation_phase1))
            .collect())
    }

    fn against_reference_nested(
        &self,
        reference: &Sequence,
        seqs: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>> {
        let hits = self.inner.against_reference_nested(reference, seqs, skip)?;
        if !self.tile_phase1 {
            return Ok(hits);
        }
        Ok(tile_then_vote(hits, self.min_row_len, self.best_orientation_phase1))
    }

    fn batches_all_vs_all(&self) -> bool {
        self.inner.batches_all_vs_all()
    }
}
