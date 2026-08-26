//! Co-linear chaining of one copy's alignments, with overlaps trimmed rather
//! than discarded.
//!
//! # Why not the masklevel rules
//!
//! Every filter in this codebase decides the fate of a whole HSP: keep it (and
//! its redundant part) or drop it (and its unique part). For a copy carrying a
//! tandem duplication the aligner emits `A+X` and `X'+B`; the second overlaps
//! the first on the consensus but continues into `B`, which nothing else
//! covers. Dropping it loses `B`; keeping it lets one copy cover `X` twice, and
//! double coverage is what distorts [`crate::lowqual::column_profile`] — it
//! averages per row, so a copy present twice votes twice on whether a column
//! looks healthy.
//!
//! Genome aligners resolve this by **trimming at a crossover**, not by
//! selection: a UCSC chain's blocks are strictly ordered and non-overlapping on
//! both axes, and overlaps between candidate blocks are cut. This module does
//! the same for one instance's HSPs.
//!
//! # What it produces
//!
//! Anchors that are strictly progressing on both axes, pairwise
//! non-overlapping after trimming. Off-chain anchors — the ones that cannot be
//! ordered consistently with the rest — are dropped, because a copy region that
//! maps out of order relative to its own neighbours is not evidence about the
//! consensus.
//!
//! # Coordinates
//!
//! An [`Alignment`]'s edit script runs in **query-forward** order and the
//! subject is reverse-complemented when the strand is minus, so walking the ops
//! advances the query from `query_start` upward and the subject either upward
//! (plus) or downward (minus). Chaining normalises the subject axis by
//! reflecting minus-strand spans, so "increasing on both axes" means the same
//! thing either way.

use aln_core::align::{Alignment, EditOp, EditScript};
#[cfg(test)]
use aln_core::seq::Strand;

/// Query bases a front trim would consume, without performing it.
///
/// The decision of whether to keep an overlap as an insertion depends on how
/// many copy bases the trim is about to take, and that is not the same as the
/// query-axis overlap: a tandem duplication has *no* query-axis overlap at all
/// — the two units are adjacent on the copy — yet trimming the subject-axis
/// overlap consumes the whole duplicated unit.
pub fn front_trim_query_cost(edits: &EditScript, n_subj: usize, n_query: usize) -> usize {
    let (mut got_s, mut got_q) = (0usize, 0usize);
    for (op, count) in edits.ops.iter().copied() {
        if got_s >= n_subj && got_q >= n_query {
            break;
        }
        let need = match op {
            EditOp::Sub => (n_subj.saturating_sub(got_s)).max(n_query.saturating_sub(got_q)),
            EditOp::GapInQuery => n_subj.saturating_sub(got_s),
            EditOp::GapInSubject => n_query.saturating_sub(got_q),
        };
        let take = need.min(count as usize);
        match op {
            EditOp::Sub => {
                got_s += take;
                got_q += take;
            }
            EditOp::GapInQuery => got_s += take,
            EditOp::GapInSubject => got_q += take,
        }
    }
    got_q
}

/// Drop the first `n_subj` subject bases and `n_query` query bases from the
/// front of an alignment, whichever requirement bites harder.
///
/// Returns the query bases actually removed, or `None` if the alignment is
/// consumed entirely. When `as_insertion` is set the removed query bases are
/// re-attached as a leading `GapInSubject` run, so the copy keeps them and they
/// occupy insertion columns.
pub fn trim_front(
    a: &mut Alignment,
    n_subj: usize,
    n_query: usize,
    as_insertion: bool,
) -> Option<usize> {
    let (mut got_s, mut got_q) = (0usize, 0usize);
    let mut ops: Vec<(EditOp, u32)> = Vec::new();
    let mut iter = a.edits.ops.iter().copied().peekable();

    while got_s < n_subj || got_q < n_query {
        let Some((op, count)) = iter.next() else { return None };
        let count = count as usize;
        // Which axes this op advances at all.
        let (adv_s, adv_q) = match op {
            EditOp::Sub => (true, true),
            EditOp::GapInQuery => (true, false),
            EditOp::GapInSubject => (false, true),
        };
        // Units of this run needed to satisfy what is still outstanding *on an
        // axis this op can advance*.
        let need = {
            let ns = if adv_s { n_subj.saturating_sub(got_s) } else { 0 };
            let nq = if adv_q { n_query.saturating_sub(got_q) } else { 0 };
            ns.max(nq)
        };
        if need == 0 {
            // This run advances nothing still outstanding — a leading
            // insertion when subject bases are wanted, say. It cannot be kept:
            // the alignment must start where the trim leaves it, so the whole
            // run goes and the walk continues.
            //
            // Treating this as a partial split instead was a bug: `take` came
            // out 0, the run was pushed back, and the loop broke before a
            // single subject base had been trimmed, returning an untrimmed
            // anchor that still overlapped its predecessor.
            if adv_s {
                got_s += count;
            }
            if adv_q {
                got_q += count;
            }
            continue;
        }
        let take = need.min(count);
        if adv_s {
            got_s += take;
        }
        if adv_q {
            got_q += take;
        }
        if count > take {
            ops.push((op, (count - take) as u32));
            // Only stop once both requirements are met; a subject-only run can
            // satisfy `n_subj` while `n_query` is still outstanding.
            if got_s >= n_subj && got_q >= n_query {
                break;
            }
        }
    }
    ops.extend(iter);
    // An alignment must not begin on a gap column.
    while matches!(ops.first(), Some((EditOp::GapInQuery, _)) | Some((EditOp::GapInSubject, _))) {
        let (op, n) = ops.remove(0);
        match op {
            EditOp::GapInQuery => got_s += n as usize,
            EditOp::GapInSubject => got_q += n as usize,
            EditOp::Sub => unreachable!(),
        }
    }
    if ops.is_empty() {
        return None;
    }

    let mut edits = EditScript::new();
    if as_insertion && got_q > 0 {
        edits.push(EditOp::GapInSubject, got_q as u32);
    }
    for (op, n) in ops {
        edits.push(op, n);
    }

    // Prorate the score by surviving columns. This is the fallback when the
    // caller has no matrix to rescore with
    // with a raw-score ratio when it does.
    let before = a.edits.align_len().max(1) as f64;
    let after = edits.align_len() as f64;
    a.score = ((a.score as f64) * (after / before)).round() as i32;

    a.edits = edits;
    if !as_insertion {
        a.query_start += got_q;
    }
    if a.strand.is_minus() {
        a.subj_end -= got_s.min(a.subj_end - a.subj_start);
    } else {
        a.subj_start += got_s;
    }
    Some(got_q)
}

#[cfg(test)]
mod tests {
    use super::*;
    use aln_core::seq;

    /// Build a plus-strand alignment from an op list.
    fn aln(qs: usize, ss: usize, ops: &[(EditOp, u32)], strand: Strand) -> Alignment {
        let mut edits = EditScript::new();
        for &(op, n) in ops {
            edits.push(op, n);
        }
        Alignment::new("q", "s", qs, ss, strand, 100, edits)
    }

    /// Trimming the front by subject bases must move `subj_start` and
    /// `query_start` by exactly what the ops consumed.
    #[test]
    fn trim_front_moves_both_axes() {
        let mut a = aln(10, 100, &[(EditOp::Sub, 20)], Strand::Plus);
        let removed = trim_front(&mut a, 5, 0, false).unwrap();
        assert_eq!(removed, 5, "5 subject bases of Sub consume 5 query bases");
        assert_eq!(a.subj_start, 105);
        assert_eq!(a.query_start, 15);
        assert_eq!(a.edits.align_len(), 15);
    }

    /// On the minus strand the subject is consumed from the far end, so the
    /// front trim must lower `subj_end` and leave `subj_start` alone.
    #[test]
    fn trim_front_on_minus_moves_the_far_end() {
        let mut a = aln(10, 100, &[(EditOp::Sub, 20)], Strand::Minus);
        trim_front(&mut a, 5, 0, false).unwrap();
        assert_eq!(a.subj_start, 100, "minus-strand start is untouched");
        assert_eq!(a.subj_end, 115);
        assert_eq!(a.query_start, 15);
    }

    /// A gap in the query advances only the subject, so trimming across it
    /// costs subject bases without costing query bases.
    #[test]
    fn trim_front_across_a_deletion() {
        let mut a = aln(
            0, 0,
            &[(EditOp::Sub, 4), (EditOp::GapInQuery, 3), (EditOp::Sub, 10)],
            Strand::Plus,
        );
        let removed = trim_front(&mut a, 7, 0, false).unwrap();
        assert_eq!(removed, 4, "only the 4 Sub columns consumed query bases");
        assert_eq!(a.subj_start, 7);
        assert_eq!(a.query_start, 4);
    }

    /// An alignment may not start on a gap column: trimming to the boundary of
    /// a deletion must eat the whole gap run.
    #[test]
    fn trim_front_never_leaves_a_leading_gap() {
        let mut a = aln(
            0, 0,
            &[(EditOp::Sub, 4), (EditOp::GapInQuery, 3), (EditOp::Sub, 10)],
            Strand::Plus,
        );
        trim_front(&mut a, 4, 0, false).unwrap();
        assert!(
            matches!(a.edits.ops.first(), Some((EditOp::Sub, _))),
            "leading gap run must be absorbed, got {:?}",
            a.edits.ops
        );
        assert_eq!(a.subj_start, 7, "the absorbed deletion counts as subject");
    }

    /// With `as_insertion` the copy keeps its bases: the query span is
    /// unchanged and the removed bases become insertion columns.
    #[test]
    fn trim_front_can_keep_the_overlap_as_an_insertion() {
        let mut a = aln(10, 100, &[(EditOp::Sub, 20)], Strand::Plus);
        let before = a.query_start;
        trim_front(&mut a, 5, 0, true).unwrap();
        assert_eq!(a.query_start, before, "no query bases are given up");
        assert_eq!(a.edits.ops.first(), Some(&(EditOp::GapInSubject, 5)));
        assert_eq!(a.edits.query_consumed(), 20, "all 20 query bases still there");
        assert_eq!(a.subj_start, 105);
    }

    /// Trimming past the end reports failure rather than producing an empty
    /// alignment.
    #[test]
    fn trim_front_past_the_end_is_none() {
        let mut a = aln(0, 0, &[(EditOp::Sub, 5)], Strand::Plus);
        assert!(trim_front(&mut a, 99, 0, false).is_none());
    }




    /// The lower-scoring row gives up the shared bases; both rows survive.
    #[test]
    fn reused_bases_are_trimmed_from_the_weaker_row() {
        let mut strong = aln(100, 0, &[(EditOp::Sub, 100)], Strand::Plus);
        strong.score = 900;
        // Shares genomic 180-200 with `strong`, and continues past it.
        let mut weak = aln(180, 500, &[(EditOp::Sub, 100)], Strand::Plus);
        weak.score = 400;
        let out = trim_reused_instance(vec![strong, weak], 10);
        assert_eq!(out.len(), 2, "both rows survive");
        assert!(
            out[0].query_end <= out[1].query_start,
            "no genomic base is claimed twice: {}-{} then {}-{}",
            out[0].query_start, out[0].query_end, out[1].query_start, out[1].query_end
        );
        assert!(out[1].query_end > 200, "the weaker row keeps its unique part");
    }

    /// A tandem expansion — one copy covering the same consensus region twice
    /// from different genomic bases — must be left alone. This is the case that
    /// separates this operation from chaining, which would trim it.
    #[test]
    fn reference_overlap_is_left_alone() {
        let mut first = aln(0, 100, &[(EditOp::Sub, 60)], Strand::Plus);
        first.score = 500;
        let mut second = aln(60, 100, &[(EditOp::Sub, 60)], Strand::Plus);
        second.score = 480;
        let out = trim_reused_instance(vec![first, second], 10);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].subj_start, out[1].subj_start, "both still cover the same consensus span");
        assert_eq!(out[0].query_end, out[1].query_start, "and neither lost a genomic base");
    }

    /// Nothing shared, nothing changed.
    #[test]
    fn disjoint_rows_are_untouched() {
        let a = aln(0, 0, &[(EditOp::Sub, 50)], Strand::Plus);
        let b = aln(200, 300, &[(EditOp::Sub, 50)], Strand::Plus);
        let out = trim_reused_instance(vec![a.clone(), b.clone()], 10);
        assert_eq!(out.len(), 2);
        assert_eq!((out[0].query_start, out[0].query_end), (a.query_start, a.query_end));
        assert_eq!((out[1].query_start, out[1].query_end), (b.query_start, b.query_end));
    }

    /// A row left too short by the trim is dropped rather than kept as a stub.
    #[test]
    fn a_row_trimmed_below_the_threshold_is_dropped() {
        let mut strong = aln(0, 0, &[(EditOp::Sub, 100)], Strand::Plus);
        strong.score = 900;
        let mut weak = aln(90, 500, &[(EditOp::Sub, 20)], Strand::Plus);
        weak.score = 100;
        let out = trim_reused_instance(vec![strong, weak], 25);
        assert_eq!(out.len(), 1, "only 10 unclaimed bases remained, under the floor");
    }

    /// Trimming the tail moves `query_end` and `subj_end`, not the starts.
    #[test]
    fn trim_back_moves_the_ends() {
        let mut a = aln(10, 100, &[(EditOp::Sub, 20)], Strand::Plus);
        trim_back(&mut a, 0, 5).unwrap();
        assert_eq!(a.query_start, 10);
        assert_eq!(a.query_end, 25);
        assert_eq!(a.subj_start, 100);
        assert_eq!(a.subj_end, 115);
    }

    /// A leading insertion must not stop the walk. Trimming subject bases from
    /// an anchor that opens with `GapInSubject` used to return the anchor
    /// untouched, because the run advanced no subject and was mistaken for a
    /// partial split — leaving two chain members covering the same reference
    /// span, the exact thing the trim exists to prevent.
    #[test]
    fn trim_front_walks_past_a_leading_insertion() {
        let mut a = aln(
            0, 100,
            &[(EditOp::GapInSubject, 6), (EditOp::Sub, 30)],
            Strand::Plus,
        );
        let removed = trim_front(&mut a, 3, 0, false).unwrap();
        assert_eq!(a.subj_start, 103, "3 subject bases really came off");
        assert_eq!(removed, 9, "6 inserted bases plus 3 aligned ones");
        assert!(matches!(a.edits.ops.first(), Some((EditOp::Sub, _))));
    }




}

/// Drop the last `n_subj` subject bases and `n_query` query bases from an
/// alignment. The mirror of [`trim_front`].
pub fn trim_back(a: &mut Alignment, n_subj: usize, n_query: usize) -> Option<usize> {
    let (mut got_s, mut got_q) = (0usize, 0usize);
    let mut ops: Vec<(EditOp, u32)> = Vec::new();
    let mut iter = a.edits.ops.iter().copied().rev();

    while got_s < n_subj || got_q < n_query {
        let Some((op, count)) = iter.next() else { return None };
        let count = count as usize;
        let (adv_s, adv_q) = match op {
            EditOp::Sub => (true, true),
            EditOp::GapInQuery => (true, false),
            EditOp::GapInSubject => (false, true),
        };
        let need = {
            let ns = if adv_s { n_subj.saturating_sub(got_s) } else { 0 };
            let nq = if adv_q { n_query.saturating_sub(got_q) } else { 0 };
            ns.max(nq)
        };
        if need == 0 {
            if adv_s { got_s += count; }
            if adv_q { got_q += count; }
            continue;
        }
        let take = need.min(count);
        if adv_s { got_s += take; }
        if adv_q { got_q += take; }
        if count > take {
            ops.push((op, (count - take) as u32));
            if got_s >= n_subj && got_q >= n_query { break; }
        }
    }
    let mut rest: Vec<(EditOp, u32)> = iter.collect();
    rest.reverse();
    ops.reverse();
    rest.extend(ops);
    // An alignment may not end on a gap column.
    while matches!(rest.last(), Some((EditOp::GapInQuery, _)) | Some((EditOp::GapInSubject, _))) {
        let (op, n) = rest.pop().expect("non-empty");
        match op {
            EditOp::GapInQuery => got_s += n as usize,
            EditOp::GapInSubject => got_q += n as usize,
            EditOp::Sub => unreachable!(),
        }
    }
    if rest.is_empty() { return None; }
    let mut edits = EditScript::new();
    for (op, n) in rest { edits.push(op, n); }
    let before = a.edits.align_len().max(1) as f64;
    let after = edits.align_len() as f64;
    a.score = ((a.score as f64) * (after / before)).round() as i32;
    a.edits = edits;
    a.query_end -= got_q.min(a.query_end - a.query_start);
    if a.strand.is_minus() {
        a.subj_start += got_s.min(a.subj_end - a.subj_start);
    } else {
        a.subj_end -= got_s.min(a.subj_end - a.subj_start);
    }
    Some(got_q)
}

/// Ensure no genomic base of a copy appears in more than one row.
///
/// An MSA row is a statement about homology: "these genomic bases align here".
/// The same bases in two rows claims they are homologous to two places at once,
/// which cannot be true of both. After the orientation rule, 99% of the
/// same-copy row pairs that survive are exactly this — the same genomic bases
/// mapped to different consensus regions — though the overlaps are small (mean
/// 6% of the shorter row, 91% under 25%).
///
/// Deliberately **not** touched: overlap on the *reference* axis. One copy
/// covering the same consensus region twice from *different* bases is a tandem
/// expansion, and two overlapping rows is a legitimate rendering of it. Folding
/// the second unit into a single row as inserted bases is a different algorithm
/// making a different claim, not a cleanup of this one.
///
/// The highest-scoring alignment claims its bases first; a lower-scoring one is
/// trimmed to the largest stretch of its own span nobody has claimed, and
/// dropped if fewer than `min_row_len` bases survive.
pub fn trim_reused_instance(mut hits: Vec<Alignment>, min_row_len: usize) -> Vec<Alignment> {
    if hits.len() < 2 {
        return hits;
    }
    hits.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.query_start.cmp(&b.query_start))
            .then(b.query_end.cmp(&a.query_end))
    });

    let mut claimed: Vec<(usize, usize)> = Vec::new();
    let mut out: Vec<Alignment> = Vec::new();
    for mut a in hits {
        let (qs, qe) = (a.query_start, a.query_end);
        // Largest sub-span of this alignment nobody has claimed. A claim lying
        // strictly inside would need the row split in two, which would create
        // the very thing being removed, so the longer free side wins instead.
        let mut free = (qs, qe);
        let mut best = 0usize;
        let mut cuts = vec![qs, qe];
        for &(cs, ce) in &claimed {
            if ce > qs && cs < qe {
                cuts.push(cs.clamp(qs, qe));
                cuts.push(ce.clamp(qs, qe));
            }
        }
        cuts.sort_unstable();
        cuts.dedup();
        for w in cuts.windows(2) {
            let (lo, hi) = (w[0], w[1]);
            if claimed.iter().any(|&(cs, ce)| lo < ce && cs < hi) {
                continue;
            }
            if hi - lo > best {
                best = hi - lo;
                free = (lo, hi);
            }
        }
        if best < min_row_len {
            continue;
        }
        if free.1 < qe && trim_back(&mut a, 0, qe - free.1).is_none() {
            continue;
        }
        if free.0 > qs && trim_front(&mut a, 0, free.0 - qs, false).is_none() {
            continue;
        }
        if a.query_end.saturating_sub(a.query_start) < min_row_len {
            continue;
        }
        claimed.push((a.query_start, a.query_end));
        claimed.sort_unstable();
        let mut merged: Vec<(usize, usize)> = Vec::with_capacity(claimed.len());
        for &(cs, ce) in &claimed {
            match merged.last_mut() {
                Some(last) if cs <= last.1 => last.1 = last.1.max(ce),
                _ => merged.push((cs, ce)),
            }
        }
        claimed = merged;
        out.push(a);
    }
    if std::env::var_os("TE_COMPOSER_TRIM_DEBUG").is_some() {
        for i in 0..out.len() {
            for j in (i + 1)..out.len() {
                if out[i].query_start < out[j].query_end
                    && out[j].query_start < out[i].query_end
                {
                    eprintln!(
                        "TRIMRESIDUAL\t{}\t{}-{}\t{}-{}",
                        out[i].query_name,
                        out[i].query_start, out[i].query_end,
                        out[j].query_start, out[j].query_end
                    );
                }
            }
        }
    }
    out.sort_by_key(|a| a.query_start);
    out
}
