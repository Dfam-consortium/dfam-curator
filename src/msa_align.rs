//! Map an externally-supplied sequence (typically a hand-built consensus) onto
//! an existing multiple alignment.
//!
//! # Why profile alignment rather than consensus alignment
//!
//! The obvious approach is to call the consensus, ungap it, align the
//! hand-built sequence to it with a plain pairwise aligner, then transfer the
//! result back onto the alignment columns.  That throws away the one piece of
//! information that matters most here: which columns are *already* gap columns.
//!
//! A real seed alignment is full of sparsely-occupied insert columns.  When the
//! hand-built sequence carries an insertion relative to the consensus, it should
//! land *in those existing columns* wherever there is room, not manufacture new
//! ones — otherwise mapping a hand-built consensus onto a seed inflates the
//! width and shreds the alignment's existing indel structure.
//!
//! Aligning directly against the column profile gets this for free.  Consuming
//! an all-gap column is scored by the substitution matrix's own gap row (`+3`
//! gap-vs-gap), so it is cheap; placing a residue into such a column costs only
//! `-6`, which is well below the affine cost of opening a brand-new column.  The
//! dynamic program therefore absorbs insertions into available gap columns on
//! its own, and widens the alignment only when it genuinely runs out of room.
//!
//! # Scoring model
//!
//! Three moves, over the query `Q` (ungapped) and the profile's `W` columns:
//!
//! | move | meaning | score |
//! |------|---------|-------|
//! | place | query residue occupies column `j` | depth-averaged `MATRIX[q][*]` over column `j` |
//! | skip | column `j` consumed, query gets `-` | depth-averaged `MATRIX['-'][*]` over column `j` |
//! | insert | query residue with no column — **widens the alignment** | affine `-(open + k·extend)` |
//!
//! Column scores are averaged over the column's depth so they stay on the same
//! scale as a single matrix entry regardless of how many instances the alignment
//! holds; that keeps one set of gap penalties meaningful across a 5-sequence and
//! a 5000-sequence seed.
//!
//! Only `insert` is affine.  `skip` deliberately is not: the cost of consuming a
//! column is already fully determined by that column's contents under the same
//! matrix the consensus caller uses, and stacking an affine penalty on top would
//! double-charge it.

use aln_core::consensus::{alpha_idx, ALPHA_LEN, MATRIX};

/// Index of the gap character in the scoring matrix.
const GAP_IDX: usize = 17;

/// Guard against pathological memory use.  One byte of traceback per cell, so
/// this caps the traceback matrix at 1 GiB.
const MAX_CELLS: u64 = 1 << 30;

// ── Parameters ────────────────────────────────────────────────────────────────

/// Tuning for [`align_to_profile`].
#[derive(Debug, Clone)]
pub struct AlignParams {
    /// Penalty (positive; subtracted) for opening a run of new columns.
    ///
    /// Calibrated against the scale of [`MATRIX`], whose entries reach `-17`
    /// for a bad mismatch.  Opening a column has to cost *more* than any single
    /// substitution, for two reasons:
    ///
    /// - It must exceed the `-6` cost of placing a residue into an existing
    ///   all-gap column, or the aligner widens the alignment instead of reusing
    ///   the room it already has.
    /// - With [`AlignParams::free_end_gaps`] on, a cheap insertion lets the
    ///   dynamic program discard a flanking column for free and shift the whole
    ///   mapping by one to dodge a single expensive mismatch — a drastic
    ///   restructuring bought for a trivial gain.
    ///
    /// A hand-built consensus is an *edit* of the called one, so the default is
    /// deliberately conservative: prefer to report a substitution rather than
    /// restructure columns.  Lower it when a genuine short indel must be forced.
    pub gap_open: f64,

    /// Per-residue penalty (positive; subtracted) for extending a run of new
    /// columns.  Charged on the first inserted column as well as later ones.
    pub gap_extend: f64,

    /// Allow the query to cover only part of the alignment's width without
    /// paying for the flanking columns it does not reach.
    ///
    /// The query is always consumed in full; only the profile gets free ends.
    /// Turn this off to force a fully global fit.
    pub free_end_gaps: bool,
}

impl Default for AlignParams {
    fn default() -> Self {
        AlignParams {
            gap_open: 30.0,
            gap_extend: 3.0,
            free_end_gaps: true,
        }
    }
}

// ── Result types ──────────────────────────────────────────────────────────────

/// One step of the mapping, in left-to-right alignment order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlnOp {
    /// Query residue (index into the ungapped query) occupies existing column `col`.
    Place { qi: usize, col: usize },
    /// Existing column `col` is consumed; the query row gets a gap there.
    Skip { col: usize },
    /// Query residue needs a new column inserted before existing column `before`.
    ///
    /// `before == width` means the new column goes after every existing column.
    Insert { qi: usize, before: usize },
}

/// The mapping of a query sequence onto a set of alignment columns.
#[derive(Debug, Clone)]
pub struct Mapping {
    /// Steps in left-to-right order.
    pub ops: Vec<AlnOp>,
    /// Alignment score under [`AlignParams`].
    pub score: f64,
    /// Width of the original profile.
    pub old_width: usize,
    /// Width after the required column insertions.
    pub new_width: usize,
}

impl Mapping {
    /// Column insertions required, as `(before_column, count)` pairs in
    /// ascending order of `before_column`.
    ///
    /// `before_column` is an index into the **original** alignment: the new
    /// columns are inserted immediately before it.  A `before_column` equal to
    /// the original width appends at the right edge.
    pub fn insertions(&self) -> Vec<(usize, usize)> {
        let mut out: Vec<(usize, usize)> = Vec::new();
        for op in &self.ops {
            if let AlnOp::Insert { before, .. } = *op {
                match out.last_mut() {
                    Some((b, n)) if *b == before => *n += 1,
                    _ => out.push((before, 1)),
                }
            }
        }
        out
    }

    /// The query rendered as an alignment row of [`Mapping::new_width`] columns.
    ///
    /// `query` must be the same ungapped sequence that was aligned.  `gap` is
    /// the gap character to emit (`b'-'` or, for Dfam-clean Stockholm, `b'.'`).
    pub fn query_row(&self, query: &[u8], gap: u8) -> Vec<u8> {
        let mut row = Vec::with_capacity(self.new_width);
        for op in &self.ops {
            match *op {
                AlnOp::Place { qi, .. } | AlnOp::Insert { qi, .. } => row.push(query[qi]),
                AlnOp::Skip { .. } => row.push(gap),
            }
        }
        // With free end gaps the ops stop at the last column the query reached;
        // pad the untouched flanks.
        debug_assert!(row.len() <= self.new_width);
        row
    }

    /// Number of columns the query row leaves untouched at the left edge.
    pub fn leading_skipped(&self) -> usize {
        match self.ops.first() {
            Some(AlnOp::Place { col, .. }) | Some(AlnOp::Skip { col }) => *col,
            Some(AlnOp::Insert { before, .. }) => *before,
            None => 0,
        }
    }
}

// ── Profile construction ──────────────────────────────────────────────────────

/// Build a per-column count profile from gapped alignment rows.
///
/// Only the region each row actually covers contributes: everything before its
/// first residue and after its last is flanking padding, not evidence of a
/// deletion, and is excluded.  Gap bytes *inside* the covered region are
/// counted as the gap symbol.  Explicit space bytes are always padding.
///
/// This matches the pre-processing in
/// [`aln_core::consensus::build_consensus_from_sequences`], and matching it is not
/// optional: a seed alignment is mostly short fragments, so counting their
/// flanks would make nearly every column look gap-dominated.  The aligner would
/// then skip the whole alignment and open a new column for every query residue.
pub fn profile_from_rows(rows: &[&[u8]]) -> Vec<[u32; ALPHA_LEN]> {
    let width = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    let mut profile = vec![[0u32; ALPHA_LEN]; width];
    for row in rows {
        let first = match row.iter().position(|b| b.is_ascii_alphabetic()) {
            Some(f) => f,
            None => continue, // row has no residues at all
        };
        let last = row
            .iter()
            .rposition(|b| b.is_ascii_alphabetic())
            .unwrap_or(first);
        for (col, &b) in row.iter().enumerate().take(last + 1).skip(first) {
            if b == b' ' {
                continue;
            }
            if let Some(idx) = alpha_idx(b.to_ascii_uppercase()) {
                profile[col][idx] += 1;
            }
        }
    }
    profile
}

/// Depth-averaged score of every alphabet symbol against every column.
///
/// `table[col][sym]` is the mean `MATRIX[sym][observed]` over the column's
/// occupants.  Empty columns score 0 for everything.
fn score_table(profile: &[[u32; ALPHA_LEN]]) -> Vec<[f64; ALPHA_LEN]> {
    profile
        .iter()
        .map(|col| {
            let depth: u32 = col.iter().sum();
            let mut row = [0.0f64; ALPHA_LEN];
            if depth == 0 {
                return row;
            }
            let inv = 1.0 / depth as f64;
            for (sym, slot) in row.iter_mut().enumerate() {
                let mut acc = 0i64;
                for obs in 0..ALPHA_LEN {
                    let c = col[obs];
                    if c > 0 {
                        acc += c as i64 * MATRIX[sym][obs] as i64;
                    }
                }
                *slot = acc as f64 * inv;
            }
            row
        })
        .collect()
}

// ── Dynamic program ───────────────────────────────────────────────────────────

// Traceback bit flags, one byte per cell.
const TB_A_IS_PLACE: u8 = 0b0001; // A[i][j] came from a place move (else skip)
const TB_A_FROM_I: u8 = 0b0010; // ...and its predecessor was in state I (else A)
const TB_I_IS_EXTEND: u8 = 0b0100; // I[i][j] extended a run (else opened one)

const NEG: f64 = f64::NEG_INFINITY;

/// Align `query` (ungapped) against a column `profile`.
///
/// The query is consumed in full.  Returns `None` only when `query` is empty.
///
/// Ties are broken deterministically — a place move is preferred over a skip,
/// and an already-open insertion run is preferred over a fresh one — and the
/// result is then canonicalised so that indels sit as far left as they can
/// without changing the score.  Identical inputs therefore always produce
/// identical column assignments, including inside homopolymer runs where the
/// placement is otherwise arbitrary.
pub fn align_to_profile(
    query: &[u8],
    profile: &[[u32; ALPHA_LEN]],
    params: &AlignParams,
) -> anyhow::Result<Mapping> {
    let m = query.len();
    let w = profile.len();
    if m == 0 {
        anyhow::bail!("cannot map an empty sequence onto an alignment");
    }

    let cells = (m as u64 + 1) * (w as u64 + 1);
    if cells > MAX_CELLS {
        anyhow::bail!(
            "alignment too large: {} query residues × {} columns exceeds the \
             {} cell limit",
            m,
            w,
            MAX_CELLS
        );
    }

    let table = score_table(profile);
    let qidx: Vec<usize> = query
        .iter()
        .map(|&b| alpha_idx(b.to_ascii_uppercase()).unwrap_or(10)) // unknown -> N
        .collect();

    let open = params.gap_open + params.gap_extend;
    let extend = params.gap_extend;

    // Score rows are streamed; only the traceback is materialised in full.
    let mut tb = vec![0u8; cells as usize];
    let mut a_prev = vec![NEG; w + 1];
    let mut i_prev = vec![NEG; w + 1];
    let mut a_cur = vec![NEG; w + 1];
    let mut i_cur = vec![NEG; w + 1];

    // Row 0: no query residues consumed, only column skips.
    a_prev[0] = 0.0;
    for j in 1..=w {
        a_prev[j] = if params.free_end_gaps {
            0.0
        } else {
            a_prev[j - 1] + table[j - 1][GAP_IDX]
        };
        // Row 0 can only have arrived by skipping, from state A.
        tb[j] = 0;
    }

    for i in 1..=m {
        let q = qidx[i - 1];
        let row_base = i * (w + 1);

        // Column 0: query residues consumed with no columns available yet —
        // reachable only as an insertion run.
        a_cur[0] = NEG;
        i_cur[0] = {
            let from_open = a_prev[0] - open;
            let from_ext = i_prev[0] - extend;
            if from_ext >= from_open {
                tb[row_base] |= TB_I_IS_EXTEND;
                from_ext
            } else {
                from_open
            }
        };

        for j in 1..=w {
            let cell = row_base + j;
            let mut flags = 0u8;

            // ── A: last move consumed column j-1 ──────────────────────────
            // place: query residue i-1 goes into column j-1
            let (place_pred, place_from_i) = better(a_prev[j - 1], i_prev[j - 1]);
            let place = place_pred + table[j - 1][q];

            // skip: column j-1 consumed, query gets a gap
            let (skip_pred, skip_from_i) = better(a_cur[j - 1], i_cur[j - 1]);
            let skip = skip_pred + table[j - 1][GAP_IDX];

            // Prefer `place` on a tie (deterministic; also biases toward
            // consuming query and column together rather than drifting).
            a_cur[j] = if place >= skip {
                flags |= TB_A_IS_PLACE;
                if place_from_i {
                    flags |= TB_A_FROM_I;
                }
                place
            } else {
                if skip_from_i {
                    flags |= TB_A_FROM_I;
                }
                skip
            };

            // ── I: last move consumed query residue i-1 with no column ────
            let from_open = a_prev[j] - open;
            let from_ext = i_prev[j] - extend;
            i_cur[j] = if from_ext >= from_open {
                flags |= TB_I_IS_EXTEND;
                from_ext
            } else {
                from_open
            };

            tb[cell] = flags;
        }

        std::mem::swap(&mut a_prev, &mut a_cur);
        std::mem::swap(&mut i_prev, &mut i_cur);
    }

    // a_prev / i_prev now hold row m.
    let (mut end_j, mut end_in_i, mut best) = (w, false, NEG);
    if params.free_end_gaps {
        // Trailing columns the query never reaches are free.
        for j in 0..=w {
            for (in_i, score) in [(false, a_prev[j]), (true, i_prev[j])] {
                if score > best {
                    best = score;
                    end_j = j;
                    end_in_i = in_i;
                }
            }
        }
    } else {
        let (score, in_i) = better(a_prev[w], i_prev[w]);
        best = score;
        end_in_i = in_i;
    }

    if !best.is_finite() {
        anyhow::bail!(
            "no valid alignment of {} residues onto {} columns",
            m,
            w
        );
    }

    // ── Traceback ────────────────────────────────────────────────────────────
    let mut ops: Vec<AlnOp> = Vec::with_capacity(m + w);
    let (mut i, mut j, mut in_i) = (m, end_j, end_in_i);
    while i > 0 || j > 0 {
        let flags = tb[i * (w + 1) + j];
        if in_i {
            // Query residue i-1 was inserted before column j.
            ops.push(AlnOp::Insert {
                qi: i - 1,
                before: j,
            });
            in_i = flags & TB_I_IS_EXTEND != 0;
            i -= 1;
        } else if j == 0 {
            // Only insertions can reach column 0 with residues left.
            debug_assert!(i > 0);
            in_i = true;
        } else if flags & TB_A_IS_PLACE != 0 {
            ops.push(AlnOp::Place {
                qi: i - 1,
                col: j - 1,
            });
            in_i = flags & TB_A_FROM_I != 0;
            i -= 1;
            j -= 1;
        } else {
            if i == 0 && params.free_end_gaps {
                // Leading columns were free; stop rather than emitting skips.
                break;
            }
            ops.push(AlnOp::Skip { col: j - 1 });
            in_i = flags & TB_A_FROM_I != 0;
            j -= 1;
        }
    }
    ops.reverse();

    canonicalise(&mut ops, &table, &qidx);

    let inserted: usize = ops
        .iter()
        .filter(|o| matches!(o, AlnOp::Insert { .. }))
        .count();

    Ok(Mapping {
        ops,
        score: best,
        old_width: w,
        new_width: w + inserted,
    })
}

/// Return the larger of `a` (state A) and `b` (state I), plus whether I won.
/// Ties go to A.
#[inline]
fn better(a: f64, b: f64) -> (f64, bool) {
    if b > a {
        (b, true)
    } else {
        (a, false)
    }
}

/// Left-align gaps: slide each gap as far left as it will go without lowering
/// the score.
///
/// The dynamic program is already optimal, so no strictly-improving swap can
/// exist; this resolves only the *ties* — most visibly inside homopolymer and
/// tandem runs, where `GG-` and `G-G` score identically and the DP's choice is
/// an artefact of its tie-break order.  Fixing a direction makes the output
/// reproducible across runs, and matches the left-alignment convention used for
/// indels everywhere else in genomics.
///
/// Runs of inserted columns are left-aligned by the same rule, in a second
/// pass: sliding a run one column left turns the residue placed just before it
/// into an inserted one and places the run's last residue instead, which is
/// score-neutral exactly when those two residues are equivalent.  This is what
/// keeps a hand-built insertion inside a homopolymer from landing at an
/// arbitrary offset within the run.
fn canonicalise(ops: &mut [AlnOp], table: &[[f64; ALPHA_LEN]], qidx: &[usize]) {
    if ops.len() < 2 {
        return;
    }
    left_align_skips(ops, table, qidx);
    left_align_inserts(ops, table, qidx);
}

/// Slide `Skip` ops left past `Place` ops while the score allows.
fn left_align_skips(ops: &mut [AlnOp], table: &[[f64; ALPHA_LEN]], qidx: &[usize]) {
    let mut k = 1;
    while k < ops.len() {
        // Look for  Place(j) , Skip(j+1)  and try to swap to  Skip(j) , Place(j+1),
        // which moves the gap one column to the left.
        let (left, right) = (ops[k - 1], ops[k]);
        if let (AlnOp::Place { qi, col: cp }, AlnOp::Skip { col: cs }) = (left, right) {
            if cs == cp + 1 {
                let q = qidx[qi];
                let before = table[cp][q] + table[cs][GAP_IDX];
                let after = table[cp][GAP_IDX] + table[cs][q];
                if after >= before {
                    ops[k - 1] = AlnOp::Skip { col: cp };
                    ops[k] = AlnOp::Place { qi, col: cs };
                    // The gap may be able to slide further left.
                    if k > 1 {
                        k -= 1;
                        continue;
                    }
                }
            }
        }
        k += 1;
    }
}

/// Slide runs of `Insert` ops left past the `Place` that precedes them.
///
/// A run `Place(p, col=b-1), Insert(p+1, b), …, Insert(p+k, b)` becomes
/// `Insert(p, b-1), …, Insert(p+k-1, b-1), Place(p+k, col=b-1)`.  Query indices
/// stay ascending and the run length is unchanged, so the only score difference
/// is which residue occupies column `b-1`.
fn left_align_inserts(ops: &mut [AlnOp], table: &[[f64; ALPHA_LEN]], qidx: &[usize]) {
    let mut i = 0usize;
    while i < ops.len() {
        // Find a maximal run of Insert ops at ops[i..end].
        if !matches!(ops[i], AlnOp::Insert { .. }) {
            i += 1;
            continue;
        }
        let mut start = i;
        let before = match ops[i] {
            AlnOp::Insert { before, .. } => before,
            _ => unreachable!(),
        };
        let mut end = start;
        while end < ops.len() {
            match ops[end] {
                AlnOp::Insert { before: b, .. } if b == before => end += 1,
                _ => break,
            }
        }

        // Try to slide the whole run left, repeatedly.
        loop {
            if start == 0 {
                break;
            }
            let (prev_qi, prev_col) = match ops[start - 1] {
                AlnOp::Place { qi, col } => (qi, col),
                _ => break, // only a Place can be traded with the run
            };
            if prev_col + 1 != before_col(ops[start]) {
                break;
            }
            let last_qi = match ops[end - 1] {
                AlnOp::Insert { qi, .. } => qi,
                _ => unreachable!(),
            };
            // Column prev_col currently holds prev_qi; after the slide it holds last_qi.
            if table[prev_col][qidx[last_qi]] < table[prev_col][qidx[prev_qi]] {
                break;
            }
            // Rewrite the window ops[start-1 ..= end-1].
            let qis: Vec<usize> = std::iter::once(prev_qi)
                .chain((start..end).map(|k| match ops[k] {
                    AlnOp::Insert { qi, .. } => qi,
                    _ => unreachable!(),
                }))
                .collect();
            let n = qis.len();
            for (slot, &qi) in (start - 1..end).zip(&qis[..n - 1]) {
                ops[slot] = AlnOp::Insert {
                    qi,
                    before: prev_col,
                };
            }
            ops[end - 1] = AlnOp::Place {
                qi: qis[n - 1],
                col: prev_col,
            };
            start -= 1;
            end -= 1;
        }
        i = end.max(start + 1);
    }
}

/// The `before` column of an `Insert` op.
#[inline]
fn before_col(op: AlnOp) -> usize {
    match op {
        AlnOp::Insert { before, .. } => before,
        _ => usize::MAX,
    }
}

// ── Mapping statistics ────────────────────────────────────────────────────────

/// Summary of how disruptive a mapping is, for sanity-checking before surgery.
#[derive(Debug, Clone)]
pub struct MappingStats {
    /// Query residues that landed in an existing column.
    pub placed: usize,
    /// Query residues that required a brand-new column.
    pub inserted: usize,
    /// Existing columns the query gapped over, within its span.
    pub skipped: usize,
    /// Of the placed residues, how many agree with the called consensus.
    pub matched: usize,
    /// `matched / placed` — 1.0 when the query reproduces the consensus exactly.
    pub identity: f64,
    /// New columns the surgery would add.
    pub columns_added: usize,
    /// `columns_added / old_width` — how much the alignment would grow.
    pub width_growth: f64,
    /// Existing columns spanned by the query (`placed + skipped`).
    pub span: usize,
    /// `span / old_width` — how much of the alignment the query covers.
    pub coverage: f64,
}

impl Mapping {
    /// Measure the mapping against the alignment's own called consensus.
    ///
    /// `consensus` must be the gapped consensus at the *original* width — the
    /// same one [`aln_core::consensus::build_consensus_from_sequences`] returns.
    /// Positions where the consensus is a gap count as mismatches, since a
    /// query residue sitting in a gap column agrees with nothing.
    pub fn stats(&self, query: &[u8], consensus: &[u8]) -> MappingStats {
        let (mut placed, mut inserted, mut skipped, mut matched) = (0, 0, 0, 0);
        for op in &self.ops {
            match *op {
                AlnOp::Place { qi, col } => {
                    placed += 1;
                    let c = consensus.get(col).copied().unwrap_or(b'-');
                    if c != b'-' && c != b'.' && c.eq_ignore_ascii_case(&query[qi]) {
                        matched += 1;
                    }
                }
                AlnOp::Insert { .. } => inserted += 1,
                AlnOp::Skip { .. } => skipped += 1,
            }
        }
        let span = placed + skipped;
        let denom = self.old_width.max(1) as f64;
        MappingStats {
            placed,
            inserted,
            skipped,
            matched,
            identity: if placed == 0 {
                0.0
            } else {
                matched as f64 / placed as f64
            },
            columns_added: self.new_width - self.old_width,
            width_growth: (self.new_width - self.old_width) as f64 / denom,
            span,
            coverage: span as f64 / denom,
        }
    }
}

// ── Column surgery ────────────────────────────────────────────────────────────

/// Splice new columns into a single alignment row.
///
/// `insertions` is the `(before_column, count)` list from
/// [`Mapping::insertions`], in ascending order.  Inserted columns are filled
/// with `gap` when they fall inside the row's occupied span and with `pad`
/// outside it, preserving the convention that flanking padding is
/// distinguishable from an internal deletion.
///
/// Pass `pad == gap` for formats such as Dfam Stockholm that make no such
/// distinction.
pub fn splice_row(row: &[u8], insertions: &[(usize, usize)], gap: u8, pad: u8) -> Vec<u8> {
    let extra: usize = insertions.iter().map(|(_, n)| n).sum();
    if extra == 0 {
        return row.to_vec();
    }
    let first = row.iter().position(|&b| b != b' ').unwrap_or(row.len());
    let last = row.iter().rposition(|&b| b != b' ');

    let mut out = Vec::with_capacity(row.len() + extra);
    let mut next = 0usize;
    for (before, count) in insertions.iter().copied() {
        let before = before.min(row.len());
        out.extend_from_slice(&row[next..before]);
        // Inside the occupied span iff there is content on both sides.
        let inside = match last {
            Some(last) => before > first && before <= last,
            None => false,
        };
        let fill = if inside { gap } else { pad };
        out.extend(std::iter::repeat_n(fill, count));
        next = before;
    }
    out.extend_from_slice(&row[next..]);
    out
}

/// Render the mapped query as a full-width row, padded to `new_width`.
///
/// The query row is placed at the offset the mapping chose; columns the query
/// never reached are filled with `pad`.
pub fn query_row_padded(mapping: &Mapping, query: &[u8], gap: u8, pad: u8) -> Vec<u8> {
    let core = mapping.query_row(query, gap);
    let lead = mapping.leading_skipped();
    let mut out = Vec::with_capacity(mapping.new_width);
    out.extend(std::iter::repeat_n(pad, lead));
    out.extend_from_slice(&core);
    while out.len() < mapping.new_width {
        out.push(pad);
    }
    out.truncate(mapping.new_width);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn prof(rows: &[&str]) -> Vec<[u32; ALPHA_LEN]> {
        let owned: Vec<Vec<u8>> = rows.iter().map(|r| r.as_bytes().to_vec()).collect();
        let refs: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
        profile_from_rows(&refs)
    }

    /// Render a mapping the way the surgery would, for eyeballing in tests.
    fn render(m: &Mapping, query: &[u8], rows: &[&str]) -> Vec<String> {
        let ins = m.insertions();
        let mut out = vec![String::from_utf8(query_row_padded(m, query, b'-', b'-')).unwrap()];
        for r in rows {
            out.push(String::from_utf8(splice_row(r.as_bytes(), &ins, b'-', b'-')).unwrap());
        }
        out
    }

    /// The motivating case: a gapless MSA, so the hand-built insertion has
    /// nowhere to go and must widen the alignment.
    #[test]
    fn insertion_widens_a_gapless_msa() {
        let rows = ["AACTGGGAT", "AACTGGGAT", "AACTGGGAT", "AACTGGGAT"];
        let hb = b"AACTTTTGGGAT";
        let m = align_to_profile(hb, &prof(&rows), &AlignParams::default()).unwrap();

        assert_eq!(m.old_width, 9);
        assert_eq!(m.new_width, 12, "three new columns for the extra TTT");
        // Left-aligned: the run of new columns sits before the existing T, not
        // after it.  Either is score-identical; the canonicaliser fixes the choice.
        assert_eq!(m.insertions(), vec![(3, 3)]);

        let r = render(&m, hb, &rows);
        assert_eq!(r[0], "AACTTTTGGGAT");
        for instance in &r[1..] {
            assert_eq!(instance, "AAC---TGGGAT");
        }
    }

    /// A length-neutral one-base indel *pair* loses to a plain substitution:
    /// opening a column plus the compensating column skip costs far more than a
    /// single mismatch.  The aligner reports substitutions and leaves the
    /// column structure alone.
    ///
    /// This is intended — it stops small hand edits from churning the alignment
    /// — but it means an equal-length hand-built consensus essentially never
    /// triggers surgery.  Surgery happens when the hand-built sequence changes
    /// *length*.
    #[test]
    fn single_substitution_beats_a_length_neutral_indel_pair() {
        let rows = ["AACTGGGAT", "AACTGGGAT", "AACTGGGAT", "AACTGGGAT"];
        let hb = b"ATCTTGGAT";
        let m = align_to_profile(hb, &prof(&rows), &AlignParams::default()).unwrap();

        assert!(m.insertions().is_empty());
        assert_eq!(m.new_width, 9);
        let r = render(&m, hb, &rows);
        assert_eq!(r[0], "ATCTTGGAT");
    }

    /// Realistic scale: a diverged seed alignment where the curator's edit is a
    /// 6 bp insertion plus a few substitutions.  The mapping must recover
    /// exactly that insertion and nothing else.
    #[test]
    fn realistic_seed_with_a_six_base_insertion() {
        let cons: &str = "GGCTAACTGCAGGATCCTTAGCAACGTTGCAATCCGGATTAGCCTTAAGGCTA\
                          CCGTTAACCGGATCAGCTTAGGCATCCGATTACGGCATTAGCCGGATTACGCA";
        // Four instances, each diverged from the consensus at scattered sites.
        let mut rows: Vec<String> = Vec::new();
        for (n, muts) in [(0usize, vec![3usize, 40, 77]),
                          (1, vec![10, 55, 90]),
                          (2, vec![7, 22, 61, 99]),
                          (3, vec![15, 33, 70])] {
            let mut s: Vec<u8> = cons.bytes().collect();
            for m in muts {
                s[m] = match s[m] { b'A' => b'G', b'G' => b'A', b'C' => b'T', _ => b'C' };
            }
            let _ = n;
            rows.push(String::from_utf8(s).unwrap());
        }
        let row_refs: Vec<&str> = rows.iter().map(|s| s.as_str()).collect();

        // Curator's hand-built version: same sequence with TTATCC inserted at 50
        // and two substitutions elsewhere.
        let mut hb: Vec<u8> = cons.bytes().collect();
        hb[20] = b'T';
        hb[80] = b'A';
        let tail = hb.split_off(50);
        hb.extend_from_slice(b"TTATCC");
        hb.extend_from_slice(&tail);

        let m = align_to_profile(&hb, &prof(&row_refs), &AlignParams::default()).unwrap();

        assert_eq!(m.insertions(), vec![(50, 6)], "exactly the curator's insertion");
        assert_eq!(m.new_width, cons.len() + 6);

        let r = render(&m, &hb, &row_refs);
        assert_eq!(r[0].len(), cons.len() + 6);
        assert_eq!(&r[0][50..56], "TTATCC");
        for instance in &r[1..] {
            assert_eq!(&instance[50..56], "------", "instances gap over the insertion");
        }

        // And the guard statistics should look unremarkable.
        let consensus = aln_core::consensus::build_consensus_from_sequences(
            &rows.iter().map(|s| s.as_bytes()).collect::<Vec<_>>(),
            &aln_core::consensus::ConsensusParams::default(),
        );
        let st = m.stats(&hb, &consensus);
        assert_eq!(st.inserted, 6);
        assert_eq!(st.skipped, 0);
        assert!(st.identity > 0.95, "identity was {}", st.identity);
        assert!(st.coverage > 0.99, "coverage was {}", st.coverage);
    }

    /// The case that matters on real seeds: the alignment already has an insert
    /// column, so the hand-built insertion is absorbed and the width is unchanged.
    #[test]
    fn insertion_absorbed_into_existing_gap_column() {
        // Column 4 is an insert column — one instance has a base, the rest gaps.
        let rows = ["AACT-GGGAT", "AACT-GGGAT", "AACTTGGGAT", "AACT-GGGAT"];
        let hb = b"AACTTGGGAT";
        let m = align_to_profile(hb, &prof(&rows), &AlignParams::default()).unwrap();

        assert!(m.insertions().is_empty(), "should reuse the existing column");
        assert_eq!(m.new_width, m.old_width);

        let r = render(&m, hb, &rows);
        assert_eq!(r[0], "AACTTGGGAT");
    }

    /// A hand-built deletion needs no surgery at all — the row just gets a gap.
    #[test]
    fn deletion_needs_no_new_columns() {
        let rows = ["AACTGGGAT", "AACTGGGAT", "AACTGGGAT"];
        let hb = b"AACTGGAT";
        let m = align_to_profile(hb, &prof(&rows), &AlignParams::default()).unwrap();

        assert!(m.insertions().is_empty());
        assert_eq!(m.new_width, 9);
        let r = render(&m, hb, &rows);
        // The gap is left-aligned to the first column of the GGG run.
        assert_eq!(r[0], "AACT-GGAT");
    }

    /// Identical sequence maps to identical columns and changes nothing.
    #[test]
    fn identity_is_a_no_op() {
        let rows = ["AACT-GGGAT", "AACTGGGGAT", "AACT-GGGAT"];
        let hb = b"AACTGGGAT";
        let m = align_to_profile(hb, &prof(&rows), &AlignParams::default()).unwrap();
        assert!(m.insertions().is_empty());
        let r = render(&m, hb, &rows);
        assert_eq!(r[0], "AACT-GGGAT");
    }

    /// Gap placement inside a homopolymer run is arbitrary on score; the
    /// canonicaliser must make it reproducible and left-shifted.
    #[test]
    fn homopolymer_gap_is_left_shifted() {
        let rows = ["AAGGGGTT", "AAGGGGTT", "AAGGGGTT"];
        let hb = b"AAGGGTT";
        let m = align_to_profile(hb, &prof(&rows), &AlignParams::default()).unwrap();
        let r = render(&m, hb, &rows);
        // The deleted G must be the leftmost of the run, not an interior one.
        assert_eq!(r[0], "AA-GGGTT");
    }

    /// The query may cover only part of the alignment's width.
    #[test]
    fn free_end_gaps_allow_a_short_query() {
        let rows = ["TTTTAACTGGGATTTTT", "TTTTAACTGGGATTTTT", "TTTTAACTGGGATTTTT"];
        let hb = b"AACTGGGAT";
        let m = align_to_profile(hb, &prof(&rows), &AlignParams::default()).unwrap();
        assert!(m.insertions().is_empty());
        assert_eq!(m.leading_skipped(), 4);
        let r = render(&m, hb, &rows);
        assert_eq!(r[0], "----AACTGGGAT----");
    }

    /// Splicing must not disturb a row's flanking-padding convention.
    #[test]
    fn splice_distinguishes_padding_from_internal_gaps() {
        let row = b"  ACGTACGT  ";
        // Insert inside the span, and beyond the right edge.
        let out = splice_row(row, &[(4, 1), (11, 2)], b'-', b' ');
        assert_eq!(String::from_utf8(out).unwrap(), "  AC-GTACGT    ");
    }

    /// A seed alignment is mostly short fragments.  Their flanking gaps must not
    /// count as evidence of a deletion, or every column looks gap-dominated and
    /// the aligner skips the whole alignment to open a new column per residue.
    ///
    /// This only shows up at scale — on a real 490-sequence LINE seed it turned
    /// a 12-column insertion into 3266 — so it is pinned here.
    #[test]
    fn flanking_padding_is_not_counted_as_deletion_evidence() {
        // Six fragments, each covering 10 of 20 columns, staggered.
        let rows = [
            "ACGTACGTAC----------",
            "ACGTACGTAC----------",
            "-----CGTACGTACG-----",
            "-----CGTACGTACG-----",
            "----------GTACGTACGT",
            "----------GTACGTACGT",
        ];
        let hb = b"ACGTACGTACGTACGTACGT";
        let m = align_to_profile(hb, &prof(&rows), &AlignParams::default()).unwrap();

        assert!(
            m.insertions().is_empty(),
            "fragments' flanks made the aligner widen the alignment: {:?}",
            m.insertions()
        );
        assert_eq!(m.new_width, 20);
        let r = render(&m, hb, &rows);
        assert_eq!(r[0], "ACGTACGTACGTACGTACGT");
    }

    /// A row consisting entirely of gaps contributes nothing and must not panic.
    #[test]
    fn all_gap_row_is_ignored() {
        let rows = ["ACGTACGT", "--------", "ACGTACGT"];
        let hb = b"ACGTACGT";
        let m = align_to_profile(hb, &prof(&rows), &AlignParams::default()).unwrap();
        assert!(m.insertions().is_empty());
        assert_eq!(render(&m, hb, &rows)[0], "ACGTACGT");
    }

    #[test]
    fn empty_query_is_an_error() {
        let rows = ["ACGT"];
        assert!(align_to_profile(b"", &prof(&rows), &AlignParams::default()).is_err());
    }
}
