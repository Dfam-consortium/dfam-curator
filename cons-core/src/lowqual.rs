//! Low-quality block detection and repair — a port of RepeatModeler `Refiner`'s
//! `resolveLowQualityBlocks`, together with the `MultAln.pm` machinery it rests
//! on (`getLowScoringAlignmentColumns` and
//! `_ruzzoTompaFindAllMaximalScoringSubsequences`).
//!
//! # What this fixes
//!
//! Refinement aligns every instance to one consensus and re-calls. Where the
//! instances genuinely disagree about an indel, that produces a *locally* bad
//! stretch of alignment: the aligner places the same event at different offsets
//! in different rows, and the column-wise caller then averages incompatible
//! registers into a consensus that matches nothing. Iterating does not fix it,
//! because each pass re-derives the same bad columns from the same bad
//! alignment. The repair identifies those stretches and re-derives the
//! consensus *for that block alone* from the instances, ignoring the column
//! structure that misled the caller.
//!
//! # The pipeline
//!
//! 1. Score every alignment column (instances against the reference row).
//! 2. Invert, and run Ruzzo-Tompa to find maximal-scoring runs of the inverted
//!    profile — i.e. the *worst*-scoring runs of the original.
//! 3. For each qualifying block, re-derive a consensus from the instances.
//! 4. Splice those into the gapped consensus, then refine again.
//!
//! # Relationship to `dfam_curator::quality`
//!
//! `dfam-curator` carries an earlier, independent take on step 1-2. It differs
//! in the scoring: it takes each column's best-fitting IUPAC score over the
//! whole column, where `Refiner` scores each *row against the reference* with
//! affine gap penalties and averages. Both find bad columns; only this one is
//! `Refiner`. They should be unified, but that module is not reachable from
//! here (`cons-core` is a dependency of `dfam-curator`, not the reverse).

use aln_core::consensus::{build_consensus_from_sequences, ConsensusParams};
use aln_core::msa::MultiAlign;
use aln_core::seq as seqmod;
use aln_core::SubstMatrix;

/// `Refiner`'s literal defaults from `getLowScoringAlignmentColumns`.
///
/// They are only meaningful alongside the matrix `Refiner` pairs them with: it
/// hardcodes `comparison.matrix` for this scoring regardless of what the
/// aligner used, and that matrix has a mean ACGT diagonal of 9.5. So the
/// penalties are really *ratios* — about -4.2 and -1.6 match-equivalents — and
/// using them verbatim with a differently-scaled matrix is a bug: against
/// xrepmask/3 (diagonal 3.0) a single gap outweighs thirteen matches, the
/// column profile goes negative nearly everywhere, and Ruzzo-Tompa returns a
/// handful of enormous blocks instead of the local defects it is looking for.
/// Measured on a Charlie1 family: 8,515 of 8,638 columns negative, 8 blocks, 7
/// of them wider than [`MAX_BLOCK`].
///
/// [`scaled_gap_penalties`] restores the ratio for any matrix, and reproduces
/// these numbers exactly on a comparison-scale one.
pub const GAP_OPEN: f64 = -40.0;
pub const GAP_EXTEND: f64 = -15.0;
/// The mean ACGT diagonal of `comparison.matrix`, which [`GAP_OPEN`] and
/// [`GAP_EXTEND`] are implicitly relative to.
const REFINER_MATCH_SCALE: f64 = 9.5;
/// A column belongs to a low-quality block when its Ruzzo-Tompa segment score
/// reaches this. `MultAln.pm` defaults to 1.
pub const THRESHOLD: f64 = 1.0;

/// Blocks narrower than this are not worth repairing (the Perl requires
/// `blockWidth > 1`).
pub const MIN_BLOCK: usize = 2;
/// `Refiner` skips wider blocks — its comment says they are "too big to align
/// with perl". The limit is kept because it also bounds the O(n^2) fallback.
pub const MAX_BLOCK: usize = 50;
/// A block needs at least this many instances to be worth re-deriving
/// (`$#{$instSeqs} >= 3`, i.e. four or more rows).
pub const MIN_INSTANCES: usize = 4;

/// A repaired block: replace columns `start..=end` of the gapped consensus with
/// `cons` (padded with gaps to the block width).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlockFix {
    pub start: usize,
    pub end: usize,
    /// Ungapped replacement sequence. May be shorter than the block.
    pub cons: Vec<u8>,
}

/// Gap penalties for `matrix`, holding `Refiner`'s penalty-to-match ratio.
///
/// Returns [`GAP_OPEN`] / [`GAP_EXTEND`] unchanged for a matrix on
/// `comparison.matrix`'s scale, and scales them down for a divided matrix such
/// as xrepmask/3 so that "a gap costs about four matches" stays true.
pub fn scaled_gap_penalties(matrix: &SubstMatrix) -> (f64, f64) {
    let diag: Vec<f64> = b"ACGT"
        .iter()
        .filter_map(|&b| matrix.score(b, b))
        .map(|v| v as f64)
        .collect();
    if diag.is_empty() {
        return (GAP_OPEN, GAP_EXTEND);
    }
    let mean = diag.iter().sum::<f64>() / diag.len() as f64;
    if mean <= 0.0 {
        return (GAP_OPEN, GAP_EXTEND);
    }
    let k = mean / REFINER_MATCH_SCALE;
    (GAP_OPEN * k, GAP_EXTEND * k)
}

/// Per-column alignment quality, averaged over the rows that reach each column.
///
/// Each instance row is walked against the reference row over the columns it
/// spans. A column where exactly one side is a gap is charged an affine gap
/// penalty; a column gapped on both sides is counted but scores nothing;
/// otherwise the substitution matrix decides. The sum is divided by the number
/// of rows contributing, so a column's score does not simply track coverage.
pub fn column_profile(
    msa: &MultiAlign,
    matrix: &SubstMatrix,
    gap_open: f64,
    gap_extend: f64,
) -> Vec<f64> {
    column_profile_weighted(msa, matrix, gap_open, gap_extend, false)
}

/// As [`column_profile`], optionally giving each **copy** a weight of one
/// rather than each row.
///
/// Under `HspPolicy::All` one copy can contribute several rows — 1.74 on
/// average over 180 simulated families, 30% of copies more than one, one copy
/// 13. Those rows are not independent observations of how well the consensus
/// explains the family; they are fragments of a single copy, and averaging over
/// rows lets a fragmented copy shout down the rest. It also dilutes the
/// profile: extra rows that agree with the consensus raise the mean and hide
/// genuinely bad columns, which is why the same families yield 16,177 flagged
/// blocks under `all` and 28,573 under `best`.
///
/// With `per_copy`, a row contributes `1 / rows_of_that_copy` to both the sum
/// and the count, so every copy carries the same weight wherever it is present.
pub fn column_profile_weighted(
    msa: &MultiAlign,
    matrix: &SubstMatrix,
    gap_open: f64,
    gap_extend: f64,
    per_copy: bool,
) -> Vec<f64> {
    let width = msa.width();
    let mut profile = vec![0.0f64; width];
    let mut counts = vec![0.0f64; width];
    let Some(reference) = msa.sequences.first() else {
        return profile;
    };

    // How many rows each copy contributed, so its weight can be split.
    let mut rows_per_copy: std::collections::HashMap<&str, usize> =
        std::collections::HashMap::new();
    if per_copy {
        for row in msa.sequences.iter().skip(1) {
            *rows_per_copy.entry(row.name.as_str()).or_insert(0) += 1;
        }
    }

    for row in msa.sequences.iter().skip(1) {
        let w = if per_copy {
            1.0 / *rows_per_copy.get(row.name.as_str()).unwrap_or(&1) as f64
        } else {
            1.0
        };
        // `in_gap` tracks affine state along this row only, so a run of gap
        // columns is charged one open and the rest extends — exactly the
        // Perl's per-sequence `$inGap` flag.
        let mut in_gap = false;
        for col in row.col_start..row.col_end.min(width) {
            let t = row.seq[col];
            let r = reference.seq[col];
            let t_gap = seqmod::is_gap(t) || t == b' ';
            let r_gap = seqmod::is_gap(r) || r == b' ';
            if t_gap != r_gap {
                profile[col] += w * if in_gap { gap_extend } else { gap_open };
                counts[col] += w;
                in_gap = true;
            } else if t_gap && r_gap {
                // Counted but unscored, as in the Perl.
                counts[col] += w;
            } else {
                profile[col] += w * matrix.score(r, t).unwrap_or(0) as f64;
                counts[col] += w;
                in_gap = false;
            }
        }
    }
    for (p, &c) in profile.iter_mut().zip(&counts) {
        if c > 0.0 {
            *p /= c;
        }
    }
    profile
}

/// Ruzzo-Tompa: all maximal-scoring contiguous subsequences, returned as a
/// per-position mask carrying each position's segment score (0 outside any
/// segment).
///
/// A direct port of `_ruzzoTompaFindAllMaximalScoringSubsequences`, including
/// its structure: intervals are kept in parallel arrays and merged leftwards
/// while a previous interval has a lower left prefix-sum and a lower right
/// prefix-sum. `k` counts *valid* intervals, and merging rolls it back — so
/// entries above `k` are stale and must not be read.
pub fn ruzzo_tompa_mask(b: &[f64]) -> Vec<f64> {
    #[derive(Clone, Copy)]
    struct Iv {
        lidx: usize,
        end: usize,
        l: f64,
        r: f64,
    }
    let mut iv: Vec<Iv> = Vec::new();
    let mut total = 0.0f64;
    let mut k = 0usize;

    for (i, &bi) in b.iter().enumerate() {
        total += bi;
        if bi <= 0.0 {
            continue;
        }
        let cur = Iv { lidx: i, end: i + 1, l: total - bi, r: total };
        if k < iv.len() {
            iv[k] = cur;
        } else {
            iv.push(cur);
        }
        loop {
            let maxj = (0..k).rev().find(|&j| iv[j].l < iv[k].l);
            match maxj {
                Some(j) if iv[j].r < iv[k].r => {
                    iv[j].end = i + 1;
                    iv[j].r = total;
                    k = j;
                }
                _ => {
                    k += 1;
                    break;
                }
            }
        }
    }

    let mut mask = vec![0.0f64; b.len()];
    for s in iv.iter().take(k) {
        let score = s.r - s.l;
        for m in mask.iter_mut().take(s.end).skip(s.lidx) {
            *m = score;
        }
    }
    mask
}

/// Column ranges whose Ruzzo-Tompa segment score reaches `threshold`.
///
/// Ranges are inclusive on both ends and in gapped (alignment) coordinates.
/// As [`low_scoring_columns`], with the per-copy profile weighting.
pub fn low_scoring_columns_weighted(
    msa: &MultiAlign,
    matrix: &SubstMatrix,
    threshold: f64,
    per_copy: bool,
) -> Vec<(usize, usize)> {
    let (go, ge) = scaled_gap_penalties(matrix);
    let profile = column_profile_weighted(msa, matrix, go, ge, per_copy);
    segments_from_profile(&profile, threshold)
}

/// Maximal low-scoring segments of an already-computed column profile.
fn segments_from_profile(profile: &[f64], threshold: f64) -> Vec<(usize, usize)> {
    // Invert: the worst-scoring runs of the profile are the best-scoring runs
    // of its negation, which is what Ruzzo-Tompa finds.
    let inverted: Vec<f64> = profile.iter().map(|p| -p).collect();
    let mask = ruzzo_tompa_mask(&inverted);

    let mut out = Vec::new();
    let mut start: Option<usize> = None;
    for (i, &v) in mask.iter().enumerate() {
        if v >= threshold {
            start.get_or_insert(i);
        } else if let Some(s) = start.take() {
            out.push((s, i - 1));
        }
    }
    if let Some(s) = start {
        out.push((s, mask.len() - 1));
    }
    out
}

pub fn low_scoring_columns(
    msa: &MultiAlign,
    matrix: &SubstMatrix,
    threshold: f64,
) -> Vec<(usize, usize)> {
    let (go, ge) = scaled_gap_penalties(matrix);
    let profile = column_profile(msa, matrix, go, ge);
    segments_from_profile(&profile, threshold)
}

/// The ungapped instance sequences spanning a block, and the reference's
/// ungapped length there.
///
/// Only rows that fully span the block contribute, matching the Perl's
/// `$start >= alignedStart && $end <= alignedEnd`. A row that spans the block
/// but is gapped throughout yields an empty string, which the Perl also keeps —
/// it is evidence that the block should be deleted.
fn block_sequences(msa: &MultiAlign, start: usize, end: usize) -> (usize, Vec<Vec<u8>>) {
    let mut inst = Vec::new();
    for row in msa.sequences.iter().skip(1) {
        if start >= row.col_start && end < row.col_end {
            let raw: Vec<u8> = row.seq[start..=end]
                .iter()
                .copied()
                .filter(|&b| !seqmod::is_gap(b) && b != b' ')
                .collect();
            inst.push(raw);
        }
    }
    let ref_len = msa
        .sequences
        .first()
        .map(|r| {
            r.seq[start..=end]
                .iter()
                .filter(|&&b| !seqmod::is_gap(b) && b != b' ')
                .count()
        })
        .unwrap_or(0);
    (ref_len, inst)
}


/// Debug: dump a candidate block's rows in alignment columns when
/// `TE_COMPOSER_BLOCK_LOG` is set. `source` names the selection, `outcome`
/// what became of the block, `repl` the replacement if any.
fn log_block(msa: &MultiAlign, source: &str, start: usize, end: usize, outcome: &str, repl: Option<&[u8]>) {
    if std::env::var_os("TE_COMPOSER_BLOCK_LOG").is_none() {
        return;
    }
    let slice = |row: &aln_core::msa::SequenceRow| -> String {
        String::from_utf8_lossy(&row.seq[start..=end.min(row.seq.len().saturating_sub(1))]).into_owned()
    };
    eprintln!(
        "BLOCK\t{source}\t{start}\t{end}\t{outcome}\t{}",
        repl.map(|r| String::from_utf8_lossy(r).into_owned()).unwrap_or_default()
    );
    if let Some(r) = msa.sequences.first() {
        eprintln!("BLOCKROW\tconsensus\t{}", slice(r));
    }
    for row in msa.sequences.iter().skip(1) {
        let spans = start >= row.col_start && end < row.col_end;
        eprintln!("BLOCKROW\t{}\t{}\t{}", row.name, slice(row), if spans { "spans" } else { "partial" });
    }
}

/// Find repairs for every low-quality block.
///
/// For each block, the instances' *ungapped* lengths are histogrammed. When the
/// most common length differs from the reference's, the block is re-derived:
///
/// * **Majority case** — that length is held by more than half the instances
///   and by at least three. Those instances are equal length by construction,
///   so a plain column-wise call over them gives the replacement directly.
/// * **Otherwise** — no length dominates, so the instances are aligned
///   all-against-all, the one scoring highest against the rest is taken as a
///   local reference, and the consensus of its alignments becomes the
///   replacement. This is the expensive branch, which is why blocks wider than
///   [`MAX_BLOCK`] are skipped.
///
/// The all-against-all branch needs a pairwise aligner; callers pass one so
/// `cons-core` does not choose a backend on their behalf. Returning `None` from
/// it (or passing `None`) skips that branch and keeps only majority repairs.
/// Why [`resolve_modal_length`] did or did not produce a replacement.
///
/// Exposed because the counts are the only way to tell a repair that found
/// nothing from one that was never given anything to work with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModalOutcome {
    /// A length held a strict majority and the block was re-derived from the
    /// instances of that length.
    Fixed(Vec<u8>),
    /// Fewer than [`MIN_INSTANCES`] instances span the block.
    TooFew,
    /// The consensus already has the length the instances vote for.
    AlreadyAgrees,
    /// No length holds a strict majority of at least three. This is the case
    /// `Refiner` hands to its all-against-all branch.
    NoMajority {
        /// The instances, for a caller that wants to run that branch.
        inst: Vec<Vec<u8>>,
    },
}

/// `Refiner`'s per-block resolution, on its own so it can be paired with a
/// block *selection* other than [`low_scoring_columns`].
///
/// The instances' ungapped lengths are histogrammed; the most common wins, ties
/// going to the longest — the Perl's 9/19/23 determinism fix. A strict majority
/// of at least three is required, and only the instances of that length vote on
/// the bases, which is what stops a mixture of lengths from being averaged into
/// a consensus matching none of them.
pub fn resolve_modal_length(
    msa: &MultiAlign,
    start: usize,
    end: usize,
    params: &ConsensusParams,
) -> ModalOutcome {
    let (ref_len, inst) = block_sequences(msa, start, end);
    if inst.len() < MIN_INSTANCES {
        return ModalOutcome::TooFew;
    }
    let mut histo: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for s in &inst {
        *histo.entry(s.len()).or_insert(0) += 1;
    }
    let mut by_count: Vec<(usize, usize)> = histo.into_iter().collect();
    by_count.sort_by(|a, b| b.1.cmp(&a.1).then(b.0.cmp(&a.0)));
    let (most_len, most_count) = by_count[0];

    if most_len == ref_len {
        return ModalOutcome::AlreadyAgrees;
    }
    if most_count * 2 > inst.len() && most_count >= 3 {
        let picked: Vec<&[u8]> = inst
            .iter()
            .filter(|s| s.len() == most_len)
            .map(|s| s.as_slice())
            .collect();
        return ModalOutcome::Fixed(seqmod::ungap(&build_consensus_from_sequences(
            &picked, params,
        )));
    }
    ModalOutcome::NoMajority { inst }
}

/// `Refiner`'s fallback when no instance length dominates a block.
///
/// `resolveLowQualityBlocks` reaches this branch whenever `mostFreqCount` fails
/// the majority test, and it is the half of the repair `cons-core` has never
/// had — `refine_with_repair` passed `None` and simply skipped those blocks.
/// The Perl runs `NeedlemanWunschGotohAlgorithm::search` all-against-all over
/// the block's instances with `linupmatrix` and -25/-5 penalties, sums each
/// sequence's scores, takes the highest as a reference, builds a `MultAln` from
/// its alignments and calls the consensus.
///
/// The same procedure here, with parasail's global alignment in place of the
/// Perl's Needleman-Wunsch: exact for the sizes involved (blocks are capped at
/// [`MAX_BLOCK`] columns) and far faster than an interpreted DP. Penalties come
/// from the caller so they stay on the matrix's scale rather than being pinned
/// to numbers that only mean anything against `comparison.matrix`.
pub fn resolve_by_all_vs_all(
    inst: &[Vec<u8>],
    matrix: &SubstMatrix,
    gap_open: u32,
    gap_extend: u32,
    params: &ConsensusParams,
) -> Option<Vec<u8>> {
    resolve_by_all_vs_all_detail(inst, matrix, gap_open, gap_extend, params).map(|d| d.0)
}

/// As [`resolve_by_all_vs_all`], also reporting which instance won and by how
/// much.
///
/// The winner is the whole basis of the re-derivation — every other instance is
/// projected onto its columns — so a diagnostic that cannot name it leaves the
/// reader unable to judge whether the block was resolved against a sensible
/// representative or an outlier that happened to score well.
pub fn resolve_by_all_vs_all_detail(
    inst: &[Vec<u8>],
    matrix: &SubstMatrix,
    gap_open: u32,
    gap_extend: u32,
    params: &ConsensusParams,
) -> Option<(Vec<u8>, usize, i64, i64)> {
    resolve_by_all_vs_all_opt(inst, matrix, gap_open, gap_extend, params, false)
}

/// As [`resolve_by_all_vs_all_detail`], optionally keeping bases the centre
/// lacks.
///
/// Projecting onto the centre's columns throws away every base a member carries
/// where the centre has a gap — measured at 6.1% of all member sequence, with
/// 15% of projections losing something. The re-derived block can then never be
/// longer than whichever instance happened to win, however many members agree
/// on an insertion it lacks. That is a hard ceiling on a routine whose purpose
/// is recovering insertions, and it is stricter than `Refiner`, whose `MultAln`
/// keeps such bases as columns and lets the consensus caller judge them.
///
/// With `keep_insertions`, each centre position gets as many extra columns as
/// the widest insertion any member places there, members are laid into them
/// left-justified, and the caller decides by occupancy as it does everywhere
/// else. An insertion a majority carries now survives; a private one still
/// resolves to a gap and disappears on ungapping.
pub fn resolve_by_all_vs_all_opt(
    inst: &[Vec<u8>],
    matrix: &SubstMatrix,
    gap_open: u32,
    gap_extend: u32,
    params: &ConsensusParams,
    keep_insertions: bool,
) -> Option<(Vec<u8>, usize, i64, i64)> {
    let (score_aligner, tb_aligner) = pack_aligners(matrix, gap_open, gap_extend)?;
    resolve_by_all_vs_all_with(inst, &score_aligner, &tb_aligner, params, keep_insertions)
}

/// The two aligners a centre-star needs: a score-only one for the
/// all-against-all that picks the centre, and one with traceback for laying
/// each member out against it. Built once per caller so a pass over thousands
/// of spans does not rebuild a parasail profile for every span.
pub fn pack_aligners(
    matrix: &SubstMatrix,
    gap_open: u32,
    gap_extend: u32,
) -> Option<(crate::FastAligner, crate::FastAligner)> {
    use aln_engine::{AlignMode, AlignParams};
    let mk = |traceback: bool| {
        let p = AlignParams {
            mode: AlignMode::Global,
            gap_open,
            gap_extend,
            // A short block can align end-to-end at a negative score and still
            // be the right answer; the floor exists to keep junk out of a
            // genome-scale search, which is not what is happening here.
            min_score: i32::MIN / 4,
            traceback,
            bandwidth: None,
        };
        crate::FastAligner::new(matrix.clone(), p).ok()
    };
    Some((mk(false)?, mk(true)?))
}

/// As [`resolve_by_all_vs_all_opt`], with the aligners supplied.
pub fn resolve_by_all_vs_all_with(
    inst: &[Vec<u8>],
    score_aligner: &crate::FastAligner,
    tb_aligner: &crate::FastAligner,
    params: &ConsensusParams,
    keep_insertions: bool,
) -> Option<(Vec<u8>, usize, i64, i64)> {
    use aln_engine::PairwiseAligner;

    if inst.len() < MIN_INSTANCES {
        return None;
    }
    let seqs: Vec<aln_core::Sequence> = inst
        .iter()
        .enumerate()
        .map(|(i, s)| aln_core::Sequence::new(format!("b{i}"), s.clone()))
        .collect();

    // All-against-all, summing each sequence's score against every other.
    let mut total = vec![0i64; seqs.len()];
    for i in 0..seqs.len() {
        for j in 0..seqs.len() {
            if i == j {
                continue;
            }
            // Score only: the centre is chosen on summed score, and a traceback
            // for every pair was the single largest cost of packing.
            if let Ok(Some(sc)) = score_aligner.score(&seqs[i], &seqs[j]) {
                total[i] += sc as i64;
            }
        }
    }
    let best = (0..seqs.len()).max_by_key(|&i| total[i])?;
    let runner_up = (0..seqs.len())
        .filter(|&i| i != best)
        .map(|i| total[i])
        .max()
        .unwrap_or(total[best]);

    let clen = seqs[best].seq.len();
    // Per member: the base sitting at each centre position, and whatever it
    // inserts *before* that position. Index `clen` holds a trailing insertion.
    let mut bases: Vec<Vec<u8>> = Vec::new();
    let mut inserts: Vec<Vec<Vec<u8>>> = Vec::new();
    for (i, q) in seqs.iter().enumerate() {
        if i == best {
            bases.push(seqs[best].seq.clone());
            inserts.push(vec![Vec::new(); clen + 1]);
            continue;
        }
        let Ok(Some(al)) = tb_aligner.align(q, &seqs[best]) else { continue };
        let Ok((gq, gs)) = al.gapped(&q.seq, &seqs[best].seq) else { continue };
        let mut b = vec![b'-'; clen];
        let mut ins = vec![Vec::new(); clen + 1];
        let mut si = 0usize;
        for (&qc, &sc) in gq.iter().zip(gs.iter()) {
            if sc == b'-' {
                if si <= clen {
                    ins[si].push(qc);
                }
            } else {
                if si < clen {
                    b[si] = qc;
                }
                si += 1;
            }
        }
        if si != clen {
            continue;
        }
        bases.push(b);
        inserts.push(ins);
        // Bases this member carries that the centre does not, and which the
        // projection therefore throws away. The ceiling of a centre-star: the
        // re-derived block can never exceed the centre's length, however many
        // members agree on an insertion the centre happens to lack.
        if std::env::var_os("TE_COMPOSER_PROJ_LOG").is_some() {
            // Every projection, lossy or not, so the discard rate is not
            // conditioned on discarding.
            let dropped = gq.iter().zip(gs.iter()).filter(|(_, &sc)| sc == b'-').count();
            eprintln!("PROJ\t{}\t{}\t{}", seqs[best].seq.len(), q.seq.len(), dropped);
        }
    }
    if bases.len() < MIN_INSTANCES {
        return None;
    }
    // Lay the members out. Without `keep_insertions` this reproduces the old
    // projection exactly: centre positions only, insertions dropped.
    let widths: Vec<usize> = if keep_insertions {
        (0..=clen)
            .map(|p| inserts.iter().map(|m| m[p].len()).max().unwrap_or(0))
            .collect()
    } else {
        vec![0; clen + 1]
    };
    // Which inserted columns are worth keeping. An insertion carried by a
    // handful of members is private and must not become consensus — but note
    // the members here are already only those with sequence in the span, so a
    // majority among *them* is not a majority of the family. Keeping every
    // column and trusting the caller produced a consensus 48 bases longer than
    // the ancestor on the first family tried.
    let keep_col: Vec<Vec<bool>> = (0..=clen)
        .map(|p| {
            (0..widths[p])
                .map(|k| {
                    let have = inserts.iter().filter(|m| m[p].len() > k).count();
                    have * 2 > bases.len()
                })
                .collect()
        })
        .collect();
    let rows: Vec<Vec<u8>> = bases
        .iter()
        .zip(inserts.iter())
        .map(|(b, ins)| {
            let mut row = Vec::with_capacity(clen + widths.iter().sum::<usize>());
            let emit = |row: &mut Vec<u8>, p: usize, seg: &[u8]| {
                for k in 0..widths[p] {
                    if !keep_col[p][k] {
                        continue;
                    }
                    row.push(seg.get(k).copied().unwrap_or(b'-'));
                }
            };
            for p in 0..clen {
                emit(&mut row, p, &ins[p]);
                row.push(b[p]);
            }
            emit(&mut row, clen, &ins[clen]);
            row
        })
        .collect();
    let refs: Vec<&[u8]> = rows.iter().map(|r| r.as_slice()).collect();
    let cons = seqmod::ungap(&build_consensus_from_sequences(&refs, params));
    if cons.is_empty() {
        None
    } else {
        Some((cons, best, total[best], runner_up))
    }
}

/// Re-derive the consensus across spans it currently has no bases for, so
/// inherited insertions can contribute.
///
/// A consensus-induced MSA aligns every copy to the reference independently, so
/// two copies' inserted bases are never aligned to *each other* and any shared
/// history among them is invisible. Measured on bootstrap-induced alignments,
/// that history is real: same-slot insertion segments share 0.72 identity at
/// 5-8 bp and 0.86 at 16+ bp, against 0.33-0.40 for unrelated segments of the
/// same length. In the simulations, where inserted bases are drawn at random,
/// such similarity can only be inherited.
///
/// Spans are merged across short stretches of called consensus before being
/// re-derived. 77% of adjacent gap runs are separated by two columns or fewer
/// (median one), so treating them separately splits single insertion events;
/// merging at five raises the inserted sequence this can reach from 47% to 95%
/// of the total. That is the problem `AutoRunBlocker`'s clustering
/// (`allowedGapDist = window/5`) and `resolveIndels`'s `-min_gap` were built
/// for, and their merge distances land in the same range.
///
/// Returns a new gapped consensus. A span is left alone unless re-deriving it
/// yields *more* sequence than is there now — this is meant to recover bases the
/// column-wise caller cannot see, not to relitigate bases it already called.
pub fn pack_insertion_spans(
    msa: &MultiAlign,
    gapped: &[u8],
    matrix: &SubstMatrix,
    call: &ConsensusParams,
    max_sep: usize,
    min_seg: usize,
    min_score: i64,
    min_occupancy: f64,
    keep_insertions: bool,
) -> Vec<u8> {
    let is_gap = |x: u8| seqmod::is_gap(x) || x == b' ';
    let width = gapped.len().min(msa.width());
    if width == 0 {
        return gapped.to_vec();
    }
    let (go, ge) = scaled_gap_penalties(matrix);
    let (go, ge) = (go.abs().round() as u32, ge.abs().round() as u32);
    let Some((score_aligner, tb_aligner)) = pack_aligners(matrix, go, ge) else {
        return gapped.to_vec();
    };

    // Maximal gap runs, then merged while the separation is short enough.
    let mut runs: Vec<(usize, usize)> = Vec::new();
    let mut col = 0usize;
    while col < width {
        if !is_gap(gapped[col]) {
            col += 1;
            continue;
        }
        let start = col;
        while col < width && is_gap(gapped[col]) {
            col += 1;
        }
        runs.push((start, col - 1));
    }
    // Merge while the span as a whole holds at most `max_sep` called columns.
    // Merging on the separation alone chained runs across the whole alignment
    // at high divergence, where the consensus has a gap column every few
    // positions: spans of 30-190 columns in which every copy had 40+ bases,
    // 6.7 million pairwise alignments per pass on one 100-copy family.
    let mut spans: Vec<(usize, usize)> = Vec::new();
    let mut called_in_span = 0usize;
    for r in runs {
        match spans.last_mut() {
            Some(prev) if called_in_span + r.0.saturating_sub(prev.1 + 1) <= max_sep => {
                called_in_span += r.0.saturating_sub(prev.1 + 1);
                prev.1 = r.1;
            }
            _ => {
                spans.push(r);
                called_in_span = 0;
            }
        }
    }

    let mut out = gapped.to_vec();
    let stats = std::env::var_os("TE_COMPOSER_PACK_LOG").is_some();
    let (mut n_spans, mut n_aligned, mut n_pairs, mut n_bases, mut n_packed) = (0usize, 0usize, 0usize, 0usize, 0usize);
    let n_spans_total = spans.len();
    for (a, b) in spans {
        n_spans += 1;
        // Copies spanning the whole span, whether or not they carry bases in it.
        let spanning = msa
            .sequences
            .iter()
            .skip(1)
            .filter(|r| r.col_start <= a && r.col_end > b)
            .count();
        // Which copies carry *inserted* bases here: bases in columns where the
        // consensus has none. Bases in the span's called columns are what every
        // copy has and say nothing about the insertion.
        let carriers = msa
            .sequences
            .iter()
            .skip(1)
            .filter(|r| r.col_start <= a && r.col_end > b)
            .map(|r| {
                (a..=b)
                    .filter(|&c| is_gap(gapped[c]) && !is_gap(r.seq[c]))
                    .count()
            })
            .filter(|&n| n > 0)
            .collect::<Vec<usize>>();
        if carriers.len() < 2 || carriers.iter().copied().max().unwrap_or(0) < min_seg {
            continue;
        }
        // Occupancy gate, measured against the copies that span the region:
        // an insertion the consensus should carry is one most spanning copies
        // carry. Without it any four copies with a shared private insertion
        // qualified, which both promoted minority insertions into the consensus
        // and made packing quadratic in copies over thousands of private spans.
        if (carriers.len() as f64) < min_occupancy * spanning as f64 {
            continue;
        }
        // What each instance carries across the span, ungapped.
        let inst: Vec<Vec<u8>> = msa
            .sequences
            .iter()
            .skip(1)
            .filter(|r| r.col_start <= a && r.col_end > b)
            .map(|r| {
                r.seq[a..=b]
                    .iter()
                    .copied()
                    .filter(|&x| !is_gap(x))
                    .collect::<Vec<u8>>()
            })
            .filter(|s: &Vec<u8>| !s.is_empty())
            .collect();
        if inst.len() < 2 {
            continue;
        }
        if inst.len() >= MIN_INSTANCES {
            n_aligned += 1;
            n_pairs += inst.len() * (inst.len() - 1) + inst.len();
            n_bases += inst.iter().map(|s| s.len()).sum::<usize>();
        }
        let Some((cons, winner, wscore, rscore)) =
            resolve_by_all_vs_all_with(&inst, &score_aligner, &tb_aligner, call, keep_insertions)
        else {
            continue;
        };
        // Only act when there is more sequence to be had, and only when it fits
        // the columns available — widening the alignment here would shift every
        // downstream column.
        // The winner's summed score against the rest is a direct test of
        // whether the block has any homology to re-derive. Inspecting the worst
        // hs1 regressions turned up winners scoring -214 against a runner-up of
        // -228: unrelated sequence forced into a common frame, where the
        // "winner" is merely the least-bad member and the column-wise call comes
        // out `N`. A non-positive score means there is nothing here to align.
        if wscore <= min_score {
            continue;
        }

        // Count only *called* bases. The old test compared raw lengths, so a
        // replacement that added nothing but `N` passed it — 23% of packed
        // spans on the worst hs1 families contained an N, and an N is a
        // mismatch against any truth while contributing no information. If the
        // instances cannot agree on a base, the honest outcome is to leave the
        // span alone.
        let called = |s: &[u8]| s.iter().filter(|c| !matches!(c, b'N' | b'n')).count();
        let here_called = gapped[a..=b]
            .iter()
            .filter(|&&x| !is_gap(x) && !matches!(x, b'N' | b'n'))
            .count();
        if called(&cons) <= here_called || cons.len() > b - a + 1 {
            continue;
        }
        if std::env::var_os("TE_COMPOSER_PACK_LOG").is_some() {
            let before: Vec<u8> =
                gapped[a..=b].iter().copied().filter(|&x| !is_gap(x)).collect();
            eprintln!(
                "PACKED\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                a,
                b,
                inst.len(),
                String::from_utf8_lossy(&before),
                String::from_utf8_lossy(&cons),
                inst.iter().map(|s| s.len().to_string()).collect::<Vec<_>>().join(","),
                winner,
                wscore,
                rscore
            );
            eprintln!("PACKROW\tconsensus\t{}", String::from_utf8_lossy(&gapped[a..=b]));
            for row in msa.sequences.iter().skip(1).filter(|r| r.col_start <= a && r.col_end > b) {
                eprintln!("PACKROW\t{}\t{}", row.name, String::from_utf8_lossy(&row.seq[a..=b]));
            }
        }
        n_packed += 1;
        for (k, col) in (a..=b).enumerate() {
            out[col] = cons.get(k).copied().unwrap_or(b'-');
        }
    }
    if stats {
        eprintln!(
            "PACKSTATS\tspans={n_spans_total}\tvisited={n_spans}\taligned={n_aligned}\tpairwise_alignments={n_pairs}\tsegment_bases={n_bases}\tpacked={n_packed}"
        );
    }
    out
}

/// Runs of columns where the *consensus* carries no base, kept when enough of
/// the instances spanning the run carry bases inside it.
///
/// This selection is neither `Refiner`'s nor `AutoRunBlocker`'s — both work
/// from column scores or fixed windows. It targets the one case those miss by
/// construction: a stretch the consensus has no positions for at all. Note that
/// an instance counts as occupying the run if it has a base *anywhere* in it,
/// so copies whose inserted bases the aligner scattered across different
/// columns still count as agreeing.
///
/// These blocks cannot be repaired by [`resolve_modal_length`]: the consensus
/// has zero bases there, and where fewer than half the spanning instances carry
/// bases the modal length is also zero, so it reports `AlreadyAgrees` and
/// declines. They need a resolution that re-derives the block from the
/// instances.
pub fn gap_run_blocks(
    msa: &MultiAlign,
    min_occupancy: f64,
    min_spanning: usize,
) -> Vec<(usize, usize)> {
    let Some(reference) = msa.sequences.first() else {
        return Vec::new();
    };
    let width = msa.width();
    let is_gap = |b: u8| seqmod::is_gap(b) || b == b' ';

    let mut runs = Vec::new();
    let mut start: Option<usize> = None;
    for col in 0..width {
        if is_gap(reference.seq[col]) {
            start.get_or_insert(col);
        } else if let Some(s) = start.take() {
            runs.push((s, col - 1));
        }
    }
    if let Some(s) = start {
        runs.push((s, width - 1));
    }

    runs.into_iter()
        .filter(|&(a, b)| {
            let spanning: Vec<_> = msa
                .sequences
                .iter()
                .skip(1)
                .filter(|r| r.col_start <= a && r.col_end > b)
                .collect();
            if spanning.len() < min_spanning {
                return false;
            }
            let occupied = spanning
                .iter()
                .filter(|r| r.seq[a..=b].iter().any(|&x| !is_gap(x)))
                .count();
            occupied as f64 / spanning.len() as f64 >= min_occupancy
        })
        .collect()
}

/// Apply [`resolve_modal_length`] to a caller-supplied block list.
///
/// `resolve_low_quality_blocks` selects *and* resolves; this resolves only, so a
/// second block source — [`length_vote_blocks`] — can feed the same vote. Blocks
/// outside `MIN_BLOCK..=max_block` are skipped, and blocks overlapping one
/// already in `existing` are dropped, so the returned fixes can be concatenated
/// and spliced in one pass.
pub fn resolve_given_blocks(
    msa: &MultiAlign,
    blocks: &[(usize, usize)],
    max_block: usize,
    params: &ConsensusParams,
    existing: &[BlockFix],
) -> Vec<BlockFix> {
    let mut out: Vec<BlockFix> = Vec::new();
    for &(start, end) in blocks {
        let width = end.saturating_sub(start) + 1;
        if end < start || width < MIN_BLOCK || width > max_block {
            continue;
        }
        let clash = |f: &BlockFix| start <= f.end && f.start <= end;
        if existing.iter().any(clash) || out.iter().any(clash) {
            continue;
        }
        match resolve_modal_length(msa, start, end, params) {
            ModalOutcome::Fixed(cons) => {
                log_block(msa, "window", start, end, "fixed", Some(&cons));
                out.push(BlockFix { start, end, cons });
            }
            other => {
                let tag = match other {
                    ModalOutcome::TooFew => "too-few",
                    ModalOutcome::AlreadyAgrees => "already-agrees",
                    _ => "no-majority",
                };
                log_block(msa, "window", start, end, tag, None);
            }
        }
    }
    out
}

pub fn resolve_low_quality_blocks<F>(
    msa: &MultiAlign,
    matrix: &SubstMatrix,
    threshold: f64,
    params: &ConsensusParams,
    mut align_block: Option<F>,
) -> Vec<BlockFix>
where
    F: FnMut(&[Vec<u8>]) -> Option<Vec<u8>>,
{
    let debug = std::env::var_os("TE_COMPOSER_REPAIR_DEBUG").is_some();
    let (mut n_narrow, mut n_wide, mut n_few, mut n_same, mut n_nomajority) = (0, 0, 0, 0, 0);
    let mut wide_widths: Vec<usize> = Vec::new();
    let mut fixes = Vec::new();
    for (start, end) in low_scoring_columns(msa, matrix, threshold) {
        let width = end - start + 1;
        if width < MIN_BLOCK {
            n_narrow += 1;
            continue;
        }
        if width > MAX_BLOCK {
            n_wide += 1;
            wide_widths.push(width);
            continue;
        }

        let cons = match resolve_modal_length(msa, start, end, params) {
            ModalOutcome::Fixed(c) => {
                log_block(msa, "ruzzo-tompa", start, end, "fixed-modal", Some(&c));
                c
            }
            ModalOutcome::TooFew => {
                n_few += 1;
                log_block(msa, "ruzzo-tompa", start, end, "too-few", None);
                continue;
            }
            ModalOutcome::AlreadyAgrees => {
                n_same += 1;
                log_block(msa, "ruzzo-tompa", start, end, "already-agrees", None);
                continue;
            }
            ModalOutcome::NoMajority { inst } => match align_block.as_mut() {
                Some(f) => match f(&inst) {
                    Some(c) => {
                        log_block(msa, "ruzzo-tompa", start, end, "fixed-centre-star", Some(&c));
                        c
                    }
                    None => continue,
                },
                None => {
                    n_nomajority += 1;
                    continue;
                }
            },
        };

        fixes.push(BlockFix { start, end, cons });
    }
    if debug {
        wide_widths.sort_unstable();
        eprintln!(
            "repair-blocks: rejected narrow={n_narrow} wide={n_wide} too-few-instances={n_few} \
             already-agrees={n_same} no-majority={n_nomajority} -> {} fixes; wide widths {:?}",
            fixes.len(),
            wide_widths
        );
    }
    fixes
}

/// Splice repairs into a gapped consensus.
///
/// Each fix overwrites its columns with the replacement, right-padded with gaps
/// to the block width — so the alignment's column count never changes, and a
/// shorter replacement simply leaves gap columns behind. Fixes are applied
/// right to left so earlier offsets stay valid.
pub fn patch_gapped_consensus(gapped: &[u8], fixes: &[BlockFix]) -> Vec<u8> {
    let mut out = gapped.to_vec();
    let mut sorted: Vec<&BlockFix> = fixes.iter().collect();
    sorted.sort_by_key(|f| std::cmp::Reverse(f.start));
    for f in sorted {
        if f.end >= out.len() {
            continue;
        }
        let width = f.end - f.start + 1;
        if f.cons.len() > width {
            // A replacement cannot widen the alignment; the Perl's substr
            // assignment would silently shift every downstream column.
            continue;
        }
        let mut repl = f.cons.clone();
        repl.resize(width, b'-');
        out[f.start..=f.end].copy_from_slice(&repl);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use aln_core::msa::SequenceRow;

    /// The mask must carry each maximal segment's own score, and zero
    /// elsewhere — the property `low_scoring_columns` thresholds on.
    #[test]
    fn ruzzo_tompa_marks_maximal_segments() {
        // One clear positive island surrounded by negatives.
        let b = [-1.0, 5.0, 3.0, -1.0, -1.0, 2.0, -9.0];
        let m = ruzzo_tompa_mask(&b);
        assert_eq!(m[0], 0.0, "negative prefix is outside any segment");
        assert_eq!(m[1], 8.0, "5+3 island scores 8");
        assert_eq!(m[2], 8.0);
        assert_eq!(m[5], 2.0, "the lone 2 is its own segment");
        assert_eq!(m[6], 0.0);
    }

    /// The penalties must reduce to Refiner's own numbers on Refiner's matrix,
    /// and scale down proportionally on a divided one.
    #[test]
    fn gap_penalties_track_the_matrix_scale() {
        let cmp = SubstMatrix::parse(
            "FREQS A 0.265 C 0.235 G 0.235 T 0.265\n\
               A   R   G   C   Y   T   N\n\
               9   1  -6 -15 -16 -17  -1\n\
               1   1   1 -15 -15 -16  -1\n\
              -6   1  10 -15 -15 -15  -1\n\
             -15 -15 -15  10   2  -6  -1\n\
             -16 -15 -15   1   1   1  -1\n\
             -17 -16 -15  -6   1   9  -1\n\
              -1  -1  -1  -1  -1  -1  -1\n").unwrap();
        let (go, ge) = scaled_gap_penalties(&cmp);
        assert!((go - GAP_OPEN).abs() < 1e-9, "got {go}");
        assert!((ge - GAP_EXTEND).abs() < 1e-9, "got {ge}");

        // A matrix at one third the scale must charge one third the penalty,
        // so "a gap costs about four matches" survives the rescaling.
        let third = SubstMatrix::parse(
            "FREQS A 0.265 C 0.235 G 0.235 T 0.265\n\
               A   R   G   C   Y   T   N\n\
               3   0  -2  -5  -5  -6  -1\n\
               0   0   0  -5  -5  -5  -1\n\
              -2   0   3  -5  -5  -5  -1\n\
              -5  -5  -5   3   1  -2  -1\n\
              -5  -5  -5   0   0   0  -1\n\
              -6  -5  -5  -2   0   3  -1\n\
              -1  -1  -1  -1  -1  -1  -1\n").unwrap();
        let (go3, _) = scaled_gap_penalties(&third);
        assert!(go3 > GAP_OPEN, "a smaller matrix must charge a smaller penalty");
        assert!((go3 / go - 3.0 / 9.5).abs() < 0.05, "ratio off: {}", go3 / go);
    }

    #[test]
    fn an_all_negative_profile_has_no_segments() {
        assert_eq!(ruzzo_tompa_mask(&[-1.0, -2.0, -3.0]), vec![0.0, 0.0, 0.0]);
    }

    /// A replacement shorter than its block is gap-padded, and the alignment
    /// width is preserved.
    #[test]
    fn patching_preserves_width() {
        let gapped = b"AAAACCCCGGGG".to_vec();
        let fixes = vec![BlockFix { start: 4, end: 7, cons: b"TT".to_vec() }];
        let out = patch_gapped_consensus(&gapped, &fixes);
        assert_eq!(out, b"AAAATT--GGGG".to_vec());
        assert_eq!(out.len(), gapped.len(), "width must not change");
    }

    /// Overlapping-free multi-block patching applies every fix.
    #[test]
    fn multiple_blocks_all_apply() {
        let gapped = b"AAAACCCCGGGGTTTT".to_vec();
        let fixes = vec![
            BlockFix { start: 0, end: 3, cons: b"GG".to_vec() },
            BlockFix { start: 12, end: 15, cons: b"CC".to_vec() },
        ];
        assert_eq!(patch_gapped_consensus(&gapped, &fixes), b"GG--CCCCGGGGCC--".to_vec());
    }

    /// A window where the copies unanimously carry more bases than the
    /// consensus must be flagged, and the delta must say how many.
    #[test]
    fn a_shared_insertion_is_detected() {
        // Consensus lacks 3 ancestral bases that every copy has. The copies'
        // insertions are deliberately placed in *different* columns, which is
        // the misregistration this detector has to see through.
        let mk = |name: &str, seq: &[u8]| SequenceRow::new(name, seq.to_vec());
        // cons:  ACGTA---CGTAC  (10 bases)
        // copies each have GGG, but at three different offsets
        let cons = mk("cons", b"ACGTA---CGTAC");
        let rows = vec![
            mk("a", b"ACGTAGGG-CGTAC".split_at(13).0.to_vec().as_slice()),
            mk("b", b"ACGTGGG--CGTAC".split_at(13).0.to_vec().as_slice()),
            mk("c", b"ACGGGG-TACGTAC".split_at(13).0.to_vec().as_slice()),
            mk("d", b"ACGTAGGG-CGTAC".split_at(13).0.to_vec().as_slice()),
            mk("e", b"ACGTAGGG-CGTAC".split_at(13).0.to_vec().as_slice()),
        ];
        let msa = MultiAlign::from_sequences(cons, rows).unwrap();
        let p = VoteParams { window: 10, min_copies: 3, min_ratio: 1.0, merge_gap: 2 };
        let blocks = length_vote_blocks(&msa, &p);
        assert!(!blocks.is_empty(), "a unanimous 3bp insertion must be flagged");
        assert!(blocks[0].length_delta() > 0,
                "delta must say the consensus is too short, got {}", blocks[0].length_delta());
    }

    /// Substitution noise alone must not trip the detector — that is the whole
    /// point of voting on length rather than on a quality score.
    #[test]
    fn substitution_noise_alone_is_ignored() {
        let mk = |name: &str, seq: &[u8]| SequenceRow::new(name, seq.to_vec());
        let msa = MultiAlign::from_sequences(
            mk("cons", b"ACGTACGTAC"),
            vec![mk("a", b"TTTTACGTAC"), mk("b", b"ACGTTTTTAC"),
                 mk("c", b"ACGTACGTTT"), mk("d", b"TTGTACGTAC")],
        ).unwrap();
        let p = VoteParams { window: 10, min_copies: 3, min_ratio: 1.0, merge_gap: 2 };
        assert!(length_vote_blocks(&msa, &p).is_empty(),
                "same length everywhere: nothing to vote on");
    }

    /// A replacement longer than its block is refused rather than shifting
    /// every downstream column.
    #[test]
    fn an_oversized_replacement_is_refused() {
        let gapped = b"AAAACCCC".to_vec();
        let fixes = vec![BlockFix { start: 0, end: 1, cons: b"GGGG".to_vec() }];
        assert_eq!(patch_gapped_consensus(&gapped, &fixes), gapped);
    }
}

// ─── Length-vote detection (AutoRunBlocker) ──────────────────────────────────

/// Tuning for [`length_vote_blocks`], from `AutoRunBlocker.pl`'s CLI.
#[derive(Debug, Clone, Copy)]
pub struct VoteParams {
    /// Window width, in **consensus positions** (`-windowSize`).
    pub window: usize,
    /// Minimum copies that must agree on the majority length (`-minCopyAgreement`).
    pub min_copies: usize,
    /// Minimum ratio of majority-length copies to copies keeping the consensus's
    /// own length (`-minRatioAgreement`). 1.0 means "at least as many copies
    /// want the change as want the status quo".
    pub min_ratio: f64,
    /// Flagged windows closer than this are merged into one region.
    /// `AutoRunBlocker` uses `window / 5`.
    pub merge_gap: usize,
}

impl Default for VoteParams {
    fn default() -> Self {
        VoteParams { window: 20, min_copies: 5, min_ratio: 1.0, merge_gap: 4 }
    }
}

impl VoteParams {
    /// `AutoRunBlocker`'s coupling of merge distance to window size.
    pub fn with_window(window: usize, min_copies: usize) -> Self {
        VoteParams { window, min_copies, min_ratio: 1.0, merge_gap: (window / 5).max(1) }
    }
}

/// A region where the copies vote for a different length than the consensus has.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VoteBlock {
    /// Inclusive alignment-column range.
    pub col_start: usize,
    pub col_end: usize,
    /// Inclusive consensus-position range (0-based, ungapped).
    pub cons_start: usize,
    pub cons_end: usize,
    /// The length the copies agree on, and how many agree.
    pub majority_len: usize,
    pub majority_count: usize,
    /// How many copies keep the consensus's own length.
    pub original_count: usize,
    /// Copies spanning the region at all.
    pub spanning: usize,
}

impl VoteBlock {
    /// Positive when the copies want the consensus *longer* here — the
    /// missing-ancestral-bases case this detector exists for.
    pub fn length_delta(&self) -> isize {
        self.majority_len as isize - (self.cons_end as isize - self.cons_start as isize + 1)
    }
}

/// Per-row prefix counts of real (non-gap, non-pad) bases, for O(1) range queries.
fn base_prefixes(msa: &MultiAlign) -> Vec<Vec<u32>> {
    msa.sequences
        .iter()
        .map(|row| {
            let mut p = Vec::with_capacity(row.seq.len() + 1);
            p.push(0);
            let mut n = 0;
            for &b in &row.seq {
                if !seqmod::is_gap(b) && b != b' ' {
                    n += 1;
                }
                p.push(n);
            }
            p
        })
        .collect()
}

/// Find regions where the copies agree on a length the consensus does not have.
///
/// A port of `AutoRunBlocker.pl`'s detection, without its per-window subprocess.
///
/// # Why this rather than a score profile
///
/// [`low_scoring_columns`] asks "is this stretch of alignment bad?", which
/// conflates substitution noise, thin coverage and misregistered indels. In a
/// reference-bootstrapped alignment the third dominates — measured MSA widths
/// run 3x the consensus length — so the profile is negative nearly everywhere
/// and Ruzzo-Tompa returns a few huge segments rather than local defects.
///
/// This asks a narrower question with a categorical answer: **do the copies
/// vote for a different length here?** A length vote cannot be swayed by
/// substitution noise, and because windows are anchored to *consensus
/// positions* rather than alignment columns, it is indifferent to how gappy the
/// alignment is — which is the property the score profile lacks.
///
/// Windows slide one consensus position at a time; those that pass are merged
/// when they lie within [`VoteParams::merge_gap`], because a single ancestral
/// insertion trips every window overlapping it and repairing them piecemeal
/// would produce inconsistent joins.
pub fn length_vote_blocks(msa: &MultiAlign, p: &VoteParams) -> Vec<VoteBlock> {
    let Some(reference) = msa.sequences.first() else {
        return Vec::new();
    };
    // Columns carrying a consensus base: the coordinate system of the windows.
    let cons_cols: Vec<usize> = reference
        .seq
        .iter()
        .enumerate()
        .filter(|(_, &b)| !seqmod::is_gap(b) && b != b' ')
        .map(|(i, _)| i)
        .collect();
    if cons_cols.len() < p.window || p.window == 0 {
        return Vec::new();
    }
    let prefixes = base_prefixes(msa);

    let mut hits: Vec<VoteBlock> = Vec::new();
    for j in 0..=(cons_cols.len() - p.window) {
        let c0 = cons_cols[j];
        let c1 = cons_cols[j + p.window - 1];

        let mut histo: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();
        let mut spanning = 0usize;
        for (ri, row) in msa.sequences.iter().enumerate().skip(1) {
            // Only copies that cover the whole window vote; a copy ending
            // inside it has no opinion about the length.
            if row.col_start > c0 || row.col_end <= c1 {
                continue;
            }
            let len = (prefixes[ri][c1 + 1] - prefixes[ri][c0]) as usize;
            *histo.entry(len).or_insert(0) += 1;
            spanning += 1;
        }
        if spanning == 0 {
            continue;
        }
        // Ties to the longer length, matching Refiner's determinism fix.
        let (maj_len, maj_count) = histo
            .iter()
            .map(|(&l, &c)| (l, c))
            .max_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)))
            .unwrap();
        let orig_count = histo.get(&p.window).copied().unwrap_or(0);

        if maj_len == p.window || maj_count < p.min_copies {
            continue;
        }
        // The ratio test is what separates "copies want a different length"
        // from "this window is merely noisy".
        if orig_count > 0 && (maj_count as f64 / orig_count as f64) < p.min_ratio {
            continue;
        }
        hits.push(VoteBlock {
            col_start: c0,
            col_end: c1,
            cons_start: j,
            cons_end: j + p.window - 1,
            majority_len: maj_len,
            majority_count: maj_count,
            original_count: orig_count,
            spanning,
        });
    }

    // Merge windows within `merge_gap`, then **recompute the vote over the
    // merged span**. Carrying a single window's statistic onto a merged region
    // is meaningless — a 20-position window's majority length says nothing
    // about a 600-position region, and comparing the two produced deltas like
    // "len 646 -> 0". `AutoRunBlocker` does the same thing for the same reason:
    // it re-runs Blocker on each clustered region rather than reusing the
    // per-window results.
    let mut spans: Vec<(usize, usize, usize, usize)> = Vec::new(); // cons0,cons1,col0,col1
    for h in hits {
        match spans.last_mut() {
            Some(prev) if h.cons_start <= prev.1 + p.merge_gap + 1 => {
                prev.1 = prev.1.max(h.cons_end);
                prev.3 = prev.3.max(h.col_end);
            }
            _ => spans.push((h.cons_start, h.cons_end, h.col_start, h.col_end)),
        }
    }

    let mut merged = Vec::new();
    for (cons0, cons1, col0, col1) in spans {
        let span_len = cons1 - cons0 + 1;
        let mut histo: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();
        let mut spanning = 0usize;
        for (ri, row) in msa.sequences.iter().enumerate().skip(1) {
            if row.col_start > col0 || row.col_end <= col1 {
                continue;
            }
            let len = (prefixes[ri][col1 + 1] - prefixes[ri][col0]) as usize;
            // A copy with no bases at all across the span is absent here, not
            // voting for "length zero" — counting it would let fragmented rows
            // that merely straddle the region dominate the vote.
            if len == 0 {
                continue;
            }
            *histo.entry(len).or_insert(0) += 1;
            spanning += 1;
        }
        if spanning == 0 {
            continue;
        }
        let (maj_len, maj_count) = histo
            .iter()
            .map(|(&l, &c)| (l, c))
            .max_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)))
            .unwrap();
        let orig_count = histo.get(&span_len).copied().unwrap_or(0);
        if maj_len == span_len || maj_count < p.min_copies {
            continue;
        }
        if orig_count > 0 && (maj_count as f64 / orig_count as f64) < p.min_ratio {
            continue;
        }
        merged.push(VoteBlock {
            col_start: col0,
            col_end: col1,
            cons_start: cons0,
            cons_end: cons1,
            majority_len: maj_len,
            majority_count: maj_count,
            original_count: orig_count,
            spanning,
        });
    }
    merged
}

/// Vote over one consensus-position span, as `resolveIndels.pl`'s
/// `evalMSABlock` does: only copies covering the whole span vote, a copy with
/// no bases at all in it is absent rather than voting zero, and ties go to the
/// longer length.
///
/// Returns `None` when nobody spans it, the majority already agrees with the
/// consensus, or the copy/ratio filters reject it.
fn vote_span(
    msa: &MultiAlign,
    prefixes: &[Vec<u32>],
    col0: usize,
    col1: usize,
    span_len: usize,
    min_copies: usize,
    min_ratio: f64,
    skip_empty: bool,
) -> Option<VoteBlock> {
    let mut histo: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut spanning = 0usize;
    for (ri, row) in msa.sequences.iter().enumerate().skip(1) {
        if row.col_start > col0 || row.col_end <= col1 {
            continue;
        }
        let len = (prefixes[ri][col1 + 1] - prefixes[ri][col0]) as usize;
        if skip_empty && len == 0 {
            continue;
        }
        *histo.entry(len).or_insert(0) += 1;
        spanning += 1;
    }
    if spanning == 0 {
        return None;
    }
    let (maj_len, maj_count) = histo
        .iter()
        .map(|(&l, &c)| (l, c))
        .max_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)))?;
    let orig_count = histo.get(&span_len).copied().unwrap_or(0);
    if maj_len == span_len || maj_count < min_copies {
        return None;
    }
    if orig_count > 0 && (maj_count as f64 / orig_count as f64) < min_ratio {
        return None;
    }
    Some(VoteBlock {
        col_start: col0,
        col_end: col1,
        cons_start: 0,
        cons_end: 0,
        majority_len: maj_len,
        majority_count: maj_count,
        original_count: orig_count,
        spanning,
    })
}

/// How strongly the copies prefer the new length over the consensus's own.
///
/// `resolveIndels.pl`'s ranking key: `numMajorityLen / numOrigLen`, with a
/// denominator of 0.1 when no copy keeps the consensus length — so a region
/// nothing agrees with sorts above every merely-contested one.
pub fn vote_ratio(v: &VoteBlock) -> f64 {
    v.majority_count as f64 / (v.original_count as f64).max(0.1)
}

/// Every window that votes for a different length, over **several window sizes
/// at once**.
///
/// `AutoRunBlocker.pl` takes a single `-windowSize`; `resolveIndels.pl` takes
/// `-discrete_windows 7,15,24,5` or a `-window_min/-window_max` range and runs
/// all of them at every consensus position, leaving the choice between a strong
/// narrow signal and a weak wide one to the aggregation step rather than to a
/// parameter. Hits are returned unmerged and unranked — feed them to
/// [`tile_hits`] or [`cluster_hits`].
pub fn multi_window_hits(
    msa: &MultiAlign,
    sizes: &[usize],
    min_copies: usize,
    min_ratio: f64,
) -> Vec<VoteBlock> {
    let Some(reference) = msa.sequences.first() else {
        return Vec::new();
    };
    let cons_cols: Vec<usize> = reference
        .seq
        .iter()
        .enumerate()
        .filter(|(_, &b)| !seqmod::is_gap(b) && b != b' ')
        .map(|(i, _)| i)
        .collect();
    let prefixes = base_prefixes(msa);

    let mut hits: Vec<VoteBlock> = Vec::new();
    for &w in sizes {
        if w == 0 || cons_cols.len() < w {
            continue;
        }
        for j in 0..=(cons_cols.len() - w) {
            let c0 = cons_cols[j];
            let c1 = cons_cols[j + w - 1];
            if let Some(mut v) =
                vote_span(msa, &prefixes, c0, c1, w, min_copies, min_ratio, false)
            {
                v.cons_start = j;
                v.cons_end = j + w - 1;
                hits.push(v);
            }
        }
    }
    hits
}

/// Greedy tiling path over candidate regions, strongest vote first.
///
/// `resolveIndels.pl -aggregation_method tile`: sort every candidate by
/// [`vote_ratio`] descending and take it unless it comes within `min_sep`
/// consensus positions of one already taken. Where the fixed-precedence
/// policies in the experiment let the *selection* decide which of two
/// overlapping proposals wins, this lets the copies decide.
pub fn tile_hits(mut hits: Vec<VoteBlock>, min_sep: usize) -> Vec<VoteBlock> {
    // Ratio descending; ties to the wider region, then to position, so the
    // result does not depend on the order the windows were generated in.
    hits.sort_by(|a, b| {
        vote_ratio(b)
            .partial_cmp(&vote_ratio(a))
            .unwrap_or(std::cmp::Ordering::Equal)
            .then((b.cons_end - b.cons_start).cmp(&(a.cons_end - a.cons_start)))
            .then(a.cons_start.cmp(&b.cons_start))
    });
    let mut kept: Vec<VoteBlock> = Vec::new();
    for h in hits {
        let lo = h.cons_start.saturating_sub(min_sep);
        let hi = h.cons_end + min_sep;
        if kept.iter().any(|k| lo <= k.cons_end && k.cons_start <= hi) {
            continue;
        }
        kept.push(h);
    }
    kept.sort_by_key(|k| k.cons_start);
    kept
}

/// Guarded left-to-right clustering of candidate regions.
///
/// `resolveIndels.pl`'s `clusterBlocks`: walk the candidates in start order and
/// extend the current cluster only while the vote taken over the *extended*
/// span still passes. When it stops passing, emit the cluster as it stood and
/// skip whatever else overlaps it.
///
/// The guard is the point. `length_vote_blocks` merges everything within
/// `merge_gap` and re-votes once at the end, which is how a family-spanning
/// region that deleted 1,467 consensus bases got proposed; here the merge stops
/// at the first extension the copies do not support, so no width cap is needed
/// to make it safe.
pub fn cluster_hits(
    msa: &MultiAlign,
    mut hits: Vec<VoteBlock>,
    gap: usize,
    min_copies: usize,
    min_ratio: f64,
) -> Vec<VoteBlock> {
    if hits.is_empty() {
        return Vec::new();
    }
    // Start ascending, longest first — the Perl's sort.
    hits.sort_by(|a, b| a.cons_start.cmp(&b.cons_start).then(b.cons_end.cmp(&a.cons_end)));
    let prefixes = base_prefixes(msa);

    let mut out: Vec<VoteBlock> = Vec::new();
    let mut i = 0usize;
    while i < hits.len() {
        let mut cluster = hits[i].clone();
        i += 1;
        while i < hits.len() && hits[i].cons_start <= cluster.cons_end + gap + 1 {
            let cand_end = hits[i].cons_end.max(cluster.cons_end);
            let span = cand_end - cluster.cons_start + 1;
            match vote_span(
                msa,
                &prefixes,
                cluster.col_start,
                hits[i].col_end.max(cluster.col_end),
                span,
                min_copies,
                min_ratio,
                true,
            ) {
                Some(mut grown) => {
                    grown.cons_start = cluster.cons_start;
                    grown.cons_end = cand_end;
                    cluster = grown;
                    i += 1;
                }
                // The extension is not supported: emit what we have and skip
                // the rest of the overlapping run, as the Perl does.
                None => {
                    while i < hits.len() && hits[i].cons_start <= cluster.cons_end + gap {
                        i += 1;
                    }
                    break;
                }
            }
        }
        out.push(cluster);
    }
    out
}
