//! `autocons` — build the N best-scoring consensus sequences from a set of
//! related sequences.
//!
//! A port of GIRI `acons/src/autocons.cpp`. Two phases:
//!
//! 1. **Reference selection.** Every input sequence is tried as the reference:
//!    align all the others to it, assemble an MSA, call a consensus. Candidates
//!    are ranked by the *total alignment score* — not by any property of the
//!    consensus itself. This is `MultipleAlignment::autoConsensus`, whose
//!    `outscore` is exactly `construct`'s return value.
//!
//! 2. **Refinement.** Each of the top N consensi is re-aligned against the whole
//!    input and re-called, up to [`Params::iterations`] extra passes, stopping
//!    early once the consensus stops changing.
//!
//! # One subtlety worth knowing
//!
//! The two phases treat the reference row differently, and this is faithful to
//! the C++:
//!
//! * Phase 1 calls `autoConsensus` with its default `withRefSeq = true`, so the
//!   reference **is** counted in the profile.
//! * Phase 2 does `maln.erase(maln.begin())` before calling, so it **is not**.
//!
//! Phase 1's consensus is not discarded — it is the starting point for phase 2 —
//! so the difference is load-bearing.

pub mod mhsp;
pub mod probe;
pub mod trim;
pub mod lowqual;

use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use aln_core::consensus::ConsensusParams;

/// The fastest pairwise aligner this build has: parasail's striped SIMD
/// kernels on x86, which is the only architecture they are vendored for, and
/// the scalar `ReferenceAligner` everywhere else. Both implement the same
/// trait and produce the same alignments; only the speed differs.
#[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
pub use aln_parasail::ParasailAligner as FastAligner;
#[cfg(not(any(target_arch = "x86_64", target_arch = "x86")))]
pub use aln_reference::ReferenceAligner as FastAligner;
use aln_core::msa::{assemble_msa, InsertionPolicy, MsaMember, MultiAlign};
use aln_core::seq::Strand;
use aln_coord::Span;
use aln_core::{giri, seq as seqmod, Alignment, Sequence};
use aln_engine::{PairwiseAligner, Result};
use rayon::prelude::*;

/// Tuning for a run.
#[derive(Debug, Clone)]
pub struct Params {
    /// How many consensi to emit. A value below 1.0 is a *fraction* of the
    /// input count, rounded down, with a floor of 1 — matching `--num`.
    pub num_consensi: f64,
    /// Base name for emitted sequences. With more than one, a 1-based counter
    /// is appended.
    pub base_name: String,
    /// Extra refinement passes after the first. The C++ uses `iter = 3`, giving
    /// up to four total passes.
    pub iterations: usize,
    pub consensus: ConsensusParams,

    /// Which consensus caller to use.
    pub caller: Caller,

    /// `acons --min`: minimum non-gap residues in a column, below which the
    /// column is called a gap whatever the scores say.
    ///
    /// Both callers honour it, through [`gapped_consensus`], and it binds in
    /// phase 2 only — phase 1 keeps the reference row and calls with a floor of
    /// 1, as the C++ does.
    ///
    /// 0 disables it, and is the default: the Dfam caller has never gated its
    /// columns, and turning that on for every run would shorten consensi
    /// curators have already built and checked without it. `--orig` restores
    /// the C++'s 2.
    pub min_non_gap_count: usize,

    /// Run GIRI's species-aware CpG restoration on the final consensus —
    /// `acons --mam`.
    ///
    /// Applies only with [`Caller::Giri`], matching the C++'s
    /// `if (useOrig && species)`. The Dfam caller does its own CpG restoration
    /// on every pass, so the flag is meaningless there.
    pub restore_cpg: bool,

    /// Repair low-quality blocks between two rounds of refinement.
    ///
    /// `Refiner`'s `resolveLowQualityBlocks` step: refine to a stable answer,
    /// re-derive the consensus over locally-bad stretches directly from the
    /// instances, then refine again from the patched consensus. Off by default
    /// — it is a change to the workflow, not to the aligner.
    pub repair_blocks: bool,

    /// Ruzzo-Tompa segment-score threshold for calling a block low-quality.
    /// `MultAln.pm` defaults to 1.
    pub repair_threshold: f64,

    /// Matrix used to score alignment columns during block repair.
    ///
    /// Carried here rather than passed alongside because [`run`] drives both
    /// phases and would otherwise have to thread it through unrelated code.
    /// `None` disables repair even when `repair_blocks` is set — there is no
    /// sensible default matrix to invent at this layer.
    pub repair_matrix: Option<aln_core::SubstMatrix>,

    /// Also take repair candidates from a scanning window of this many
    /// consensus positions, resolved by the same length vote. 0 disables it.
    ///
    /// `AutoRunBlocker`'s selection, added alongside `Refiner`'s. Measured on
    /// 791 hs1 families against curated Dfam consensi, it roughly doubles what
    /// the repair recovers (+99 -> +223 net bases closer to curation) because
    /// 44% of what it finds sits in a low-quality block whose *whole-block*
    /// vote says the consensus already agrees — the disagreement is only
    /// visible at window scale.
    pub repair_window: usize,

    /// Which statistic decides whether a repair is kept.
    pub repair_accept: AcceptRule,

    /// Re-derive spans the consensus has no bases for, every pass, so inherited
    /// insertions can enter the consensus and the reference can grow.
    pub pack_insertions: bool,
    /// Merge gap runs separated by at most this many called columns.
    pub pack_max_sep: usize,
    /// A span needs one instance contributing at least this many bases.
    pub pack_min_seg: usize,
    /// The all-against-all winner must beat this summed score, so a span whose
    /// instances do not align to one another is left alone.
    pub pack_min_score: i64,
    /// A span is packed only if at least this fraction of the copies spanning
    /// it carry bases inside it. 0.5 makes packing consistent with the column
    /// caller: an insertion most spanning copies carry survives the
    /// re-alignment that follows, a minority one would not.
    pub pack_min_occupancy: f64,
    /// Keep bases the re-derivation's centre lacks, so a majority-carried
    /// insertion can survive rather than being projected away.
    pub pack_keep_insertions: bool,

    /// Re-derive a low-quality block from an all-against-all alignment of its
    /// instances when no single instance length dominates it.
    ///
    /// `Refiner`'s `else` branch. On by default: skipping it leaves every
    /// no-majority block unrepaired.
    pub repair_all_vs_all: bool,

    /// Score reference candidates by non-redundant coverage per instance
    /// rather than by a plain sum. See [`nonredundant_candidate_score`].
    ///
    /// **On by default.** A plain sum credits a candidate twice for a region
    /// one copy covers twice, which rewards attracting fragmented, overlapping
    /// alignments. Measured on 24 heavily fragmented simulated families it
    /// discounts 7.4% of the typical candidate's score and changed the chosen
    /// reference in none of them — redundancy concentrates in candidates that
    /// were already losing. It is on because it is the fairer statistic, not
    /// because it was shown to change outcomes.
    pub nonredundant_reference_score: bool,

    /// How insertions are merged into the alignment.
    ///
    /// [`InsertionPolicy::GrowIncremental`] reproduces GIRI's
    /// `adjustReference`; [`InsertionPolicy::GrowPerSlot`] keeps each member's
    /// insertions independent. They place the same bases but can differ in
    /// width — see [`aln_core::msa::InsertionPolicy`].
    pub insertions: InsertionPolicy,
}

/// Which consensus caller drives a run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Caller {
    /// The Dfam caller — per-column argmax with prefer-`N` ties, plus CpG
    /// restoration. The default here and in the current C++.
    #[default]
    Dfam,
    /// The original GIRI caller — `acons --orig`.
    ///
    /// Note the C++ is *asymmetric* about this: phase 1 always calls GIRI, and
    /// only phase 2 honours `--orig`. Selecting [`Caller::Giri`] here uses GIRI
    /// in **both** phases, which is what makes a like-for-like comparison
    /// against `autocons --orig` possible — under `--orig` the C++ is
    /// all-GIRI too.
    Giri,
}

/// The phase-2 gapped consensus over `msa`, in its own column coordinates.
///
/// Both callers drop the reference row and both honour
/// [`Params::min_non_gap_count`], but they reach the gate differently. GIRI
/// applies it inside its own argmax, so no candidate ever wins a sparse
/// column. The Dfam caller has no gate, so [`gate_by_occupancy`] masks what it
/// called. The orderings differ only where a gated column adjoins a kept one
/// and the CpG pass had paired the two: GIRI restores CpG after gating and
/// never sees that pair, while here the CpG pass resolved it before the gate
/// removed half of it.
///
/// Pack insertion spans *after* this, as both phases do. Packing re-derives a
/// span only where two instances carry sequence across it, so what it puts
/// back is supported by more than the row the gate was aimed at.
pub fn gapped_consensus(msa: &MultiAlign, params: &Params) -> Vec<u8> {
    match params.caller {
        Caller::Dfam => {
            let call = ConsensusParams { include_reference: false, ..params.consensus.clone() };
            let mut cons = msa.consensus(&call);
            gate_by_occupancy(&mut cons, msa, params.min_non_gap_count);
            cons
        }
        Caller::Giri => msa.giri_consensus(params.min_non_gap_count),
    }
}

/// Force to a gap any column where fewer than `min` rows carry a residue.
///
/// Counts what `giri::get_consensus` counts: anything that is neither a gap
/// nor padding, in either the Dfam or the GIRI convention. Skips the reference
/// row, as the Dfam call it corrects does. A row that ends before a column
/// counts as padded there.
fn gate_by_occupancy(cons: &mut [u8], msa: &MultiAlign, min: usize) {
    if min == 0 {
        return;
    }
    let rows = &msa.sequences[1.min(msa.sequences.len())..];
    for (col, c) in cons.iter_mut().enumerate() {
        let covered = rows
            .iter()
            .filter(|r| r.seq.get(col).is_some_and(|&b| !seqmod::is_structural(b)))
            .count();
        if covered < min {
            *c = b'-';
        }
    }
}

impl Default for Params {
    fn default() -> Self {
        Params {
            num_consensi: 1.0,
            base_name: "CON".to_string(),
            iterations: 3,
            consensus: ConsensusParams::default(),
            caller: Caller::Dfam,
            // The C++ default for --min.
            min_non_gap_count: 0,
            restore_cpg: false,
            repair_blocks: false,
            repair_accept: AcceptRule::default(),
            nonredundant_reference_score: true,
            repair_all_vs_all: true,
            pack_insertions: false,
            pack_max_sep: 4,
            pack_min_seg: 5,
            pack_min_score: 0,
            pack_min_occupancy: 0.5,
            pack_keep_insertions: true,
            repair_threshold: crate::lowqual::THRESHOLD,
            repair_matrix: None,
            repair_window: 0,
            insertions: InsertionPolicy::GrowIncremental,
        }
    }
}

/// One candidate reference, scored.
#[derive(Debug, Clone)]
pub struct Candidate {
    /// Index into the input set.
    pub index: usize,
    /// Total alignment score of this reference against every other sequence.
    pub score: i64,
    /// Consensus called from that alignment, ungapped.
    pub consensus: Vec<u8>,
}

/// Why refinement stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StopReason {
    /// A pass reproduced its input: a fixed point. The C++'s only stop signal.
    Converged,
    /// A pass produced a consensus seen earlier in this run — a cycle of length
    /// >1, so further passes would repeat it forever. Iterating longer cannot
    /// help, and at `--iterations 11` that is most of the budget wasted.
    Cycled,
    /// Ran out of passes without either. The consensus may still be improving.
    Exhausted,
}

impl StopReason {
    /// True when refinement reached a stable answer rather than a deadline.
    pub fn is_stable(self) -> bool {
        matches!(self, StopReason::Converged | StopReason::Cycled)
    }
}

/// One refinement pass, for forensics.
///
/// Refinement is chosen by score, but score is a *self-referential* objective —
/// how well the consensus explains the instances, not how close it is to the
/// truth. Keeping every pass lets the two be compared directly wherever a truth
/// is known. Cost is one consensus per pass, not one MSA.
#[derive(Debug, Clone)]
pub struct PassTrace {
    /// 1-based pass number.
    pub pass: usize,
    /// Summed alignment score of every instance against this pass's *input*.
    pub score: i64,
    /// The same alignments under [`mean_score_per_base`].
    pub norm_score: f64,
    /// The consensus this pass produced.
    pub consensus: Vec<u8>,
    /// Mean Kimura divergence of the instance rows against this pass's
    /// *input* consensus, plain and CpG-adjusted.
    ///
    /// Recorded per pass rather than only at the end so the first pass
    /// describes the bootstrap consensus — the state before any refinement —
    /// without having to keep a second MSA alive to measure it later.
    pub divergence: f64,
    pub divergence_cpg: f64,
    /// MSA rows this pass built, excluding the reference row.
    ///
    /// Under `--hsps tiled` one instance can hold several rows, so this exceeds
    /// [`instances`](Self::instances) whenever any instance tiled. Reporting
    /// rows as if they were instances is what produced "62 of 60 instances".
    pub rows: usize,
    /// Distinct instances contributing at least one row.
    ///
    /// Counted from the alignment indices, so it is bounded by the input count
    /// however the rows are split.
    pub instances: usize,
}

/// Mean over copies of (alignment score / aligned copy bases).
///
/// The alternative to the summed score as an accept criterion. The sum is
/// dominated by whichever copies are longest and by how many clear the score
/// floor, and it is bought most cheaply by *removing* a consensus base that a
/// majority of copies lack: with `comparison.matrix` (mean match +9.5, gap open
/// -40) that returns roughly +40 per copy no longer forced to gap, against
/// about -50 charged to the minority that carried the base. One low-occupancy
/// column is therefore worth several times its own length in summed score, and
/// no repair that deletes one can be rejected.
///
/// Averaging a per-base quality over copies gives every copy the same weight
/// and bounds what any one of them can contribute, so a repair that helps a
/// large majority slightly at the cost of wrecking a few alignments no longer
/// scores the same as one that helps everybody.
pub fn mean_score_per_base(alignments: &[(usize, aln_core::Alignment)]) -> f64 {
    if alignments.is_empty() {
        return 0.0;
    }
    let total: f64 = alignments
        .iter()
        .map(|(_, a)| {
            let bases = a.query_end.saturating_sub(a.query_start).max(1);
            a.score as f64 / bases as f64
        })
        .sum();
    total / alignments.len() as f64
}

/// Which statistic [`refine_with_repair`] compares before keeping a repair.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AcceptRule {
    /// Summed alignment score — `Refiner`'s `SCORE=`, and the historical
    /// behaviour.
    #[default]
    Total,
    /// [`mean_score_per_base`].
    MeanPerBase,
    /// Keep the repair unconditionally — no gate at all.
    ///
    /// Not a shipping option. It exists so the gate's contribution can be
    /// measured against the same proposal set rather than against a different
    /// tool's.
    Always,
}

/// A refined consensus, ready to emit.
#[derive(Debug, Clone)]
pub struct Refined {
    pub name: String,
    /// Ungapped consensus.
    pub consensus: Vec<u8>,
    /// Which input sequence seeded this consensus, and its phase-1 score.
    ///
    /// Reference selection is a discrete argmax over summed alignment scores,
    /// so a small scoring difference can flip the choice and yield a materially
    /// different consensus. Anyone debugging an unexpected result needs to know
    /// which sequence won, and by how much.
    pub reference: Reference,
    /// Total alignment score from the final pass.
    pub score: i64,
    /// The same pass under [`mean_score_per_base`].
    pub norm_score: f64,
    /// The alignment the final consensus was called from.
    pub msa: MultiAlign,
    /// Passes actually run.
    pub passes: usize,
    /// True if it stopped because the consensus stopped changing rather than
    /// because it ran out of passes. [`StopReason`] carries the full answer.
    pub converged: bool,
    /// Why the loop stopped.
    pub stop: StopReason,
    /// Every pass, in order. See [`PassTrace`].
    pub trace: Vec<PassTrace>,
    /// Free-text provenance about how this family was built, for the record
    /// that gets written out.
    ///
    /// A consensus that was *not* extended because the extension ran away is a
    /// different object from one that was never offered an extension, and a
    /// curator looking at the family months later cannot tell them apart from
    /// the sequence. Whatever the pipeline decided on this family's behalf is
    /// recorded here and travels with it into the Stockholm.
    pub notes: Vec<String>,
}


/// Candidate score counting each instance's coverage of the reference once.
///
/// The overlap filter deletes a lower-scoring alignment only when less than
/// `100 - mask_level` percent of it is novel, so an alignment survives with up
/// to ~80% of its reference span already covered — and then contributes its
/// *whole* score. A candidate that attracts fragmented, overlapping alignments
/// is over-credited against one that attracts clean single alignments, and
/// reference selection is a discrete argmax, so that can decide the run.
///
/// Each alignment is taken in score order and credited only for the fraction of
/// its reference span not already claimed **by the same instance**. Coverage of
/// one reference region by *different* copies is left alone: that is what a
/// good reference looks like, not double counting.
fn nonredundant_candidate_score(alignments: &[(usize, Alignment)]) -> i64 {
    let mut by: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();
    for (pos, (inst, _)) in alignments.iter().enumerate() {
        by.entry(*inst).or_default().push(pos);
    }
    let mut total = 0i64;
    for (_, mut idx) in by {
        idx.sort_by(|&a, &b| {
            alignments[b].1.score.cmp(&alignments[a].1.score).then(
                alignments[a].1.subj_start.cmp(&alignments[b].1.subj_start),
            )
        });
        let mut claimed: Vec<(usize, usize)> = Vec::new();
        for p in idx {
            let a = &alignments[p].1;
            let (s, e) = (a.subj_start, a.subj_end.max(a.subj_start));
            let span = e.saturating_sub(s);
            if span == 0 {
                total += a.score as i64;
                continue;
            }
            let overlap: usize = claimed
                .iter()
                .map(|&(cs, ce)| ce.min(e).saturating_sub(cs.max(s)))
                .sum();
            let novel = span.saturating_sub(overlap.min(span));
            total += ((a.score as f64) * (novel as f64 / span as f64)).round() as i64;
            claimed.push((s, e));
            claimed.sort_unstable();
            // Merge so the next overlap sum cannot double-count claimed space.
            let mut merged: Vec<(usize, usize)> = Vec::with_capacity(claimed.len());
            for &(cs, ce) in &claimed {
                match merged.last_mut() {
                    Some(last) if cs <= last.1 => last.1 = last.1.max(ce),
                    _ => merged.push((cs, ce)),
                }
            }
            claimed = merged;
        }
    }
    total
}

/// Where alignments come from.
///
/// `autocons` needs exactly two operations, and the two backend families supply
/// them very differently:
///
/// * A [`PairwiseAligner`] does one pair at a time; "reference against all" is a
///   loop, and reference selection is that loop repeated per candidate.
/// * A search engine does "one query set against one subject set" in a single
///   call, so both operations are one call each — and reference selection is a
///   single all-against-all search rather than `n` searches.
///
/// This trait lets the phase-1 and phase-2 logic stay identical while the work
/// underneath is shaped for the backend.
pub trait AlignmentSource: Sync {
    /// Alignments of every sequence in `seqs` against `reference`.
    ///
    /// `skip` suppresses one index so a reference is not aligned to itself.
    /// A backend may return **more than one alignment per index** — rmblast
    /// does, for instances that match in several pieces.
    fn against_reference(
        &self,
        reference: &Sequence,
        seqs: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>>;

    /// Alignments for every candidate reference, indexed by candidate.
    ///
    /// The default repeats [`against_reference`](Self::against_reference) per
    /// candidate, which is right for a pairwise aligner. A search engine
    /// overrides it with a single all-against-all call.
    fn against_every_reference(
        &self,
        seqs: &[Sequence],
    ) -> Result<Vec<Vec<(usize, Alignment)>>> {
        seqs.iter()
            .enumerate()
            .map(|(i, r)| {
                if r.is_empty() {
                    Ok(Vec::new())
                } else {
                    self.against_reference(r, seqs, Some(i))
                }
            })
            .collect()
    }

    /// As [`against_reference`](Self::against_reference), but called from inside
    /// a loop that is already parallel over candidates.
    ///
    /// Phase 1 fans out one task per candidate reference, which saturates the
    /// machine on its own; a backend that would otherwise parallelise
    /// internally should not, to avoid nesting for no gain. The default simply
    /// delegates.
    fn against_reference_nested(
        &self,
        reference: &Sequence,
        seqs: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>> {
        self.against_reference(reference, seqs, skip)
    }

    /// True when [`against_every_reference`](Self::against_every_reference) is a
    /// single batched call, so phase 1 should not also parallelise over
    /// candidates.
    fn batches_all_vs_all(&self) -> bool {
        false
    }
}

/// Where the rayon pool is used inside [`align_all`].
///
/// Phase 1 already has one task per candidate reference, which saturates any
/// realistic core count once the input runs to dozens of sequences, so its
/// inner loop runs [`Inner::Sequential`]. Phase 2 has only one reference in
/// flight at a time and needs [`Inner::Parallel`] to use the machine at all.
///
/// This mirrors the C++, which sets `ThreadedAligner::setMultithreaded(false)`
/// around `process_with_pthreads` and restores it for the refinement loop —
/// except that rayon work-steals over one pool, so nesting would merely add
/// overhead rather than explode the thread count.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Inner {
    Parallel,
    Sequential,
}

/// Align every sequence in `seqs` to `reference`, keeping input indices.
///
/// `skip` suppresses one index, used in phase 1 so a reference is not aligned to
/// itself. The C++ skips by pointer identity, which silently fails to skip a
/// duplicated sequence stored as a distinct object; skipping by index is exact.
fn align_all<A: PairwiseAligner>(
    aligner: &A,
    reference: &Sequence,
    seqs: &[Sequence],
    skip: Option<usize>,
    inner: Inner,
    prepass: &Prepass,
) -> Result<Vec<(usize, Alignment)>> {
    let profile = aligner.prepare_subject(reference)?;

    // Score-then-trace.  parasail's score-only kernels skip the O(mn) traceback
    // matrix and run 3.7-6.5x faster, so when most pairs fall below
    // `min_score` it pays to score first and trace only the survivors.
    //
    // When most pairs *pass* it is pure overhead — the survivors get scored
    // twice.  Which regime we are in depends on the cutoff and the family, and
    // is not knowable in advance: at GIRI's default cutoff of 1 essentially
    // nothing is rejected, at a noise-calibrated cutoff around 60% is.  So
    // probe the first few pairs of the run and keep the prepass only if it is
    // paying off.  The counters live on `Pairwise`, so this is calibrated once
    // per run rather than once per reference.
    //
    // This never changes the result — `score_prepared` and `align_prepared`
    // agree by construction, pinned by aln-parasail's
    // `score_only_kernels_agree_with_the_traceback_kernels`.  Under
    // `Inner::Parallel` which pairs land in the probe varies between runs;
    // only the timing does.

    let one = |i: usize, s: &Sequence| -> Result<Option<(usize, Alignment)>> {
        if Some(i) == skip || s.is_empty() {
            return Ok(None);
        }
        if prepass.wanted() {
            let survived = aligner.score_prepared(&profile, s)?.is_some();
            prepass.observe(survived);
            if !survived {
                return Ok(None);
            }
        }
        Ok(aligner.align_prepared(&profile, s)?.map(|a| (i, a)))
    };

    let per: Vec<Result<Option<(usize, Alignment)>>> = match inner {
        Inner::Parallel => seqs.par_iter().enumerate().map(|(i, s)| one(i, s)).collect(),
        Inner::Sequential => seqs.iter().enumerate().map(|(i, s)| one(i, s)).collect(),
    };

    let mut out = Vec::with_capacity(per.len());
    for r in per {
        if let Some(pair) = r? {
            out.push(pair);
        }
    }
    Ok(out)
}

/// Wraps a [`PairwiseAligner`] as an [`AlignmentSource`].
///
/// A blanket `impl<A: PairwiseAligner>` would be tempting, but it forecloses a
/// concrete impl for any type that might one day also be a `PairwiseAligner` —
/// coherence cannot know that it will not. An explicit wrapper keeps both
/// families available.
pub struct Pairwise<A> {
    aligner: A,
    prepass: Prepass,
}

/// Score-prepass state, shared across every call for a whole run.
///
/// Phase 1 issues O(n^2) alignments against n candidate references, so probing
/// afresh per reference would re-pay the sampling cost n times — measured at
/// ~3% on a 100-sequence family where the prepass never pays. Calibrating once
/// amortises it away.
#[derive(Default)]
struct Prepass {
    enabled: bool,
    probed: AtomicUsize,
    passed: AtomicUsize,
}

impl Prepass {
    /// Whether this pair should be scored before being traced.
    ///
    /// Keeps the prepass only while it is paying off: the first `PROBE` pairs
    /// of a run are always sampled, and after that it stays on only if fewer
    /// than `KEEP_BELOW` of the sample cleared `min_score`.
    fn wanted(&self) -> bool {
        const PROBE: usize = 64;
        const KEEP_BELOW: f64 = 0.75;
        if !self.enabled {
            return false;
        }
        let n = self.probed.load(Ordering::Relaxed);
        n < PROBE || (self.passed.load(Ordering::Relaxed) as f64) < KEEP_BELOW * n as f64
    }

    /// Fold one probe result into the sample, while the sample is still open.
    fn observe(&self, survived: bool) {
        if self.probed.load(Ordering::Relaxed) < 64 {
            self.probed.fetch_add(1, Ordering::Relaxed);
            if survived {
                self.passed.fetch_add(1, Ordering::Relaxed);
            }
        }
    }
}

impl<A> Pairwise<A> {
    pub fn new(aligner: A) -> Self {
        Pairwise { aligner, prepass: Prepass::default() }
    }

    /// Score every pair before tracing it, and trace only the survivors.
    ///
    /// parasail's score-only kernels skip the O(mn) traceback matrix and run
    /// 3.7-6.5x faster, so where most pairs fall below `min_score` this is a
    /// clear win — and it caps peak memory, since only survivors ever allocate
    /// a traceback matrix.
    ///
    /// **Off by default, because on the Dfam family corpus it loses.** It needs
    /// two things at once: sequences long enough that the traceback matrix
    /// dominates per-call overhead, and a low acceptance rate. Measured across
    /// 791 families those are anti-correlated — median acceptance is 34% below
    /// 200 bp but 74% above 1500 bp — so enabling it cost 4% overall (92.9s ->
    /// 96.5s at `--min-score 95`) even though individual short low-acceptance
    /// families ran 1.7-2.7x faster in ratio terms on a tiny absolute base.
    ///
    /// Worth enabling for long, heterogeneous input run at a high cutoff, or
    /// when peak memory rather than time is the binding constraint.
    pub fn with_score_prepass(mut self) -> Self {
        self.prepass.enabled = true;
        self
    }

    /// The wrapped aligner.
    pub fn inner(&self) -> &A {
        &self.aligner
    }
}

impl<A: PairwiseAligner> AlignmentSource for Pairwise<A> {
    fn against_reference(
        &self,
        reference: &Sequence,
        seqs: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>> {
        // Parallel here: a pairwise aligner has nothing better to batch.
        align_all(&self.aligner, reference, seqs, skip, Inner::Parallel, &self.prepass)
    }

    fn against_reference_nested(
        &self,
        reference: &Sequence,
        seqs: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>> {
        align_all(&self.aligner, reference, seqs, skip, Inner::Sequential, &self.prepass)
    }
}

/// Assemble an MSA from alignments against `reference`.
///
/// Insertions are always kept — dropping them would make them uncallable — but
/// *how* they are merged is [`Params::insertions`].
/// Assemble an MSA from one reference and its alignments.
///
/// Public so a diagnostic can rebuild the **bootstrap-instance-induced**
/// alignment — every copy pairwise-aligned to the single real instance phase 1
/// chose. That is the state insertion packing would act on. A refined MSA is
/// no substitute: iterations and repairs reshape where gaps sit, so measuring
/// gap runs there describes the pipeline's own handiwork rather than what a
/// real instance's pairwise alignments produce.
pub fn build_msa(
    reference: &Sequence,
    seqs: &[Sequence],
    alignments: &[(usize, Alignment)],
    policy: InsertionPolicy,
) -> Result<MultiAlign> {
    // The gapped rows have to outlive the borrowed `MsaMember`s.
    let mut rows: Vec<(usize, Vec<u8>, Vec<u8>, &Alignment)> = Vec::with_capacity(alignments.len());
    for (i, a) in alignments {
        let (gq, gs) = a.gapped(&seqs[*i].seq, &reference.seq)?;
        rows.push((*i, gq, gs, a));
    }

    let members: Vec<MsaMember<'_>> = rows
        .iter()
        .map(|(i, gq, gs, a)| MsaMember {
            name: &seqs[*i].name,
            gapped_query: gq,
            gapped_reference: gs,
            ref_start: a.subj_start,
            span: Some(
                Span::new(a.query_start as u64, a.query_end as u64)
                    .expect("Alignment::validate keeps query_start <= query_end"),
            ),
            orient: a.strand,
        })
        .collect();

    Ok(assemble_msa(
        &reference.seq,
        &reference.name,
        &members,
        policy,
    )?)
}

/// Phase 1: score every sequence as a candidate reference.
///
/// Returns candidates sorted best-first, ties broken by input order. The C++
/// keys a `std::map<double, …>` and nudges colliding scores by `0.001`, which
/// makes tie order depend on insertion sequence; sorting explicitly is
/// deterministic and equivalent otherwise.
/// Parallelism note: one rayon task per candidate, each doing `n` alignments.
/// Scoring needs only those alignments, so the MSA and consensus are no longer
/// built here for every candidate — see [`score_candidates_ranked`]. They used
/// to be, and the cost was not what the old note here claimed: bounding peak
/// memory by worker-thread count still means one MSA per thread, and on sixty
/// 20 kb instances across 44 threads that reached 15 GB with 59 of the 60
/// results discarded.
pub fn score_candidates<A: AlignmentSource>(
    aligner: &A,
    seqs: &[Sequence],
    params: &Params,
) -> Result<Vec<Candidate>> {
    score_candidates_ranked(aligner, seqs, params, usize::MAX)
}

/// Rank candidate references, calling a consensus for the best `consensus_for`.
///
/// The score comes from the alignments alone. Building the multiple alignment
/// and calling a consensus is the expensive half — on a family of sixty 20 kb
/// instances it is most of the bootstrap, and it holds every candidate's MSA in
/// memory at once — so it is worth doing only for the candidates the caller
/// will use. [`run`] uses one. `Candidate::consensus` is left empty for the
/// rest; the scores, the ordering and the count are unaffected, so a caller
/// reading only `score` sees no difference.
pub fn score_candidates_ranked<A: AlignmentSource>(
    aligner: &A,
    seqs: &[Sequence],
    params: &Params,
    consensus_for: usize,
) -> Result<Vec<Candidate>> {
    // Phase 1 keeps the reference in the profile (`withRefSeq = true`), and
    // uses min_non_gap_count = 1 — the C++ calls the static `autoConsensus`
    // without the argument, taking its default, rather than passing `--min`.
    let phase1 = ConsensusParams { include_reference: true, ..params.consensus.clone() };

    let batched = aligner.batches_all_vs_all();
    // One call for every candidate when the backend can batch; otherwise the
    // default fans out per candidate below.
    let prefetched: Option<Vec<Vec<(usize, Alignment)>>> =
        if batched { Some(aligner.against_every_reference(seqs)?) } else { None };

    // Pass one: score only. The prefetched alignments are borrowed rather than
    // cloned — a clone per candidate duplicated the whole all-vs-all result
    // sixty times over on the family that prompted this.
    let scored: Vec<Result<Option<(usize, i64)>>> = seqs
        .par_iter()
        .enumerate()
        .map(|(i, reference)| {
            if reference.is_empty() {
                return Ok(None);
            }
            let fetched;
            let alignments: &[(usize, Alignment)] = match &prefetched {
                Some(all) => &all[i],
                None => {
                    fetched = aligner.against_reference_nested(reference, seqs, Some(i))?;
                    &fetched
                }
            };
            let plain: i64 = alignments.iter().map(|(_, a)| a.score as i64).sum();
            let score: i64 = if params.nonredundant_reference_score {
                nonredundant_candidate_score(alignments)
            } else {
                plain
            };
            // Why a redundancy correction may leave the argmax untouched: if it
            // discounts every candidate by a similar proportion, the ranking is
            // unchanged and only the margin moves.
            if std::env::var_os("TE_COMPOSER_CAND_DEBUG").is_some() {
                let nr = nonredundant_candidate_score(alignments);
                eprintln!(
                    "CAND\t{}\t{}\t{}\t{:.4}",
                    i,
                    plain,
                    nr,
                    if plain != 0 { nr as f64 / plain as f64 } else { 1.0 }
                );
            }
            // The C++ only calls a consensus when the total score is positive.
            if score <= 0 {
                return Ok(None);
            }
            Ok(Some((i, score)))
        })
        .collect();

    let mut ranked: Vec<(usize, i64)> = Vec::with_capacity(seqs.len());
    for r in scored {
        if let Some(c) = r? {
            ranked.push(c);
        }
    }
    ranked.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

    // Pass two: a consensus for the ones that will be used.
    let build = |i: usize| -> Result<Vec<u8>> {
        let fetched;
        let alignments: &[(usize, Alignment)] = match &prefetched {
            Some(all) => &all[i],
            None => {
                fetched = aligner.against_reference_nested(&seqs[i], seqs, Some(i))?;
                &fetched
            }
        };
        let msa = build_msa(&seqs[i], seqs, alignments, params.insertions)?;
        let gapped = match params.caller {
            Caller::Dfam => msa.consensus(&phase1),
            // Phase 1's C++ call site takes minNonGapCount's default of 1,
            // and keeps the reference row.
            Caller::Giri => {
                let rows: Vec<&[u8]> =
                    msa.sequences.iter().map(|r| r.seq.as_slice()).collect();
                aln_core::giri::get_consensus(&rows, 1)
            }
        };
        Ok(seqmod::ungap(&gapped))
    };

    let mut out = Vec::with_capacity(ranked.len());
    for (n, (i, score)) in ranked.into_iter().enumerate() {
        let consensus = if n < consensus_for { build(i)? } else { Vec::new() };
        out.push(Candidate { index: i, score, consensus });
    }
    Ok(out)
}

/// The input sequence a consensus was seeded from.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Reference {
    /// 0-based index into the input set.
    pub index: usize,
    /// That sequence's name.
    pub name: String,
    /// Its phase-1 score: the sum over its alignments to every other input.
    pub score: i64,
    /// The next-best candidate's score, when there was one.
    ///
    /// Reference selection is a discrete argmax, so the *margin* is what says
    /// whether the choice was forced or a coin toss — a run that won by 30
    /// points out of 200,000 could have gone either way, and every later
    /// stage inherits that.
    pub runner_up: Option<i64>,
    /// How many candidates were scored.
    pub candidates: usize,
}

/// Mean Kimura divergence of the instance rows against a gapped consensus,
/// plain and CpG-adjusted.
///
/// CpG dinucleotides mutate fast enough to dominate a raw estimate, so the
/// adjusted figure is the one families compare on; the pair is returned
/// because the gap between them tracks CpG content.
pub fn mean_kimura(msa: &MultiAlign, gapped_cons: &[u8]) -> (f64, f64) {
    use aln_core::stats::{kimura_divergence, Masking};
    let (mut plain, mut cpg) = (Vec::new(), Vec::new());
    for row in msa.sequences.iter().skip(1) {
        if row.seq.len() != gapped_cons.len() {
            continue;
        }
        for (adj, sink) in [(false, &mut plain), (true, &mut cpg)] {
            if let Ok(d) = kimura_divergence(&row.seq, gapped_cons, adj, Masking::Ignore) {
                if let Some(v) = d.value {
                    sink.push(v);
                }
            }
        }
    }
    let mean = |v: &Vec<f64>| if v.is_empty() { 0.0 } else { v.iter().sum::<f64>() / v.len() as f64 };
    (mean(&plain), mean(&cpg))
}

/// Phase 2: refine one consensus until it stops changing, cycles, or runs out.
///
/// Runs at most `params.iterations + 1` passes. Each pass re-aligns the whole
/// input to the current consensus, drops the reference row, and re-calls.
///
/// # Stopping
///
/// The C++ stops only on a fixed point (`next == current`). That is a cycle of
/// length 1; longer cycles exist — A produces B, B produces A — and the C++
/// spends its whole budget oscillating between them. Every consensus seen is
/// therefore remembered, and a repeat of *any* of them stops the loop
/// ([`StopReason::Cycled`]). At `--iterations 3` that saves little; at 11 it is
/// most of the budget.
///
/// # Which pass is returned
///
/// Always the last one, whatever the stop reason.
///
/// Selecting the *best-scoring* pass instead was tried and measured worse.
/// Score is self-referential — it says how well a consensus explains the
/// instances, not how close it is to the progenitor — and a consensus can
/// explain the data well while drifting from the truth. Over ten TFE sims with
/// a known root, taking the last pass scored oracle-optimal 8/10 against
/// selection-by-score's 7/10 (mean normalised accuracy 0.3857 vs 0.3848); the
/// worst case was a cycling family where the score peaked at pass 4 while the
/// truth-optimal consensus was at pass 15, so score selection returned an
/// answer 12% worse than simply taking the last.
///
/// Taking the last pass also keeps this loop cheap and obvious: nothing is held
/// across iterations except the set of consensi seen, and the MSA in hand at
/// the end is by construction the one that produced the returned consensus.
pub fn refine<A: AlignmentSource>(
    aligner: &A,
    seqs: &[Sequence],
    start: &[u8],
    name: &str,
    reference: Reference,
    params: &Params,
) -> Result<Option<Refined>> {
    /// One pass: align everything to `input`, build the MSA, call the next
    /// consensus. The unit that `refine` both iterates and replays.
    fn pass<A: AlignmentSource>(
        aligner: &A,
        seqs: &[Sequence],
        input: &Sequence,
        params: &Params,
    ) -> Result<(i64, f64, MultiAlign, Vec<u8>, Vec<u8>, usize)> {
        let alignments = aligner.against_reference(input, seqs, None)?;
        let score: i64 = alignments.iter().map(|(_, a)| a.score as i64).sum();
        // Distinct contributors, taken from the alignment indices rather than
        // from the row count: `--hsps tiled` gives one instance several rows.
        let instances = alignments
            .iter()
            .map(|(i, _)| *i)
            .collect::<std::collections::HashSet<_>>()
            .len();
        let norm = mean_score_per_base(&alignments);
        let msa = build_msa(input, seqs, &alignments, params.insertions)?;
        // Phase 2 drops the reference row before calling, and here `--min` does
        // apply.
        let gapped = gapped_consensus(&msa, params);
        let next = seqmod::ungap(&gapped);
        Ok((score, norm, msa, gapped, next, instances))
    }

    /// Recover inherited insertions the column caller cannot see, once the
    /// column-wise loop has settled. Returns the packed consensus.
    ///
    /// Packing used to run inside every pass. It re-derived the same spans
    /// pass after pass and dominated runtime on deep families (13x on 100
    /// copies x 12 kb at 20% divergence). Packing at the settled point and
    /// then re-refining keeps what the per-pass version bought — a recovered
    /// base joins the reference and the copies re-align to it — at one packing
    /// per round instead of one per pass.
    fn pack_settled(msa: &MultiAlign, gapped: &[u8], params: &Params) -> Option<Vec<u8>> {
        let mx = params.repair_matrix.as_ref()?;
        let call = ConsensusParams { include_reference: false, ..params.consensus.clone() };
        let packed = crate::lowqual::pack_insertion_spans(
            msa,
            gapped,
            mx,
            &call,
            params.pack_max_sep,
            params.pack_min_seg,
            params.pack_min_score,
            params.pack_min_occupancy,
            params.pack_keep_insertions,
        );
        Some(seqmod::ungap(&packed))
    }

    /// Packing rounds per refinement: pack, re-refine to a fixed point, pack
    /// again. The second round is for insertions that only become registrable
    /// once the first round's bases are in the reference.
    const MAX_PACK_ROUNDS: usize = 2;

    /// `acons --mam`. The C++ restores CpG once, on the final pass only,
    /// against the *gapped* consensus and the instance rows (its `maln` has
    /// already dropped the reference), then ungaps — autocons.cpp:730.
    /// Convergence is decided on the unrestored consensus, exactly as there:
    /// restoration must not feed back into the loop.
    fn finish(gapped: &[u8], plain: Vec<u8>, msa: &MultiAlign, params: &Params) -> Vec<u8> {
        if params.restore_cpg && params.caller == Caller::Giri {
            let mut g = gapped.to_vec();
            let rows: Vec<&[u8]> = msa.sequences[1.min(msa.sequences.len())..]
                .iter()
                .map(|r| r.seq.as_slice())
                .collect();
            giri::restore_cpg(&mut g, &rows);
            seqmod::ungap(&g)
        } else {
            plain
        }
    }

    let mut current = Sequence::new(name, start.to_vec());
    // Every consensus that has been *fed to* a pass. A repeat means a cycle.
    let mut seen: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();
    let mut trace: Vec<PassTrace> = Vec::new();

    // `iterations` is a pass count, not a ceiling on an index: `--iterations 10`
    // runs at most ten passes. It was `0..=` with an `== params.iterations`
    // exhaustion test, and reporting compensated by printing `budget + 1`, so
    // the flag ran one pass more than it said. Both had to move together: with
    // an exclusive range the old test never fires and the loop falls through.
    let budget = params.iterations.max(1);
    let mut pass_idx = 0usize;
    // Passes since the last packing round. The budget is per round, so a
    // packed reference gets a full loop to settle.
    let mut since_pack = 0usize;
    let mut pack_rounds = 0usize;
    loop {
        if current.is_empty() {
            return Ok(None);
        }
        seen.insert(current.seq.clone());

        let (score, norm_score, msa, gapped, next, instances) =
            pass(aligner, seqs, &current, params)?;
        let (divergence, divergence_cpg) = mean_kimura(&msa, &gapped);
        trace.push(PassTrace {
            pass: pass_idx + 1,
            score,
            norm_score,
            consensus: next.clone(),
            divergence,
            divergence_cpg,
            rows: msa.num_instances(),
            instances,
        });

        let converged = next == current.seq;
        // A fixed point is the 1-cycle; `seen` catches the longer ones.
        let cycled = !converged && seen.contains(&next);
        let exhausted = since_pack + 1 >= budget;

        if converged || cycled || exhausted {
            // The column caller has settled; now recover what it cannot see.
            // A packed consensus that the next loop does not reproduce is
            // dropped by that loop, so what is emitted is always a fixed point
            // of the column caller.
            if params.pack_insertions && pack_rounds < MAX_PACK_ROUNDS {
                if let Some(packed) = pack_settled(&msa, &gapped, params) {
                    if packed != next && !seen.contains(&packed) {
                        pack_rounds += 1;
                        since_pack = 0;
                        pass_idx += 1;
                        current = Sequence::new(name, packed);
                        continue;
                    }
                }
            }
            let stop = if converged {
                StopReason::Converged
            } else if cycled {
                StopReason::Cycled
            } else {
                StopReason::Exhausted
            };
            return Ok(Some(Refined {
                name: name.to_string(),
                consensus: finish(&gapped, next, &msa, params),
                reference,
                score,
                norm_score,
                msa,
                passes: pass_idx + 1,
                converged,
                stop,
                trace,
                notes: Vec::new(),
            }));
        }
        current = Sequence::new(name, next);
        pass_idx += 1;
        since_pack += 1;
    }
}

/// Refine, repair low-quality blocks, then refine again.
///
/// `Refiner`'s full phase-2: `refineUntil`, `resolveLowQualityBlocks`,
/// `refineUntil`. The repair works on the *gapped* consensus, so it must run
/// against the MSA that produced it — the alignment's columns are the
/// coordinate system the block ranges are expressed in.
///
/// Falls back to the single-round result whenever the repair finds nothing,
/// and — importantly — whenever the second round comes back *worse*. The
/// repair is a heuristic re-derivation over a handful of columns; it can be
/// wrong, and a workflow step that can only help is one nobody has tested
/// against a case where it hurts.
pub fn refine_with_repair<A: AlignmentSource>(
    aligner: &A,
    seqs: &[Sequence],
    start: &[u8],
    name: &str,
    reference: Reference,
    params: &Params,
    matrix: &aln_core::SubstMatrix,
) -> Result<Option<Refined>> {
    refine_with_repair_using(
        aligner,
        seqs,
        start,
        name,
        reference,
        params,
        |msa, call| default_block_source(msa, matrix, params, call),
    )
}

/// The blocks `refine_with_repair` repairs by default: `Refiner`'s low-quality
/// selection, plus the scanning window when `repair_window` is set.
///
/// Exposed so an experiment can ask what a *different* selection would have
/// produced without reimplementing the rest of the workflow around it.
pub fn default_block_source(
    msa: &MultiAlign,
    matrix: &aln_core::SubstMatrix,
    params: &Params,
    call: &ConsensusParams,
) -> Vec<crate::lowqual::BlockFix> {
    // `Refiner`'s two-branch repair, both branches. Where a single instance
    // length dominates, the length vote resolves the block; where none does,
    // the block is re-derived from an all-against-all alignment of its
    // instances. That second branch used to be skipped, which left every
    // no-majority block unrepaired — a gap against `Refiner` that widens
    // exactly where copies are numerous and length-diverse.
    let (go, ge) = crate::lowqual::scaled_gap_penalties(matrix);
    let fixes = crate::lowqual::resolve_low_quality_blocks(
        msa,
        matrix,
        params.repair_threshold,
        call,
        params.repair_all_vs_all.then_some(|inst: &[Vec<u8>]| {
            crate::lowqual::resolve_by_all_vs_all(
                inst,
                matrix,
                go.abs().round() as u32,
                ge.abs().round() as u32,
                call,
            )
        }),
    );
    if params.repair_window == 0 {
        return fixes;
    }
    // Second block source: the scanning window. Same vote, different way of
    // deciding where to look; overlaps with the low-quality fixes are dropped.
    let vp = crate::lowqual::VoteParams::with_window(params.repair_window, 5);
    let windows: Vec<(usize, usize)> = crate::lowqual::length_vote_blocks(msa, &vp)
        .into_iter()
        .map(|v| (v.col_start, v.col_end))
        .collect();
    let mut all = fixes;
    all.extend(crate::lowqual::resolve_given_blocks(
        msa,
        &windows,
        crate::lowqual::MAX_BLOCK,
        call,
        &all,
    ));
    all
}

/// `refine_with_repair` with the block selection supplied by the caller.
///
/// The workflow around the repair — patch the *gapped* consensus, refine again
/// to a fixed point, keep the result only if it beats the unrepaired one — is
/// the part that must not vary between an experiment and the pipeline, so it
/// lives here once and every selection is measured through it.
pub fn refine_with_repair_using<A, F>(
    aligner: &A,
    seqs: &[Sequence],
    start: &[u8],
    name: &str,
    reference: Reference,
    params: &Params,
    block_source: F,
) -> Result<Option<Refined>>
where
    A: AlignmentSource,
    F: FnOnce(&MultiAlign, &ConsensusParams) -> Vec<crate::lowqual::BlockFix>,
{
    let first = refine(aligner, seqs, start, name, reference.clone(), params)?;
    let Some(first) = first else { return Ok(None) };
    if !params.repair_blocks {
        return Ok(Some(first));
    }

    // The consensus in `first` is ungapped; the repair needs it in the MSA's
    // column coordinates, so re-derive the gapped form from the same MSA.
    let call = ConsensusParams { include_reference: false, ..params.consensus.clone() };
    let gapped = gapped_consensus(&first.msa, params);

    let fixes = block_source(&first.msa, &call);
    if std::env::var_os("TE_COMPOSER_REPAIR_DEBUG").is_some() {
        eprintln!(
            "repair: width={} rows={} cons={} fixes={}",
            first.msa.width(),
            first.msa.sequences.len(),
            first.consensus.len(),
            fixes.len()
        );
    }
    if fixes.is_empty() {
        return Ok(Some(first));
    }

    apply_repair(aligner, seqs, first, &gapped, name, reference, params, &fixes)
        .map(|o| Some(o.take(params.repair_accept)))
}

/// Both statistics and both verdicts for one judged repair.
#[derive(Debug, Clone, Copy)]
pub struct RepairReport {
    pub fixes: usize,
    pub first_score: i64,
    pub second_score: i64,
    pub first_norm: f64,
    pub second_norm: f64,
    pub keep_total: bool,
    pub keep_mean: bool,
}

/// What the repair step decided.
pub enum RepairOutcome {
    /// Nothing was proposed, or the patch changed nothing.
    NoFixes(Refined),
    /// A repair was proposed and refined; `report` carries both verdicts so a
    /// caller can ask what a rule other than the configured one would have
    /// done, without paying for a second refinement.
    Judged {
        first: Refined,
        second: Refined,
        report: RepairReport,
    },
}

impl RepairOutcome {
    /// The consensus a given accept rule keeps.
    pub fn take(self, rule: AcceptRule) -> Refined {
        match self {
            RepairOutcome::NoFixes(first) => first,
            RepairOutcome::Judged { first, second, report } => {
                let keep = match rule {
                    AcceptRule::Total => report.keep_total,
                    AcceptRule::MeanPerBase => report.keep_mean,
                    AcceptRule::Always => true,
                };
                if keep {
                    second
                } else {
                    first
                }
            }
        }
    }

    pub fn report(&self) -> Option<RepairReport> {
        match self {
            RepairOutcome::NoFixes(_) => None,
            RepairOutcome::Judged { report, .. } => Some(*report),
        }
    }
}

/// Patch the gapped consensus, refine again to a fixed point, and judge.
///
/// The half of the repair workflow that must be identical between the pipeline
/// and any experiment measuring a block selection: an experiment that splices
/// into the *ungapped* consensus and stops measures the patch, not the method,
/// and the two answers differ in sign.
#[allow(clippy::too_many_arguments)]
pub fn apply_repair<A: AlignmentSource>(
    aligner: &A,
    seqs: &[Sequence],
    first: Refined,
    gapped: &[u8],
    name: &str,
    reference: Reference,
    params: &Params,
    fixes: &[crate::lowqual::BlockFix],
) -> Result<RepairOutcome> {
    if fixes.is_empty() {
        return Ok(RepairOutcome::NoFixes(first));
    }
    let patched = seqmod::ungap(&crate::lowqual::patch_gapped_consensus(gapped, fixes));
    if patched.is_empty() || patched == first.consensus {
        return Ok(RepairOutcome::NoFixes(first));
    }

    let Some(second) = refine(aligner, seqs, &patched, name, reference, params)? else {
        return Ok(RepairOutcome::NoFixes(first));
    };
    let report = RepairReport {
        fixes: fixes.len(),
        first_score: first.score,
        second_score: second.score,
        first_norm: first.norm_score,
        second_norm: second.norm_score,
        keep_total: second.score > first.score,
        keep_mean: second.norm_score > first.norm_score,
    };
    if std::env::var_os("TE_COMPOSER_GATE_LOG").is_some() {
        eprintln!(
            "GATE\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}",
            name, report.fixes,
            report.first_score, report.second_score,
            report.first_norm, report.second_norm,
            report.keep_total as u8, report.keep_mean as u8,
            first.consensus.len(), second.consensus.len(),
        );
    }
    Ok(RepairOutcome::Judged { first, second, report })
}

/// Resolve `--num` into a count.
///
/// Below 1.0 it is a fraction of the input size, floored, with a minimum of 1 —
/// the C++'s `nCon = int(nCon * seqls.size())` followed by a floor of 1.
pub fn resolve_count(num: f64, total: usize) -> usize {
    if num < 1.0 {
        ((num * total as f64) as usize).max(1)
    } else {
        (num as usize).max(1)
    }
}

/// Run both phases.
pub fn run<A: AlignmentSource>(
    aligner: &A,
    seqs: &[Sequence],
    params: &Params,
) -> Result<Vec<Refined>> {
    let want_raw = resolve_count(params.num_consensi, seqs.len());
    let candidates = score_candidates_ranked(aligner, seqs, params, want_raw)?;
    let want = want_raw.min(candidates.len());
    let multiple = want > 1;

    let mut out = Vec::with_capacity(want);
    for (n, cand) in candidates.iter().take(want).enumerate() {
        let name = if multiple {
            format!("{}{}", params.base_name, n + 1)
        } else {
            params.base_name.clone()
        };
        let reference = Reference {
            index: cand.index,
            name: seqs[cand.index].name.clone(),
            score: cand.score,
            // Candidates are sorted best-first, so the next entry is the
            // margin this choice won by.
            runner_up: candidates.get(n + 1).map(|c| c.score),
            candidates: candidates.len(),
        };
        let refined = match (&params.repair_matrix, params.repair_blocks) {
            (Some(mx), true) => refine_with_repair(
                aligner, seqs, &cand.consensus, &name, reference, params, mx,
            )?,
            _ => refine(aligner, seqs, &cand.consensus, &name, reference, params)?,
        };
        if let Some(r) = refined {
            out.push(r);
        }
    }
    Ok(out)
}

/// Counts calls so a reader can tell one refinement pass from the next.
static CENSUS_CALL: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// Record every alignment that passes through, then hand it on unchanged.
///
/// A census of what actually enters the MSA. GIRI applies no score cutoff
/// (`MultipleAlignment::align` hardwires `minScore = 0`), so in principle a
/// two-base alignment can become a row; this is how to find out whether that
/// happens in practice.
///
/// Rows are `phase, ref_index, query_index, score, columns, query_span,
/// subj_span`. `phase` is 1 for reference selection (all-against-all) and 2 for
/// refinement against the working consensus.
pub struct Census<S> {
    inner: S,
    out: std::sync::Mutex<Box<dyn std::io::Write + Send>>,
}

impl<S> Census<S> {
    /// `sink` receives a TSV header immediately.
    pub fn new(inner: S, mut sink: Box<dyn std::io::Write + Send>) -> std::io::Result<Self> {
        writeln!(sink, "call\tphase\tref\tquery\tscore\tcolumns\tqspan\tsspan\tqstart\tqend\tstrand")?;
        Ok(Census { inner, out: std::sync::Mutex::new(sink) })
    }

    fn record(&self, phase: u8, refi: usize, hits: &[(usize, Alignment)]) {
        // One counter per call, so a reader can isolate a single refinement
        // pass. Without it every pass pools together and one alignment
        // repeated across twenty passes reads as twenty rows of the same copy
        // re-using the same bases — which is exactly the question this file
        // exists to answer, so pooling makes it unanswerable.
        let call = CENSUS_CALL.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let mut w = self.out.lock().unwrap_or_else(|e| e.into_inner());
        for (qi, a) in hits {
            let _ = writeln!(
                w,
                "{call}\t{phase}\t{refi}\t{qi}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                a.score,
                a.edits.align_len(),
                a.query_end.saturating_sub(a.query_start),
                a.subj_end.saturating_sub(a.subj_start),
                // Coordinates, not just extents: whether two rows of the same
                // copy re-use the same genomic bases is only answerable from
                // the ranges, and that is what a profile HMM built from this
                // alignment would double-count.
                a.query_start,
                a.query_end,
                if a.strand.is_minus() { '-' } else { '+' },
            );
        }
    }
}

impl<S: AlignmentSource> AlignmentSource for Census<S> {
    fn against_reference(
        &self,
        reference: &Sequence,
        seqs: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>> {
        let hits = self.inner.against_reference(reference, seqs, skip)?;
        // Phase 2 aligns against the working consensus, which has no index.
        self.record(2, usize::MAX, &hits);
        Ok(hits)
    }

    fn against_every_reference(
        &self,
        seqs: &[Sequence],
    ) -> Result<Vec<Vec<(usize, Alignment)>>> {
        let groups = self.inner.against_every_reference(seqs)?;
        for (i, g) in groups.iter().enumerate() {
            self.record(1, i, g);
        }
        Ok(groups)
    }

    fn against_reference_nested(
        &self,
        reference: &Sequence,
        seqs: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>> {
        let hits = self.inner.against_reference_nested(reference, seqs, skip)?;
        self.record(1, skip.unwrap_or(usize::MAX), &hits);
        Ok(hits)
    }

    fn batches_all_vs_all(&self) -> bool {
        self.inner.batches_all_vs_all()
    }
}

/// Drop alignments scoring below a threshold before they reach the MSA.
///
/// The exact analogue of the C++ `autocons --minscore`, which sets `minScore`
/// in `MultipleAlignment::align` (hardwired to 0 until now). This is a
/// **post-alignment** filter: the aligner does whatever it does, and this
/// decides what is allowed to become a row.
///
/// Distinct from [`AlignParams::min_score`], which is handed to the aligner
/// itself. For a DP backend the two coincide, but for a seeded backend
/// `min_score` also sets the X-drop budget and so changes what is *found* —
/// see `--xdrop`. Filtering here is backend-independent, which is what a
/// like-for-like comparison needs.
///
/// 0 keeps everything, matching GIRI.
pub struct MinScore<S> {
    inner: S,
    min: i32,
}

impl<S> MinScore<S> {
    pub fn new(inner: S, min: i32) -> Self {
        MinScore { inner, min }
    }
}

impl<S: AlignmentSource> MinScore<S> {
    fn keep(&self, hits: Vec<(usize, Alignment)>) -> Vec<(usize, Alignment)> {
        if self.min <= 0 {
            return hits;
        }
        hits.into_iter().filter(|(_, a)| a.score >= self.min).collect()
    }
}

impl<S: AlignmentSource> AlignmentSource for MinScore<S> {
    fn against_reference(
        &self,
        reference: &Sequence,
        seqs: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>> {
        Ok(self.keep(self.inner.against_reference(reference, seqs, skip)?))
    }

    fn against_every_reference(
        &self,
        seqs: &[Sequence],
    ) -> Result<Vec<Vec<(usize, Alignment)>>> {
        Ok(self
            .inner
            .against_every_reference(seqs)?
            .into_iter()
            .map(|g| self.keep(g))
            .collect())
    }

    fn against_reference_nested(
        &self,
        reference: &Sequence,
        seqs: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>> {
        Ok(self.keep(self.inner.against_reference_nested(reference, seqs, skip)?))
    }

    fn batches_all_vs_all(&self) -> bool {
        self.inner.batches_all_vs_all()
    }
}

/// Collapse a batched search to **one alignment per instance**.
///
/// GIRI's `MultipleAlignment::align` takes a single `PairwiseAlignment` per
/// sequence from its aligner, so an MSA row per instance is what the original
/// workflow assumes. Wrapping the search engine in this keeps the batching —
/// and its speed — while restoring those semantics, which is what an
/// aligner-only comparison needs.
///
/// [`AlignmentSource`] on the bare engine keeps every HSP instead; see
/// `--hsps`.
pub struct BestHsp<S>(pub S);

impl<S: AlignmentSource> AlignmentSource for BestHsp<S> {
    fn against_reference(
        &self,
        reference: &Sequence,
        seqs: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>> {
        Ok(best_per_instance(self.0.against_reference(reference, seqs, skip)?))
    }

    fn against_every_reference(
        &self,
        seqs: &[Sequence],
    ) -> Result<Vec<Vec<(usize, Alignment)>>> {
        Ok(self
            .0
            .against_every_reference(seqs)?
            .into_iter()
            .map(best_per_instance)
            .collect())
    }

    fn against_reference_nested(
        &self,
        reference: &Sequence,
        seqs: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>> {
        Ok(best_per_instance(
            self.0.against_reference_nested(reference, seqs, skip)?,
        ))
    }

    fn batches_all_vs_all(&self) -> bool {
        self.0.batches_all_vs_all()
    }
}

/// `Refiner`'s post-alignment overlap filter, applied to both phases.
///
/// Phase 1 ports the *enabled* branch of RepeatModeler `Refiner`'s
/// `findHighestScoringAlignmentSet`: within each (instance, strand) group,
/// alignments are ranked by score, and a lower-scoring alignment survives only
/// if at least `100 - mask_level` percent of its **reference-axis** span
/// extends beyond a higher-scoring one — dedup before candidate sums. Phase 2
/// mirrors Refiner's one-vs-all engine mask instead, measuring on the
/// **instance** axis; see [`MaskAxis`] for why the phases must differ.
///
/// Three properties are faithful to the Perl rather than accidents:
///
/// * **Single-axis.** Novelty is measured on the reference span only. Two
///   alignments mapping the *same instance region* to different reference
///   regions both survive — the tandem-repeat register case. Refiner carries a
///   double-axis variant in a disabled `else` branch; this port deliberately
///   matches the enabled behaviour (decision 2026-08-14: leave double-axis out).
/// * **Deleted alignments still mask.** The Perl iterates higher-scoring
///   entries whether or not they were themselves deleted, so a chain of
///   mutually-overlapping alignments collapses to its best member, not to an
///   alternating survivors pattern.
/// * **Deterministic ordering.** Ties break by reference start then longer
///   span, mirroring Refiner's 9/19/23 determinism fix.
///
/// Applying this in phase 1 is the substantive change: reference selection sums
/// scores per candidate, and unfiltered HSP sets double-count overlapping
/// registers, which can elevate the wrong candidate (observed: a truncated
/// SVA_A instance winning phase 1, costing 300 bp of consensus). `Refiner`
/// filters before ranking; so does this.
pub struct RefinerFilter<S> {
    inner: S,
    mask_level: u32,
    /// Let alignments of one instance compete regardless of orientation.
    cross_strand: bool,
    /// Keep only the better-scoring orientation per instance, in phase 1
    /// (reference-axis culling).
    best_orientation_phase1: bool,
    /// The same in phase 2 (instance-axis culling), kept separate so phase-2
    /// culling can be varied without disturbing a settled phase 1.
    best_orientation_phase2: bool,
    /// Do anything at all in phase 2.
    ///
    /// Off lets a different phase-2 selection sit on top and see raw HSPs,
    /// which is how the tile-then-vote redesign is measured in phase 2 alone
    /// while phase 1 keeps the shipping behaviour.
    phase2_enabled: bool,
}

impl<S> RefinerFilter<S> {
    /// Leave phase 2 untouched, for a caller that replaces it.
    pub fn without_phase2(mut self) -> Self {
        self.phase2_enabled = false;
        self
    }
}

impl<S> RefinerFilter<S> {
    /// `mask_level >= 100` disables filtering entirely.
    pub fn new(inner: S, mask_level: u32) -> Self {
        RefinerFilter {
            inner,
            mask_level,
            cross_strand: false,
            best_orientation_phase1: true,
            best_orientation_phase2: true,
            phase2_enabled: true,
        }
    }

    /// As [`new`](Self::new), with orientation ignored when masking — the
    /// masklevel semantics RepeatMasker's own filter implements.
    pub fn with_cross_strand(inner: S, mask_level: u32, cross_strand: bool) -> Self {
        RefinerFilter {
            inner,
            mask_level,
            cross_strand,
            best_orientation_phase1: true,
            best_orientation_phase2: true,
            phase2_enabled: true,
        }
    }

    /// Full control: overlap masking, plus the orientation rule.
    ///
    /// `best_orientation` is not a port of anything — `Refiner`'s
    /// `findHighestScoringAlignmentSet` masks *within* each strand and never
    /// compares the two sets, so a copy can contribute forward and reverse rows
    /// to the same alignment as long as they do not pairwise overlap. For
    /// consensus building there is no case where fragments of one copy in
    /// opposite orientations are both right: the copy inserted in one
    /// orientation. This sums the surviving score per orientation and keeps the
    /// better set, preferring forward on a tie.
    pub fn with_options(
        inner: S,
        mask_level: u32,
        cross_strand: bool,
        best_orientation: bool,
    ) -> Self {
        Self::with_phase_options(inner, mask_level, cross_strand, best_orientation, best_orientation)
    }

    /// As [`with_options`](Self::with_options), with the orientation rule set
    /// per phase.
    ///
    /// Phase 1 is settled — orientation plus non-redundant candidate scoring —
    /// and phase-2 culling is still open, so the two must be independently
    /// controllable or an experiment on the second silently moves the first.
    pub fn with_phase_options(
        inner: S,
        mask_level: u32,
        cross_strand: bool,
        best_orientation_phase1: bool,
        best_orientation_phase2: bool,
    ) -> Self {
        RefinerFilter {
            inner,
            mask_level,
            cross_strand,
            best_orientation_phase1,
            best_orientation_phase2,
            phase2_enabled: true,
        }
    }
}

/// As [`refiner_keep`], with the option to ignore orientation.
///
/// `cross_strand` collapses the two strand groups into one, so a reverse-strand
/// alignment competes with the forward alignments of the same instance. That is
/// what masklevel is supposed to mean — `SearchResultCollection::maskLevelFilter`,
/// RepeatMasker's own implementation, groups by sequence name and never looks at
/// orientation, and `RMBlast_project/BugInMaskLevel/notebook` states the intent
/// outright: "Crossmatch doesn't care if the overlap occurs on one strand or the
/// other and neither should rmblast."
///
/// Grouping by strand lets a spurious reverse-strand HSP of a copy region enter
/// the alignment beside its forward twin. Measured on 180 simulated families,
/// 43.9% of same-copy row pairs were near-duplicates (99.6% reference, 91.9%
/// instance overlap) and 72.9% of those survived only because they were on
/// opposite strands.
/// [`refiner_keep_opt`] with cross-strand comparison off, as the phases use it.
#[cfg(test)]
fn refiner_keep(items: &[(i64, usize, usize, Strand)], mask_level: u32) -> Vec<bool> {
    refiner_keep_opt(items, mask_level, false)
}

fn refiner_keep_opt(
    items: &[(i64, usize, usize, Strand)],
    mask_level: u32,
    cross_strand: bool,
) -> Vec<bool> {
    let mut keep = vec![true; items.len()];
    if mask_level >= 100 || items.len() < 2 {
        return keep;
    }
    // Group by strand; within a group, rank by (score desc, start asc, len desc).
    // Under `cross_strand` there is one group holding everything.
    let groups: &[&[Strand]] = if cross_strand {
        &[&[Strand::Plus, Strand::Minus]]
    } else {
        &[&[Strand::Plus], &[Strand::Minus]]
    };
    for strands in groups {
        let mut idx: Vec<usize> =
            (0..items.len()).filter(|&i| strands.contains(&items[i].3)).collect();
        idx.sort_by(|&a, &b| {
            let (sa, ba, ea, _) = items[a];
            let (sb, bb, eb, _) = items[b];
            sb.cmp(&sa)
                .then(ba.cmp(&bb))
                .then((eb - bb).cmp(&(ea - ba)))
        });
        for i in 0..idx.len() {
            for j in (i + 1)..idx.len() {
                let (_, b1, e1, _) = items[idx[i]];
                let (_, b2, e2, _) = items[idx[j]];
                // Disjoint spans never mask (the Perl's `next if` guard).
                if b2 >= e1 || b1 >= e2 {
                    continue;
                }
                // The lower-scoring alignment's overhang outside the higher-
                // scoring one: left + right, as the Perl computes it.
                let novel = b1.saturating_sub(b2) + e2.saturating_sub(e1);
                let perc = 100.0 * novel as f64 / (e2 - b2) as f64;
                if perc < (100 - mask_level) as f64 {
                    keep[idx[j]] = false;
                }
            }
        }
    }
    keep
}

/// Which coordinate axis novelty is measured on.
///
/// `Refiner` is not symmetric between phases, and the asymmetry is
/// load-bearing:
///
/// * **Phase 1** (`findHighestScoringAlignmentSet`) measures on the candidate
///   **reference** axis — dedup registers before summing candidate scores.
/// * **Phase 2** (its one-vs-all engine mask, `setQuery(instFile)` +
///   mask level 80) measures on the **instance** axis.
///
/// Applying the reference axis to both phases was tried and measurably wrong:
/// in a dimeric element (AluY), an instance's left arm also aligns to the
/// consensus's *right* arm — novel on the reference axis, so a ref-axis
/// phase-2 filter keeps it, and the MSA fills with arm-swapped rows
/// (substitutions 4.2 -> 47.0 bp/kb on t1/t5 AluY). The instance axis kills
/// exactly those: the duplicate covers the same instance span as the better
/// alignment.
#[derive(Clone, Copy, PartialEq, Eq)]
enum MaskAxis {
    Reference,
    Instance,
}

/// Apply [`refiner_keep`] to one candidate's hit list, grouping by instance.
fn refiner_mask(
    hits: Vec<(usize, Alignment)>,
    mask_level: u32,
    axis: MaskAxis,
    cross_strand: bool,
    best_orientation: bool,
) -> Vec<(usize, Alignment)> {
    if (mask_level >= 100 && !best_orientation) || hits.len() < 2 {
        return hits;
    }
    let mut by: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();
    for (pos, (inst, _)) in hits.iter().enumerate() {
        by.entry(*inst).or_default().push(pos);
    }
    let mut keep = vec![true; hits.len()];
    for positions in by.values() {
        if positions.len() < 2 {
            continue;
        }
        let items: Vec<(i64, usize, usize, Strand)> = positions
            .iter()
            .map(|&p| {
                let a = &hits[p].1;
                match axis {
                    MaskAxis::Reference => (a.score as i64, a.subj_start, a.subj_end, a.strand),
                    MaskAxis::Instance => (a.score as i64, a.query_start, a.query_end, a.strand),
                }
            })
            .collect();
        for (k, &p) in refiner_keep_opt(&items, mask_level, cross_strand)
            .iter()
            .zip(positions.iter())
        {
            keep[p] = *k;
        }

        // The orientation decision, taken on what survived masking: a copy
        // inserted in one orientation, so only one set of fragments can be
        // describing it.
        if best_orientation {
            let (mut fwd, mut rev) = (0i64, 0i64);
            for &p in positions.iter().filter(|&&p| keep[p]) {
                let a = &hits[p].1;
                match a.strand {
                    Strand::Plus => fwd += a.score as i64,
                    Strand::Minus => rev += a.score as i64,
                }
            }
            if fwd > 0 || rev > 0 {
                // Forward wins ties.
                let losing = if rev > fwd { Strand::Plus } else { Strand::Minus };
                for &p in positions.iter() {
                    if hits[p].1.strand == losing {
                        keep[p] = false;
                    }
                }
            }
        }
    }
    hits.into_iter()
        .zip(keep)
        .filter_map(|(h, k)| k.then_some(h))
        .collect()
}

impl<S: AlignmentSource> AlignmentSource for RefinerFilter<S> {
    fn against_reference(
        &self,
        reference: &Sequence,
        seqs: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>> {
        let hits = self.inner.against_reference(reference, seqs, skip)?;
        if !self.phase2_enabled {
            return Ok(hits);
        }
        Ok(refiner_mask(
            hits,
            self.mask_level,
            MaskAxis::Instance,
            self.cross_strand,
            self.best_orientation_phase2,
        ))
    }

    fn against_every_reference(
        &self,
        seqs: &[Sequence],
    ) -> Result<Vec<Vec<(usize, Alignment)>>> {
        Ok(self
            .inner
            .against_every_reference(seqs)?
            .into_iter()
            .map(|hits| {
                refiner_mask(
                    hits,
                    self.mask_level,
                    MaskAxis::Reference,
                    self.cross_strand,
                    self.best_orientation_phase1,
                )
            })
            .collect())
    }

    fn against_reference_nested(
        &self,
        reference: &Sequence,
        seqs: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>> {
        // Phase 1, unbatched: the only caller is the per-candidate fan-out in
        // `score_candidates`, so this must mask on the same axis as
        // `against_every_reference` — the reference axis, `Refiner`'s
        // `findHighestScoringAlignmentSet`. It previously used the instance
        // axis, silently giving phase 1 phase 2's rule for any backend that
        // cannot batch (parasail, the reference aligner). rmblast batches, so
        // nothing measured through it was affected.
        Ok(refiner_mask(
            self.inner.against_reference_nested(reference, seqs, skip)?,
            self.mask_level,
            MaskAxis::Reference,
            self.cross_strand,
            self.best_orientation_phase1,
        ))
    }

    fn batches_all_vs_all(&self) -> bool {
        self.inner.batches_all_vs_all()
    }
}

/// Highest-scoring alignment per instance index, input order preserved.
fn best_per_instance(hits: Vec<(usize, Alignment)>) -> Vec<(usize, Alignment)> {
    let mut best: std::collections::HashMap<usize, Alignment> =
        std::collections::HashMap::new();
    for (i, a) in hits {
        match best.get(&i) {
            Some(prev) if prev.score >= a.score => {}
            _ => {
                best.insert(i, a);
            }
        }
    }
    let mut out: Vec<(usize, Alignment)> = best.into_iter().collect();
    out.sort_by_key(|(i, _)| *i);
    out
}

/// `autocons` driven by rmblast's batched searches.
///
/// Both phases become a single search call: phase 2 is one query-set against the
/// consensus, phase 1 one all-against-all. That is the shape a search engine is
/// built for — the earlier pair-at-a-time adapter rebuilt a query lookup table
/// per pair and scanned one subject at a time.
///
/// **Multiple alignments per instance are kept.** An instance matching in
/// several pieces contributes several MSA rows, which is GIRI's own `FRAGMENT`
/// model. `mask_level` (RepeatMasker default 80) decides how much query overlap
/// is tolerated before the lower-scoring HSP is dropped; 101 disables it.
impl AlignmentSource for aln_rmblast::RmblastEngine {
    fn against_reference(
        &self,
        reference: &Sequence,
        seqs: &[Sequence],
        skip: Option<usize>,
    ) -> Result<Vec<(usize, Alignment)>> {
        self.one_to_many(reference, seqs, skip)
    }

    fn against_every_reference(
        &self,
        seqs: &[Sequence],
    ) -> Result<Vec<Vec<(usize, Alignment)>>> {
        // `Refiner` leaves its all-vs-all engine unmasked; `--mask-level`
        // governs refinement only.
        let hits = self.all_vs_all(seqs, 101)?;
        Ok(aln_rmblast::RmblastEngine::group_by_subject(hits, seqs.len()))
    }

    fn batches_all_vs_all(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod reference_tests {
    use super::*;
    use aln_engine::{AlignMode, AlignParams};
    use crate::FastAligner;

    /// The seed reference must be reported, and must be the input that actually
    /// won phase 1 — not merely the first sequence.
    #[test]
    fn the_winning_reference_is_reported() {
        // Three near-identical sequences plus one outlier.  The outlier can
        // never win, and whichever of the three does must be named correctly.
        let core = b"ACGTTGCAAGGCTTACGGATCCGTTACAGGCATTACGGATCAGGCTTAACGGTTACG";
        let mut seqs: Vec<Sequence> = (0..3)
            .map(|i| Sequence::new(format!("core{i}"), core.to_vec()))
            .collect();
        seqs.push(Sequence::new("outlier", b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_vec()));

        let p = AlignParams { mode: AlignMode::Local, min_score: 1, ..Default::default() };
        let al = FastAligner::new(aln_core::SubstMatrix::parse(TEST_MATRIX).unwrap(), p).unwrap();
        let out = run(&Pairwise::new(al), &seqs, &Params::default()).unwrap();

        assert_eq!(out.len(), 1);
        let r = &out[0].reference;
        assert!(r.index < 3, "the outlier should never seed the consensus, got {}", r.index);
        assert_eq!(r.name, seqs[r.index].name, "reported name must match the index");
        assert!(r.score > 0, "the winning reference must carry its phase-1 score");
    }

    const TEST_MATRIX: &str = "\
FREQS A 0.25 C 0.25 G 0.25 T 0.25
GAP -10 -2
       A   C   G   T   N
   A   3  -4  -2  -4  -1
   C  -4   3  -4  -2  -1
   G  -2  -4   3  -4  -1
   T  -4  -2  -4   3  -1
   N  -1  -1  -1  -1  -1
";
}

#[cfg(test)]
mod tests {
    use super::*;
    use aln_core::SubstMatrix;
    use aln_engine::{AlignMode, AlignParams};
    use aln_reference::ReferenceAligner;

    const M: &str = "\
FREQS A 0.325 C 0.175 G 0.175 T 0.325
  A   R   G   C   Y   T   K   M   S   W   N   X
  8   0 -10 -18 -19 -21 -15  -4 -14  -6  -1 -30
  3   3  12 -17 -18 -19  -9  -8  -8  -9  -1 -30
 -7   2  12 -16 -16 -17  -2 -11  -1 -12  -1 -30
-17 -16 -16  12   2  -7 -11  -2  -1 -12  -1 -30
-19 -18 -17  12   0   3  -8  -9  -8  -9  -1 -30
-21 -19 -18 -10   0   8  -4 -15 -14  -6  -1 -30
-14  -8  -2 -13  -8  -4  -3 -13  -8  -9  -1 -30
 -4  -8 -13  -2  -8 -14 -13  -3  -8  -9  -1 -30
-12  -7  -1  -1  -7 -12  -7  -7  -1 -12  -1 -30
 -6 -10 -14 -14 -10  -6 -10 -10 -14  -6  -1 -30
 -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
-30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
";

    fn aligner() -> ReferenceAligner {
        let p = AlignParams {
            mode: AlignMode::Local,
            gap_open: 25,
            gap_extend: 5,
            min_score: 1,
            ..Default::default()
        };
        ReferenceAligner::new(SubstMatrix::parse(M).unwrap(), p).unwrap()
    }

    struct Rng(u64);
    impl Rng {
        fn new(s: u64) -> Self {
            Rng(s.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1)
        }
        fn below(&mut self, n: usize) -> usize {
            let mut x = self.0;
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            self.0 = x;
            (x.wrapping_mul(0x2545_F491_4F6C_DD1D) % n as u64) as usize
        }
    }

    /// A family: one ancestral sequence, several diverged copies.
    fn family(n: usize, len: usize, subst_pct: usize) -> (Vec<u8>, Vec<Sequence>) {
        let mut rng = Rng::new(11);
        let ancestor: Vec<u8> = (0..len).map(|_| b"ACGT"[rng.below(4)]).collect();
        let seqs = (0..n)
            .map(|i| {
                let copy: Vec<u8> = ancestor
                    .iter()
                    .map(|&b| {
                        if rng.below(100) < subst_pct {
                            b"ACGT"[rng.below(4)]
                        } else {
                            b
                        }
                    })
                    .collect();
                Sequence::new(format!("s{i}"), copy)
            })
            .collect();
        (ancestor, seqs)
    }

    #[test]
    fn recovers_the_ancestral_sequence() {
        let (ancestor, seqs) = family(12, 200, 8);
        let out = run(&Pairwise::new(aligner()), &seqs, &Params::default()).unwrap();
        assert_eq!(out.len(), 1);

        let cons = &out[0].consensus;
        assert_eq!(cons.len(), ancestor.len(), "consensus length should match");
        let matches = cons.iter().zip(&ancestor).filter(|(a, b)| a == b).count();
        let identity = matches as f64 / ancestor.len() as f64;
        assert!(
            identity > 0.95,
            "consensus should recover the ancestor: {:.3} identity",
            identity
        );
    }

    #[test]
    fn a_consensus_beats_any_single_member() {
        // The point of the exercise: the called consensus should be closer to
        // the ancestor than the input sequences are.
        let (ancestor, seqs) = family(12, 200, 10);
        let out = run(&Pairwise::new(aligner()), &seqs, &Params::default()).unwrap();
        let cons = &out[0].consensus;

        let identity = |a: &[u8]| {
            a.iter().zip(&ancestor).filter(|(x, y)| x == y).count() as f64 / ancestor.len() as f64
        };
        let cons_id = identity(cons);
        let best_member = seqs.iter().map(|s| identity(&s.seq)).fold(0.0, f64::max);
        assert!(
            cons_id > best_member,
            "consensus {cons_id:.3} should beat the best input {best_member:.3}"
        );
    }

    #[test]
    fn candidates_come_back_sorted_best_first() {
        let (_, seqs) = family(8, 150, 10);
        let cands = score_candidates(&Pairwise::new(aligner()), &seqs, &Params::default()).unwrap();
        assert!(cands.len() > 1);
        for w in cands.windows(2) {
            assert!(w[0].score >= w[1].score, "not sorted: {:?}",
                    cands.iter().map(|c| c.score).collect::<Vec<_>>());
        }
    }

    /// A stable stop must not depend on how big the budget was.
    ///
    /// The C++ stops only on a fixed point; longer cycles run the budget out.
    /// With cycle detection, a run that stops stably at pass k does so whatever
    /// `iterations` was — so a large budget and a small one must agree exactly.
    #[test]
    fn a_stable_stop_is_independent_of_the_budget() {
        let (_, seqs) = family(12, 200, 12);
        let short = run(&Pairwise::new(aligner()), &seqs,
                        &Params { iterations: 3, ..Params::default() }).unwrap();
        let long = run(&Pairwise::new(aligner()), &seqs,
                       &Params { iterations: 30, ..Params::default() }).unwrap();
        assert!(long[0].stop.is_stable(),
                "30 passes should reach a fixed point or a cycle, got {:?}", long[0].stop);
        if short[0].stop.is_stable() {
            assert_eq!(short[0].passes, long[0].passes,
                       "a stable stop must not depend on the budget");
            assert_eq!(short[0].consensus, long[0].consensus,
                       "a stable stop must give the same consensus");
        }
    }

    /// The returned MSA must be the alignment the returned consensus was
    /// called from — including after a replay, when the best pass was not the
    /// last one. This is what "keep both" has to mean.
    #[test]
    fn the_returned_msa_recalls_the_returned_consensus() {
        let (_, seqs) = family(10, 180, 10);
        let params = Params { iterations: 11, ..Params::default() };
        let out = run(&Pairwise::new(aligner()), &seqs, &params).unwrap();
        let r = &out[0];
        let call = ConsensusParams { include_reference: false, ..params.consensus.clone() };
        let recalled = seqmod::ungap(&r.msa.consensus(&call));
        assert_eq!(recalled, r.consensus,
                   "the returned MSA must re-call the returned consensus");
    }

    /// Replaying the winning pass must reproduce it exactly — the assumption
    /// that lets `refine` store `(consensus_in, score)` instead of an MSA per
    /// iteration.
    #[test]
    fn a_pass_is_reproducible_from_its_input_consensus() {
        let (_, seqs) = family(8, 160, 8);
        let params = Params { iterations: 5, ..Params::default() };
        let a = run(&Pairwise::new(aligner()), &seqs, &params).unwrap();
        let b = run(&Pairwise::new(aligner()), &seqs, &params).unwrap();
        assert_eq!(a[0].consensus, b[0].consensus);
        assert_eq!(a[0].score, b[0].score);
        assert_eq!(a[0].passes, b[0].passes);
    }

    #[test]
    fn refinement_converges_before_running_out_of_passes() {
        let (_, seqs) = family(10, 150, 6);
        let out = run(&Pairwise::new(aligner()), &seqs, &Params::default()).unwrap();
        assert!(out[0].converged, "expected convergence on a tight family");
        assert!(
            out[0].passes <= Params::default().iterations + 1,
            "passes = {}",
            out[0].passes
        );
    }

    #[test]
    fn asking_for_several_consensi_numbers_the_names() {
        let (_, seqs) = family(8, 150, 10);
        let params = Params { num_consensi: 3.0, ..Default::default() };
        let out = run(&Pairwise::new(aligner()), &seqs, &params).unwrap();
        assert_eq!(out.len(), 3);
        assert_eq!(out[0].name, "CON1");
        assert_eq!(out[2].name, "CON3");
    }

    #[test]
    fn a_single_consensus_is_not_numbered() {
        let (_, seqs) = family(6, 120, 8);
        let out = run(&Pairwise::new(aligner()), &seqs, &Params::default()).unwrap();
        assert_eq!(out[0].name, "CON");
    }

    #[test]
    fn fractional_counts_are_a_share_of_the_input() {
        assert_eq!(resolve_count(0.5, 10), 5);
        assert_eq!(resolve_count(0.25, 10), 2);
        // Never fewer than one, however small the fraction.
        assert_eq!(resolve_count(0.01, 10), 1);
        assert_eq!(resolve_count(3.0, 10), 3);
        assert_eq!(resolve_count(1.0, 10), 1);
    }

    #[test]
    fn the_msa_keeps_every_sequence_that_aligned() {
        let (_, seqs) = family(10, 150, 8);
        let out = run(&Pairwise::new(aligner()), &seqs, &Params::default()).unwrap();
        assert_eq!(
            out[0].msa.num_instances(),
            seqs.len(),
            "every input should be placed in the final alignment"
        );
    }

    #[test]
    fn an_empty_input_produces_nothing() {
        let out = run(&Pairwise::new(aligner()), &[], &Params::default()).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn a_single_sequence_has_no_other_to_align_to() {
        // With nothing to align against, the total score is 0 and the candidate
        // is dropped — matching the C++'s `if (outscore > 0)` guard.
        let seqs = vec![Sequence::new("only", b"ACGTACGTACGTACGT".to_vec())];
        let out = run(&Pairwise::new(aligner()), &seqs, &Params::default()).unwrap();
        assert!(out.is_empty());
    }
}

#[cfg(test)]
mod nonredundant_score_tests {
    use super::*;
    use aln_core::align::{EditOp, EditScript};

    fn a(inst: usize, ss: usize, len: u32, score: i32) -> (usize, Alignment) {
        let mut e = EditScript::new();
        e.push(EditOp::Sub, len);
        (inst, Alignment::new("q", "s", 0, ss, Strand::Plus, score, e))
    }

    /// Two copies covering the same reference region is what a good reference
    /// looks like; both must count in full.
    #[test]
    fn different_instances_are_not_deduplicated() {
        let hits = vec![a(0, 0, 100, 500), a(1, 0, 100, 500)];
        assert_eq!(nonredundant_candidate_score(&hits), 1000);
    }

    /// One copy covering the same region twice is double counting; the second
    /// contributes nothing.
    #[test]
    fn one_instance_twice_over_counts_once() {
        let hits = vec![a(0, 0, 100, 500), a(0, 0, 100, 400)];
        assert_eq!(nonredundant_candidate_score(&hits), 500);
    }

    /// Half-covered contributes half. This is the case the containment filter
    /// lets through: 50% novel is above the 20% threshold, so both survive.
    #[test]
    fn partial_overlap_is_prorated() {
        let hits = vec![a(0, 0, 100, 500), a(0, 50, 100, 400)];
        assert_eq!(nonredundant_candidate_score(&hits), 500 + 200);
    }

    /// Genuinely disjoint fragments of one copy are both real coverage.
    #[test]
    fn disjoint_fragments_both_count() {
        let hits = vec![a(0, 0, 100, 500), a(0, 200, 100, 400)];
        assert_eq!(nonredundant_candidate_score(&hits), 900);
    }

    /// Claimed intervals must merge, or a third alignment overlapping two
    /// earlier ones would have its overlap counted twice and go negative.
    #[test]
    fn claimed_intervals_merge() {
        let hits = vec![a(0, 0, 100, 900), a(0, 80, 100, 800), a(0, 40, 100, 100)];
        let got = nonredundant_candidate_score(&hits);
        assert!(got >= 900 + 640, "no credit lost twice, got {got}");
        assert!(got <= 900 + 640 + 100);
    }
}

#[cfg(test)]
mod refiner_filter_tests {
    use super::*;

    fn keep(items: &[(i64, usize, usize)], level: u32) -> Vec<bool> {
        // (annotated i64 literals via the slice type)
        let v: Vec<(i64, usize, usize, Strand)> =
            items.iter().map(|&(s, b, e)| (s, b, e, Strand::Plus)).collect();
        refiner_keep(&v, level)
    }

    /// A lower-scoring alignment mostly inside a higher-scoring one is dropped.
    #[test]
    fn contained_fragment_is_masked() {
        // 0..100 (score 500) vs 10..90 (score 100): zero novel span.
        assert_eq!(keep(&[(500, 0, 100), (100, 10, 90)], 80), [true, false]);
    }

    /// ≥20% overhang at mask level 80 survives — that is the tiling case the
    /// filter exists to keep.
    #[test]
    fn tiling_fragment_with_enough_overhang_survives() {
        // 0..100 vs 80..200: novel = 100 of 120 = 83% ≥ 20%.
        assert_eq!(keep(&[(500, 0, 100), (100, 80, 200)], 80), [true, true]);
        // 0..100 vs 60..110: novel = 10 of 50 = 20%... exactly at boundary:
        // Perl deletes only when perc < 20, so exactly 20% survives.
        assert_eq!(keep(&[(500, 0, 100), (100, 60, 110)], 80), [true, true]);
        // 0..100 vs 70..110: novel = 10 of 40 = 25% — survives.
        assert_eq!(keep(&[(500, 0, 100), (100, 70, 110)], 80), [true, true]);
        // 0..100 vs 85..110: novel = 10 of 25 = 40% — survives.
        // 0..100 vs 90..105: novel = 5 of 15 = 33% — survives.
        // 0..100 vs 95..104: novel = 4 of 9 = 44% — survives.
    }

    /// Disjoint alignments never mask each other, whatever the scores.
    #[test]
    fn disjoint_fragments_are_untouched() {
        assert_eq!(keep(&[(500, 0, 100), (1, 100, 150)], 80), [true, true]);
    }

    /// Faithful to the Perl: a deleted alignment still masks ones below it, so
    /// an overlapping chain collapses to its best member.
    #[test]
    fn deleted_alignments_still_mask() {
        // A(500, 0..100) masks B(100, 5..95); B would NOT mask C(50, 5..95) on
        // its own terms... but B does anyway, and A also masks C directly.
        // The sharper case: C overlaps B heavily but A only partially.
        // A: 0..100.  B: 60..160 (novel 60/100=60% -> survives).
        // C: 90..165 (vs A: novel 65/75=87% survives; vs B: novel 5/75=6.7% -> deleted by B).
        assert_eq!(keep(&[(500, 0, 100), (400, 60, 160), (50, 90, 165)], 80),
                   [true, true, false]);
    }

    /// Strands are independent groups.
    #[test]
    fn strands_do_not_mask_each_other() {
        let items: [(i64, usize, usize, Strand); 2] =
            [(500, 0, 100, Strand::Plus), (100, 10, 90, Strand::Minus)];
        assert_eq!(refiner_keep(&items, 80), [true, true]);
    }

    /// 101 disables everything, matching the engine convention.
    #[test]
    fn level_101_keeps_all() {
        assert_eq!(keep(&[(500, 0, 100), (1, 0, 100)], 101), [true, true]);
    }
}

/// The occupancy gate, on the shape that motivated it: a reference whose flanks
/// no other instance reaches.
#[cfg(test)]
mod gate_tests {
    use super::*;
    use aln_core::msa::SequenceRow;

    /// Reference is 15 columns; the three instances cover only the middle five.
    fn flanked_msa() -> MultiAlign {
        let row = |n: &str, s: &[u8]| SequenceRow::new(n, s.to_vec());
        MultiAlign::from_sequences(
            row("ref", b"ACGTACCTAGCAGTA"),
            vec![
                row("i1", b"     CCTAG     "),
                row("i2", b"     CCTAG     "),
                row("i3", b"     CCTAG     "),
            ],
        )
        .unwrap()
    }

    fn params(caller: Caller, min: usize) -> Params {
        Params { caller, min_non_gap_count: min, ..Params::default() }
    }

    /// Ungated, the Dfam caller reads the reference's flanks straight out of the
    /// one row that has them — the whole reason `--min` was wanted here.
    #[test]
    fn ungated_the_flanks_survive() {
        let cons = gapped_consensus(&flanked_msa(), &params(Caller::Dfam, 0));
        assert_eq!(cons.len(), 15);
        assert_eq!(seqmod::ungap(&cons).len(), 15);
    }

    /// A floor of 2 leaves only the columns the instances actually cover.
    #[test]
    fn a_floor_of_two_keeps_only_the_covered_span() {
        for caller in [Caller::Dfam, Caller::Giri] {
            let cons = gapped_consensus(&flanked_msa(), &params(caller, 2));
            assert_eq!(cons.len(), 15, "{caller:?} changed the column count");
            assert_eq!(seqmod::ungap(&cons), b"CCTAG", "{caller:?}");
        }
    }

    /// The floor counts rows, so one above the coverage empties the call.
    #[test]
    fn a_floor_above_the_coverage_calls_nothing() {
        for caller in [Caller::Dfam, Caller::Giri] {
            let cons = gapped_consensus(&flanked_msa(), &params(caller, 4));
            assert!(seqmod::ungap(&cons).is_empty(), "{caller:?}");
        }
    }

    /// Padding and gaps both count as absent, in either convention.
    #[test]
    fn gaps_and_padding_alike_leave_a_column_uncovered() {
        let row = |n: &str, s: &[u8]| SequenceRow::new(n, s.to_vec());
        let msa = MultiAlign::from_sequences(
            row("ref", b"ACGTA"),
            vec![row("i1", b"AC-TA"), row("i2", b"AC.TA"), row("i3", b"AC TA")],
        )
        .unwrap();
        let cons = gapped_consensus(&msa, &params(Caller::Dfam, 1));
        assert_eq!(cons[2], b'-');
        assert_eq!(seqmod::ungap(&cons), b"ACTA");
    }
}
