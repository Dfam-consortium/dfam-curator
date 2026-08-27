//! `te-composer` — TE Composer.
//!
//! Consensus building beyond the GIRI `autocons` workflow. Where `autocons`
//! exists to reproduce the C++ exactly and must not change, this is where
//! changes go: score cutoffs, unscaled matrices, tunable iteration, alternative
//! search backends, Stockholm output.
//!
//! It started as a copy of the `autocons` CLI and still carries most of it;
//! options will diverge as improvements land.

use std::io::Write;
use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};

use aln_core::consensus::ConsensusParams;
use aln_core::io;
use aln_core::msa::InsertionPolicy;
use aln_core::SubstMatrix;
use aln_engine::{AlignMode, AlignParams};
use aln_engine::engine::{ScoreMode, SearchParams};
use cons_core::FastAligner;
use aln_reference::ReferenceAligner;
use aln_rmblast::{RmblastEngine, RmblastOptions};
use cons_core::{BestHsp, Census, MinScore, RefinerFilter, run, Caller, Pairwise, Params};

mod extend;

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum HspPolicy {
    /// One alignment per instance — matches GIRI's workflow.
    ///
    /// An instance that aligns in several fragments contributes only its
    /// highest-scoring one; the rest of its coverage is discarded.
    Best,
    /// Tile each instance's alignments, one strand per instance. **The default.**
    ///
    /// Per instance, in both phases: trim each strand separately so that
    /// strand's instance ranges do not overlap, then sum the trimmed scores and
    /// keep the higher-scoring tiling. A trimmed alignment carries a
    /// proportionally reduced score, so the strand vote counts only evidence
    /// that becomes rows.
    ///
    /// An instance contributes one or more rows, all on the same strand,
    /// together covering its sequence without overlap. Overlap on the
    /// *reference* axis is deliberately left alone — one instance covering the
    /// same consensus region twice from different bases is a tandem expansion,
    /// not a defect.
    Tiled,
    /// Every HSP surviving `--mask-level`, as separate MSA rows — Refiner's
    /// behaviour.
    ///
    /// Keeps all the coverage, but the same genomic bases can appear in two
    /// rows, which claims they are homologous to two places at once.
    All,
}

/// CLI spelling of [`cons_core::AcceptRule`].
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum AcceptArg {
    /// Summed alignment score.
    Total,
    /// Mean over instances of score per aligned base.
    MeanPerBase,
    /// No gate: keep every repair. Diagnostic only.
    Always,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum Backend {
    /// parasail striped SIMD. Exact Smith-Waterman, and memory-hungry.
    ///
    /// Peak memory is about `4 bytes x min(threads, sequences) x length^2`, so
    /// the thread count is chosen to fit `--max-memory` unless `--threads`
    /// overrides it. Against rmblast it buys a fraction of a point of consensus
    /// identity for 6-36x the runtime and up to 150x the memory.
    Parasail,
    /// Plain O(mn) scalar aligner — slower, used to arbitrate between the
    /// other two when they disagree. Hidden: it exists to referee, not to run
    /// production work.
    #[value(hide = true)]
    Reference,
    /// Seeded search via the rmblast port, best HSP per pair. **The default.**
    ///
    /// Not equivalent to the DP backends: seeding can miss weak relationships
    /// that full dynamic programming finds. Measured against parasail on
    /// simulated families, the two agree to within 0.3 points of consensus
    /// identity up to 20% divergence; above that parasail recovers roughly 3
    /// points more, at 6-36x the runtime and up to 150x the memory.
    Rmblast,
}

/// How insertions are merged into the alignment.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum Insertions {
    /// GIRI's incremental merge (`adjustReference`) — the C++ behaviour, and
    /// the default so output is comparable with `bin/autocons`.
    Incremental,
    /// Per-slot maximum: each member's insertions stay independent, so a shared
    /// column always implies positional homology. Slightly wider.
    PerSlot,
    /// Discard insertions; width stays equal to the reference.
    Drop,
}

impl From<Insertions> for InsertionPolicy {
    fn from(i: Insertions) -> Self {
        match i {
            Insertions::Incremental => InsertionPolicy::GrowIncremental,
            Insertions::PerSlot => InsertionPolicy::GrowPerSlot,
            Insertions::Drop => InsertionPolicy::Drop,
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum OutFormat {
    /// GIRI `.ig` (the C++ default).
    Ig,
    Fasta,
}

/// Build the N best-scoring consensus sequences from a set of related sequences.
#[derive(Parser, Debug)]
#[command(
    name = "te-composer",
    about = "Build TE consensus sequences, with improvements beyond GIRI autocons.",
    version
)]
struct Cli {
    /// Input FASTA. Use `-` for stdin.
    input: String,

    /// Supply the bootstrap consensus instead of deriving one.
    ///
    /// The bootstrap phase exists only to pick a starting consensus from the
    /// instances; given one, it is skipped and the run goes straight to
    /// refinement. Use it to improve an existing family against a new instance
    /// set, or to re-refine a curated consensus, rather than rediscovering a
    /// reference the instances may not agree on. The first record of the FASTA
    /// is used.
    #[arg(long = "bootstrap-cons", value_name = "FASTA", conflicts_with = "num")]
    bootstrap_cons: Option<PathBuf>,

    /// Suppress the progress report on stderr.
    #[arg(long, default_value_t = false)]
    silent: bool,

    /// Stockholm output: the final alignment plus its consensus as `#=GC RF`,
    /// one record per consensus. Defaults to stdout.
    ///
    /// This is the primary output — it carries the alignment *and* the
    /// consensus, so nothing is lost. Use `--consensus` when a bare consensus
    /// file is also wanted.
    output: Option<PathBuf>,

    /// Number of consensi to emit. Below 1.0 it is a fraction of the input count.
    #[arg(short = 'n', long, default_value_t = 1.0)]
    num: f64,

    /// Base name for emitted sequences.
    #[arg(long, default_value = "CON")]
    name: String,

    /// Refinement passes after the first.
    ///
    /// Defaults to **3 for the DP backends** (the C++ `autocons` value) and
    /// **20 for rmblast**. A pass is cheap for a seeded search and expensive
    /// for full dynamic programming, so the budget that is reasonable differs
    /// by roughly the cost ratio between them.
    ///
    /// 20 is measured, not guessed: over 486 TFE simulations a budget of 3
    /// left 55% of families still refining, 11 left 15%, and 20 leaves 7%.
    /// The tail is concentrated in long, structured families at mid-to-high
    /// divergence (Charlie1 needs a median of 12 passes; AluY needs 2), which
    /// are exactly the ones that most need the refinement. Raising it is
    /// close to free because 63% of families stop by pass 5 and never see the
    /// extra budget.
    ///
    /// Raising it is safe: refinement stops as soon as the consensus reaches a
    /// fixed point *or* repeats one seen earlier (a cycle), so a larger budget
    /// buys more passes only where they are still changing the answer. Runs
    /// that stop stably give the same result at 3 as at 30.
    #[arg(long)]
    iterations: Option<usize>,

    /// Alignment matrix, in crossmatch layout.
    #[arg(long)]
    matrix: Option<PathBuf>,

    /// Gap-open penalty, as a positive magnitude.
    ///
    /// Defaults to the matrix's own `GAP` line when it has one, so a matrix and
    /// its penalties cannot drift apart. Falls back to 25.
    #[arg(long)]
    gap_open: Option<u32>,

    /// Gap-extension penalty, as a positive magnitude. Defaults to the matrix's
    /// `GAP` line, else 5.
    #[arg(long)]
    gap_extend: Option<u32>,

    /// Discard alignments below this score. Default 150, for every backend.
    ///
    /// 150 is Refiner's own setMinScore(150), in units of the default matrix
    /// (mean diagonal 9.5), so it means the same thing to either backend.
    /// RepeatMasker's 200 drops short families outright.
    ///
    /// It matters more to the seeded search than to full DP: rmblast derives
    /// its X-drop cutoffs from this floor, so a low value collapses gapped
    /// extension and returns roughly a hundred junk HSPs per pair — 11.0 s and
    /// 2.7 GB at a floor of 1 against 1.5 s and 88 MB at 150, on the same
    /// input. A DP backend reports one optimal alignment per pair either way,
    /// so for it this is purely about what enters the MSA.
    #[arg(long)]
    min_score: Option<i32>,

    /// Seed word size. `--backend rmblast` only.
    #[arg(long, default_value_t = 7)]
    word_size: u32,

    /// Turn off DUST low-complexity query masking.
    ///
    /// DUST is **on** by default here, which differs from RepeatMasker and
    /// `Refiner`: both emit `-dust no` unconditionally
    /// (`NCBIBlastSearchEngine.pm`, and note its `setUseDustSeg(1)` is dead —
    /// `getUseDustSeg` is never read). They lean on complexity-adjusted scoring
    /// instead, which the Rust rmblast port currently applies as a no-op, so
    /// without DUST a low-complexity family has no defence at all: 60 GC-rich
    /// 20 kb instances produced 8.7M HSPs and 25 GB of peak memory.
    ///
    /// Measured cost on a well-behaved family: 0.7% of HSPs.
    #[arg(long, hide = true)]
    no_dust: bool,

    /// How an instance's alignments become MSA rows. `--backend rmblast` only.
    ///
    /// A seeded search reports several HSPs per instance, and the three
    /// policies differ in what they do with them. `best` keeps one and
    /// discards the instance's other coverage; `all` keeps them all and lets
    /// the same genomic bases appear in two rows; `tiled`, the default, keeps
    /// them all from the winning strand but tiles them so no base is used
    /// twice.
    ///
    /// The DP backends report one alignment per pair, so this does not apply
    /// to them — they are always effectively `best`.
    ///
    /// `tiled` trims each strand into a tiling and then keeps the
    /// higher-scoring one, so every score it reports or ranks on is a trimmed
    /// score and one convention serves both phases.
    #[arg(long, value_enum, default_value_t = HspPolicy::Tiled)]
    hsps: HspPolicy,

    /// Repair low-quality blocks between two rounds of refinement.
    ///
    /// A port of `Refiner`'s `resolveLowQualityBlocks`. Refinement aligns every
    /// instance to one consensus and re-calls column-wise; where instances
    /// disagree about an indel, the aligner places the same event at different
    /// offsets in different rows and the caller averages incompatible registers
    /// into a consensus matching none of them. Iterating cannot fix that — each
    /// pass re-derives the same columns from the same alignment. This finds
    /// those stretches (Ruzzo-Tompa over an inverted column-quality profile),
    /// re-derives the consensus for each from the instances directly, splices
    /// the results in, and refines again.
    ///
    /// The second round is kept only if it scores better than the first, so
    /// this cannot lower the reported score.
    ///
    /// **On by default**, as policy D — the low-quality selection plus the
    /// scanning window, both resolved by the length vote. Measured through this
    /// exact path on both benchmarks, it was the only policy that beat the
    /// low-quality selection alone on each.
    #[arg(long, default_value_t = false)]
    no_repair_blocks: bool,

    /// Do not re-derive spans the consensus has no bases for.
    ///
    /// Packing is on by default. A consensus-induced MSA never aligns one
    /// instance's inserted bases to another's, so an insertion inherited by many
    /// instances is invisible to the consensus caller however many passes run.
    /// Re-deriving those spans inside each pass lets a recovered base join the
    /// reference for the next one, which compounds. Measured to help against
    /// known ancestors and to be roughly neutral on curated hs1 families.
    #[arg(long, default_value_t = false)]
    no_pack_insertions: bool,

    /// Merge gap runs separated by at most this many called columns.
    #[arg(long, default_value_t = 2, hide = true)]
    pack_max_sep: usize,

    /// A span needs one instance contributing at least this many bases.
    #[arg(long, default_value_t = 5, hide = true)]
    pack_min_seg: usize,

    /// The all-against-all winner must beat this summed score for a span to be
    /// re-derived. 0 means it must align positively to the other instances.
    #[arg(long, default_value_t = 0, hide = true)]
    pack_min_score: i64,

    /// Genome the input sequences were taken from, as a 2bit. Given one,
    /// the finished consensus is extended past its edges with RAMExtend and
    /// then refined again against the widened instances.
    ///
    /// Input names must be Smitten identifiers — `chr1:1000-2000_+`, or
    /// `hg38:chr1:1000-2000_+` — since that is what locates an instance in the
    /// genome. Copies that cannot be located, or whose bases do not match the
    /// genome there, are reported and left out of the extension; they stay in
    /// the family.
    #[arg(long, value_name = "2BIT")]
    genome: Option<PathBuf>,

    /// Assembly the genome represents. When set, an input whose identifier
    /// carries a different assembly prefix is not extended.
    #[arg(long, value_name = "ID", requires = "genome")]
    assembly: Option<String>,

    /// Maximum extension per side.
    #[arg(long, default_value_t = 20_000, value_name = "BP", requires = "genome")]
    extend_max: i32,

    /// Band half-width for the extension alignment.
    #[arg(long, default_value_t = 40, requires = "genome")]
    extend_bandwidth: i32,

    /// Copies that must reach an edge before extension is attempted; also the
    /// multiplier in the per-column score improvement the extension demands.
    #[arg(long, default_value_t = 3, requires = "genome")]
    extend_min_seqs: i32,

    /// Extension matrix. Default picks one from the family's Kimura
    /// divergence, as `extend-stk.pl` does.
    #[arg(long, value_name = "NAME", requires = "genome")]
    extend_matrix: Option<String>,

    /// Score improvement per consensus column the extension must sustain.
    /// Default is `--extend-min-seqs` x the matrix diagonal average.
    #[arg(long, value_name = "N", requires = "genome")]
    extend_min_improvement: Option<i32>,

    /// Columns of slack allowed between an instance's alignment edge and the
    /// consensus edge before it stops driving the extension on that side.
    #[arg(long, default_value_t = 10, requires = "genome")]
    extend_edge_slop: usize,

    /// Refuse an extension whose two sides together add more than this.
    ///
    /// An extension is also refused outright if either side reaches
    /// `--extend-max`, since a side that was stopped has not found an edge.
    #[arg(long, default_value_t = 25_000, value_name = "BP", requires = "genome", hide = true)]
    extend_max_total: i32,

    /// Fraction of comparable bases an input must share with the genome at
    /// its stated coordinates to take part in the extension.
    #[arg(long, default_value_t = 0.95, value_name = "F", requires = "genome")]
    extend_min_identity: f64,

    /// Extend but do not refine afterwards: the extension is spliced onto the
    /// consensus and reported as-is.
    #[arg(long, default_value_t = false, requires = "genome", hide = true)]
    extend_no_refine: bool,

    /// Ruzzo-Tompa segment-score threshold for calling a block low-quality.
    /// Higher is more conservative. `--repair-blocks` only.
    #[arg(long, default_value_t = 1.0, hide = true)]
    repair_threshold: f64,

    /// Also take repair candidates from a scanning window of this many
    /// consensus positions (`AutoRunBlocker`'s selection), resolved by the same
    /// length vote. 0 disables it. `--repair-blocks` only.
    ///
    /// Measured on 791 hs1 families against curated Dfam consensi: the
    /// low-quality selection alone recovers +99 net bases, adding this window
    /// takes it to +223. 44% of what the window finds sits inside a low-quality
    /// block whose whole-block vote reports the consensus already agrees — the
    /// disagreement is only visible at window scale.
    ///
    /// 10 is the default and is what policy D means; 0 falls back to the
    /// low-quality selection alone (policy A).
    #[arg(long, default_value_t = 10)]
    repair_window: usize,

    /// Shortest row worth keeping after `--trim-reused`.
    #[arg(long, default_value_t = 25, hide = true)]
    min_row_len: usize,

    /// Write every refinement pass's consensus to `<FILE>` as FASTA.
    ///
    /// One record per pass, headed with that pass's score. Refinement selects
    /// by score, which measures how well a consensus explains the instances —
    /// not how close it is to the truth. Where a truth is known, this makes the
    /// two comparable.
    #[arg(long, value_name = "FASTA", hide = true)]
    dump_refinement: Option<PathBuf>,

    /// Write a TSV census of every alignment entering the MSA.
    ///
    /// One row per alignment: phase, reference index, query index, score,
    /// alignment columns, query span, subject span. Use it to see what a low
    /// `--min-score` actually admits.
    #[arg(long, value_name = "TSV", hide = true)]
    dump_alignments: Option<PathBuf>,

    /// Disable Phil Green's complexity-adjusted scoring.
    /// `--backend rmblast` only.
    ///
    /// Adjustment is ON by default, matching `Refiner`, which sets
    /// `complexityAdjustedScoreMode` on both its engines — so score floors and
    /// reference ranking are all in adjusted units, and low-complexity matches
    /// cannot inflate a candidate's phase-1 sum. (Before 2026-08-14 the default
    /// was off, as `--complexity-adjust`, to keep the GIRI-aligner comparison
    /// varying one thing at a time; that comparison is done.)
    #[arg(long = "no-complexity-adjust", default_value_t = true,
          action = clap::ArgAction::SetFalse, hide = true)]
    complexity_adjust: bool,

    /// Mammalian sequences: run GIRI's species-aware CpG restoration on the
    /// final consensus. `--orig` only; the Dfam caller restores CpG on every
    /// pass regardless.
    #[arg(long, hide = true)]
    mam: bool,

    /// Skip the CpG restoration pass. Dfam caller only.
    #[arg(long)]
    no_cpg: bool,

    /// Use the original GIRI consensus caller instead of the Dfam one.
    ///
    /// Note this applies it to *both* phases. The C++ `--orig` is the same in
    /// effect, since its phase 1 always used the GIRI caller anyway.
    #[arg(long)]
    orig: bool,

    /// Minimum non-gap residues in an alignment column. Applies only with
    /// `--orig`, matching the C++.
    #[arg(long, default_value_t = 2, hide = true)]
    min: usize,

    /// Alignment backend. Defaults to the seeded search, which is adequate
    /// across the whole size range measured and orders of magnitude cheaper.
    ///
    /// Full dynamic programming buys about 3 points of consensus identity, but
    /// only above ~20% divergence; below that the two agree to within 0.3
    /// points at equal coverage.
    #[arg(long, value_enum, default_value_t = Backend::Rmblast)]
    backend: Backend,

    /// How insertions are merged into the multiple alignment.
    #[arg(long, value_enum, default_value_t = Insertions::Incremental, hide = true)]
    insertions: Insertions,

    /// Also write the bare consensus here, in `--format`.
    ///
    /// The Stockholm output already contains the consensus; this is for
    /// downstream tools that want it on its own.
    #[arg(long, value_name = "FILE")]
    consensus: Option<PathBuf>,

    /// Format for `--consensus`.
    #[arg(long, value_enum, default_value_t = OutFormat::Ig)]
    format: OutFormat,

    /// Write each final multiple alignment to `<PREFIX>.<name>` as aligned
    /// FASTA. Superseded by the Stockholm output; kept for scripts that read
    /// aligned FASTA.
    #[arg(long, hide = true)]
    aln: Option<String>,

    /// Run with exactly N threads, BYPASSING the memory guard.
    ///
    /// Without this flag, and with `--parasail`, the thread count is chosen to
    /// fit `--max-memory` — which is the recommended way to run. Passing
    /// `--threads` takes that protection off: peak memory is roughly
    /// `4 bytes x threads x length^2`, so 64 threads on a 16 kb element wants
    /// about 63 GB, and exceeding what the machine has will get the run killed
    /// rather than slowed.
    ///
    /// Capped at the number of input sequences either way — reference
    /// selection cannot use more concurrency than there are sequences to try.
    ///
    /// 0 (the default) means "decide for me".
    #[arg(long, default_value_t = 0)]
    threads: usize,

    /// Memory budget for `--parasail`, as a size (`48G`, `512M`) or a
    /// percentage of MemAvailable (`80%`).
    ///
    /// The thread count is derived from it. Ignored when `--threads` is given,
    /// and irrelevant to the seeded-search backend, whose memory is a small
    /// fraction of this and is not driven by thread count.
    ///
    /// The default is a percentage rather than a fixed size because it is
    /// measured against MemAvailable, which already excludes memory in use —
    /// so it is not the same quantity as the "half of RAM" that long-running
    /// daemons reserve to protect the page cache.
    #[arg(long, default_value = "80%")]
    max_memory: String,

    /// Line width for sequence output; 0 disables wrapping.
    #[arg(long, default_value_t = 60)]
    width: usize,
}

/// The default matrix, used when `--matrix` is not given.
///
/// RepeatModeler's `comparison.matrix` — the matrix `Refiner` hands to rmblastn
/// for exactly this job, all-against-all among the instances of one family.
///
/// # Why not a more permissive matrix, given the benchmark says so
///
/// Mismatch tolerance is the largest quality lever measured on the TFE
/// benchmark: GIRI's xrepmask (mism/match ~ -0.92) beats this one by up to
/// +0.077 normalised score at 30% divergence (~30 bp/kb more root coverage,
/// ~12/kb fewer wrong bases; p ~ 1e-19), because a strict matrix drains score
/// through diverged regions faster and trips x-drop truncation sooner.
///
/// It is kept as the default anyway, deliberately, pending a decision on which
/// permissive matrix to adopt. 25p43g was tried and rejected on measurement:
/// it sits between the two in tolerance (-1.22) and recovers only about a
/// fifth of the gap, while being significantly *worse* than xrepmask from 20%
/// divergence up. The live candidate is xrepmask **unscaled** (every entry and
/// gap cost x3), which is the same matrix as xrepmask/3 — identical optimal
/// alignments — but with a mean diagonal of 9.0, so a score threshold means
/// what it looks like.
///
/// No `GAP` line, so gap costs fall back to 25/5 — matching `Refiner`'s
/// `setGapInit(-25)` / `setInsGapExt(-5)`.
const DEFAULT_MATRIX: &str = "\
FREQS A 0.265 C 0.235 G 0.235 T 0.265
  A   R   G   C   Y   T   K   M   S   W   N
  9   1  -6 -15 -16 -17 -12  -2 -10  -4  -1
  1   1   1 -15 -15 -16  -6  -6  -6  -7  -1
 -6   1  10 -15 -15 -15  -2 -10  -2 -10  -1
-15 -15 -15  10   2  -6  -9  -2  -2  -9  -1
-16 -15 -15   1   1   1  -6  -7  -7  -7  -1
-17 -16 -15  -6   1   9  -2 -12 -11  -4  -1
-12  -6  -2 -11  -6  -2  -2 -11  -7  -7  -1
 -2  -6 -10  -2  -7 -12 -11  -2  -7  -7  -1
-10  -6  -2  -2  -7 -11  -7  -7  -2 -10  -1
 -4  -7 -10 -11  -7  -4  -7  -7 -10  -4  -1
 -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1
";

/// Bytes of scratch a single in-flight alignment needs, as a function of the two
/// sequence lengths.
///
/// parasail's striped **traceback** kernels allocate a full `m * n` traceback
/// matrix per alignment — score-only kernels do not, but the consensus
/// pipeline needs the path.
///
/// # Why 4 bytes per cell, and why not 2
///
/// The real figure is **bimodal**, because parasail escalates lane width when a
/// score will not fit: 16-bit lanes cost ~2 bytes/cell, and once the best
/// alignment score exceeds 32,767 it re-runs the whole matrix at 32 bits for
/// ~4. Measured across 41 runs spanning 8-32 kb, `peak / (min(threads, seqs) *
/// L^2)` lands at 1.6-2.0 or 3.3-3.9 with nothing in between, and which one it
/// is tracks the predicted lane width in 37 of them.
///
/// This deliberately always assumes the expensive branch. Guessing 16-bit and
/// being wrong means using **double** the predicted memory, and no budget
/// fraction protects against a 2x miss; guessing 32-bit and being wrong only
/// costs some parallelism. Escalation also happens *inside* an alignment, so by
/// the time it is observable the memory is already committed — there is nothing
/// to re-plan.
///
/// Note the 32,767 boundary is in units of the matrix in play. A larger
/// diagonal escalates at shorter lengths, so this is an upper bound rather than
/// a prediction for any particular matrix.
const BYTES_PER_CELL: u64 = 4;

fn alignment_scratch_bytes(m: usize, n: usize) -> u64 {
    (m as u64)
        .saturating_mul(n as u64)
        .saturating_mul(BYTES_PER_CELL)
}

/// Concurrency the reference-selection phase can actually use.
///
/// Phase 1 parallelises over *candidate references*, so a family of 25
/// sequences never has more than 25 alignments in flight however many cores
/// exist. Budgeting on raw thread count over-predicts memory by up to 2.6x on
/// shallow input and would refuse runs that comfortably fit.
fn effective_threads(threads: usize, n_seqs: usize) -> usize {
    threads.min(n_seqs).max(1)
}

/// Parse `--max-memory`: `48G`, `512M`, `80%`, or a bare byte count.
fn parse_memory_budget(spec: &str, available: Option<u64>) -> Result<u64> {
    let s = spec.trim();
    if let Some(pct) = s.strip_suffix('%') {
        let pct: f64 = pct
            .trim()
            .parse()
            .with_context(|| format!("--max-memory {spec}: not a percentage"))?;
        if !(0.0..=100.0).contains(&pct) {
            bail!("--max-memory {spec}: percentage must be between 0 and 100");
        }
        let avail = available.ok_or_else(|| {
            anyhow::anyhow!(
                "--max-memory {spec} is a percentage of MemAvailable, which is \
                 only readable on Linux; give an absolute size instead"
            )
        })?;
        return Ok((avail as f64 * pct / 100.0) as u64);
    }
    let (num, mult) = match s.chars().last() {
        Some('G') | Some('g') => (&s[..s.len() - 1], 1024u64 * 1024 * 1024),
        Some('M') | Some('m') => (&s[..s.len() - 1], 1024u64 * 1024),
        Some('K') | Some('k') => (&s[..s.len() - 1], 1024u64),
        _ => (s, 1),
    };
    let v: f64 = num
        .trim()
        .parse()
        .with_context(|| format!("--max-memory {spec}: not a size"))?;
    Ok((v * mult as f64) as u64)
}

fn gb(bytes: u64) -> f64 {
    bytes as f64 / (1024.0 * 1024.0 * 1024.0)
}

/// How many threads parasail should run with, and what to tell the user.
///
/// Never returns an error and never suggests refusing: if the budget cannot fit
/// even one thread, it returns one and warns. One thread is the floor of what
/// there is to decide, and a machine with swap — or a user watching — may cope
/// perfectly well. A run that was warned about teaches more than one that
/// refused to start.
fn plan_threads(
    longest: usize,
    n_seqs: usize,
    requested: Option<usize>,
    budget: u64,
    available: Option<u64>,
) -> (usize, Vec<String>) {
    let per = alignment_scratch_bytes(longest, longest);
    let mut notes = Vec::new();

    if let Some(n) = requested {
        let eff = effective_threads(n, n_seqs);
        if eff < n {
            notes.push(format!(
                "note: --threads {n} reduced to {eff} — reference selection cannot \
                 use more concurrency than there are input sequences"
            ));
        }
        notes.push("note: --threads given; --max-memory not applied".to_string());
        let want = per.saturating_mul(eff as u64);
        if available.is_some_and(|a| want > a) {
            notes.push(format!(
                "warning: {eff} threads x {longest} bp needs roughly {:.1} GB, \
                 MemAvailable is {:.1} GB — the run may be killed",
                gb(want),
                gb(available.unwrap())
            ));
        }
        return (eff, notes);
    }

    let cap = (budget / per.max(1)) as usize;
    let eff = effective_threads(cap.min(num_cpus()), n_seqs);
    if cap == 0 {
        notes.push(format!(
            "warning: 1 thread x {longest} bp needs roughly {:.1} GB, over the \
             {:.1} GB budget",
            gb(per),
            gb(budget)
        ));
        notes.push(
            "warning: running at 1 thread anyway — the budget cannot be met. \
             Raise --max-memory, or use the default seeded-search backend."
                .to_string(),
        );
    } else if eff < num_cpus().min(n_seqs) {
        notes.push(format!(
            "note: using {eff} of {} threads to stay inside the {:.1} GB budget \
             (about {:.1} GB); --max-memory raises it, --threads overrides it",
            num_cpus().min(n_seqs),
            gb(budget),
            gb(per.saturating_mul(eff as u64))
        ));
    }
    (eff, notes)
}

fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

/// `MemAvailable` from `/proc/meminfo`, in bytes. `None` off Linux.
fn available_memory() -> Option<u64> {
    let text = std::fs::read_to_string("/proc/meminfo").ok()?;
    let line = text.lines().find(|l| l.starts_with("MemAvailable:"))?;
    let kb: u64 = line.split_whitespace().nth(1)?.parse().ok()?;
    Some(kb * 1024)
}

/// Report the plan for a parasail run, and return the thread count to use.
///
/// Only parasail is planned this way. The seeded-search backend holds hits
/// rather than a matrix, so its memory is a small fraction of this and is
/// driven by sequence count rather than thread count — there is nothing here
/// for a budget to control.
fn plan_parasail(seqs: &[aln_core::Sequence], cli: &Cli) -> Result<usize> {
    let longest = seqs.iter().map(|s| s.len()).max().unwrap_or(0);
    let available = available_memory();
    let budget = parse_memory_budget(&cli.max_memory, available)?;
    let requested = (cli.threads > 0).then_some(cli.threads);
    let (threads, notes) = plan_threads(longest, seqs.len(), requested, budget, available);
    for n in notes {
        eprintln!("{n}");
    }
    Ok(threads)
}

/// Run the pipeline, optionally recording every alignment on the way through.
/// Turn an instance's surviving HSPs into MSA rows, per `--hsps`.
///
/// One axis, deliberately. These policies constrain one another: reused-base
/// trimming cannot do anything once `Best` has collapsed an instance to a single
/// row, so expressing them as independent flags admitted combinations that
/// parsed, ran, and were structurally incapable of having an effect. A 425-
/// family run was lost to exactly that before this was folded into an enum.
fn drive_rows<A: cons_core::AlignmentSource>(
    a: A,
    cli: &Cli,
    seqs: &[aln_core::Sequence],
    params: &Params,
    ext: Option<&ExtendOpts>,
) -> anyhow::Result<Vec<cons_core::Refined>> {
    match cli.hsps {
        HspPolicy::Best => go_boxed(BestHsp(a), cli, seqs, params, ext),
        HspPolicy::All => go_boxed(a, cli, seqs, params, ext),
        // `TileVoteFilter` has already done the whole job upstream — trimming
        // and the strand vote both — so nothing wraps it here.
        HspPolicy::Tiled => go_boxed(a, cli, seqs, params, ext),
    }
}

/// `drive` with the generic argument erased, so the match arms above can each
/// hand it a differently-wrapped source.
fn go_boxed<A: cons_core::AlignmentSource>(
    a: A,
    cli: &Cli,
    seqs: &[aln_core::Sequence],
    params: &Params,
    ext: Option<&ExtendOpts>,
) -> anyhow::Result<Vec<cons_core::Refined>> {
    drive(a, cli, seqs, params, ext)
}

/// Write one Stockholm record: the alignment, its consensus as `#=GC RF`, and
/// the run's provenance as `#=GF` lines.
///
/// The body is produced by `dfam_stk_io::msa::write` — the same code that
/// backs the rest of the Dfam toolchain, so identifier suffixes, the `.`
/// gap convention and the single-block layout cannot drift from it here. That
/// writer takes no annotations, so the record is parsed straight back and the
/// provenance spliced in. The round trip is cheap next to a refinement pass and
/// it means this function cannot emit something the parser would reject.
fn write_stockholm<W: Write>(
    out: &mut W,
    r: &cons_core::Refined,
    params: &Params,
    comment: &str,
) -> anyhow::Result<()> {
    // `#=GC RF` needs the consensus in the alignment's own column coordinates,
    // not the ungapped form the caller reports.
    let call = ConsensusParams { include_reference: false, ..params.consensus.clone() };
    let gapped = match params.caller {
        cons_core::Caller::Dfam => r.msa.consensus(&call),
        cons_core::Caller::Giri => r.msa.giri_consensus(params.min_non_gap_count),
    };

    let mut buf: Vec<u8> = Vec::new();
    dfam_stk_io::msa::write(&r.msa, &mut buf, Some(&r.name), Some(&gapped), false)?;

    let mut rec = dfam_stk_io::iter_records(std::io::Cursor::new(&buf))
        .next()
        .ok_or_else(|| anyhow::anyhow!("Stockholm writer produced no record"))?
        .map_err(|e| anyhow::anyhow!("re-reading generated Stockholm: {e}"))?;

    // `msa::write` labels each row by appending its aligned extent to the
    // input name, which on a Smitten-named instance yields a recursive
    // identifier: `chr22:33101-33204_-:1-104_+`. That is well-formed but no
    // downstream tool should have to compose it, so collapse each one to the
    // single absolute range it denotes. Names that are not Smitten identifiers
    // (bare `gi|4`, say) are left exactly as the writer produced them.
    for row in rec.sequences.iter_mut() {
        if let Ok((id, _)) = smitten::Identifier::from_unknown_format(&row.original_id, false, false)
        {
            if id.ranges.len() > 1 {
                if let Ok(flat) = id.normalize() {
                    row.original_id = flat.to_string();
                }
            }
        }
    }

    // Provenance after ID, before SQ, so the record reads top-down.
    let at = rec.gf.iter().position(|(t, _)| t == "SQ").unwrap_or(rec.gf.len());
    let mut n = 0;
    for field in comment.split(' ') {
        if let Some((k, v)) = field.split_once('=') {
            rec.gf.insert(at + n, ("**".to_string(), format!("{k}: {v}")));
            n += 1;
        }
    }
    for note in &r.notes {
        rec.gf.insert(at + n, ("**".to_string(), note.clone()));
        n += 1;
    }
    rec.write_to(out)?;
    Ok(())
}

/// Print each instance at its core boundary, before the extension runs.
///
/// Ported from `printCoreEdges` in the C RAMExtend's `report.c`; the Rust port
/// had not carried it across. Flanks are ten bases in consensus orientation,
/// `*` marking one cut short by a boundary rather than by the window; the core
/// is 24 columns, centred when short and elided in the middle when not.
///
/// It is worth the space because it answers a question no summary statistic
/// does: if every instance shows the *same* flanking sequence, the family is a
/// segmental duplication rather than a transposable element, and the extension
/// is about to run away into shared context.
fn report_core_edges(cli: &Cli, edges: &[ram_core::library::CoreEdge]) {
    if cli.silent || edges.is_empty() {
        return;
    }
    eprintln!("  pre-extension core boundaries:");
    let idw = edges.iter().map(|e| e.identifier.len()).max().unwrap_or(5).max(5);
    eprintln!(
        "  {:>4} {:<idw$} {:<21} {:<6} {:<4}  {:>11} {:<26} {}",
        "Seq", "Ident", "Range", "Orient", "L/R?", "Left-Flank", "Core", "Right-Flank",
    );
    for e in edges {
        let blank = "";
        eprintln!(
            "  {:>4} {:<idw$} {:<21} {:<6} {}/{}   {:>11} [{:<24}] {}",
            e.index,
            e.identifier,
            // The C report prints the core 0-based fully closed; kept so
            // the two tools' output can be compared line for line.
            format!("{}-{}", e.span.start(), e.span.end() - 1),
            if e.minus { "-" } else { "+" },
            u8::from(e.left_extendable),
            u8::from(e.right_extendable),
            if e.left_extendable { e.left_flank.as_str() } else { blank },
            e.core,
            if e.right_extendable { e.right_flank.as_str() } else { blank },
        );
    }
}

// ── Run report ────────────────────────────────────────────────────────────────
//
// Everything here goes to stderr, so stdout stays a clean Stockholm stream.
// `--silent` turns the lot off.

/// The four stages a run passes through, named consistently everywhere.
const PHASE_BOOTSTRAP: &str = "consensus bootstrap";
const PHASE_REFINE: &str = "iterative refinement";
const PHASE_EXTEND: &str = "extension";
const PHASE_REPAIR: &str = "block repair";

fn banner(cli: &Cli) {
    if cli.silent {
        return;
    }
    eprintln!("te-composer {}", env!("CARGO_PKG_VERSION"));
    // The invocation verbatim: a report that does not say what produced it is
    // hard to trust six months later.
    let args: Vec<String> = std::env::args().skip(1).collect();
    eprintln!("  invocation: te-composer {}", args.join(" "));
}

/// What came in: how many instances, and how their lengths are spread.
///
/// The spread matters more than the count. A family of uniform full-length
/// instances and one of mostly fragments behave differently at every later
/// stage, and the median against the extremes says which you have.
fn report_input(cli: &Cli, seqs: &[aln_core::Sequence]) {
    if cli.silent || seqs.is_empty() {
        return;
    }
    let mut lens: Vec<usize> = seqs.iter().map(|s| s.len()).collect();
    lens.sort_unstable();
    let n = lens.len();
    let med = if n % 2 == 0 { (lens[n / 2 - 1] + lens[n / 2]) / 2 } else { lens[n / 2] };
    let total: usize = lens.iter().sum();
    // GC against the scoring matrix's own background. `comparison.matrix`
    // declares FREQS A .265 C .235 G .235 T .265, i.e. 47% GC, and the score
    // floor assumes it: family-49 at 69% GC produced 8.7M HSPs because
    // GC-rich sequence matches itself far above the background the floor was
    // set for. Worth seeing before a run rather than inferring afterwards.
    let (mut gc, mut acgt) = (0usize, 0usize);
    for sq in seqs {
        for &b in &sq.seq {
            match b.to_ascii_uppercase() {
                b'G' | b'C' => { gc += 1; acgt += 1; }
                b'A' | b'T' => acgt += 1,
                _ => {}
            }
        }
    }
    let gc_pct = if acgt > 0 { 100.0 * gc as f64 / acgt as f64 } else { 0.0 };
    let flag = if !(35.0..=60.0).contains(&gc_pct) { "  <- far from the 47% matrix background" } else { "" };
    eprintln!(
        "  input: {n} instances, {} bp total; lengths {}-{} bp (median {}, q1 {}, q3 {})",
        total, lens[0], lens[n - 1], med, lens[n / 4], lens[(3 * n) / 4].min(lens[n - 1]),
    );
    eprintln!("  composition: {:.0}% GC{}", gc_pct, flag);
}

fn phase(cli: &Cli, name: &str) {
    if !cli.silent {
        eprintln!("\n── {name} ──");
    }
}

/// Mean Kimura divergence of the instance rows against the consensus, plain and
/// CpG-adjusted, plus the consensus length.
///
/// CpG sites mutate fast enough to dominate a raw divergence estimate, so the
/// adjusted figure is the one to compare families on; both are reported because
/// their gap is itself informative about CpG content.
/// Returns `(divergence, divergence_cpg, consensus_len, rows, instances)`.
///
/// `rows` and `instances` differ under `--hsps tiled`, which can give one
/// instance several rows. Row names are the instance names, unmodified by MSA
/// assembly, so distinct names count distinct contributors.
fn family_stats(
    msa: &aln_core::msa::MultiAlign,
    gapped_cons: &[u8],
) -> (f64, f64, usize, usize, usize) {
    let (plain, cpg) = cons_core::mean_kimura(msa, gapped_cons);
    let ungapped = gapped_cons.iter().filter(|&&b| b != b'-' && b != b' ').count();
    let instances = msa
        .sequences
        .iter()
        .skip(1) // row 0 is the reference
        .map(|r| r.name.as_str())
        .collect::<std::collections::HashSet<_>>()
        .len();
    (plain, cpg, ungapped, msa.num_instances(), instances)
}

/// Overlap threshold for the `--hsps all|best` cull, matching `Refiner`'s
/// `setMaskLevel(80)`. The default `tiled` policy has no cull — trimming
/// subsumes it.
const MASK_LEVEL: u32 = 80;

/// `N rows (X.XX per instance)`, or just `N rows` when they are one to one.
fn rows_note(rows: usize, instances: usize) -> String {
    format!("{} rows{}", rows, rows_ratio(rows, instances))
}

/// The bare count plus the ratio, for a line that already says "rows:".
fn rows_count(rows: usize, instances: usize) -> String {
    format!("{}{}", rows, rows_ratio(rows, instances))
}

fn rows_ratio(rows: usize, instances: usize) -> String {
    if instances == 0 || rows == instances {
        String::new()
    } else {
        format!(" ({:.2} per instance)", rows as f64 / instances as f64)
    }
}

/// The bootstrap consensus as pass 1 measured it, before any refinement.
fn report_bootstrap_state(cli: &Cli, r: &cons_core::Refined, total: usize) {
    if cli.silent {
        return;
    }
    let Some(t) = r.trace.first() else { return };
    eprintln!(
        "  consensus {} bp, score {}, {} of {total} instances, {}, \
Kimura {:.2}% (CpG-adjusted {:.2}%)",
        t.consensus.len(),
        t.score,
        t.instances,
        rows_note(t.rows, t.instances),
        t.divergence,
        t.divergence_cpg,
    );
}

/// The gapped consensus for a finished refinement, in the MSA's own columns.
fn gapped_consensus(r: &cons_core::Refined, params: &Params) -> Vec<u8> {
    let call = ConsensusParams { include_reference: false, ..params.consensus.clone() };
    match params.caller {
        cons_core::Caller::Dfam => r.msa.consensus(&call),
        cons_core::Caller::Giri => r.msa.giri_consensus(params.min_non_gap_count),
    }
}

/// One line of family state, printed at the end of every stage so the stages
/// can be compared directly.
fn report_state(
    cli: &Cli,
    r: &cons_core::Refined,
    params: &Params,
    total_instances: usize,
    stage: Option<&str>,
) {
    if cli.silent {
        return;
    }
    let g = gapped_consensus(r, params);
    let (div, divcpg, len, rows, instances) = family_stats(&r.msa, &g);
    // How many copies the consensus actually rests on. A family where a third
    // of the input never aligns is a different object from one where all of it
    // does, and the bare row count does not say which you have.
    let head = match stage {
        Some(t) => format!("{t}: "),
        None => String::new(),
    };
    eprintln!(
        "  {head}consensus {len} bp, score {}, {instances} of {total_instances} instances, \
{}, Kimura {div:.2}% (CpG-adjusted {divcpg:.2}%)",
        r.score,
        rows_note(rows, instances),
    );
}

/// The finished family, restated in one place.
///
/// Redundant with the per-stage lines by design: a reader who skipped the
/// middle should not have to reconstruct where it ended up.
fn report_summary(
    cli: &Cli,
    r: &cons_core::Refined,
    params: &Params,
    total_instances: usize,
    elapsed: std::time::Duration,
) {
    if cli.silent {
        return;
    }
    let g = gapped_consensus(r, params);
    let (div, divcpg, len, rows, instances) = family_stats(&r.msa, &g);
    eprintln!("\n── summary ──");
    eprintln!("  {}: consensus {len} bp", r.name);
    eprintln!("  instances:  {instances} of {total_instances} participating");
    eprintln!("  rows:       {}", rows_count(rows, instances));
    eprintln!("  divergence: Kimura {div:.2}%, CpG-adjusted {divcpg:.2}%");
    eprintln!("  score:      {}", r.score);
    eprintln!("  runtime:    {:.1}s", elapsed.as_secs_f64());
}

/// Which sequence seeded the consensus, and by how much it won.
fn report_reference(cli: &Cli, r: &cons_core::Refined, seqs: &[aln_core::Sequence]) {
    if cli.silent {
        return;
    }
    let f = &r.reference;
    let len = seqs.get(f.index).map(|s| s.len()).unwrap_or(0);
    eprint!("  reference: {} ({} bp, index {}), score {}", f.name, len, f.index, f.score);
    match f.runner_up {
        Some(next) if f.score > 0 => {
            let margin = f.score - next;
            eprintln!(
                ", {:.1}% ahead of the runner-up ({next}) of {} candidates",
                100.0 * margin as f64 / f.score as f64,
                f.candidates
            );
        }
        Some(next) => eprintln!(", runner-up {next} of {} candidates", f.candidates),
        None => eprintln!(" (only candidate)"),
    }
}

/// Report how a phase-2 loop terminated.
///
/// A run contains up to three of these — the initial refinement, one after a
/// genome extension, and one after the block repair — and they behave quite
/// differently. Printing the pass count and the stop reason for each makes the
/// iteration budget visible: `exhausted` means the consensus was still moving
/// when the budget ran out, which is the case worth noticing.
fn report_passes(cli: &Cli, label: Option<&str>, r: &cons_core::Refined, budget: usize) {
    if cli.silent {
        return;
    }
    let why = match r.stop {
        cons_core::StopReason::Converged => "converged",
        cons_core::StopReason::Cycled => "cycled",
        cons_core::StopReason::Exhausted => "EXHAUSTED — still changing at the limit",
    };
    let plural = if r.passes == 1 { "pass" } else { "passes" };
    let head = match label {
        Some(t) => format!("{t}: "),
        None => String::new(),
    };
    eprintln!("  {head}{} {plural} of {budget}, {why}", r.passes);
}

/// The genome, the extension settings, and whether to refine afterwards.
///
/// Held together because the 2bit is opened once and read by every family.
struct ExtendOpts {
    genome: aln_core::twobit::TwoBitReader,
    cfg: extend::Config,
    refine: bool,
}

fn drive<A: cons_core::AlignmentSource>(
    a: A,
    cli: &Cli,
    seqs: &[aln_core::Sequence],
    params: &Params,
    ext: Option<&ExtendOpts>,
) -> anyhow::Result<Vec<cons_core::Refined>> {
    let a = MinScore::new(a, 0);
    match &cli.dump_alignments {
        Some(path) => {
            let f = std::fs::File::create(path)
                .with_context(|| format!("creating {}", path.display()))?;
            let census = Census::new(a, Box::new(std::io::BufWriter::new(f)))
                .context("starting the alignment census")?;
            compose(&census, cli, seqs, params, ext)
        }
        None => compose(&a, cli, seqs, params, ext),
    }
}

/// Run the pipeline, then extend each consensus against the genome if asked.
/// Run the pipeline, extend if asked, and repair blocks **last**.
///
/// Block repair is an evaluate-and-maybe-keep step: it patches the consensus,
/// re-refines, and keeps the result only if the score improved. That makes it
/// the wrong thing to run in the middle of a pipeline that goes on to modify
/// the consensus again — `Refiner` repairs before extending, so the bases
/// RAMExtend adds are never offered to the repair at all, and the verdict is
/// reached against a consensus that no longer exists by the time the run ends.
///
/// Here it runs after every other modification, against the final instance set
/// (widened, when the family was extended). The loop count is unchanged: the
/// refinement the repair triggers replaces the one that used to follow it.
fn compose<A: cons_core::AlignmentSource>(
    a: &A,
    cli: &Cli,
    seqs: &[aln_core::Sequence],
    params: &Params,
    ext: Option<&ExtendOpts>,
) -> anyhow::Result<Vec<cons_core::Refined>> {
    // The main loop runs without repair; it is reinstated at the end.
    let mut deferred = params.clone();
    deferred.repair_blocks = false;

    let supplied = cli.bootstrap_cons.is_some();
    let results = match &cli.bootstrap_cons {
        // Bootstrap exists only to invent a starting consensus. Given one,
        // there is nothing for it to do.
        Some(path) => {
            let recs = io::read_fasta_file(path)
                .with_context(|| format!("reading bootstrap consensus {}", path.display()))?;
            let seed = recs
                .into_iter()
                .next()
                .ok_or_else(|| anyhow::anyhow!("{} contains no sequence", path.display()))?;
            phase(cli, PHASE_BOOTSTRAP);
            if !cli.silent {
                eprintln!(
                    "  skipped — starting from supplied consensus {} ({} bp)",
                    seed.name,
                    seed.seq.len()
                );
            }
            // No candidate was chosen, so there is no index, score or margin to
            // report; `candidates: 0` is what marks the reference as supplied.
            let reference = cons_core::Reference {
                index: usize::MAX,
                name: seed.name.clone(),
                score: 0,
                runner_up: None,
                candidates: 0,
            };
            match cons_core::refine(a, seqs, &seed.seq, &params.base_name, reference, &deferred)? {
                Some(r) => vec![r],
                None => Vec::new(),
            }
        }
        None => {
            phase(cli, PHASE_BOOTSTRAP);
            run(a, seqs, &deferred)?
        }
    };
    let mut out = Vec::with_capacity(results.len());

    for mut refined in results {
        if !supplied {
            report_reference(cli, &refined, seqs);
            report_bootstrap_state(cli, &refined, seqs.len());
        }
        phase(cli, PHASE_REFINE);
        report_passes(cli, None, &refined, params.iterations);
        report_state(cli, &refined, params, seqs.len(), None);
        let widened = match ext {
            Some(opts) => {
                phase(cli, PHASE_EXTEND);
                let w = extend_one(a, cli, seqs, &deferred, opts, &mut refined)?;
                if w.is_some() {
                    // Probe the *extended consensus*, which is where a
                    // satellite becomes visible. The pre-search probe cannot
                    // see this case: `hs1/rnd-5_family-3330` has 71 bp
                    // instances, far too short to self-align into anything, and
                    // extension grew them into 15 kb that is 50% `(GGAAT)n` —
                    // HSATII. Measured on the extended consensus it scores 55
                    // self-HSPs against 4 and 6 for genuine 18 kb and 28 kb
                    // extensions, so the same threshold separates them.
                    //
                    // This only ever warns. The consensus is already built, the
                    // curator is the one who can judge whether a satellite
                    // family is wanted, and refusing here would discard work
                    // rather than avoid it.
                    phase(cli, PHASE_REFINE);
                    report_passes(cli, None, &refined, params.iterations);
                    report_state(cli, &refined, params, seqs.len(), None);
                }
                w
            }
            None => None,
        };
        let final_seqs: &[aln_core::Sequence] = widened.as_deref().unwrap_or(seqs);

        out.push(match (&params.repair_matrix, params.repair_blocks) {
            (Some(mx), true) => {
                phase(cli, PHASE_REPAIR);
                let r = repair_last(a, cli, final_seqs, params, mx, refined)?;
                report_state(cli, &r, params, seqs.len(), None);
                r
            }
            _ => refined,
        });
    }
    Ok(out)
}

/// The block repair, applied to a finished consensus.
///
/// Mirrors `cons_core::refine_with_repair`'s second half: derive the gapped
/// consensus from the MSA that produced it, select blocks, patch, re-refine,
/// and let the accept rule choose between the two fully-refined candidates.
fn repair_last<A: cons_core::AlignmentSource>(
    a: &A,
    cli: &Cli,
    seqs: &[aln_core::Sequence],
    params: &Params,
    matrix: &aln_core::SubstMatrix,
    refined: cons_core::Refined,
) -> anyhow::Result<cons_core::Refined> {
    let call = ConsensusParams { include_reference: false, ..params.consensus.clone() };
    let gapped = match params.caller {
        cons_core::Caller::Dfam => refined.msa.consensus(&call),
        cons_core::Caller::Giri => refined.msa.giri_consensus(params.min_non_gap_count),
    };
    let fixes = cons_core::default_block_source(&refined.msa, matrix, params, &call);
    if fixes.is_empty() {
        if !cli.silent {
            eprintln!("  no low-quality blocks selected");
        }
        return Ok(refined);
    }
    let reference = cons_core::Reference {
        index: refined.reference.index,
        name: refined.reference.name.clone(),
        score: refined.reference.score,
        runner_up: refined.reference.runner_up,
        candidates: refined.reference.candidates,
    };
    let name = refined.name.clone();
    // `apply_repair` may return a freshly refined Refined, which starts with no
    // provenance; carry what the family has accumulated so far across it.
    let carried = refined.notes.clone();
    let outcome =
        cons_core::apply_repair(a, seqs, refined, &gapped, &name, reference, params, &fixes)?;
    if let cons_core::RepairOutcome::Judged { second, report, .. } = &outcome {
        report_passes(cli, None, second, params.iterations);
        let kept = match params.repair_accept {
            cons_core::AcceptRule::Total => report.keep_total,
            cons_core::AcceptRule::MeanPerBase => report.keep_mean,
            cons_core::AcceptRule::Always => true,
        };
        if !cli.silent {
        eprintln!(
            "  {} block{} patched, score {} -> {} — {}",
            report.fixes,
            if report.fixes == 1 { "" } else { "s" },
            report.first_score,
            report.second_score,
            if kept { "kept" } else { "discarded, keeping the unrepaired consensus" },
        );
        }
    }
    let mut kept = outcome.take(params.repair_accept);
    if kept.notes.is_empty() {
        kept.notes = carried;
    }
    Ok(kept)
}

/// Extend one finished consensus past its edges, then refine what came back.
///
/// Splicing the extension on and stopping would leave the new bases as
/// RAMExtend called them, never having been seen by the consensus caller that
/// produced the rest. Refining again puts the whole sequence — old and new —
/// through one code path, which is also what the Perl `Refiner` did with the
/// extension it got back.
fn extend_one<A: cons_core::AlignmentSource>(
    a: &A,
    cli: &Cli,
    seqs: &[aln_core::Sequence],
    params: &Params,
    opts: &ExtendOpts,
    refined: &mut cons_core::Refined,
) -> anyhow::Result<Option<Vec<aln_core::Sequence>>> {
    let (outcome, warnings) = extend::extend(&opts.genome, &opts.cfg, &refined.msa, seqs)
        .with_context(|| format!("extending {}", refined.name))?;
    if !cli.silent {
        for w in &warnings {
            eprintln!("  warning: {w}");
        }
    }
    // Whatever the extension decided, the family carries the record of it.
    // A refusal is the interesting case: the consensus looks exactly like an
    // unextended one, so without this there is nothing to say it was
    // considered and declined.
    for w in &warnings {
        refined.notes.push(format!("Extension: {w}"));
    }
    let Some(outcome) = outcome else {
        refined
            .notes
            .push("Extension: not applied".to_string());
        return Ok(None);
    };
    // Pre-extension: what each copy looks like at its core boundary, and
    // therefore what the extension is about to reason over. Printed before
    // the result so it reads as input, not outcome.
    refined.notes.push(format!(
        "Extension: anchored {} of {} copies, divergence {:.1}%, matrix {}",
        outcome.contributors,
        seqs.len(),
        outcome.divergence,
        outcome.matrix_name,
    ));
    if outcome.extended() {
        // Say when only one side was taken. A curator reading `left 0 bp,
        // right 4210 bp` cannot tell a side that found no extension from a
        // side that was dropped for running away, and those mean opposite
        // things about how far to trust the edge.
        let note = match (outcome.left.len(), outcome.right.len()) {
            (0, r) => format!("Extension: right {r} bp only (left side dropped or empty)"),
            (l, 0) => format!("Extension: left {l} bp only (right side dropped or empty)"),
            (l, r) => format!("Extension: left {l} bp, right {r} bp"),
        };
        refined.notes.push(note);
    }
    report_core_edges(cli, &outcome.edges);
    if !cli.silent {
    eprintln!(
        "  anchored {} copies, divergence {:.1}%, matrix {}; \
extended left {} bp, right {} bp",
        outcome.contributors,
        outcome.divergence,
        outcome.matrix_name,
        outcome.left.len(),
        outcome.right.len(),
    );
    }
    if !outcome.extended() {
        return Ok(None);
    }

    let start = outcome.splice(&refined.consensus);
    if !opts.refine {
        refined.consensus = start;
        return Ok(Some(outcome.seqs));
    }
    // Reference selection belongs to phase 1 and is not revisited; the record
    // is carried forward so the report still names the sequence that seeded
    // this consensus.
    let reference = cons_core::Reference {
        index: refined.reference.index,
        name: refined.reference.name.clone(),
        score: refined.reference.score,
        runner_up: refined.reference.runner_up,
        candidates: refined.reference.candidates,
    };
    match cons_core::refine(a, &outcome.seqs, &start, &refined.name, reference, params)? {
        Some(mut again) => {
            // `refine` returns a fresh Refined; the provenance accumulated so
            // far would otherwise be dropped here.
            again.notes = std::mem::take(&mut refined.notes);
            *refined = again;
        }
        None => {
            if !cli.silent {
                eprintln!("  refinement after extension produced nothing; keeping the spliced consensus");
            }
            refined.consensus = start;
        }
    }
    // The widened instances become the input for everything downstream — the
    // final block repair has to align against the same sequences the
    // consensus was just refined from, not the un-extended originals.
    Ok(Some(outcome.seqs))
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    banner(&cli);
    let started = std::time::Instant::now();
    let backend = cli.backend;

    // Sequences are read before the thread pool is built: the pool size depends
    // on how long and how many they are, and reading does no parallel work.
    let seqs = if cli.input == "-" {
        io::read_fasta(std::io::BufReader::new(std::io::stdin().lock()))
            .context("reading FASTA from stdin")?
    } else {
        io::read_fasta_file(&cli.input)
            .with_context(|| format!("reading {}", cli.input))?
    };
    if seqs.is_empty() {
        bail!("no sequences read from {}", cli.input);
    }

    let threads = match backend {
        Backend::Parasail | Backend::Reference => plan_parasail(&seqs, &cli)?,
        // The seeded search is not memory-planned; it takes what it is given.
        Backend::Rmblast => {
            if cli.threads > 0 { cli.threads } else { num_cpus() }
        }
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .context("configuring the thread pool")?;

    let matrix = match &cli.matrix {
        Some(p) => SubstMatrix::from_file(p)
            .with_context(|| format!("reading matrix {}", p.display()))?,
        None => SubstMatrix::parse(DEFAULT_MATRIX).expect("built-in matrix must parse"),
    };

    // Penalties come from the matrix's GAP line unless overridden.  Pairing a
    // matrix with penalties on a different scale silently wrecks alignments —
    // a max-match-3 matrix with a gap cost of 30 can never extend through an
    // indel — so the default keeps them together.
    let from_matrix = AlignParams::from_matrix(&matrix);
    // GIRI applies no score cutoff at all: `MultipleAlignment::align` hardwires
    // `minScore = 0` (the adaptive `0.3 * mean` is commented out) and neither SW
    // aligner has a floor.  So 1 is the comparable default *for a DP backend*,
    // which reports one optimal alignment per pair and pays nothing for a low
    // floor.
    //
    // A seeded search is different in kind, and now that it is the default
    // backend the difference matters.  rmblast reports every extension clearing
    // the floor *and* derives its X-drop cutoffs from it
    // (`min_score*2 / min_score/2 / min_score`), so a floor of 1 collapses the
    // extension budget to 2/0/1 and turns gapped extension off, while returning
    // roughly a hundred junk HSPs per pair.  Measured on 25 x 8 kb: 11.0 s and
    // 2.7 GB at a floor of 1, against 1.5 s and 88 MB at 150.
    //
    // 150 is Refiner's own `setMinScore(150)`.  The help text has claimed this
    // default for some time; it was never implemented.
    // 150 for every backend now — `Refiner`'s own `setMinScore(150)`, in units
    // of the default matrix, which both backends share.
    //
    // The DP backends previously defaulted to 1, inherited from the C++, where
    // `MultipleAlignment::align` hardwires `minScore = 0` with GIRI's adaptive
    // `0.3 * mean` commented out just above it. That floor admits a great deal
    // of noise: against a dinucleotide-shuffled null, ~68% of
    // reference-selection alignments score at or below the noise ceiling, and
    // the shortest alignment observed was two matched bases.
    //
    // Reproducing the C++ is `autocons`' job, not this one.
    let min_score = cli.min_score.unwrap_or(150);
    let align_params = AlignParams {
        mode: AlignMode::Local,
        gap_open: cli.gap_open.unwrap_or(from_matrix.gap_open),
        gap_extend: cli.gap_extend.unwrap_or(from_matrix.gap_extend),
        min_score,
        traceback: true,
        bandwidth: None,
    };
    if align_params.gap_open < 3 * align_params.gap_extend.max(1) {
        eprintln!(
            "warning: gap_open {} is small relative to gap_extend {} — check the \
             penalties match the matrix's scale",
            align_params.gap_open, align_params.gap_extend
        );
    }

    // A pass costs a full re-align of every instance, so the sensible budget
    // tracks the backend's per-pass cost: the C++'s 3 for DP, 11 for the
    // seeded search. Cycle detection makes the larger value safe.
    //
    // Long families once got half that budget, on the reasoning that a pass is
    // quadratic-ish in length. Removed after measuring it. Across 4,823
    // families it fired on one, mm39/rnd-1_family-93 (median instance 16,046
    // bp), whose three rounds converged in 2, 3 and 2 passes and whose
    // consensus is byte-identical either way.
    //
    // It does bind on family-49 (608,515 bp of input, median instance 11,893
    // bp), and there it only makes the stop worse: capped, round 1 is cut off
    // at 10 passes still changing; uncapped it reaches its own cycle detector
    // at 14. Both land on the same 12,154 bp consensus at score 5,059,173,
    // because block repair converges them. Cycle detection is the real bound.
    let iterations = cli.iterations.unwrap_or(match backend {
        Backend::Rmblast => 20,
        Backend::Parasail | Backend::Reference => 3,
    });

    let params = Params {
        num_consensi: cli.num,
        base_name: cli.name.clone(),
        iterations,
        consensus: ConsensusParams {
            enable_cpg: !cli.no_cpg,
            ..Default::default()
        },
        caller: if cli.orig { Caller::Giri } else { Caller::Dfam },
        min_non_gap_count: cli.min,
        restore_cpg: cli.mam,
        repair_blocks: !cli.no_repair_blocks,
        repair_all_vs_all: true,
        pack_insertions: !cli.no_pack_insertions,
        pack_max_sep: cli.pack_max_sep,
        pack_min_seg: cli.pack_min_seg,
        pack_min_score: cli.pack_min_score,
        pack_keep_insertions: true,
        repair_threshold: cli.repair_threshold,
        repair_window: cli.repair_window,
        nonredundant_reference_score: false,
        repair_accept: cons_core::AcceptRule::default(),
        // The repair scores columns with the same matrix the aligner uses, so a
        // run cannot end up judging alignment quality on a different scale from
        // the one that produced the alignment.
        // Packing needs the matrix too, and asks for it independently of the
        // block repair. Gating it on `--repair-blocks` alone made
        // `--pack-insertions --no-repair-blocks` a silent no-op: the packing
        // code checked for a matrix that had just been set to `None`.
        repair_matrix: (!cli.no_repair_blocks || !cli.no_pack_insertions)
            .then(|| matrix.clone()),
        insertions: cli.insertions.into(),
    };

    // The genome is opened once, before any family runs: a bad path or an
    // unreadable 2bit should fail here rather than after the whole pipeline
    // has been paid for.
    let ext = match &cli.genome {
        Some(path) => {
            let genome = aln_core::twobit::TwoBitReader::open(path)
                .with_context(|| format!("opening genome {}", path.display()))?;
            Some(ExtendOpts {
                genome,
                cfg: extend::Config {
                    assembly: cli.assembly.clone(),
                    l_max: cli.extend_max,
                    bandwidth: cli.extend_bandwidth,
                    min_aligning_seqs: cli.extend_min_seqs,
                    matrix: cli.extend_matrix.clone(),
                    min_improvement: cli.extend_min_improvement,
                    edge_slop: cli.extend_edge_slop,
                    min_identity: cli.extend_min_identity,
                    max_total: cli.extend_max_total,
                    ..extend::Config::default()
                },
                refine: !cli.extend_no_refine,
            })
        }
        None => None,
    };

    report_input(&cli, &seqs);

    let mut probe_engine: Option<RmblastEngine> = None;
    let mut results = match backend {
        Backend::Parasail => {
            #[cfg(not(any(target_arch = "x86_64", target_arch = "x86")))]
            if !cli.silent {
                eprintln!(
                    "note: parasail's SIMD kernels are x86-only; \
                     --backend parasail runs the scalar reference aligner here"
                );
            }
            let a = FastAligner::new(matrix, align_params)
                .context("building the parasail aligner")?;
            drive(Pairwise::new(a), &cli, &seqs, &params, ext.as_ref())
        }
        Backend::Reference => {
            let a = ReferenceAligner::new(matrix, align_params)
                .context("building the reference aligner")?;
            drive(Pairwise::new(a), &cli, &seqs, &params, ext.as_ref())
        }
        Backend::Rmblast => {
            // rmblast takes signed penalties and does its own NCBI gap-cost
            // conversion internally; hand it the same magnitudes, negated.
            let search = SearchParams {
                matrix: Some(matrix),
                gap_init: -(align_params.gap_open as i32),
                ins_gap_ext: -(align_params.gap_extend as i32),
                del_gap_ext: -(align_params.gap_extend as i32),
                min_match: cli.word_size,
                min_score: align_params.min_score.max(0),
                // Engine-level masking is retired: overlap filtering happens
                // in cons-core's RefinerFilter below, which applies Refiner's
                // post-alignment semantics (reference-axis novelty) to BOTH
                // phases — including reference selection, which the engine's
                // own per-pair query-axis mask never covered. 101 = off.
                mask_level: 101,
                cores: Some(threads.max(1)),
                // `Refiner` sets this; GIRI does not.  Off by default so the
                // aligner comparison varies one thing at a time.
                score_mode: if cli.complexity_adjust {
                    ScoreMode::ComplexityAdjusted
                } else {
                    ScoreMode::Basic
                },
                ..Default::default()
            };
            let mut opts = RmblastOptions::default();
            opts.dust = !cli.no_dust;
            // Pre-flight: is the input internally repetitive? Run before the
            // all-vs-all, because that search is what explodes. DUST is
            // deliberately off here whatever the shipping default is — it masks
            // low-complexity query sequence, so with it on a pure tandem repeat
            // self-aligns to nothing and the probe reports 0.
            {
                // Kept for the post-extension check too — same requirements
                // (raw hits, DUST off), so building it twice would be waste.
                let pe = RmblastEngine::new(
                    search.clone(),
                    RmblastOptions { dust: false, ..RmblastOptions::default() },
                )
                .context("building the repeat-check engine")?;
                let pcfg = cons_core::probe::ProbeConfig::default();
                let rep = cons_core::probe::probe_repetitive(&pe, &seqs, &pcfg)
                    .context("probing the input for internal repeats")?;
                if !rep.flagged.is_empty() {
                    let worst = &rep.flagged[0];
                    let detail = format!(
                        "{} of {} probed instances are internally repetitive \
(worst: {} with {} self-alignments in a {} bp window, limit {})",
                        rep.flagged.len(),
                        rep.probed.len(),
                        worst.name,
                        worst.hsps,
                        pcfg.window,
                        pcfg.max_hsps,
                    );
                    if !cli.silent {
                        eprintln!("  WARNING: {detail}");
                        eprintln!("           if these are a satellite or tandem array, the consensus \
will be a tiling of the unit rather than one element");
                    }
                }
                probe_engine = Some(pe);
            }

            let a = RmblastEngine::new(search, opts)
                .context("building the rmblast engine")?;
            // Refiner's overlap filter, then (optionally) best-per-instance.
            // Order matters only in principle — the chain's best member always
            // survives the filter, so Best sees the same winner either way —
            // but filtering first means phase-1 candidate sums are deduped in
            // both modes, which is the point.
            // The orientation rule runs in both phases; there is no switch.
            let (p1, p2) = (true, true);
            if cli.hsps == HspPolicy::Tiled {
                // Tile-then-vote replaces the cull, the orientation rule and the
                // trim in one pass, so it wraps the raw engine and `RefinerFilter`
                // is skipped entirely.
                let a = cons_core::mhsp::TileVoteFilter::new(a, cli.min_row_len, p1, p2);
                drive_rows(a, &cli, &seqs, &params, ext.as_ref())
            } else {
                let a = RefinerFilter::with_phase_options(a, MASK_LEVEL, false, p1, p2);
                drive_rows(a, &cli, &seqs, &params, ext.as_ref())
            }
        }
    }
    .context("building consensus sequences")?;

    // A satellite the input was too short to reveal only becomes visible after
    // extension: `hs1/rnd-5_family-3330` has 71 bp instances — far too short to
    // self-align into anything — and extension grew them into 15 kb that is 50%
    // `(GGAAT)n`, HSATII. Measured on the extended consensus it scores 55
    // self-alignments against 4 and 6 for genuine 18 kb and 28 kb extensions,
    // so the input threshold separates them unchanged.
    //
    // This only ever warns, whatever `--repeat-check` says. The consensus is
    // already built; refusing here would discard work rather than avoid it, and
    // whether a satellite family is wanted is the curator's call.
    if ext.is_some() {
        if let Some(pe) = probe_engine.as_ref() {
            let cfg = cons_core::probe::ProbeConfig::default();
            for r in results.iter_mut() {
                let n = cons_core::probe::self_hsp_count(pe, &r.consensus, &cfg)
                    .unwrap_or(0);
                if n > cfg.max_hsps {
                    let msg = format!(
                        "the extended consensus is internally repetitive \
({n} self-alignments in a {} bp window, limit {}) — likely a satellite or \
tandem array rather than a transposable element",
                        cfg.window, cfg.max_hsps
                    );
                    if !cli.silent {
                        eprintln!("  WARNING: {msg}");
                    }
                    r.notes.push(msg);
                }
            }
        }
    }

    if results.is_empty() {
        bail!(
            "no consensus could be built from {} sequence(s) — nothing aligned above \
             the minimum score of {}",
            seqs.len(),
            min_score
        );
    }

    let mut out: Box<dyn Write> = match &cli.output {
        Some(p) => Box::new(std::io::BufWriter::new(
            std::fs::File::create(p).with_context(|| format!("creating {}", p.display()))?,
        )),
        None => Box::new(std::io::BufWriter::new(std::io::stdout().lock())),
    };

    if let Some(path) = &cli.dump_refinement {
        let mut f = std::io::BufWriter::new(
            std::fs::File::create(path)
                .with_context(|| format!("creating {}", path.display()))?,
        );
        for r in &results {
            for t in &r.trace {
                io::write_fasta(
                    &mut f,
                    &format!("{}_pass{}", r.name, t.pass),
                    Some(&format!("SCORE={} PASS={}", t.score, t.pass)),
                    &t.consensus,
                    cli.width,
                )?;
            }
        }
    }

    let mut cons_out: Option<Box<dyn Write>> = match &cli.consensus {
        Some(p) => Some(Box::new(std::io::BufWriter::new(
            std::fs::File::create(p).with_context(|| format!("creating {}", p.display()))?,
        ))),
        None => None,
    };

    for r in &results {
        // `SCORE=` matches the C++.  `REF=` and `REFSCORE=` are additions:
        // reference selection is a discrete argmax, so knowing which input
        // seeded a consensus — and by what margin — is the first thing anyone
        // debugging an unexpected result needs.
        let stop = match r.stop {
            cons_core::StopReason::Converged => "converged",
            cons_core::StopReason::Cycled => "cycled",
            cons_core::StopReason::Exhausted => "exhausted",
        };
        let comment = format!(
            "SCORE={:.2} REF={} REFIDX={} REFSCORE={} PASSES={} CONVERGED={} STOP={}",
            r.score as f64, r.reference.name, r.reference.index, r.reference.score,
            r.passes, r.converged, stop,
        );

        write_stockholm(&mut out, r, &params, &comment)
            .with_context(|| format!("writing Stockholm for {}", r.name))?;

        if let Some(w) = cons_out.as_mut() {
            match cli.format {
                OutFormat::Fasta => {
                    io::write_fasta(w, &r.name, Some(&comment), &r.consensus, cli.width)?
                }
                OutFormat::Ig => {
                    io::write_ig(w, &r.name, Some(&comment), &r.consensus, cli.width)?
                }
            }
        }

        if let Some(prefix) = &cli.aln {
            let path = format!("{prefix}.{}", r.name);
            let mut f = std::io::BufWriter::new(
                std::fs::File::create(&path)
                    .with_context(|| format!("creating {path}"))?,
            );
            // Rows are already full width, padded with spaces outside their own
            // aligned span (`row.start`..=`row.end`). Write them verbatim: the
            // padding is what distinguishes "this copy does not reach here" from
            // "this copy has a deletion here", and a reader that collapses the
            // two will count absent copies as evidence of deletion.
            for (i, row) in r.msa.sequences.iter().enumerate() {
                let tag = (i == 0).then_some("reference");
                io::write_fasta(&mut f, &row.name, tag, &row.seq, cli.width)?;
            }
        }
    }
    if let Some(mut w) = cons_out {
        w.flush()?;
    }
    for r in &results {
        report_summary(&cli, r, &params, seqs.len(), started.elapsed());
    }
    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_builtin_matrix_parses() {
        let m = SubstMatrix::parse(DEFAULT_MATRIX).unwrap();
        // comparison.matrix carries 11 symbols (A R G C Y T K M S W N).
        assert_eq!(m.size(), 11);
        assert!(m.lambda().is_some());
        // The scale matters more than the size: thresholds like --min-score 150
        // are in these units, so a silently rescaled matrix would silently
        // change what they mean. 25p43g's mean diagonal is 8.5 — the same
        // scale family as 25p43g's 8.5, NOT xrepmask/3's 3.0.
        for (sym, want) in [(b'A', 9i32), (b'C', 10), (b'G', 10), (b'T', 9)] {
            assert_eq!(m.score(sym, sym), Some(want), "diagonal for {}", sym as char);
        }
    }

    #[test]
    fn cli_parses_a_minimal_invocation() {
        let cli = Cli::try_parse_from(["autocons", "in.fa"]).unwrap();
        assert_eq!(cli.input, "in.fa");
        assert_eq!(cli.num, 1.0);
        assert_eq!(cli.backend, Backend::Rmblast);
        assert_eq!(cli.format, OutFormat::Ig);
    }

    #[test]
    fn cli_accepts_a_fractional_count() {
        let cli = Cli::try_parse_from(["autocons", "in.fa", "-n", "0.25"]).unwrap();
        assert_eq!(cli.num, 0.25);
    }

    #[test]
    fn scratch_estimate_is_quadratic_in_length() {
        let ten_kb = alignment_scratch_bytes(10_000, 10_000);
        // 4 bytes/cell: the 32-bit-lane branch, which is what the budget always
        // assumes.  10 kb x 10 kb is 1e8 cells, so ~381 MiB.
        let mb = ten_kb as f64 / (1024.0 * 1024.0);
        assert!(
            (370.0..390.0).contains(&mb),
            "estimate for 10kb x 10kb is {mb:.0} MB, not the 4 bytes/cell the \
             32-bit lane measurements pin it to"
        );
        assert_eq!(alignment_scratch_bytes(20_000, 20_000), ten_kb * 4);
    }

    #[test]
    fn scratch_estimate_does_not_overflow_on_absurd_input() {
        let huge = alignment_scratch_bytes(usize::MAX, usize::MAX);
        assert!(huge > 0, "saturating arithmetic should not wrap to zero");
    }

    const GB: u64 = 1024 * 1024 * 1024;

    #[test]
    fn concurrency_is_capped_by_sequence_count() {
        // A 25-sequence family cannot use 64 threads however much memory exists:
        // phase 1 parallelises over candidate references.
        assert_eq!(effective_threads(64, 25), 25);
        assert_eq!(effective_threads(8, 25), 8);
        assert_eq!(effective_threads(64, 0), 1, "never returns zero");
    }

    #[test]
    fn memory_budget_parses_sizes_and_percentages() {
        assert_eq!(parse_memory_budget("48G", None).unwrap(), 48 * GB);
        assert_eq!(parse_memory_budget("512M", None).unwrap(), 512 * 1024 * 1024);
        assert_eq!(parse_memory_budget("80%", Some(100 * GB)).unwrap(), 80 * GB);
        // A percentage with no MemAvailable to measure against must say so
        // rather than silently pick a number.
        assert!(parse_memory_budget("80%", None).is_err());
        assert!(parse_memory_budget("200%", Some(GB)).is_err());
        assert!(parse_memory_budget("banana", None).is_err());
    }

    #[test]
    fn a_comfortable_run_uses_every_thread_it_can() {
        // 2 kb x 100 seqs: 4 bytes/cell x 64 x 2000^2 is ~1 GB, well inside.
        let (t, notes) = plan_threads(2_000, 100, None, 64 * GB, Some(111 * GB));
        assert_eq!(t, num_cpus().min(100));
        assert!(notes.is_empty(), "should not editorialise on a run that fits: {notes:?}");
    }

    #[test]
    fn an_oversubscribed_run_is_throttled_not_refused() {
        // 30 kb against a 16 GB budget: 4 x 30000^2 is ~3.35 GB per thread.
        let (t, notes) = plan_threads(30_000, 100, None, 16 * GB, Some(64 * GB));
        assert!(t >= 1 && t <= 4, "expected a handful of threads, got {t}");
        assert!(
            notes.iter().any(|n| n.contains("--max-memory")),
            "should name the knob that would raise it: {notes:?}"
        );
        // The plan must actually fit the budget it was given.
        assert!(alignment_scratch_bytes(30_000, 30_000) * t as u64 <= 16 * GB);
    }

    #[test]
    fn an_impossible_budget_still_runs_at_one_thread() {
        // Deliberately not an error: the machine may have swap, and a warned
        // run tells the user more than a refusal.
        let (t, notes) = plan_threads(1_000_000, 50, None, 2 * GB, Some(4 * GB));
        assert_eq!(t, 1);
        assert!(
            notes.iter().any(|n| n.contains("running at 1 thread anyway")),
            "{notes:?}"
        );
    }

    #[test]
    fn explicit_threads_bypass_the_budget_and_say_so() {
        // 16 kb x 64 threads is ~63 GB, far over a 4 GB budget — honoured anyway.
        let (t, notes) = plan_threads(16_000, 100, Some(64), 4 * GB, Some(41 * GB));
        assert_eq!(t, 64, "--threads must win over the budget");
        assert!(
            notes.iter().any(|n| n.contains("--max-memory not applied")),
            "the budget must never be silently inert: {notes:?}"
        );
        assert!(
            notes.iter().any(|n| n.contains("may be killed")),
            "should still predict the consequence: {notes:?}"
        );
    }

    #[test]
    fn explicit_threads_are_still_capped_by_sequence_count() {
        let (t, notes) = plan_threads(2_000, 25, Some(64), 64 * GB, Some(111 * GB));
        assert_eq!(t, 25);
        assert!(notes.iter().any(|n| n.contains("reduced to 25")), "{notes:?}");
    }

    #[test]
    fn rmblast_is_the_default_backend() {
        let cli = Cli::try_parse_from(["te-composer", "in.fa"]).unwrap();
        assert_eq!(cli.backend, Backend::Rmblast);
        assert_eq!(cli.max_memory, "80%");
    }

    /// `--backend reference` is hidden but must stay reachable: it is the
    /// arbiter the other two backends are checked against.
    #[test]
    fn the_reference_backend_is_hidden_but_selectable() {
        let cli = Cli::try_parse_from(["te-composer", "in.fa",
                                       "--backend", "reference"]).unwrap();
        assert_eq!(cli.backend, Backend::Reference);
    }
}
