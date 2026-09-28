# The te-composer default pipeline

What `te-composer instances.fa` does when you give it no flags, and where it
differs from the Perl `Refiner`. Written against te-composer as of 2026-08-21;
Refiner references are to `RepeatModeler/Refiner` and `MultAln.pm`.

## Summary


te-composer takes a FASTA of instance sequences and produces one consensus.
The default run is four stages, preceded by a cheap sanity check:

0. **Check the input.** Report the instances' GC against the matrix background,
   and self-align the middle 2 kb of each of the ten longest instances to see
   whether they are tandemly repetitive. Costs about 34 ms, and only warns:
   the run continues either way.
1. **Pick a reference.** Align every input sequence against every other one.
   Score each sequence as a candidate. The winner seeds the first consensus.
2. **Refine, packing insertions as it goes.** Align all instances to the
   current consensus, recover inherited insertions the alignment has hidden,
   rebuild the consensus, repeat. Stop when it stops changing, when it starts
   cycling, or after 20 passes.
3. **Extend, then refine again** — *only with `--genome <2bit>`*. RAMExtend
   pushes the consensus past both edges using flanking genomic sequence. What
   it hands back is **not refined**, so stage 2 runs again over the widened
   instances.
4. **Repair bad blocks, then refine again, then judge — last.** Find stretches
   where the consensus disagrees with its own instances and patch them.
   Patching changes the consensus, so the MSA no longer matches it — a **full
   stage-2 loop runs again** on top of the patch to rebuild it. Then compare:
   **keep the repaired version only if the total alignment score went up**,
   otherwise throw it away and keep the unrepaired one. The block *selection*
   happens once; the refinement it triggers is a whole loop.
5. **Write out** the consensus.

So the full order is:

```
check input (GC, tandem repeats) — warns only
  → pick reference
  → refine+pack (loops)
  → [ extend → check consensus for tandem repeats → refine+pack (loops) ]
                                            only with --genome
  → patch blocks → refine+pack (loops) → judge: keep whichever scored better
  → out
```

So a default run contains **two** full refinement loops, and three when
extending. Both candidates the gate compares are fully refined, so the
comparison is like-for-like and whichever is kept has an MSA consistent with
its own consensus — a patched-but-unrefined consensus is never returned.

Stage 3 is bracketed because it needs an input you may not have: the genome the
instances came from, and Smitten-named sequences to locate them in it. It is
not a judgement call — if you have the genome, you almost certainly want it.

### Why repair runs last

**This is a deliberate departure from Refiner.** Refiner repairs *before* it
extends (`refineUntil → resolveLowQualityBlocks → refineUntil →
extendAlignment → refineUntil`), which has two consequences: the bases
RAMExtend adds are never offered to the repair at all, and the repair's
keep-or-discard verdict is reached against a consensus that no longer exists
by the time the run finishes.

Block repair is an evaluate-and-maybe-keep step, so it belongs after every
other modification — it should judge the consensus that is actually going to
be emitted. te-composer therefore runs the main loop without repair, extends
if asked, and applies the repair at the end against the final instance set
(widened, when the family was extended).

The loop count does not change: the refinement the repair triggers replaces
the one that used to follow it. Runs without `--genome` are byte-identical to
the old ordering — verified on simulated families — because with no extension
in between, "repair then refine" and "refine then repair" are the same two
steps in the same order.

> **Two things fix a consensus here, and they run at different rates.**
> *Insertion packing* runs inside **every** refinement pass — and there are two
> or three whole loops of those. *Block repair* selects blocks **once**, at
> the end, though it then triggers a full refinement loop of its own. Both are on by default and they are
> easy to confuse — see
> [Two repair mechanisms](#two-repair-mechanisms-dont-confuse-them).

Refiner's shape is the same — bootstrap, refine, repair, refine, optionally
extend — because te-composer is a port of it. The differences that matter:

| # | te-composer | Refiner |
|---|---|---|
| 1 | **Flags internally repetitive input**, twice: before the search, on the instances; and after extension, on the extended consensus. Also reports the instances' GC against the matrix background. All of it warns; nothing is refused. The post-extension warning is recorded in the Stockholm | No equivalent. Runs the search and returns whatever comes out |
| 2 | **DUST low-complexity masking is on.** Measured on 647 simulated families with known ancestors: mean distance to truth improves from 188.3 to 182.0, and the paired sign test against Refiner flips from losing (p = 0.0002) to winning (p = 0.0058). Costs 0.7% of HSPs on a well-behaved family | Emits `-dust no` unconditionally. Its `NCBIBlastSearchEngine` constructor calls `setUseDustSeg(1)`, but `getUseDustSeg` is never read when building the command line, so the setting is dead and masking is always off |
| 3 | **Phase 1 (reference selection).** No cull. The aligner runs an all-vs-all search with mask level off; te-composer splits each instance's HSPs by strand, trims each strand so that strand's **instance** ranges do not overlap (dropping anything left under `--min-row-len`, 25 bp), sums the **trimmed** scores per strand, and keeps the higher-scoring tiling. A candidate's score is then the plain sum of those surviving trimmed scores — the same quantity phase 2 sums, so one convention serves both. That same surviving set becomes the MSA rows | **Phase 1 (reference selection).** The all-vs-all search sets no mask level (rmblastn's `-1`, off); `findHighestScoringAlignmentSet` culls at 80% on the **reference** axis within each strand. **One surviving set then does both jobs** — summed raw to score the candidate, and handed to `MultAln` as the rows. Nothing reconsiders strand afterwards, so an instance that aligned in both orientations contributes both; and because the cull is pairwise, its overlapping fragments can credit the same reference bases more than once |
| 4 | **Phase 2 (MHSP selection).** The same rule as #3, applied to the instance-vs-reference search: no cull, trim each strand into a tiling, sum the trimmed scores, keep the better tiling. There is no candidate scoring here — the surviving set *is* the MSA. Reference overlap is deliberately left unbounded. Every row of an instance is on one strand, and no instance base appears in more than one row | **Phase 2 (MHSP selection).** rmblastn culls during the search (`-mask_level 80`): an HSP is dropped when more than 80% of its **instance** range is covered by a single higher-scoring HSP, strand ignored. Nothing filters afterwards. Survivors' instance bases may overlap by up to 80%, and one instance may produce rows on both strands |
| 5 | **Insertion packing, once the refinement loop has settled.** A consensus-induced MSA never aligns one instance's inserted bases to another's, so an insertion inherited by many instances stays invisible to the consensus caller however many passes run. When the column-wise loop converges, te-composer finds each run of consensus-gap columns, merges neighbouring runs while the span holds at most four called columns (`--pack-max-sep`), and keeps the spans where at least half of the spanning copies carry bases (`--pack-min-occupancy`). For each, the copies' bases are scored **all-against-all**; the one scoring highest against the rest becomes the centre, the others are aligned to it, and a consensus is called from that small MSA. It replaces the span only if it has **more called bases** than what is there now — this recovers bases the column-wise caller cannot see rather than relitigating ones it already called. The recovered bases join the reference and the loop runs again, at most twice. On by default (`--no-pack-insertions` disables) | **No equivalent.** `MultAln::_alignFromSearchResultCollection` reserves the widest insertion seen at each reference position, then **left-justifies** each instance's inserted bases in that block and pads the remainder with gaps on the right. The bases are never compared across instances, so whatever shares a column does so because the insertions happen to start at the same offset, not because they were aligned |
| 6 | Stops on any repeated consensus (cycle detection) | Stops only when a pass reproduces the previous consensus exactly |
| 7 | Up to **20** passes per refinement round; the loop also stops on a cycle (#6) | **11**, not 10: `$maxIterations = 10` but the loop is `for ($i = 0; $i <= $maxIterations; $i++)` |
| 8 | Extension is an in-process library call | Shells out to the `RAMExtend` binary |
| 9 | **Extension guards and flank measurement.** A side that is *capped* — it exhausted `--extend-max` (20 kb) without converging — is discarded on its own, and the other side kept. If what remains would still add more than `--extend-max-total` (25 kb), the whole extension is refused. Fixed the bug that caused same orientation neighbor seeds to limit extension to their midpoint rather than their edges | Refuses only when *both* sides cap (at `L=10000`), so a one-sided runaway is accepted whole and there is no ceiling on the total. And its `RAMExtend` tests same-strandedness with `s->strand == neighbor->strand`, comparing `char*` pointers from separate `cloneString` allocations — always false — so it midpoints same-strand neighbours too and under-extends |
| 10 | Block repair runs **last**, after the extension, so it judges the consensus that is actually emitted | Repairs before extending, so the extended bases are never repaired and the verdict is reached against a consensus that no longer exists |
| 11 | Finds bad blocks two ways: low-scoring segments **and** a 10-column scanning window | Low-scoring segments only |
| 12 | A repair is kept only if it improves the alignment score | Every repair is kept, unconditionally |

Everything else — the matrix, the score floor of 150, gap penalties of
-25/-5, word size 7, mask level 80, complexity-adjusted scoring — is the same
on both sides, deliberately.

---

## How many passes does it actually take?


Each phase-2 run reports its pass count and why it stopped:

```
te-composer: alu: phase 2 [initial]: 1 pass of 21, converged
te-composer: alu: extension anchored 30 copies, ...
te-composer: alu: phase 2 [after extension]: 2 passes of 21, converged
te-composer: alu: block repair: no low-quality blocks selected
```

`converged` and `cycled` mean the loop finished thinking. `EXHAUSTED — still
changing at the limit` means it ran out of budget mid-improvement, and is the
one case worth acting on. All of it goes to stderr, so stdout stays clean.

Measured across all 648 simulated families at the default budget of 21 passes:

| phase | families | mean passes | max | converged | cycled | exhausted |
|---|---|---|---|---|---|---|
| initial | 648 | **3.01** | 21 | 529 | 118 | **1** |
| after block repair | 310 | **2.19** | 8 | 249 | 61 | 0 |

35% of families converge in a single pass and 61% within two. Exactly **one
family out of 648** ever reached the 21-pass limit. So a run with three phase-2
loops is not three lots of 20 iterations — it is typically about five passes in
total.

Two things this makes visible:

- **Cycle detection earns its place.** 118 families (18%) stop because they
  started oscillating, not because they converged. Refiner has no cycle
  detection, so every one of those would run to its iteration limit swapping
  between two answers. This is what makes te-composer's larger budget safe.
- **The accept gate is doing real work.** Of the 310 families where blocks were
  selected, 152 repairs were kept and **158 discarded** — the gate rejects 51%
  of what the block selection proposes, consistent with the 45–53% measured
  earlier on other corpora. 305 families had no low-quality blocks at all.

## Stage 0: settings


These are the defaults, and where they come from.

| Setting | Default | Notes |
|---|---|---|
| `--backend` | `rmblast` | Seeded search. `parasail` and `reference` are full dynamic programming. |
| `--matrix` | built-in | Equivalent to `comparison.matrix`; mean diagonal 9.5. |
| `--gap-open` / `--gap-extend` | 25 / 5 | Read from the matrix's own `GAP` line, so the two cannot drift apart. Matches Refiner's `-25` / `-5`. |
| `--min-score` | 150 | Refiner's `setMinScore(150)`. |
| `--word-size` | 7 | Refiner's `setMinMatch(7)`. |
| complexity-adjusted scoring | on | Refiner's `complexityAdjustedScoreMode`. |
| X-drop cutoffs | derived | rmblast derives them from the score floor the way `NCBIBlastSearchEngine` does: `2x/÷2/x` = 300/75/150. Refiner never sets them either. |
| `--iterations` | 20 (rmblast), 3 (DP) | A pass is cheap for a seeded search and expensive for full DP. Not adjusted for family size. A length-based cap on long families fired on 1 of 4,823 families, changed nothing there, and on the one pathological family it did bind it only stopped the loop early with the consensus still changing. Removed; cycle detection (#6) is what bounds these runs. |
| `--hsps` | `tiled` | rmblast only — the DP backends return one alignment per pair anyway. **Differs from Refiner.** See Stage 2. |
| repeat check | always on | Warns about internally repetitive input, before the search and after extension. **No Refiner equivalent.** See below. |
| DUST masking | on | **Differs from Refiner**, which emits `-dust no`. See divergence #2. |
| `-n` / `--num` | 1 | One consensus out. |

### Repetitive-input checks

Two cheap checks, both of which only warn.

An earlier version refused. Across 4,823 families it fired four times, and all
four ran to completion in 0.3 to 1.9 s of wall clock at 13 to 30 MB of peak
RSS, on 8,113 to 14,895 bp of input each. There was no expensive run to avoid
there, and a refusal would have withheld four consensi that cost nothing to
produce. Whether a satellite family is worth keeping is a curator's
decision, and it is cheaper to make it by looking at the consensus than by
reading a refusal.

The post-extension warning is written into the Stockholm output, so it travels
with the family. The pre-search one goes to stderr only.

#### Before the search, on the instances

An all-vs-all search costs time proportional to the square of the instance
count, and on internally repetitive sequence the aligner also returns far more
HSPs per pair. One real family — 60 instances at 69% GC, **608,515 bp of input
in total** (median instance 11,893 bp, longest 20,004) — produces **8.7 million
HSPs**, takes 401 s and 14,311 MB of peak RSS, and yields a 12,154 bp consensus
that RepeatMasker resolves into a 1,268 bp unit tiled five times. It is a tandem
array, not an element.

te-composer self-aligns the middle **2 kb** of each of the **10 longest**
instances and counts the HSPs. The window is deterministic, not sampled, so a
rerun gives the same verdict. A sequence that is not internally repetitive
aligns to itself as one diagonal HSP plus a little noise; a tandemly repetitive
one matches itself at every period.

| input | self-HSPs in a 2 kb window |
|---|---|
| random 20 kb, 10% divergence | 1 |
| worst curated hs1 family (426 tested) | 20 |
| worst simulated family (648 tested) | 27 |
| pure `(GGGAGG)n` | 57 |
| the 8.7M-HSP family above | 78 |

The threshold is **40**, which flags neither corpus and catches the
pathological cases.

#### After extension, on the consensus

The instance-level probe cannot see a satellite whose instances are short. In
`hs1/rnd-5_family-3330` the instances are **71 bp** — far too short to
self-align into anything — and extension grew them into 15 kb that is 50%
`(GGAAT)n`, HSATII. RepeatMasker annotates 681 bp of that 15,024 (**5%**), all
scattered Alu fragments. The extension walked out through chromosome 14 rather
than along an element.

So the same probe runs again on the extended consensus, where the repeat is
plainly visible. It scores 58 self-alignments there, against 4 and 6 for
genuine 18 kb and 28 kb extensions, so the threshold carries over unchanged.

By this point the consensus already exists, so declining it would save nothing.

#### Composition

The instances' GC is reported at the start, against the scoring matrix's own
background (`FREQS A .265 C .235 G .235 T .265`, i.e. 47% GC) and flagged when
it falls outside 35–60%:

```
  composition: 69% GC  <- far from the 47% matrix background
```

That 69% is the property that made the family above produce 8.7M HSPs: the
score floor assumes a 47% background, and GC-rich purine-biased sequence
matches itself far above it, so it is worth reporting up front.

#### Three details for anyone changing the probe

* **It runs with DUST off**, whatever the shipping default is. DUST masks
  low-complexity query sequence, so with it on a pure tandem repeat self-aligns
  to nothing and the probe reports 0 — the one input it most needs to catch.
* **Truncation is what makes it affordable.** Probe cost tracks the HSP count
  it is hunting: self-aligning a full 20 kb pathological instance takes 6.07 s
  against 0.14 s for a clean one. At 2 kb the same instance costs 0.04 s, 150x
  less, and still separates by 100x. Across 648 families the check adds 22 s,
  about 34 ms each.
* **A flat threshold only works because the window is fixed.** Raw self-HSP
  counts scale with length — a random, non-repetitive 20 kb sequence yields
  21–29 — so a threshold set against short sequences would flag every large
  family. Probing a constant 2 kb removes the length term.

Refiner has no equivalent to any of this; it runs the search and returns
whatever comes out.

## Stage 1: pick a reference


te-composer has no consensus to start from, so it borrows one. Every input
sequence is aligned against every other, and each is scored as a candidate to
be the starting point. The best-scoring sequence becomes the reference, and the
alignment of everything else to it gives the first consensus.

Two things happen to the alignments before scoring:

1. **Trim each strand into a tiling.** Within the forward set and within the
   reverse set separately, the highest-scoring alignment claims its instance
   bases first; each lower-scoring one is trimmed to the largest stretch of its
   own range nobody has claimed, and dropped if fewer than `--min-row-len`
   bases survive. A trimmed alignment carries a proportionally reduced score.

   *Refiner does not trim in phase 1.* It culls instead — dropping an alignment
   whose **reference** span is more than 80% covered by a higher-scoring one
   from the same instance and strand.

2. **Keep the better strand.** The trimmed scores are summed per strand and the
   losing strand is dropped entirely.

   *Refiner keeps both strands and sums them.* An instance that aligns to the
   candidate in both orientations contributes twice there, once here.

The candidate's score is then the plain sum of what survives. Coverage of one
region by *different* instances counts in full from each — that is real support
from separate copies. Coverage of one region twice by the *same* instance also
counts in full, because trimming works on the instance axis and leaves
reference overlap alone; an earlier version discounted that with a novelty
weight, but over 648 simulated families the weight left 645 consensi
byte-identical, so the weight was removed.

Both tools keep the reference row in the bootstrap consensus — te-composer's
`include_reference: true`, Refiner's `consensus( inclRef => 1 )`, whose comment
gives the reason: it is a real input sequence. Stage 2 drops it and calls the
consensus from the instance rows alone.

One more difference in bookkeeping: te-composer builds an MSA and a consensus
for *every* candidate and then sorts, where Refiner picks the winner on score
first and builds a single MSA. Same answer, more work.

Note the axis: stage 1 filters on the **reference** axis (which parts of the
candidate are covered). Stage 2 filters on the **instance** axis. They are
different questions and the code uses different settings for each.

## Stage 2: refine


Now there is a consensus, so the work becomes: align every instance to it,
rebuild the consensus, repeat.

Each pass:

1. Align all instances against the current consensus, with the aligner's own
   masking **off**.
2. Filter overlapping alignments on the **instance** axis at mask level 80,
   forward and reverse separately. Refiner gets a similar reduction from
   rmblastn's `-mask_level 80` instead, but that one is strand-blind — see
   [Overlap culling](#overlap-culling-what-is-shared-and-what-is-not).
3. Keep only the better strand per instance, judged on the summed score of
   what survived step 2.
4. Trim each strand into a tiling, then keep the higher-scoring one (`tiled`).
5. Build the multiple alignment, drop the reference row, and call a new
   consensus from the instance rows alone.

> **Measurements pending.** The comparative figures that used to sit here were
> taken before the pre-release cleanup, on a different default policy and
> without DUST. They are being re-run against the release build on the
> simulated, hs1 and RepeatModeler raw-family corpora, and will be restored
> once measured.

**How rows are counted (`--hsps`).** An instance can align to the consensus in
several pieces, and the three policies differ in what becomes of them:

| value | an instance contributes | |
|---|---|---|
| `best` | 1 row | its highest-scoring HSP; the rest of its coverage is discarded |
| `all` | n rows | Refiner's shape — every survivor, overlaps and all |
| **`tiled`** | n rows | every survivor from the winning strand, trimmed so no genomic base appears twice |

The flag applies to `--backend rmblast` only; the DP backends report a single
optimal alignment per pair, so there is nothing to tile.

**When it stops.** Refiner stops when a pass produces a consensus identical to
the one it started from. That catches a fixed point but not a cycle: if pass A
produces B and pass B produces A, Refiner oscillates until it runs out of its
10 iterations. te-composer remembers every consensus it has seen and stops on
any repeat, which makes the larger budget of 20 passes safe.

**Consensus caller.** te-composer uses the Dfam caller (per-column argmax,
ties prefer `N`, plus CpG restoration) by default. `--orig` selects the older
GIRI caller. Refiner calls `MultAln::consensus` with `linupmatrix`.

## Stage 3: extend against a genome

*Only with `--genome <2bit>`.* A consensus can only describe the part of
the element its instances cover; this pushes it past both edges.

A consensus can only describe the part of the element its instances cover. If
you supply the genome the instances came from, te-composer extends past both
edges with RAMExtend and then **refines again** against the widened instances,
which is what Refiner does with the extension it gets back.

Input names must be Smitten identifiers (`chr1:1000-2000_+`, or with an
assembly prefix `hg38:chr1:1000-2000_+`) because that is what locates an instance in
the genome. Two separate things can stop an instance contributing to the extension:

- **It does not reach the consensus edge.** Its alignment stops more than 10
  columns short, so it knows nothing about what lies beyond. It stays in the
  alignment and still contributes on its other side.
- **It cannot be located in the genome.** Unparseable name, wrong assembly,
  sequence missing from the 2bit, coordinates past the end, or bases that do
  not match the genome there (below 95% identity). It is dropped from the
  extension with a warning naming the instance and the reason, and stays in the
  family.

Matrix and minimum improvement are chosen from the family's Kimura divergence
using the same buckets as `ram-extend -stk` and `extend-stk.pl`
(`<16 → 14p43g, <19 → 18p43g, <22.5 → 20p43g, else 25p43g`).

Refiner shells out to the `RAMExtend` binary and parses its files back;
te-composer calls the library in-process.

Extension runs **once**, not iteratively.

### Overlap culling: what is shared, and what is not

It is easy to read divergence #4 as "Refiner keeps everything, te-composer
filters". That is wrong: both tools start from the same raw HSP set and both
reduce it at 80% overlap. Culling can happen in two places — inside the aligner
(`rmblastn -mask_level`, default off) and after it (Refiner's
`findHighestScoringAlignmentSet`, ported by te-composer as `RefinerFilter`) —
and the two tools put the work in different places in different phases. So the
comparison has to be per phase.

#### Phase 1 — reference selection

| step | te-composer | Refiner |
|---|---|---|
| search | every instance against every other, self-hits skipped | the instance file against itself; self-hits deleted by the filter |
| in-aligner mask | `mask_level: 101` → off; every HSP is returned | **none set** → rmblastn's `-1`; every HSP is returned |
| cull | **none** — trimming subsumes it | `findHighestScoringAlignmentSet`: drop an HSP when more than 80% of its **reference** range is covered by a higher-scoring one from the same instance and strand |
| trimming | each strand trimmed so its **instance** ranges do not overlap; anything under `--min-row-len` (25 bp) dropped | none |
| strand | the **trimmed** scores are summed per strand and the losing strand is discarded (#3) | survivors of both strands are **pooled** into the candidate's score |
| candidate score | plain sum of the surviving **trimmed** scores | plain sum of surviving scores, untrimmed |
| tie-break | lower input index | lower sequence ID (`$a cmp $b`) |
| consensus | built for the candidates that will be used, reference row included | built once, for the winner only, reference row included (`inclRef => 1`) |

Note what phase 1 does *not* do: nothing bounds how much of the reference an
instance may cover twice. Trimming works on the instance axis only, so an
instance covering one consensus region twice from *different* bases — a tandem
expansion — contributes both alignments to the candidate's score. An earlier
version discounted that with a novelty weight; measured over 648 simulated
families the weight left 645 consensi byte-identical, and of the three that
moved, two were worse with it on. It was removed.

Refiner has neither. Its cull is pairwise on the reference axis, so an
instance's overlapping fragments can each clear the 80% test against every
individual competitor and then each contribute their full score for reference
bases already counted.

#### Phase 2 — MHSP selection

The two tools no longer share a mechanism here, so take them one at a time.

**te-composer applies no cull.** Per instance it splits the HSPs by strand,
trims each strand so that strand's instance ranges do not overlap, drops
anything left under `--min-row-len` (25 bp), sums the **trimmed** scores per
strand, and keeps the higher-scoring tiling. A trimmed HSP carries a
proportionally reduced score (`score × surviving ÷ original`), so the strand
vote counts only evidence that becomes rows.

The cull can go because trimming already covers it. When one HSP is
contained in another, both rules discard it. They part company only in a narrow
band — an HSP mostly but not wholly covered, where the cull deletes it outright
and the trim keeps its unclaimed remainder — and there the trim is the better
rule: it keeps novel bases as evidence, and it thresholds in absolute bp rather
than as a percentage of a span.

**Refiner's whole filter is `-mask_level 80`**, applied by rmblastn during the
search: an HSP is dropped when more than 80% of its instance range is covered by
a single higher-scoring one, strand ignored. Nothing filters afterwards. That is
worth stating plainly because it is surprising: **nothing else in Refiner's
phase 2 reduces the HSP set.** Verified by reading the path end to end —
`refineUntil` (line 576) is a bare loop; `refineConsensus` (1230–1387) contains
no filtering call other than `setMaskLevel(80)` at 1252; and the row builder it
hands the collection to, `MultAln::_alignFromSearchResultCollection`, has no
overlap, strand or dedup logic at all — its only `next` is a gap-pattern guard,
so every `SearchResult` becomes exactly one row.

`findHighestScoringAlignmentSet` — the per-strand cull that measures overlap on
the reference only — is real, but it is called exactly once, at line 1101,
inside `bootstrapConsensus`. It is phase 1 and never runs in phase 2.

Two details of Refiner's rule are what make the difference bite.

**It compares against one HSP at a time, not against everything kept so far.**
Build a subject holding three segments of the same query — an exact copy of
query 150–300, an exact copy of 300–450, and a 15%-diverged copy of 200–400 —
and all three survive at `-mask_level 80`:

```
150 300     1 151  plus  1423
300 450   212 362  plus  1415
200 400   423 623  plus  1229     <- kept
```

The third HSP has **no** instance range of its own: it is 100% covered by the
union of the other two. But it is only ~50% covered by either one alone, so
neither individually clears the threshold and it stays. This is why "keep only
HSPs with 20% unique coverage" is not quite what either tool does — *unique*
would have dropped it.

**Refiner's version ignores strand; te-composer's does not.** rmblastn measures
overlap on the query range and never looks at orientation, so a reverse HSP is
compared against a forward one. Subject built as `S + revcomp(S)`, so one query
range matches both ways:

```
-mask_level 101   1 600      1  600  plus   5685
                  1 600   1250  651  minus  5685      <- both returned
-mask_level 80    1 600      1  600  plus   5685      <- reverse culled
```

But that is still a test between two HSPs, not a decision about the instance. Two
*non-overlapping* instance ranges matching on opposite strands both survive:

```
-mask_level 80  301 600  600  301  minus  2841
                  1 300    1  300  plus   2824
```

So Refiner collapses cross-strand overlap where the instance ranges overlap,
and nowhere else. An instance whose left half matches forward and right half matches
reverse still contributes rows on both strands, and same-strand survivors
overlapping by 80% or less both become rows with the shared bases counted in
each.

**te-composer's order, per instance:**

1. **Trim each strand into a tiling.** Within the forward set and within the
   reverse set separately, the highest-scoring HSP claims its instance bases
   first; each lower-scoring one is trimmed to the largest stretch of its own
   range nobody has claimed, and dropped if fewer than `--min-row-len` bases
   survive. On the three-HSP example above, the third gets nothing and is
   dropped.
2. **Vote on the trimmed scores.** Sum each strand's trimmed scores and keep the
   higher-scoring tiling — forward wins ties. An instance was inserted in one
   orientation, so at most one of the two sets can be right.
3. **Each survivor becomes one row.** An instance contributes one or more rows,
   all on the same strand, together covering its sequence without overlap.

Trimming first is the point: the winning strand is the one contributing more
*rows*, not the one with more raw HSP score. An earlier version
voted first and trimmed only the winner, so a strand could win on evidence that
was then cut away. Measured end to end the change was neutral, so it is a
consistency fix rather than an accuracy one.

Three limits worth being explicit about:

* **Trimming has no reference-axis limit** — there is no maximum-overlap
  setting to tune. It looks only at instance ranges, and two rows of one
  instance may overlap freely on the consensus. That is deliberate, and guarded
  by a test (`reference_overlap_is_left_alone`): an instance covering the same
  consensus region twice from *different* bases is a tandem expansion, and two
  overlapping rows is a legitimate rendering of it.
* **The rows need not be co-linear.** Nothing requires an instance's rows to
  advance together on both axes, so an inversion or rearrangement inside an
  element survives as separate rows rather than being resolved into one.
  Measured on 79 families this matters less than it sounds: **0.53%** of
  simulated instances and **0.78%** of curated human ones end up with rows that
  are not co-linear. A co-linear chaining policy was built and measured — it was
  indistinguishable from tiling and is not shipped.
* **The trim order uses untrimmed scores.** HSPs are ranked by the score the
  aligner reported, so a long mediocre HSP claims bases ahead of a short
  excellent one. This is circular by nature — a trimmed score cannot be known
  before trimming — and ranking by score density instead is the obvious
  alternative, unmeasured.

**How much of this engages at all.** Measured over 79 families, the median
instance contributes exactly **one** row (mean 1.12 on simulations, 1.07 on
curated human), and only **7.6%** / **8.1%** of instances have more than one
HSP to reconcile. Trimming, the strand vote and the co-linearity gap all operate on that
minority. That bounds how much any of them can matter on these corpora; it is
not a claim about TE biology, where inversions and rearrangements are real.

**Which axis the overlap is measured on changes with the phase**, in both
tools. Phase 1 measures on the **reference**; phase 2 measures on the
**instance**. Using the reference axis in phase 2 was tried and is measurably
wrong: in a dimeric element an instance's left arm also aligns to the
consensus's *right* arm, which is novel on the reference axis, so a
reference-axis filter keeps it and the alignment fills with arm-swapped rows —
substitutions rose by an order of magnitude on simulated AluY families.

### What the extension will refuse

An extension is only as good as the edge it found. Each side is judged
separately:

* **A side that reaches `--extend-max` (20,000 bp) is dropped.** Reaching the
  cap means the DP was *stopped*, not that it converged — nothing about the
  sequence justifies where it ended. At 20 kb it is well past any real element
  (a full-length ERV runs to ~9–10 kb), so a side still going at the cap is
  tracking an array or a duplication, not a boundary.
* **The other side is kept.** A side that converged found an edge whatever its
  neighbour did. Refusing the whole extension because one side ran away cost
  190 kb of real element across 35 mouse families before this was split.
* **If neither side converges, or the survivors together exceed
  `--extend-max-total` (25,000 bp), the whole extension is refused** and the
  family keeps its unextended consensus.

Every one of those decisions is written into the Stockholm record, because a
refused extension leaves a consensus indistinguishable from one that was never
offered an extension:

```
#=GF ** TE-COMPOSER: Extension: left extension hit the 20000 bp limit — no edge found on that side, dropping it
#=GF ** TE-COMPOSER: Extension: anchored 16 of 16 copies, divergence 1.8%, matrix 14p43g
#=GF ** TE-COMPOSER: Extension: right 13383 bp only (left side dropped or empty)
```

**Refiner refuses only when both sides cap**, at `L=10000`, testing the
reported lengths against `> 9990`. That misses the one-sided runaway: a 45 bp
human family extended 20,000 bp left and 5,102 right, was never "capped both
ways", and became a 25 kb consensus matching four different TE classes.

### Why te-composer extends further than Refiner

Not a difference in policy — a **bug in the C tool that the Rust port fixed**.
Where two cores could extend toward each other the gap is split at the
midpoint, but a *same-strand* neighbour should cede the full distance. The C
tests this with `s->strand == neighbor->strand`, comparing `char*` pointers
from separate allocations, which is always false, so it halves the flank for
same-strand neighbours too.

The midpoint split itself is deliberate and is kept: two cores extending toward each other in one directional pass each take half the gap so their extensions meet without overlapping. What the fix changes is only *which* pairs get it. Two same-strand neighbours travel the same genomic direction in a given pass, so they never collide and never needed splitting. Neither branch extends into a neighbour's core — the distance being split is the gap *between* them — so restoring the full gap does not reintroduce overlapping extensions. A test pins both halves: `same_strand_neighbors_split_full_vs_midpoint` asserts the full gap for same-strand pairs and the midpoint for an opposite-strand neighbour whose facing flag is set.

The signature is unmistakable: identical extension on one side, very different
on the other. On `hs1 rnd-5_family-3330`, Refiner gets left 7120 / right 5115
and te-composer left 7120 / right 7779. Running the port with `-ccompat`, which
reproduces the C bug, gives right 4860 against 7779 — left unchanged.

So **Refiner under-extends wherever a same-strand neighbour sits on that
side**, and it will keep doing so until RepeatModeler is pointed at the Rust
`ram-extend` rather than the C binary.

## Two repair mechanisms, don't confuse them


te-composer has two separate ways of correcting a consensus against its
instances. Both were tested heavily, which is exactly why they get mixed up.
Only one of them is on by default.

| | Block repair | Insertion packing |
|---|---|---|
| **Runs** | **Once**, between two refinement runs | **When the refinement loop settles**, then the loop runs again; at most two rounds |
| **Fixes** | Stretches where the consensus disagrees with its instances, usually a length disagreement | Insertions carried by many instances that a consensus-induced MSA never aligns to one another |
| **Judged by** | Accept gate — kept only if the total alignment score improves | No score gate; a span is packed only if at least half of the spanning copies carry bases in it (`--pack-min-occupancy`), the centre aligns positively to the rest (`--pack-min-score`), and the result has more called bases than the span holds. The loop that follows keeps a packed base only if the re-aligned copies support it |
| **Comes from** | Refiner's `resolveLowQualityBlocks`, plus `AutoRunBlocker`'s window | New. No Refiner equivalent |
| **Default** | **On** | **On** — disable with `--no-pack-insertions` |

Packing was built to test a specific idea: that a consensus-induced MSA never
aligns one instance's inserted bases to another's, so an insertion inherited by
many instances is invisible to the consensus caller no matter how many passes you
run. Recovering those spans *inside* each pass means a recovered base joins the
reference for the next pass and compounds — a longer reference gives the
following round's alignments more to anchor on. If that worked, one refinement
round would be enough and block repair would become redundant.

Measured on 648 simulations where the true ancestor is known exactly, with all
four arms produced by the same binary (2026-08-21, `benchmarks/score_pack_arms.py`):

The arms were measured against each other on simulated and hs1 families;
the figures are being re-measured against the release build.
Both are on by default, and **running both beats running either alone.**

Packing was originally designed to *replace* block repair —
`run_pack_arms.sh` says so in its header: *"with the re-alignment inside each
iteration there is no call for a second refinement round, and keeping the old
repair alongside it would only obscure which of the two did the work."* That
turned out to be wrong. The original experiment never ran a
packing-plus-repair arm, so the question went untested until it was measured
directly:

The two mechanisms were measured against each other and are complementary
rather than redundant; the figures are being re-measured.
The two do not overlap as much as expected because they fix different things:
packing recovers bases the alignment hid, block repair fixes length
disagreements the consensus already had. The accept gate keeps repair honest —
it is kept only when it improves the total alignment score.

The sims/hs1 disagreement is **unresolved, not settled**. One live question is
whether the hs1 metric can detect the improvement packing makes at all, given
its "truth" is a curated consensus produced with Refiner's help.

## Stage 4: repair bad blocks


This stage runs **once**, after refinement has settled — not on every
iteration. After refinement settles, some stretches of the consensus still
disagree with the instances underneath them — usually a length disagreement, where the
consensus has picked up or lost a base or two relative to what most instances
actually have.

**Finding candidate blocks.** Two selections run, and their results are
combined:

- *Low-scoring segments.* Score each alignment column against the matrix and
  run Ruzzo–Tompa to find maximal runs of poor columns (threshold 1). This is
  Refiner's `getLowScoringAlignmentColumns`.
- *Scanning window.* Slide a 10-column window across the consensus and test
  each position the same way. This is `AutoRunBlocker`'s selection, and
  **Refiner does not have it**.

The window matters because roughly 44% of what it finds sits inside a
low-quality block whose *whole-block* vote says the consensus is already
right — the disagreement is only visible at window scale. Measured on 791 hs1
families against curated Dfam consensi, adding it roughly doubled what the
repair recovers (+99 → +223 net bases closer to curation).

**Resolving a block.** For each candidate block, look at the instance
sequences spanning it:

- If one length is held by a majority of the instances and differs from the
  consensus's length there, rebuild that block from just the instances of that
  length. (Refiner additionally requires the block to be 2–50 columns wide,
  with at least 4 instances and at least 3 agreeing.)
- Otherwise, align the block's instances against each other, pick the one with
  the best summed score, and use it. Refiner uses Needleman–Wunsch–Gotoh with
  `linupmatrix`; te-composer uses a parasail global alignment with the same
  penalties. That branch is always taken now; skipping it leaves every
  no-majority block unrepaired.

**The accept gate — the important difference.** Refiner patches the consensus
with whatever the blocks produced and moves on. te-composer patches, runs a
full extra refinement from the patched consensus, and then compares the total
alignment score against the unpatched result. If it did not go up, the whole
repair is discarded and the unpatched consensus is kept.

This is not a formality: the gate rejects roughly 45–53% of proposed repairs.
Without it, the repair is a net loss on the simulated benchmark; with it, it
is a clear win. `--repair-accept mean-per-base` scores by mean score per
aligned base instead of the total, and `--repair-accept always` disables the
gate (which is Refiner's behaviour).

Turn the whole stage off with `--no-repair-blocks`.

## Stage 5: output

**Stockholm is the primary output** — it carries the alignment *and* the
consensus, so nothing is lost:

```
te-composer instances.fa out.stk                              # MSA + consensus
te-composer instances.fa out.stk --consensus c.fa --format fasta
```

The positional output is the Stockholm record (stdout if omitted);
`--consensus` writes the bare consensus in `--format ig|fasta` for tools that
want it on its own.

```
# STOCKHOLM 1.0
#=GF ID MyFam
#=GF ** TE-COMPOSER: BOOTSTRAP-REF: chr22:33101-33204_-
#=GF ** TE-COMPOSER: SCORE: 3038
#=GF ** TE-COMPOSER: DIVERGENCE: Kimura 7.82%, CpG-adjusted 5.91%
#=GF ** TE-COMPOSER: Extension: anchored 30 of 30 copies, divergence 7.8%, matrix 14p43g
#=GF ** TE-COMPOSER: Extension: left 101 bp, right 100 bp
#=GF SQ 30
#=GC RF       TTATGGAGGCTGAGAAGTCCCACGAT...
chr22:33101-33204_-  .TATGCAGGCAGAGAAGTCCCATGATC...
//
```

Row identifiers are **normalized Smitten ranges**. The alignment writer labels
each row by appending its aligned extent to the input name, which on a
Smitten-named instance yields a recursive identifier
(`chr22:33101-33204_-:1-104_+`); te-composer collapses that to the single
absolute range it denotes. Names that are not Smitten identifiers are left
alone.

Everything the pipeline decided on the family's behalf goes into `#=GF **` —
the reference it chose, how it stopped, and what the extension did or refused
to do. A curator reading the family months later cannot recover that from the
sequence.

The record passes `stk lint` with no structural errors; the only complaints are
missing curation fields (`DE`, `AU`, `TP`, `OC`), which is correct for a family
that has not been curated yet.

## The run report

Progress goes to **stderr**, so stdout stays a clean Stockholm stream.
`--silent` turns it off entirely.

```
te-composer 0.1.0
  invocation: te-composer --genome hs1.2bit instances.fa out.stk
  input: 30 instances, 3065 bp total; lengths 91-115 bp (median 102, q1 98, q3 107)

── consensus bootstrap ──
  reference: chr22:33101-33204_- (104 bp, index 0), score 20689, 2.5% ahead of the runner-up (20166) of 30 candidates
  consensus 104 bp, score 22736, 30 of 30 instances participating, Kimura 7.83% (CpG-adjusted 4.56%)

── iterative refinement ──
  1 pass of 21, converged
  consensus 104 bp, score 22736, 30 of 30 instances participating, Kimura 7.83% (CpG-adjusted 4.56%)

── extension ──
  pre-extension core boundaries:
   Seq Ident Range        Orient L/R?   Left-Flank Core                       Right-Flank
     0 chr22 33100-33203  -      1/1    CTAACATGGT [GAAACCCCGT....GGCGTGAACC] AAGGAGGCAG
  anchored 30 copies, divergence 7.8%, matrix 14p43g; extended left 101 bp, right 100 bp

── iterative refinement ──
  2 passes of 21, converged
  consensus 305 bp, score 63815, 30 of 30 instances participating, Kimura 9.56% (CpG-adjusted 5.20%)

── block repair ──
  no low-quality blocks selected
  consensus 305 bp, score 63815, 30 of 30 instances participating, Kimura 9.56% (CpG-adjusted 5.20%)

── summary ──
  AluY_test: consensus 305 bp
  instances:  30 of 30 participating
  divergence: Kimura 9.56%, CpG-adjusted 5.20%
  score:      63815
  runtime:    0.1s
```

The same three numbers appear at every stage so the stages compare directly.
Three things on those lines are worth knowing:

* **`EXHAUSTED — still changing at the limit`** replaces `converged`/`cycled`
  when a loop runs out of budget mid-improvement. It is rare (1 family in 648
  simulated) and is the one stopping condition worth acting on.
* **`N of M instances participating`** is the fraction of the input the
  consensus actually rests on. A family where a third of the input never
  aligns is a different object from one where all of it does.
* **The pre-extension core view** is ported from the C tool's `printCoreEdges`,
  which the Rust port had not carried across. It shows what the extension is
  about to reason over: if every instance shows the *same* flanking sequence, the
  family is a segmental duplication and the extension is about to run into
  shared context.

## Starting from an existing consensus

```
te-composer instances.fa out.stk --bootstrap-cons existing.fa
```

The bootstrap phase exists only to pick a starting consensus from the
instances. Given one, it is skipped:

```
── consensus bootstrap ──
  skipped — starting from supplied consensus alu (104 bp)
```

Everything downstream — refinement, extension, block repair — is unchanged.
Use it to improve an existing family against a new instance set, or to
re-refine a curated consensus, rather than rediscovering a reference the
instances may not agree on.

## Insertion packing

See [Two repair mechanisms](#two-repair-mechanisms-dont-confuse-them) for what
this does and why it ships. Disable it with `--no-pack-insertions`.
`--pack-insertions` is still accepted and ignored, so scripts that switched it
on explicitly keep working.

The knobs are `--pack-max-sep` (merge neighbouring gap runs while the span
holds at most this many called columns in total, default 4), `--pack-min-seg` (a span needs one instance
contributing at least this many bases, default 5), `--pack-min-occupancy` (the
fraction of copies spanning the region that must carry bases in it, default
0.5), `--pack-min-score` (the all-against-all winner must beat this summed
score, default 0).

**Why packing waits for the loop to settle, and why the occupancy gate.**
Packing originally ran inside every pass with no occupancy gate: any four
copies with a shared insertion qualified, and each qualifying span cost an
all-against-all alignment with traceback. Two things followed. Any consensus-gap
span that four copies happened to fill was re-derived, so an insertion carried
by 5 of 15 spanning copies was written into the consensus that the column
caller had correctly declined. And the work grew with the number of private
insertions times the square of the copy number, pass after pass: on a
100-copy, 12 kb family at 20% divergence packing took about 60 s per pass and
the run 374 s against 27 s without it. Measured 2026-09-12; the earlier
estimate of a 16% cost came from shallower families.

Packing now runs when the column-wise loop has converged, cycled, or spent its
budget. te-composer keeps a span only when at least half of the spanning copies
carry bases in it, the same majority the column caller demands once the copies
are re-aligned to the packed reference. It picks the centre on score-only
alignments and runs one traceback per member. If packing changed the
consensus, the loop runs again from it with a fresh pass budget; what is
emitted is therefore always a fixed point of the column caller, and a packed
base survives only if the re-aligned copies support it. At most two packing
rounds run per refinement.

---

## Quick flag reference

| Flag | Effect |
|---|---|
| `--genome <2bit>` | Extend against a genome, then refine again |
| `--bootstrap-cons <FASTA>` | Supply the starting consensus; skip the bootstrap phase |
| `--consensus <FILE>` | Also write the bare consensus, in `--format ig\|fasta` |
| `--hsps <best\|all\|tiled>` | How an instance's HSPs become MSA rows. Default `tiled` |
| `--silent` | Suppress the progress report on stderr |
| `--no-repair-blocks` | Skip block repair entirely |
| `--no-pack-insertions` | Turn off per-pass insertion packing |
| `--extend-max <BP>` | Cap per side, default 20000 |
| `--iterations <N>` | Refinement passes after the first; default 20 (rmblast), 3 (DP) |

Sixteen further flags exist for experiments and are hidden from `--help`; they
still work when given explicitly. They cover block-selection tuning, packing
parameters, orientation and scoring switches, x-drop, and mask level.

To run te-composer as close to Refiner as its options allow:

```
te-composer --hsps all --repair-accept always --repair-window 0 \
            --no-best-orientation --no-pack-insertions --no-dust \
            --iterations 10 instances.fa
```

Three rows of the table cannot be flagged away. **Stage ordering** (#10):
repair runs last and there is no switch to move it back — though without
`--genome` the two orderings are equivalent and the outputs byte-identical, so
this recipe *is* exact for runs that do not extend. **Cycle detection** (#6):
the loop still stops on a repeated consensus rather than only on a fixed point,
which changes the stopping point for 18% of families. And **#9** contributes
two: the extension guards have no off switch, only wider thresholds, and the
same-strand flank fix is a property of the port rather than a setting.
