# recover-seeds

Builds seed alignments and consensi from a RepeatMasker `.align` file under two
sampling policies, and reports how each consensus compares with the library
consensus it was searched with.

This tool never realigns copies. The `.align` file is its only source of
alignments, so a difference in output comes from the sampling rather than from
a different aligner. parasail runs in one place only, aligning one consensus to
another for the ConsCmp report.

## Policies

`population` reproduces `generateSeedAlignments.pl`: drop the most diverged
quartile of copies, take the longest first up to `--target-copies`, then keep
adding any copy that covers a part of the consensus still below `--min-depth`.

`closest` takes the copies nearest the library consensus, reconstructing the
small, closely related sample RepBase consensi look like they were built from.
It takes `--closest-copies` of them, ten by default, or fewer when the genome
holds fewer. Span ranks copies rather than excluding them, the way the
population rule puts its divergence outliers last instead of dropping them:
copies covering `--min-span-fraction` of the consensus come first, and shorter
ones fill the target when too few long copies exist.

Set `--closest-copies 0` to try each count in `--closest-grid` and keep
whichever lands nearest the library consensus. That search scores a sample
against the consensus it is meant to explain, and on the Aedes set it settled
on one or two copies for half the families, too thin to call a consensus from.

Each policy gets its own Stockholm file, carrying the sample and the called
consensus as `#=GC RF`. The summary TSV holds both side by side.

## The source-sample scan

Those scores are worth having with the size fixed too, so recover-seeds records
them at every `--closest-grid` size: the identity each sample's consensus
reaches against the library, and how far the sample it kept falls short of the
best one. `--no-source-scan` turns this off and saves most of the policy's
runtime.

Transib5_AA in Aedes aegypti reads

    SourceScan: 1:99.22 2:75.17 3:74.16 5:71.47 7:70.72 | best=1 copies at
    99.22%, kept 7 copies at 70.72%, drop=28.51, hint=narrow

One genome copy reproduces the RepBase consensus to 99.2%, and adding the
second drops it to 75%. RepBase most likely built this family from that one
insertion, worth knowing before reading the ten-copy seed as a disagreement.
The `hint` reads `flat` when the best size lands within a point of the kept
sample, `narrow` when five or fewer copies beat it by five points or more, and
`mixed` in between; `scan_best_k`, `scan_best_identity`, `scan_drop`
and `source_hint` carry the same numbers in the summary TSV.

### How a grid point is scored

recover-seeds runs a full pass per size: take the first n of the ranked list,
build the MSA, call the consensus, align it to the library consensus with
parasail, and keep the `ConsCmp` identity over the whole reference. Every
sample is a prefix of the same ranking, so the profile follows that one order
instead of drawing afresh at each size. A cliff like Transib5_AA's is
therefore attributable to one copy: the second instance in the ranking moves
the consensus 24 points.

Two kinds of grid point disappear without a note. A size the family cannot
reach repeats the previous sample and is skipped, which is why a ten-value
grid often prints fewer than ten entries. A sample covering under 80% of the
consensus is not scored at all, which is why AeBuster3 starts at two copies
and why 29 of the 3,561 Aedes families carry no scan.

### What `drop` does not account for

`--cover-gap` tops up the kept sample and never a scanned one. Where a family
was topped up, `drop` measures a topped-up consensus against consensi built
without those extra copies, and `hint` can turn on that difference by itself.
It reaches 172 of 3,561 families in Aedes and 176 of 2,701 in Xenopus. For the
rest the kept sample is the ten-copy grid point and the two numbers are
measured the same way.

## Unoccupied columns

A reference position that no sampled copy reaches carries no evidence either
way. Scoring it as disagreement is wrong, and so is dropping it silently.

`generateSeedAlignments.pl` handles them in two different ways, and this tool
matches both by default:

- Outside the copies, before the first and after the last, no column exists.
  Those reference bases are deletions.
- Between copies, the column exists but every row is padding. `MultAln`'s
  consensus caller skips padding, every base ties at zero, and the tie rule
  calls N.

So an N in these seeds means one of two things. `ConsCmp` therefore reports
`amb` for columns a copy covered and `unocc` for columns none did, and gives
`id` over the whole reference alongside `id_occupied` over covered columns
only. An `Occupancy` line lists the unoccupied regions in reference
coordinates.

Gypsy-263_AA-I in Aedes aegypti is the case that motivated this, and makes a
good test: 3 copies cover 239 of 4,892 bases, its identity reads 5.16% across
the reference and 79.50% over covered columns, and the Perl seed carries a
2,171 base run of N.

`--min-non-gap 1` calls unoccupied columns as gaps instead, which shortens the
consensus and moves them into the trim fields. That departs from
`generateSeedAlignments.pl`, so it is off by default.

## Validation

Against the Perl on Aedes aegypti, RepBase 31.08, 3,561 shared families:
median identity 98.95% against 98.99%, differing by 0.02 points as a median,
within one point for 94% of families, with identical length change in 77%.
The closest policy reproduces an independent Python reconstruction to within
one point for all 3,528 families the two share.

Both policies over 3,563 families take about 8 minutes on 24 threads. The Perl
took 3 hours 12 minutes for one policy on 20 threads.

## Caveats

`--target-copies` counts alignments, as the Perl does. Fragments sharing a
RepeatMasker element id stay separate.

dfam-lib 0.2.1's crossmatch parser read RepeatMasker's per-alignment trailer as
sequence, appending text like `20p39g.matrix` to the query row. 0.2.2 fixes it,
and this crate requires that version.

## Usage

```
recover-seeds --align genome.fa.align --library lib.fa --out-prefix out/aedes \
            --assembly-id GCA_002204515.1 --taxon "Aedes aegypti" --threads 24
```

Writes `out/aedes-population.stk`, `out/aedes-closest.stk` and
`out/aedes-summary.tsv`. `--family` builds named families only; `--index-only`
reports the index and stops.

The crate carries its own `[workspace]` table, so it builds without being a
member of the dfam-curator workspace. Fold it in by adding `recover-seeds` to
that manifest's members when it belongs there.
