//! Genome-anchored extension of a finished consensus, via RAMExtend.
//!
//! Phase 2 refines a consensus against a fixed set of instance sequences, so
//! it can only ever describe the part of the element those sequences already
//! cover. The element usually continues past them. RAMExtend answers the
//! separate question — what lies beyond the edges — by going back to the
//! genome the instances came from and extending every copy simultaneously
//! under a fit-preferred banded model, which is how the Perl `Refiner`
//! invoked it.
//!
//! The flow here mirrors that: run the normal pipeline, extend the result
//! against the genome, then hand the extended consensus and the widened
//! instances back to phase 2 so the new bases are refined rather than merely
//! appended.
//!
//! # Which copies get a say
//!
//! Two different reasons disqualify a copy, and they are not interchangeable:
//!
//!  - **Doesn't reach the edge.** A copy whose alignment stops well inside
//!    the consensus knows nothing about what lies beyond it. RAMExtend has
//!    its own mechanism for this — `left_extendable` / `right_extendable` —
//!    and the copy stays in the alignment, contributing to the other side and
//!    bounded by [`ExtendParams::cap_penalty`]. The threshold is the Perl's:
//!    within [`Config::edge_slop`] columns of the edge.
//!
//!  - **Can't be located in the genome.** A name that isn't a Smitten
//!    identifier, an assembly that doesn't match, a sequence absent from the
//!    2bit, coordinates past its end, or bases that don't match what the
//!    genome holds there. There is no range to anchor, so the copy cannot
//!    take part in the extension at all. It is dropped from the extension and
//!    a warning names it and says why — never silently, and never fatally:
//!    the copy keeps its place in the family and in the refinement that
//!    follows.
//!
//! Both are non-contributions to the *extension*, not exclusions from the
//! family.

use std::collections::HashMap;

use aln_core::msa::MultiAlign;
use aln_core::seq::{Sequence, Strand};
use aln_core::stats::{kimura_divergence, Masking};
use aln_core::twobit::TwoBitReader;
use aln_coord::Span;
use anyhow::{Context, Result};
use ram_core::alphabet::num_to_char;
use ram_core::engine::{
    apply_overlap_avoidance, extend_alignment, Direction, ExtendParams, ScoreArena,
};
use ram_core::library::{load_sequence_subset_minimal, RangeRecord};
use ram_core::matrix::ScoringSystem;
use smitten::Identifier;

/// How to run the extension.
#[derive(Debug, Clone)]
pub struct Config {
    /// Expected assembly identifier. When set, an input whose Smitten
    /// identifier names a different assembly is rejected; when unset, any
    /// assembly prefix is accepted and only the sequence name is matched.
    pub assembly: Option<String>,
    /// Maximum extension per side (the RAMExtend `L`).
    pub l_max: i32,
    /// Refuse an extension whose two sides together exceed this many bases.
    ///
    /// A cap on how much an extension may add regardless of how confident the
    /// model is. Reaching it means the copies agree over a span far larger
    /// than the element that was discovered, which is a satellite, a
    /// segmental duplication, or a nested insertion — not a boundary.
    pub max_total: i32,
    pub bandwidth: i32,
    /// Extendable copies needed before extension is attempted, and the
    /// multiplier in the adaptive `min_improvement`.
    pub min_aligning_seqs: i32,
    /// Overrides the divergence-adaptive matrix choice.
    pub matrix: Option<String>,
    /// Overrides the adaptive `min_aligning_seqs x diagonal_average`.
    pub min_improvement: Option<i32>,
    pub cap_penalty: i32,
    pub when_to_stop: i32,
    pub max_occurrences: usize,
    /// Columns of slack allowed between the alignment edge and the consensus
    /// edge before a copy stops being extendable on that side.
    pub edge_slop: usize,
    /// Fraction of comparable bases that must match the genome for an input
    /// sequence to be accepted as genuinely from these coordinates.
    pub min_identity: f64,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            assembly: None,
            // The Stockholm-mode defaults from `ram-extend -stk`, which is the
            // closest existing caller to this one.
            l_max: 20_000,
            max_total: 25_000,
            bandwidth: 40,
            min_aligning_seqs: 3,
            matrix: None,
            min_improvement: None,
            cap_penalty: -15,
            when_to_stop: 100,
            max_occurrences: 10_000,
            edge_slop: 10,
            min_identity: 0.95,
        }
    }
}

/// What the extension produced.
pub struct Outcome {
    /// Bases to prepend to the consensus, in consensus orientation.
    pub left: Vec<u8>,
    /// Bases to append.
    pub right: Vec<u8>,
    /// One entry per input sequence, in input order: re-fetched at the
    /// extended range for copies that took part, unchanged for the rest.
    pub seqs: Vec<Sequence>,
    /// Copies that anchored to the genome and were offered to the extension.
    pub contributors: usize,
    /// Divergence and matrix actually used, for the run report.
    pub divergence: f64,
    pub matrix_name: String,
    /// Each contributing copy at its core boundary, captured *before* the
    /// extension runs — the view of what the extension is about to reason
    /// over.
    pub edges: Vec<ram_core::library::CoreEdge>,
}

impl Outcome {
    pub fn extended(&self) -> bool {
        !self.left.is_empty() || !self.right.is_empty()
    }

    /// Left extension + old consensus + right extension.
    pub fn splice(&self, consensus: &[u8]) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.left.len() + consensus.len() + self.right.len());
        out.extend_from_slice(&self.left);
        out.extend_from_slice(consensus);
        out.extend_from_slice(&self.right);
        out
    }
}

/// The Perl's divergence buckets: matrix name and its diagonal average.
///
/// Kept identical to `ram-extend -stk` (which ports `extend-stk.pl`) so a
/// family extended through either path gets the same scoring system.
fn matrix_for_divergence(div: f64) -> (&'static str, i32) {
    if div >= 22.5 {
        ("25p43g", 9)
    } else if div >= 19.0 {
        ("20p43g", 10)
    } else if div >= 16.0 {
        ("18p43g", 10)
    } else {
        ("14p43g", 10)
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            other => other,
        })
        .collect()
}

/// Fraction of positions where both sides hold an unambiguous base and they
/// agree. `None` when the lengths differ or nothing is comparable.
///
/// Ns are skipped rather than counted as mismatches: an instance drawn from an
/// assembly gap would otherwise fail validation for a reason that says nothing
/// about whether the coordinates are right.
fn identity(a: &[u8], b: &[u8]) -> Option<f64> {
    if a.len() != b.len() || a.is_empty() {
        return None;
    }
    let (mut comparable, mut same) = (0usize, 0usize);
    for (&x, &y) in a.iter().zip(b.iter()) {
        let (x, y) = (x.to_ascii_uppercase(), y.to_ascii_uppercase());
        if !matches!(x, b'A' | b'C' | b'G' | b'T') || !matches!(y, b'A' | b'C' | b'G' | b'T') {
            continue;
        }
        comparable += 1;
        if x == y {
            same += 1;
        }
    }
    if comparable == 0 {
        None
    } else {
        Some(same as f64 / comparable as f64)
    }
}

/// A copy that anchored to the genome, with its core range and edge flags.
struct Anchor {
    /// Index into the caller's input sequences.
    seq_index: usize,
    /// Sequence name within the assembly.
    chrom: String,
    /// Genomic core range (the aligned extent).
    span: Span,
    minus: bool,
    left_extendable: bool,
    right_extendable: bool,
    /// Assembly prefix to put back on the re-fetched name, if the input had one.
    assembly: Option<String>,
    /// The input identifier, for diagnostics.
    source: String,
}

/// Aligned extent and edge reach of one input sequence across all its rows.
struct RowSpan {
    /// Envelope of the copy's rows on the input sequence. `None` when any of
    /// its rows carried no coordinates, so the copy cannot be placed.
    span: Option<Span>,
    orient: Strand,
    reaches_left: bool,
    reaches_right: bool,
    /// Width of the widest row, used to pick the orientation when a copy
    /// contributed rows on both strands.
    best_span: u64,
}

/// Collapse the MSA's instance rows to one span per input sequence.
///
/// With `--hsp all` one copy can contribute several rows. The core is their
/// envelope, and the copy reaches an edge if any of its rows does.
fn row_spans(
    msa: &MultiAlign,
    by_name: &HashMap<&str, usize>,
    edge_slop: usize,
) -> HashMap<usize, RowSpan> {
    let width = msa.width();
    let mut out: HashMap<usize, RowSpan> = HashMap::new();
    // Row 0 is the reference (the consensus), not a copy.
    for row in msa.sequences.iter().skip(1) {
        let Some(&idx) = by_name.get(row.name.as_str()) else {
            continue;
        };
        let reaches_left = row.col_start <= edge_slop;
        let reaches_right = row.col_end + edge_slop >= width;
        let span = row.span.map_or(0, |s| s.len());
        out.entry(idx)
            .and_modify(|e| {
                e.span = match (e.span, row.span) {
                    (Some(a), Some(b)) => Some(
                        Span::new(a.start().min(b.start()), a.end().max(b.end()))
                            .expect("an envelope of ascending spans is ascending"),
                    ),
                    _ => None,
                };
                e.reaches_left |= reaches_left;
                e.reaches_right |= reaches_right;
                if span > e.best_span {
                    e.best_span = span;
                    e.orient = row.orient;
                }
            })
            .or_insert(RowSpan {
                span: row.span,
                orient: row.orient,
                reaches_left,
                reaches_right,
                best_span: span,
            });
    }
    out
}

/// Place a sub-range of an oriented parent range back on the genome.
///
/// This is exactly what a recursive Smitten identifier means, so it is left to
/// Smitten: `chr1:1000-2000_-` with an aligned extent of `10-50_+` is written
/// `chr1:1000-2000_-:10-50_+` and normalized to a single absolute range. The
/// coordinates flip end-for-end on a reverse parent and the strand is the XOR
/// of the two.
///
/// Requires Smitten >= 1.0.2. v1.0.1-alpha returned `-` for every `+` parent
/// regardless of the sub-range's strand — right coordinates, wrong strand, on
/// the commonest composition there is — which pointed the extension at the
/// opposite end of the element. [`tests::smitten_composition_truth_table`]
/// pins the behaviour so a dependency bump cannot quietly reintroduce it.
fn compose_range(
    parent: (usize, usize, char),
    sub: (usize, usize, char),
) -> std::result::Result<(usize, usize, char), String> {
    let composed = format!(
        "seq:{}-{}_{}:{}-{}_{}",
        parent.0, parent.1, parent.2, sub.0, sub.1, sub.2
    );
    let normalized = Identifier::from_v2(&composed)
        .and_then(|i| i.normalize())
        .map_err(|e| format!("could not resolve {composed:?}: {e}"))?;
    let r = normalized
        .ranges
        .first()
        .ok_or_else(|| format!("{composed:?} normalized to no range"))?;
    Ok((r.start, r.end, r.orientation))
}

/// Resolve one input sequence to a genomic core range, or say why not.
///
/// The identifier is composed with the aligned sub-range and normalized by
/// Smitten rather than by hand: `chr1:1000-2000_-` aligned over its own
/// positions 10-50 is `chr1:1951-1991_+`, and that arithmetic is exactly what
/// Smitten's recursive-range normalization exists to get right.
fn anchor(
    seq_index: usize,
    seq: &Sequence,
    span: &RowSpan,
    genome: &TwoBitReader,
    cfg: &Config,
) -> std::result::Result<Anchor, String> {
    let (ident, _version) = Identifier::from_unknown_format(&seq.name, false, false)
        .map_err(|e| format!("not a Smitten identifier ({e})"))?;
    let ident = ident
        .normalize()
        .map_err(|e| format!("identifier ranges could not be normalized ({e})"))?;

    match ident.ranges.len() {
        1 => {}
        0 => {
            return Err(
                "identifier carries no genomic range; extension needs \
                 `sequence:start-end_orient`"
                    .to_string(),
            )
        }
        n => {
            return Err(format!(
                "identifier resolves to {n} ranges after normalization; a copy \
                 must be one contiguous genomic range"
            ))
        }
    }
    if let (Some(want), Some(got)) = (cfg.assembly.as_deref(), ident.assembly_id.as_deref()) {
        if want != got {
            return Err(format!("assembly {got:?} does not match {want:?}"));
        }
    }
    let chrom = ident.sequence_id.clone();
    let seq_len = genome
        .seq_len(&chrom)
        .ok_or_else(|| format!("sequence {chrom:?} is not in the genome"))?
        as i64;

    // Whole-sequence validation, before any sub-range arithmetic: does this
    // input actually hold what the genome holds at these coordinates?
    let whole = &ident.ranges[0];
    let (w_start, w_end) = (whole.start as i64 - 1, whole.end as i64);
    if w_start < 0 || w_end > seq_len || w_start >= w_end {
        return Err(format!(
            "range {}-{} lies outside {chrom} (length {seq_len})",
            whole.start, whole.end
        ));
    }
    let mut genomic = genome
        .fetch(&chrom, w_start as u64, w_end as u64)
        .map_err(|e| format!("could not read {chrom}:{}-{}: {e}", whole.start, whole.end))?;
    if whole.orientation == '-' {
        genomic = revcomp(&genomic);
    }
    match identity(&genomic, &seq.seq) {
        None => {
            return Err(format!(
                "sequence is {} bp but its range {}-{} is {} bp",
                seq.seq.len(),
                whole.start,
                whole.end,
                genomic.len()
            ))
        }
        Some(id) if id < cfg.min_identity => {
            return Err(format!(
                "only {:.1}% identical to {chrom}:{}-{}_{} in the genome \
                 (need {:.1}%)",
                id * 100.0,
                whole.start,
                whole.end,
                whole.orientation,
                cfg.min_identity * 100.0
            ))
        }
        Some(_) => {}
    }

    // Core = the aligned extent, placed back on the genome.
    //
    // The strand that comes out is the copy's orientation relative to the
    // consensus, which is what the extension needs: it walks every core
    // outward in consensus space, so `-` has to mean "this copy runs the
    // other way", not "this copy is on the genome's minus strand".
    let orient_ch = if span.orient == Strand::Minus { '-' } else { '+' };
    let Some((start_1b, end_1b)) = span.span.and_then(|s| s.as_1b_closed()) else {
        return Err("aligned extent is unknown: a row carried no coordinates".to_string());
    };
    if end_1b as usize > seq.seq.len() {
        return Err(format!(
            "aligned extent ends at {} but the sequence is {} bp",
            end_1b,
            seq.seq.len()
        ));
    }
    // compose_range works in Smitten's 1-based closed terms on both sides.
    let (c_start, c_end, c_orient) = compose_range(
        (whole.start, whole.end, whole.orientation),
        (start_1b as usize, end_1b as usize, orient_ch),
    )?;

    Ok(Anchor {
        seq_index,
        chrom,
        span: Span::from_1b_closed(c_start as u64, c_end as u64)
            .map_err(|e| format!("composed range {c_start}-{c_end}: {e}"))?,
        minus: c_orient == '-',
        left_extendable: span.reaches_left,
        right_extendable: span.reaches_right,
        assembly: ident.assembly_id.clone(),
        source: seq.name.clone(),
    })
}

/// Mean unadjusted Kimura divergence of the instance rows against the
/// reference row, which selects the scoring matrix.
fn msa_divergence(msa: &MultiAlign) -> f64 {
    let Some(reference) = msa.reference_seq() else {
        return 0.0;
    };
    let mut divs = Vec::new();
    for row in msa.sequences.iter().skip(1) {
        if row.seq.len() != reference.len() {
            continue;
        }
        if let Ok(d) = kimura_divergence(&row.seq, reference, false, Masking::Ignore) {
            if let Some(v) = d.value {
                divs.push(v);
            }
        }
    }
    if divs.is_empty() {
        0.0
    } else {
        divs.iter().sum::<f64>() / divs.len() as f64
    }
}

/// Extend `msa`'s consensus against the genome.
///
/// `Ok(None)` means the extension was not attempted or produced nothing —
/// too few copies could be anchored, or the model declined to extend either
/// way. Warnings are returned rather than printed so the caller controls
/// where they go.
pub fn extend(
    genome: &TwoBitReader,
    cfg: &Config,
    msa: &MultiAlign,
    seqs: &[Sequence],
) -> Result<(Option<Outcome>, Vec<String>)> {
    let mut warnings = Vec::new();

    let mut by_name: HashMap<&str, usize> = HashMap::new();
    for (i, s) in seqs.iter().enumerate() {
        if by_name.insert(s.name.as_str(), i).is_some() {
            warnings.push(format!(
                "duplicate sequence name {:?}; extension will anchor only the last",
                s.name
            ));
        }
    }

    let spans = row_spans(msa, &by_name, cfg.edge_slop);
    let mut anchors: Vec<Anchor> = Vec::new();
    // Input order, so warnings and the extension are reproducible.
    for (i, seq) in seqs.iter().enumerate() {
        let Some(span) = spans.get(&i) else { continue };
        if !(span.reaches_left || span.reaches_right) {
            continue; // reaches neither edge: nothing to say about the flanks
        }
        match anchor(i, seq, span, genome, cfg) {
            Ok(a) => anchors.push(a),
            Err(why) => warnings.push(format!(
                "{}: not contributing to the extension — {why}",
                seq.name
            )),
        }
    }

    if anchors.len() <= cfg.min_aligning_seqs as usize {
        warnings.push(format!(
            "only {} copies could be anchored to the genome ({} needed); \
             skipping extension",
            anchors.len(),
            cfg.min_aligning_seqs as usize + 1
        ));
        return Ok((None, warnings));
    }
    if anchors.len() > cfg.max_occurrences {
        warnings.push(format!(
            "{} anchored copies exceeds --extend-max-occurrences {}; \
             skipping extension",
            anchors.len(),
            cfg.max_occurrences
        ));
        return Ok((None, warnings));
    }

    let divergence = msa_divergence(msa);
    let (adaptive_matrix, diag_avg) = matrix_for_divergence(divergence);
    let matrix_name = cfg.matrix.clone().unwrap_or_else(|| adaptive_matrix.into());
    let scoring = ScoringSystem::by_name(&matrix_name)
        .map_err(|e| anyhow::anyhow!("unknown extension matrix {matrix_name:?}: {e}"))?;
    let min_improvement = cfg
        .min_improvement
        .unwrap_or(cfg.min_aligning_seqs * diag_avg);

    // `load_sequence_subset_minimal` sorts internally (name asc, start asc,
    // end desc) and returns cores in that order. Sorting to match first — the
    // sort is stable — keeps `cores[i]` paired with `anchors[i]`, which is how
    // each copy's extension length finds its way back to its sequence.
    anchors.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.span.start().cmp(&b.span.start()))
            .then(b.span.end().cmp(&a.span.end()))
    });
    if std::env::var_os("TE_COMPOSER_EXTEND_DEBUG").is_some() {
        for a in &anchors {
            eprintln!(
                "ANCHOR\t{}\t{}\t{}\t{}\tleft={}\tright={}\tfrom={}",
                a.chrom, a.span.start(), a.span.end(),
                if a.minus { '-' } else { '+' },
                a.left_extendable, a.right_extendable, a.source
            );
        }
    }
    let ranges: Vec<RangeRecord> = anchors
        .iter()
        .map(|a| RangeRecord {
            name: a.chrom.clone(),
            span: a.span,
            left_flag: i32::from(a.left_extendable),
            right_flag: i32::from(a.right_extendable),
            minus: a.minus,
        })
        .collect();

    let (lib, mut cores) = load_sequence_subset_minimal(
        genome,
        &ranges,
        (cfg.l_max + cfg.bandwidth) as i64,
        false,
    )
    .map_err(|e| anyhow::anyhow!("loading flanking sequence for extension: {e}"))?;
    anyhow::ensure!(
        cores.len() == anchors.len(),
        "extension loader returned {} cores for {} ranges",
        cores.len(),
        anchors.len()
    );

    // Taken before extension: afterwards the cores carry their extension
    // lengths and no longer describe the starting position.
    let edges = ram_core::library::core_edges(&lib, &cores);

    let l = cfg.l_max as usize;
    let mut master = vec![0u8; 2 * l + 2];
    master[l] = ram_core::alphabet::SYM_N;
    let params = ExtendParams {
        bandwidth: cfg.bandwidth,
        cap_penalty: cfg.cap_penalty,
        min_improvement,
        l_max: cfg.l_max,
        when_to_stop: cfg.when_to_stop,
    };
    let mut arena = ScoreArena::new(cores.len(), cfg.bandwidth);

    // Right first, then overlap avoidance, then left — the order RAMExtend
    // requires: the left pass has to see the bounds the right pass consumed,
    // or two copies of a tandem array extend into each other.
    let right = extend_alignment(
        Direction::Right,
        &mut cores,
        &mut arena,
        &lib,
        &mut master,
        &params,
        &scoring,
    );
    if right.hit_limit {
        warnings.push(format!("extended right to the limit (L={})", cfg.l_max));
    }
    for _ in apply_overlap_avoidance(&mut cores, &lib) {
        // Events are per-pair bookkeeping; the interesting part is the bound
        // change, which the left pass already respects.
    }
    let left = extend_alignment(
        Direction::Left,
        &mut cores,
        &mut arena,
        &lib,
        &mut master,
        &params,
        &scoring,
    );
    if left.hit_limit {
        warnings.push(format!("extended left to the limit (L={})", cfg.l_max));
    }

    // Hitting the cap in BOTH directions means the model never found an edge:
    // the copies go on agreeing for as far as we let them look. That is a
    // satellite array or a segmental duplication, not a transposable element
    // with boundaries, and the "extension" is an arbitrary `l_max` bases of
    // whatever the copies share. Refiner refuses it for the same reason
    // ("Extension hit limits in both directions...probably a segmental
    // duplication, keeping unextended") and so do we.
    //
    // It is also a runaway: an accepted double-cap extension turns a 204 bp
    // family into a 40 kb consensus and then asks phase 2 to refine every
    // instance against it, which on tandem-repeat sequence costs hours and
    // more than a gigabyte for an answer that was wrong to begin with.
    // Refuse the extension if EITHER side ran out of room, or if the two
    // together add more than `max_total`.
    //
    // A capped side exhausted `l_max` rows rather than converging: it was
    // stopped, not finished, so it never found an edge. Nothing about the
    // sequence justifies where it happens to end, and at 20 kb it is well past
    // any real element — full-length ERVs run to ~9-10 kb — so a side still
    // going at the cap is tracking an array or a duplication, not a boundary.
    //
    // Each side is judged on its own. Refusing the whole extension because one
    // side ran away also discards the other side's answer, and on mouse that
    // cost 190 kb of extension across 35 families which annotates at 61% —
    // i.e. real element, thrown away for its neighbour's sins. A side that
    // converged found an edge whatever the other side did.
    let drop_left = left.capped;
    let drop_right = right.capped;
    let kept_left = if drop_left { 0 } else { left.bp };
    let kept_right = if drop_right { 0 } else { right.bp };

    if drop_left {
        warnings.push(format!(
            "left extension hit the {} bp limit — no edge found on that side, \
dropping it",
            cfg.l_max
        ));
    }
    if drop_right {
        warnings.push(format!(
            "right extension hit the {} bp limit — no edge found on that side, \
dropping it",
            cfg.l_max
        ));
    }
    if kept_left + kept_right > cfg.max_total {
        warnings.push(format!(
            "extension would add {} bp, over the {} bp limit; keeping the \
family unextended",
            kept_left + kept_right,
            cfg.max_total
        ));
        return Ok((None, warnings));
    }
    if kept_left == 0 && kept_right == 0 {
        warnings.push(
            "neither side of the extension found an edge; keeping the family \
unextended"
                .to_string(),
        );
        return Ok((None, warnings));
    }

    // Zero the dropped side on every core as well: the per-copy re-fetch below
    // widens each instance by its own extension lengths, and widening into
    // sequence the consensus no longer carries would leave every row hanging
    // off the end of the alignment.
    for c in cores.iter_mut() {
        if drop_left {
            c.left_extension_len = 0;
        }
        if drop_right {
            c.right_extension_len = 0;
        }
    }

    let left_bases: Vec<u8> = master[l - kept_left.max(0) as usize..l]
        .iter()
        .map(|&z| num_to_char(z))
        .collect();
    let right_bases: Vec<u8> = master[l + 1..l + 1 + kept_right.max(0) as usize]
        .iter()
        .map(|&z| num_to_char(z))
        .collect();

    // Re-fetch every contributing copy at its own extended range. Copies
    // extend by different amounts — a neighbour or a scaffold end can stop
    // one short — so this is per-copy, not the consensus-wide totals.
    let mut out_seqs: Vec<Sequence> = seqs.to_vec();
    for (a, core) in anchors.iter().zip(cores.iter()) {
        let (before, after) = if a.minus {
            (core.right_extension_len as i64, core.left_extension_len as i64)
        } else {
            (core.left_extension_len as i64, core.right_extension_len as i64)
        };
        let seq_len = genome.seq_len(&a.chrom).unwrap_or(0) as i64;
        let start = (a.span.start() as i64 - before).max(0);
        let end = (a.span.end() as i64 + after).min(seq_len);
        if start >= end {
            continue;
        }
        let mut bases = genome
            .fetch(&a.chrom, start as u64, end as u64)
            .with_context(|| format!("re-reading {}:{start}-{end}", a.chrom))?;
        if a.minus {
            bases = revcomp(&bases);
        }
        let name = match &a.assembly {
            Some(asm) => format!(
                "{asm}:{}:{}-{}_{}",
                a.chrom,
                start + 1,
                end,
                if a.minus { '-' } else { '+' }
            ),
            None => format!(
                "{}:{}-{}_{}",
                a.chrom,
                start + 1,
                end,
                if a.minus { '-' } else { '+' }
            ),
        };
        out_seqs[a.seq_index] = Sequence::new(name, bases);
    }

    let outcome = Outcome {
        left: left_bases,
        right: right_bases,
        seqs: out_seqs,
        contributors: anchors.len(),
        divergence,
        matrix_name,
        edges,
    };
    Ok((Some(outcome), warnings))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every parent/sub strand combination, plus the offsets on each side.
    ///
    /// The strand column is an XOR and the coordinates flip end-for-end on a
    /// `-` parent; getting either wrong points the extension at the opposite
    /// end of the element, which is how the Smitten bug was found. These now
    /// assert on Smitten itself, so they are a guard on the dependency: a
    /// version that regresses fails here rather than silently in the field.
    #[test]
    fn smitten_composition_truth_table() {
        // (parent, sub) -> expected
        let cases = [
            // Forward parent: sub offsets add from the low end.
            (((1001, 1100, '+'), (1, 100, '+')), (1001, 1100, '+')),
            (((1001, 1100, '+'), (11, 20, '+')), (1011, 1020, '+')),
            (((1001, 1100, '+'), (11, 20, '-')), (1011, 1020, '-')),
            (((1001, 1100, '+'), (1, 10, '+')), (1001, 1010, '+')),
            (((1001, 1100, '+'), (91, 100, '+')), (1091, 1100, '+')),
            // Reverse parent: position 1 is the high genomic coordinate.
            (((1001, 1100, '-'), (11, 20, '+')), (1081, 1090, '-')),
            (((1001, 1100, '-'), (11, 20, '-')), (1081, 1090, '+')),
            (((1001, 1100, '-'), (1, 10, '+')), (1091, 1100, '-')),
            (((1001, 1100, '-'), (91, 100, '+')), (1001, 1010, '-')),
        ];
        for ((parent, sub), want) in cases {
            let got = compose_range(parent, sub).expect("composition failed");
            assert_eq!(
                got, want,
                "{}:{}-{}_{} composed with {}-{}_{}",
                "chr", parent.0, parent.1, parent.2, sub.0, sub.1, sub.2
            );
        }
    }

    /// Smitten's own documented multi-range example, which its `normalize`
    /// gets right: `chr1:100-200_+:10-50_-:1-5_+` is `chr1:145-149_-`.
    /// Composing left-to-right has to reach the same answer.
    #[test]
    fn compose_range_matches_smitten_documented_example() {
        let first = compose_range((100, 200, '+'), (10, 50, '-')).unwrap();
        assert_eq!(first, (109, 149, '-'));
        let second = compose_range((first.0, first.1, first.2), (1, 5, '+')).unwrap();
        assert_eq!(second, (145, 149, '-'));
    }

    /// A whole-span sub-range is the identity on a forward parent. This is
    /// the case Smitten v1.0.1-alpha returns `-` for.
    #[test]
    fn compose_range_forward_identity_stays_forward() {
        assert_eq!(
            compose_range((48389, 48491, '+'), (1, 103, '+')).unwrap(),
            (48389, 48491, '+')
        );
    }

    #[test]
    fn identity_skips_ambiguous_positions() {
        // Ns on either side are not evidence either way.
        assert_eq!(identity(b"ACGT", b"ACGT"), Some(1.0));
        assert_eq!(identity(b"ACNT", b"ACGT"), Some(1.0));
        assert_eq!(identity(b"ACGT", b"ACGA"), Some(0.75));
        assert_eq!(identity(b"NNNN", b"ACGT"), None);
        assert_eq!(identity(b"ACG", b"ACGT"), None);
    }

    #[test]
    fn revcomp_round_trips() {
        assert_eq!(revcomp(b"ACGTN"), b"NACGT".to_vec());
        assert_eq!(revcomp(&revcomp(b"ACGTTGCA")), b"ACGTTGCA".to_vec());
    }

    /// The divergence buckets that pick the extension matrix, kept in step
    /// with `ram-extend -stk`.
    #[test]
    fn matrix_buckets() {
        assert_eq!(matrix_for_divergence(0.0).0, "14p43g");
        assert_eq!(matrix_for_divergence(15.9).0, "14p43g");
        assert_eq!(matrix_for_divergence(16.0).0, "18p43g");
        assert_eq!(matrix_for_divergence(19.0).0, "20p43g");
        assert_eq!(matrix_for_divergence(22.5).0, "25p43g");
        assert_eq!(matrix_for_divergence(22.5).1, 9);
    }
}
