//! Per-family work: scrub the alignments, pick a sample, build the MSA.
//!
//! The unit is one alignment, as in `generateSeedAlignments.pl`.  Fragments
//! sharing a RepeatMasker element id stay separate, because the Perl treats
//! them separately and joining them would change what a sample holds.

use std::collections::HashSet;
use std::path::Path;

use aln_core::crossmatch::{self, PairwiseHit};
use aln_core::msa::{assemble_msa, InsertionPolicy, MsaMember, MultiAlign};
use aln_core::seq::{revcomp, Strand};
use aln_core::twobit::TwoBitReader;

use anyhow::{Context, Result};

use crate::align_index::RecordRef;

/// One alignment, with its rows read back from the file.
pub struct Instance {
    pub rec: RecordRef,
    pub hit: PairwiseHit,
    /// Rows in the consensus's orientation.
    pub gapped_query: Vec<u8>,
    pub gapped_reference: Vec<u8>,
    pub name: String,
}

impl Instance {
    /// Consensus bases this alignment covers.
    pub fn span(&self) -> u64 {
        self.rec.cons_end - self.rec.cons_start + 1
    }
}

/// Why alignments were dropped before sampling.
#[derive(Debug, Default, Clone)]
pub struct Scrubbed {
    /// Records whose consensus length is not the library's.  RepeatMasker's
    /// own simple-repeat annotations look like this: they carry a per-locus
    /// consensus under the same name as a library entry.
    pub wrong_frame: usize,
    pub too_short: usize,
    pub duplicate: usize,
    pub masked_run: usize,
    pub failed_assembly_check: usize,
    pub unusable_rows: usize,
    /// Records the crossmatch parser rejects.  RepeatMasker writes a few
    /// reverse-strand fragments that align bases to no consensus at all
    /// (consensus end one below start); the parser refuses the empty span.
    pub unparsable: usize,
    /// The parser's message for the first unparsable record.
    pub first_parse_error: Option<String>,
}

impl Scrubbed {
    pub fn total(&self) -> usize {
        self.wrong_frame
            + self.too_short
            + self.duplicate
            + self.masked_run
            + self.failed_assembly_check
            + self.unusable_rows
            + self.unparsable
    }
}

pub struct ScrubParams<'a> {
    /// `-minAlignedLength`: an instance shorter than this is dropped.
    pub min_aligned_length: u64,
    /// Drop an instance carrying this many masked bases in a row.
    pub max_masked_run: usize,
    /// Check every instance against the assembly, as the Perl does.
    pub assembly: Option<&'a TwoBitReader>,
}

/// Read a family's alignments and apply the Perl's scrubbing rules.
pub fn load_instances(
    align_path: &Path,
    family: &str,
    consensus: &[u8],
    records: &[RecordRef],
    assembly_id: Option<&str>,
    p: &ScrubParams,
) -> Result<(Vec<Instance>, Scrubbed)> {
    let text = crate::align_index::read_records(align_path, records)
        .with_context(|| format!("reading alignments for {family}"))?;
    let mut by_offset: Vec<&RecordRef> = records.iter().collect();
    by_offset.sort_by_key(|r| r.offset);

    // read_records concatenates the records in offset order, so each record
    // is the next rec.len bytes.  Parsing them one at a time costs a family
    // one record when the parser rejects it, not all of them, and keeps every
    // hit paired with its own index entry.
    let mut out = Vec::new();
    let mut drop = Scrubbed::default();
    let mut seen: HashSet<(String, u64, u64)> = HashSet::new();
    let mut pos = 0;
    for rec in by_offset {
        let chunk = &text.as_bytes()[pos..pos + rec.len];
        pos += rec.len;
        let mut parsed: Vec<PairwiseHit> = match crossmatch::parse(std::io::Cursor::new(chunk)) {
            Ok(v) => v,
            Err(e) => {
                drop.unparsable += 1;
                if drop.first_parse_error.is_none() {
                    drop.first_parse_error = Some(e.to_string());
                }
                continue;
            }
        };
        if parsed.len() != 1 {
            drop.unusable_rows += 1;
            continue;
        }
        let h = &parsed.pop().unwrap();
        if rec.cons_len as usize != consensus.len() {
            drop.wrong_frame += 1;
            continue;
        }
        if h.subj_seq.is_empty() || h.query_seq.len() != h.subj_seq.len() {
            drop.unusable_rows += 1;
            continue;
        }
        if h.query.len() < p.min_aligned_length {
            drop.too_short += 1;
            continue;
        }
        if !seen.insert((h.query_name.clone(), h.query.start(), h.query.end())) {
            drop.duplicate += 1;
            continue;
        }
        let bases: Vec<u8> = h
            .query_seq
            .iter()
            .filter(|b| b.is_ascii_alphabetic())
            .map(|b| b.to_ascii_uppercase())
            .collect();
        if p.max_masked_run > 0
            && bases
                .windows(p.max_masked_run)
                .any(|w| w.iter().all(|&b| b == b'X'))
        {
            drop.masked_run += 1;
            continue;
        }
        if let Some(tb) = p.assembly {
            // The instance has to match the assembly it came from.  A 2bit
            // file holds only A, C, G, T and N, and faToTwoBit writes every
            // other IUB code as N, so an N there matches any IUB code in the
            // alignment.  Without this, copies from assemblies that carry IUB
            // codes failed the check; in one RepBase 31.08 genome that cost
            // 293 of 314 families at least one of their closest ten copies,
            // often the longest and least diverged.
            let ok = match tb.fetch(&h.query_name, h.query.start(), h.query.end()) {
                Ok(g) => {
                    g.len() == bases.len()
                        && g.iter().zip(&bases).all(|(a, b)| {
                            let a = a.to_ascii_uppercase();
                            a == *b || (a == b'N' && b"RYKMSWBDHV".contains(b))
                        })
                }
                Err(_) => false,
            };
            if !ok {
                drop.failed_assembly_check += 1;
                continue;
            }
        }
        let (mut q, s) = if h.orientation == Strand::Minus {
            (revcomp(&h.query_seq), revcomp(&h.subj_seq))
        } else {
            (h.query_seq.clone(), h.subj_seq.clone())
        };
        // A seed row carries the assembly's bases as its 2bit holds them:
        // an IUB code other than N becomes N.
        let mut h = h.clone();
        for seq in [&mut q, &mut h.query_seq] {
            for b in seq.iter_mut() {
                if b"RYKMSWBDHVrykmswbdhv".contains(b) {
                    *b = if b.is_ascii_lowercase() { b'n' } else { b'N' };
                }
            }
        }
        if h.subj.start() as usize + aln_core::seq::ungapped_len(&s) > consensus.len() {
            drop.unusable_rows += 1;
            continue;
        }
        let strand = if h.orientation == Strand::Minus { '-' } else { '+' };
        let name = match assembly_id {
            Some(a) => format!(
                "{a}:{}:{}-{}_{strand}",
                h.query_name,
                h.query.start() + 1,
                h.query.end()
            ),
            None => format!(
                "{}:{}-{}_{strand}",
                h.query_name,
                h.query.start() + 1,
                h.query.end()
            ),
        };
        out.push(Instance {
            rec: (*rec).clone(),
            hit: h.clone(),
            gapped_query: q,
            gapped_reference: s,
            name,
        });
    }
    Ok((out, drop))
}

/// How to choose the instances that go into a seed alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Policy {
    /// generateSeedAlignments.pl's rule.
    Population,
    /// The instances closest to the library consensus.
    Closest,
}

pub struct SampleParams {
    pub target_copies: usize,
    pub min_depth: u32,
    /// Closest policy: an instance must cover this fraction of the consensus.
    pub min_span_fraction: f64,
}

/// Every instance under the closest rule, nearest first.
///
/// Kimura divergence ascending orders the list, and the sort is stable, so
/// ties keep the order the `.align` file gave them.  Span then splits the
/// list in two rather than sorting it: instances covering `min_span_fraction`
/// of the consensus come first, the shorter ones after, each part still in
/// divergence order.  A short copy is therefore used only once the long ones
/// run out.
pub fn rank_closest(inst: &[Instance], cons_len: u64, min_span_fraction: f64) -> Vec<usize> {
    let mut order: Vec<usize> = (0..inst.len()).collect();
    order.sort_by(|&a, &b| {
        inst[a]
            .rec
            .kimura
            .partial_cmp(&inst[b].rec.kimura)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let need = (min_span_fraction * cons_len as f64) as u64;
    let (mut full, short): (Vec<usize>, Vec<usize>) =
        order.into_iter().partition(|&i| inst[i].span() >= need);
    full.extend(short);
    full
}

/// Choose instances under `policy`, returning indices into `inst`.
///
/// The population rule, from `buildFamily`: rank by Kimura divergence, move
/// the top quartile to the end of the list, order each part by consensus span
/// descending, then walk until the budget is met and every tenth consensus
/// position carries `min_depth` instances.  An instance past the budget is
/// still taken when it covers a position still short of that depth.
pub fn sample(inst: &[Instance], cons_len: u64, policy: Policy, p: &SampleParams) -> Vec<usize> {
    let mut order: Vec<usize> = (0..inst.len()).collect();
    order.sort_by(|&a, &b| {
        inst[a]
            .rec
            .kimura
            .partial_cmp(&inst[b].rec.kimura)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    match policy {
        Policy::Closest => {
            let mut ranked = rank_closest(inst, cons_len, p.min_span_fraction);
            ranked.truncate(p.target_copies);
            ranked
        }
        Policy::Population => {
            let n = order.len();
            let cut = (3 * n) / 4;
            let mut main: Vec<usize> = order[..cut].to_vec();
            let mut outliers: Vec<usize> = order[cut..].to_vec();
            let by_span =
                |v: &mut Vec<usize>| v.sort_by_key(|&i| std::cmp::Reverse(inst[i].span()));
            by_span(&mut main);
            by_span(&mut outliers);
            let budget_end = main.len();
            main.extend(outliers);

            let bins = ((cons_len + 9) / 10) as usize;
            let mut depth = vec![0u32; bins.max(1)];
            let mut low = depth.len();
            let mut chosen = Vec::new();
            for (rank, &i) in main.iter().enumerate() {
                if chosen.len() >= p.target_copies && low == 0 {
                    break;
                }
                let lo = ((inst[i].rec.cons_start + 8) / 10) as usize;
                let hi = (((inst[i].rec.cons_end - 1) / 10) as usize).min(depth.len() - 1);
                if lo > hi {
                    continue;
                }
                let fills_gap = depth[lo..=hi].iter().any(|&d| d < p.min_depth);
                let within_budget = rank < budget_end && chosen.len() < p.target_copies;
                if !(within_budget || fills_gap) {
                    continue;
                }
                for d in depth[lo..=hi].iter_mut() {
                    *d += 1;
                    if *d == p.min_depth {
                        low -= 1;
                    }
                }
                chosen.push(i);
            }
            chosen
        }
    }
}

/// Add instances until no more than `min_gap` consensus bases are left
/// uncovered, taking the closest ones that reach into a gap.
///
/// The closest policy often settles on one or two instances, which on a long
/// family can leave hundreds of bases with nothing over them.  This fills
/// those in without disturbing a sample that already covers the consensus:
/// a family under the threshold gets nothing added.  With `ends_only`, gaps
/// before the first and after the last covered base are filled and interior
/// ones are left alone.
/// Extend `chosen` until no more than `min_gap` consensus bases stay
/// uncovered.  `ranked` carries every instance in the order the policy
/// prefers them, closest first; candidates come from there, so a gap is
/// filled by the nearest copy that reaches it rather than the one that
/// happens to sit earliest in the .align file.
pub fn top_up_coverage(
    inst: &[Instance],
    ranked: &[usize],
    chosen: &[usize],
    cons_len: u64,
    min_gap: u64,
    ends_only: bool,
) -> Vec<usize> {
    let mut covered = vec![false; cons_len as usize];
    let mut mark = |i: usize, covered: &mut Vec<bool>| {
        for p in inst[i].rec.cons_start..=inst[i].rec.cons_end {
            if let Some(c) = covered.get_mut((p - 1) as usize) {
                *c = true;
            }
        }
    };
    let mut out = chosen.to_vec();
    for &i in chosen {
        mark(i, &mut covered);
    }
    // Which uncovered bases count: all of them, or only the two flanks.
    let wanted = |covered: &[bool]| -> Vec<bool> {
        let mut want = covered.iter().map(|c| !c).collect::<Vec<bool>>();
        if ends_only {
            let first = covered.iter().position(|c| *c);
            let last = covered.iter().rposition(|c| *c);
            if let (Some(f), Some(l)) = (first, last) {
                for w in want[f..=l].iter_mut() {
                    *w = false;
                }
            }
        }
        want
    };
    let taken: std::collections::HashSet<usize> = chosen.iter().copied().collect();
    loop {
        let want = wanted(&covered);
        if want.iter().filter(|w| **w).count() as u64 <= min_gap {
            break;
        }
        let next = ranked.iter().copied().find(|i| {
            !taken.contains(i)
                && !out.contains(i)
                && (inst[*i].rec.cons_start..=inst[*i].rec.cons_end)
                    .any(|p| want.get((p - 1) as usize).copied().unwrap_or(false))
        });
        match next {
            Some(i) => {
                mark(i, &mut covered);
                out.push(i);
            }
            None => break,
        }
    }
    out
}

/// Add instances until every consensus position has `min_depth` of them, or
/// no unchosen instance reaches a position still short.
///
/// The closest sample and the coverage fill leave stretches covered by one or
/// two copies, and a consensus called there follows whichever copy happens to
/// cover it, so ties and single-copy calls go against the library even where
/// most copies in the genome carry the library's base.  This walks the same
/// ranking.  Among the instances within `window` Kimura points of the first
/// one that reaches a short position, it takes the one reaching the most short
/// positions (earlier in the ranking on a tie): a long copy rather than a
/// fragment, which gives fewer rows and fewer fragment ends in the alignment.
/// `window` 0 keeps the ranking's order.
pub fn top_up_depth(
    inst: &[Instance],
    ranked: &[usize],
    chosen: &[usize],
    cons_len: u64,
    min_depth: u32,
    window: f64,
) -> Vec<usize> {
    let n = cons_len as usize;
    let span = |i: usize| -> (usize, usize) {
        let a = (inst[i].rec.cons_start.max(1) - 1) as usize;
        let b = (inst[i].rec.cons_end as usize).min(n);
        (a.min(b), b)
    };
    let mut depth = vec![0u32; n];
    let mut out = chosen.to_vec();
    let mut taken = vec![false; inst.len()];
    for &i in chosen {
        let (a, b) = span(i);
        depth[a..b].iter_mut().for_each(|d| *d += 1);
        taken[i] = true;
    }
    loop {
        // short[p] counts the positions before p still under min_depth
        let mut short = vec![0u32; n + 1];
        for p in 0..n {
            short[p + 1] = short[p] + u32::from(depth[p] < min_depth);
        }
        if short[n] == 0 {
            break;
        }
        let reach = |i: usize| -> u32 {
            let (a, b) = span(i);
            short[b] - short[a]
        };
        let Some(first) = ranked.iter().copied().find(|&i| !taken[i] && reach(i) > 0) else {
            break;
        };
        let mut best = first;
        if window > 0.0 {
            let limit = inst[first].rec.kimura + window;
            let mut best_reach = reach(first);
            for &i in ranked {
                if !taken[i] && inst[i].rec.kimura <= limit {
                    let r = reach(i);
                    if r > best_reach {
                        best = i;
                        best_reach = r;
                    }
                }
            }
        }
        let (a, b) = span(best);
        depth[a..b].iter_mut().for_each(|d| *d += 1);
        taken[best] = true;
        out.push(best);
    }
    out
}

/// Merge the chosen instances onto the consensus.
pub fn build_msa(
    family: &str,
    consensus: &[u8],
    inst: &[Instance],
    chosen: &[usize],
    policy: InsertionPolicy,
) -> Result<MultiAlign> {
    let members: Vec<MsaMember> = chosen
        .iter()
        .map(|&i| {
            let h = &inst[i].hit;
            MsaMember {
                name: &inst[i].name,
                gapped_query: &inst[i].gapped_query,
                gapped_reference: &inst[i].gapped_reference,
                ref_start: h.subj.start() as usize,
                span: Some(h.query),
                orient: h.orientation,
            }
        })
        .collect();
    assemble_msa(consensus, family, &members, policy)
        .map_err(|e| anyhow::anyhow!("{family}: {e}"))
}
