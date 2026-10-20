//! Build seed alignments from a RepeatMasker .align file under two sampling
//! policies, and report how each consensus compares with the library's.
//!
//! This tool never realigns copies. The .align file is its only source of
//! alignments, so a change in output comes from the sampling rather than from
//! a different aligner.

mod subfamily;
mod search;
mod align_index;
mod compare;
mod family;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use rayon::prelude::*;

use aln_core::matrix::SubstMatrix;
use cons_core::Params as ConsParams;
use dfam_stk_io::{SeqRow, StkRecord};

use aln_core::msa::InsertionPolicy;
use family::{Policy, SampleParams, ScrubParams};

#[derive(Parser, Debug)]
#[command(about = "Sample copies from a RepeatMasker .align file and build seed alignments")]
struct Args {
    /// RepeatMasker .align file
    #[arg(long)]
    align: PathBuf,
    /// Library FASTA that RepeatMasker searched with
    #[arg(long)]
    library: PathBuf,
    /// Prefix for the output files
    #[arg(long)]
    out_prefix: PathBuf,
    /// Assembly identifier to prefix instance names with
    #[arg(long)]
    assembly_id: Option<String>,
    /// Assembly 2bit; every instance is checked against it, as the Perl does
    #[arg(long)]
    assembly: Option<PathBuf>,
    /// Drop instances shorter than this (-minAlignedLength)
    #[arg(long, default_value_t = 50)]
    min_aligned_length: u64,
    /// Drop instances with this many masked bases in a row
    #[arg(long, default_value_t = 10)]
    max_masked_run: usize,
    /// Taxon written to each record's OC line
    #[arg(long)]
    taxon: Option<String>,
    /// Copies to aim for under the population policy
    #[arg(long, default_value_t = 200)]
    target_copies: usize,
    /// Copies to take under the closest policy, or fewer when fewer instances
    /// qualify.  0 searches --closest-grid for the count that reproduces the
    /// library consensus best, which lands on one or two copies for about half
    /// the families and leaves a thin sample.
    #[arg(long, default_value_t = 10)]
    closest_copies: usize,
    /// Copy counts the search tries
    #[arg(long, value_delimiter = ',', default_value = "1,2,3,5,8,10,15,20,30,50")]
    closest_grid: Vec<usize>,

    /// Skip the source-sample scan, which scores every --closest-grid size
    /// against the library consensus to show how small a sample explains it.
    /// The scan costs most of the closest policy's runtime.
    #[arg(long)]
    no_source_scan: bool,
    /// Keep adding copies until every sampled position has this depth
    #[arg(long, default_value_t = 10)]
    min_depth: u32,
    /// A closest-policy copy covering this fraction of the consensus ranks
    /// ahead of shorter ones, which fill the target when too few reach it
    #[arg(long, default_value_t = 0.5)]
    min_span_fraction: f64,
    /// Only these families (repeatable)
    #[arg(long)]
    family: Vec<String>,
    /// Worker threads; 0 uses every core
    #[arg(long, default_value_t = 0)]
    threads: usize,
    /// A consensus column needs this many copies with a residue, or it is a
    /// gap.  0 is generateSeedAlignments' behaviour: columns no copy covers
    /// keep the reference's frame and are called N.
    #[arg(long, default_value_t = 0)]
    min_non_gap: usize,
    /// Closest policy: add instances until no more than this many consensus
    /// bases are uncovered.  0 leaves the sample as chosen.
    #[arg(long, default_value_t = 20)]
    cover_gap: u64,
    /// Fill only the gaps before the first and after the last covered base
    #[arg(long)]
    cover_ends_only: bool,
    /// Closest policy: after the coverage fill, keep adding copies in the
    /// closest ranking until every position has this many, or no copy left
    /// reaches a short position.  0 turns it off.
    #[arg(long, default_value_t = 0)]
    closest_min_depth: u32,
    /// With --closest-min-depth: among copies within this many Kimura points of
    /// the next in the ranking, take the one reaching the most short
    /// positions.  0 keeps the ranking's order.
    #[arg(long, default_value_t = 0.0)]
    fill_kimura_window: f64,
    /// Closest policy: keep the plain closest sample even when it leaves
    /// disputed sites unsupported.  By default such a family also gets a
    /// "closest-by-disputed-sites" sample (see disputed_fallback), which
    /// replaces the plain one when it supports more of the library.
    #[arg(long)]
    no_disputed_fallback: bool,
    /// Closest policy: also try a "closest-by-lineage" sample (see
    /// lineage_fallback) and keep it when its consensus is nearer the
    /// library's.  Off by default while it is being evaluated.
    #[arg(long)]
    lineage_fallback: bool,
    /// Lineage fallback: the copies' agreement with the library at the
    /// informative sites must spread at least this much more than chance
    /// (observed variance over binomial variance)
    #[arg(long, default_value_t = 2.0)]
    lineage_min_dispersion: f64,
    /// Lineage fallback: at least this many copies must carry the library's
    /// base at 80% or more of the informative sites they reach
    #[arg(long, default_value_t = 5)]
    lineage_min_copies: usize,
    /// Use the closest-by-lineage sample whenever one can be built (the
    /// gates above still apply), not only when its consensus is nearer the
    /// library's: lineage selection as a policy of its own.  Implies
    /// --lineage-fallback.
    #[arg(long)]
    lineage_always: bool,
    /// Closest policy: also try a "closest-by-repbase-sites" sample (see
    /// rbsites_fallback) and keep it when its consensus is nearer the
    /// library's.  Off by default while it is being evaluated.
    #[arg(long)]
    rbsites_fallback: bool,
    /// RepBase-sites fallback: weight each site by -ln of the share of copies
    /// carrying the library's base there, so rare matches count most
    #[arg(long)]
    rbsites_weighted: bool,
    /// RepBase-sites fallback: copies that must reach a disputed site for it
    /// to count
    #[arg(long, default_value_t = 3)]
    rbsites_min_cover: usize,
    /// RepBase-sites fallback: at least this many copies must score 0.8 or more
    #[arg(long, default_value_t = 3)]
    rbsites_min_copies: usize,
    /// RepBase-sites fallback: take the copies scoring 0.8 or more in the
    /// closest ranking's order (spanning half the consensus first, then lower
    /// Kimura divergence), then the rest by score, instead of all by score.
    /// Ranking by score alone favours short copies that match every site they
    /// reach over long ones that match most.
    #[arg(long)]
    rbsites_rank_order: bool,
    /// With --rbsites-rank-order: when fewer copies are RepBase-like than the
    /// sample takes, fill the rest in the closest ranking's order, not by
    /// score
    #[arg(long)]
    rbsites_pad_rank: bool,
    /// With --rbsites-rank-order: when at least as many copies are
    /// RepBase-like as the sample takes, take the highest-scoring of those
    /// spanning at least --min-span-fraction of the consensus (ties to the
    /// earlier in the closest ranking), then any other RepBase-like copies in
    /// ranking order.  With many to choose from, ranking order alone drifts
    /// back toward the plain closest sample.
    #[arg(long)]
    rbsites_strict_many: bool,
    /// Use the RepBase-sites sample whenever one can be built, not only when
    /// its consensus is nearer the library's.  Implies --rbsites-fallback.
    #[arg(long)]
    rbsites_always: bool,
    /// Closest policy: every fallback sample competes under one rule.  The
    /// sample kept is the one whose consensus leaves the fewest differences
    /// from the library as edit identity counts them (substitutions,
    /// ambiguity codes, insertions, deletions, uncovered library bases; see
    /// edits()) among those leaving no more library bases uncovered than the
    /// plain closest sample.  Without it the disputed-sites sample is
    /// judged on unsupported substitutions and each later sample against the
    /// one kept before it.
    #[arg(long)]
    best_of: bool,
    /// Consensus caller for every consensus recover-seeds calls (the seed's
    /// own, the fallback candidates', the source scan's): dfam (the Dfam
    /// caller with its CpG pass, as before), dfam-nocpg (its CpG pass off),
    /// giri (GIRI's acons caller, which restores nothing) or giri-cpg (GIRI
    /// plus its statistical CpG restoration, acons --mam).  RepBase was built
    /// with GIRI's tools, so scoring a seed against RepBase may call for
    /// GIRI's conventions.  The seed rows are the same whichever is chosen.
    #[arg(long, default_value = "dfam", value_parser = ["dfam", "dfam-nocpg", "giri", "giri-cpg"])]
    caller: String,
    /// Closest policy: build these families' seeds from exactly the copies
    /// listed, with no coverage fill and no fallback.  Tab-separated lines of
    /// family and copy, a copy written seq:start-end as in the .align record
    /// header (1-based).  For testing a chosen set with the real alignment and
    /// caller.
    #[arg(long)]
    copy_list: Option<PathBuf>,
    /// Closest policy: also search the closest copies for the set whose
    /// column calls best reproduce the library (see search.rs), and keep it
    /// when its real consensus is within --search-margin-per-kb of the best
    /// candidate's edit count without more uncovered bases
    #[arg(long)]
    search: bool,
    /// Search: how many of the closest ranking may be chosen
    #[arg(long, default_value_t = 30)]
    search_pool: usize,
    /// Search: once the most library positions are won, copies are added
    /// while the set loses no more than this many positions per kb beyond
    /// that; deeper seeds are preferred when sets come this close
    #[arg(long, default_value_t = 1.0)]
    search_margin_per_kb: f64,
    /// Search: most copies in a set
    #[arg(long, default_value_t = 50)]
    search_max_copies: usize,
    /// Search: random subsets of the pool to climb from as well; any value
    /// above 0 also adds the top 5, 10, 20 and 30 of the ranking and the
    /// RepBase-sites sample as starts
    #[arg(long, default_value_t = 0)]
    search_restarts: usize,
    /// Search: when single moves stop helping, try two-copy moves
    #[arg(long)]
    search_pairs: bool,
    /// Search: the margin stage may add copies from this many of the closest
    /// ranking, beyond the search pool; those copies take no part in the
    /// search's moves
    #[arg(long, default_value_t = 60)]
    search_depth_pool: usize,
    /// How insertions are merged: per-slot matches generateSeedAlignments,
    /// which opens the widest insertion seen at each slot
    #[arg(long, default_value = "per-slot")]
    insertions: String,
    /// Skip the subfamily test, which splits a family's copies on how much of
    /// the library's version they carry and reports whether RepBase described
    /// one group or their average.  It costs two extra consensus calls per
    /// family.
    #[arg(long)]
    no_subfamily: bool,

    /// Also write the subfamily test as a table here
    #[arg(long)]
    subfamily_out: Option<PathBuf>,

    /// Share of a family's disputed sites a copy must reach to be scored
    #[arg(long, default_value_t = 0.25)]
    subfamily_min_frac: f64,

    /// Copies each subfamily group may carry into its consensus
    #[arg(long, default_value_t = 200)]
    subfamily_cap: usize,

    /// Report the index and stop
    #[arg(long)]
    index_only: bool,
}

struct Built {
    family: String,
    policy: &'static str,
    record: StkRecord,
    cmp: compare::ConsCmp,
    copies_used: usize,
    copies_available: usize,
    dropped: family::Scrubbed,
    coverage: f64,
    /// Sample size and the identity it reached, one entry per scanned size.
    scan: Vec<(usize, f64)>,
    /// Subfamily structure, population policy only.
    sub: Option<SubReport>,
}

struct SubReport {
    scored: usize,
    n_a: usize,
    n_b: usize,
    between: usize,
    id_a: f64,
    id_b: f64,
    cov_a: f64,
    cov_b: f64,
    diagnostic: usize,
    sites: usize,
    span: Option<(usize, usize, f64)>,
    call: &'static str,
}

/// Close an insertion or deletion run that reaches the end of a ConsCAF line.
///
/// dfam-lib's `to_caf` writes the closing '+' or '-' when the next column is
/// aligned or the other kind of run starts, and writes nothing when the
/// alignment ends with a run still open.  A consensus that stops short of the
/// library's end therefore trails one long unterminated deletion, which breaks
/// any reader that pairs delimiters left to right, and breaks it most often on
/// the samples that cover least.  See dfam-lib TODO.md; this closes the run on
/// the way out and does nothing once the library stops leaving it open.
fn close_caf(caf: &str) -> String {
    let Some((head, payload)) = caf.rsplit_once(',') else {
        return caf.to_string();
    };
    let (mut ins, mut del) = (false, false);
    for c in payload.chars() {
        match c {
            '+' if !del => ins = !ins,
            '-' if !ins => del = !del,
            _ => {}
        }
    }
    let mut out = format!("{head},{payload}");
    if del {
        out.push('-');
    } else if ins {
        out.push('+');
    }
    out
}

/// The scanned sample size that came nearest the library consensus.
fn best_scan(scan: &[(usize, f64)]) -> Option<(usize, f64)> {
    scan.iter().copied().max_by(|a, b| a.1.total_cmp(&b.1))
}

/// How the scan reads: `flat` when the best size lands within a point of the
/// sample we kept, `narrow` when five or fewer copies beat it by five points
/// or more, `mixed` in between.  Empty without a scan.
fn source_hint(scan: &[(usize, f64)], kept_identity: f64) -> &'static str {
    let Some((best_k, best_id)) = best_scan(scan) else {
        return "";
    };
    let drop = best_id - kept_identity;
    if drop < 1.0 {
        "flat"
    } else if best_k <= 5 && drop >= 5.0 {
        "narrow"
    } else {
        "mixed"
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }

    let t = std::time::Instant::now();
    let idx = align_index::build(&args.align)?;
    eprintln!(
        "indexed {} records over {} families ({} simple/low/short skipped) in {:.1}s",
        idx.records,
        idx.families.len(),
        idx.skipped_simple,
        t.elapsed().as_secs_f64()
    );
    if args.index_only {
        return Ok(());
    }

    let copy_list: std::collections::HashMap<String, std::collections::HashSet<String>> = match &args.copy_list {
        Some(p) => {
            let mut m: std::collections::HashMap<String, std::collections::HashSet<String>> = Default::default();
            for line in std::fs::read_to_string(p).with_context(|| format!("cannot read {}", p.display()))?.lines() {
                if let Some((f, c)) = line.split_once('\t') {
                    m.entry(f.to_string()).or_default().insert(c.trim().to_string());
                }
            }
            m
        }
        None => Default::default(),
    };
    let assembly = match &args.assembly {
        Some(p) => Some(
            aln_core::twobit::TwoBitReader::open(p)
                .with_context(|| format!("cannot open {}", p.display()))?,
        ),
        None => None,
    };

    // RepeatMasker drops X from a -lib consensus, so a library that still
    // contains X leaves its alignments in an X-free frame.  Which frame a run
    // used is recorded in the .align headers, so reconcile rather than assume:
    // take the library as written when its length matches, and strip X only
    // when that is what agrees.
    let mut stripped = 0usize;
    let mut mismatched: Vec<String> = Vec::new();
    let library: HashMap<String, Vec<u8>> = aln_core::io::read_fasta_file(&args.library)
        .with_context(|| format!("cannot read {}", args.library.display()))?
        .into_iter()
        .map(|s| {
            let name = s.name.split('#').next().unwrap_or(&s.name).to_string();
            let as_written = s.seq.clone();
            let no_x: Vec<u8> = s.seq.iter().filter(|b| !b.eq_ignore_ascii_case(&b'X')).copied().collect();
            // Records can disagree on the consensus length, so count votes
            // rather than trusting the first one.
            let records = idx.families.get(&name).map(Vec::as_slice).unwrap_or(&[]);
            let (mut as_is, mut without_x) = (0usize, 0usize);
            for r in records {
                let n = r.cons_len as usize;
                as_is += usize::from(n == as_written.len());
                without_x += usize::from(n == no_x.len());
            }
            let seq = if as_is >= without_x && as_is > 0 {
                as_written
            } else if without_x > 0 {
                stripped += 1;
                no_x
            } else {
                // Nothing to reconcile when the family has no alignments.
                if !records.is_empty() {
                    mismatched.push(name.clone());
                }
                as_written
            };
            (name, seq)
        })
        .collect();
    if stripped > 0 {
        eprintln!("{stripped} consensi carry X; matched the alignments in their X-free frame");
    }
    if !mismatched.is_empty() {
        eprintln!(
            "{} consensi match no alignment's consensus length: {}",
            mismatched.len(),
            mismatched.join(", ")
        );
    }

    let wanted: Vec<String> = if args.family.is_empty() {
        let mut v: Vec<String> = idx
            .families
            .keys()
            .filter(|f| library.contains_key(*f))
            .cloned()
            .collect();
        v.sort();
        v
    } else {
        args.family.clone()
    };
    eprintln!("{} families to build", wanted.len());

    let insertions = match args.insertions.as_str() {
        "per-slot" => InsertionPolicy::GrowPerSlot,
        "incremental" => InsertionPolicy::GrowIncremental,
        "drop" => InsertionPolicy::Drop,
        other => anyhow::bail!("unknown --insertions {other}"),
    };

    let built: Vec<Built> = wanted
        .par_iter()
        .flat_map(|fam| {
            let mut out = Vec::new();
            let (Some(recs), Some(consensus)) = (idx.families.get(fam), library.get(fam)) else {
                return out;
            };
            let aligner = match compare::global_aligner(lineup_matrix(), 25, 5) {
                Ok(a) => a,
                Err(e) => {
                    eprintln!("{fam}: {e:#}");
                    return out;
                }
            };
            let scrub = ScrubParams {
                min_aligned_length: args.min_aligned_length,
                max_masked_run: args.max_masked_run,
                assembly: assembly.as_ref(),
            };
            let (inst, dropped) = match family::load_instances(
                &args.align,
                fam,
                consensus,
                recs,
                args.assembly_id.as_deref(),
                &scrub,
            ) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("{fam}: {e:#}");
                    return out;
                }
            };
            if let Some(why) = &dropped.first_parse_error {
                eprintln!("{fam}: skipped {} unparsable record(s); first: {why}", dropped.unparsable);
            }
            if inst.is_empty() {
                return out;
            }
            // The population sample's disputed sites, kept for the closest
            // policy's fallback, which runs second.
            let mut pop_disputed: Vec<(usize, u8, u8)> = Vec::new();
            for (policy, label) in [
                (Policy::Population, "population"),
                (Policy::Closest, "closest"),
            ] {
                // The closest policy can search for the copy count that lands
                // nearest the library consensus, which is what reconstructing
                // RepBase's own sample calls for.
                // Source-sample scan: how nearly each sample size reproduces
                // the library consensus.  A family scoring best at one or two
                // copies and falling away from there was most likely built
                // from that handful, which the curator should see even when
                // the sample we keep is the larger one.
                //
                // Each point is the first n of the ranked list, so the samples
                // nest.  The cover-gap top-up below applies to the kept sample
                // only, which leaves `drop` comparing unlike consensi for the
                // few families that get topped up.  The README says so under
                // "What `drop` does not account for".
                let mut scan: Vec<(usize, f64)> = Vec::new();
                if policy == Policy::Closest && !args.no_source_scan {
                    let mut last_len = 0usize;
                    for &n in &args.closest_grid {
                        let params = SampleParams {
                            target_copies: n,
                            min_depth: args.min_depth,
                            min_span_fraction: args.min_span_fraction,
                        };
                        let picked =
                            family::sample(&inst, consensus.len() as u64, policy, &params);
                        if picked.is_empty() || picked.len() == last_len {
                            continue;
                        }
                        last_len = picked.len();
                        if let Some(score) = trial_identity(
                            fam,
                            consensus,
                            &inst,
                            &picked,
                            &aligner,
                            args.min_non_gap,
                            insertions,
                            &args.caller,
                        ) {
                            scan.push((picked.len(), score));
                        }
                    }
                }
                let target = match policy {
                    Policy::Population => args.target_copies,
                    Policy::Closest if args.closest_copies > 0 => args.closest_copies,
                    // Without a fixed target the scan picks the sample size,
                    // which is the old search.
                    Policy::Closest => best_scan(&scan)
                        .map(|(k, _)| k)
                        .unwrap_or(*args.closest_grid.last().unwrap_or(&1)),
                };
                let params = SampleParams {
                    target_copies: target,
                    min_depth: args.min_depth,
                    min_span_fraction: args.min_span_fraction,
                };
                let chosen = family::sample(&inst, consensus.len() as u64, policy, &params);
                if chosen.is_empty() {
                    continue;
                }
                let chosen = if policy == Policy::Closest && args.cover_gap > 0 {
                    let ranked = family::rank_closest(
                        &inst,
                        consensus.len() as u64,
                        args.min_span_fraction,
                    );
                    let filled = family::top_up_coverage(
                        &inst,
                        &ranked,
                        &chosen,
                        consensus.len() as u64,
                        args.cover_gap,
                        args.cover_ends_only,
                    );
                    depth_fill(&inst, &ranked, filled, consensus.len() as u64, &args)
                } else if policy == Policy::Closest {
                    let ranked = family::rank_closest(
                        &inst,
                        consensus.len() as u64,
                        args.min_span_fraction,
                    );
                    depth_fill(&inst, &ranked, chosen, consensus.len() as u64, &args)
                } else {
                    chosen
                };
                // A listed family gets exactly its listed copies.
                let given = policy == Policy::Closest && copy_list.contains_key(fam.as_str());
                let chosen = if given {
                    let want = &copy_list[fam.as_str()];
                    let idx: Vec<usize> = (0..inst.len())
                        .filter(|&i| {
                            let h = &inst[i].hit;
                            want.contains(&format!("{}:{}-{}", h.query_name, h.query.start() + 1, h.query.end()))
                        })
                        .collect();
                    if idx.is_empty() {
                        eprintln!("{fam}: none of the {} listed copies found", want.len());
                        continue;
                    }
                    idx
                } else {
                    chosen
                };
                // Alignment, consensus and comparison to the library for one
                // sample of copies.
                let evaluate = |chosen: &[usize]| -> Result<(
                    aln_core::msa::MultiAlign,
                    Vec<u8>,
                    Vec<bool>,
                    f64,
                    compare::ConsCmp,
                )> {
                    let mut msa = family::build_msa(fam, consensus, &inst, chosen, insertions)?;
                    // generateSeedAlignments' alignment spans only its
                    // instances, so reference positions outside them never
                    // become columns.  Trim to the covered extent before
                    // calling, or the flanks would collect N and the consensus
                    // would keep the reference's full length.
                    let covered = covered_positions(consensus.len(), &inst, chosen);
                    if let Some((lo, hi)) = covered_extent_columns(&msa, &covered) {
                        msa.slice_columns(lo, hi + 1);
                    }
                    let gapped = call_consensus(&msa, args.min_non_gap, &args.caller);
                    let called: Vec<u8> =
                        gapped.iter().filter(|b| b.is_ascii_alphabetic()).copied().collect();
                    let coverage =
                        covered.iter().filter(|c| **c).count() as f64 / consensus.len() as f64;
                    let cmp = compare::compare(&aligner, consensus, &called, &covered)?;
                    Ok((msa, gapped, covered, coverage, cmp))
                };
                let (mut msa, mut gapped, covered, mut coverage, mut cmp) = match evaluate(&chosen) {
                    Ok(v) => v,
                    Err(e) => {
                        eprintln!("{fam} ({label}): {e:#}");
                        continue;
                    }
                };
                let mut chosen = chosen;
                // the plain closest sample, one of the search's starts
                let plain = chosen.clone();
                let mut selection: Option<String> = None;
                if given {
                    selection = Some(format!("Selection: given copy list, {} copies", chosen.len()));
                }
                if policy == Policy::Population {
                    pop_disputed = cmp.disputed.clone();
                }
                // The plain closest sample's uncovered bases, the bound every
                // fallback is held to under --best-of.
                let unc0 = covered.iter().filter(|c| !**c).count();
                if policy == Policy::Closest && !given && !args.no_disputed_fallback && !pop_disputed.is_empty() {
                    let lib_len = consensus.len();
                    let before = unsupported(&inst, &chosen, &cmp.disputed, lib_len);
                    let e1 = edits(&cmp);
                    if before > 0 || (args.best_of && e1 > 0) {
                        let (cand, n_like) = disputed_fallback(
                            &inst,
                            &pop_disputed,
                            lib_len,
                            target,
                            &args,
                        );
                        match evaluate(&cand) {
                            Ok((m2, g2, cv2, co2, c2)) => {
                                let after = unsupported(&inst, &cand, &c2.disputed, lib_len);
                                let unc1 = covered.iter().filter(|c| !**c).count();
                                let unc2 = cv2.iter().filter(|c| !**c).count();
                                let e2 = edits(&c2);
                                let take = if args.best_of {
                                    e2 < e1 && unc2 <= unc0
                                } else {
                                    after < before && unc2 <= unc1
                                };
                                if take {
                                    selection = Some(format!(
                                        "Selection: closest-by-disputed-sites, unsupported {before} -> {after}, \
                                         differences {e1} -> {e2}, uncovered {unc1} -> {unc2}, library-like copies {n_like}"
                                    ));
                                    msa = m2;
                                    gapped = g2;
                                    coverage = co2;
                                    cmp = c2;
                                    chosen = cand;
                                } else {
                                    selection = Some(format!(
                                        "Selection: closest; closest-by-disputed-sites tried, unsupported \
                                         {before} -> {after}, differences {e1} -> {e2}, uncovered {unc1} -> {unc2}, \
                                         library-like copies {n_like}"
                                    ));
                                }
                            }
                            Err(e) => eprintln!("{fam} (closest-by-disputed-sites): {e:#}"),
                        }
                    }
                }
                // Further samples, each tried against the seed kept so far and
                // kept when its consensus leaves fewer differences from the
                // library without more uncovered bases.
                if policy == Policy::Closest && !given && !pop_disputed.is_empty() {
                    let lib_len = consensus.len();
                    let mut tries: Vec<(&str, Vec<usize>, String)> = Vec::new();
                    if (args.lineage_fallback || args.lineage_always)
                        && (args.lineage_always || edits(&cmp) > 0)
                    {
                        if let Some((c, st)) = lineage_fallback(&inst, &pop_disputed, lib_len, target, &args) {
                            let label = if args.lineage_always { "lineage" } else { "closest-by-lineage" };
                            tries.push((label, c, format!(
                                "informative sites {}, dispersion {:.2}, lineage copies {}",
                                st.sites, st.dispersion, st.lineage_copies
                            )));
                        }
                    }
                    if (args.rbsites_fallback || args.rbsites_always)
                        && (args.rbsites_always || edits(&cmp) > 0)
                    {
                        if let Some((c, st)) = rbsites_fallback(&inst, &pop_disputed, lib_len, target, &args) {
                            tries.push(("closest-by-repbase-sites", c, format!(
                                "RepBase sites {}, dispersion {:.2}, RepBase-like copies {}, {}",
                                st.sites, st.dispersion, st.lineage_copies,
                                if args.rbsites_weighted { "rarity-weighted" } else { "unweighted" }
                            )));
                        }
                    }
                    let mut kept = match &selection {
                        Some(s) if s.starts_with("Selection: closest-by-disputed-sites") => "closest-by-disputed-sites".to_string(),
                        _ => "closest".to_string(),
                    };
                    let mut notes: Vec<String> = selection
                        .take()
                        .map(|s| {
                            // "closest; x tried" or "closest-by-disputed-sites, ..." -> the note alone
                            let s = s.trim_start_matches("Selection: ");
                            vec![s.strip_prefix("closest; ").unwrap_or(s).to_string()]
                        })
                        .unwrap_or_default();
                    for (label, cand, stats) in tries {
                        match evaluate(&cand) {
                            Ok((m2, g2, cv2, co2, c2)) => {
                                let e1 = edits(&cmp);
                                let e2 = edits(&c2);
                                let unc1 = if args.best_of {
                                    unc0
                                } else {
                                    covered_positions(lib_len, &inst, &chosen).iter().filter(|c| !**c).count()
                                };
                                let unc2 = cv2.iter().filter(|c| !**c).count();
                                let take = (args.lineage_always && label == "lineage")
                                    || (args.rbsites_always && label == "closest-by-repbase-sites")
                                    || (e2 < e1 && unc2 <= unc1);
                                notes.push(format!(
                                    "{label} {}, differences {e1} -> {e2}, uncovered {unc1} -> {unc2}, {stats}",
                                    if take { "kept" } else { "tried" }
                                ));
                                if take {
                                    kept = label.to_string();
                                    msa = m2;
                                    gapped = g2;
                                    coverage = co2;
                                    cmp = c2;
                                    chosen = cand;
                                }
                            }
                            Err(e) => eprintln!("{fam} ({label}): {e:#}"),
                        }
                    }
                    if !notes.is_empty() {
                        selection = Some(format!("Selection: {kept}; {}", notes.join("; ")));
                    }
                }
                // The search: from the seed kept so far and from the plain
                // closest sample, climb on the caller's column rule over the
                // closest copies; build the result for real and keep it when
                // its edits are within the margin of the kept seed's.
                if policy == Policy::Closest && !given && args.search {
                    let lib_len = consensus.len();
                    let ranked = family::rank_closest(&inst, lib_len as u64, args.min_span_fraction);
                    let rule = if args.caller.starts_with("giri") { search::Rule::giri() } else { search::Rule::dfam() };
                    let sp = search::Params {
                        pool: args.search_pool,
                        max_copies: args.search_max_copies,
                        margin_per_kb: args.search_margin_per_kb,
                        restarts: args.search_restarts,
                        pairs: args.search_pairs,
                        depth_pool: args.search_depth_pool,
                    };
                    let mut starts = vec![chosen.clone(), plain.clone()];
                    if args.search_restarts > 0 && !pop_disputed.is_empty() {
                        if let Some((c, _)) = rbsites_fallback(&inst, &pop_disputed, lib_len, target, &args) {
                            starts.push(c);
                        }
                    }
                    if let Some((set, lost)) = search::search(&inst, &ranked, &starts, consensus, &rule, unc0, &sp) {
                        if set != chosen {
                            match evaluate(&set) {
                                Ok((m2, g2, cv2, co2, c2)) => {
                                    let (e1, e2) = (edits(&cmp), edits(&c2));
                                    let unc2 = cv2.iter().filter(|c| !**c).count();
                                    let slack = (args.search_margin_per_kb * lib_len as f64 / 1000.0) as usize;
                                    let take = e2 <= e1 + slack && unc2 <= unc0;
                                    let note = format!(
                                        "closest-by-search {}, differences {e1} -> {e2}, uncovered {unc0} -> {unc2}, \
                                         copies {} -> {}, column-rule losses {lost}",
                                        if take { "kept" } else { "tried" },
                                        chosen.len(),
                                        set.len()
                                    );
                                    let prior = selection
                                        .take()
                                        .map(|s| s.trim_start_matches("Selection: ").to_string());
                                    let prior_kept = prior
                                        .as_deref()
                                        .map(|p| p.split([';', ',']).next().unwrap_or("closest").trim().to_string())
                                        .unwrap_or_else(|| "closest".to_string());
                                    let rest = prior
                                        .map(|p| {
                                            let p = p.strip_prefix(&format!("{prior_kept}; ")).unwrap_or(&p).to_string();
                                            format!("; {p}")
                                        })
                                        .unwrap_or_default();
                                    if take {
                                        selection = Some(format!("Selection: closest-by-search; {note}{rest}"));
                                        msa = m2;
                                        gapped = g2;
                                        coverage = co2;
                                        cmp = c2;
                                        chosen = set;
                                    } else {
                                        selection = Some(format!("Selection: {prior_kept}; {note}{rest}"));
                                    }
                                }
                                Err(e) => eprintln!("{fam} (closest-by-search): {e:#}"),
                            }
                        }
                    }
                }
                let sub = if policy == Policy::Population && !args.no_subfamily {
                    let pool: Vec<&family::Instance> = inst.iter().collect();
                    subfamily::groups(
                        &pool,
                        &cmp.disputed,
                        consensus.len(),
                        args.subfamily_min_frac,
                        args.subfamily_cap,
                    )
                    .and_then(|g| {
                        // Each side gets its own consensus, and comparing both
                        // to the library says which side the library is.
                        let side = |idx: &[usize]| -> Option<(f64, f64)> {
                            if idx.len() < 5 {
                                return None;
                            }
                            let msa =
                                family::build_msa(fam, consensus, &inst, idx, insertions).ok()?;
                            let gapped = call_consensus(&msa, args.min_non_gap, &args.caller);
                            let called: Vec<u8> = gapped
                                .iter()
                                .filter(|b| b.is_ascii_alphabetic())
                                .copied()
                                .collect();
                            let covered = covered_positions(consensus.len(), &inst, idx);
                            let frac = covered.iter().filter(|c| **c).count() as f64
                                / consensus.len().max(1) as f64;
                            compare::compare(&aligner, consensus, &called, &covered)
                                .ok()
                                .map(|c| (c.identity_occupied(), frac))
                        };
                        let (id_a, cov_a) = side(&g.a).unwrap_or((f64::NAN, f64::NAN));
                        let (id_b, cov_b) = side(&g.b).unwrap_or((f64::NAN, f64::NAN));
                        Some(SubReport {
                            scored: g.scored,
                            n_a: g.a.len(),
                            n_b: g.b.len(),
                            between: g.between,
                            id_a,
                            id_b,
                            cov_a,
                            cov_b,
                            diagnostic: g.diagnostic.len(),
                            sites: cmp.disputed.len(),
                            span: g.diagnostic_span(consensus.len()),
                            call: subfamily::call(&g, id_a, id_b),
                        })
                    })
                } else {
                    None
                };
                cmp.agreed = Vec::new();   // large, and nothing downstream reads it
                let mut record = to_stk(
                    fam,
                    label,
                    &msa,
                    &gapped,
                    &cmp,
                    chosen.len(),
                    args.taxon.as_deref(),
                    args.assembly_id.as_deref(),
                    &scan,
                    sub.as_ref(),
                );
                if let Some(line) = selection {
                    let at = record
                        .gf
                        .iter()
                        .position(|(t, v)| t == "**" && v.starts_with("Source:"))
                        .map(|i| i + 1)
                        .unwrap_or(record.gf.len());
                    record.gf.insert(at, ("**".into(), line));
                }
                out.push(Built {
                    family: fam.clone(),
                    policy: label,
                    record,
                    cmp,
                    copies_used: chosen.len(),
                    copies_available: inst.len(),
                    dropped: dropped.clone(),
                    coverage,
                    scan,
                    sub,
                });
            }
            out
        })
        .collect();

    write_outputs(&args, built)?;
    eprintln!("done in {:.1}s", t.elapsed().as_secs_f64());
    Ok(())
}

/// Disputed sites (positions where a sample's consensus calls another base
/// than the library) at which `MIN_UNSUPPORTED_DEPTH` or more of the sampled
/// copies reach and none carries the library's base.  Indel differences are
/// not counted: `compare` lists substitutions only.
fn unsupported(
    inst: &[family::Instance],
    chosen: &[usize],
    sites: &[(usize, u8, u8)],
    lib_len: usize,
) -> usize {
    const MIN_UNSUPPORTED_DEPTH: usize = 3;
    if sites.is_empty() || chosen.is_empty() {
        return 0;
    }
    let refs: Vec<&family::Instance> = chosen.iter().map(|&i| &inst[i]).collect();
    let rows = subfamily::rows(&refs, sites, lib_len);
    (0..sites.len())
        .filter(|&j| {
            let reach = rows.iter().filter(|r| r[j] != 0).count();
            let agree = rows.iter().filter(|r| r[j] != 0 && r[j] == sites[j].1).count();
            reach >= MIN_UNSUPPORTED_DEPTH && agree == 0
        })
        .count()
}

/// What the lineage fallback saw: informative sites, how far the copies'
/// agreement with the library at them spreads beyond chance, and how many
/// copies carry the library's base at most of them.
struct LineageStats {
    sites: usize,
    dispersion: f64,
    lineage_copies: usize,
}

/// The "closest-by-lineage" sample.  A copy's Kimura divergence is mostly its
/// own decay since it inserted, so ranking by it barely separates the lineage
/// the library describes from a sister lineage that differs at a few
/// positions.  This ranks copies at those positions only.
///
/// Informative sites are the population sample's disputed sites (where the
/// family's majority departs from the library) at which, over every copy,
/// 10 copies or more reach the site and 10-50% carry the library's base.  A
/// copy reaching at least min(5, sites) of them scores the share where it
/// carries the library's base.  If copies carried the library's base at each
/// site independently, scores would spread as a binomial; the dispersion is
/// their observed variance over that, and a lineage shows as well above 1.
/// With dispersion of at least --lineage-min-dispersion and at least
/// --lineage-min-copies copies scoring 0.8 or more, the sample is the
/// highest-scoring copies (the closest ranking breaks ties), then the
/// coverage and depth fills.
///
/// In a prototype on 244 RepBase 31.08 families whose closest consensus was
/// under 99% identical to the library (majority vote over RepeatMasker's
/// pairwise alignments, no fills), the top 10 by this score left 28% fewer
/// differences from the library than the closest 10, and families at 99% or
/// more went from 16 to 34; where copies did not group (dispersion under 2)
/// it did worse, hence the gate.
fn lineage_fallback(
    inst: &[family::Instance],
    pop_disputed: &[(usize, u8, u8)],
    lib_len: usize,
    target: usize,
    args: &Args,
) -> Option<(Vec<usize>, LineageStats)> {
    const MIN_COVER: usize = 10;
    let refs: Vec<&family::Instance> = inst.iter().collect();
    let rows = subfamily::rows(&refs, pop_disputed, lib_len);
    let mut freq: Vec<Option<f64>> = vec![None; pop_disputed.len()];
    for (j, (_, lib, _)) in pop_disputed.iter().enumerate() {
        let reach = rows.iter().filter(|r| r[j] != 0).count();
        let agree = rows.iter().filter(|r| r[j] != 0 && r[j] == *lib).count();
        if reach >= MIN_COVER {
            let f = agree as f64 / reach as f64;
            if (0.10..=0.50).contains(&f) {
                freq[j] = Some(f);
            }
        }
    }
    let sites = freq.iter().filter(|f| f.is_some()).count();
    if sites < 3 {
        return None;
    }
    let need = sites.min(5);
    let mut score = vec![f64::NAN; inst.len()];
    let (mut dev2, mut var) = (Vec::new(), Vec::new());
    for (i, row) in rows.iter().enumerate() {
        let mut n = 0usize;
        let (mut m, mut e, mut v) = (0usize, 0.0, 0.0);
        for (j, f) in freq.iter().enumerate() {
            if let (Some(f), b) = (f, row[j]) {
                if b != 0 {
                    n += 1;
                    e += f;
                    v += f * (1.0 - f);
                    if b == pop_disputed[j].1 {
                        m += 1;
                    }
                }
            }
        }
        if n >= need {
            let s = m as f64 / n as f64;
            score[i] = s;
            dev2.push((s - e / n as f64).powi(2));
            var.push(v / (n * n) as f64);
        }
    }
    if dev2.len() < 5 {
        return None;
    }
    let mean = |x: &[f64]| x.iter().sum::<f64>() / x.len() as f64;
    let dispersion = mean(&dev2) / mean(&var).max(1e-12);
    let lineage_copies = score.iter().filter(|s| **s >= 0.8).count();
    if dispersion < args.lineage_min_dispersion || lineage_copies < args.lineage_min_copies {
        return None;
    }
    let ranked = family::rank_closest(inst, lib_len as u64, args.min_span_fraction);
    let mut order: Vec<usize> = ranked.iter().copied().filter(|&i| !score[i].is_nan()).collect();
    // stable: equal scores keep the closest ranking's order
    order.sort_by(|&a, &b| score[b].partial_cmp(&score[a]).unwrap_or(std::cmp::Ordering::Equal));
    let mut cand: Vec<usize> = order.into_iter().take(target).collect();
    if args.cover_gap > 0 {
        cand = family::top_up_coverage(inst, &ranked, &cand, lib_len as u64, args.cover_gap, args.cover_ends_only);
    }
    let cand = depth_fill(inst, &ranked, cand, lib_len as u64, args);
    Some((cand, LineageStats { sites, dispersion, lineage_copies }))
}

/// The "closest-by-repbase-sites" sample: the copies that carry the
/// library's bases where the family's majority does not.  It generalises the
/// disputed-sites and lineage samples for the question of whether the genome
/// holds copies that support the library, whatever their number.
///
/// Sites are the population sample's disputed sites that --rbsites-min-cover
/// copies (over every copy) reach and at least one of them carries the
/// library's base there.  There is no upper or lower bound on how many carry
/// it: a curator who built from a handful of copies left the library's bases
/// at sites only those copies share, and those rare sites are the strongest
/// evidence of which copies they were.  A copy reaching at least min(3,
/// sites) of them scores the share where it carries the library's base, or
/// with --rbsites-weighted the share weighted by -ln of each site's
/// frequency.  At least 3 sites and --rbsites-min-copies copies scoring 0.8
/// or more are needed.  The dispersion (see lineage_fallback) is reported,
/// not used: whether the copies form a lineage is a question about the
/// family, not about whether they support the library.
///
/// The sample is the highest-scoring copies (ties to the earlier copy in the
/// closest ranking), or with --rbsites-rank-order the copies scoring 0.8 or
/// more in the closest ranking's order and then the rest by score; then the
/// coverage and depth fills.
fn rbsites_fallback(
    inst: &[family::Instance],
    pop_disputed: &[(usize, u8, u8)],
    lib_len: usize,
    target: usize,
    args: &Args,
) -> Option<(Vec<usize>, LineageStats)> {
    let refs: Vec<&family::Instance> = inst.iter().collect();
    let rows = subfamily::rows(&refs, pop_disputed, lib_len);
    let mut freq: Vec<Option<f64>> = vec![None; pop_disputed.len()];
    for (j, (_, lib, _)) in pop_disputed.iter().enumerate() {
        let reach = rows.iter().filter(|r| r[j] != 0).count();
        let agree = rows.iter().filter(|r| r[j] != 0 && r[j] == *lib).count();
        if reach >= args.rbsites_min_cover && agree > 0 {
            freq[j] = Some(agree as f64 / reach as f64);
        }
    }
    let sites = freq.iter().filter(|f| f.is_some()).count();
    if sites < 3 {
        return None;
    }
    let need = sites.min(3);
    let mut score = vec![f64::NAN; inst.len()];
    let (mut dev2, mut var) = (Vec::new(), Vec::new());
    for (i, row) in rows.iter().enumerate() {
        let (mut n, mut m) = (0usize, 0usize);
        let (mut w, mut wm, mut e, mut v) = (0.0, 0.0, 0.0, 0.0);
        for (j, f) in freq.iter().enumerate() {
            if let (Some(f), b) = (f, row[j]) {
                if b != 0 {
                    let wt = if args.rbsites_weighted { -f.ln() } else { 1.0 };
                    n += 1;
                    w += wt;
                    e += f;
                    v += f * (1.0 - f);
                    if b == pop_disputed[j].1 {
                        m += 1;
                        wm += wt;
                    }
                }
            }
        }
        if n >= need {
            score[i] = if w > 0.0 { wm / w } else { m as f64 / n as f64 };
            dev2.push((m as f64 / n as f64 - e / n as f64).powi(2));
            var.push(v / (n * n) as f64);
        }
    }
    let liked = score.iter().filter(|s| **s >= 0.8).count();
    if liked < args.rbsites_min_copies {
        return None;
    }
    let mean = |x: &[f64]| if x.is_empty() { 0.0 } else { x.iter().sum::<f64>() / x.len() as f64 };
    let dispersion = mean(&dev2) / mean(&var).max(1e-12);
    let ranked = family::rank_closest(inst, lib_len as u64, args.min_span_fraction);
    let mut order: Vec<usize> = ranked.iter().copied().filter(|&i| !score[i].is_nan()).collect();
    if args.rbsites_rank_order {
        // RepBase-like copies in ranking order, then the rest by score (or in
        // ranking order with --rbsites-pad-rank)
        let (mut like, mut rest): (Vec<usize>, Vec<usize>) = order.into_iter().partition(|&i| score[i] >= 0.8);
        if args.rbsites_strict_many && like.len() >= target {
            let need = (args.min_span_fraction * lib_len as f64) as u64;
            let (mut long, short): (Vec<usize>, Vec<usize>) = like.into_iter().partition(|&i| inst[i].span() >= need);
            long.sort_by(|&a, &b| score[b].partial_cmp(&score[a]).unwrap_or(std::cmp::Ordering::Equal));
            like = long.into_iter().chain(short).collect();
        }
        if !args.rbsites_pad_rank {
            rest.sort_by(|&a, &b| score[b].partial_cmp(&score[a]).unwrap_or(std::cmp::Ordering::Equal));
        }
        order = like.into_iter().chain(rest).collect();
    } else {
        order.sort_by(|&a, &b| score[b].partial_cmp(&score[a]).unwrap_or(std::cmp::Ordering::Equal));
    }
    let mut cand: Vec<usize> = order.into_iter().take(target).collect();
    if args.cover_gap > 0 {
        cand = family::top_up_coverage(inst, &ranked, &cand, lib_len as u64, args.cover_gap, args.cover_ends_only);
    }
    let cand = depth_fill(inst, &ranked, cand, lib_len as u64, args);
    Some((cand, LineageStats { sites, dispersion, lineage_copies: liked }))
}

/// The "closest-by-disputed-sites" sample.  The population sample's disputed
/// sites are where the family as a whole departs from the library, so they
/// are where the library's subfamily shows.  Copies carrying the library's
/// base at 80% or more of the disputed sites they reach (and reaching at
/// least three, or all of them when there are fewer) are library-like; they
/// come first, in the closest ranking, then the rest, and the cover-gap
/// top-up fills what they leave uncovered.  Returns the sample and how many
/// library-like copies there were.
///
/// In a test on 522 RepBase 31.08 families whose closest sample left 1 to 10
/// sites unsupported, using this sample only when it left fewer unsupported
/// sites and no more uncovered bases raised the families with none from 65
/// to 104; it switched 1 of 149 families the closest sample already supported.
fn disputed_fallback(
    inst: &[family::Instance],
    pop_disputed: &[(usize, u8, u8)],
    lib_len: usize,
    target: usize,
    args: &Args,
) -> (Vec<usize>, usize) {
    let refs: Vec<&family::Instance> = inst.iter().collect();
    let rows = subfamily::rows(&refs, pop_disputed, lib_len);
    let need = pop_disputed.len().min(3);
    let like: Vec<bool> = rows
        .iter()
        .map(|row| {
            let reach = row.iter().filter(|b| **b != 0).count();
            let agree = row
                .iter()
                .zip(pop_disputed)
                .filter(|(b, (_, lib, _))| **b != 0 && **b == *lib)
                .count();
            reach >= need && agree * 5 >= reach * 4
        })
        .collect();
    let n_like = like.iter().filter(|l| **l).count();
    let ranked = family::rank_closest(inst, lib_len as u64, args.min_span_fraction);
    let mut cand: Vec<usize> = ranked.iter().copied().filter(|&i| like[i]).take(target).collect();
    for &i in &ranked {
        if cand.len() >= target {
            break;
        }
        if !like[i] {
            cand.push(i);
        }
    }
    if args.cover_gap > 0 {
        cand = family::top_up_coverage(
            inst,
            &ranked,
            &cand,
            lib_len as u64,
            args.cover_gap,
            args.cover_ends_only,
        );
    }
    let cand = depth_fill(inst, &ranked, cand, lib_len as u64, args);
    (cand, n_like)
}

/// The closest policy's depth fill, when --closest-min-depth asks for one.
fn depth_fill(
    inst: &[family::Instance],
    ranked: &[usize],
    chosen: Vec<usize>,
    lib_len: u64,
    args: &Args,
) -> Vec<usize> {
    if args.closest_min_depth == 0 {
        return chosen;
    }
    family::top_up_depth(inst, ranked, &chosen, lib_len, args.closest_min_depth, args.fill_kimura_window)
}

/// Build a candidate sample's consensus and report its identity to the
/// library, or None if it covers too little to judge.
#[allow(clippy::too_many_arguments)]
fn trial_identity(
    family: &str,
    consensus: &[u8],
    inst: &[family::Instance],
    chosen: &[usize],
    aligner: &aln_parasail::ParasailAligner,
    min_non_gap: usize,
    insertions: InsertionPolicy,
    caller: &str,
) -> Option<f64> {
    let msa = family::build_msa(family, consensus, inst, chosen, insertions).ok()?;
    let gapped = call_consensus(&msa, min_non_gap, caller);
    let called: Vec<u8> = gapped.iter().filter(|b| b.is_ascii_alphabetic()).copied().collect();
    let covered = covered_positions(consensus.len(), inst, chosen);
    if covered.iter().filter(|c| **c).count() * 5 < consensus.len() * 4 {
        return None;                      // under 80% of the consensus reached
    }
    let cmp = compare::compare(aligner, consensus, &called, &covered).ok()?;
    let id = cmp.identity();
    id.is_finite().then_some(id)
}

/// First and last column, inclusive, that the instances cover.  Columns are
/// found through the reference row, so insertions inside the span come along.
fn covered_extent_columns(
    msa: &aln_core::msa::MultiAlign,
    covered: &[bool],
) -> Option<(usize, usize)> {
    let reference = msa.reference()?;
    let mut ref_pos = 0usize;
    let (mut lo, mut hi) = (None, 0usize);
    for (col, &b) in reference.seq.iter().enumerate() {
        if aln_core::seq::is_gap(b) || aln_core::seq::is_pad(b) {
            continue;
        }
        if covered.get(ref_pos).copied().unwrap_or(false) {
            lo.get_or_insert(col);
            hi = col;
        }
        ref_pos += 1;
    }
    lo.map(|l| (l, hi))
}

/// Which library positions any chosen copy reached.
fn covered_positions(len: usize, inst: &[family::Instance], chosen: &[usize]) -> Vec<bool> {
    let mut covered = vec![false; len];
    for &i in chosen {
        for p in inst[i].rec.cons_start..=inst[i].rec.cons_end {
            if let Some(slot) = covered.get_mut((p - 1) as usize) {
                *slot = true;
            }
        }
    }
    covered
}

/// Every difference between a called consensus and the library that edit
/// identity counts: substitutions, ambiguity codes, insertions, deletions and
/// library bases no copy covers.  Edit identity is 1 - edits / library
/// length; the fallbacks compete on it.
fn edits(c: &compare::ConsCmp) -> usize {
    c.sub + c.amb + c.ins + c.del + c.unoccupied
}

/// The gapped consensus of `msa` under --caller.
fn call_consensus(msa: &aln_core::msa::MultiAlign, min_non_gap: usize, caller: &str) -> Vec<u8> {
    let mut params = ConsParams { min_non_gap_count: min_non_gap, ..ConsParams::default() };
    match caller {
        "giri" | "giri-cpg" => params.caller = cons_core::Caller::Giri,
        "dfam-nocpg" => params.consensus.enable_cpg = false,
        _ => {}
    }
    let mut gapped = cons_core::gapped_consensus(msa, &params);
    if caller == "giri-cpg" {
        // as cons-core's refinement does: the reference row is left out
        let rows: Vec<&[u8]> = msa.sequences[1.min(msa.sequences.len())..]
            .iter()
            .map(|r| r.seq.as_slice())
            .collect();
        aln_core::giri::restore_cpg(&mut gapped, &rows);
    }
    gapped
}

/// The lineup matrix, taken from the consensus caller itself so the
/// consensus-to-consensus comparison scores exactly what the caller scored.
/// The gap symbol is dropped: parasail takes gap costs as parameters.
fn lineup_matrix() -> SubstMatrix {
    const ALPHABET: &[u8] = b"ARGCYTKMSWNXZVHDB";
    let n = ALPHABET.len();
    let mut scores = Vec::with_capacity(n * n);
    for r in 0..n {
        for c in 0..n {
            scores.push(aln_core::consensus::MATRIX[r][c]);
        }
    }
    SubstMatrix::from_parts("linup", ALPHABET, scores, None).expect("lineup matrix is well formed")
}

#[allow(clippy::too_many_arguments)]
fn to_stk(
    family: &str,
    policy: &str,
    msa: &aln_core::msa::MultiAlign,
    gapped_consensus: &[u8],
    cmp: &compare::ConsCmp,
    copies_used: usize,
    taxon: Option<&str>,
    assembly_id: Option<&str>,
    scan: &[(usize, f64)],
    sub: Option<&SubReport>,
) -> StkRecord {
    let mut gf: Vec<(String, String)> = vec![("ID".into(), family.into())];
    if let Some(t) = taxon {
        gf.push(("OC".into(), t.into()));
    }
    gf.push((
        "CC".into(),
        format!(
            "Seed alignment built by recover-seeds from RepeatMasker alignments, \
             {policy} sampling, {copies_used} copies."
        ),
    ));
    gf.push((
        "**".into(),
        format!(
            "Source:recover-seeds, policy={policy}, copies={copies_used}{}",
            assembly_id
                .map(|a| format!(", assembly={a}"))
                .unwrap_or_default()
        ),
    ));
    gf.push((
        "**".into(),
        format!(
            "ConsCmp: src=library, sub={}, ins={}, del={} ({} unoccupied), amb={}, unocc={}, \
             ext5={}, ext3={}, trim5={} ({} unoccupied), trim3={} ({} unoccupied), \
             lens={}/{}, change={:+}, id={:.2}%, id_occupied={:.2}%",
            cmp.sub,
            cmp.ins,
            cmp.del,
            cmp.del_uncovered,
            cmp.amb,
            cmp.unoccupied,
            cmp.ext5,
            cmp.ext3,
            cmp.trim5,
            cmp.trim5_uncovered,
            cmp.trim3,
            cmp.trim3_uncovered,
            cmp.called_len,
            cmp.library_len,
            cmp.change(),
            cmp.identity(),
            cmp.identity_occupied()
        ),
    ));
    if !scan.is_empty() {
        let (best_k, best_id) = best_scan(scan).unwrap();
        let profile: Vec<String> = scan.iter().map(|(k, id)| format!("{k}:{id:.2}")).collect();
        gf.push((
            "**".into(),
            format!(
                "SourceScan: {} | best={} copies at {:.2}%, kept {} copies at {:.2}%, \
                 drop={:.2}, hint={}",
                profile.join(" "),
                best_k,
                best_id,
                copies_used,
                cmp.identity(),
                best_id - cmp.identity(),
                source_hint(scan, cmp.identity()),
            ),
        ));
    }
    let ranges: Vec<String> = cmp
        .unoccupied_ranges
        .iter()
        .take(8)
        .map(|(a, b)| if a == b { a.to_string() } else { format!("{a}-{b}") })
        .collect();
    gf.push((
        "**".into(),
        format!(
            "Occupancy: reference={} bp, covered={:.1}%, unoccupied={} bp in {} region(s), \
             longest={} bp{}{}",
            cmp.library_len,
            100.0 * (cmp.library_len - cmp.unoccupied_ranges.iter().map(|(a, b)| b - a + 1).sum::<usize>()) as f64
                / cmp.library_len.max(1) as f64,
            cmp.unoccupied_ranges.iter().map(|(a, b)| b - a + 1).sum::<usize>(),
            cmp.unoccupied_ranges.len(),
            cmp.longest_unoccupied,
            if ranges.is_empty() { String::new() } else { format!(", at {}", ranges.join(" ")) },
            if cmp.unoccupied_ranges.len() > 8 { ", ..." } else { "" }
        ),
    ));
    if let Some(x) = sub {
        let span = match x.span {
            Some((lo, hi, frac)) => {
                format!(", at {lo}-{hi} ({:.0}% of the consensus)", frac * 100.0)
            }
            None => String::new(),
        };
        let pct = |v: f64| if v.is_nan() { "na".to_string() } else { format!("{v:.2}%") };
        let cov = |v: f64| {
            if v.is_nan() {
                "na".to_string()
            } else {
                format!("{:.0}%", v * 100.0)
            }
        };
        gf.push((
            "**".into(),
            format!(
                "Subfamily: call={}, scored={}, library-like={}, other={}, between={}, \
                 idLibraryLike={} (covers {}), idOther={} (covers {}), diagnostic={}/{}{}",
                x.call,
                x.scored,
                x.n_a,
                x.n_b,
                x.between,
                pct(x.id_a),
                cov(x.cov_a),
                pct(x.id_b),
                cov(x.cov_b),
                x.diagnostic,
                x.sites,
                span,
            ),
        ));
    }
    gf.push(("**".into(), format!("ConsCAF: {}", close_caf(&cmp.caf))));

    // Dfam writes every gap and every flanking pad as '.', so an MSA's '-'
    // and ' ' both become '.' on the way out.  A space inside the alignment
    // field would not survive a parser in any case, and RF takes the same rule.
    let dotted = |seq: &[u8]| -> String {
        seq.iter()
            .map(|&b| if b.is_ascii_alphabetic() { b as char } else { '.' })
            .collect()
    };
    let sequences: Vec<SeqRow> = msa.sequences[1.min(msa.sequences.len())..]
        .iter()
        .map(|r| SeqRow::new_raw(r.name.clone(), dotted(&r.seq)))
        .collect();
    let mut gc = HashMap::new();
    gc.insert(
        "RF".to_string(),
        dotted(gapped_consensus),
    );
    StkRecord {
        record_num: 0,
        start_line: 0,
        header: "# STOCKHOLM 1.0".into(),
        gf,
        gc,
        sequences,
        unknown_annotations: Vec::new(),
        terminated: true,
        block_separator_line: None,
    }
}

fn write_outputs(args: &Args, built: Vec<Built>) -> Result<()> {
    let prefix = args.out_prefix.display().to_string();
    let mut files: HashMap<&str, BufWriter<File>> = HashMap::new();
    for policy in ["population", "closest"] {
        let path = format!("{prefix}-{policy}.stk");
        files.insert(
            policy,
            BufWriter::new(File::create(&path).with_context(|| path.clone())?),
        );
    }
    let mut by_family: HashMap<String, HashMap<&str, &Built>> = HashMap::new();
    for b in &built {
        by_family
            .entry(b.family.clone())
            .or_default()
            .insert(b.policy, b);
        if let Some(f) = files.get_mut(b.policy) {
            b.record.write_to(f)?;
        }
    }
    for (_, mut f) in files {
        f.flush()?;
    }

    let path = format!("{prefix}-summary.tsv");
    let mut tsv = BufWriter::new(File::create(&path).with_context(|| path.clone())?);
    writeln!(
        tsv,
        "family\tlibrary_len\tcopies_available\t\
         pop_copies\tpop_identity\tpop_change\tpop_ext\tpop_trim\tpop_uncovered\tpop_coverage\t\
         near_copies\tnear_identity\tnear_change\tnear_ext\tnear_trim\tnear_uncovered\tnear_coverage\t\
         scan_best_k\tscan_best_identity\tscan_drop\tsource_hint"
    )?;
    let mut families: Vec<&String> = by_family.keys().collect();
    families.sort();
    for fam in families {
        let m = &by_family[fam];
        let cell = |b: Option<&&Built>| match b {
            Some(b) => format!(
                "{}\t{:.4}\t{}\t{}\t{}\t{}\t{:.3}",
                b.copies_used,
                b.cmp.identity(),
                b.cmp.change(),
                b.cmp.ext5 + b.cmp.ext3,
                b.cmp.trim5 + b.cmp.trim3,
                b.cmp.trim5_uncovered + b.cmp.trim3_uncovered + b.cmp.del_uncovered,
                b.coverage
            ),
            None => "\t\t\t\t\t\t".to_string(),
        };
        let any = m.values().next().unwrap();
        let scan_cells = match m.get("closest").filter(|b| !b.scan.is_empty()) {
            Some(b) => {
                let (k, id) = best_scan(&b.scan).unwrap();
                format!(
                    "{k}\t{id:.4}\t{:.4}\t{}",
                    id - b.cmp.identity(),
                    source_hint(&b.scan, b.cmp.identity())
                )
            }
            None => "\t\t\t".to_string(),
        };
        writeln!(
            tsv,
            "{fam}\t{}\t{}\t{}\t{}\t{scan_cells}",
            any.cmp.library_len,
            any.copies_available,
            cell(m.get("population")),
            cell(m.get("closest"))
        )?;
    }
    tsv.flush()?;
    eprintln!("wrote {prefix}-population.stk, {prefix}-closest.stk, {prefix}-summary.tsv");

    if let Some(path) = &args.subfamily_out {
        let mut fh = BufWriter::new(File::create(path).with_context(|| path.display().to_string())?);
        writeln!(
            fh,
            "family\tcall\tscored\tlibrary_like\tother\tbetween\tbetween_share\t\
             id_library_like\tcov_library_like\tid_other\tcov_other\tid_population\t\
             disputed\tdiagnostic\tdiag_start\tdiag_end\tdiag_span_frac"
        )?;
        let mut n = 0;
        for b in &built {
            let Some(x) = &b.sub else { continue };
            n += 1;
            let (lo, hi, frac) = x
                .span
                .map(|(l, h, f)| (l.to_string(), h.to_string(), format!("{f:.3}")))
                .unwrap_or_default();
            let v = |x: f64| if x.is_nan() { String::new() } else { format!("{x:.4}") };
            writeln!(
                fh,
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{lo}\t{hi}\t{frac}",
                b.family, x.call, x.scored, x.n_a, x.n_b, x.between,
                x.between as f64 / x.scored.max(1) as f64,
                v(x.id_a), v(x.cov_a), v(x.id_b), v(x.cov_b),
                b.cmp.identity(), x.sites, x.diagnostic,
            )?;
        }
        fh.flush()?;
        eprintln!("subfamily: {n} families -> {}", path.display());
    }
    Ok(())
}
