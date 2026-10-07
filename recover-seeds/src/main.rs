//! Build seed alignments from a RepeatMasker .align file under two sampling
//! policies, and report how each consensus compares with the library's.
//!
//! This tool never realigns copies. The .align file is its only source of
//! alignments, so a change in output comes from the sampling rather than from
//! a different aligner.

mod subfamily;
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
                    family::top_up_coverage(
                        &inst,
                        &ranked,
                        &chosen,
                        consensus.len() as u64,
                        args.cover_gap,
                        args.cover_ends_only,
                    )
                } else {
                    chosen
                };
                let msa = match family::build_msa(fam, consensus, &inst, &chosen, insertions) {
                    Ok(m) => m,
                    Err(e) => {
                        eprintln!("{fam} ({label}): {e:#}");
                        continue;
                    }
                };
                let cons_params = ConsParams {
                    min_non_gap_count: args.min_non_gap,
                    ..ConsParams::default()
                };
                // generateSeedAlignments' alignment spans only its instances,
                // so reference positions outside them never become columns.
                // Trim to the covered extent before calling, or the flanks
                // would collect N and the consensus would keep the reference's
                // full length.
                let covered = covered_positions(consensus.len(), &inst, &chosen);
                let mut msa = msa;
                if let Some((lo, hi)) = covered_extent_columns(&msa, &covered) {
                    msa.slice_columns(lo, hi + 1);
                }
                let gapped = cons_core::gapped_consensus(&msa, &cons_params);
                let called: Vec<u8> = gapped
                    .iter()
                    .filter(|b| b.is_ascii_alphabetic())
                    .copied()
                    .collect();

                let coverage =
                    covered.iter().filter(|c| **c).count() as f64 / consensus.len() as f64;
                let cmp = match compare::compare(&aligner, consensus, &called, &covered) {
                    Ok(c) => c,
                    Err(e) => {
                        eprintln!("{fam} ({label}): {e:#}");
                        continue;
                    }
                };
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
                            let gapped = cons_core::gapped_consensus(&msa, &cons_params);
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
                let mut cmp = cmp;
                cmp.agreed = Vec::new();   // large, and nothing downstream reads it
                out.push(Built {
                    family: fam.clone(),
                    policy: label,
                    record: to_stk(
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
                    ),
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
) -> Option<f64> {
    let msa = family::build_msa(family, consensus, inst, chosen, insertions).ok()?;
    let params = ConsParams { min_non_gap_count: min_non_gap, ..ConsParams::default() };
    let gapped = cons_core::gapped_consensus(&msa, &params);
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
