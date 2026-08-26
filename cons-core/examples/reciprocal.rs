//! Is te-composer's phase-1 all-vs-all doing redundant reciprocal work?
//!
//! `all_vs_all` passes the same indexed set as both query and subject, so
//! rmblast computes every *ordered* pair — (i,j) and (j,i) — plus N self-pairs
//! that are discarded after the fact. If the two directions always agree, the
//! upper triangle plus a transpose would do, halving the work.
//!
//! Agreement is not guaranteed even with a symmetric matrix and equal
//! insertion/deletion penalties: rmblast is a seed-extend heuristic, and
//! seeding scans the query against an indexed subject, which is not a symmetric
//! operation.
//!
//! This reports, per unordered pair, whether (j,i) is the transpose of (i,j):
//! same HSP count, same scores, and same spans with the axes swapped.
//!
//! ```sh
//! cargo run -p cons-core --release --example reciprocal -- family.fa [matrix]
//! ```

use aln_core::{io, SubstMatrix};
use aln_engine::engine::ScoreMode;
use aln_engine::{AlignMode, AlignParams, SearchParams};
use aln_rmblast::{RmblastEngine, RmblastOptions};

/// One HSP reduced to what a transpose must preserve.
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone)]
struct Hsp {
    score: i64,
    /// Span on the sequence that was the *query* of this pair.
    q: (usize, usize),
    /// Span on the sequence that was the *subject* of this pair.
    s: (usize, usize),
    minus: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let fasta = args.next().ok_or("usage: reciprocal <family.fa> [matrix]")?;
    let matrix = match args.next() {
        Some(p) => SubstMatrix::from_file(&p)?,
        None => return Err("a matrix path is required".into()),
    };
    // Third arg "ca" turns on complexity-adjusted scoring, which is what
    // `Refiner` always uses and te-composer leaves off by default.
    // Remaining args are flags: "ca" = complexity-adjusted scoring (Refiner's
    // setting), "dust" = low-complexity query masking (which RepeatMasker turns
    // off with `-dust no`, so this is an experiment, not a fidelity fix).
    let flags: Vec<String> = args.collect();
    let complexity = flags.iter().any(|f| f == "ca");
    let dust = flags.iter().any(|f| f == "dust");

    let seqs = io::read_fasta_file(&fasta)?;
    let matrix_for_check = matrix.clone();
    let from_matrix = AlignParams::from_matrix(&matrix);
    let align = AlignParams {
        mode: AlignMode::Local,
        min_score: 150,
        traceback: true,
        ..from_matrix
    };

    // The same settings te-composer's rmblast backend uses.
    let search = SearchParams {
        matrix: Some(matrix),
        gap_init: -(align.gap_open as i32),
        ins_gap_ext: -(align.gap_extend as i32),
        del_gap_ext: -(align.gap_extend as i32),
        min_match: 7,
        min_score: align.min_score.max(0),
        mask_level: 101,
        cores: Some(std::thread::available_parallelism()?.get()),
        score_mode: if complexity {
            ScoreMode::ComplexityAdjusted
        } else {
            ScoreMode::Basic
        },
        ..Default::default()
    };
    // Diagnostics: complexity adjustment short-circuits to identity when the
    // matrix carries no background frequencies (`lambda <= 0.0`), so confirm
    // both the flag and lambda rather than inferring from output.
    let rm_matrix = aln_rmblast::matrix::to_rmblast(&matrix_for_check)?;
    eprintln!("matrix lambda = {}", rm_matrix.lambda);
    let engine = RmblastEngine::new(
        search,
        RmblastOptions { dust, ..RmblastOptions::default() },
    )?;
    eprintln!(
        "engine complexity_adjust = {}",
        engine.rmblast_params().complexity_adjust
    );

    // "selfprobe": align each sequence against itself and count HSPs.
    //
    // A sequence that is not internally repetitive aligns to itself as a single
    // diagonal HSP. A tandemly repetitive one matches itself off-diagonal at
    // every period, so the count explodes. `all_vs_all` already computes every
    // self-pair and throws it away (`if qi == si { continue }`), so this signal
    // is free to collect and is a candidate pre-flight check: the all-vs-all is
    // what blows up, and this predicts it from N alignments instead of N^2.
    if flags.iter().any(|f| f == "selfprobe") {
        let mut rows: Vec<(usize, String, usize)> = Vec::new();
        let t = std::time::Instant::now();
        for (i, s) in seqs.iter().enumerate() {
            let one = vec![s.clone()];
            let hits = engine.one_to_many(s, &one, None)?;
            rows.push((s.len(), s.name.clone(), hits.len()));
            let _ = i;
        }
        let elapsed = t.elapsed().as_secs_f64();
        rows.sort_by(|a, b| b.2.cmp(&a.2));
        let counts: Vec<usize> = rows.iter().map(|r| r.2).collect();
        let n = counts.len();
        let mut sorted = counts.clone();
        sorted.sort_unstable();
        println!("self-alignment HSPs per sequence ({n} seqs, {elapsed:.1}s)");
        println!("  median {}  max {}  mean {:.1}",
                 sorted[n / 2], sorted[n - 1],
                 counts.iter().sum::<usize>() as f64 / n as f64);
        println!("  worst 5:");
        for (len, name, c) in rows.iter().take(5) {
            println!("    {name:<20} {len:7} bp  {c} self-HSPs");
        }
        return Ok(());
    }

    eprintln!("{} sequences; complexity-adjusted={complexity} dust={dust}; running all-vs-all...", seqs.len());
    let t0 = std::time::Instant::now();
    let hits = engine.all_vs_all(&seqs, 101)?;
    eprintln!(
        "all-vs-all: {} HSPs over {} ordered pairs in {:.1}s",
        hits.len(),
        seqs.len() * (seqs.len() - 1),
        t0.elapsed().as_secs_f64()
    );

    // Bucket HSPs by ordered pair.
    let mut by: std::collections::HashMap<(usize, usize), Vec<Hsp>> =
        std::collections::HashMap::new();
    for (qi, si, a) in hits {
        by.entry((qi, si)).or_default().push(Hsp {
            score: a.score as i64,
            q: (a.query_start, a.query_end),
            s: (a.subj_start, a.subj_end),
            minus: a.strand.is_minus(),
        });
    }
    for v in by.values_mut() {
        v.sort();
    }

    let (mut both, mut one_only, mut identical, mut same_count, mut same_scores) =
        (0usize, 0usize, 0usize, 0usize, 0usize);
    let mut score_delta_bp: Vec<i64> = Vec::new();

    let n = seqs.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let f = by.get(&(i, j));
            let r = by.get(&(j, i));
            match (f, r) {
                (None, None) => {}
                (Some(_), None) | (None, Some(_)) => one_only += 1,
                (Some(f), Some(r)) => {
                    both += 1;
                    if f.len() == r.len() {
                        same_count += 1;
                    }
                    // Transpose the reverse direction back into forward frame.
                    let mut rt: Vec<Hsp> = r
                        .iter()
                        .map(|h| Hsp { score: h.score, q: h.s, s: h.q, minus: h.minus })
                        .collect();
                    rt.sort();
                    if *f == rt {
                        identical += 1;
                    }
                    let fs: i64 = f.iter().map(|h| h.score).sum();
                    let rs: i64 = r.iter().map(|h| h.score).sum();
                    if fs == rs {
                        same_scores += 1;
                    }
                    score_delta_bp.push((fs - rs).abs());
                }
            }
        }
    }

    let pairs = both + one_only;
    println!("\nunordered pairs with any HSP: {pairs}");
    println!("  both directions present : {both}");
    println!("  one direction only      : {one_only}");
    if both > 0 {
        let pct = |x: usize| 100.0 * x as f64 / both as f64;
        println!("  same HSP count          : {same_count} ({:.1}%)", pct(same_count));
        println!("  identical after transpose: {identical} ({:.1}%)", pct(identical));
        println!("  same summed score       : {same_scores} ({:.1}%)", pct(same_scores));
        score_delta_bp.sort_unstable();
        let med = score_delta_bp[score_delta_bp.len() / 2];
        let max = score_delta_bp.last().copied().unwrap_or(0);
        println!("  |score difference|: median {med}, max {max}");
    }
    Ok(())
}
