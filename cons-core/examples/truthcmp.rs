//! Score a consensus against a known progenitor.
//!
//! Built for the `dfam-curator` `cpg-topology` forward-simulation benchmark,
//! where each simulated family carries its own root sequence (`int-node-0`), so
//! the right answer is known. That makes this an *absolute* accuracy measure,
//! unlike the sum-score used elsewhere in this comparison — which is
//! self-referential and rewards a longer consensus.
//!
//! Prints one TSV line:
//! `identity_pct  coverage_pct  cons_len  truth_len  score  self_score  norm`
//!
//!   * identity — matched columns / aligned columns
//!   * coverage — progenitor positions aligned to a consensus base, over its length
//!   * norm     — `score / self_score`, the fraction of the achievable score
//!
//! # Why the normalised score
//!
//! Identity and coverage have to be read together, and they trade off
//! invisibly: a short pristine fragment is 100% identical at low coverage, a
//! sprawling consensus covers everything at poor identity, and two arms can
//! split the pair with no way to say which is better. Under `global` the
//! alignment score folds both into one number — surplus or missing consensus
//! length costs gap penalties, base errors cost mismatches.
//!
//! Raw score scales with length, so comparing a 311 bp AluY against a 3 kb L2
//! needs normalising. `self_score` is the progenitor aligned to itself — the
//! sum of its own diagonal entries, which is the most any consensus could
//! score — so `norm` is 1.0 for a perfect reconstruction and falls from there.
//! That makes it comparable across families of different length and base
//! composition.
//!
//! `norm` can go negative: a consensus bad enough that mismatch and gap
//! penalties outweigh its matches has a negative global score, which is
//! meaningful rather than an error.
//!
//! # Mode
//!
//! `local` (the default) keeps the historical behaviour and ignores unaligned
//! flanks entirely — which is exactly why it cannot see length differences.
//! `global` is the one to use with `norm`.
//!
//! cargo run --release -p cons-core --example truthcmp -- <truth.fa> <cons.fa> <matrix> [mode]

use aln_core::{io, SubstMatrix};
use aln_engine::{AlignMode, AlignParams, PairwiseAligner};
use cons_core::FastAligner;

fn read_one(path: &str) -> Result<Option<aln_core::Sequence>, Box<dyn std::error::Error>> {
    let f = std::fs::File::open(path)?;
    let mut seqs = io::read_fasta(std::io::BufReader::new(f))?;
    seqs.retain(|s| !s.is_empty());
    Ok(seqs.into_iter().next())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = std::env::args().skip(1);
    let truth_p = a.next().expect("truth.fa");
    let cons_p = a.next().expect("cons.fa");
    let mx_p = a.next().expect("matrix");

    let truth = read_one(&truth_p)?.expect("truth sequence must be present");
    // A backend that produced nothing scores zero rather than erroring, so an
    // empty run stays comparable with a populated one.
    let mode = match a.next().as_deref().unwrap_or("local") {
        "local" => AlignMode::Local,
        "global" => AlignMode::Global,
        other => return Err(format!("unknown mode {other}; want local or global").into()),
    };

    let matrix = SubstMatrix::parse(&std::fs::read_to_string(&mx_p)?)?;

    // The most any consensus could score: the progenitor aligned to itself,
    // which under any sane matrix is its own diagonal summed. Computed
    // analytically rather than by aligning truth to truth, so it cannot be
    // perturbed by the aligner's tie-breaking.
    let self_score: i64 = truth
        .seq
        .iter()
        .map(|&b| matrix.score(b, b).unwrap_or(0) as i64)
        .sum();

    let Some(cons) = read_one(&cons_p)? else {
        println!("0.00\t0.00\t0\t{}\t0\t{self_score}\t0.0000\t0", truth.len());
        return Ok(());
    };

    let mut p = AlignParams::from_matrix(&matrix);
    p.mode = mode;
    // In global mode the score of a poor consensus is legitimately negative —
    // the doc above promises `norm` can go negative, and a min_score of 1 was
    // silently converting those alignments to None, recorded as 0/0/0.0000.
    // Local mode keeps the historical floor: a local alignment's score is
    // non-negative by construction, so 1 only suppresses the empty alignment.
    p.min_score = match mode {
        AlignMode::Global => i32::MIN / 4,
        _ => 1,
    };
    p.traceback = true;
    let al = FastAligner::new(matrix, p)?;

    // Truth as subject so coverage is measured against the truth.
    match al.align(&cons, &truth)? {
        Some(r) => {
            // Count directly from the gapped strings rather than through
            // identity_counts, because the mismatch bucket needs splitting:
            // a consensus `N` over a real root base is declared uncertainty,
            // not a wrong call, and the two must be reportable separately —
            // the callers under comparison emit N at different rates, so
            // pooling them biases the substitution column by N policy.
            let (gq, gs) = r.gapped(&cons.seq, &truth.seq)?;
            let mut matches = 0u32;
            let mut wrong = 0u32;
            let mut n_cols = 0u32;
            for (&qb, &sb) in gq.iter().zip(&gs) {
                if aln_core::seq::is_gap(qb) || aln_core::seq::is_gap(sb) {
                    continue;
                }
                if qb.eq_ignore_ascii_case(&sb) {
                    matches += 1;
                } else if qb == b'N' || qb == b'n' {
                    n_cols += 1;
                } else {
                    wrong += 1;
                }
            }
            let aligned = matches + wrong + n_cols;
            let ident = if aligned == 0 {
                0.0
            } else {
                100.0 * matches as f64 / aligned as f64
            };
            // Truth positions carrying an aligned consensus base. Not the
            // subject *span*, which counts deletions inside the alignment as
            // covered — under `global` the span is the whole sequence by
            // construction and would always read 100%.
            let cov = 100.0 * aligned as f64 / truth.len() as f64;
            let norm = if self_score != 0 {
                r.score as f64 / self_score as f64
            } else {
                0.0
            };
            println!(
                "{ident:.2}\t{cov:.2}\t{}\t{}\t{}\t{self_score}\t{norm:.4}\t{n_cols}",
                cons.len(),
                truth.len(),
                r.score
            );
        }
        None => println!(
            "0.00\t0.00\t{}\t{}\t0\t{self_score}\t0.0000\t0",
            cons.len(),
            truth.len()
        ),
    }
    Ok(())
}
