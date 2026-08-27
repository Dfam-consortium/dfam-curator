//! Score a consensus against the family it came from.
//!
//! Aligns the consensus to every input sequence and sums the best score per
//! sequence. A consensus that represents the family well pulls a high total;
//! one that is truncated, chimeric, or padded with noise does not.
//!
//! The scorer is fixed (parasail, local, one matrix) regardless of which
//! backend produced the consensus — otherwise this would compare scoring
//! systems rather than consensus quality.
//!
//! cargo run --release -p cons-core --example sumscore -- <cons.fa> <family.fa> <matrix>

use aln_core::{io, SubstMatrix};
use aln_engine::{AlignMode, AlignParams, PairwiseAligner};
use cons_core::FastAligner;

fn read(p: &str) -> Result<Vec<aln_core::Sequence>, Box<dyn std::error::Error>> {
    Ok(io::read_fasta(std::io::BufReader::new(std::fs::File::open(p)?))?)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = std::env::args().skip(1);
    let (cons_p, fam_p, mx_p) = (
        a.next().expect("consensus"),
        a.next().expect("family"),
        a.next().expect("matrix"),
    );

    let cons = read(&cons_p)?;
    let fam = read(&fam_p)?;
    let matrix = SubstMatrix::parse(&std::fs::read_to_string(&mx_p)?)?;

    // An absent or empty consensus scores zero rather than erroring, so a
    // backend that produced nothing is comparable with one that did.
    let Some(cons) = cons.into_iter().find(|s| !s.is_empty()) else {
        println!("0\t0\t0\t0");
        return Ok(());
    };

    let mut p = AlignParams::from_matrix(&matrix);
    p.mode = AlignMode::Local;
    p.min_score = 1;
    p.traceback = false; // score-only kernels; nothing here needs the path
    let al = FastAligner::new(matrix, p)?;

    let mut total: i64 = 0;
    let mut hit = 0usize;
    for s in &fam {
        if let Some(r) = al.align(&cons, s)? {
            total += r.score as i64;
            hit += 1;
        }
    }
    // total, mean over all inputs, sequences aligned, consensus length
    println!(
        "{total}\t{:.1}\t{hit}\t{}",
        total as f64 / fam.len().max(1) as f64,
        cons.len()
    );
    Ok(())
}
