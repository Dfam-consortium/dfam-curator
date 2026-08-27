//! Align sequence 0 of a FASTA against every other, printing the same records
//! as the C++ `pairaln` harness so the two can be diffed directly.
//!
//! ```sh
//! cargo run -p cons-core --release --example pairaln -- fam.fa matrix.txt [mode]
//! ```
//!
//! `mode` is `local` (default), `sg` (both ends free), `fit` (subject ends free
//! — fit the whole instance into the reference), or `global`.

use aln_core::{io, Sequence, SubstMatrix};
use aln_engine::{AlignMode, AlignParams, PairwiseAligner};
use cons_core::FastAligner;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let fasta = args.next().ok_or("usage: pairaln <fasta> <matrix> [mode]")?;
    let matrix_path = args.next().ok_or("usage: pairaln <fasta> <matrix> [mode]")?;
    let mode = match args.next().as_deref().unwrap_or("local") {
        "local" => AlignMode::Local,
        "global" => AlignMode::Global,
        "sg" => AlignMode::SemiGlobal { free_query_ends: true, free_subject_ends: true },
        // Fit the whole instance into the reference: the instance must be used
        // in full, the reference's flanks are free.
        "fit" => AlignMode::SemiGlobal { free_query_ends: false, free_subject_ends: true },
        other => return Err(format!("unknown mode {other}").into()),
    };

    let seqs = io::read_fasta_file(&fasta)?;
    if seqs.len() < 2 {
        return Err("need >=2 sequences".into());
    }
    let matrix = SubstMatrix::from_file(&matrix_path)?;
    // Take the gap penalties from the matrix's own GAP line, as the C++ does —
    // hardcoding them is how this harness first went wrong, pairing a
    // max-match-3 matrix with a 30/6 gap cost.
    let params = AlignParams { mode, min_score: 1, ..AlignParams::from_matrix(&matrix) };
    eprintln!(
        "matrix {} gap_open={} gap_extend={} mode={:?}",
        matrix_path, params.gap_open, params.gap_extend, params.mode
    );
    let aligner = FastAligner::new(matrix, params)?;

    // The C++ makes the reference the *subject* (`bot`) and the instance the
    // query (`top`); mirror that.
    let reference: &Sequence = &seqs[0];
    let profile = aligner.prepare_subject(reference)?;

    for q in &seqs[1..] {
        match aligner.align_prepared(&profile, q)? {
            Some(a) => {
                let (gq, gs) = a.gapped(&q.seq, &reference.seq)?;
                println!("PAIR\t{}\t{}\t{:.1}\t{}", q.name, q.len(), a.score, gq.len());
                println!("Q\t{}", String::from_utf8_lossy(&gq));
                println!("S\t{}", String::from_utf8_lossy(&gs));
            }
            None => {
                println!("PAIR\t{}\t{}\t0.0\t0", q.name, q.len());
                println!("Q\t");
                println!("S\t");
            }
        }
    }
    Ok(())
}
