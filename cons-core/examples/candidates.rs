//! Dump the phase-1 candidate ranking, for comparison against the C++.
//!
//! The C++ can be instrumented to emit the same shape by adding a line to
//! `process_locally` after its `autoConsensus` call:
//!
//! ```c
//! fprintf(stderr, "CAND\t%d\t%s\t%zu\t%.3f\t%zu\n", cnt,
//!         (*slsptr1).getName().c_str(), (*slsptr1).size(), score, sCON.size());
//! ```
//!
//! Run it with `--nothreads` so the candidate order — and hence the C++'s
//! `score += 0.001` tie-bumping — is deterministic.
//!
//! ```sh
//! cargo run -p cons-core --release --example candidates -- family.fa matrix.txt
//! ```

use aln_core::consensus::ConsensusParams;
use aln_core::{io, SubstMatrix};
use aln_engine::{AlignMode, AlignParams};
use aln_parasail::ParasailAligner;
use cons_core::{score_candidates, Caller, Pairwise, Params};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let fasta = args.next().ok_or("usage: candidates <family.fa> <matrix>")?;
    let matrix_path = args.next().ok_or("usage: candidates <family.fa> <matrix>")?;

    let seqs = io::read_fasta_file(&fasta)?;
    let matrix = SubstMatrix::from_file(&matrix_path)?;

    // Penalties from the matrix's own GAP line — hardcoding them here once
    // produced a bogus 3x score gap against the C++.
    let params = AlignParams {
        mode: AlignMode::Local,
        min_score: 1,
        ..AlignParams::from_matrix(&matrix)
    };
    eprintln!("gap_open={} gap_extend={}", params.gap_open, params.gap_extend);
    let aligner = ParasailAligner::new(matrix, params)?;

    let params = Params {
        caller: Caller::Giri,
        consensus: ConsensusParams::default(),
        ..Default::default()
    };

    let cands = score_candidates(&Pairwise::new(aligner), &seqs, &params)?;

    // Same columns as the instrumented C++: index, name, ref length, score,
    // consensus length.  Emitted in input order so the two dumps line up.
    let mut by_index: Vec<_> = cands.iter().collect();
    by_index.sort_by_key(|c| c.index);
    for c in &by_index {
        println!(
            "CAND\t{}\t{}\t{}\t{}\t{}",
            c.index,
            seqs[c.index].name,
            seqs[c.index].len(),
            c.score,
            c.consensus.len()
        );
    }

    if let Some(best) = cands.first() {
        eprintln!(
            "WINNER index={} name={} ref_len={} score={} cons_len={}",
            best.index,
            seqs[best.index].name,
            seqs[best.index].len(),
            best.score,
            best.consensus.len()
        );
    }
    Ok(())
}
