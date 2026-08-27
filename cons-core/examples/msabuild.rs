//! Build the star MSA from sequence 0 of a FASTA and dump it, in the same
//! format as the C++ `msabuild` harness so the two can be diffed.
//!
//! This isolates MSA *assembly* from consensus *calling*: the pairwise
//! alignments are already known to agree exactly (see `pairaln`), so any
//! difference here is `assemble_msa` against GIRI's `adjustReference`.
//!
//! ```sh
//! cargo run -p cons-core --release --example msabuild -- fam.fa matrix.txt
//! ```

use aln_core::msa::{assemble_msa, InsertionPolicy, MsaMember};
use aln_core::{io, SubstMatrix};
use aln_engine::{AlignMode, AlignParams, PairwiseAligner};
use cons_core::FastAligner;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let fasta = args.next().ok_or("usage: msabuild <fasta> <matrix> [policy]")?;
    let matrix_path = args.next().ok_or("usage: msabuild <fasta> <matrix> [policy]")?;
    let policy = match args.next().as_deref().unwrap_or("incremental") {
        "incremental" => InsertionPolicy::GrowIncremental,
        "per-slot" => InsertionPolicy::GrowPerSlot,
        "drop" => InsertionPolicy::Drop,
        o => return Err(format!("unknown policy {o}").into()),
    };
    // Incremental merging is order-dependent; let the caller vary it.
    let order = args.next().unwrap_or_else(|| "input".to_string());

    let seqs = io::read_fasta_file(&fasta)?;
    let matrix = SubstMatrix::from_file(&matrix_path)?;
    let params = AlignParams {
        mode: AlignMode::Local,
        min_score: 1,
        ..AlignParams::from_matrix(&matrix)
    };
    let aligner = FastAligner::new(matrix, params)?;

    let reference = &seqs[0];
    let profile = aligner.prepare_subject(reference)?;

    let mut score = 0i64;
    let mut rows = Vec::new();
    for (i, q) in seqs.iter().enumerate() {
        if i == 0 {
            continue; // the C++ skips aligning the reference to itself
        }
        if let Some(a) = aligner.align_prepared(&profile, q)? {
            score += a.score as i64;
            let (gq, gs) = a.gapped(&q.seq, &reference.seq)?;
            rows.push((i, gq, gs, a));
        }
    }

    match order.as_str() {
        "input" => {}
        "reverse" => rows.reverse(),
        // Longest first, which is a plausible stable ordering for a star merge.
        "longest" => rows.sort_by_key(|(i, _, _, _)| std::cmp::Reverse(seqs[*i].len())),
        "shortest" => rows.sort_by_key(|(i, _, _, _)| seqs[*i].len()),
        o => return Err(format!("unknown order {o}").into()),
    }

    let members: Vec<MsaMember<'_>> = rows
        .iter()
        .map(|(i, gq, gs, a)| MsaMember {
            name: &seqs[*i].name,
            gapped_query: gq,
            gapped_reference: gs,
            ref_start: a.subj_start,
            span: Some(aln_coord::Span::new(a.query_start as u64, a.query_end as u64).unwrap()),
            orient: a.strand,
        })
        .collect();

    let msa = assemble_msa(
        &reference.seq,
        &reference.name,
        &members,
        policy,
    )?;

    eprintln!(
        "score={score} rows={} width={}",
        msa.sequences.len(),
        msa.width()
    );
    for row in &msa.sequences {
        println!("ROW\t{}\t{}", row.name, String::from_utf8_lossy(&row.seq));
    }
    Ok(())
}
