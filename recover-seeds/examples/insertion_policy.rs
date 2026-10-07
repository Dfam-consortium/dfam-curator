//! Show how the two insertion policies lay out the same three alignments.
//!
//! Reference ACGTACGT.  Three members, each aligned to part of it:
//!   A: inserts TT after reference base 4
//!   B: inserts G  after reference base 4
//!   C: no insertion, and it stops before base 4

use aln_core::msa::{assemble_msa, InsertionPolicy, MsaMember};
use aln_core::seq::Strand;

fn show(policy: InsertionPolicy, label: &str, reference: &[u8], members: &[MsaMember]) {
    let msa = assemble_msa(reference, "ref", members, policy).expect("assembles");
    println!("\n{label}  (width {})", msa.width());
    for row in &msa.sequences {
        println!("  {:<4} {}", row.name, String::from_utf8_lossy(&row.seq).replace(' ', "·"));
    }
}

fn main() {
    let reference = b"ACGTACGT";
    //            ref: A C G T - - A C G T
    let a_query = b"ACGTTTACGT";
    let a_ref = b"ACGT--ACGT";
    //            ref: A C G T - A C G T
    let b_query = b"ACGTGACGT";
    let b_ref = b"ACGT-ACGT";
    // C covers only the first three bases.
    let c_query = b"ACG";
    let c_ref = b"ACG";

    let members = vec![
        MsaMember { name: "A", gapped_query: a_query, gapped_reference: a_ref,
                    ref_start: 0, span: None, orient: Strand::Plus },
        MsaMember { name: "B", gapped_query: b_query, gapped_reference: b_ref,
                    ref_start: 0, span: None, orient: Strand::Plus },
        MsaMember { name: "C", gapped_query: c_query, gapped_reference: c_ref,
                    ref_start: 0, span: None, orient: Strand::Plus },
    ];
    println!("reference {}", String::from_utf8_lossy(reference));
    println!("(· marks padding, - marks a gap)");
    show(InsertionPolicy::GrowPerSlot, "GrowPerSlot (what the Perl does)", reference, &members);
    show(InsertionPolicy::GrowIncremental, "GrowIncremental (GIRI's merge)", reference, &members);

    // Order matters to the incremental merge, so run it reversed as well.
    let reversed: Vec<MsaMember> = members.iter().rev().cloned().collect();
    show(InsertionPolicy::GrowIncremental, "GrowIncremental, members reversed", reference, &reversed);
    show(InsertionPolicy::GrowPerSlot, "GrowPerSlot, members reversed", reference, &reversed);
}
