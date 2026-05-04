/// Build a `MultiAlign` from a collection of pairwise alignments via transitive
/// alignment.
///
/// This replicates the Perl `MultAln::_alignFromSearchResultCollection()` method.
///
/// **Algorithm overview:**
///
/// Given N pairwise alignments of a *reference* sequence against N *instance*
/// sequences (each alignment covers a potentially different region of the
/// reference), we:
///
/// 1. Collect all gapped reference sub-sequences and determine the full reference
///    span (`ref_start` .. `ref_end`).
/// 2. Reconstruct the ungapped combined reference over that span.
/// 3. For each alignment position in the reference, find the maximum number of
///    inserted gaps that any alignment introduces at that position — the **gap
///    pattern**.
/// 4. Apply the merged gap pattern to the combined reference to produce the final
///    gapped reference row.
/// 5. For each instance, replay its alignment against the gapped reference to
///    produce its gapped instance row (inserting extra dashes wherever the merged
///    gap pattern is wider than that instance's gap pattern).
///
/// The result is a `MultiAlign` where every row has the same width.
use crate::alignment::{MultiAlign, Orientation, SequenceRow};
use crate::io::crossmatch::PairwiseHit;

/// Which side of each `PairwiseHit` is the shared reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reference {
    /// The *subject* field of each hit is the reference (typical for RepeatMasker
    /// output where the repeat consensus is the subject).
    Subject,
    /// The *query* field of each hit is the reference (less common).
    Query,
}

/// Build a `MultiAlign` from a set of pairwise hits.
///
/// `ref_side`: which column of each hit is the shared reference sequence.
/// `provided_ref_seq`: if `Some`, use this as the reference sequence instead of
///   reconstructing it from the hits.  It must already be ungapped and span
///   positions 1..=len.
pub fn build_from_pairwise(
    hits: &[PairwiseHit],
    ref_side: Reference,
    provided_ref_seq: Option<&[u8]>,
) -> Result<MultiAlign, BuildError> {
    if hits.is_empty() {
        return Err(BuildError::NoHits);
    }

    // ── Step 1: determine reference span and reconstruct combined reference ───
    let (ref_min, _ref_max, combined_ref) =
        reconstruct_reference(hits, ref_side, provided_ref_seq)?;

    let ref_len = combined_ref.len(); // ungapped length

    // ── Step 2: compute per-alignment gap patterns for the reference side ─────
    //
    // gap_pattern[i][j] = number of gap characters inserted in alignment i
    //                     *before* (to the left of) reference position j.
    //
    // Reference position j here is 0-based within the combined reference
    // (i.e., j == 0 is ref_min).
    let mut gap_patterns: Vec<Vec<usize>> = Vec::with_capacity(hits.len());

    for hit in hits {
        let ref_seq_gapped = ref_gapped_seq(hit, ref_side);
        let ref_start_offset = ref_start(hit, ref_side) as usize - ref_min;

        // Compute gap count before each ungapped reference base.
        // split_on_non_gap yields the gap runs between reference bases.
        let gaps = gaps_before_each_base(ref_seq_gapped);

        // Pad the front with zeros for the offset into the combined reference.
        let mut pattern = vec![0usize; ref_start_offset];
        pattern.extend_from_slice(&gaps);

        gap_patterns.push(pattern);
    }

    // ── Step 3: merge gap patterns — take the maximum at each position ────────
    let mut merged_gaps = vec![0usize; ref_len + 1]; // +1 for trailing gaps
    for pattern in &gap_patterns {
        for (j, &g) in pattern.iter().enumerate() {
            if j < merged_gaps.len() && g > merged_gaps[j] {
                merged_gaps[j] = g;
            }
        }
    }

    // ── Step 4: build the gapped reference row ────────────────────────────────
    let gapped_ref = build_gapped_seq(&combined_ref, &merged_gaps);
    let ref_width = gapped_ref.len();

    // Cumulative gap counts — used to compute gapped start positions.
    // total_gaps[j] = sum of merged_gaps[0..=j]
    let mut total_gaps = vec![0usize; ref_len + 1];
    total_gaps[0] = merged_gaps[0];
    for j in 1..=ref_len {
        total_gaps[j] = total_gaps[j - 1] + merged_gaps[j];
    }

    // ── Step 5: build each instance row ──────────────────────────────────────
    let mut instance_rows: Vec<SequenceRow> = Vec::with_capacity(hits.len());

    for (i, hit) in hits.iter().enumerate() {
        let inst_name = inst_name(hit, ref_side).to_string();
        let inst_seq_gapped = inst_gapped_seq(hit, ref_side);
        let ref_seq_gapped = ref_gapped_seq(hit, ref_side);
        let orient = if ref_side == Reference::Subject {
            hit.orientation
        } else {
            Orientation::Forward
        };

        let ref_start_offset = ref_start(hit, ref_side) as usize - ref_min;

        // Gapped column where this instance's aligned region begins.
        // We want the column just before the merged_gaps[ref_start_offset] insertion
        // columns at that position; the walk will then emit those gap chars.
        // total_gaps[j] = sum(merged_gaps[0..=j]), so we need sum(merged_gaps[0..j]),
        // which is total_gaps[j-1].  For j=0 there are no leading columns.
        let gapped_start = if ref_start_offset > 0 {
            ref_start_offset + total_gaps[ref_start_offset - 1]
        } else {
            0
        };

        // Build the gapped instance sequence by walking both the alignment's
        // reference and instance strings in lockstep.
        let mut inst_out: Vec<u8> = vec![b' '; gapped_start]; // leading padding
        let mut k = ref_start_offset; // position in ungapped combined reference
        let instance_gap_pattern = &gap_patterns[i];

        let ref_bytes = ref_seq_gapped.as_bytes();
        let inst_bytes = inst_seq_gapped.as_bytes();
        let aln_len = ref_bytes.len().min(inst_bytes.len());

        let mut ai = 0usize;
        while ai < aln_len {
            let ref_char = ref_bytes[ai];
            let inst_char = inst_bytes[ai];

            if ref_char != b'-' {
                // Reference has a base at this alignment column.
                // Insert any *extra* gaps required by the merged pattern.
                let local_gaps = instance_gap_pattern.get(k).copied().unwrap_or(0);
                let extra = merged_gaps[k].saturating_sub(local_gaps);
                for _ in 0..extra {
                    inst_out.push(b'-');
                }
                k += 1;
            }
            inst_out.push(inst_char);
            ai += 1;
        }

        // Handle trailing gaps in merged pattern (after the last reference base).
        let extra_trailing = merged_gaps
            .get(k)
            .copied()
            .unwrap_or(0)
            .saturating_sub(instance_gap_pattern.get(k).copied().unwrap_or(0));
        for _ in 0..extra_trailing {
            inst_out.push(b'-');
        }

        // Pad to full width with spaces.
        while inst_out.len() < ref_width {
            inst_out.push(b' ');
        }

        let inst_start_abs = inst_start(hit, ref_side);
        let inst_end_abs = inst_end(hit, ref_side);

        let mut row = SequenceRow::new(inst_name, inst_out);
        row.seq_start = inst_start_abs;
        row.seq_end = inst_end_abs;
        row.orient = orient;
        row.div = Some(hit.pct_div);
        row.src_div = Some(hit.pct_div);
        instance_rows.push(row);
    }

    // ── Assemble MultiAlign ───────────────────────────────────────────────────
    let ref_name = ref_name(hits.first().unwrap(), ref_side).to_string();
    let mut ref_row = SequenceRow::new(ref_name, gapped_ref);
    ref_row.seq_start = ref_min as u64;
    Ok(MultiAlign::from_sequences(ref_row, instance_rows))
}

#[derive(Debug, thiserror::Error)]
pub enum BuildError {
    #[error("no pairwise hits provided")]
    NoHits,
    #[error("gap in reference coverage: position {0} has no alignment data")]
    CoverageGap(usize),
    #[error("reference sequence length mismatch")]
    RefLenMismatch,
}

// ── Reference reconstruction ──────────────────────────────────────────────────

/// Determine (ref_min, ref_max) and reconstruct the combined ungapped reference.
fn reconstruct_reference(
    hits: &[PairwiseHit],
    ref_side: Reference,
    provided: Option<&[u8]>,
) -> Result<(usize, usize, Vec<u8>), BuildError> {
    // Find the span of the reference covered by all hits.
    let ref_min = hits
        .iter()
        .map(|h| ref_start(h, ref_side) as usize)
        .min()
        .unwrap();
    let ref_max = hits
        .iter()
        .map(|h| ref_end(h, ref_side) as usize)
        .max()
        .unwrap();

    if let Some(seq) = provided {
        return Ok((1, seq.len(), seq.to_vec()));
    }

    let span = ref_max - ref_min + 1;
    let mut combined = vec![b' '; span];

    for hit in hits {
        let gapped = ref_gapped_seq(hit, ref_side);
        let ungapped: Vec<u8> = gapped.bytes().filter(|&b| b != b'-').collect();
        let offset = ref_start(hit, ref_side) as usize - ref_min;
        let end = offset + ungapped.len();
        if end > span {
            return Err(BuildError::RefLenMismatch);
        }
        combined[offset..end].copy_from_slice(&ungapped);
    }

    // Check for gaps in coverage.
    for (i, &b) in combined.iter().enumerate() {
        if b == b' ' {
            return Err(BuildError::CoverageGap(ref_min + i));
        }
    }

    Ok((ref_min, ref_max, combined))
}

/// Number of gap characters inserted before each ungapped base in `gapped_seq`.
/// Returns a Vec of length == ungapped_length.
fn gaps_before_each_base(gapped_seq: &str) -> Vec<usize> {
    let mut result = Vec::new();
    let mut gap_run = 0usize;
    for b in gapped_seq.bytes() {
        if b == b'-' {
            gap_run += 1;
        } else {
            result.push(gap_run);
            gap_run = 0;
        }
    }
    // Trailing gaps (after last base) stored at position == ungapped_length.
    result.push(gap_run);
    result
}

/// Build a gapped sequence from an ungapped sequence and a gap pattern.
fn build_gapped_seq(ungapped: &[u8], gaps: &[usize]) -> Vec<u8> {
    let mut out = Vec::with_capacity(ungapped.len() + gaps.iter().sum::<usize>());
    for (i, &b) in ungapped.iter().enumerate() {
        let g = gaps.get(i).copied().unwrap_or(0);
        for _ in 0..g {
            out.push(b'-');
        }
        out.push(b);
    }
    // Trailing gaps after last base.
    if let Some(&trailing) = gaps.last() {
        // Only append if gaps has an extra element beyond ungapped.len().
        if gaps.len() > ungapped.len() {
            for _ in 0..trailing {
                out.push(b'-');
            }
        }
    }
    out
}

// ── Field accessors ───────────────────────────────────────────────────────────

fn ref_gapped_seq<'a>(hit: &'a PairwiseHit, side: Reference) -> &'a str {
    let bytes = match side {
        Reference::Subject => &hit.subj_seq,
        Reference::Query   => &hit.query_seq,
    };
    std::str::from_utf8(bytes).unwrap_or("")
}

fn inst_gapped_seq<'a>(hit: &'a PairwiseHit, side: Reference) -> &'a str {
    let bytes = match side {
        Reference::Subject => &hit.query_seq,
        Reference::Query   => &hit.subj_seq,
    };
    std::str::from_utf8(bytes).unwrap_or("")
}

fn ref_start(hit: &PairwiseHit, side: Reference) -> u64 {
    match side {
        Reference::Subject => hit.subj_start,
        Reference::Query   => hit.query_start,
    }
}

fn ref_end(hit: &PairwiseHit, side: Reference) -> u64 {
    match side {
        Reference::Subject => hit.subj_end,
        Reference::Query   => hit.query_end,
    }
}

fn inst_start(hit: &PairwiseHit, side: Reference) -> u64 {
    match side {
        Reference::Subject => hit.query_start,
        Reference::Query   => hit.subj_start,
    }
}

fn inst_end(hit: &PairwiseHit, side: Reference) -> u64 {
    match side {
        Reference::Subject => hit.query_end,
        Reference::Query   => hit.subj_end,
    }
}

fn ref_name<'a>(hit: &'a PairwiseHit, side: Reference) -> &'a str {
    match side {
        Reference::Subject => &hit.subj_name,
        Reference::Query   => &hit.query_name,
    }
}

fn inst_name<'a>(hit: &'a PairwiseHit, side: Reference) -> &'a str {
    match side {
        Reference::Subject => &hit.query_name,
        Reference::Query   => &hit.subj_name,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::Orientation;
    use crate::io::crossmatch::PairwiseHit;

    fn make_hit(
        subj_start: u64, subj_end: u64, subj_seq: &[u8],
        query_start: u64, query_end: u64, query_seq: &[u8],
    ) -> PairwiseHit {
        PairwiseHit {
            sw_score: 100,
            pct_div: 5.0,
            pct_del: 0.0,
            pct_ins: 0.0,
            query_name: "inst".to_string(),
            query_start,
            query_end,
            query_remaining: 0,
            query_seq: query_seq.to_vec(),
            subj_name: "ref".to_string(),
            subj_start,
            subj_end,
            subj_remaining: 0,
            subj_seq: subj_seq.to_vec(),
            orientation: Orientation::Forward,
            id: None,
        }
    }

    #[test]
    fn trivial_single_hit_no_gaps() {
        // Single hit, no gaps in either sequence — trivial alignment.
        let hit = make_hit(1, 4, b"ACGT", 1, 4, b"ACGT");
        let msa = build_from_pairwise(&[hit], Reference::Subject, None).unwrap();
        assert_eq!(msa.width(), 4);
        assert_eq!(msa.reference().unwrap().seq, b"ACGT");
        assert_eq!(msa.instance(0).unwrap().seq, b"ACGT");
    }

    #[test]
    fn single_hit_with_gap_in_instance() {
        // Reference: ACGT (no gap), instance: AC-T (1 deletion).
        let hit = make_hit(1, 4, b"ACGT", 1, 3, b"AC-T");
        let msa = build_from_pairwise(&[hit], Reference::Subject, None).unwrap();
        // Reference should stay ungapped (no other alignments to force gaps).
        assert_eq!(msa.reference().unwrap().seq, b"ACGT");
        assert_eq!(msa.instance(0).unwrap().seq, b"AC-T");
    }

    #[test]
    fn two_hits_merged_gaps() {
        // Reference ACGT.
        // Hit 1: ref AC-GT, inst ACCGT (insertion at pos 2 in ref).
        // Hit 2: ref ACGT,  inst ACGT  (no gaps).
        // Merged ref must accommodate the insertion in hit 1 → AC-GT.
        let h1 = make_hit(1, 4, b"AC-GT", 1, 5, b"ACCGT");
        let h2 = make_hit(1, 4, b"ACGT",  1, 4, b"ACGT");
        let msa = build_from_pairwise(&[h1, h2], Reference::Subject, None).unwrap();
        // The gapped reference must include a gap column.
        assert!(msa.width() > 4, "merged ref must be wider than ungapped");
        // Instance 2 should have a gap inserted to match the merged width.
        let inst2 = &msa.instance(1).unwrap().seq;
        assert!(inst2.contains(&b'-'), "instance 2 must have a gap at merged column");
    }
}
