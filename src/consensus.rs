/// Parameters for the consensus-calling algorithm.
///
/// Defaults match the Perl MultAln.pm buildConsensusFromArray defaults.
#[derive(Debug, Clone)]
pub struct ConsensusParams {
    /// Include the reference sequence in the profile (inclRef in Perl).
    pub include_reference: bool,

    /// Bonus applied when a CpG dinucleotide hypothesis wins over the
    /// plain dinucleotide score.  Observation TG or CA each adds this.
    /// Default: 12.
    pub cg_param: i32,

    /// Penalty applied when the observed dinucleotide is TA (two-step CpG
    /// deamination product).  Default: -5.
    pub ta_param: i32,

    /// Transition bonus applied when the observed dinucleotide is a
    /// transition pair (TC/TT on forward or AA/GA on reverse) that could
    /// have arisen from a CpG.  Default: 2.
    pub cg_trans_param: i32,
}

impl Default for ConsensusParams {
    fn default() -> Self {
        ConsensusParams {
            include_reference: false,
            cg_param: 12,
            ta_param: -5,
            cg_trans_param: 2,
        }
    }
}

/// Build a gapped consensus from raw (gapped) aligned sequences.
///
/// This is the primary entry point.  It replicates the Perl
/// `buildConsensusFromArray` function exactly, including the per-sequence
/// CpG dinucleotide correction pass.
///
/// Each element of `sequences` is a byte slice of the same length (gapped
/// alignment width).  Space bytes (b' ') denote flanking padding and are
/// excluded from scoring.  Gap bytes (b'-') within the covered region are
/// included and can win at sparsely-covered columns.
pub fn build_consensus_from_sequences(
    sequences: &[&[u8]],
    params: &ConsensusParams,
) -> Vec<u8> {
    use crate::matrix::{alpha_byte, alpha_idx, ALPHA_LEN, MATRIX};

    if sequences.is_empty() {
        return Vec::new();
    }

    let width = sequences.iter().map(|s| s.len()).max().unwrap_or(0);

    // Pre-process: convert leading and trailing gap characters to spaces.
    // This matches Perl's MultAln::consensus() which pads each sequence with
    // spaces before the first base and does NOT extend past the last base.
    // Spaces are skipped when building the profile; internal gaps are counted.
    let processed: Vec<Vec<u8>> = sequences.iter().map(|&seq| {
        let first = seq.iter().position(|&b| b.is_ascii_alphabetic()).unwrap_or(seq.len());
        let last  = seq.iter().rposition(|&b| b.is_ascii_alphabetic()).unwrap_or(0);
        let mut out = seq.to_vec();
        for b in &mut out[..first] { if *b == b'-' { *b = b' '; } }
        if first < seq.len() {
            for b in &mut out[last + 1..] { if *b == b'-' { *b = b' '; } }
        }
        out
    }).collect();
    let seqs: Vec<&[u8]> = processed.iter().map(|v| v.as_slice()).collect();

    // ── Build position-frequency profile ─────────────────────────────────────
    // Spaces (flanking padding) are excluded.  Gap characters within the
    // covered region ARE counted — they can win at gap-dominant columns.
    let mut profile = vec![[0u32; ALPHA_LEN]; width];
    for &seq in &seqs {
        for (col, &b) in seq.iter().enumerate() {
            if b == b' ' { continue; }
            if let Some(idx) = alpha_idx(b.to_ascii_uppercase()) {
                profile[col][idx] += 1;
            }
        }
    }

    // ── First pass: per-column best IUPAC letter ──────────────────────────────
    let mut consensus: Vec<u8> = Vec::with_capacity(width);
    for col in &profile {
        let mut max_score = i64::MIN;
        let mut best_idx = 10usize; // N
        let mut n_score = i64::MIN;

        for cand in 0..ALPHA_LEN { // '-' (gap_idx) is a valid candidate
            let mut score: i64 = 0;
            for obs in 0..ALPHA_LEN {
                let cnt = col[obs] as i64;
                if cnt > 0 {
                    score += cnt * MATRIX[cand][obs] as i64;
                }
            }
            if cand == 10 { n_score = score; }
            if score > max_score {
                max_score = score;
                best_idx = cand;
            }
        }
        // Prefer N when it ties any other candidate (Perl behaviour).
        if best_idx != 10 && n_score == max_score {
            best_idx = 10;
        }
        consensus.push(alpha_byte(best_idx));
    }

    // ── Second pass: CpG dinucleotide correction ─────────────────────────────
    // For each adjacent pair of non-gap consensus positions (i, k), compare:
    //   dnScore  — how well the current consensus dinucleotide fits the data
    //   cgScore  — how well a CG dinucleotide (with CpG deamination model) fits
    // If cgScore > dnScore overwrite both positions with C and G.
    //
    // Only space characters (flanking padding) are skipped; gap characters
    // within the covered region are included, matching Perl's behaviour.
    let c_idx = alpha_idx(b'C').unwrap();
    let g_idx = alpha_idx(b'G').unwrap();

    let mut i = 0usize;
    'outer: while i + 1 < consensus.len() {
        if consensus[i] == b'-' { i += 1; continue; }

        let mut k = i + 1;
        loop {
            if k >= consensus.len() { break 'outer; }
            if consensus[k] != b'-' { break; }
            k += 1;
        }

        let cl_idx = match alpha_idx(consensus[i]) { Some(x) => x, None => { i += 1; continue; } };
        let cr_idx = match alpha_idx(consensus[k]) { Some(x) => x, None => { i += 1; continue; } };

        let mut dn_score: i64 = 0;
        let mut cg_score: i64 = 0;

        for &seq in &seqs {
            if i >= seq.len() { continue; }
            let hl_raw = seq[i];
            if hl_raw == b' ' { continue; }
            let hr_raw = if k < seq.len() { seq[k] } else { b' ' };
            if hr_raw == b' ' { continue; }

            let hl = hl_raw.to_ascii_uppercase();
            let hr = hr_raw.to_ascii_uppercase();

            dn_score += MATRIX[cl_idx][alpha_idx(hl).unwrap_or(10)] as i64;
            dn_score += MATRIX[cr_idx][alpha_idx(hr).unwrap_or(10)] as i64;

            match (hl, hr) {
                (b'C', b'A') | (b'T', b'G') => { cg_score += params.cg_param as i64; }
                (b'T', b'A') => { cg_score += params.ta_param as i64; }
                (b'T', b'C') | (b'T', b'T') => {
                    cg_score += params.cg_trans_param as i64
                        + MATRIX[g_idx][alpha_idx(hr).unwrap_or(10)] as i64;
                }
                (b'A', b'A') | (b'G', b'A') => {
                    cg_score += params.cg_trans_param as i64
                        + MATRIX[c_idx][alpha_idx(hl).unwrap_or(10)] as i64;
                }
                _ => {
                    cg_score += MATRIX[c_idx][alpha_idx(hl).unwrap_or(10)] as i64;
                    cg_score += MATRIX[g_idx][alpha_idx(hr).unwrap_or(10)] as i64;
                }
            }
        }

        if cg_score > dn_score {
            consensus[i] = b'C';
            consensus[k] = b'G';
        }

        i += 1;
    }

    consensus
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conserved_cg_stays() {
        let seqs: Vec<&[u8]> = vec![b"CG", b"CG", b"CG"];
        let params = ConsensusParams::default();
        let cons = build_consensus_from_sequences(&seqs, &params);
        assert_eq!(&cons, b"CG");
    }

    /// When TG and CA dominate, the deamination model should restore CG.
    #[test]
    fn cpg_correction_tg_ca() {
        let seqs: Vec<&[u8]> = vec![
            b"TG", b"TG", b"TG", b"CA", b"CA", b"CA", b"CG", b"CG",
        ];
        let params = ConsensusParams::default();
        let cons = build_consensus_from_sequences(&seqs, &params);
        assert_eq!(&cons, b"CG", "CpG correction should restore CG");
    }

    #[test]
    fn no_false_cpg_correction() {
        let seqs: Vec<&[u8]> = vec![b"AC", b"AC", b"AC", b"AC"];
        let params = ConsensusParams::default();
        let cons = build_consensus_from_sequences(&seqs, &params);
        assert_ne!(&cons, b"CG");
    }

    #[test]
    fn gaps_skipped_in_cpg_correction() {
        // C-G with an interior gap — should still be seen as a CG dinucleotide.
        let seqs: Vec<&[u8]> = vec![b"T-G", b"T-G", b"T-G", b"C-A", b"C-A", b"C-A"];
        let params = ConsensusParams::default();
        let cons = build_consensus_from_sequences(&seqs, &params);
        // Column 1 is all '-' so consensus[1] should be '-'.
        assert_eq!(cons[1], b'-');
        // Positions 0 and 2 should be corrected to CG.
        assert_eq!(cons[0], b'C');
        assert_eq!(cons[2], b'G');
    }
}
