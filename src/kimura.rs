/// Kimura two-parameter divergence statistics between two sequences.
#[derive(Debug, Default, Clone)]
pub struct KimuraStats {
    /// Classic Kimura 2-parameter divergence (%).
    pub kimura: f64,

    /// CpG-adjusted Kimura divergence (%).
    pub kimura_adjusted: f64,

    pub transitions: u32,
    /// CpG-adjusted transition count (fractional, matching the Perl algorithm).
    pub transitions_adjusted: f64,
    pub transversions: u32,

    /// Number of well-characterised (non-N/gap) compared positions.
    pub well_characterised: u32,

    /// Number of CpG sites detected.
    pub cpg_sites: u32,

    /// Number of positions with high divergence.
    pub high_div: u32,
}

/// Classify a base substitution as transition (Ti) or transversion (Tv).
///
/// Returns `None` for identical bases, gaps, or ambiguity codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubstClass {
    Transition,
    Transversion,
}

pub fn classify_subst(ref_base: u8, obs_base: u8) -> Option<SubstClass> {
    let r = ref_base.to_ascii_uppercase();
    let o = obs_base.to_ascii_uppercase();
    if r == o { return None; }
    // Only score unambiguous bases.
    if !matches!(r, b'A' | b'C' | b'G' | b'T') { return None; }
    if !matches!(o, b'A' | b'C' | b'G' | b'T') { return None; }

    // Purines: A, G  |  Pyrimidines: C, T
    let r_pur = matches!(r, b'A' | b'G');
    let o_pur = matches!(o, b'A' | b'G');

    if r_pur == o_pur {
        Some(SubstClass::Transition)
    } else {
        Some(SubstClass::Transversion)
    }
}

/// Compute Kimura 2-parameter divergence for a single gapped sequence pair
/// (consensus vs. one instance), matching the Perl MultAln::calculateKimuraDivergence
/// algorithm exactly.
///
/// CpG adjustment uses a look-back state machine:
/// - Each consecutive non-gap consensus C→G pair is a CpG site.
/// - A transition at one of the two positions → 0.1 adjusted transitions.
/// - Transitions at both positions → 1.0 adjusted transitions.
/// - Transitions outside CpG context contribute 1.0 each (unchanged).
///
/// `consensus` and `seq` must have the same length (alignment width).
pub fn kimura_pair(consensus: &[u8], seq: &[u8]) -> KimuraStats {
    let mut ti = 0u32;
    let mut tv = 0u32;
    let mut well = 0u32;
    let mut cpg = 0u32;

    let len = consensus.len().min(seq.len());

    // Perl state: previous non-gap consensus base, and pending transition credit.
    let mut prev_cons: u8 = 0;
    let mut prev_trans: f64 = 0.0;
    let mut ti_adj: f64 = 0.0;

    for i in 0..len {
        let r = consensus[i].to_ascii_uppercase();
        let o = seq[i].to_ascii_uppercase();

        // Skip consensus gaps (insertions in the instance relative to the consensus).
        // Perl: `next if $sBases[$i] eq "-"`.  prev_cons is NOT updated.
        if r == b'-' {
            continue;
        }

        // Skip positions where the instance is absent (flanking padding).
        // Outside the aligned region there is nothing to compare.
        if o == b' ' {
            continue;
        }

        // Well-characterised: both sides are unambiguous ACGT (not gap, not IUPAC code).
        if matches!(r, b'A' | b'C' | b'G' | b'T') && matches!(o, b'A' | b'C' | b'G' | b'T') {
            well += 1;
        }

        // CpG detection: look back at the previous non-gap consensus base.
        // Perl: `if ($pSBase eq "C" && $sBases[$i] eq "G")`.
        let in_cpg = prev_cons == b'C' && r == b'G';

        if in_cpg {
            cpg += 1;
            // Determine mutation type at the G position.
            let mt = classify_subst(r, o);
            if let Some(SubstClass::Transition) = mt {
                prev_trans += 1.0; // includes C transition if any
                ti += 1;
            } else if let Some(SubstClass::Transversion) = mt {
                tv += 1;
            }
            // Apply fractional CpG weighting.
            // 2 transitions (C and G both mutated) → count as 1.
            // 1 transition (C or G mutated) → count as 0.1.
            // 0 transitions → count as 0 (no change).
            if prev_trans >= 2.0 {
                prev_trans = 1.0;
            } else if prev_trans >= 1.0 {
                prev_trans = 0.1;
            }
        } else {
            // Non-CpG position: flush the pending credit from the previous position.
            ti_adj += prev_trans;
            prev_trans = 0.0;

            let mt = classify_subst(r, o);
            match mt {
                Some(SubstClass::Transition) => {
                    ti += 1;
                    prev_trans = 1.0; // hold — next position might form a CpG pair
                }
                Some(SubstClass::Transversion) => {
                    tv += 1;
                }
                None => {}
            }
        }

        // Update the look-back state.  Perl: `$pSBase = $sBases[$i]`.
        prev_cons = r;
    }
    // Flush any pending credit at end-of-sequence.
    ti_adj += prev_trans;

    let kimura = k2p(ti as f64, tv as f64, well);
    let kimura_adjusted = k2p(ti_adj, tv as f64, well);

    KimuraStats {
        kimura,
        kimura_adjusted,
        transitions: ti,
        transitions_adjusted: ti_adj,
        transversions: tv,
        well_characterised: well,
        cpg_sites: cpg,
        high_div: 0,
    }
}

/// Kimura 2-parameter formula.
///
/// Accepts float transition count to support fractional CpG adjustment.
/// Returns divergence in percent (0–100), or NaN if the formula breaks down.
fn k2p(ti: f64, tv: f64, n: u32) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let p = ti / n as f64; // transition fraction
    let q = tv / n as f64; // transversion fraction

    let a = 1.0 - 2.0 * p - q;
    let b = 1.0 - 2.0 * q;

    if a <= 0.0 || b <= 0.0 {
        return f64::NAN;
    }

    // Guard against IEEE-754 negative zero when p=q=0.
    (-0.5 * a.ln() - 0.25 * b.ln()).max(0.0) * 100.0
}

/// Compute per-sequence Kimura divergence for a whole alignment and return the
/// mean over all instances.
pub fn mean_kimura(
    consensus: &[u8],
    instances: &[&[u8]],
    cpg_adjusted: bool,
) -> f64 {
    if instances.is_empty() {
        return 0.0;
    }
    let sum: f64 = instances
        .iter()
        .map(|seq| {
            let s = kimura_pair(consensus, seq);
            if cpg_adjusted { s.kimura_adjusted } else { s.kimura }
        })
        .filter(|v| v.is_finite())
        .sum();
    sum / instances.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_seqs_zero_divergence() {
        let cons = b"ACGT";
        let seq  = b"ACGT";
        let s = kimura_pair(cons, seq);
        assert_eq!(s.transitions, 0);
        assert_eq!(s.transversions, 0);
        assert!((s.kimura - 0.0).abs() < 1e-9);
    }

    #[test]
    fn some_transitions() {
        // One A→G transition out of four sites.
        let cons = b"AAAA";
        let seq  = b"GAAA";
        let s = kimura_pair(cons, seq);
        assert_eq!(s.transitions, 1);
        assert_eq!(s.transversions, 0);
        assert!(s.kimura > 0.0 && s.kimura.is_finite());
    }

    #[test]
    fn all_transitions_is_undefined() {
        // 100% transition rate saturates K2P (returns NaN).
        let cons = b"AAAA";
        let seq  = b"GGGG";
        let s = kimura_pair(cons, seq);
        assert_eq!(s.transitions, 4);
        assert!(s.kimura.is_nan());
    }

    #[test]
    fn classify_transition_transversion() {
        assert_eq!(classify_subst(b'A', b'G'), Some(SubstClass::Transition));
        assert_eq!(classify_subst(b'C', b'T'), Some(SubstClass::Transition));
        assert_eq!(classify_subst(b'A', b'C'), Some(SubstClass::Transversion));
        assert_eq!(classify_subst(b'G', b'T'), Some(SubstClass::Transversion));
        assert_eq!(classify_subst(b'A', b'A'), None);
    }

    #[test]
    fn cpg_both_transitions_count_as_one() {
        // Consensus: CG, instance: TA (C→T transition, G→A transition).
        // Both CpG positions mutated → ti_adj = 1.0 (not 2.0).
        let cons = b"ACGX";
        let seq  = b"ATAX";
        let s = kimura_pair(cons, seq);
        assert_eq!(s.transitions, 2); // raw transitions = 2
        assert!((s.transitions_adjusted - 1.0).abs() < 1e-9); // adjusted = 1
        assert_eq!(s.cpg_sites, 1);
    }

    #[test]
    fn cpg_single_transition_counts_as_tenth() {
        // Consensus: CG, instance: TG (C→T transition only).
        // One CpG position mutated → ti_adj = 0.1.
        let cons = b"ACGX";
        let seq  = b"ATGX";
        let s = kimura_pair(cons, seq);
        assert_eq!(s.transitions, 1); // raw = 1
        assert!((s.transitions_adjusted - 0.1).abs() < 1e-9); // adjusted = 0.1
        assert_eq!(s.cpg_sites, 1);
    }

    #[test]
    fn non_cpg_transition_unaffected() {
        // Consensus: AG (no CpG), instance: GG (A→G transition).
        let cons = b"AG";
        let seq  = b"GG";
        let s = kimura_pair(cons, seq);
        assert_eq!(s.transitions, 1);
        assert!((s.transitions_adjusted - 1.0).abs() < 1e-9);
        assert_eq!(s.cpg_sites, 0);
    }
}
