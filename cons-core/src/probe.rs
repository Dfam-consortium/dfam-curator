//! Cheap pre-flight check for internally repetitive input.
//!
//! An all-vs-all search is quadratic in the number of instances and, for
//! tandemly repetitive sequence, explodes in the number of HSPs per pair: 60
//! GC-rich 20 kb instances produced 8.7M HSPs, 250 s and 25 GB of peak memory,
//! and a consensus with a nonsensical 94% Kimura divergence. The work is wasted
//! — the answer is junk either way — so it is worth a fraction of a second to
//! find out first.
//!
//! The test is self-alignment. A sequence that is not internally repetitive
//! aligns to itself as one diagonal HSP plus a little spurious noise; a
//! tandemly repetitive one matches itself at every period, so the count
//! explodes. Measured on a 2 kb window:
//!
//! | input                        | self-HSPs |
//! |------------------------------|-----------|
//! | random 20 kb, 10% divergence | 1         |
//! | AluY simulation              | ~2        |
//! | GC-rich pathological family  | 133       |
//! | pure `(GGGAGG)n`             | ~57       |
//!
//! Two details keep it cheap. The probe cost tracks the HSP count it is
//! looking for — 6.07 s to self-align a full 20 kb pathological instance
//! against 0.14 s for a clean one — so it is *truncated* to a window: at 2 kb
//! the same instance costs 0.04 s, 150x less, and still separates by 100x.
//! And because every probe is the same length, one threshold serves all
//! inputs rather than scaling with sequence length.
//!
//! The window is the middle of the sequence, chosen so the probe sees the bulk
//! of the element rather than whichever end the assembly happened to start at.
//! It is deterministic, not sampled, so a rerun on the same input gives the same
//! verdict. The cost is that a tandem array confined to one end is missed — the
//! probe is a guard against families that are repetitive throughout, which is
//! what makes the search explode, not a general repeat detector.

use crate::AlignmentSource;
use aln_core::Sequence;
use aln_engine::Result;

/// How much to probe, and what counts as repetitive.
#[derive(Debug, Clone)]
pub struct ProbeConfig {
    /// Bases probed per instance. Sequences shorter than this are probed whole.
    pub window: usize,
    /// How many instances to probe, longest first.
    ///
    /// One instance catches a uniformly repetitive family, which is the common
    /// case, but would miss a family where only a few instances are
    /// contaminated — the realistic `RepeatModeler` failure. At 0.04 s each,
    /// sampling costs nothing.
    pub sample: usize,
    /// Self-HSPs above which an instance is called repetitive.
    ///
    /// Calibrated on 1,074 real families with DUST off and a 2 kb window. The
    /// worst instance of any curated hs1 family scored 20 and of any simulated
    /// family 27; the pathological cases score 57 (`(GGGAGG)n`) and 142
    /// (a GC-rich family that took 250 s and 25 GB). 40 sits in that gap.
    pub max_hsps: usize,
}

impl Default for ProbeConfig {
    fn default() -> Self {
        ProbeConfig { window: 2000, sample: 10, max_hsps: 40 }
    }
}

/// One instance's probe result.
#[derive(Debug, Clone)]
pub struct ProbeHit {
    pub index: usize,
    pub name: String,
    pub len: usize,
    pub hsps: usize,
}

#[derive(Debug, Clone, Default)]
pub struct ProbeReport {
    /// Every instance probed, worst first.
    pub probed: Vec<ProbeHit>,
    /// Those over [`ProbeConfig::max_hsps`].
    pub flagged: Vec<ProbeHit>,
}

impl ProbeReport {
    pub fn worst(&self) -> usize {
        self.probed.first().map(|p| p.hsps).unwrap_or(0)
    }
    pub fn median(&self) -> usize {
        if self.probed.is_empty() {
            return 0;
        }
        // `probed` is sorted worst-first, so the midpoint is the median.
        self.probed[self.probed.len() / 2].hsps
    }
}

/// Self-HSP count for a single sequence, using the same window and rules as
/// [`probe_repetitive`].
///
/// Used on an *extended consensus*, where a satellite becomes visible that the
/// input instances were too short to reveal. Same caller obligations: a raw
/// engine with DUST off.
pub fn self_hsp_count<A: AlignmentSource>(
    aligner: &A,
    seq: &[u8],
    cfg: &ProbeConfig,
) -> Result<usize> {
    if seq.is_empty() {
        return Ok(0);
    }
    let start = seq.len().saturating_sub(cfg.window) / 2;
    let end = (start + cfg.window).min(seq.len());
    let win = Sequence::new("consensus", seq[start..end].to_vec());
    let one = vec![win.clone()];
    Ok(aligner.against_reference(&win, &one, None)?.len())
}

/// Self-align a sample of instances and report how repetitive they look.
///
/// `aligner` should be the **raw** engine, not a filtered stack: the count of
/// raw HSPs is the signal, and a cull would erase it.
///
/// It must also have **DUST off**, whatever the shipping default is. DUST masks
/// low-complexity query sequence, so a pure tandem repeat self-aligns to
/// nothing and the probe reports 0 — the one input it most needs to catch. The
/// `(GGGAGG)n` synthetic scores 0 with DUST on and 133 with it off.
pub fn probe_repetitive<A: AlignmentSource>(
    aligner: &A,
    seqs: &[Sequence],
    cfg: &ProbeConfig,
) -> Result<ProbeReport> {
    let mut order: Vec<usize> = (0..seqs.len()).filter(|&i| !seqs[i].is_empty()).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(seqs[i].len()));
    order.truncate(cfg.sample.max(1));

    let mut probed = Vec::with_capacity(order.len());
    for i in order {
        let s = &seqs[i];
        let start = s.len().saturating_sub(cfg.window) / 2;
        let end = (start + cfg.window).min(s.len());
        let win = Sequence::new(s.name.clone(), s.seq[start..end].to_vec());
        // Self-alignment: the window is both the reference and the only
        // instance, with no `skip`, so it aligns against itself.
        let one = vec![win.clone()];
        let hits = aligner.against_reference(&win, &one, None)?;
        probed.push(ProbeHit {
            index: i,
            name: s.name.clone(),
            len: s.len(),
            hsps: hits.len(),
        });
    }
    probed.sort_by_key(|p| std::cmp::Reverse(p.hsps));
    let flagged = probed.iter().filter(|p| p.hsps > cfg.max_hsps).cloned().collect();
    Ok(ProbeReport { probed, flagged })
}
