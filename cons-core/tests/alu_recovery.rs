//! End-to-end check on real sequence: can `autocons` reconstruct a known
//! consensus from a diverged family?
//!
//! The AluY consensus below is the 311 bp sequence from `rmblast-port/alu.fa`.
//! The test diverges it into 20 copies at ~15% substitution and ~3% indel, buries
//! each in random flanking sequence so hits are embedded rather than
//! edge-aligned, and asks `autocons` to recover the original.
//!
//! This exercises the whole stack in one shot: alignment (parasail striped SIMD),
//! MSA assembly under `InsertionPolicy::GrowReference`, the Dfam consensus
//! caller with CpG restoration, and refinement to convergence.

use aln_core::{Sequence, SubstMatrix};
use aln_engine::{AlignMode, AlignParams};
use cons_core::FastAligner;
use cons_core::{run, Pairwise, Params};

/// AluY, 311 bp.
const ALUY: &str = concat!(
    "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGA",
    "TCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAA",
    "AAATACAAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGG",
    "CTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGC",
    "CACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAA",
    "AAAAAAAAAAA",
);

const MATRIX: &str = "\
FREQS A 0.325 C 0.175 G 0.175 T 0.325
  A   R   G   C   Y   T   K   M   S   W   N   X
  8   0 -10 -18 -19 -21 -15  -4 -14  -6  -1 -30
  3   3  12 -17 -18 -19  -9  -8  -8  -9  -1 -30
 -7   2  12 -16 -16 -17  -2 -11  -1 -12  -1 -30
-17 -16 -16  12   2  -7 -11  -2  -1 -12  -1 -30
-19 -18 -17  12   0   3  -8  -9  -8  -9  -1 -30
-21 -19 -18 -10   0   8  -4 -15 -14  -6  -1 -30
-14  -8  -2 -13  -8  -4  -3 -13  -8  -9  -1 -30
 -4  -8 -13  -2  -8 -14 -13  -3  -8  -9  -1 -30
-12  -7  -1  -1  -7 -12  -7  -7  -1 -12  -1 -30
 -6 -10 -14 -14 -10  -6 -10 -10 -14  -6  -1 -30
 -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
-30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
";

struct Rng(u64);
impl Rng {
    fn new(s: u64) -> Self {
        Rng(s.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1)
    }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
    /// Percent roll in [0, 100).
    fn pct(&mut self) -> usize {
        self.below(100)
    }
}

const BASES: &[u8; 4] = b"ACGT";

fn diverge(rng: &mut Rng, src: &[u8], subst_pct: usize, indel_pct: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(src.len());
    for &b in src {
        let roll = rng.pct();
        if roll < indel_pct {
            if rng.below(2) == 0 {
                continue; // deletion
            }
            out.push(BASES[rng.below(4)]); // insertion, then keep the base
            out.push(b);
        } else if roll < indel_pct + subst_pct {
            out.push(BASES[rng.below(4)]);
        } else {
            out.push(b);
        }
    }
    out
}

fn random_seq(rng: &mut Rng, len: usize) -> Vec<u8> {
    (0..len).map(|_| BASES[rng.below(4)]).collect()
}

/// Build the diverged family, each copy buried in random flanks.
fn family(seed: u64, n: usize, subst_pct: usize, indel_pct: usize) -> Vec<Sequence> {
    let mut rng = Rng::new(seed);
    (0..n)
        .map(|i| {
            let core = diverge(&mut rng, ALUY.as_bytes(), subst_pct, indel_pct);
            let lf_len = rng.below(41);
            let mut s = random_seq(&mut rng, lf_len);
            s.extend_from_slice(&core);
            let rf_len = rng.below(41);
            s.extend_from_slice(&random_seq(&mut rng, rf_len));
            Sequence::new(format!("copy{i}"), s)
        })
        .collect()
}

fn aligner() -> FastAligner {
    let p = AlignParams {
        mode: AlignMode::Local,
        gap_open: 25,
        gap_extend: 5,
        min_score: 200,
        traceback: true,
        bandwidth: None,
    };
    FastAligner::new(SubstMatrix::parse(MATRIX).unwrap(), p).unwrap()
}

/// Smith-Waterman (+1 / -1 / -2, linear gaps) — enough to score the recovery
/// without dragging the reference aligner in as a dependency.
fn local_identity(a: &[u8], b: &[u8]) -> (usize, usize) {
    let (n, m) = (a.len(), b.len());
    let mut prev = vec![0i32; m + 1];
    let mut ptr = vec![vec![0u8; m + 1]; n + 1];
    let mut best = (0i32, 0usize, 0usize);

    for i in 1..=n {
        let mut cur = vec![0i32; m + 1];
        for j in 1..=m {
            let d = prev[j - 1] + if a[i - 1] == b[j - 1] { 1 } else { -1 };
            let up = prev[j] - 2;
            let left = cur[j - 1] - 2;
            let v = 0.max(d).max(up).max(left);
            cur[j] = v;
            ptr[i][j] = if v == 0 {
                0
            } else if v == d {
                1
            } else if v == up {
                2
            } else {
                3
            };
            if v > best.0 {
                best = (v, i, j);
            }
        }
        prev = cur;
    }

    let (_, mut i, mut j) = best;
    let (mut cols, mut ident) = (0usize, 0usize);
    while i > 0 && j > 0 && ptr[i][j] != 0 {
        match ptr[i][j] {
            1 => {
                if a[i - 1] == b[j - 1] {
                    ident += 1;
                }
                i -= 1;
                j -= 1;
            }
            2 => i -= 1,
            _ => j -= 1,
        }
        cols += 1;
    }
    (ident, cols)
}

#[test]
fn recovers_aluy_from_a_diverged_family() {
    let seqs = family(20_260_807, 20, 15, 3);
    let out = run(&Pairwise::new(aligner()), &seqs, &Params::default()).unwrap();
    assert_eq!(out.len(), 1, "expected a single consensus");

    let cons = &out[0].consensus;
    let (ident, cols) = local_identity(ALUY.as_bytes(), cons);
    let identity = ident as f64 / cols as f64;
    let coverage = cols as f64 / ALUY.len() as f64;

    assert!(
        identity > 0.98,
        "recovered consensus is only {:.1}% identical to AluY over {cols} columns",
        identity * 100.0
    );
    assert!(
        coverage > 0.95,
        "recovered consensus covers only {:.1}% of AluY",
        coverage * 100.0
    );
    assert!(out[0].converged, "refinement should converge");
    assert_eq!(
        out[0].msa.num_instances(),
        seqs.len(),
        "every copy should be placed in the final alignment"
    );
}

/// Divergence high enough to be interesting but still within reach.
#[test]
fn still_recovers_aluy_at_higher_divergence() {
    let seqs = family(99, 30, 22, 4);
    let out = run(&Pairwise::new(aligner()), &seqs, &Params::default()).unwrap();
    let cons = &out[0].consensus;
    let (ident, cols) = local_identity(ALUY.as_bytes(), cons);
    let identity = ident as f64 / cols as f64;
    assert!(
        identity > 0.95,
        "at 22% divergence the consensus is {:.1}% identical over {cols} columns",
        identity * 100.0
    );
}

/// The consensus must be a better reconstruction than any single input copy —
/// otherwise the exercise achieved nothing.
#[test]
fn the_consensus_beats_every_input_copy() {
    let seqs = family(7, 20, 15, 3);
    let out = run(&Pairwise::new(aligner()), &seqs, &Params::default()).unwrap();

    let score = |s: &[u8]| {
        let (ident, cols) = local_identity(ALUY.as_bytes(), s);
        ident as f64 / cols.max(1) as f64
    };
    let cons_id = score(&out[0].consensus);
    let best_copy = seqs.iter().map(|s| score(&s.seq)).fold(0.0, f64::max);

    assert!(
        cons_id > best_copy,
        "consensus {cons_id:.3} should beat the best single copy {best_copy:.3}"
    );
}
