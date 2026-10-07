//! Does a family hold more than one subfamily, and did RepBase describe one of
//! them rather than their average?
//!
//! Score each copy by the share of its covered disputed sites carrying the
//! library's base.  A family built from one homogeneous set puts every copy in
//! one hump.  A family holding two subfamilies, one of which RepBase described,
//! puts copies at both ends with few between: the copies carrying the library's
//! bases, and the copies carrying the other version of them.
//!
//! Splitting the copies is only half of it.  A consensus called from each side
//! and compared to the library says which of the two the library is.  If the
//! library matches one side closely and the other poorly, it is a subfamily
//! consensus and the other subfamily has no entry of its own.  If it sits
//! between them, it is an average and nothing was missed.
//!
//! Two limits worth carrying: the groups come from a threshold on a continuous
//! score, so a family near the cut can fall either way, which is why the raw
//! counts go out alongside any call.  And every disputed site comes from the
//! population consensus disagreeing with the library, so a subfamily split that
//! the population consensus happens to follow leaves no trace here at all.

use crate::family::Instance;

pub struct Groups {
    /// Copies covering enough disputed sites to be scored.
    pub scored: usize,
    /// Copies carrying the library's base at 80% or more of the sites they
    /// cover, most library-like first.
    pub a: Vec<usize>,
    /// Copies carrying it at 20% or less, least library-like first.
    pub b: Vec<usize>,
    pub between: usize,
    /// Library positions where group A almost always carries the library base
    /// and group B almost never does.
    pub diagnostic: Vec<usize>,
}

impl Groups {
    /// Share of scored copies that sit between the two groups.  A low share is
    /// what a two-subfamily family looks like.
    pub fn between_share(&self) -> f64 {
        self.between as f64 / self.scored.max(1) as f64
    }

    /// First and last diagnostic position, and what share of the consensus
    /// they span.  A split confined to one end reads differently from one
    /// spread over the whole element.
    pub fn diagnostic_span(&self, lib_len: usize) -> Option<(usize, usize, f64)> {
        let lo = *self.diagnostic.first()?;
        let hi = *self.diagnostic.last()?;
        Some((lo + 1, hi + 1, (hi - lo + 1) as f64 / lib_len.max(1) as f64))
    }
}

/// A copy's base at each disputed position, 0 where it does not reach.
fn rows(inst: &[&Instance], sites: &[(usize, u8, u8)], lib_len: usize) -> Vec<Vec<u8>> {
    let mut slot = vec![u16::MAX; lib_len];
    for (j, (pos, _, _)) in sites.iter().enumerate() {
        if *pos < lib_len {
            slot[*pos] = j as u16;
        }
    }
    inst.iter()
        .map(|x| {
            let mut row = vec![0u8; sites.len()];
            let mut pos = x.rec.cons_start as usize;
            for (r, q) in x.gapped_reference.iter().zip(&x.gapped_query) {
                if !aln_core::seq::is_gap(*r) && !aln_core::seq::is_pad(*r) {
                    if pos >= 1 && pos - 1 < lib_len {
                        let s = slot[pos - 1];
                        if s != u16::MAX && q.is_ascii_alphabetic() {
                            row[s as usize] = q.to_ascii_uppercase();
                        }
                    }
                    pos += 1;
                }
            }
            row
        })
        .collect()
}

/// Split a family's copies on how much of the library's version they carry.
///
/// `min_frac` is the share of a family's disputed sites a copy must reach
/// before it counts.  An absolute floor does not work: a copy covering five of
/// eighty-four sites gets its score from almost no evidence and lands between
/// the two groups, which buries the very structure the test looks for.  `cap`
/// bounds each group so a satellite does not make the consensus calls
/// expensive.
pub fn groups(
    inst: &[&Instance],
    disputed: &[(usize, u8, u8)],
    lib_len: usize,
    min_frac: f64,
    cap: usize,
) -> Option<Groups> {
    if disputed.len() < 10 {
        return None;
    }
    let m = rows(inst, disputed, lib_len);
    let mut scored: Vec<(f64, usize)> = Vec::new();
    for (i, row) in m.iter().enumerate() {
        let (mut hit, mut seen) = (0usize, 0usize);
        for (j, b) in row.iter().enumerate() {
            if *b != 0 {
                seen += 1;
                hit += usize::from(*b == disputed[j].1);
            }
        }
        if seen as f64 >= (min_frac * disputed.len() as f64).ceil().max(5.0) {
            scored.push((hit as f64 / seen as f64, i));
        }
    }
    if scored.len() < 20 {
        return None;
    }
    let mut a: Vec<(f64, usize)> = scored.iter().copied().filter(|(f, _)| *f >= 0.8).collect();
    let mut b: Vec<(f64, usize)> = scored.iter().copied().filter(|(f, _)| *f <= 0.2).collect();
    let between = scored.len() - a.len() - b.len();
    // Each group's consensus wants the members that span the element, not the
    // ones that score most extremely.  The extremes are often short fragments
    // that disagree everywhere they touch, and a consensus called from those
    // covers a fifth of the element and reads as a diverged subfamily when it
    // is only a thin one.
    let by_span = |v: &mut Vec<(f64, usize)>| {
        v.sort_by_key(|(_, i)| std::cmp::Reverse(inst[*i].span()));
    };
    by_span(&mut a);
    by_span(&mut b);

    // A site is diagnostic when group A carries the library's base there and
    // group B does not, which is what makes the two groups distinguishable.
    let mut diagnostic = Vec::new();
    for (j, (pos, lib, _)) in disputed.iter().enumerate() {
        let side = |g: &[(f64, usize)]| {
            let v: Vec<u8> = g.iter().map(|(_, i)| m[*i][j]).filter(|c| *c != 0).collect();
            (v.iter().filter(|c| *c == lib).count(), v.len())
        };
        let (ha, na) = side(&a);
        let (hb, nb) = side(&b);
        if na >= 5 && nb >= 5 && ha * 10 >= na * 9 && hb * 10 <= nb {
            diagnostic.push(*pos);
        }
    }

    Some(Groups {
        scored: scored.len(),
        a: a.into_iter().take(cap).map(|(_, i)| i).collect(),
        b: b.into_iter().take(cap).map(|(_, i)| i).collect(),
        between,
        diagnostic,
    })
}

/// What the numbers say a curator is looking at.
///
/// `id_a` and `id_b` are each group's consensus against the library.  The call
/// is a summary of the columns beside it, not a substitute for reading them.
pub fn call(g: &Groups, id_a: f64, id_b: f64) -> &'static str {
    let both_populated = g.a.len() >= 5 && g.b.len() >= 5;
    let bimodal = both_populated && g.between_share() <= 0.15;
    if !bimodal {
        if g.a.len() < 5 && g.b.len() >= 5 {
            // Almost every copy carries the other version: the library's bases
            // come from a handful of insertions, or from none in this genome.
            return "library-thin";
        }
        return "single";
    }
    if id_a - id_b >= 2.0 {
        "split"
    } else {
        "mixed"
    }
}
