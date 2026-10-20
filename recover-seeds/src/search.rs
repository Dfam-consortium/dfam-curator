//! Local search for the copies whose consensus reproduces the library's.
//!
//! The closest ranking, and the samples built on it, land near the best
//! subset of a family's closest copies but rarely on it.  An integer program
//! over the closest 30 copies, scoring each library position by whether the
//! consensus caller's column rule calls the library's base there, found sets
//! whose real consensi reached 99% edit identity in 189 of 300 RepBase 31.08
//! families that failed it (69 under best-of).  This searches the same
//! objective by hill climbing, fast enough to run on every family.
//!
//! The objective is the column rule alone: each copy's state at each library
//! position (its base, a deletion, or not reaching it), the caller's matrix
//! and its tie order.  Insertions relative to the library, the CpG pass and
//! how a multiple alignment places gaps are left to the real build, which
//! the caller of [`search`] does afterwards.

use crate::family::Instance;

/// One consensus caller's column rule.
pub struct Rule {
    /// Candidate symbols; observed states use the same indices.
    syms: &'static [u8],
    /// score[candidate][observed]
    score: Vec<Vec<i32>>,
    /// For the library's symbol r and a rival c: must r score strictly
    /// higher (true) or is a tie enough (false)?
    strict: Vec<Vec<bool>>,
}

const GIRI_SYMS: &[u8] = b"ARGCYTKMSWNX-";
const GIRI_M: [[i32; 12]; 12] = [
    [9, 0, -8, -15, -16, -17, -13, -3, -11, -4, -5, -7],
    [2, 1, 1, -15, -15, -16, -7, -6, -6, -7, -5, -7],
    [-4, 3, 10, -14, -14, -15, -2, -9, -2, -9, -5, -7],
    [-15, -14, -14, 10, 3, -4, -9, -2, -2, -9, -5, -7],
    [-16, -15, -15, 1, 1, 2, -6, -7, -6, -7, -5, -7],
    [-17, -16, -15, -8, 0, 9, -3, -13, -11, -4, -5, -7],
    [-11, -6, -2, -11, -7, -3, -2, -11, -6, -7, -5, -7],
    [-3, -7, -11, -2, -6, -11, -11, -2, -6, -7, -5, -7],
    [-9, -5, -2, -2, -5, -9, -5, -5, -2, -9, -5, -7],
    [-4, -8, -11, -11, -8, -4, -8, -8, -11, -4, -5, -7],
    [-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -1, -7],
    [-7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7],
];
const DFAM_SYMS: &[u8] = b"ARGCYTKMSWNXZVHDB-";

impl Rule {
    /// GIRI's acons caller (aln-core giri.rs): a gap scores -5 against a
    /// base and +2 against a gap; candidates are scanned from the gap down to
    /// A, a later one taking over only on a strictly higher score, so the gap
    /// wins any tie.
    pub fn giri() -> Rule {
        let n = GIRI_SYMS.len();
        let gap = n - 1;
        let score = (0..n)
            .map(|c| {
                (0..n)
                    .map(|o| {
                        if c == gap || o == gap {
                            if c == o { 2 } else { -5 }
                        } else {
                            GIRI_M[c][o]
                        }
                    })
                    .collect()
            })
            .collect();
        let strict = (0..n).map(|r| (0..n).map(|c| c > r).collect()).collect();
        Rule { syms: GIRI_SYMS, score, strict }
    }

    /// The Dfam caller's column rule without its CpG pass (aln-core
    /// consensus.rs): scanned from A up, an earlier candidate keeps a tie,
    /// and N then takes any tie for the top score.
    pub fn dfam() -> Rule {
        let n = DFAM_SYMS.len();
        let score = (0..n).map(|c| (0..n).map(|o| aln_core::consensus::MATRIX[c][o]).collect()).collect();
        let nidx = DFAM_SYMS.iter().position(|&b| b == b'N').unwrap();
        let strict = (0..n).map(|r| (0..n).map(|c| c < r || c == nidx).collect()).collect();
        Rule { syms: DFAM_SYMS, score, strict }
    }

    fn index(&self, b: u8) -> usize {
        let b = b.to_ascii_uppercase();
        self.syms.iter().position(|&s| s == b).unwrap_or_else(|| self.syms.iter().position(|&s| s == b'N').unwrap())
    }
}

/// Tuning for [`search`].
pub struct Params {
    /// Copies most closely ranked that may be chosen, besides the starts'
    pub pool: usize,
    /// Most copies in a set
    pub max_copies: usize,
    /// After the best column count is reached, copies are added while the
    /// set stays within this many lost library positions per kb of it
    pub margin_per_kb: f64,
    /// Random subsets of the pool to climb from, besides the given starts and
    /// the top 5, 10, 20 and 30 of the ranking (which are added when this is
    /// above 0)
    pub restarts: usize,
    /// When single moves stop helping, try adding two copies, dropping two,
    /// or swapping one out for two in
    pub pairs: bool,
    /// Closest-ranked copies the margin stage may add from; those beyond the
    /// pool take no part in the search's moves
    pub depth_pool: usize,
    /// Count a deaminated CpG call (T at a library CG's C, A at its G) as won
    pub cpg_tolerant: bool,
}

/// A copy's states at the library positions it reaches: (position index,
/// observed symbol index), positions 0-based.
struct Copy {
    lo: usize,
    hi: usize,
    states: Vec<(usize, u8)>,
}

struct State<'a> {
    rule: &'a Rule,
    lib: Vec<u8>,          // library symbol index per position, or u8::MAX where not A/C/G/T
    alt: Vec<u8>,          // a call also counted as won here (CpG deamination product), or u8::MAX
    tot: Vec<Vec<i32>>,    // per position, total score of each candidate
    depth: Vec<u32>,
    won: Vec<bool>,
    wins: i64,
    uncovered: usize,
}

impl<'a> State<'a> {
    fn new(rule: &'a Rule, library: &[u8], cpg_tolerant: bool) -> State<'a> {
        let n = rule.syms.len();
        let lib: Vec<u8> = library
            .iter()
            .map(|&b| if b"ACGTacgt".contains(&b) { rule.index(b) as u8 } else { u8::MAX })
            .collect();
        let len = library.len();
        // At a library CG, a call of T at the C or A at the G is what a
        // consensus of deaminated copies gives before CpG restoration, so it
        // counts as won: the library may have been built with restoration.
        let up = |i: usize| library.get(i).map(|b| b.to_ascii_uppercase());
        let alt: Vec<u8> = (0..len)
            .map(|i| {
                if !cpg_tolerant {
                    u8::MAX
                } else if up(i) == Some(b'C') && up(i + 1) == Some(b'G') {
                    rule.index(b'T') as u8
                } else if up(i) == Some(b'G') && i > 0 && up(i - 1) == Some(b'C') {
                    rule.index(b'A') as u8
                } else {
                    u8::MAX
                }
            })
            .collect();
        State { rule, lib, alt, tot: vec![vec![0; n]; len], depth: vec![0; len], won: vec![false; len], wins: 0, uncovered: len }
    }

    fn wins_at(&self, p: usize) -> bool {
        let r = self.lib[p];
        if r == u8::MAX || self.depth[p] == 0 {
            return false;
        }
        let t = &self.tot[p];
        let wins = |r: usize| (0..t.len()).all(|c| c == r || if self.rule.strict[r][c] { t[r] > t[c] } else { t[r] >= t[c] });
        wins(r as usize) || (self.alt[p] != u8::MAX && wins(self.alt[p] as usize))
    }

    /// Add (sign 1) or remove (sign -1) a copy, keeping wins and uncovered.
    fn apply(&mut self, cp: &Copy, sign: i32) {
        for &(p, o) in &cp.states {
            let o = o as usize;
            for c in 0..self.tot[p].len() {
                self.tot[p][c] += sign * self.rule.score[c][o];
            }
        }
        for p in cp.lo..=cp.hi {
            let before = self.depth[p];
            if sign > 0 { self.depth[p] += 1 } else { self.depth[p] -= 1 }
            if before == 0 && self.depth[p] > 0 {
                self.uncovered -= 1;
            } else if before > 0 && self.depth[p] == 0 {
                self.uncovered += 1;
            }
            let w = self.wins_at(p);
            if w != self.won[p] {
                self.wins += if w { 1 } else { -1 };
                self.won[p] = w;
            }
        }
    }
}

fn copy_states(inst: &Instance, rule: &Rule, len: usize) -> Option<Copy> {
    let mut states = Vec::new();
    let mut pos = inst.rec.cons_start as usize;
    for (r, q) in inst.gapped_reference.iter().zip(&inst.gapped_query) {
        if aln_core::seq::is_gap(*r) || aln_core::seq::is_pad(*r) {
            continue;
        }
        if pos >= 1 && pos <= len {
            let o = if aln_core::seq::is_gap(*q) { rule.syms.len() - 1 } else { rule.index(*q) };
            states.push((pos - 1, o as u8));
        }
        pos += 1;
    }
    let lo = states.first()?.0;
    let hi = states.last()?.0;
    Some(Copy { lo, hi, states })
}

/// Search for the set of copies whose column calls match `library` best.
///
/// `ranked` is the closest ranking; `starts` are the sets to climb from (the
/// first is usually the seed kept so far).  No set may leave more library
/// positions uncovered than `max_uncovered`.  Returns the chosen copy indices
/// and the library positions the column rule loses with them.
pub fn search(
    inst: &[Instance],
    ranked: &[usize],
    starts: &[Vec<usize>],
    library: &[u8],
    rule: &Rule,
    max_uncovered: usize,
    p: &Params,
) -> Option<(Vec<usize>, i64)> {
    let len = library.len();
    let mut pool: Vec<usize> = ranked.iter().copied().take(p.pool).collect();
    for s in starts {
        for &i in s {
            if !pool.contains(&i) {
                pool.push(i);
            }
        }
    }
    // copies past this index are for the margin stage only
    let n_search = pool.len();
    for &i in ranked.iter().take(p.depth_pool) {
        if !pool.contains(&i) {
            pool.push(i);
        }
    }
    let rank_of: std::collections::HashMap<usize, usize> = ranked.iter().enumerate().map(|(k, &i)| (i, k)).collect();
    let copies: Vec<Option<Copy>> = pool.iter().map(|&i| copy_states(&inst[i], rule, len)).collect();
    let judged = library.iter().filter(|b| b"ACGTacgt".contains(b)).count() as i64;
    // Copies added for depth, or on a tie in columns won, must be long: span
    // at least half the library, or the median span of the closest-ranked
    // copies, whichever is less, so a fragmentary family still gets depth
    // from its best copies while fragments enter only by winning columns
    // (or through the coverage fill, before the search).
    let cspan = |k: usize| copies[k].as_ref().map(|c| c.hi - c.lo + 1).unwrap_or(0);
    let mut top: Vec<usize> = ranked
        .iter()
        .take(p.depth_pool.max(p.pool))
        .filter_map(|i| pool.iter().position(|x| x == i))
        .map(cspan)
        .filter(|&s| s > 0)
        .collect();
    top.sort_unstable();
    let min_len = if top.is_empty() { len / 2 } else { (len / 2).min(top[top.len() / 2]) };
    let long = |k: usize| cspan(k) >= min_len;
    let mut starts: Vec<Vec<usize>> = starts.to_vec();
    if p.restarts > 0 {
        for k in [5, 10, 20, 30] {
            starts.push(ranked.iter().copied().take(k.min(p.pool)).collect());
        }
        // a fixed-seed generator, so a run repeats exactly
        let mut x: u64 = 0x9E37_79B9_7F4A_7C15 ^ (len as u64) ^ ((pool.len() as u64) << 32);
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let top: Vec<usize> = ranked.iter().copied().take(p.pool).collect();
        for _ in 0..p.restarts {
            let size = 3 + (next() % 13) as usize;
            let mut set: Vec<usize> = Vec::new();
            for _ in 0..size * 3 {
                if set.len() >= size.min(top.len()) {
                    break;
                }
                let c = top[(next() % top.len().max(1) as u64) as usize];
                if !set.contains(&c) {
                    set.push(c);
                }
            }
            if !set.is_empty() {
                starts.push(set);
            }
        }
    }
    // per start: (columns won after phase 1, final set, final lost)
    let mut best: Option<(i64, Vec<usize>, i64)> = None;

    for start in starts {
        let mut st = State::new(rule, library, p.cpg_tolerant);
        let mut chosen = vec![false; pool.len()];
        // Start from the start set's long copies only; fragments come back
        // only where coverage needs them, closest-ranked first.  A start of
        // nothing but fragments is kept whole.
        let in_start: Vec<usize> = (0..pool.len()).filter(|&k| start.contains(&pool[k]) && copies[k].is_some()).collect();
        let longs: Vec<usize> = in_start.iter().copied().filter(|&k| long(k)).collect();
        for &k in if longs.is_empty() { &in_start } else { &longs } {
            st.apply(copies[k].as_ref().unwrap(), 1);
            chosen[k] = true;
        }
        let n_chosen = |c: &[bool]| c.iter().filter(|x| **x).count();
        let ok = |st: &State| st.uncovered <= max_uncovered;
        if !ok(&st) {
            let mut order: Vec<usize> = (0..pool.len()).filter(|&k| !chosen[k] && copies[k].is_some()).collect();
            order.sort_by_key(|&k| rank_of.get(&pool[k]).copied().unwrap_or(usize::MAX));
            for k in order {
                if ok(&st) {
                    break;
                }
                let before = st.uncovered;
                st.apply(copies[k].as_ref().unwrap(), 1);
                if st.uncovered < before {
                    chosen[k] = true;
                } else {
                    st.apply(copies[k].as_ref().unwrap(), -1);
                }
            }
        }
        // a copy's aligned length on the library, and a set's depth: their sum
        let span = |k: usize| copies[k].as_ref().map(|c| c.hi - c.lo + 1).unwrap_or(0);
        let depth_of = |c: &[bool]| (0..c.len()).filter(|&k| c[k]).map(span).sum::<usize>();
        let frag = |k: usize| usize::from(!long(k));
        let frags_of = |c: &[bool]| (0..c.len()).filter(|&k| c[k]).map(frag).sum::<usize>();
        // phase 1: climb on columns won; a tie goes to the set with fewer
        // fragments, then to the deeper one (more aligned bases), so a
        // fragment that wins nothing is dropped rather than kept
        type Key = (i64, i64, usize);
        let key = |w: i64, f: usize, d: usize| -> Key { (w, -(f as i64), d) };
        loop {
            let (w0, k0, d0, f0) = (st.wins, n_chosen(&chosen), depth_of(&chosen), frags_of(&chosen));
            let k_now = key(w0, f0, d0);
            let mut best_move: Option<(Key, Option<usize>, Option<usize>)> = None;
            let better = |k: Key, bm: &Option<(Key, Option<usize>, Option<usize>)>| {
                k > k_now && bm.as_ref().is_none_or(|b| k > b.0)
            };
            for a in 0..n_search {
                let Some(ca) = &copies[a] else { continue };
                if chosen[a] {
                    // drop a
                    if k0 > 1 {
                        st.apply(ca, -1);
                        let k = key(st.wins, f0 - frag(a), d0 - span(a));
                        if ok(&st) && better(k, &best_move) {
                            best_move = Some((k, Some(a), None));
                        }
                        // swap a for b
                        for b in 0..n_search {
                            if chosen[b] {
                                continue;
                            }
                            let Some(cb) = &copies[b] else { continue };
                            st.apply(cb, 1);
                            let k = key(st.wins, f0 + frag(b) - frag(a), d0 + span(b) - span(a));
                            if ok(&st) && (st.wins > w0 || long(b)) && better(k, &best_move) {
                                best_move = Some((k, Some(a), Some(b)));
                            }
                            st.apply(cb, -1);
                        }
                        st.apply(ca, 1);
                    }
                } else if k0 < p.max_copies {
                    // add a
                    st.apply(ca, 1);
                    let k = key(st.wins, f0 + frag(a), d0 + span(a));
                    if ok(&st) && (st.wins > w0 || long(a)) && better(k, &best_move) {
                        best_move = Some((k, None, Some(a)));
                    }
                    st.apply(ca, -1);
                }
            }
            if let Some((_, out, inn)) = best_move {
                if let Some(a) = out {
                    st.apply(copies[a].as_ref().unwrap(), -1);
                    chosen[a] = false;
                }
                if let Some(b) = inn {
                    st.apply(copies[b].as_ref().unwrap(), 1);
                    chosen[b] = true;
                }
                continue;
            }
            if !p.pairs {
                break;
            }
            // stuck: compound moves, the best improving one taken
            let ins: Vec<usize> = (0..n_search).filter(|&k| !chosen[k] && copies[k].is_some()).collect();
            let outs: Vec<usize> = (0..pool.len()).filter(|&k| chosen[k]).collect();
            let mut best2: Option<(Key, Vec<usize>, Vec<usize>)> = None;
            let consider = |st: &State, o: Vec<usize>, i: Vec<usize>, best2: &mut Option<(Key, Vec<usize>, Vec<usize>)>| {
                let w = st.wins;
                let d = d0 + i.iter().map(|&k| span(k)).sum::<usize>() - o.iter().map(|&k| span(k)).sum::<usize>();
                let f = f0 + i.iter().map(|&k| frag(k)).sum::<usize>() - o.iter().map(|&k| frag(k)).sum::<usize>();
                let k = key(w, f, d);
                if ok(st) && k > k_now && (w > w0 || i.iter().all(|&x| long(x)))
                    && best2.as_ref().is_none_or(|b| k > b.0)
                {
                    *best2 = Some((k, o, i));
                }
            };
            // add two
            if k0 + 2 <= p.max_copies {
                for (x, &a) in ins.iter().enumerate() {
                    st.apply(copies[a].as_ref().unwrap(), 1);
                    for &b in &ins[x + 1..] {
                        st.apply(copies[b].as_ref().unwrap(), 1);
                        consider(&st, vec![], vec![a, b], &mut best2);
                        st.apply(copies[b].as_ref().unwrap(), -1);
                    }
                    st.apply(copies[a].as_ref().unwrap(), -1);
                }
            }
            // drop two
            if k0 > 2 {
                for (x, &a) in outs.iter().enumerate() {
                    st.apply(copies[a].as_ref().unwrap(), -1);
                    for &b in &outs[x + 1..] {
                        st.apply(copies[b].as_ref().unwrap(), -1);
                        consider(&st, vec![a, b], vec![], &mut best2);
                        st.apply(copies[b].as_ref().unwrap(), 1);
                    }
                    st.apply(copies[a].as_ref().unwrap(), 1);
                }
            }
            // one out, two in
            if k0 + 1 <= p.max_copies {
                for &a in &outs {
                    st.apply(copies[a].as_ref().unwrap(), -1);
                    for (x, &b) in ins.iter().enumerate() {
                        st.apply(copies[b].as_ref().unwrap(), 1);
                        for &c in &ins[x + 1..] {
                            st.apply(copies[c].as_ref().unwrap(), 1);
                            consider(&st, vec![a], vec![b, c], &mut best2);
                            st.apply(copies[c].as_ref().unwrap(), -1);
                        }
                        st.apply(copies[b].as_ref().unwrap(), -1);
                    }
                    st.apply(copies[a].as_ref().unwrap(), 1);
                }
            }
            let Some((_, o, i)) = best2 else { break };
            for a in o {
                st.apply(copies[a].as_ref().unwrap(), -1);
                chosen[a] = false;
            }
            for b in i {
                st.apply(copies[b].as_ref().unwrap(), 1);
                chosen[b] = true;
            }
        }
        // phase 2: add copies, closest-ranked first, within the margin
        let phase1 = st.wins;
        let floor = st.wins - (p.margin_per_kb * len as f64 / 1000.0) as i64;
        let mut order: Vec<usize> = (0..pool.len()).filter(|&k| !chosen[k] && copies[k].is_some()).collect();
        order.sort_by_key(|&k| rank_of.get(&pool[k]).copied().unwrap_or(usize::MAX));
        // A copy added for depth may lose no more library positions than the
        // margin allows for its own length (5 for a 5 kb copy at 1 per kb, none
        // for a 100 bp fragment), and the set no more than the margin allows
        // for the library's length in all.
        for k in order {
            if n_chosen(&chosen) >= p.max_copies {
                break;
            }
            if !long(k) {
                continue;
            }
            let cp = copies[k].as_ref().unwrap();
            let before = st.wins;
            st.apply(cp, 1);
            let allow = (p.margin_per_kb * span(k) as f64 / 1000.0) as i64;
            if before - st.wins <= allow && st.wins >= floor && ok(&st) {
                chosen[k] = true;
            } else {
                st.apply(cp, -1);
            }
        }
        // Two copies give no majority: every column they disagree on is a tie
        // the caller breaks by its own rule.  Two stay only when each is needed
        // for coverage; otherwise a third copy is added within the margin, or
        // the set drops to the one copy that keeps coverage and wins most.
        if n_chosen(&chosen) == 2 {
            let two: Vec<usize> = (0..pool.len()).filter(|&k| chosen[k]).collect();
            let mut removable = Vec::new();
            for &k in &two {
                st.apply(copies[k].as_ref().unwrap(), -1);
                if ok(&st) {
                    removable.push((st.wins, k));
                }
                st.apply(copies[k].as_ref().unwrap(), 1);
            }
            if !removable.is_empty() {
                let mut third: Option<(i64, usize, usize)> = None;
                for k in 0..pool.len() {
                    if chosen[k] || copies[k].is_none() || !long(k) {
                        continue;
                    }
                    st.apply(copies[k].as_ref().unwrap(), 1);
                    if st.wins >= floor && ok(&st)
                        && third.is_none_or(|t| st.wins > t.0 || (st.wins == t.0 && span(k) > t.1))
                    {
                        third = Some((st.wins, span(k), k));
                    }
                    st.apply(copies[k].as_ref().unwrap(), -1);
                }
                if let Some((_, _, k)) = third {
                    st.apply(copies[k].as_ref().unwrap(), 1);
                    chosen[k] = true;
                } else {
                    let (_, k) = *removable.iter().max().unwrap();
                    st.apply(copies[k].as_ref().unwrap(), -1);
                    chosen[k] = false;
                }
            }
        }
        let set: Vec<usize> = (0..pool.len()).filter(|&k| chosen[k]).map(|k| pool[k]).collect();
        let lost = judged - st.wins;
        if best.as_ref().is_none_or(|b| phase1 > b.0 || (phase1 == b.0 && set.len() > b.1.len())) {
            best = Some((phase1, set, lost));
        }
    }
    best.map(|(_, set, lost)| (set, lost))
}
