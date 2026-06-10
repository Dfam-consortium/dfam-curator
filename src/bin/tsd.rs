/// tsd — Target Site Duplication finder for TE families.
///
/// Reads a multiple sequence alignment in Stockholm, FASTA, or crossmatch
/// format, identifies full-length instances, extracts genomic flanking
/// sequence, and scores left/right flank pairs at each candidate TSD length.
use anyhow::Context;
use clap::Parser;
use std::collections::HashMap;
use std::path::PathBuf;

use dfam_curator::{
    alignment::Orientation,
    io::read_alignment,
    matrix,
};

// ── TSD scoring matrix ────────────────────────────────────────────────────────

/// Rows/cols indexed by matrix::alpha_idx(): A R G C Y T K M S W N X
/// Values from the Dfam `tsdmatrix` file (FREQS A=0.295 C=0.205 G=0.205 T=0.295).
#[rustfmt::skip]
static TSD_MATRIX: [[i32; 12]; 12] = [
    //    A    R    G    C    Y    T    K    M    S    W    N    X
    [    9,   0,  -8, -15, -16, -17, -13,  -3, -11,  -4,  -2,  -7], // A
    [    3,   3,   3, -15, -15, -16,  -7,  -6,  -6,  -7,  -2,  -7], // R
    [   -5,   3,  10, -14, -14, -15,  -2,  -9,  -2,  -9,  -2,  -7], // G
    [  -15, -14, -14,  10,   3,  -5,  -9,  -2,  -2,  -9,  -2,  -7], // C
    [  -16, -15, -15,   3,   3,   3,  -6,  -7,  -6,  -7,  -2,  -7], // Y
    [  -17, -16, -15,  -8,   0,   9,  -3, -13, -11,  -4,  -2,  -7], // T
    [  -11,  -6,  -2, -11,  -7,  -2,  -2, -11,  -6,  -7,  -2,  -7], // K
    [   -2,  -7, -11,  -2,  -6, -11, -11,  -2,  -6,  -7,  -2,  -7], // M
    [   -9,  -5,  -2,  -2,  -5,  -9,  -5,  -5,  -2,  -9,  -2,  -7], // S
    [   -3,  -8, -11, -11,  -8,  -3,  -8,  -8, -11,  -4,  -2,  -7], // W
    [   -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -7], // N
    [   -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7], // X
];

/// Mirrors Perl `calc()` from TSD.pl: scores one left-flank vs right-flank base.
/// N or gap vs anything → 0.  Match: AT=9, GC=10.
/// Mismatch: purine-purine or pyrimidine-pyrimidine=-6, AT transversion=-17, other=-15.
fn tsd_pair_score(l: u8, r: u8) -> i32 {
    let l = l.to_ascii_uppercase();
    let r = r.to_ascii_uppercase();
    if l == b'N' || r == b'N' || l == b'-' || r == b'-' {
        return 0;
    }
    if l == r {
        return if l == b'A' || l == b'T' { 9 } else { 10 };
    }
    let is_purine = |b: u8| b == b'A' || b == b'G';
    let is_pyrim  = |b: u8| b == b'C' || b == b'T';
    let is_at     = |b: u8| b == b'A' || b == b'T';
    if (is_purine(l) && is_purine(r)) || (is_pyrim(l) && is_pyrim(r)) {
        -6
    } else if is_at(l) && is_at(r) {
        -17
    } else {
        -15
    }
}

// ── CLI ───────────────────────────────────────────────────────────────────────

#[derive(Parser, Debug)]
#[command(
    name = "tsd",
    about = "Detect target site duplications from an MSA + reference genome",
    after_help = "\
Identifies full-length TE copies in the MSA, extracts genomic flanking
sequence from the reference, and scores left/right flank pairs at each
candidate TSD length.

Score table columns: LENGTH  SCORE  [OFFSET with --all-offsets]
OFFSET = number of TSD bases falling inside the alignment boundaries
(0 = TSD is entirely external; OFFSET=LENGTH = fully internal).

The edge view printed after the consensus shows each contributing instance:
  lowercase  = non-TSD genomic flank
  UPPERCASE in flank region = TSD bases
  UPPERCASE in MSA region   = first/last alignment columns

PREREQUISITES: sequence identifiers must be in chrom:start-end_+/- format.
Run `stk lint --genome` to validate coordinates before using this command.

Examples:
  tsd --genome assembly.2bit family.stk
  tsd --genome assembly.2bit family.fa
  tsd --genome assembly.2bit --max-flank 15 --tolerance 5 family.stk
  tsd --genome assembly.2bit --scores-only family.stk
  tsd --genome assembly.2bit --all-offsets family.stk",
    version,
)]
struct Args {
    /// Input alignment (Stockholm, FASTA/A2M, or crossmatch .align).
    input: PathBuf,

    /// Reference genome (UCSC .2bit preferred; FASTA also accepted).
    #[arg(long, value_name = "FILE")]
    genome: PathBuf,

    /// Maximum flank length to test (bp).
    #[arg(long, default_value = "12")]
    max_flank: usize,

    /// Allow up to this many alignment columns of truncation from either end
    /// when deciding whether an instance is full-length.
    #[arg(long, default_value = "3")]
    tolerance: usize,

    /// Print only the per-length score table; skip TSD pairs, consensus, and edge view.
    #[arg(long)]
    scores_only: bool,

    /// Also search non-zero offsets (0 < k ≤ L) to detect TSDs embedded inside
    /// the alignment boundaries.  Off by default because for DNA transposons the
    /// internal comparison can produce false-positive signal from terminal
    /// inverted repeats (TIRs).
    #[arg(long)]
    all_offsets: bool,

    /// Number of MSA columns to display on each side in the edge view.
    #[arg(long, default_value = "10")]
    edge_cols: usize,
}

// ── Genome access ─────────────────────────────────────────────────────────────

enum Genome {
    TwoBit(dfam_curator::io::twobit::TwoBitReader),
    InMem(HashMap<String, Vec<u8>>),
}

impl Genome {
    fn fetch(&self, chrom: &str, start: u64, end: u64) -> Option<Vec<u8>> {
        match self {
            Genome::TwoBit(r) => r.fetch(chrom, start, end).ok(),
            Genome::InMem(m)  => {
                let seq = m.get(chrom)?;
                let (s, e) = (start as usize, end as usize);
                if e <= seq.len() { Some(seq[s..e].to_vec()) } else { None }
            }
        }
    }
    fn contains(&self, chrom: &str) -> bool {
        match self {
            Genome::TwoBit(r) => r.contains(chrom),
            Genome::InMem(m)  => m.contains_key(chrom),
        }
    }
}

// ── Window ────────────────────────────────────────────────────────────────────

/// Genomic windows around one full-length instance plus its MSA row.
struct Window {
    /// genome[lb - max_flank .. lb + max_flank]  (lb = seq_start-1, 0-based)
    lw: Vec<u8>,
    /// genome[rb - max_flank .. rb + max_flank]  (rb = seq_end, 0-based exclusive)
    rw: Vec<u8>,
    /// Aligned row from the MSA (gap characters preserved).
    seq: Vec<u8>,
    /// Display label: chrom:start-end_orient
    label: String,
}

// ── main ──────────────────────────────────────────────────────────────────────

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    let msa = read_alignment(&args.input)
        .with_context(|| format!("cannot read {:?}", args.input))?;

    let width     = msa.width();
    let max_flank = args.max_flank;
    let tol       = args.tolerance;
    let mf        = max_flank as u64;

    // Open the genome.
    let genome_str = args.genome.to_str()
        .ok_or_else(|| anyhow::anyhow!("genome path is not valid UTF-8"))?;
    let genome = if genome_str.ends_with(".2bit") {
        eprintln!("Indexing 2bit file: {}", genome_str);
        Genome::TwoBit(
            dfam_curator::io::twobit::TwoBitReader::open(&args.genome)
                .with_context(|| format!("cannot open {:?}", args.genome))?,
        )
    } else {
        eprintln!("Loading reference genome: {}", genome_str);
        Genome::InMem(
            dfam_coord::load_reference(genome_str)
                .with_context(|| format!("cannot load {:?}", args.genome))?,
        )
    };

    // Collect windows for full-length instances with valid coordinates.
    let mut windows: Vec<Window> = Vec::new();
    let mut skipped_partial   = 0usize;
    let mut skipped_no_coord  = 0usize;
    let mut skipped_no_chrom  = 0usize;
    let mut skipped_bounds    = 0usize;

    for inst in &msa.sequences[1..] {
        // Use the position of the first/last actual base (not gap or space) so
        // that partial copies whose truncated ends are encoded as '-' rather
        // than ' ' are correctly excluded.
        let first_base = inst.seq.iter().position(|&b| b != b'-' && b != b' ')
            .unwrap_or(width);
        let last_base  = inst.seq.iter().rposition(|&b| b != b'-' && b != b' ')
            .unwrap_or(0);
        if first_base > tol || last_base + tol + 1 < width {
            skipped_partial += 1;
            continue;
        }
        if inst.seq_start >= inst.seq_end {
            skipped_no_coord += 1;
            continue;
        }
        let lb = inst.seq_start - 1; // 0-based left boundary (exclusive)
        let rb = inst.seq_end;       // 0-based right boundary (exclusive)

        if !genome.contains(&inst.name) {
            skipped_no_chrom += 1;
            continue;
        }
        if lb < mf || rb < mf {
            skipped_bounds += 1;
            continue;
        }
        let lw = genome.fetch(&inst.name, lb - mf, lb + mf);
        let rw = genome.fetch(&inst.name, rb - mf, rb + mf);
        match (lw, rw) {
            (Some(l), Some(r)) if l.len() == 2 * max_flank && r.len() == 2 * max_flank => {
                let orient = match inst.orient {
                    Orientation::Forward => '+',
                    Orientation::Reverse => '-',
                };
                let label = format!("{}:{}-{}_{}", inst.name, inst.seq_start, inst.seq_end, orient);
                windows.push(Window { lw: l, rw: r, seq: inst.seq.clone(), label });
            }
            _ => { skipped_bounds += 1; }
        }
    }

    eprintln!(
        "{} full-length instances ({} partial, {} no-coords, {} chrom-not-found, {} near-edge skipped)",
        windows.len(), skipped_partial, skipped_no_coord, skipped_no_chrom, skipped_bounds
    );

    if windows.is_empty() {
        anyhow::bail!("no full-length instances with valid genomic coordinates found");
    }

    // Phase 1: score every (length L, offset k) pair.
    //
    // For offset k and length L, TSD candidate bases are:
    //   left:  lw[mf - L + k .. mf + k)
    //   right: rw[mf - k     .. mf - k + L)
    // Background correction: +6.9 per base (expected score at ~42-45% GC).
    let mut score_at = vec![vec![0.0f64; max_flank + 1]; max_flank + 1];

    for w in &windows {
        for l in 1..=max_flank {
            let k_max = if args.all_offsets { l } else { 0 };
            for k in 0..=k_max {
                let mut raw: i32 = 0;
                for j in 0..l {
                    raw += tsd_pair_score(w.lw[max_flank - l + k + j], w.rw[max_flank - k + j]);
                }
                score_at[l][k] += raw as f64 + l as f64 * 6.9;
            }
        }
    }

    let mut best_len    = 1usize;
    let mut best_offset = 0usize;
    let mut best_avg    = i64::MIN;

    for l in 1..=max_flank {
        let k_max = if args.all_offsets { l } else { 0 };
        let (k_opt, avg_opt) = (0..=k_max)
            .map(|k| (k, (score_at[l][k] / l as f64 + 0.5) as i64))
            .max_by_key(|&(_, s)| s)
            .unwrap();
        if args.all_offsets {
            println!("{} {} {}", l, avg_opt, k_opt);
        } else {
            println!("{} {}", l, avg_opt);
        }
        if avg_opt >= best_avg {
            best_avg    = avg_opt;
            best_len    = l;
            best_offset = k_opt;
        }
    }

    if args.all_offsets {
        let desc = if best_offset == 0 {
            "TSD is external to the alignment".to_string()
        } else if best_offset == best_len {
            "TSD appears to be fully internal — may reflect TIR structure \
             for DNA transposons rather than a genuine TSD".to_string()
        } else {
            format!("TSD partially overlaps alignment boundary ({}/{} bases internal)",
                    best_offset, best_len)
        };
        println!("\nBest TSD: length={}, offset={} ({})\n", best_len, best_offset, desc);
    } else {
        println!();
    }

    if args.scores_only {
        return Ok(());
    }

    // Phase 2: re-score at (best_len, best_offset), collect pairs, build IUPAC consensus.
    let mut count = vec![[0u32; 12]; best_len];
    let mut passing: Vec<(&Window, Vec<u8>, Vec<u8>)> = Vec::new();

    for w in &windows {
        let mut raw: i32 = 0;
        for j in 0..best_len {
            raw += tsd_pair_score(
                w.lw[max_flank - best_len + best_offset + j],
                w.rw[max_flank - best_offset + j],
            );
        }
        // Require a positive raw score (genuine net matching) rather than
        // just above the background-corrected zero, which is far too permissive.
        if raw <= 0 {
            continue;
        }

        let lslice = w.lw[max_flank - best_len + best_offset .. max_flank + best_offset].to_vec();
        let rslice = w.rw[max_flank - best_offset .. max_flank - best_offset + best_len].to_vec();

        let ls: String = lslice.iter().map(|&b| b as char).collect();
        let rs: String = rslice.iter().map(|&b| b as char).collect();
        println!("{}\n{}", ls, rs);

        for i in 0..best_len {
            for b in [lslice[i].to_ascii_uppercase(), rslice[i].to_ascii_uppercase()] {
                if let Some(idx) = matrix::alpha_idx(b) {
                    if idx < 12 { count[i][idx] += 1; }
                }
            }
        }
        passing.push((w, lslice, rslice));
    }
    println!();

    // Consensus.
    let mut cons = String::with_capacity(best_len);
    for i in 0..best_len {
        let mut best_char = b'N';
        let mut max_score = i32::MIN;
        for row in 0..12usize {
            let score: i32 = (0..12)
                .map(|col| TSD_MATRIX[row][col] * count[i][col] as i32)
                .sum();
            if score > max_score { max_score = score; best_char = matrix::alpha_byte(row); }
        }
        cons.push(best_char as char);
    }
    println!("{}\n", cons);

    // Phase 3: edge view — one line per passing instance.
    //
    // Layout (lowercase=flank, UPPERCASE=TSD or MSA):
    //   [left_non_tsd][LEFT_TSD][left_msa_edge]  ...  [right_msa_edge][RIGHT_TSD][right_non_tsd]
    //
    // lw window: [lb-mf .. lb+mf].  Non-TSD left flank = lw[0 .. mf-L+k].  TSD = lw[mf-L+k .. mf+k].
    // rw window: [rb-mf .. rb+mf].  TSD = rw[mf-k .. mf-k+L].  Non-TSD right flank = rw[mf-k+L ..].
    let ec        = args.edge_cols;
    let fnl       = max_flank - best_len + best_offset; // non-TSD flank length each side
    let name_w    = passing.iter().map(|(w, _, _)| w.label.len()).max().unwrap_or(8);

    // Header: two lines.  The +1 accounts for the separator space between TSD and MSA.
    let left_w  = fnl + best_len + 1 + ec;
    let right_w = ec + 1 + best_len + fnl;
    let sep     = " ... ";

    {
        // Line 1: "Left" and "Right" labels positioned over their halves.
        let left_label  = "Left";
        let right_label = "Right";
        let left_pad    = (left_w.saturating_sub(left_label.len())) / 2;
        let right_pad   = (right_w.saturating_sub(right_label.len())) / 2;
        println!(
            "{:<name_w$}  {:<left_w$}{}{:>right_w$}",
            "",
            format!("{:>pad$}{}", "", left_label, pad = left_pad),
            sep,
            format!("{:<pad$}{}", right_label, "", pad = right_pad),
            name_w = name_w,
            left_w = left_w,
            right_w = right_w,
        );

        // Line 2: component breakdown.
        let l_flank = format!("{:^width$}", "flank", width = fnl.max(5));
        let l_tsd   = format!("{:^width$}", "TSD",   width = best_len.max(3));
        let l_msa   = format!("{:^width$}", "MSA",   width = ec.max(3));
        let r_msa   = format!("{:^width$}", "MSA",   width = ec.max(3));
        let r_tsd   = format!("{:^width$}", "TSD",   width = best_len.max(3));
        let r_flank = format!("{:^width$}", "flank", width = fnl.max(5));
        println!(
            "{:<name_w$}  {}{} {}{}{} {}{}",
            "",
            &l_flank[..fnl.min(l_flank.len())],
            &l_tsd[..best_len.min(l_tsd.len())],
            &l_msa[..ec.min(l_msa.len())],
            sep,
            &r_msa[..ec.min(r_msa.len())],
            &r_tsd[..best_len.min(r_tsd.len())],
            &r_flank[..fnl.min(r_flank.len())],
            name_w = name_w,
        );
    }

    for (w, lslice, rslice) in &passing {
        // Left non-TSD flank (lowercase).
        let left_flank: String = w.lw[..fnl].iter().map(|&b| b.to_ascii_lowercase() as char).collect();
        // Left TSD (uppercase).
        let left_tsd: String  = lslice.iter().map(|&b| b.to_ascii_uppercase() as char).collect();
        // Left MSA edge: first ec alignment columns.
        let left_msa: String  = w.seq[..ec.min(w.seq.len())]
            .iter().map(|&b| if b == b'-' { '-' } else { b.to_ascii_uppercase() as char }).collect();

        // Right MSA edge: last ec alignment columns.
        let seq_len = w.seq.len();
        let right_msa: String = w.seq[seq_len.saturating_sub(ec)..]
            .iter().map(|&b| if b == b'-' { '-' } else { b.to_ascii_uppercase() as char }).collect();
        // Right TSD (uppercase).
        let right_tsd: String = rslice.iter().map(|&b| b.to_ascii_uppercase() as char).collect();
        // Right non-TSD flank (lowercase).
        let rw_start = max_flank - best_offset + best_len;
        let right_flank: String = w.rw[rw_start..].iter().map(|&b| b.to_ascii_lowercase() as char).collect();

        println!(
            "{:<name_w$}  {}{} {}{}{} {}{}",
            w.label,
            left_flank, left_tsd, left_msa,
            sep,
            right_msa, right_tsd, right_flank,
            name_w = name_w,
        );
    }

    Ok(())
}
