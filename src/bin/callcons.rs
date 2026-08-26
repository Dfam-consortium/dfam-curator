//! callcons — call a consensus sequence from an IG or Stockholm alignment.
//!
//! Two callers are exposed as subcommands so they can be compared directly on
//! the same input:
//!
//! * `acons` — faithful port of the GIRI `acons` caller ([`dfam_curator::giri`]).
//! * `dfam`  — the Dfam / `MultAln.pm` caller ([`dfam_curator::consensus`]).
//!
//! Output is always FASTA on stdout.
//!
//! # Why this reads files itself
//!
//! The crate's `io::read_alignment` normalises flanking gaps to padding, which
//! is equivalent to `acons -t`.  Reproducing `acons` exactly requires the file
//! bytes verbatim, so this binary parses IG and Stockholm directly.

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use std::io::Write;
use std::path::{Path, PathBuf};

use aln_core::consensus::ConsensusParams;
use aln_core::giri;

#[derive(Parser, Debug)]
#[command(
    name = "callcons",
    about = "Call consensus sequences from IG or Stockholm alignments",
    version
)]
struct Args {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand, Debug)]
enum Cmd {
    /// GIRI acons-compatible consensus caller.
    Acons(AconsArgs),
    /// Dfam (MultAln.pm lineage) consensus caller.
    Dfam(DfamArgs),
}

#[derive(Parser, Debug)]
struct AconsArgs {
    /// Input alignment (IG or Stockholm).
    input: PathBuf,

    /// Base name for the consensus (acons -N).
    #[arg(short = 'N', long, default_value = "CON")]
    name: String,

    /// Mammalian sequences: restore CpG doublets where possible (acons --mam).
    #[arg(long)]
    mam: bool,

    /// Use the bug-fixed CpG restorer ([`giri::fixed`]) instead of the faithful
    /// acons port: corrects the inert `first_base` lookup and the masked-CpG
    /// over-skip.  Only affects `--mam`.
    #[arg(long)]
    fixed: bool,

    /// Trim flanking gap symbols padding aligned sequences (acons -t).
    #[arg(short = 't', long)]
    trim: bool,

    /// Minimum number of non-gap symbols in an alignment column (acons --min).
    #[arg(long, default_value_t = 2, value_name = "N")]
    min: usize,

    /// Emit the full-length gapped consensus (one char per input column)
    /// instead of the gap-stripped sequence.  For per-column evaluation.
    #[arg(long)]
    gapped: bool,
}

#[derive(Parser, Debug)]
struct DfamArgs {
    /// Input alignment (IG or Stockholm).
    input: PathBuf,

    /// Do not infer CpG sites (skip the CpG correction pass).
    #[arg(long)]
    no_cpg: bool,

    /// Emit the full-length gapped consensus (one char per input column)
    /// instead of the gap-stripped sequence.  For per-column evaluation.
    #[arg(long)]
    gapped: bool,
}

fn main() -> Result<()> {
    match Args::parse().cmd {
        Cmd::Acons(a) => {
            let rows = read_rows(&a.input)?;
            let seqs: Vec<Vec<u8>> = if a.trim {
                rows.iter().map(|r| giri::trim_flanking_gaps(&r.seq)).collect()
            } else {
                rows.iter().map(|r| r.seq.clone()).collect()
            };
            let refs: Vec<&[u8]> = seqs.iter().map(|v| v.as_slice()).collect();

            let mut cons = giri::get_consensus(&refs, a.min);
            if a.mam {
                // acons restores CpG on the still-gapped consensus, before
                // gaps are stripped.
                if a.fixed {
                    giri::fixed::restore_cpg(&mut cons, &refs);
                } else {
                    giri::restore_cpg(&mut cons, &refs);
                }
            }
            let cons = if a.gapped { cons } else { giri::strip_gaps(&cons) };
            write_fasta(&a.name, &cons)
        }
        Cmd::Dfam(d) => {
            let rows = read_rows(&d.input)?;
            let seqs: Vec<&[u8]> = rows.iter().map(|r| r.seq.as_slice()).collect();
            let params = ConsensusParams {
                // acons has no reference row; use every sequence so the two
                // subcommands are directly comparable on one file.
                include_reference: true,
                enable_cpg: !d.no_cpg,
                ..Default::default()
            };
            let cons = aln_core::consensus::build_consensus_from_sequences(&seqs, &params);
            let cons: Vec<u8> = if d.gapped {
                cons
            } else {
                cons.into_iter().filter(|&b| b != b'-').collect()
            };
            write_fasta("CON", &cons)
        }
    }
}

// ── Output ───────────────────────────────────────────────────────────────────

fn write_fasta(name: &str, seq: &[u8]) -> Result<()> {
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    writeln!(out, ">{}", name)?;
    for chunk in seq.chunks(60) {
        out.write_all(chunk)?;
        writeln!(out)?;
    }
    if seq.is_empty() {
        writeln!(out)?;
    }
    Ok(())
}

// ── Verbatim readers ─────────────────────────────────────────────────────────

struct Row {
    #[allow(dead_code)]
    name: String,
    seq: Vec<u8>,
}

fn read_rows(path: &Path) -> Result<Vec<Row>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading {}", path.display()))?;
    let rows = if text
        .lines()
        .find(|l| !l.trim().is_empty())
        .is_some_and(|l| l.starts_with("# STOCKHOLM"))
    {
        read_stockholm(&text)
    } else {
        read_ig(&text)
    };
    if rows.is_empty() {
        bail!("no sequences found in {}", path.display());
    }
    let w = rows[0].seq.len();
    if let Some(bad) = rows.iter().find(|r| r.seq.len() != w) {
        eprintln!(
            "warning: ragged alignment — '{}' is {} columns, first row is {}",
            bad.name,
            bad.seq.len(),
            w
        );
    }
    Ok(rows)
}

/// Parse IG: `;`-prefixed comment lines, then a bare identifier line, then the
/// gapped sequence.  Bytes are preserved verbatim — in particular `-` stays a
/// gap (GIRI `GAPCHAR`) and unknown symbols such as `.` are left alone, since
/// `acons` scores them as 0 but still counts them toward column coverage.
fn read_ig(text: &str) -> Vec<Row> {
    let mut rows = Vec::new();
    let mut name: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();
    let mut expect_name = true;

    for line in text.lines() {
        let t = line.trim();
        if t.starts_with(';') {
            if let Some(n) = name.take() {
                rows.push(Row { name: n, seq: std::mem::take(&mut seq) });
            }
            expect_name = true;
            continue;
        }
        if t.is_empty() {
            continue;
        }
        if expect_name {
            name = Some(t.to_string());
            expect_name = false;
        } else {
            seq.extend_from_slice(t.as_bytes());
        }
    }
    if let Some(n) = name.take() {
        rows.push(Row { name: n, seq });
    }
    rows
}

/// Parse Stockholm sequence lines.  `.` denotes columns outside a sequence's
/// aligned region, so it maps to padding (`b' '`); `-` stays an interior gap.
fn read_stockholm(text: &str) -> Vec<Row> {
    let mut order: Vec<String> = Vec::new();
    let mut map: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();

    for line in text.lines() {
        let t = line.trim_end();
        if t.is_empty() || t.starts_with('#') || t.starts_with("//") {
            continue;
        }
        let mut it = t.split_whitespace();
        let (Some(n), Some(s)) = (it.next(), it.next()) else { continue };
        let bytes: Vec<u8> = s
            .bytes()
            .map(|b| if b == b'.' { giri::PAD } else { b })
            .collect();
        if !map.contains_key(n) {
            order.push(n.to_string());
        }
        map.entry(n.to_string()).or_default().extend_from_slice(&bytes);
    }
    order
        .into_iter()
        .map(|n| {
            let seq = map.remove(&n).unwrap_or_default();
            Row { name: n, seq }
        })
        .collect()
}
