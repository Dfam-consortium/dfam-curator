/// linup — MSA viewer (Rust port of RepeatModeler/util/Linup)
///
/// Reads a multiple sequence alignment in Stockholm, FASTA/A2M, or crossmatch
/// format and renders it as a pretty-printed alignment or Stockholm output.
///
/// When --blast-tab is set the input file is treated as rmblastn tabular output
/// (outfmt: `6 score qseqid qstart qend qlen sstrand sseqid sstart send slen qseq sseq`)
/// and --ref-seq must point to the FASTA file that was used as the BLAST subject.
use anyhow::Context;
use clap::{Parser, ValueEnum};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use dfam_curator::{
    alignment::MultiAlign,
    consensus::ConsensusParams,
    io::read_alignment,
};

#[derive(Debug, Clone, Copy, ValueEnum)]
enum OutFormat {
    /// Perl Linup–compatible pretty-print blocks (default).
    Linup,
    /// Stockholm 1.0 with #=GC RF consensus line.
    Stockholm,
    /// Print only the consensus sequence.
    Consensus,
    /// Print per-sequence statistics.
    Stats,
}

#[derive(Parser, Debug)]
#[command(
    name = "linup",
    about = "MSA viewer / converter for repeat element alignments",
    version
)]
struct Args {
    /// Input alignment file (Stockholm, FASTA/A2M, crossmatch .align, or BLAST tabular
    /// when --blast-tab is set).
    input: PathBuf,

    /// Treat input as rmblastn tabular output rather than an alignment file.
    /// Requires --ref-seq.
    #[arg(long)]
    blast_tab: bool,

    /// Reference sequence FASTA used as the BLAST subject (required with --blast-tab).
    /// The first sequence in this file becomes the MSA reference row.
    #[arg(long)]
    ref_seq: Option<PathBuf>,

    /// Output format (default: pretty-print to stdout).
    #[arg(long, value_enum)]
    format: Option<OutFormat>,

    /// Include the reference sequence when calling the consensus.
    #[arg(short = 'i', long)]
    include_ref: bool,

    /// Trim this many ungapped reference bp from the left.
    #[arg(long, default_value_t = 0)]
    trim_left: usize,

    /// Trim this many ungapped reference bp from the right.
    #[arg(long, default_value_t = 0)]
    trim_right: usize,

    /// Trim ambiguous (non-ACGT) bases from both ends of the consensus.
    /// Mutually exclusive with --trim-left/--trim-right, --sub-align, and --revcomp.
    #[arg(long)]
    trim_ambig: bool,

    /// Slice the alignment to the given 1-based, fully-closed consensus coordinate
    /// range (e.g. --sub-align 10-200).  Requires a stable RF/consensus.
    /// Mutually exclusive with --trim-left/--trim-right, --trim-ambig, and --revcomp.
    #[arg(long, value_name = "START-END")]
    sub_align: Option<String>,

    /// Minimum number of bases a sequence must retain after --sub-align to be kept.
    /// Only valid with --sub-align.
    #[arg(long, value_name = "N")]
    min_len: Option<usize>,

    /// Reverse-complement the entire alignment before output.
    /// Mutually exclusive with all trim/slice options.
    #[arg(long)]
    revcomp: bool,

    /// Override the family name / ID in output.
    #[arg(long)]
    name: Option<String>,

    /// With --format consensus: include gap characters in output (default: strip gaps).
    #[arg(long)]
    include_gaps: bool,

    /// Select a specific record from a multi-record Stockholm file.
    /// Only applies to STK-format input; ignored for all other formats.
    /// A purely numeric value selects by 1-based record number (e.g. --select 3).
    /// Any non-numeric value selects by exact #=GF ID match (e.g. --select MyFam).
    #[arg(long, value_name = "SELECT")]
    select: Option<String>,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    // ── Validate mutually exclusive options ───────────────────────────────────
    let has_trim_lr = args.trim_left > 0 || args.trim_right > 0;
    if has_trim_lr && (args.trim_ambig || args.sub_align.is_some()) {
        anyhow::bail!(
            "--trim-left/--trim-right and --trim-ambig/--sub-align are mutually exclusive"
        );
    }
    if args.trim_ambig && args.sub_align.is_some() {
        anyhow::bail!("--trim-ambig and --sub-align are mutually exclusive");
    }
    if args.revcomp && (has_trim_lr || args.trim_ambig || args.sub_align.is_some()) {
        anyhow::bail!(
            "--revcomp cannot be combined with --trim-left/--trim-right/--trim-ambig/--sub-align"
        );
    }
    if args.min_len.is_some() && args.sub_align.is_none() {
        anyhow::bail!("--min-len only applies to --sub-align");
    }
    if args.select.is_some() && args.blast_tab {
        anyhow::bail!("--select cannot be combined with --blast-tab");
    }

    // ── Load alignment ────────────────────────────────────────────────────────
    let mut msa: MultiAlign = if args.blast_tab {
        let ref_path = args.ref_seq.as_ref()
            .context("--ref-seq is required when --blast-tab is set")?;
        let (ref_name, ref_seq) = read_first_fasta_seq(ref_path)?;
        let text = std::fs::read_to_string(&args.input)
            .with_context(|| format!("failed to read {}", args.input.display()))?;
        let hits = dfam_curator::blast::parse_tab(&text)?;
        dfam_curator::blast::hits_to_multialign(&ref_seq, &ref_name, &hits)
    } else if let Some(ref select) = args.select {
        let fmt = dfam_curator::io::detect_format(&args.input)
            .with_context(|| format!("failed to detect format of {}", args.input.display()))?;
        if fmt != dfam_curator::io::Format::Stockholm {
            anyhow::bail!("--select only applies to Stockholm (.stk) files");
        }
        dfam_curator::io::stockholm::read_select(&args.input, select)
            .with_context(|| format!("failed to read {}", args.input.display()))?
    } else {
        read_alignment(&args.input)
            .with_context(|| format!("failed to read {}", args.input.display()))?
    };

    // ── Trimming / slicing ────────────────────────────────────────────────────
    if has_trim_lr {
        msa.trim(args.trim_left, args.trim_right);
    }

    if args.trim_ambig {
        let cons = build_consensus(&msa, false);
        let width = msa.width();
        let left_col = cons.iter()
            .position(|&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
            .unwrap_or(width);
        let right_end = cons.iter()
            .rposition(|&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
            .map(|p| p + 1)
            .unwrap_or(0);
        let right_trim = width.saturating_sub(right_end);
        if left_col > 0 || right_trim > 0 {
            eprintln!("Ambiguous edges trimmed: left = {}, right = {}", left_col, right_trim);
            msa.slice_columns(left_col, right_end);
        }
    }

    if let Some(ref range_str) = args.sub_align {
        let (start_bp, end_bp) = parse_sub_align_range(range_str)?;

        // Check that RF and consensus agree before slicing.
        let check_cons = build_consensus(&msa, args.include_ref);
        let rf = msa.sequences.first()
            .map(|r| r.seq.as_slice())
            .unwrap_or(&[]);
        let is_occupancy_rf = rf.iter().all(|&b| matches!(b, b'x' | b'X' | b'-' | b'.' | b' '));
        let is_stable = if is_occupancy_rf {
            let cons_occ: Vec<u8> = check_cons.iter().map(|&b| {
                if b == b'-' || b == b'.' || b == b' ' { b } else { b'x' }
            }).collect();
            rf == cons_occ.as_slice()
        } else {
            rf == check_cons.as_slice()
        };
        if !is_stable {
            anyhow::bail!(
                "--sub-align can only be used on MSAs with a stable consensus; \
                 please update the RF line in {} before using --sub-align",
                args.input.display()
            );
        }

        // Use inclRef=false for locating columns, matching Perl's MultAln::slice().
        let slice_cons = build_consensus(&msa, false);
        let (start_col, end_col) = consensus_bp_to_cols(&slice_cons, start_bp, end_bp)
            .with_context(|| format!("--sub-align {}-{} out of range", start_bp, end_bp))?;
        msa.slice_columns(start_col, end_col + 1);

        if let Some(min_len) = args.min_len {
            if min_len > 0 {
                let mut i = 1;
                while i < msa.sequences.len() {
                    let n = msa.sequences[i].seq.iter()
                        .filter(|&&b| b.is_ascii_alphabetic())
                        .count();
                    if n < min_len {
                        msa.sequences.remove(i);
                    } else {
                        i += 1;
                    }
                }
                msa.invalidate_consensus();
            }
        }
    }

    // ── Orientation ───────────────────────────────────────────────────────────
    if args.revcomp {
        msa.reverse_complement();
    }

    // ── Consensus ─────────────────────────────────────────────────────────────
    let params = ConsensusParams {
        include_reference: args.include_ref,
        ..Default::default()
    };

    let raw_seqs: Vec<&[u8]> = {
        let start = if params.include_reference { 0 } else { 1 };
        msa.sequences[start..].iter().map(|s| s.seq.as_slice()).collect()
    };
    let consensus = dfam_curator::consensus::build_consensus_from_sequences(&raw_seqs, &params);

    // ── Output ────────────────────────────────────────────────────────────────
    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    let family_id = args.name.as_deref().or_else(|| {
        msa.reference().map(|r| r.name.as_str())
    });

    match args.format {
        None | Some(OutFormat::Linup) => {
            dfam_curator::io::linup_fmt::write(&msa, &consensus, &mut out, 100)?;
        }
        Some(OutFormat::Stockholm) => {
            dfam_curator::io::stockholm::write(
                &msa,
                &mut out,
                family_id,
                Some(&consensus),
                true,
            )?;
        }
        Some(OutFormat::Consensus) => {
            let name = family_id.unwrap_or("consensus");
            writeln!(out, ">{}", name)?;
            if args.include_gaps {
                out.write_all(&consensus)?;
            } else {
                let ungapped: Vec<u8> = consensus.iter().copied().filter(|&b| b != b'-').collect();
                out.write_all(&ungapped)?;
            }
            writeln!(out)?;
            writeln!(out)?;
        }
        Some(OutFormat::Stats) => {
            print_stats(&msa, &consensus, &mut out)?;
        }
    }

    Ok(())
}

/// Build a consensus from `msa` using the default CpG parameters.
/// When `include_ref` is false the reference row (index 0) is excluded.
fn build_consensus(msa: &MultiAlign, include_ref: bool) -> Vec<u8> {
    let params = ConsensusParams {
        include_reference: include_ref,
        ..Default::default()
    };
    let start = if include_ref { 0 } else { 1 };
    let raw: Vec<&[u8]> = msa.sequences[start..]
        .iter()
        .map(|s| s.seq.as_slice())
        .collect();
    dfam_curator::consensus::build_consensus_from_sequences(&raw, &params)
}

/// Parse "START-END" into (start, end) as 1-based integers.
fn parse_sub_align_range(s: &str) -> anyhow::Result<(usize, usize)> {
    let parts: Vec<&str> = s.splitn(2, '-').collect();
    if parts.len() != 2 {
        anyhow::bail!("--sub-align must be in the form START-END (e.g. 10-200)");
    }
    let start: usize = parts[0].parse()
        .context("--sub-align START is not a valid integer")?;
    let end: usize = parts[1].parse()
        .context("--sub-align END is not a valid integer")?;
    if start < 1 {
        anyhow::bail!("--sub-align start must be ≥ 1 (1-based coordinates)");
    }
    if start > end {
        anyhow::bail!("--sub-align start ({}) must be ≤ end ({})", start, end);
    }
    Ok((start, end))
}

/// Map 1-based, fully-closed consensus bp positions to (start_col, end_col) column
/// indices (both inclusive) by counting non-gap, non-space characters in `cons`.
fn consensus_bp_to_cols(cons: &[u8], start_bp: usize, end_bp: usize)
    -> anyhow::Result<(usize, usize)>
{
    let mut bp = 0usize;
    let mut start_col = None;
    let mut end_col   = None;
    for (i, &b) in cons.iter().enumerate() {
        if b == b'-' || b == b'.' || b == b' ' { continue; }
        bp += 1;
        if bp == start_bp { start_col = Some(i); }
        if bp == end_bp   { end_col   = Some(i); break; }
    }
    match (start_col, end_col) {
        (Some(s), Some(e)) => Ok((s, e)),
        _ => anyhow::bail!(
            "range {}-{} is beyond the alignment length ({} ungapped positions)",
            start_bp, end_bp, bp
        ),
    }
}

/// Read the first sequence from a FASTA file, returning (name, sequence).
fn read_first_fasta_seq(path: &PathBuf) -> anyhow::Result<(String, Vec<u8>)> {
    use std::io::BufRead;
    let f = std::io::BufReader::new(
        std::fs::File::open(path)
            .with_context(|| format!("cannot open {}", path.display()))?,
    );
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::new();
    for line in f.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() { continue; }
        if let Some(rest) = trimmed.strip_prefix('>') {
            if !name.is_empty() { break; }
            name = rest.split_whitespace().next().unwrap_or("ref").to_string();
        } else if !name.is_empty() {
            seq.extend_from_slice(trimmed.as_bytes());
        }
    }
    if name.is_empty() {
        anyhow::bail!("no sequences found in {}", path.display());
    }
    Ok((name, seq))
}

fn print_stats(msa: &MultiAlign, consensus: &[u8], out: &mut dyn Write) -> std::io::Result<()> {
    use dfam_curator::alignment::Orientation;
    use dfam_curator::kimura::kimura_pair;

    let reference = match msa.reference() {
        Some(r) => r,
        None => return Ok(()),
    };

    writeln!(out, "Alignment Stats:")?;
    writeln!(out, "\tReference Sequence\t\t\t\t\t\t\t\tConsensus Sequence")?;
    writeln!(out, "seq_id\tTranS\tTranSMod\tTranV\tKimura\tKimuraMod\tCpGSites\tWellCharBases\t\tTranS\tTranSMod\tTranV\tKimura\tKimuraMod\tCpGSites\tWellCharBases")?;

    let mut count = 0usize;
    let mut sum_r_kim = 0.0f64;
    let mut sum_r_kim_adj = 0.0f64;
    let mut sum_c_kim = 0.0f64;
    let mut sum_c_kim_adj = 0.0f64;
    let mut num_high = 0usize;

    for inst in &msa.sequences[1..] {
        let rs = kimura_pair(&reference.seq, &inst.seq);
        let cs = kimura_pair(consensus, &inst.seq);

        let rk  = if rs.kimura.is_finite()          { rs.kimura }          else { 100.0 };
        let rkm = if rs.kimura_adjusted.is_finite()  { rs.kimura_adjusted } else { 100.0 };
        let ck  = if cs.kimura.is_finite()           { cs.kimura }          else { 100.0 };
        let ckm = if cs.kimura_adjusted.is_finite()  { cs.kimura_adjusted } else { 100.0 };

        let is_high = rk >= 90.0;
        if is_high { num_high += 1; }

        let (coords, prefix) = if inst.orient == Orientation::Reverse {
            (format!("{}-{}", inst.seq_end, inst.seq_start), is_high)
        } else {
            (format!("{}-{}", inst.seq_start, inst.seq_end), is_high)
        };
        let label = if prefix {
            format!("**{}:{}", inst.name, coords)
        } else {
            format!("{}:{}", inst.name, coords)
        };

        writeln!(out,
            "{}\t{}\t{:.1}\t{}\t{:.2}\t{:.2}\t{}\t{}\t\t{}\t{:.1}\t{}\t{:.2}\t{:.2}\t{}\t{}",
            label,
            rs.transitions, rs.transitions_adjusted, rs.transversions, rk, rkm, rs.cpg_sites, rs.well_characterised,
            cs.transitions, cs.transitions_adjusted, cs.transversions, ck, ckm, cs.cpg_sites, cs.well_characterised,
        )?;

        sum_r_kim     += rk;
        sum_r_kim_adj += rkm;
        sum_c_kim     += ck;
        sum_c_kim_adj += ckm;
        count += 1;
    }

    writeln!(out)?;
    writeln!(out)?;

    if num_high > 0 {
        writeln!(out, "WARNING: There were {} alignments which had > 90% divergenge relative to the consensus!\n", num_high)?;
    }

    let avg_c_kim     = if count > 0 { sum_c_kim     / count as f64 } else { 0.0 };
    let avg_c_kim_adj = if count > 0 { sum_c_kim_adj / count as f64 } else { 0.0 };
    let avg_r_kim     = if count > 0 { sum_r_kim     / count as f64 } else { 0.0 };
    let avg_r_kim_adj = if count > 0 { sum_r_kim_adj / count as f64 } else { 0.0 };

    writeln!(out, "Total sequences: {}", count)?;
    writeln!(out, "Avg Kimura Div: {:.2}", avg_c_kim / 100.0)?;
    writeln!(out, "Avg Kimura Div (CpG adjusted): {:.2}", avg_c_kim_adj / 100.0)?;
    writeln!(out, "Relative to Reference Sequence (if available):")?;
    writeln!(out, "\tKimura Divergence:\t{:.2}", avg_r_kim)?;
    writeln!(out, "\tKimura Divergence (CpG Mod):\t{:.2}\n", avg_r_kim_adj)?;
    writeln!(out, "Relative to Consensus Sequence:")?;
    writeln!(out, "\tKimura Divergence:\t{:.2}", avg_c_kim)?;
    writeln!(out, "\tKimura Divergence (CpG Mod):\t{:.2}\n", avg_c_kim_adj)?;

    Ok(())
}
