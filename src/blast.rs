/// BLAST integration: external process wrappers, tabular hit parsing, and
/// conversion of pairwise hits into a MultiAlign.
///
/// The rmblastn search uses a custom outfmt-6 field list that includes the
/// gapped alignment strings (`qseq`/`sseq`), which are needed to place each
/// instance into the shared reference coordinate space.
use std::io::Write;
use std::path::Path;
use std::process::Command;

use anyhow::{bail, Context};

use crate::alignment::{MultiAlign, Orientation, SequenceRow};

// ── Search parameters ─────────────────────────────────────────────────────────

/// Parameters forwarded to makeblastdb and rmblastn.
#[derive(Debug, Clone)]
pub struct BlastParams {
    /// Path (or name on $PATH) of the makeblastdb executable.
    pub makeblastdb: String,
    /// Path (or name on $PATH) of the rmblastn executable.
    pub rmblastn: String,
    /// Optional substitution matrix file (rmblastn -matrix).
    pub matrix: Option<std::path::PathBuf>,
    /// Minimum raw gapped score for a hit to be reported (default 150).
    pub min_score: u32,
    /// rmblastn -gapopen value (default 20; derived from Perl gap_init=-25, gap_ext=-5).
    pub gap_open: u32,
    /// rmblastn -gapextend value (default 5).
    pub gap_extend: u32,
    /// rmblastn -word_size (default 7).
    pub word_size: u32,
    /// rmblastn -num_threads (default 4).
    pub num_threads: u32,
    /// rmblastn -mask_level (default 80; rmblastn-specific flag).
    pub mask_level: u32,
}

impl Default for BlastParams {
    fn default() -> Self {
        BlastParams {
            makeblastdb: "makeblastdb".to_string(),
            rmblastn: "rmblastn".to_string(),
            matrix: None,
            min_score: 150,
            gap_open: 20,
            gap_extend: 5,
            word_size: 7,
            num_threads: 4,
            mask_level: 80,
        }
    }
}

// ── Hit ───────────────────────────────────────────────────────────────────────

/// One pairwise alignment record parsed from rmblastn tabular output.
#[derive(Debug, Clone)]
pub struct BlastHit {
    /// Raw alignment score.
    pub score: u32,
    /// Query sequence name.
    pub query_name: String,
    /// Query start (1-based, closed).
    pub query_start: u64,
    /// Query end (1-based, closed).
    pub query_end: u64,
    /// Total query sequence length.
    pub query_len: u64,
    /// Subject sequence name.
    pub subj_name: String,
    /// Subject start (1-based; always <= subj_end in BLAST tabular output).
    pub subj_start: u64,
    /// Subject end (1-based).
    pub subj_end: u64,
    /// Total subject sequence length.
    pub subj_len: u64,
    /// Forward = plus strand of subject; Reverse = minus strand.
    pub orientation: Orientation,
    /// Gapped query sequence from the alignment (gaps = b'-').
    pub query_seq: Vec<u8>,
    /// Gapped subject sequence from the alignment (gaps = b'-').
    pub subj_seq: Vec<u8>,
}

// ── External process wrappers ────────────────────────────────────────────────

/// Build a nucleotide BLAST database from `fasta`, storing files at `db_prefix`.
pub fn make_blast_db(params: &BlastParams, fasta: &Path, db_prefix: &Path) -> anyhow::Result<()> {
    let status = Command::new(&params.makeblastdb)
        .args([
            "-blastdb_version", "4",
            "-dbtype", "nucl",
            "-in",  fasta.to_str().context("non-UTF8 path in makeblastdb -in")?,
            "-out", db_prefix.to_str().context("non-UTF8 path in makeblastdb -out")?,
        ])
        .stderr(std::process::Stdio::null())
        .status()
        .with_context(|| format!("failed to run {}", params.makeblastdb))?;

    if !status.success() {
        bail!("makeblastdb exited with {}", status);
    }
    Ok(())
}

/// Run rmblastn (query vs db) and return parsed hits.
///
/// Output format:
/// `6 score qseqid qstart qend qlen sstrand sseqid sstart send slen qseq sseq`
pub fn run_rmblastn(
    params: &BlastParams,
    query: &Path,
    db: &Path,
) -> anyhow::Result<Vec<BlastHit>> {
    const OUTFMT: &str =
        "6 score qseqid qstart qend qlen sstrand sseqid sstart send slen qseq sseq";

    let xdrop_ungap  = params.min_score * 2;
    let xdrop_gap    = params.min_score / 2;
    let xdrop_final  = params.min_score;

    let mut cmd = Command::new(&params.rmblastn);
    cmd.args([
        "-query",               query.to_str().context("non-UTF8 query path")?,
        "-db",                  db.to_str().context("non-UTF8 db path")?,
        "-outfmt",              OUTFMT,
        "-gapopen",             &params.gap_open.to_string(),
        "-gapextend",           &params.gap_extend.to_string(),
        "-word_size",           &params.word_size.to_string(),
        "-xdrop_ungap",         &xdrop_ungap.to_string(),
        "-xdrop_gap",           &xdrop_gap.to_string(),
        "-xdrop_gap_final",     &xdrop_final.to_string(),
        "-min_raw_gapped_score",&params.min_score.to_string(),
        "-dust",                "no",
        "-num_threads",         &params.num_threads.to_string(),
        "-mask_level",          &params.mask_level.to_string(),
        "-complexity_adjust",
    ]);

    if let Some(m) = &params.matrix {
        cmd.args(["-matrix", m.to_str().context("non-UTF8 matrix path")?]);
    }

    let output = cmd
        .output()
        .with_context(|| format!("failed to run {}", params.rmblastn))?;

    if !output.status.success() {
        bail!(
            "rmblastn failed ({}): {}",
            output.status,
            String::from_utf8_lossy(&output.stderr).trim()
        );
    }

    let text = String::from_utf8_lossy(&output.stdout);
    parse_hits(&text)
}

/// Parse BLAST tabular output produced with the custom field list:
/// `6 score qseqid qstart qend qlen sstrand sseqid sstart send slen qseq sseq`
///
/// This is the public API for building a `MultiAlign` from saved tabular output
/// without going through `run_rmblastn()`.
pub fn parse_tab(text: &str) -> anyhow::Result<Vec<BlastHit>> {
    parse_hits(text)
}

fn parse_hits(text: &str) -> anyhow::Result<Vec<BlastHit>> {
    let mut hits = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.splitn(13, '\t').collect();
        if f.len() < 12 {
            continue;
        }
        // score qseqid qstart qend qlen sstrand sseqid sstart send slen qseq sseq
        let score: u32       = f[0].parse().context("blast score")?;
        let query_name       = f[1].to_string();
        let query_start: u64 = f[2].parse().context("qstart")?;
        let query_end: u64   = f[3].parse().context("qend")?;
        let query_len: u64   = f[4].parse().context("qlen")?;
        let orientation      = if f[5] == "minus" { Orientation::Reverse } else { Orientation::Forward };
        let subj_name        = f[6].to_string();
        let subj_start: u64  = f[7].parse().context("sstart")?;
        let subj_end: u64    = f[8].parse().context("send")?;
        let subj_len: u64    = f[9].parse().context("slen")?;
        let query_seq        = f[10].as_bytes().to_vec();
        let subj_seq         = f[11].as_bytes().to_vec();

        hits.push(BlastHit {
            score,
            query_name,
            query_start,
            query_end,
            query_len,
            subj_name,
            subj_start,
            subj_end,
            subj_len,
            orientation,
            query_seq,
            subj_seq,
        });
    }
    Ok(hits)
}

// ── MSA builder ───────────────────────────────────────────────────────────────

/// Build a `MultiAlign` from pairwise BLAST hits against a single reference.
///
/// `ref_seq` is the **ungapped** consensus/reference that was used as the BLAST
/// subject.  The MSA width equals `ref_seq.len()`.  Insertions relative to the
/// reference (gaps in `subj_seq`) are silently dropped so every row has the
/// same width.  `ref_seq` is stored at index 0; each hit becomes one instance
/// row.
///
/// For reverse-strand hits the gapped alignment strings are reverse-
/// complemented before mapping so that all rows are in the same strand
/// orientation as the reference.
pub fn hits_to_multialign(ref_seq: &[u8], ref_name: &str, hits: &[BlastHit]) -> MultiAlign {
    let width = ref_seq.len();
    let reference = SequenceRow::new(ref_name, ref_seq.to_vec());

    let instances: Vec<SequenceRow> = hits
        .iter()
        .map(|hit| {
            let mut row = vec![b' '; width];
            let ref_start = (hit.subj_start as usize).saturating_sub(1); // 0-based

            // For a minus-strand hit, BLAST reports:
            //   sseq  = RC of the plus-strand subject (read 5'→3' on the minus strand)
            //   qseq  = query aligned against that RC subject
            //   sstart/send are still plus-strand coords (sstart <= send)
            //
            // To place in the forward MSA we RC both strings; after reversal the
            // pair reads left-to-right in plus-strand order starting at sstart.
            let (sseq, qseq): (Vec<u8>, Vec<u8>) = if hit.orientation == Orientation::Reverse {
                let s: Vec<u8> = hit.subj_seq.iter().rev().map(|&b| iupac_complement(b)).collect();
                let q: Vec<u8> = hit.query_seq.iter().rev().map(|&b| iupac_complement(b)).collect();
                (s, q)
            } else {
                (hit.subj_seq.clone(), hit.query_seq.clone())
            };

            let mut ref_col = ref_start;
            for (&s, &q) in sseq.iter().zip(qseq.iter()) {
                if ref_col >= width {
                    break;
                }
                if s == b'-' {
                    // Insertion in query relative to reference — drop.
                    continue;
                }
                row[ref_col] = q;
                ref_col += 1;
            }

            let mut seq_row = SequenceRow::new(hit.query_name.clone(), row);
            seq_row.seq_start = hit.query_start;
            seq_row.seq_end   = hit.query_end;
            seq_row.orient    = hit.orientation;
            seq_row
        })
        .collect();

    MultiAlign::from_sequences(reference, instances)
}

/// Write a single named sequence to a FASTA file (utility used by the refiner).
pub fn write_fasta(path: &Path, name: &str, seq: &[u8]) -> anyhow::Result<()> {
    let mut f = std::io::BufWriter::new(
        std::fs::File::create(path)
            .with_context(|| format!("cannot create {}", path.display()))?,
    );
    writeln!(f, ">{}", name)?;
    f.write_all(seq)?;
    writeln!(f)?;
    Ok(())
}

// ── Internal helpers ──────────────────────────────────────────────────────────

fn iupac_complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T', b'T' => b'A', b'G' => b'C', b'C' => b'G',
        b'R' => b'Y', b'Y' => b'R', b'K' => b'M', b'M' => b'K',
        b'S' => b'S', b'W' => b'W', b'B' => b'V', b'V' => b'B',
        b'D' => b'H', b'H' => b'D', b'N' => b'N', b'X' => b'X',
        b'-' => b'-', b' ' => b' ',
        _ => b'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_forward_hit() {
        let line = "500\tseq1\t10\t200\t300\tplus\tcons\t1\t190\t800\tACGT-ACG\tACGTTACG";
        let hits = parse_hits(line).unwrap();
        assert_eq!(hits.len(), 1);
        let h = &hits[0];
        assert_eq!(h.score, 500);
        assert_eq!(h.query_name, "seq1");
        assert_eq!(h.orientation, Orientation::Forward);
        assert_eq!(h.subj_start, 1);
        assert_eq!(&h.query_seq, b"ACGT-ACG");
    }

    #[test]
    fn parse_reverse_hit() {
        let line = "300\tseq2\t50\t100\t200\tminus\tcons\t10\t60\t800\tACGT\tACGT";
        let hits = parse_hits(line).unwrap();
        assert_eq!(hits[0].orientation, Orientation::Reverse);
    }

    #[test]
    fn hits_to_msa_basic() {
        // Reference: ACGT (width 4)
        // One hit covering the full reference, no gaps.
        let ref_seq = b"ACGT";
        let hit = BlastHit {
            score: 100,
            query_name: "q1".to_string(),
            query_start: 1, query_end: 4, query_len: 4,
            subj_name: "ref".to_string(),
            subj_start: 1, subj_end: 4, subj_len: 4,
            orientation: Orientation::Forward,
            query_seq: b"ACGT".to_vec(),
            subj_seq:  b"ACGT".to_vec(),
        };
        let msa = hits_to_multialign(ref_seq, "ref", &[hit]);
        assert_eq!(msa.num_instances(), 1);
        assert_eq!(msa.instance(0).unwrap().seq, b"ACGT");
    }

    #[test]
    fn hits_to_msa_partial_and_deletion() {
        // Reference: ACGTACGT (width 8)
        // Hit covers positions 3-6 (1-based), with a deletion at pos 5.
        //   subj_seq: CGTA  (no gaps — all 4 ref positions covered)
        //   query_seq: CG-A  (gap at position 5 → deletion in query)
        let ref_seq = b"ACGTACGT";
        let hit = BlastHit {
            score: 80,
            query_name: "q1".to_string(),
            query_start: 1, query_end: 3, query_len: 3,
            subj_name: "ref".to_string(),
            subj_start: 3, subj_end: 6, subj_len: 8,
            orientation: Orientation::Forward,
            query_seq: b"CG-A".to_vec(),
            subj_seq:  b"CGTA".to_vec(),
        };
        let msa = hits_to_multialign(ref_seq, "ref", &[hit]);
        let row = &msa.instance(0).unwrap().seq;
        // Columns 0-1: space padding; 2: C; 3: G; 4: -; 5: A; 6-7: space padding
        assert_eq!(row[0], b' ');
        assert_eq!(row[2], b'C');
        assert_eq!(row[3], b'G');
        assert_eq!(row[4], b'-');
        assert_eq!(row[5], b'A');
        assert_eq!(row[7], b' ');
    }

    #[test]
    fn hits_to_msa_insertion_dropped() {
        // An insertion in the query (gap in subj_seq) should not expand the MSA.
        // subj_seq: AC-GT  (gap = insertion in query)
        // query_seq: ACXGT (X is the inserted base)
        // Only 4 ref positions covered; X is dropped.
        let ref_seq = b"ACGT";
        let hit = BlastHit {
            score: 60,
            query_name: "q1".to_string(),
            query_start: 1, query_end: 5, query_len: 5,
            subj_name: "ref".to_string(),
            subj_start: 1, subj_end: 4, subj_len: 4,
            orientation: Orientation::Forward,
            query_seq: b"ACXGT".to_vec(),
            subj_seq:  b"AC-GT".to_vec(),
        };
        let msa = hits_to_multialign(ref_seq, "ref", &[hit]);
        let row = &msa.instance(0).unwrap().seq;
        assert_eq!(row.len(), 4);
        assert_eq!(row, b"ACGT"); // X dropped; A,C,G,T placed at cols 0-3
    }
}
