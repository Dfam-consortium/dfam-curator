/// Parser for RepeatMasker / crossmatch `.align` pairwise alignment files.
///
/// Format overview (from CrossmatchSearchEngine.pm documentation):
///
/// ```text
/// <score> <pctDiv> <pctDel> <pctIns> <queryName> <queryStart> <queryEnd> (<queryLeft>)
///     [C] <subjName> [<subjType>] [(<subjLeft>)] <subjEnd> <subjStart>  [id] [overlap]
///
/// [C] <queryName> <qStart> <gappedQuerySeq> <qEnd>
///          <matchString>
/// [C] <subjName>  <sStart> <gappedSubjSeq>  <sEnd>
///
/// Gap_init rate = ...
/// ```
///
/// The `C` prefix on a line indicates the complement (reverse-complement) strand.
/// Fields marked `[...]` are optional depending on the variant (`.out`, `.align`, standard).
///
/// This parser reads the **`.align`** variant produced by RepeatMasker, which always
/// has the gapped alignment strings interleaved after the score header.
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use crate::alignment::Orientation;

/// One pairwise alignment record from a crossmatch `.align` file.
#[derive(Debug, Clone)]
pub struct PairwiseHit {
    /// Smith-Waterman score.
    pub sw_score: u32,
    /// Percent divergence.
    pub pct_div: f64,
    /// Percent deletion (gaps in subject relative to query).
    pub pct_del: f64,
    /// Percent insertion (gaps in query relative to subject).
    pub pct_ins: f64,

    /// Query sequence name.
    pub query_name: String,
    /// Query start (1-based, fully-closed).
    pub query_start: u64,
    /// Query end (1-based, fully-closed).
    pub query_end: u64,
    /// Bases remaining in query after query_end.
    pub query_remaining: u64,
    /// Gapped query sequence (may be empty if file has no alignment strings).
    pub query_seq: Vec<u8>,

    /// Subject sequence name.
    pub subj_name: String,
    /// Subject start (1-based, fully-closed; always ≤ subj_end regardless of strand).
    pub subj_start: u64,
    /// Subject end (1-based, fully-closed).
    pub subj_end: u64,
    /// Bases remaining in subject.
    pub subj_remaining: u64,
    /// Gapped subject sequence.
    pub subj_seq: Vec<u8>,

    /// Strand of the **query** relative to the subject.
    /// `Forward` = same strand; `Reverse` = reverse complement.
    pub orientation: Orientation,

    /// Optional unique alignment identifier field.
    pub id: Option<String>,
}

/// Read all pairwise hits from a crossmatch `.align` file.
pub fn read(path: &Path) -> io::Result<Vec<PairwiseHit>> {
    let f = BufReader::new(std::fs::File::open(path)?);
    parse(f)
}

/// Parse crossmatch alignment records from any `BufRead` source.
pub fn parse<R: BufRead>(reader: R) -> io::Result<Vec<PairwiseHit>> {
    let mut hits: Vec<PairwiseHit> = Vec::new();
    let mut current: Option<PairwiseHit> = None;
    // 0 = expecting header or blank; 1 = query seq line; 2 = match line; 3 = subj seq line
    let mut align_pos: u8 = 0;

    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();

        if trimmed.is_empty() {
            // A blank line at align_pos==2 IS the match-indicator line (happens when
            // query and subject are identical — all spaces, no v/i characters).
            // Advance through it; ignore blank lines everywhere else.
            if align_pos == 2 {
                align_pos = 3;
            }
            continue;
        }

        // ── Record terminator ────────────────────────────────────────────────
        // Standalone:  "Gap_init rate = ..."
        // Combined:    "Transitions / transversions = ...; Gap_init rate = ..."
        if trimmed.starts_with("Gap_init") || trimmed.starts_with("Transitions /") {
            if let Some(hit) = current.take() {
                hits.push(hit);
            }
            align_pos = 0;
            continue;
        }

        // ── Score header detection ────────────────────────────────────────────
        // Starts with an integer score (possibly with leading spaces).
        // Pattern: <int> <float> <float> <float> <name> <int> <int> (<int>) ...
        if is_score_line(trimmed) {
            // Push any pending hit.
            if let Some(hit) = current.take() {
                hits.push(hit);
            }
            current = parse_score_line(trimmed)?;
            align_pos = 1; // next non-blank line should be query alignment
            continue;
        }

        // ── Alignment sequence lines ──────────────────────────────────────────
        if let Some(ref mut hit) = current {
            match align_pos {
                1 => {
                    // Query sequence line — append (multi-block alignments concatenate).
                    if let Some(seq) = parse_align_seq_line(trimmed) {
                        hit.query_seq.extend_from_slice(&seq);
                    }
                    align_pos = 2;
                }
                2 => {
                    // Match-indicator line (v/i/space characters) — discard content,
                    // just advance state.  Blank indicator lines are handled above.
                    align_pos = 3;
                }
                3 => {
                    // Subject sequence line — append.
                    if let Some(seq) = parse_align_seq_line(trimmed) {
                        hit.subj_seq.extend_from_slice(&seq);
                    }
                    align_pos = 1; // cycle back for next block
                }
                _ => {}
            }
        }
    }

    // Flush last record.
    if let Some(hit) = current.take() {
        hits.push(hit);
    }

    Ok(hits)
}

// ── Internal helpers ──────────────────────────────────────────────────────────

/// True if `line` looks like a crossmatch score header.
///
/// Requires the first four whitespace-delimited tokens to be numeric:
/// an integer SW score followed by three floating-point percent values.
/// This rejects preamble lines like "4 distinct alphabetic chars..." that
/// happen to start with a digit.
fn is_score_line(line: &str) -> bool {
    let mut fields = line.split_whitespace();
    let f0 = fields.next().unwrap_or("");
    if f0.is_empty() || !f0.chars().all(|c| c.is_ascii_digit()) {
        return false;
    }
    // The next three fields must parse as floats (pct_div, pct_del, pct_ins).
    for _ in 0..3 {
        if !fields.next().map(|f| f.parse::<f64>().is_ok()).unwrap_or(false) {
            return false;
        }
    }
    true
}

/// Parse a score header line into a PairwiseHit (without alignment sequences).
///
/// Two field layouts are supported:
///
/// **Forward strand** (11–12 fields):
/// ```text
/// score div del ins qName qStart qEnd (qLeft) subjName sStart sEnd (sLeft) [id]
/// ```
///
/// **Reverse-complement strand** (12–13 fields):
/// ```text
/// score div del ins qName qStart qEnd (qLeft) C subjName (sLeft) sEnd sStart [id]
/// ```
fn parse_score_line(line: &str) -> io::Result<Option<PairwiseHit>> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < 9 {
        return Ok(None); // Not enough fields — skip (preamble / header comment).
    }

    let sw_score = parse_u32(fields[0])?;
    let pct_div  = parse_f64(fields[1])?;
    let pct_del  = parse_f64(fields[2])?;
    let pct_ins  = parse_f64(fields[3])?;

    let query_name  = fields[4].to_string();
    let query_start = parse_u64(fields[5])?;
    let query_end   = parse_u64(fields[6])?;
    let query_remaining = parse_paren_u64(fields[7])?;

    // Field 8 is either "C" (complement flag) or the subject name.
    let (orientation, subj_field_offset) = if fields[8] == "C" {
        (Orientation::Reverse, 1)
    } else {
        (Orientation::Forward, 0)
    };

    let subj_name_field = 8 + subj_field_offset;
    if subj_name_field >= fields.len() {
        return Ok(None);
    }
    let subj_name = fields[subj_name_field].to_string();

    // Skip optional subject "type" field (e.g. "SINE/Alu") if present.
    // Heuristic: the next field that is NOT in parentheses and IS numeric is subj_start.
    let mut fi = subj_name_field + 1;

    // The remaining/start/end triplet order depends on orientation:
    //   Forward: sStart sEnd (sLeft)
    //   Reverse: (sLeft) sEnd sStart
    let (subj_start, subj_end, subj_remaining);

    if orientation == Orientation::Reverse {
        // (sLeft) sEnd sStart
        // Skip any non-parenthesised, non-numeric fields (subject type).
        while fi < fields.len() && !fields[fi].starts_with('(') && parse_u64(fields[fi]).is_err() {
            fi += 1;
        }
        subj_remaining = if fi < fields.len() { parse_paren_u64(fields[fi]).unwrap_or(0) } else { 0 };
        fi += 1;
        subj_end   = if fi < fields.len() { parse_u64(fields[fi]).unwrap_or(0) } else { 0 };
        fi += 1;
        subj_start = if fi < fields.len() { parse_u64(fields[fi]).unwrap_or(0) } else { 0 };
        fi += 1;
    } else {
        // Skip optional subject type fields.
        while fi < fields.len() && !fields[fi].starts_with('(') && parse_u64(fields[fi]).is_err() {
            fi += 1;
        }
        subj_start = if fi < fields.len() { parse_u64(fields[fi]).unwrap_or(0) } else { 0 };
        fi += 1;
        subj_end   = if fi < fields.len() { parse_u64(fields[fi]).unwrap_or(0) } else { 0 };
        fi += 1;
        subj_remaining = if fi < fields.len() { parse_paren_u64(fields[fi]).unwrap_or(0) } else { 0 };
        fi += 1;
    }

    // Optional ID field.
    let id = fields.get(fi).map(|s| s.to_string());

    Ok(Some(PairwiseHit {
        sw_score,
        pct_div,
        pct_del,
        pct_ins,
        query_name,
        query_start,
        query_end,
        query_remaining,
        query_seq: Vec::new(),
        subj_name,
        subj_start,
        subj_end,
        subj_remaining,
        subj_seq: Vec::new(),
        orientation,
        id,
    }))
}

/// Parse an alignment sequence line of the form:
/// `[C] <name>  <pos1>  <gappedSeq>  <pos2>`
///
/// The gapped sequence is the third whitespace-delimited field (or second if no
/// "C" prefix), surrounded by position numbers.
fn parse_align_seq_line(line: &str) -> Option<Vec<u8>> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    // Layout: [C] name pos1 GAPPEDSEQ pos2
    // With C: fields = ["C", name, pos1, seq, pos2] -> seq at index 3
    // Without: fields = [name, pos1, seq, pos2] -> seq at index 2
    let seq_idx = if fields.first().map(|f| *f == "C").unwrap_or(false) {
        3
    } else {
        2
    };
    fields.get(seq_idx).map(|s| s.as_bytes().to_vec())
}

// ── Number parsing helpers ────────────────────────────────────────────────────

fn parse_u32(s: &str) -> io::Result<u32> {
    s.parse::<u32>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn parse_u64(s: &str) -> io::Result<u64> {
    s.parse::<u64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn parse_f64(s: &str) -> io::Result<f64> {
    s.parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Parse a parenthesised integer like `(1234)`.
fn parse_paren_u64(s: &str) -> io::Result<u64> {
    let inner = s.trim_matches(|c| c == '(' || c == ')');
    inner.parse::<u64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const SAMPLE: &str = "\
2334 8.44 0.00 3.25 Human 127 737 (8222) AluSx#SINE/Alu 1 298 (14)

  Human       127 TAAAGTCCCTGCTCGCCCCCGCTCCAGCTCGCC 159
                   vvvi  vv  v v   v   v      v v  v
  AluSx#SINE  1   TAAGATCCCTGCTCACCCCCGCTCCAGGTCACC 33

Gap_init rate = 25.0

254 28.0 4.0 9.0 seq-13 3873751 3873963 (114) C L2d 3331 3075 (6)

C seq-13   3873963 CTGCAACATGCGACACAACACGCGT 3939
                    v  i  ii   v ii   i   i
  L2d      3331    CTACAACATACGCCACAACACGCGT 3305

Gap_init rate = 22.0
";

    #[test]
    fn parse_two_records() {
        let hits = parse(Cursor::new(SAMPLE)).unwrap();
        assert_eq!(hits.len(), 2);

        let h0 = &hits[0];
        assert_eq!(h0.sw_score, 2334);
        assert_eq!(h0.query_name, "Human");
        assert_eq!(h0.query_start, 127);
        assert_eq!(h0.query_end, 737);
        assert_eq!(h0.orientation, Orientation::Forward);
        assert!(!h0.query_seq.is_empty());
        assert!(!h0.subj_seq.is_empty());

        let h1 = &hits[1];
        assert_eq!(h1.sw_score, 254);
        assert_eq!(h1.query_name, "seq-13");
        assert_eq!(h1.orientation, Orientation::Reverse);
    }

    #[test]
    fn is_score_line_check() {
        assert!(is_score_line("2334 8.44 0.00 3.25 Human 127 737 (8222) AluSx 1 298 (14)"));
        assert!(!is_score_line("  Human       127 TAAAGT 159"));
        assert!(!is_score_line("Gap_init rate = 25.0"));
    }

    /// RepeatMasker produces a combined "Transitions / transversions = ...; Gap_init rate = ..."
    /// terminator instead of a standalone "Gap_init" line.  Verify both sequences are clean.
    #[test]
    fn combined_transitions_terminator() {
        const COMBINED: &str = "\
2334 8.44 0.00 3.25 Human 127 159 (8222) AluSx 1 33 (14)

  Human       127 TAAAGTCCCTGCTCGCCCCCGCTCCAGCTCGCC 159
                   vvvi  vv  v v   v   v      v v  v
  AluSx       1   TAAGATCCCTGCTCACCCCCGCTCCAGGTCACC 33

Transitions / transversions = 2.17 (26 / 12); Gap_init rate = 0.00 (2 / 421), avg. gap size = 2.50 (5 / 2)

254 28.0 4.0 9.0 seq-13 3873751 3873963 (114) C L2d 3331 3075 (6)

C seq-13   3873963 CTGCAACATGCGACACAACACGCGT 3939
                    v  i  ii   v ii   i   i
  L2d      3331    CTACAACATACGCCACAACACGCGT 3305

Transitions / transversions = 1.00 (4 / 4); Gap_init rate = 0.00 (0 / 24), avg. gap size = 0.00 (0 / 0)
";
        let hits = parse(Cursor::new(COMBINED)).unwrap();
        assert_eq!(hits.len(), 2);
        // Sequences must not contain "transversions", "entries", or other stats text.
        let q0 = std::str::from_utf8(&hits[0].query_seq).unwrap();
        let s0 = std::str::from_utf8(&hits[0].subj_seq).unwrap();
        assert!(!q0.contains("transversions"), "query_seq contaminated: {q0}");
        assert!(!s0.contains("transversions"), "subj_seq contaminated: {s0}");
        assert_eq!(hits[1].orientation, Orientation::Reverse);
    }
}
