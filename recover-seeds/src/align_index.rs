//! Index a RepeatMasker `.align` file by family.
//!
//! RepeatMasker orders the file by genome position, so one family's copies are
//! scattered through it.  A single scan records where each record starts and
//! ends, so a later pass can read one family on its own, in any order and in
//! parallel.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use anyhow::{Context, Result};

/// Where one alignment record sits in the file, and what it says in its header.
#[derive(Debug, Clone)]
pub struct RecordRef {
    pub offset: u64,
    pub len: usize,
    /// Percent divergence from the header line.
    pub pct_div: f64,
    /// Kimura divergence with the CpG modification, from the record's trailer.
    /// generateSeedAlignments.pl ranks instances by this, not by pct_div.
    pub kimura: f64,
    /// Genome span, 1-based fully closed.
    pub query_start: u64,
    pub query_end: u64,
    /// Consensus span, 1-based fully closed, low end first.
    pub cons_start: u64,
    pub cons_end: u64,
    /// Consensus length as RepeatMasker saw it (span plus what it left over).
    pub cons_len: u64,
    /// Genome sequence name.
    pub seq_name: String,
    /// RepeatMasker's element id: fragments of one copy share it.
    pub element_id: u32,
}

#[derive(Debug, Default)]
pub struct AlignIndex {
    pub families: HashMap<String, Vec<RecordRef>>,
    pub records: usize,
    /// Records skipped as simple, low-complexity or short, as
    /// generateSeedAlignments.pl skips them.  RepeatMasker's own tandem-repeat
    /// annotations arrive this way, carrying a per-locus consensus that
    /// matches no library entry.
    pub skipped_simple: usize,
    pub unparsed_headers: usize,
}

/// The families generateSeedAlignments.pl drops on sight.
fn is_simple_or_low(header: &str) -> bool {
    header.contains("#Simple") || header.contains("#Low") || header.contains("short")
}

/// True for a record header: an integer score then three percentages.
fn is_header(line: &str) -> bool {
    let mut f = line.split_whitespace();
    match f.next() {
        Some(t) if !t.is_empty() && t.bytes().all(|b| b.is_ascii_digit()) => {}
        _ => return false,
    }
    (0..3).all(|_| f.next().is_some_and(|t| t.parse::<f64>().is_ok()))
}

/// Family name, consensus span and element id from a header line.
fn parse_header(line: &str) -> Option<(String, f64, u64, u64, u64, String, u32, u64, u64)> {
    let f: Vec<&str> = line.split_whitespace().collect();
    if f.len() < 12 {
        return None;
    }
    let pct_div: f64 = f[1].parse().ok()?;
    let seq_name = f[4].to_string();
    let qs: u64 = f[5].parse().ok()?;
    let qe: u64 = f[6].parse().ok()?;
    let paren = |s: &str| -> Option<u64> { s.trim_matches(['(', ')']).parse().ok() };
    let (fam, lo, hi, left, rest) = if f[8] == "C" {
        // score div del ins qName qStart qEnd (qLeft) C subjName (sLeft) sEnd sStart [id]
        let left = paren(f[10])?;
        let a: u64 = f[11].parse().ok()?;
        let b: u64 = f[12].parse().ok()?;
        (f[9], b.min(a), a.max(b), left, &f[13..])
    } else {
        // score div del ins qName qStart qEnd (qLeft) subjName sStart sEnd (sLeft) [id]
        let a: u64 = f[9].parse().ok()?;
        let b: u64 = f[10].parse().ok()?;
        (f[8], a.min(b), a.max(b), paren(f[11])?, &f[12..])
    };
    // The element id is the last integer field on the line.
    let element_id = rest
        .iter()
        .rev()
        .find_map(|t| t.parse::<u32>().ok())
        .unwrap_or(0);
    let fam = fam.split('#').next().unwrap_or(fam).to_string();
    Some((fam, pct_div, lo, hi, hi + left, seq_name, element_id, qs, qe))
}

/// Scan the file once, recording each record's extent and header facts.
pub fn build(path: &Path) -> Result<AlignIndex> {
    let file = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let mut reader = BufReader::with_capacity(1 << 22, file);
    let mut idx = AlignIndex::default();
    let mut line = String::new();
    let mut offset: u64 = 0;
    let mut open: Option<(String, RecordRef)> = None;

    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        let start = offset;
        offset += n as u64;
        // The Kimura line belongs to the record still open.
        if let Some(rest) = line.trim_start().strip_prefix("Kimura (with divCpGMod) = ") {
            if let (Some((_, rec)), Ok(v)) = (open.as_mut(), rest.trim().parse::<f64>()) {
                rec.kimura = v;
            }
            continue;
        }
        if !line.as_bytes().first().is_some_and(u8::is_ascii_digit) || !is_header(line.trim_end()) {
            continue;
        }
        if let Some((fam, mut rec)) = open.take() {
            rec.len = (start - rec.offset) as usize;
            idx.families.entry(fam).or_default().push(rec);
        }
        if is_simple_or_low(&line) {
            idx.skipped_simple += 1;
            continue;
        }
        match parse_header(line.trim_end()) {
            Some((fam, pct_div, lo, hi, clen, seq_name, element_id, qs, qe)) => {
                idx.records += 1;
                open = Some((
                    fam,
                    RecordRef {
                        offset: start,
                        len: 0,
                        pct_div,
                        kimura: pct_div,     // replaced when the trailer is read
                        query_start: qs,
                        query_end: qe,
                        cons_start: lo,
                        cons_end: hi,
                        cons_len: clen,
                        seq_name,
                        element_id,
                    },
                ));
            }
            None => idx.unparsed_headers += 1,
        }
    }
    if let Some((fam, mut rec)) = open.take() {
        rec.len = (offset - rec.offset) as usize;
        idx.families.entry(fam).or_default().push(rec);
    }
    Ok(idx)
}

/// Read the given records back as text, in file order.
pub fn read_records(path: &Path, recs: &[RecordRef]) -> Result<String> {
    let mut file = File::open(path)?;
    let mut order: Vec<&RecordRef> = recs.iter().collect();
    order.sort_by_key(|r| r.offset);
    let mut out = String::new();
    let mut buf = Vec::new();
    for r in order {
        file.seek(SeekFrom::Start(r.offset))?;
        buf.resize(r.len, 0);
        file.read_exact(&mut buf)?;
        out.push_str(std::str::from_utf8(&buf)?);
    }
    Ok(out)
}
