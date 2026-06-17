/// Read and write aligned FASTA / A2M format.
///
/// In A2M format every sequence has the same length when gaps are included.
/// The first sequence is treated as the reference.
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

use crate::alignment::{Orientation, SequenceRow, MultiAlign};

/// Read an aligned FASTA / A2M file into a MultiAlign.
///
/// The first sequence becomes the reference; all others are instances.
/// Lines beginning with `>` are header lines; everything else is sequence data.
pub fn read(path: &Path) -> io::Result<MultiAlign> {
    let f = BufReader::new(std::fs::File::open(path)?);
    let mut entries: Vec<(String, Vec<u8>)> = Vec::new();
    let mut current_name = String::new();
    let mut current_seq: Vec<u8> = Vec::new();

    for line in f.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('>') {
            if !current_name.is_empty() || !current_seq.is_empty() {
                entries.push((current_name.clone(), current_seq.clone()));
                current_seq.clear();
            }
            current_name = rest.trim().to_string();
        } else {
            current_seq.extend_from_slice(trimmed.as_bytes());
        }
    }
    if !current_name.is_empty() || !current_seq.is_empty() {
        entries.push((current_name, current_seq));
    }

    if entries.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "empty FASTA file"));
    }

    let (ref_name, ref_seq) = entries.remove(0);
    let reference = SequenceRow::new(ref_name, ref_seq);
    let instances = entries
        .into_iter()
        .map(|(name, seq)| make_instance_row(name, seq))
        .collect();

    Ok(MultiAlign::from_sequences(reference, instances))
}

/// Write a MultiAlign as aligned FASTA (A2M).
///
/// Only instance sequences (not the reference) are written, matching
/// Perl `MultAln::toFASTA()` default behaviour.  Leading/trailing space
/// padding is converted to `-`.  The optional `consensus_seq` is
/// prepended as `>consensus` if provided.  The optional `description` is
/// appended to each sequence header after two spaces (FASTA description
/// field), e.g. `"families.stk:3 id=MamSINE1"`.
pub fn write(
    msa: &MultiAlign,
    out: &mut dyn Write,
    consensus_seq: Option<&[u8]>,
    description: Option<&str>,
) -> io::Result<()> {
    let desc_suffix = description.map(|d| format!("  {}", d)).unwrap_or_default();
    if let Some(cons) = consensus_seq {
        writeln!(out, ">consensus{}", desc_suffix)?;
        let s: Vec<u8> = cons.iter().map(|&b| if b == b' ' { b'-' } else { b }).collect();
        out.write_all(&s)?;
        writeln!(out)?;
    }
    for seq in &msa.sequences[1..] {
        writeln!(out, ">{}{}", seq_label(&seq.name, seq.seq_start, seq.seq_end, seq.orient), desc_suffix)?;
        let s: Vec<u8> = seq.seq.iter().map(|&b| if b == b' ' { b'-' } else { b }).collect();
        out.write_all(&s)?;
        writeln!(out)?;
    }
    Ok(())
}

/// Write a MultiAlign as unaligned FASTA (gaps stripped).
///
/// Only instance sequences are written, matching Perl `toFASTA(seqOnly => 1)`.
pub fn write_ungapped(
    msa: &MultiAlign,
    out: &mut dyn Write,
    consensus_seq: Option<&[u8]>,
) -> io::Result<()> {
    if let Some(cons) = consensus_seq {
        let ungapped: Vec<u8> = cons.iter().filter(|&&b| b != b'-').copied().collect();
        writeln!(out, ">consensus")?;
        out.write_all(&ungapped)?;
        writeln!(out)?;
    }
    for seq in &msa.sequences[1..] {
        let ungapped: Vec<u8> = seq.seq.iter()
            .filter(|&&b| b != b'-' && b != b' ')
            .copied()
            .collect();
        writeln!(out, ">{}", seq_label(&seq.name, seq.seq_start, seq.seq_end, seq.orient))?;
        out.write_all(&ungapped)?;
        writeln!(out)?;
    }
    Ok(())
}

fn make_instance_row(orig_name: String, seq: Vec<u8>) -> SequenceRow {
    let (name, seq_start, seq_end, orient) = parse_seq_name_coords(&orig_name);
    let mut row = SequenceRow::new(name, seq);
    row.seq_start = seq_start;
    row.seq_end = seq_end;
    row.orient = orient;
    row
}

fn parse_seq_name_coords(name: &str) -> (String, u64, u64, Orientation) {
    if let Some(colon) = name.rfind(':') {
        let prefix = &name[..colon];
        let coords = &name[colon + 1..];
        if let Some(dash) = coords.find('-') {
            let (s, e_raw) = (&coords[..dash], &coords[dash + 1..]);
            let e = e_raw
                .trim_end_matches(|c: char| c == '+' || c == '-')
                .trim_end_matches('_');
            if let (Ok(a), Ok(b)) = (s.parse::<u64>(), e.parse::<u64>()) {
                let orient = if e_raw.ends_with("_-") {
                    Orientation::Reverse
                } else if e_raw.ends_with("_+") {
                    Orientation::Forward
                } else if a > b {
                    Orientation::Reverse
                } else {
                    Orientation::Forward
                };
                let (seq_start, seq_end) = if a <= b { (a, b) } else { (b, a) };
                return (prefix.to_string(), seq_start, seq_end, orient);
            }
        }
    }
    (name.to_string(), 0, 0, Orientation::Forward)
}

/// Build the FASTA label matching Perl's toFASTA id convention:
/// forward → `name:seq_start-seq_end`; reverse → `name:seq_end-seq_start`.
/// Falls back to bare `name` when coordinates are both zero.
fn seq_label(name: &str, seq_start: u64, seq_end: u64, orient: Orientation) -> String {
    if seq_start == 0 && seq_end == 0 {
        return name.to_string();
    }
    match orient {
        Orientation::Forward => format!("{}:{}-{}", name, seq_start, seq_end),
        Orientation::Reverse => format!("{}:{}-{}", name, seq_end, seq_start),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn roundtrip(msa: &MultiAlign) -> Vec<u8> {
        let mut buf = Vec::new();
        write(msa, &mut buf, None, None).unwrap();
        buf
    }

    #[test]
    fn write_basic() {
        let ref_seq = SequenceRow::new("ref", b"AC-GT".to_vec());
        let inst = SequenceRow::new("s1", b"AC-GT".to_vec());
        let msa = MultiAlign::from_sequences(ref_seq, vec![inst]);
        let out = roundtrip(&msa);
        let s = String::from_utf8(out).unwrap();
        // Reference is excluded from MSA output (Perl toFASTA behaviour).
        assert!(!s.contains(">ref"), "reference should not appear in MSA output");
        assert!(s.contains(">s1"));
        assert!(s.contains("AC-GT"));
    }
}
