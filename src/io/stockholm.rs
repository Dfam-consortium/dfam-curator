/// Read and write Stockholm 1.0 format (Dfam / Pfam / Rfam standard).
///
/// Key features handled:
/// - `#=GF` (per-file annotation) — ID, DE, AU, etc.
/// - `#=GC RF` — reference/consensus annotation line
/// - `#=GS` (per-sequence annotation) — currently ignored
/// - `//` — record terminator
/// - Multi-block alignments (sequences interleaved across blocks)
use std::collections::HashMap;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

use crate::alignment::{Orientation, SequenceRow, MultiAlign};

/// Read a Stockholm file and return the first record as a MultiAlign.
///
/// The `RF` annotation line (if present) is used as the reference sequence.
/// All other sequence lines become instances.
pub fn read(path: &Path) -> io::Result<MultiAlign> {
    let f = BufReader::new(std::fs::File::open(path)?);
    parse_record(f)
}

/// Read a specific record from a (possibly multi-record) Stockholm file.
///
/// `select` is interpreted as:
/// - A decimal integer → 1-based record number within the file.
/// - Any other string → exact match against the record's `#=GF ID` field.
///
/// Returns an error if no matching record is found.
pub fn read_select(path: &Path, select: &str) -> io::Result<MultiAlign> {
    use dfam_stk_io::iter_records;

    let f = BufReader::new(std::fs::File::open(path)?);

    for result in iter_records(f) {
        let record = result?;
        let matched = if let Ok(n) = select.parse::<usize>() {
            record.record_num == n
        } else {
            record.gf_first("ID").map(str::trim) == Some(select)
        };
        if matched {
            let mut buf = Vec::new();
            record.write_to(&mut buf)
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
            return parse_record(io::Cursor::new(buf));
        }
    }

    Err(io::Error::new(
        io::ErrorKind::NotFound,
        format!("no STK record matching {:?}", select),
    ))
}

/// Parse genomic coordinates from a RepeatModeler-style sequence identifier.
///
/// Recognises `prefix:start-end` (e.g. `gi|57:120437225-120436960`).
/// When `raw_start > raw_end` the sequence is reverse-strand; seq_start/end are
/// normalised so that seq_start ≤ seq_end regardless of orientation.
/// Returns the display name (prefix only), seq_start, seq_end, and orientation.
/// Falls back to the original name with start=end=1, Forward for unparseable IDs.
fn parse_seq_name_coords(name: &str) -> (String, u64, u64, Orientation) {
    if let Some(colon) = name.rfind(':') {
        let prefix = &name[..colon];
        let coords = &name[colon + 1..];
        if let Some(dash) = coords.find('-') {
            let (s, e_raw) = (&coords[..dash], &coords[dash + 1..]);
            // Optional _+ / _- strand suffix (standard Smitten/Dfam format).
            // When present, coords are always low-to-high so we read orientation
            // from the suffix.  Without the suffix we fall back to the old
            // RepeatModeler convention where start > end signals reverse strand.
            let e = e_raw.trim_end_matches(|c: char| c == '+' || c == '-')
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
                return (prefix.to_string(), a.min(b), a.max(b), orient);
            }
        }
    }
    (name.to_string(), 1, 1, Orientation::Forward)
}

/// Build a `SequenceRow` from an original Stockholm sequence name + normalised bytes.
///
/// Coordinates and strand are parsed from the name when it matches the
/// `prefix:start-end` pattern used by RepeatModeler.
fn make_instance_row(orig_name: String, seq: Vec<u8>) -> SequenceRow {
    let (display_name, seq_start, seq_end, orient) = parse_seq_name_coords(&orig_name);
    let mut row = SequenceRow::new(display_name, seq);
    row.seq_start = seq_start;
    row.seq_end   = seq_end;
    row.orient    = orient;
    row
}

fn parse_record<R: BufRead>(reader: R) -> io::Result<MultiAlign> {
    // Sequences accumulated across interleaved blocks: original name -> bytes.
    let mut seq_order: Vec<String> = Vec::new();
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
    let mut rf_line: Option<Vec<u8>> = None;
    let mut gf_id: Option<String> = None;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();

        if trimmed == "//" {
            break; // end of record
        }
        if trimmed.starts_with('#') {
            // Annotation lines.
            if let Some(rest) = trimmed.strip_prefix("#=GC RF") {
                // Reference annotation — may be split across blocks.
                // Normalize '.' (Stockholm insert-column gap) -> '-'.
                let data: Vec<u8> = rest.trim().bytes()
                    .map(|b| if b == b'.' { b'-' } else { b })
                    .collect();
                rf_line.get_or_insert_with(Vec::new).extend_from_slice(&data);
            } else if let Some(rest) = trimmed.strip_prefix("#=GF ID") {
                gf_id = Some(rest.trim().to_string());
            }
            // Other GF/GS/GR lines are silently ignored for now.
            continue;
        }
        if trimmed.is_empty() || trimmed.starts_with("# STOCKHOLM") {
            continue;
        }

        // Sequence line: "name  ACGT-..."
        // Normalize '.' -> '-' so the consensus builder and display treat them
        // uniformly as gap characters.  RepeatModeler writes STK files using '.'
        // for all gap positions; Perl MultAln converts them back to '-' on read.
        let mut parts = trimmed.splitn(2, char::is_whitespace);
        let name = match parts.next() { Some(n) => n.to_string(), None => continue };
        let data: Vec<u8> = match parts.next() {
            Some(d) => d.trim().bytes()
                .map(|b| if b == b'.' { b'-' } else { b })
                .collect(),
            None => continue,
        };
        if !seqs.contains_key(&name) {
            seq_order.push(name.clone());
        }
        seqs.entry(name).or_default().extend_from_slice(&data);
    }

    if seqs.is_empty() && rf_line.is_none() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "empty Stockholm record"));
    }

    // Build the MultiAlign.
    // If there is an RF line, use it as the reference; otherwise use the first sequence.
    let (reference, instances) = match rf_line {
        Some(rf) => {
            let ref_name = gf_id.as_deref().unwrap_or("consensus").to_string();
            let reference = SequenceRow::new(ref_name, rf);
            let instances = seq_order
                .into_iter()
                .map(|orig_name| {
                    let seq = seqs.remove(&orig_name).unwrap();
                    make_instance_row(orig_name, seq)
                })
                .collect();
            (reference, instances)
        }
        None => {
            let first_name = seq_order.remove(0);
            let first_seq = seqs.remove(&first_name).unwrap();
            let reference = SequenceRow::new(first_name, first_seq);
            let instances = seq_order
                .into_iter()
                .map(|orig_name| {
                    let seq = seqs.remove(&orig_name).unwrap();
                    make_instance_row(orig_name, seq)
                })
                .collect();
            (reference, instances)
        }
    };

    Ok(MultiAlign::from_sequences(reference, instances))
}

/// Write a MultiAlign to Stockholm 1.0 format.
///
/// Matches Perl `MultAln::toSTK(consRF => 1, idFormat => 2)`:
/// - Single uninterleaved block (no column wrapping).
/// - All gap characters (`-`) and space padding are replaced with `.`.
/// - Consensus written as `#=GC RF` before instance rows.
/// - Instance names get a `:[seq_start]-[seq_end]_[orient]` suffix.
/// - When `include_template` is true, a Dfam-style placeholder header is
///   prepended (matches Perl's `includeTemplate => 1` default).
pub fn write(
    msa: &MultiAlign,
    out: &mut dyn Write,
    family_id: Option<&str>,
    consensus_seq: Option<&[u8]>,
    include_template: bool,
) -> io::Result<()> {
    use crate::alignment::Orientation;

    let id = family_id.unwrap_or("consensus");

    writeln!(out, "# STOCKHOLM 1.0")?;

    if include_template {
        writeln!(out, "#=GF ID {}", id)?;
        writeln!(out, "#=GF DE My favorite ERVL ~:Title")?;
        writeln!(out, "#=GF AU Foobar Jones ~:Author")?;
        writeln!(out, "#=GF TP LTR/ERVL ~:Classification")?;
        writeln!(out, "#=GF OC Muridae ~:Clade1")?;
        writeln!(out, "#=GF OC Drosophila melanogaster ~:multiple OC lines allowed")?;
        writeln!(out, "#=GF TD CATATAC ~:TSD")?;
        writeln!(out, "#=GF RN [1]")?;
        writeln!(out, "#=GF RM 12343244 ~:PubMed ID")?;
        writeln!(out, "#=GF RN [2]")?;
        writeln!(out, "#=GF RM 289283 ~: Another PubMed ID")?;
        writeln!(out, "#=GF DR RepBase;{}; ~:Database Reference", id)?;
        writeln!(out, "#=GF CC ~:Public description with more details of the family")?;
        writeln!(out, "#=GF ** ~:Curation details and metdata go in the ** field")?;
        writeln!(out, "#=GF **")?;
        writeln!(out, "#=GF ** SearchStages: 3, 5, 10  ~: stages, separated by commas")?;
        writeln!(out, "#=GF ** BufferStages: 3[1-2], 5[3-6], 5 ~: 'stage'[start-end] or just 'stage'")?;
        writeln!(out, "#=GF ** HC: CACTACCCCC ~: Handbuilt consensus")?;
        writeln!(out, "#=GF BM RepeatModeler/MultAln")?;
    } else {
        writeln!(out, "#=GF ID {}", id)?;
    }

    writeln!(out, "#=GF SQ {}", msa.num_instances())?;

    // Build per-instance labels: name:seq_start-seq_end_orient (idFormat=2).
    let labels: Vec<String> = msa.sequences[1..].iter().map(|s| {
        if s.seq_start == 0 && s.seq_end == 0 {
            s.name.clone()
        } else {
            let orient = match s.orient {
                Orientation::Forward => "+",
                Orientation::Reverse => "-",
            };
            format!("{}:{}-{}_{}", s.name, s.seq_start, s.seq_end, orient)
        }
    }).collect();

    // Max name width (at least as wide as "#=GC RF ").
    let gc_rf_tag = "#=GC RF ";
    let max_name = labels.iter().map(|l| l.len()).max().unwrap_or(0)
        .max(gc_rf_tag.len());

    // #=GC RF line: consensus with gaps/spaces as '.'.
    if let Some(cons) = consensus_seq {
        let rf: Vec<u8> = cons.iter().map(|&b| if b == b'-' || b == b' ' { b'.' } else { b }).collect();
        write!(out, "{:<width$}  ", gc_rf_tag, width = max_name)?;
        out.write_all(&rf)?;
        writeln!(out)?;
    }

    // Instance rows.
    for (i, seq) in msa.sequences[1..].iter().enumerate() {
        let s: Vec<u8> = seq.seq.iter()
            .map(|&b| if b == b'-' || b == b' ' { b'.' } else { b })
            .collect();
        write!(out, "{:<width$}  ", labels[i], width = max_name)?;
        out.write_all(&s)?;
        writeln!(out)?;
    }

    writeln!(out, "//")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const SIMPLE_STK: &str = "# STOCKHOLM 1.0
#=GF ID   TestFam
seq1       ACGT
seq2       AC-T
#=GC RF    ACGT
//
";

    #[test]
    fn parse_simple() {
        let msa = parse_record(Cursor::new(SIMPLE_STK)).unwrap();
        assert_eq!(msa.num_instances(), 2);
        // GF ID "TestFam" is used as the reference name when an RF line exists.
        assert_eq!(msa.reference().unwrap().name, "TestFam");
        assert_eq!(msa.reference().unwrap().seq, b"ACGT");
        assert_eq!(msa.instance(0).unwrap().name, "seq1");
    }

    #[test]
    fn roundtrip() {
        let msa = parse_record(Cursor::new(SIMPLE_STK)).unwrap();
        let mut out = Vec::new();
        write(&msa, &mut out, Some("TestFam"), Some(b"ACGT"), false).unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(s.contains("# STOCKHOLM 1.0"));
        assert!(s.contains("#=GC RF"));
        assert!(s.contains("//"));
    }

    #[test]
    fn dot_normalized_to_dash() {
        let stk = "# STOCKHOLM 1.0\n#=GF ID Fam\nseq1  A.C\nseq2  A-C\n#=GC RF xxx\n//\n";
        let msa = parse_record(Cursor::new(stk)).unwrap();
        assert_eq!(msa.instance(0).unwrap().seq, b"A-C");
        assert_eq!(msa.instance(1).unwrap().seq, b"A-C");
    }
}
