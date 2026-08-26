/// Read and write Clustal ALN format.
///
/// Clustal ALN is an interleaved block format: sequences are written in 60-column
/// blocks, each block followed by a conservation line.  The conservation line uses:
///   `*`  — identical residue in every sequence at that column
///   ` `  — not conserved (or gap-only column)
///
/// For nucleotide alignments ClustalW only uses `*` and ` `; the `:` and `.`
/// symbols are protein-specific (Gonnet PAM 250 conserved groups).
use std::collections::HashMap;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

use aln_core::msa::{MultiAlign, SequenceRow};
use aln_core::Strand;
use dfam_stk_io::msa::parse_seq_name_coords;

const BLOCK_WIDTH: usize = 60;

/// Read a Clustal ALN file into a `MultiAlign`.
///
/// The first sequence becomes the reference; all others are instances.
/// Trailing per-line position counts (integers appended after the alignment
/// characters) are silently discarded.
pub fn read(path: &Path) -> io::Result<MultiAlign> {
    let f = BufReader::new(std::fs::File::open(path)?);
    let mut found_header = false;
    let mut seq_order: Vec<String> = Vec::new();
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();

    for line in f.lines() {
        let line = line?;

        if !found_header {
            if line.trim().to_ascii_uppercase().starts_with("CLUSTAL") {
                found_header = true;
            }
            continue;
        }

        if line.trim().is_empty() {
            continue;
        }

        // Conservation lines start with whitespace — skip.
        if line.starts_with(' ') || line.starts_with('\t') {
            continue;
        }

        // Sequence line: first whitespace-delimited token is the name.
        let mut iter = line.splitn(2, char::is_whitespace);
        let name = match iter.next() {
            Some(n) if !n.is_empty() => n.to_string(),
            _ => continue,
        };
        let rest = iter.next().unwrap_or("").trim();

        // Retain only valid alignment characters.  Filtering out digits also
        // discards any trailing position-count integer some writers append.
        let data: Vec<u8> = rest
            .bytes()
            .filter(|&b| b.is_ascii_alphabetic() || b == b'-')
            .map(|b| b.to_ascii_uppercase())
            .collect();

        if data.is_empty() {
            continue;
        }

        if !seqs.contains_key(&name) {
            seq_order.push(name.clone());
        }
        seqs.entry(name).or_default().extend_from_slice(&data);
    }

    if !found_header {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "no CLUSTAL header found",
        ));
    }
    if seq_order.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "empty Clustal alignment",
        ));
    }

    let first_name = seq_order.remove(0);
    let first_seq = seqs.remove(&first_name).unwrap();
    // First sequence is the reference; use its name as-is (no coord parsing).
    let reference = SequenceRow::new(first_name, first_seq);
    // Instances: parse genomic coordinates from the name so that round-trips
    // through STK format reconstruct the correct per-instance labels.
    let instances = seq_order
        .into_iter()
        .map(|name| {
            let seq = seqs.remove(&name).unwrap();
            make_instance_row(name, seq)
        })
        .collect();

    MultiAlign::from_sequences(reference, instances)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))
}

/// Write a `MultiAlign` as a Clustal ALN file.
///
/// All sequences (including the reference at index 0) are written so the
/// consensus/RF row appears alongside the instances.  The name column is wide
/// enough to accommodate the longest label plus six spaces of padding.
pub fn write(msa: &MultiAlign, out: &mut dyn Write) -> io::Result<()> {
    if msa.sequences.is_empty() {
        return Ok(());
    }

    let labels: Vec<String> = msa.sequences.iter().map(seq_label).collect();
    let seqs: Vec<&[u8]> = msa.sequences.iter().map(|s| s.seq.as_slice()).collect();
    let aln_len = seqs[0].len();
    let name_col = labels.iter().map(|l| l.len()).max().unwrap_or(0).max(10) + 6;

    writeln!(out, "CLUSTAL W (1.81) multiple sequence alignment")?;
    writeln!(out)?;

    let mut pos = 0;
    while pos < aln_len {
        let end = (pos + BLOCK_WIDTH).min(aln_len);

        for (label, seq) in labels.iter().zip(seqs.iter()) {
            let block = std::str::from_utf8(&seq[pos..end]).unwrap_or("?");
            writeln!(out, "{:<width$}{}", label, block, width = name_col)?;
        }

        // Conservation: '*' when every non-gap base in the column is identical.
        let conservation: String = (pos..end)
            .map(|col| {
                let bases: Vec<u8> = seqs
                    .iter()
                    .map(|s| s[col])
                    .filter(|&b| b != b'-' && b != b'.' && b != b' ')
                    .map(|b| b.to_ascii_uppercase())
                    .collect();
                if !bases.is_empty() && bases.iter().all(|&b| b == bases[0]) {
                    '*'
                } else {
                    ' '
                }
            })
            .collect();

        writeln!(out, "{:<width$}{}", "", conservation, width = name_col)?;
        writeln!(out)?;

        pos = end;
    }

    Ok(())
}

/// Build a `SequenceRow` from a raw name+sequence, parsing genomic coordinates
/// from the name when it matches RepeatModeler's `prefix:start-end[_strand]`
/// format. Uses the same parser as the Stockholm reader so round-trips agree.
fn make_instance_row(orig_name: String, seq: Vec<u8>) -> SequenceRow {
    let (name, span, orient) = parse_seq_name_coords(&orig_name);
    let mut row = SequenceRow::new(name, seq);
    row.span = span;
    row.orient = orient;
    row
}

/// `name:start-end`, 1-based closed and written descending on the minus
/// strand as Perl's `toFASTA` does; the bare name when there are no
/// coordinates.
fn seq_label(s: &SequenceRow) -> String {
    match s.span.and_then(|sp| sp.as_1b_closed()) {
        None => s.name.clone(),
        Some((start, end)) => match s.orient {
            Strand::Plus => format!("{}:{}-{}", s.name, start, end),
            Strand::Minus => format!("{}:{}-{}", s.name, end, start),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_msa() -> MultiAlign {
        let ref_seq = SequenceRow::new("consensus", b"ACGT-ACGT".to_vec());
        let inst1 = SequenceRow::new("seq1", b"ACGT-ACGT".to_vec());
        let inst2 = SequenceRow::new("seq2", b"ACGT-TTTT".to_vec());
        MultiAlign::from_sequences(ref_seq, vec![inst1, inst2]).expect("fixture rows are equal width")
    }

    #[test]
    fn write_produces_clustal_header() {
        let msa = make_msa();
        let mut buf = Vec::new();
        write(&msa, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("CLUSTAL W"));
        assert!(s.contains("consensus"));
        assert!(s.contains("seq1"));
        assert!(s.contains("seq2"));
    }

    #[test]
    fn conservation_line_correct() {
        let msa = make_msa();
        let mut buf = Vec::new();
        write(&msa, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        // The conservation line is the first line that starts with whitespace
        // and contains '*' or ' ' (not blank).
        let cons_line = s
            .lines()
            .find(|l| l.starts_with(' ') && !l.trim().is_empty())
            .expect("no conservation line found");
        let aln_part = cons_line.trim_start();
        // ref/seq1=ACGT-ACGT, seq2=ACGT-TTTT
        // Col 0-3: ACGT identical → "****"
        assert_eq!(&aln_part[..4], "****");
        // Col 4: all gaps → " "
        assert_eq!(&aln_part[4..5], " ");
        // Cols 5-8: A,A,T / C,C,T / G,G,T / T,T,T → "   *"
        assert_eq!(&aln_part[5..8], "   ");
        assert_eq!(&aln_part[8..9], "*");
    }
}
