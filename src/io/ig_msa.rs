/// Read Repbase's IG-derived aligned-FASTA (MSA) format.
///
/// The format is loosely based on IntelliGenetics/Stanford (IG) format: `;`-prefixed
/// comment lines, then a bare identifier line, then the (gapped) sequence.  Repbase
/// precedes each row with a `; FRAGMENT <start> -> <end>` comment giving the row's
/// coordinates, followed by a blank `;` comment:
///
/// ```text
/// ; FRAGMENT 1 -> 383
/// ;
/// Mariner-N5_CyaStr
/// CTGGATAATTTCGACC....----TAAGCC....ATTATCCAG
/// ; FRAGMENT 1 -> 383
/// ;
/// JAOVFP010000019.1_1
/// CTGGATAATTTCGACC....----TAAGCC....ATTATTCAG
/// ```
///
/// The first record is treated as the consensus (the MSA reference); the rest are
/// instances.  A `FRAGMENT a -> b` populates the row's `seq_start`/`seq_end` (and
/// orientation when `a > b`).  These coordinates are taken as-is and are not assumed
/// accurate — they are expected to be verified/repaired later by DisCoord.  There is
/// no `1`/`2` sequence terminator; sequence data runs until the next `;` comment or
/// end of file.
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use crate::alignment::{MultiAlign, Orientation, SequenceRow};

/// Read an IG MSA file into a `MultiAlign` (first record = reference/consensus).
pub fn read(path: &Path) -> io::Result<MultiAlign> {
    let f = BufReader::new(std::fs::File::open(path)?);

    // (identifier, gapped sequence, FRAGMENT coordinates)
    let mut entries: Vec<(String, Vec<u8>, Option<(u64, u64)>)> = Vec::new();
    let mut pending_frag: Option<(u64, u64)> = None;
    let mut cur_name: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();
    let mut cur_frag: Option<(u64, u64)> = None;

    for line in f.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        if let Some(comment) = trimmed.strip_prefix(';') {
            // A comment line delimits records: close out any row in progress, then
            // capture FRAGMENT coordinates to attach to the next identifier line.
            if let Some(name) = cur_name.take() {
                entries.push((name, std::mem::take(&mut cur_seq), cur_frag.take()));
            }
            if let Some(frag) = parse_fragment(comment.trim()) {
                pending_frag = Some(frag);
            }
            continue;
        }

        // First bare line after the comment block is the identifier; the rest is
        // sequence data (possibly wrapped across several lines).
        if cur_name.is_none() {
            cur_name = Some(trimmed.to_string());
            cur_frag = pending_frag.take();
        } else {
            cur_seq.extend_from_slice(trimmed.as_bytes());
        }
    }
    if let Some(name) = cur_name.take() {
        entries.push((name, cur_seq, cur_frag));
    }

    if entries.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "empty IG MSA file"));
    }

    let (ref_name, ref_seq, _) = entries.remove(0);
    let reference = SequenceRow::new(ref_name, ref_seq);
    let instances = entries
        .into_iter()
        .map(|(name, seq, frag)| make_instance_row(name, seq, frag))
        .collect();

    Ok(MultiAlign::from_sequences(reference, instances))
}

/// Parse the coordinates from a `FRAGMENT a -> b` comment body (the text after `;`).
/// Returns `None` for any other comment.
fn parse_fragment(comment: &str) -> Option<(u64, u64)> {
    let mut it = comment.split_whitespace();
    if it.next()? != "FRAGMENT" {
        return None;
    }
    // Tolerant of the `->` separator: just take the first two integer tokens.
    let nums: Vec<u64> = it.filter_map(|t| t.parse::<u64>().ok()).collect();
    match nums.as_slice() {
        [a, b, ..] => Some((*a, *b)),
        _ => None,
    }
}

/// Build an instance row, applying FRAGMENT coordinates and orientation.
fn make_instance_row(name: String, seq: Vec<u8>, frag: Option<(u64, u64)>) -> SequenceRow {
    let mut row = SequenceRow::new(name, seq);
    if let Some((a, b)) = frag {
        let (seq_start, seq_end, orient) = if a <= b {
            (a, b, Orientation::Forward)
        } else {
            (b, a, Orientation::Reverse)
        };
        row.seq_start = seq_start;
        row.seq_end = seq_end;
        row.orient = orient;
    }
    row
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_tmp(name: &str, body: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir();
        let path = dir.join(name);
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        path
    }

    const SAMPLE: &str = "\
; FRAGMENT 1 -> 12
;
Mariner-N5_CyaStr
CTGGA----TAAT
; FRAGMENT 1 -> 9
;
JAOVFP010000019.1_1
CTGGA----TAAT
; FRAGMENT 9 -> 1
;
JAOVFP010000020.1_1
CTGGA----TAAT
";

    #[test]
    fn reads_consensus_as_reference() {
        let path = write_tmp("ig_msa_ref.ig", SAMPLE);
        let msa = read(&path).unwrap();
        // reference (index 0) is the consensus; two instances follow.
        assert_eq!(msa.sequences[0].name, "Mariner-N5_CyaStr");
        assert_eq!(msa.sequences.len(), 3);
        assert_eq!(msa.sequences[1].name, "JAOVFP010000019.1_1");
    }

    #[test]
    fn fragment_sets_instance_coords_and_orientation() {
        let path = write_tmp("ig_msa_frag.ig", SAMPLE);
        let msa = read(&path).unwrap();
        // Forward fragment 1 -> 9.
        assert_eq!(msa.sequences[1].seq_start, 1);
        assert_eq!(msa.sequences[1].seq_end, 9);
        assert_eq!(msa.sequences[1].orient, Orientation::Forward);
        // Reverse fragment 9 -> 1 is normalized to start<=end, orient reverse.
        assert_eq!(msa.sequences[2].seq_start, 1);
        assert_eq!(msa.sequences[2].seq_end, 9);
        assert_eq!(msa.sequences[2].orient, Orientation::Reverse);
    }

    #[test]
    fn gapped_sequence_preserved() {
        let path = write_tmp("ig_msa_gap.ig", SAMPLE);
        let msa = read(&path).unwrap();
        assert_eq!(msa.sequences[0].seq, b"CTGGA----TAAT");
    }

    #[test]
    fn wrapped_sequence_lines_are_concatenated() {
        let body = "; FRAGMENT 1 -> 8\n;\ncons\nACGT\nACGT\n";
        let path = write_tmp("ig_msa_wrap.ig", body);
        let msa = read(&path).unwrap();
        assert_eq!(msa.sequences[0].seq, b"ACGTACGT");
    }
}
