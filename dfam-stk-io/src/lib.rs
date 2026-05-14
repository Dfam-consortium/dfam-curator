/// Stockholm 1.0 parser with Smitten identifier enrichment.
///
/// Provides `StkRecord` (one `//`-terminated block), `SeqRow` (one sequence
/// row with parsed coordinates), and `iter_records` (streaming iterator).
pub use smitten::IDVersion;

use smitten::Identifier;
use std::collections::HashMap;
use std::io::{self, BufRead, Write};

// ── SeqRow ────────────────────────────────────────────────────────────────────

/// One sequence row from a Stockholm alignment block.
#[derive(Debug, Clone)]
pub struct SeqRow {
    /// Exact identifier string as written in the file.
    pub original_id: String,
    /// Assembly accession parsed by Smitten, if recognised.
    pub assembly_id: Option<String>,
    /// Chromosome / scaffold name.  `None` when the identifier cannot be
    /// parsed (e.g. bare consensus labels).
    pub sequence_id: Option<String>,
    /// 0-based start coordinate from the Smitten identifier, if present.
    pub seq_start: Option<u64>,
    /// 0-based end coordinate from the Smitten identifier, if present.
    pub seq_end: Option<u64>,
    /// Strand orientation (`'+'` or `'-'`), if present.
    pub orient: Option<char>,
    /// Smitten identifier version inferred during parsing.
    pub inferred_version: Option<IDVersion>,
    /// Aligned sequence as stored in the file (gap characters preserved).
    pub aligned_seq: String,
}

impl SeqRow {
    /// Construct a `SeqRow` from a raw name/sequence pair.
    ///
    /// Tries Smitten parsing first; falls back gracefully so that sequences
    /// with unparseable identifiers (consensus lines, etc.) are still stored.
    pub fn from_name_seq(name: &str, aligned: &str) -> Self {
        let parsed = Identifier::from_unknown_format(name, false, true)
            .ok()
            .and_then(|(id, version)| id.normalize().ok().map(|n| (n, version)));

        match parsed {
            Some((norm, version)) => SeqRow {
                original_id: name.to_string(),
                assembly_id: norm.assembly_id,
                sequence_id: Some(norm.sequence_id),
                seq_start: norm.ranges.first().map(|r| r.start as u64),
                seq_end: norm.ranges.first().map(|r| r.end as u64),
                orient: norm.ranges.first().map(|r| r.orientation),
                inferred_version: Some(version),
                aligned_seq: aligned.to_string(),
            },
            None => SeqRow {
                original_id: name.to_string(),
                assembly_id: None,
                sequence_id: None,
                seq_start: None,
                seq_end: None,
                orient: None,
                inferred_version: None,
                aligned_seq: aligned.to_string(),
            },
        }
    }

    /// Construct a `SeqRow` without attempting Smitten identifier parsing.
    ///
    /// Useful in tests or when only the raw name and sequence are needed.
    pub fn new_raw(name: impl Into<String>, aligned: impl Into<String>) -> Self {
        SeqRow {
            original_id: name.into(),
            assembly_id: None,
            sequence_id: None,
            seq_start: None,
            seq_end: None,
            orient: None,
            inferred_version: None,
            aligned_seq: aligned.into(),
        }
    }
}

// ── StkRecord ─────────────────────────────────────────────────────────────────

/// A single record parsed from a Dfam Stockholm file.
#[derive(Debug, Default)]
pub struct StkRecord {
    /// 1-based record number within the file.
    pub record_num: usize,
    /// Line number (1-based) of the `# STOCKHOLM 1.0` header.
    pub start_line: usize,
    /// The raw header line (e.g. `# STOCKHOLM 1.0`) as read from the file.
    pub header: String,
    /// All `#=GF` annotations in document order: `(tag, value)`.
    pub gf: Vec<(String, String)>,
    /// All `#=GC` annotations, concatenated per tag across interleaved blocks.
    pub gc: HashMap<String, String>,
    /// Sequence rows in document order; multi-block rows are concatenated.
    pub sequences: Vec<SeqRow>,
    /// Lines beginning with `#` that are not a valid Stockholm annotation.
    pub unknown_annotations: Vec<String>,
    /// `true` when the record was closed by a `//` line; `false` when EOF
    /// was reached without a terminator.
    pub terminated: bool,
}

impl StkRecord {
    /// All values for a given GF tag (e.g. multiple `OC` clade lines).
    pub fn gf_all(&self, tag: &str) -> Vec<&str> {
        self.gf
            .iter()
            .filter(|(t, _)| t == tag)
            .map(|(_, v)| v.as_str())
            .collect()
    }

    /// First value for a GF tag, or `None`.
    pub fn gf_first(&self, tag: &str) -> Option<&str> {
        self.gf
            .iter()
            .find(|(t, _)| t == tag)
            .map(|(_, v)| v.as_str())
    }

    /// Whether a GF tag appears at least once.
    pub fn gf_has(&self, tag: &str) -> bool {
        self.gf.iter().any(|(t, _)| t == tag)
    }

    /// Display label used in diagnostics: `record:N/ID` when ID is present,
    /// otherwise `record:N`.
    pub fn label(&self) -> String {
        match self.gf_first("ID") {
            Some(id) if !id.trim().is_empty() => {
                format!("record:{}/{}", self.record_num, id.trim())
            }
            _ => format!("record:{}", self.record_num),
        }
    }

    /// Serialise the record back to Stockholm format.
    ///
    /// Always writes a canonical `# STOCKHOLM 1.0` header.  The `#=GC RF`
    /// line (when present) is emitted immediately before the sequence rows as
    /// the last header annotation, regardless of where it appeared in the
    /// input.  All other `#=GC` lines follow the sequence rows, sorted
    /// alphabetically for deterministic output.
    pub fn write_to<W: Write>(&self, w: &mut W) -> io::Result<()> {
        writeln!(w, "# STOCKHOLM 1.0")?;

        for (tag, value) in &self.gf {
            writeln!(w, "#=GF {} {}", tag, value)?;
        }

        if !self.sequences.is_empty() || !self.gc.is_empty() {
            let col = self
                .sequences
                .iter()
                .map(|row| row.original_id.len())
                .chain(self.gc.keys().map(|k| "#=GC ".len() + k.len()))
                .max()
                .unwrap_or(0)
                + 2;

            // RF goes before sequences as the final header line.
            if let Some(rf) = self.gc.get("RF") {
                writeln!(w, "{:<width$}{}", "#=GC RF", rf, width = col)?;
            }

            for row in &self.sequences {
                writeln!(w, "{:<width$}{}", row.original_id, row.aligned_seq, width = col)?;
            }

            // All other #=GC tags follow the sequences, sorted alphabetically.
            let mut gc_tags: Vec<&String> = self.gc.keys().filter(|k| k.as_str() != "RF").collect();
            gc_tags.sort();
            for tag in gc_tags {
                writeln!(
                    w,
                    "{:<width$}{}",
                    format!("#=GC {}", tag),
                    self.gc[tag],
                    width = col
                )?;
            }
        }

        writeln!(w, "//")?;
        Ok(())
    }
}

// ── Iterator ──────────────────────────────────────────────────────────────────

/// Streaming iterator that yields one `StkRecord` per `//`-terminated block.
pub struct StkRecordIter<R: BufRead> {
    reader: R,
    line_num: usize,
    record_num: usize,
    done: bool,
    /// When `false`, sequence rows are stored via `SeqRow::new_raw` — skipping
    /// Smitten identifier parsing.  Use this for operations that only touch
    /// `#=GF` fields and never inspect parsed coordinates.
    parse_sequences: bool,
}

impl<R: BufRead> StkRecordIter<R> {
    pub fn new(reader: R) -> Self {
        StkRecordIter { reader, line_num: 0, record_num: 0, done: false, parse_sequences: true }
    }

    pub fn new_raw(reader: R) -> Self {
        StkRecordIter { reader, line_num: 0, record_num: 0, done: false, parse_sequences: false }
    }
}

impl<R: BufRead> Iterator for StkRecordIter<R> {
    type Item = io::Result<StkRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let mut in_record = false;
        let mut record = StkRecord::default();

        loop {
            let mut line = String::new();
            match self.reader.read_line(&mut line) {
                Err(e) => return Some(Err(e)),
                Ok(0) => {
                    self.done = true;
                    return if in_record { Some(Ok(record)) } else { None };
                }
                Ok(_) => {}
            }
            self.line_num += 1;
            let trimmed = line.trim_end();

            if trimmed.starts_with("# STOCKHOLM") {
                if !in_record {
                    self.record_num += 1;
                    record.record_num = self.record_num;
                    record.start_line = self.line_num;
                    record.header = trimmed.to_string();
                    in_record = true;
                }
                continue;
            }

            if !in_record {
                continue;
            }

            if trimmed == "//" {
                record.terminated = true;
                return Some(Ok(record));
            }

            if trimmed.is_empty() {
                continue;
            }

            if trimmed.starts_with("#=GF") {
                let rest = trimmed[4..].trim_start();
                let mut it = rest.splitn(2, char::is_whitespace);
                let tag = it.next().unwrap_or("").to_string();
                let value = it.next().map(|s| s.trim().to_string()).unwrap_or_default();
                if !tag.is_empty() {
                    record.gf.push((tag, value));
                }
                continue;
            }

            if trimmed.starts_with("#=GC") {
                let rest = trimmed[4..].trim_start();
                let mut it = rest.splitn(2, char::is_whitespace);
                let tag = it.next().unwrap_or("").to_string();
                let value = it.next().map(|s| s.trim().to_string()).unwrap_or_default();
                if !tag.is_empty() {
                    record.gc.entry(tag).or_default().push_str(&value);
                }
                continue;
            }

            if trimmed.starts_with('#') {
                if !trimmed.starts_with("#=GS") && !trimmed.starts_with("#=GR") {
                    record.unknown_annotations.push(trimmed.to_string());
                }
                continue;
            }

            // Sequence line: "name   ACGT..."
            // Each row is always a complete, independent sequence — Stockholm
            // has no interleaved/wrapped-block format.  Duplicate identifiers
            // (tandem repeats, split alignments) are valid and must be kept as
            // separate rows, never concatenated.
            let mut it = trimmed.splitn(2, char::is_whitespace);
            if let (Some(name), Some(seq)) = (it.next(), it.next()) {
                let seq = seq.trim();
                let row = if self.parse_sequences {
                    SeqRow::from_name_seq(name, seq)
                } else {
                    SeqRow::new_raw(name, seq)
                };
                record.sequences.push(row);
            }
        }
    }
}

/// Create a streaming record iterator from any `BufRead`.
pub fn iter_records<R: BufRead>(reader: R) -> StkRecordIter<R> {
    StkRecordIter::new(reader)
}

/// Like `iter_records` but skips Smitten identifier parsing on sequence rows.
///
/// Sequence rows are stored with their original identifier and aligned sequence
/// intact, but the parsed coordinate fields are left empty.  Use this for
/// operations that only read or write `#=GF` annotations and never inspect
/// sequence coordinates.
pub fn iter_records_raw<R: BufRead>(reader: R) -> StkRecordIter<R> {
    StkRecordIter::new_raw(reader)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const MULTI: &str = "\
# STOCKHOLM 1.0
#=GF ID    Fam1
#=GF DE    First family
#=GF AU    Smith J
#=GF TP    Interspersed_Repeat;Unknown
#=GF OC    Mus musculus
#=GF SQ    2
#=GC RF    xxxx
seq1        ACGT
seq2        ACGT
//
# STOCKHOLM 1.0
#=GF ID    Fam2
#=GF DE    Second family
#=GF AU    Jones B
#=GF TP    Interspersed_Repeat;Unknown
#=GF OC    Homo sapiens
#=GF SQ    1
#=GC RF    xx
seq3        AC
//
";

    #[test]
    fn parses_two_records() {
        let records: Vec<_> = iter_records(Cursor::new(MULTI))
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].gf_first("ID"), Some("Fam1"));
        assert_eq!(records[1].gf_first("ID"), Some("Fam2"));
        assert_eq!(records[0].sequences.len(), 2);
        assert_eq!(records[1].sequences.len(), 1);
    }

    #[test]
    fn multi_oc_preserved() {
        const STK: &str = "\
# STOCKHOLM 1.0
#=GF ID    Fam1
#=GF DE    Test
#=GF AU    A B
#=GF TP    X
#=GF OC    Mammalia
#=GF OC    Aves
#=GF SQ    1
#=GC RF    xx
s1          AC
//
";
        let records: Vec<_> = iter_records(Cursor::new(STK))
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(records[0].gf_all("OC"), vec!["Mammalia", "Aves"]);
    }

    #[test]
    fn label_with_id() {
        let mut r = StkRecord::default();
        r.record_num = 3;
        r.gf.push(("ID".to_string(), "MyFam".to_string()));
        assert_eq!(r.label(), "record:3/MyFam");
    }

    #[test]
    fn label_without_id() {
        let mut r = StkRecord::default();
        r.record_num = 7;
        assert_eq!(r.label(), "record:7");
    }

    #[test]
    fn seq_row_id_and_seq_preserved() {
        const STK: &str = "\
# STOCKHOLM 1.0
#=GF SQ    1
#=GC RF    ACGT
consensus   ACGT
//
";
        let records: Vec<_> = iter_records(Cursor::new(STK))
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(records[0].sequences.len(), 1);
        let row = &records[0].sequences[0];
        assert_eq!(row.original_id, "consensus");
        assert_eq!(row.aligned_seq, "ACGT");
    }

    #[test]
    fn duplicate_id_rows_kept_separate() {
        // Two rows with the same identifier (e.g. a tandem repeat) must be
        // stored as separate SeqRows, never concatenated.
        const STK: &str = "\
# STOCKHOLM 1.0
#=GF SQ    3
#=GC RF    xxxx
seq1        ACGT
seq2        TGCA
seq1        AAAA
//
";
        let records: Vec<_> = iter_records(Cursor::new(STK))
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(records[0].sequences.len(), 3);
        assert_eq!(records[0].sequences[0].aligned_seq, "ACGT");
        assert_eq!(records[0].sequences[1].aligned_seq, "TGCA");
        assert_eq!(records[0].sequences[2].aligned_seq, "AAAA");
    }
}
