//! Dfam Stockholm JSON: a structured JSON form of a Dfam seed alignment file.
//!
//! The format is specified in `docs/dfam-stockholm-json.md`, with a JSON Schema
//! in `docs/dfam-stockholm-json.schema.json`.  In short, a document is
//!
//! ```json
//! { "format": "dfam-stockholm-json", "version": "1.0", "records": [ ... ] }
//! ```
//!
//! and each record names the `#=GF` fields instead of using their two-letter
//! tags.  Where Stockholm encodes structure through line order — a reference
//! block opened by `RN`, an `FT` feature followed by its `/qualifier` lines,
//! semicolon-separated `AU` authors — the JSON makes that structure explicit.
//! Text fields stay text.
//!
//! The conversion is lossless for the content of a record but not for its
//! layout.  Going JSON → Stockholm writes the `#=GF` lines in a canonical order
//! (see [`to_record`]), one `AU` line, and an `SQ` computed from the rows.
//! A record whose structure the JSON cannot hold — two `DE` lines, an `RM`
//! before any `RN`, an `SQ` that disagrees with the rows — is refused by
//! [`from_record`] rather than silently reshaped; `stk lint` reports all of
//! these as errors, so a lint-clean file always converts.

use std::collections::BTreeMap;
use std::io::{self, Read, Write};

use dfam_stk_io::{SeqRow, StkRecord};
use serde::{Deserialize, Serialize};

/// Value of the document's `format` member.
pub const FORMAT_NAME: &str = "dfam-stockholm-json";

/// Version written by this module.  Readers accept any `1.x`.
pub const FORMAT_VERSION: &str = "1.0";

/// A whole file: the records in file order.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JsonDocument {
    pub format: String,
    pub version: String,
    pub records: Vec<JsonRecord>,
}

/// One Stockholm record (`# STOCKHOLM 1.0` … `//`).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JsonRecord {
    /// `ID`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    /// `AC`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub accession: Option<String>,
    /// `DE`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    /// `AU`, every author of every line in order.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub authors: Vec<JsonAuthor>,
    /// `SE`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
    /// `TP`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub classification: Option<String>,
    /// `OC`, one entry per line.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub clades: Vec<String>,
    /// `TD`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_site_duplication: Option<String>,
    /// `CT`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub consensus_type: Option<String>,
    /// `KD`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kimura_divergence: Option<f64>,
    /// `BM`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub build_method: Option<String>,
    /// `RN` blocks with their `RM`/`RD`/`RT`/`RA`/`RL` lines.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub references: Vec<JsonReference>,
    /// `DR`, one entry per line.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub database_references: Vec<String>,
    /// `CC`, one entry per line.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub comments: Vec<String>,
    /// `FT` features with their qualifiers.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub features: Vec<JsonFeature>,
    /// `**`, one entry per line.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub curator_notes: Vec<String>,
    /// Any other `#=GF` line, in document order.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub other_gf: Vec<JsonTagValue>,
    /// `#=GC RF`, gapped: one character per alignment column.  Named for HMMER's
    /// term, not "consensus": in Dfam the consensus is this line ungapped.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference_annotation: Option<String>,
    /// `#=GC MM`, gapped: one character per alignment column.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub model_mask: Option<String>,
    /// Any other `#=GC` line, keyed by tag.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub other_gc: BTreeMap<String, String>,
    /// Sequence rows in document order.
    #[serde(default)]
    pub sequences: Vec<JsonSequence>,
}

/// One author from an `AU` line, with its optional `ORCID:` prefix split off.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JsonAuthor {
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub orcid: Option<String>,
}

/// One literature reference block.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JsonReference {
    /// `RN [n]`
    pub number: u32,
    /// `RM`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pubmed_id: Option<String>,
    /// `RD`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub doi: Option<String>,
    /// `RT`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub title: Option<String>,
    /// `RA`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub authors: Option<String>,
    /// `RL`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub location: Option<String>,
}

/// One `FT` feature: a key + location line and the qualifier lines under it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JsonFeature {
    pub key: String,
    /// As written, e.g. `490..4569` or `join(312..602,1346..2053)`; 1-based on the
    /// ungapped consensus, descending coordinates on the minus strand.
    pub location: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub qualifiers: Vec<JsonQualifier>,
}

/// One `/name="value"` qualifier.  A valueless qualifier (`/pseudo`) has `value: null`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JsonQualifier {
    pub name: String,
    #[serde(default)]
    pub value: Option<String>,
}

/// A `#=GF` line with a tag this format has no named field for.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JsonTagValue {
    pub tag: String,
    pub value: String,
}

/// One sequence row.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JsonSequence {
    /// The row name exactly as written.  Authoritative: `coordinates` is derived from it.
    pub id: String,
    /// Aligned sequence, gaps included, as wide as the alignment.
    pub aligned: String,
    /// The Smitten identifier parsed out of `id`, when it parses.  Written for the
    /// reader's convenience; ignored on input.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coordinates: Option<JsonCoordinates>,
}

/// Genomic location parsed from a Smitten identifier (its outermost range).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JsonCoordinates {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assembly: Option<String>,
    pub sequence: String,
    /// 1-based, fully closed, `start <= end` regardless of strand.
    pub start: u64,
    pub end: u64,
    /// `"+"` or `"-"`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strand: Option<String>,
}

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

// ── Stockholm → JSON ──────────────────────────────────────────────────────────

/// Set a single-valued field, refusing a second occurrence of its tag.
fn set_once(slot: &mut Option<String>, tag: &str, value: &str) -> io::Result<()> {
    if slot.is_some() {
        return Err(invalid(format!("#=GF {} appears more than once", tag)));
    }
    *slot = Some(value.to_string());
    Ok(())
}

/// Split one `AU` line into its semicolon-separated authors.
fn parse_authors(line: &str) -> impl Iterator<Item = JsonAuthor> + '_ {
    line.split(';').map(str::trim).filter(|t| !t.is_empty()).map(|token| {
        match token.strip_prefix("ORCID:") {
            Some(rest) => {
                let (orcid, name) = rest.split_once(char::is_whitespace).unwrap_or((rest, ""));
                JsonAuthor { name: name.trim().to_string(), orcid: Some(orcid.to_string()) }
            }
            None => JsonAuthor { name: token.to_string(), orcid: None },
        }
    })
}

/// Parse an `RN` value, `[n]`.
fn parse_ref_number(value: &str) -> Option<u32> {
    value.trim().strip_prefix('[')?.strip_suffix(']')?.trim().parse().ok()
}

/// Parse an `FT` qualifier line, `/name="value"` or `/name`.
fn parse_qualifier(line: &str) -> JsonQualifier {
    let body = &line[1..];
    match body.split_once('=') {
        Some((name, value)) => {
            let value = value
                .strip_prefix('"')
                .and_then(|v| v.strip_suffix('"'))
                .unwrap_or(value);
            JsonQualifier { name: name.to_string(), value: Some(value.to_string()) }
        }
        None => JsonQualifier { name: body.to_string(), value: None },
    }
}

/// Convert a parsed Stockholm record to its JSON form.
///
/// Fails on structure the JSON cannot hold; see the module documentation.
pub fn from_record(r: &StkRecord) -> io::Result<JsonRecord> {
    from_record_inner(r).map_err(|e| invalid(format!("{}: {}", r.label(), e)))
}

fn from_record_inner(r: &StkRecord) -> io::Result<JsonRecord> {
    let mut j = JsonRecord::default();
    let mut kd: Option<String> = None;
    let mut sq: Option<String> = None;

    for (tag, value) in &r.gf {
        let v = value.as_str();
        match tag.as_str() {
            "ID" => set_once(&mut j.name, tag, v)?,
            "AC" => set_once(&mut j.accession, tag, v)?,
            "DE" => set_once(&mut j.description, tag, v)?,
            "SE" => set_once(&mut j.source, tag, v)?,
            "TP" => set_once(&mut j.classification, tag, v)?,
            "TD" => set_once(&mut j.target_site_duplication, tag, v)?,
            "CT" => set_once(&mut j.consensus_type, tag, v)?,
            "BM" => set_once(&mut j.build_method, tag, v)?,
            "KD" => set_once(&mut kd, tag, v)?,
            "SQ" => set_once(&mut sq, tag, v)?,
            "AU" => j.authors.extend(parse_authors(v)),
            "OC" => j.clades.push(v.to_string()),
            "DR" => j.database_references.push(v.to_string()),
            "CC" => j.comments.push(v.to_string()),
            "**" => j.curator_notes.push(v.to_string()),
            "RN" => {
                let number = parse_ref_number(v)
                    .ok_or_else(|| invalid(format!("RN {:?} is not of the form [n]", v)))?;
                j.references.push(JsonReference { number, ..Default::default() });
            }
            "RM" | "RD" | "RT" | "RA" | "RL" => {
                let reference = j
                    .references
                    .last_mut()
                    .ok_or_else(|| invalid(format!("{} appears before any RN", tag)))?;
                let slot = match tag.as_str() {
                    "RM" => &mut reference.pubmed_id,
                    "RD" => &mut reference.doi,
                    "RT" => &mut reference.title,
                    "RA" => &mut reference.authors,
                    _ => &mut reference.location,
                };
                set_once(slot, tag, v).map_err(|_| {
                    invalid(format!("{} appears more than once in reference [{}]", tag, reference.number))
                })?;
            }
            "FT" if v.starts_with('/') => {
                let feature = j.features.last_mut().ok_or_else(|| {
                    invalid(format!("FT qualifier {:?} has no feature above it", v))
                })?;
                feature.qualifiers.push(parse_qualifier(v));
            }
            "FT" => {
                let (key, location) = v
                    .split_once(char::is_whitespace)
                    .ok_or_else(|| invalid(format!("FT feature {:?} has no location", v)))?;
                j.features.push(JsonFeature {
                    key: key.to_string(),
                    location: location.trim().to_string(),
                    qualifiers: Vec::new(),
                });
            }
            _ => j.other_gf.push(JsonTagValue { tag: tag.clone(), value: value.clone() }),
        }
    }

    if let Some(kd) = kd {
        j.kimura_divergence = Some(
            kd.trim()
                .parse()
                .map_err(|_| invalid(format!("KD {:?} is not a number", kd)))?,
        );
    }
    // SQ is not stored: it is the length of `sequences`, and is written from it.
    if let Some(sq) = sq {
        match sq.trim().parse::<usize>() {
            Ok(n) if n == r.sequences.len() => {}
            _ => {
                return Err(invalid(format!(
                    "SQ {:?} does not match the {} sequence rows",
                    sq,
                    r.sequences.len()
                )))
            }
        }
    }

    for (tag, value) in &r.gc {
        match tag.as_str() {
            "RF" => j.reference_annotation = Some(value.clone()),
            "MM" => j.model_mask = Some(value.clone()),
            _ => {
                j.other_gc.insert(tag.clone(), value.clone());
            }
        }
    }

    j.sequences = r.sequences.iter().map(sequence_from_row).collect();
    Ok(j)
}

fn sequence_from_row(row: &SeqRow) -> JsonSequence {
    // Smitten accepts a bare name such as `consensus` as a sequence ID with no
    // range; that locates nothing, so only a row with a range gets coordinates.
    let closed = row.span.and_then(|s| s.as_1b_closed());
    let coordinates = row.sequence_id.as_ref().zip(closed).map(|(sequence, (start, end))| {
        JsonCoordinates {
            assembly: row.assembly_id.clone(),
            sequence: sequence.clone(),
            start,
            end,
            strand: row.orient.map(String::from),
        }
    });
    JsonSequence { id: row.original_id.clone(), aligned: row.aligned_seq.clone(), coordinates }
}

// ── JSON → Stockholm ──────────────────────────────────────────────────────────

/// Convert a JSON record back to a Stockholm record.
///
/// `#=GF` lines are written in this order: `ID AC DE AU SE TP OC SQ TD CT KD BM`,
/// the reference blocks (`RN RM RD RT RA RL` each), `DR CC FT **`, then
/// `other_gf` in its own order.  All authors go on one `AU` line, and `SQ` is the
/// number of sequence rows.
///
/// Fails on a value Stockholm cannot hold: a line break anywhere, or whitespace
/// in a tag or sequence ID, either of which would change the record's lines.
pub fn to_record(j: &JsonRecord, record_num: usize) -> io::Result<StkRecord> {
    check_writable(j).map_err(|e| invalid(format!("record {}: {}", record_num, e)))?;

    let mut gf: Vec<(String, String)> = Vec::new();
    let mut put = |tag: &str, value: &str| gf.push((tag.to_string(), value.to_string()));

    if let Some(v) = &j.name { put("ID", v); }
    if let Some(v) = &j.accession { put("AC", v); }
    if let Some(v) = &j.description { put("DE", v); }
    if !j.authors.is_empty() {
        let line: Vec<String> = j
            .authors
            .iter()
            .map(|a| match &a.orcid {
                Some(orcid) => format!("ORCID:{} {}", orcid, a.name).trim_end().to_string(),
                None => a.name.clone(),
            })
            .collect();
        put("AU", &line.join("; "));
    }
    if let Some(v) = &j.source { put("SE", v); }
    if let Some(v) = &j.classification { put("TP", v); }
    for v in &j.clades { put("OC", v); }
    put("SQ", &j.sequences.len().to_string());
    if let Some(v) = &j.target_site_duplication { put("TD", v); }
    if let Some(v) = &j.consensus_type { put("CT", v); }
    if let Some(kd) = j.kimura_divergence { put("KD", &kd.to_string()); }
    if let Some(v) = &j.build_method { put("BM", v); }
    for r in &j.references {
        put("RN", &format!("[{}]", r.number));
        for (tag, v) in [
            ("RM", &r.pubmed_id),
            ("RD", &r.doi),
            ("RT", &r.title),
            ("RA", &r.authors),
            ("RL", &r.location),
        ] {
            if let Some(v) = v { put(tag, v); }
        }
    }
    for v in &j.database_references { put("DR", v); }
    for v in &j.comments { put("CC", v); }
    for f in &j.features {
        put("FT", &format!("{} {}", f.key, f.location));
        for q in &f.qualifiers {
            match &q.value {
                Some(v) => put("FT", &format!("/{}=\"{}\"", q.name, v)),
                None => put("FT", &format!("/{}", q.name)),
            }
        }
    }
    for v in &j.curator_notes { put("**", v); }
    for tv in &j.other_gf { put(&tv.tag, &tv.value); }

    let mut gc = std::collections::HashMap::new();
    if let Some(rf) = &j.reference_annotation { gc.insert("RF".to_string(), rf.clone()); }
    if let Some(mm) = &j.model_mask { gc.insert("MM".to_string(), mm.clone()); }
    for (tag, v) in &j.other_gc { gc.insert(tag.clone(), v.clone()); }

    Ok(StkRecord {
        record_num,
        header: "# STOCKHOLM 1.0".to_string(),
        gf,
        gc,
        sequences: j.sequences.iter().map(|s| SeqRow::from_name_seq(&s.id, &s.aligned)).collect(),
        terminated: true,
        ..Default::default()
    })
}

/// Check that every value of `j` fits on its Stockholm line.
fn check_writable(j: &JsonRecord) -> io::Result<()> {
    fn walk(v: &serde_json::Value) -> io::Result<()> {
        match v {
            serde_json::Value::String(s) if s.contains(['\n', '\r']) => {
                Err(invalid(format!("value {:?} contains a line break", s)))
            }
            serde_json::Value::Array(a) => a.iter().try_for_each(walk),
            serde_json::Value::Object(o) => o.values().try_for_each(walk),
            _ => Ok(()),
        }
    }
    walk(&serde_json::to_value(j)?)?;

    let tags = j.other_gf.iter().map(|t| &t.tag).chain(j.other_gc.keys());
    let ids = j.sequences.iter().map(|s| &s.id);
    for (what, name) in tags.map(|t| ("tag", t)).chain(ids.map(|i| ("sequence id", i))) {
        if name.is_empty() || name.contains(char::is_whitespace) {
            return Err(invalid(format!("{} {:?} is empty or contains whitespace", what, name)));
        }
    }
    Ok(())
}

// ── Document I/O ──────────────────────────────────────────────────────────────

/// Streaming writer for a JSON document: records are written as they arrive, so
/// a large Stockholm file never has to be held in memory whole.
pub struct JsonWriter<W: Write> {
    w: W,
    count: usize,
}

impl<W: Write> JsonWriter<W> {
    /// Write the document header.
    pub fn new(mut w: W) -> io::Result<Self> {
        write!(
            w,
            "{{\n  \"format\": \"{}\",\n  \"version\": \"{}\",\n  \"records\": [",
            FORMAT_NAME, FORMAT_VERSION
        )?;
        Ok(JsonWriter { w, count: 0 })
    }

    /// Append one record.
    pub fn write_record(&mut self, record: &JsonRecord) -> io::Result<()> {
        let text = serde_json::to_string_pretty(record)?;
        if self.count > 0 {
            write!(self.w, ",")?;
        }
        // Serialized JSON strings never contain a raw newline, so indenting
        // every line only touches the layout.
        for line in text.lines() {
            write!(self.w, "\n    {}", line)?;
        }
        self.count += 1;
        Ok(())
    }

    /// Close the document and return the underlying writer.
    pub fn finish(mut self) -> io::Result<W> {
        if self.count > 0 {
            write!(self.w, "\n  ")?;
        }
        writeln!(self.w, "]\n}}")?;
        Ok(self.w)
    }
}

/// Read a whole JSON document, checking its `format` and `version`.
pub fn read_document<R: Read>(r: R) -> io::Result<JsonDocument> {
    let doc: JsonDocument = serde_json::from_reader(r).map_err(|e| invalid(e.to_string()))?;
    if doc.format != FORMAT_NAME {
        return Err(invalid(format!(
            "format is {:?}, expected {:?}",
            doc.format, FORMAT_NAME
        )));
    }
    if doc.version.split('.').next() != Some("1") {
        return Err(invalid(format!(
            "unsupported {} version {:?} (this stk reads 1.x)",
            FORMAT_NAME, doc.version
        )));
    }
    Ok(doc)
}

/// `true` if `head` (the start of a file) looks like a JSON document: its first
/// non-whitespace byte is `{`.  Stockholm starts with `#`, so the two never collide.
pub fn looks_like_json(head: &[u8]) -> bool {
    head.iter().find(|b| !b.is_ascii_whitespace()) == Some(&b'{')
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use dfam_stk_io::iter_records;
    use std::io::Cursor;

    const FULL: &str = "\
# STOCKHOLM 1.0
#=GF ID    MER4-int
#=GF AC    DF000000123.4
#=GF DE    MER4 internal sequence
#=GF AU    ORCID:0000-0001-2345-6789 Barbara McClintock; Weidong Bao
#=GF SE    RepeatModeler 2.0
#=GF TP    Interspersed_Repeat;LTR/ERV1
#=GF OC    Bos taurus
#=GF OC    Homo sapiens
#=GF SQ    2
#=GF TD    NNNNN
#=GF CT    handbuilt
#=GF KD    12.5
#=GF BM    hand
#=GF RN    [1]
#=GF RM    12490638
#=GF RN    [2]
#=GF RD    10.1093/nar/gkl1049
#=GF DR    Repbase; MER4-int
#=GF CC    Line one.
#=GF CC    Line two.
#=GF FT    CDS 1..6
#=GF FT    /product=\"Gag\"
#=GF FT    /pseudo
#=GF **    private note
#=GF XX    unknown tag
#=GC RF    ACG.TAC
GCA_000001405.15:chr1:101-106_+    ACGTTAC
consensus                          ACG.TAC
#=GC MM    ...m...
//
";

    fn parse(text: &str) -> StkRecord {
        iter_records(Cursor::new(text)).next().unwrap().unwrap()
    }

    fn stk_text(r: &StkRecord) -> String {
        let mut out = Vec::new();
        r.write_to(&mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn full_record_fields() {
        let j = from_record(&parse(FULL)).unwrap();
        assert_eq!(j.name.as_deref(), Some("MER4-int"));
        assert_eq!(j.authors.len(), 2);
        assert_eq!(j.authors[0].orcid.as_deref(), Some("0000-0001-2345-6789"));
        assert_eq!(j.authors[0].name, "Barbara McClintock");
        assert_eq!(j.authors[1], JsonAuthor { name: "Weidong Bao".into(), orcid: None });
        assert_eq!(j.clades, vec!["Bos taurus", "Homo sapiens"]);
        assert_eq!(j.kimura_divergence, Some(12.5));
        assert_eq!(j.references.len(), 2);
        assert_eq!(j.references[0].pubmed_id.as_deref(), Some("12490638"));
        assert_eq!(j.references[1].number, 2);
        assert_eq!(j.references[1].doi.as_deref(), Some("10.1093/nar/gkl1049"));
        assert_eq!(j.features[0].key, "CDS");
        assert_eq!(j.features[0].location, "1..6");
        assert_eq!(j.features[0].qualifiers[0].value.as_deref(), Some("Gag"));
        assert_eq!(j.features[0].qualifiers[1].value, None);
        assert_eq!(j.curator_notes, vec!["private note"]);
        assert_eq!(j.other_gf, vec![JsonTagValue { tag: "XX".into(), value: "unknown tag".into() }]);
        assert_eq!(j.reference_annotation.as_deref(), Some("ACG.TAC"));
        assert_eq!(j.model_mask.as_deref(), Some("...m..."));
        let c = j.sequences[0].coordinates.as_ref().unwrap();
        assert_eq!(c.assembly.as_deref(), Some("GCA_000001405.15"));
        assert_eq!((c.start, c.end, c.strand.as_deref()), (101, 106, Some("+")));
        assert!(j.sequences[1].coordinates.is_none());
    }

    #[test]
    fn canonical_record_round_trips_exactly() {
        // FULL is already in canonical order, so STK → JSON → STK is the identity.
        let r = parse(FULL);
        let back = to_record(&from_record(&r).unwrap(), 1).unwrap();
        assert_eq!(stk_text(&back), stk_text(&r));
    }

    #[test]
    fn json_round_trips_through_document() {
        let j = from_record(&parse(FULL)).unwrap();
        let mut w = JsonWriter::new(Vec::new()).unwrap();
        w.write_record(&j).unwrap();
        w.write_record(&j).unwrap();
        let bytes = w.finish().unwrap();
        assert!(looks_like_json(&bytes));
        let doc = read_document(Cursor::new(bytes)).unwrap();
        assert_eq!(doc.records, vec![j.clone(), j]);
    }

    #[test]
    fn empty_document() {
        let bytes = JsonWriter::new(Vec::new()).unwrap().finish().unwrap();
        assert!(read_document(Cursor::new(bytes)).unwrap().records.is_empty());
    }

    #[test]
    fn unrepresentable_structure_is_refused() {
        for (bad, why) in [
            ("#=GF DE    a\n#=GF DE    b\n", "more than once"),
            ("#=GF RM    123\n", "before any RN"),
            ("#=GF RN    [1]\n#=GF RM    1\n#=GF RM    2\n", "more than once in reference"),
            ("#=GF FT    /note=\"x\"\n", "no feature above"),
            ("#=GF FT    CDS\n", "no location"),
            ("#=GF KD    high\n", "not a number"),
            ("#=GF SQ    5\n", "does not match"),
        ] {
            let text = format!("# STOCKHOLM 1.0\n{}seq1  ACGT\n//\n", bad);
            let err = from_record(&parse(&text)).unwrap_err().to_string();
            assert!(err.contains(why), "{:?}: {}", bad, err);
        }
    }

    #[test]
    fn values_that_would_break_stockholm_lines_are_refused() {
        let mut j = from_record(&parse(FULL)).unwrap();
        j.description = Some("ok\n#=GF AC DF000000001".into());
        assert!(to_record(&j, 1).unwrap_err().to_string().contains("line break"));

        let mut j = from_record(&parse(FULL)).unwrap();
        j.sequences[0].id = "two words".into();
        assert!(to_record(&j, 1).unwrap_err().to_string().contains("whitespace"));
    }

    #[test]
    fn rejects_wrong_format_and_major_version() {
        let wrong = r#"{"format":"other","version":"1.0","records":[]}"#;
        assert!(read_document(Cursor::new(wrong)).is_err());
        let future = r#"{"format":"dfam-stockholm-json","version":"2.0","records":[]}"#;
        assert!(read_document(Cursor::new(future)).is_err());
    }

    #[test]
    fn rejects_unknown_record_field() {
        let typo = r#"{"format":"dfam-stockholm-json","version":"1.0",
                       "records":[{"classifcation":"x","sequences":[]}]}"#;
        assert!(read_document(Cursor::new(typo)).is_err());
    }
}
