/// Tier-1 and tier-2 lint checks for Dfam Stockholm records.
use std::collections::{HashMap, HashSet};

use crate::consensus::{build_consensus_from_sequences, ConsensusParams};
use crate::dfam::cache::Cache;
use crate::dfam::record::RawDfamRecord;

// ── Character sets ────────────────────────────────────────────────────────────

/// IUPAC/IUB nucleotide codes (upper- and lower-case).
const IUB: &[u8] = b"ACGTRYSWKMBDHVNacgtrymkswhbvdn";

/// Valid characters in sequence rows: IUB codes and `.`.
const SEQ_VALID: &[u8] = b"ACGTRYSWKMBDHVNacgtrymkswhbvdn.";

/// Valid characters in `#=GC RF`: IUB codes, `.`, and `X`/`x`.
const RF_VALID: &[u8] = b"ACGTRYSWKMBDHVNacgtrymkswhbvdnXx.";

/// All recognised `#=GF` tags (for unknown-tag detection).
const KNOWN_GF_TAGS: &[&str] = &[
    "AC", "ID", "DE", "AU", "SE", "TP", "OC", "SQ",
    "TD", "RN", "RT", "RA", "RM", "RL", "DR", "CC", "**", "KD", "BM",
];

// ── Diagnostic types ──────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Severity {
    Info,
    Warn,
    Error,
}

impl std::fmt::Display for Severity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Severity::Error => write!(f, "ERROR"),
            Severity::Warn  => write!(f, "WARN"),
            Severity::Info  => write!(f, "INFO"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Diagnostic {
    pub severity: Severity,
    pub check: &'static str,
    pub message: String,
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn err(check: &'static str, msg: impl Into<String>) -> Diagnostic {
    Diagnostic { severity: Severity::Error, check, message: msg.into() }
}

fn warn(check: &'static str, msg: impl Into<String>) -> Diagnostic {
    Diagnostic { severity: Severity::Warn, check, message: msg.into() }
}

fn info(check: &'static str, msg: impl Into<String>) -> Diagnostic {
    Diagnostic { severity: Severity::Info, check, message: msg.into() }
}

/// Return the first byte in `s` that is not in `valid`, as a char.
fn first_invalid(s: &str, valid: &[u8]) -> Option<char> {
    s.bytes().find(|b| !valid.contains(b)).map(|b| b as char)
}

/// `true` if `ac` matches `DF|DR` + 7 or 9 digits + optional `.N` version.
fn valid_ac(ac: &str) -> bool {
    let ac = ac.trim();
    let base = match ac.rfind('.') {
        Some(pos) if ac[pos + 1..].chars().all(|c| c.is_ascii_digit()) => &ac[..pos],
        _ => ac,
    };
    if base.len() < 4 { return false; }
    let (prefix, digits) = base.split_at(2);
    matches!(prefix, "DF" | "DR")
        && (digits.len() == 7 || digits.len() == 9)
        && digits.chars().all(|c| c.is_ascii_digit())
}

// ── Per-record lint checks ────────────────────────────────────────────────────

/// Run all applicable lint checks on a single record.
///
/// Tier-2 checks (TP/OC/ID against external databases) only run when
/// `cache` is `Some` and the relevant data file was loaded.
pub fn lint_record(record: &RawDfamRecord, cache: Option<&Cache>) -> Vec<Diagnostic> {
    let mut d: Vec<Diagnostic> = Vec::new();

    check_header(record, &mut d);
    check_terminator(record, &mut d);
    check_required_fields(record, &mut d);
    check_ac(record, &mut d);
    check_id(record, &mut d);
    check_de(record, &mut d);
    check_se(record, &mut d);
    check_au(record, &mut d);
    check_tp(record, &mut d);
    check_td(record, &mut d);
    check_kd(record, &mut d);
    check_sq(record, &mut d);
    check_ref_blocks(record, &mut d);
    check_rf(record, &mut d);
    check_rf_consensus(record, &mut d);
    check_sequences(record, &mut d);
    check_unknown_tags(record, &mut d);
    check_unknown_annotations(record, &mut d);

    if let Some(cache) = cache {
        tier2_tp(record, cache, &mut d);
        tier2_oc(record, cache, &mut d);
        tier2_id(record, cache, &mut d);
    }

    d
}

/// Cross-record check: report duplicate IDs within a file (case-insensitive).
///
/// Records that share an ID but carry an AC field are treated as update records
/// (INFO); records that share an ID without AC are flagged as errors.
///
/// Returns file-level diagnostics (the caller prints them with label `FILE`).
pub fn check_duplicate_ids(records: &[RawDfamRecord]) -> Vec<Diagnostic> {
    // key = lowercased ID → Vec<(original ID, record label, has AC)>
    let mut seen: HashMap<String, Vec<(String, String, bool)>> = HashMap::new();
    for r in records {
        if let Some(id) = r.gf_first("ID") {
            let id = id.trim().to_string();
            if !id.is_empty() {
                let has_ac = r.gf_has("AC");
                seen
                    .entry(id.to_lowercase())
                    .or_default()
                    .push((id, r.label(), has_ac));
            }
        }
    }
    let mut out = Vec::new();
    for (_, entries) in &seen {
        if entries.len() > 1 {
            for (orig_id, label, has_ac) in entries {
                if *has_ac {
                    out.push(info(
                        "duplicate_id_update",
                        format!(
                            "record {} has ID {:?} which also appears in other records; \
                             AC present — treating as an update record",
                            label, orig_id
                        ),
                    ));
                } else {
                    out.push(err(
                        "duplicate_id",
                        format!(
                            "record {} has ID {:?} which also appears in other records; \
                             add AC to mark as an update, or use a unique ID",
                            label, orig_id
                        ),
                    ));
                }
            }
        }
    }
    out
}

// ── Individual check functions ────────────────────────────────────────────────

fn check_terminator(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    if !r.terminated {
        d.push(err(
            "missing_terminator",
            "record is not closed by a '//' line (unexpected end of file)",
        ));
    }
}

/// Validate the `# STOCKHOLM <major>.<minor>` header line.
fn check_header(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    // Parser stores the raw header; we validate the version token here.
    // Expected: "# STOCKHOLM <digits>.<digits>"
    let ok = r.header.starts_with("# STOCKHOLM ")
        && r.header["# STOCKHOLM ".len()..]
            .split_once('.')
            .map(|(maj, min)| {
                !maj.is_empty()
                    && !min.is_empty()
                    && maj.chars().all(|c| c.is_ascii_digit())
                    && min.chars().all(|c| c.is_ascii_digit())
            })
            .unwrap_or(false);
    if !ok {
        d.push(err(
            "invalid_header",
            format!(
                "invalid Stockholm header {:?}; expected '# STOCKHOLM <N>.<N>'",
                r.header
            ),
        ));
    }
}

fn check_required_fields(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    for tag in &["DE", "AU", "TP", "OC", "SQ"] {
        if !r.gf_has(tag) {
            d.push(err("missing_required_field", format!("#{} field is absent", tag)));
        }
    }
}

fn check_ac(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    for ac in r.gf_all("AC") {
        if ac.trim().is_empty() {
            d.push(err("empty_field", "AC field is present but empty"));
        } else if !valid_ac(ac) {
            d.push(err(
                "ac_format",
                format!(
                    "AC {:?} does not match DF/DR + 7 or 9 digits (e.g. DF0000001 or DF000000001)",
                    ac.trim()
                ),
            ));
        }
    }
}

fn check_id(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    let ids = r.gf_all("ID");
    if ids.len() > 1 {
        d.push(err("id_multi_line", format!("ID appears {} times; must be a single value", ids.len())));
    }
    if let Some(id) = ids.first() {
        let id = id.trim();
        if id.is_empty() {
            d.push(err("empty_field", "ID field is present but empty"));
        } else if id.len() > 45 {
            d.push(err(
                "id_too_long",
                format!("ID {:?} is {} characters (max 45)", id, id.len()),
            ));
        } else if id.chars().all(|c| c.is_ascii_digit()) {
            d.push(err(
                "id_numeric",
                format!(
                    "ID {:?} is a purely numeric string; \
                     `stk edit --select` interprets numeric values as record numbers, \
                     making this ID unreachable by name",
                    id
                ),
            ));
        }
    }
}

fn check_de(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    let lines = r.gf_all("DE");
    if lines.len() > 1 {
        d.push(err(
            "de_multi_line",
            format!("DE appears {} times; must be a single line", lines.len()),
        ));
    }
    if let Some(de) = lines.first() {
        if de.trim().is_empty() {
            d.push(err("empty_field", "DE field is present but empty"));
        } else if de.trim().len() > 80 {
            d.push(err(
                "de_too_long",
                format!("DE is {} characters (max 80)", de.trim().len()),
            ));
        }
    }
}

fn check_se(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    if let Some(se) = r.gf_first("SE") {
        if se.trim().len() > 80 {
            d.push(err(
                "se_too_long",
                format!("SE is {} characters (max 80)", se.trim().len()),
            ));
        }
    }
}

fn check_au(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    if let Some(au) = r.gf_first("AU") {
        if au.trim().is_empty() {
            d.push(err("empty_field", "AU field is present but empty"));
        }
    }
}

fn check_tp(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    if let Some(tp) = r.gf_first("TP") {
        if tp.trim().is_empty() {
            d.push(err("empty_field", "TP field is present but empty"));
        }
    }
}

fn check_td(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    if let Some(td) = r.gf_first("TD") {
        if let Some(bad) = first_invalid(td.trim(), IUB) {
            d.push(err(
                "td_invalid_chars",
                format!("TD contains invalid character {:?} (only IUB codes allowed)", bad),
            ));
        }
    }
}

fn check_kd(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    if let Some(kd) = r.gf_first("KD") {
        if kd.trim().parse::<f64>().is_err() {
            d.push(err(
                "kd_not_numeric",
                format!("KD value {:?} is not a number", kd.trim()),
            ));
        }
    }
}

fn check_sq(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    if let Some(sq_str) = r.gf_first("SQ") {
        match sq_str.trim().parse::<usize>() {
            Ok(sq) => {
                let actual = r.sequences.len();
                if sq != actual {
                    d.push(err(
                        "sq_mismatch",
                        format!("SQ={} but {} sequence rows found", sq, actual),
                    ));
                }
            }
            Err(_) => {
                d.push(err(
                    "sq_not_numeric",
                    format!("SQ value {:?} is not a non-negative integer", sq_str.trim()),
                ));
            }
        }
    }
}

fn check_ref_blocks(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    let has_rn = r.gf_has("RN");
    let has_rm = r.gf_has("RM");
    if has_rn && !has_rm {
        d.push(warn("ref_block_incomplete", "RN is present but no RM (PubMed ID) found"));
    }
    if has_rm && !has_rn {
        d.push(warn("ref_block_incomplete", "RM is present but no RN (reference number) found"));
    }
}

fn check_rf(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    match r.gc.get("RF") {
        None => {
            d.push(err("rf_missing", "#=GC RF line is absent (required)"));
        }
        Some(rf) => {
            if let Some(bad) = first_invalid(rf, RF_VALID) {
                d.push(err(
                    "rf_invalid_chars",
                    format!("#=GC RF contains invalid character {:?} (allowed: IUB codes, X, .)", bad),
                ));
            }
            let rf_len = rf.len();
            for row in &r.sequences {
                if row.aligned_seq.len() != rf_len {
                    d.push(err(
                        "rf_length_mismatch",
                        format!(
                            "#=GC RF length {} does not match sequence {:?} length {}",
                            rf_len, row.original_id, row.aligned_seq.len()
                        ),
                    ));
                }
            }
        }
    }
}

fn check_rf_consensus(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    let rf = match r.gc.get("RF") {
        Some(rf) if !rf.is_empty() => rf,
        _ => return,
    };
    if r.sequences.is_empty() {
        return;
    }

    // Sequences in STK files use '.' for gaps; the consensus builder uses '-'.
    let converted: Vec<Vec<u8>> = r.sequences.iter()
        .map(|row| row.aligned_seq.bytes().map(|b| if b == b'.' { b'-' } else { b }).collect())
        .collect();
    let raw_seqs: Vec<&[u8]> = converted.iter().map(|v| v.as_slice()).collect();
    let called = build_consensus_from_sequences(&raw_seqs, &ConsensusParams::default());

    // Match the gap character used in the RF line ('.' or '-').
    let rf_gap: u8 = if rf.as_bytes().contains(&b'.') { b'.' } else { b'-' };
    let called_rf: Vec<u8> = called.iter()
        .map(|&b| if b == b'-' { rf_gap } else { b })
        .collect();

    // Case-insensitive comparison so lowercase RF letters are handled.
    let mismatch = called_rf.len() != rf.len()
        || called_rf.iter().zip(rf.bytes())
            .any(|(&c, r)| c.to_ascii_uppercase() != r.to_ascii_uppercase());

    if mismatch {
        d.push(warn(
            "rf_consensus_mismatch",
            "the #=GC RF line does not match the consensus called from the alignment sequences",
        ));
    }
}

fn check_sequences(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    let mut first_len: Option<usize> = None;

    for row in &r.sequences {
        // Character validation: report the first bad character per sequence.
        for (pos, &b) in row.aligned_seq.as_bytes().iter().enumerate() {
            if !SEQ_VALID.contains(&b) {
                let bad = b as char;
                let hint = if bad == '-' || bad == '~' {
                    format!(" ({:?} is not valid; use '.' for gaps)", bad)
                } else {
                    String::new()
                };
                d.push(err(
                    "seq_invalid_chars",
                    format!(
                        "sequence {:?} contains invalid character {:?} at position {}{}",
                        row.original_id,
                        bad,
                        pos + 1,
                        hint,
                    ),
                ));
                break;
            }
        }

        // All sequences must be the same aligned length.  When RF is present,
        // rf_length_mismatch (from check_rf) already flags any outlier against
        // the canonical column width, so skip the redundant first-seq comparison.
        let len = row.aligned_seq.len();
        if r.gc.get("RF").is_none() {
            match first_len {
                None => first_len = Some(len),
                Some(fl) if len != fl => {
                    d.push(err(
                        "seq_length_mismatch",
                        format!(
                            "sequence {:?} length {} differs from first sequence length {}",
                            row.original_id, len, fl
                        ),
                    ));
                }
                _ => {}
            }
        }
    }
}

fn check_unknown_tags(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    for (tag, _) in &r.gf {
        if !KNOWN_GF_TAGS.contains(&tag.as_str()) {
            d.push(info(
                "unknown_gf_tag",
                format!("unrecognised #=GF tag {:?} (possible typo?)", tag),
            ));
        }
    }
}

fn check_unknown_annotations(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    for line in &r.unknown_annotations {
        let prefix = line.split_whitespace().next().unwrap_or(line.as_str());
        // Offer a "did you mean #=GF X" hint when the bare token after # or #=
        // is a known GF tag (catches "#AU ...", "#=AU ...", etc.).
        let candidate = prefix.trim_start_matches('#').trim_start_matches('=').to_uppercase();
        let hint = if KNOWN_GF_TAGS.contains(&candidate.as_str()) {
            format!(" — did you mean '#=GF {}'?", candidate)
        } else {
            String::new()
        };
        d.push(err(
            "unknown_annotation",
            format!(
                "invalid Stockholm line {:?}; Stockholm has no comments — \
                 all annotation lines must use #=GF, #=GC, #=GS, or #=GR{}",
                prefix, hint
            ),
        ));
    }
}

// ── Tier-2 checks ────────────────────────────────────────────────────────────

fn tier2_tp(r: &RawDfamRecord, cache: &Cache, d: &mut Vec<Diagnostic>) {
    if let Some(ref cls) = cache.classification {
        if let Some(tp) = r.gf_first("TP") {
            if !tp.trim().is_empty() && !cls.contains(tp.trim()) {
                d.push(err(
                    "tp_unknown",
                    format!("TP {:?} is not a known Dfam classification string", tp.trim()),
                ));
            }
        }
    }
}

fn tier2_oc(r: &RawDfamRecord, cache: &Cache, d: &mut Vec<Diagnostic>) {
    let Some(ref tax) = cache.taxonomy else { return };

    for oc in r.gf_all("OC") {
        let oc = oc.trim();
        if oc.is_empty() || tax.contains(oc) {
            continue;
        }

        let hint = suggest_taxon(oc, cache);
        let msg = if hint.is_empty() {
            format!("OC {:?} is not a recognised NCBI scientific taxon name", oc)
        } else {
            format!(
                "OC {:?} is not a recognised NCBI scientific taxon name; {}",
                oc, hint
            )
        };
        d.push(err("oc_unknown", msg));
    }
}

/// Build a suggestion string for an unrecognised OC value.
///
/// Checks in order:
///   1. Exact case-insensitive match in the common-name/synonym table.
///   2. Fuzzy Jaro-Winkler match against scientific names (≥ 0.85 similarity).
///
/// Returns an empty string when no useful suggestion can be made.
fn suggest_taxon(query: &str, cache: &Cache) -> String {
    // ── 1. Common-name lookup ─────────────────────────────────────────────────
    if let Some(ref common) = cache.taxonomy_common {
        let q = query.to_lowercase();
        if let Some(sci) = common.get(&q) {
            return format!("did you mean {:?} (common name)?", sci);
        }
    }

    // ── 2. Fuzzy match against scientific names ───────────────────────────────
    if let Some(ref tax) = cache.taxonomy {
        let hits = fuzzy_taxon_matches(query, tax);
        if !hits.is_empty() {
            let list = hits
                .iter()
                .map(|(name, score)| format!("{:?} ({:.0}%)", name, score * 100.0))
                .collect::<Vec<_>>()
                .join(", ");
            return format!("did you mean: {}", list);
        }
    }

    String::new()
}

/// Return the top scientific-name candidates for `query` by Jaro-Winkler similarity.
///
/// Pre-filters to names that share the first character and are within ±8 characters
/// of the query length, then applies a 0.85 threshold.  Returns at most 3 results,
/// sorted by descending score.
fn fuzzy_taxon_matches(query: &str, tax: &HashSet<String>) -> Vec<(String, f64)> {
    const THRESHOLD: f64 = 0.85;
    const MAX_LEN_DIFF: usize = 8;
    const MAX_RESULTS: usize = 3;

    if query.len() < 4 {
        return vec![];
    }

    let q_lower = query.to_lowercase();
    let first   = q_lower.chars().next().unwrap_or('\0');
    let q_len   = query.len();

    let mut hits: Vec<(String, f64)> = tax
        .iter()
        .filter(|name| {
            let n = name.to_lowercase();
            n.starts_with(first) && n.len().abs_diff(q_len) <= MAX_LEN_DIFF
        })
        .filter_map(|name| {
            let score = strsim::jaro_winkler(&q_lower, &name.to_lowercase());
            if score >= THRESHOLD { Some((name.clone(), score)) } else { None }
        })
        .collect();

    hits.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    hits.truncate(MAX_RESULTS);
    hits
}

fn tier2_id(r: &RawDfamRecord, cache: &Cache, d: &mut Vec<Diagnostic>) {
    if let Some(ref names) = cache.dfam_names {
        if let Some(id) = r.gf_first("ID") {
            let id = id.trim();
            if !id.is_empty() && names.contains(&id.to_lowercase()) {
                if r.gf_has("AC") {
                    d.push(info(
                        "id_in_dfam_update",
                        format!("ID {:?} already exists in Dfam; AC present — treating as an update record", id),
                    ));
                } else {
                    d.push(err(
                        "id_in_dfam",
                        format!("ID {:?} already exists in Dfam; add AC to mark as an update, or use a unique ID", id),
                    ));
                }
            }
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dfam::cache::Cache;
    use crate::dfam::record::{RawDfamRecord, SeqRow};
    use std::collections::{HashMap, HashSet};

    fn make_record(gf: &[(&str, &str)], gc_rf: Option<&str>, seqs: &[(&str, &str)]) -> RawDfamRecord {
        let mut r = RawDfamRecord::default();
        r.record_num = 1;
        r.header = "# STOCKHOLM 1.0".to_string();
        r.terminated = true;
        for (tag, val) in gf {
            r.gf.push((tag.to_string(), val.to_string()));
        }
        if let Some(rf) = gc_rf {
            r.gc.insert("RF".to_string(), rf.to_string());
        }
        for (name, seq) in seqs {
            r.sequences.push(SeqRow::new_raw(*name, *seq));
        }
        r
    }

    fn has_check(diags: &[Diagnostic], check: &str) -> bool {
        diags.iter().any(|d| d.check == check)
    }

    fn errors(diags: &[Diagnostic]) -> Vec<&str> {
        diags.iter()
            .filter(|d| d.severity == Severity::Error)
            .map(|d| d.check)
            .collect()
    }

    #[test]
    fn clean_record_no_errors() {
        let r = make_record(
            &[("DE","A test family"),("AU","Smith J"),("TP","Interspersed_Repeat;Unknown"),
              ("OC","Mus musculus"),("SQ","1")],
            Some("xxxx"),
            &[("s1", "ACGT")],
        );
        let diags = lint_record(&r, None);
        assert!(errors(&diags).is_empty(), "unexpected errors: {:?}", errors(&diags));
    }

    #[test]
    fn missing_required_fields() {
        let r = make_record(&[], None, &[]);
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "missing_required_field"));
        assert!(has_check(&diags, "rf_missing"));
    }

    #[test]
    fn de_too_long() {
        let long_de = "A".repeat(81);
        let r = make_record(
            &[("DE", &long_de),("AU","X"),("TP","Y"),("OC","Z"),("SQ","0")],
            Some(""),
            &[],
        );
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "de_too_long"));
    }

    #[test]
    fn id_too_long() {
        let r = make_record(
            &[("ID","ThisNameIsWayTooLongForDfamAndExceedsTheFortyFiveCharacterLimit"),
              ("DE","desc"),("AU","X"),("TP","Y"),("OC","Z"),("SQ","0")],
            Some(""),
            &[],
        );
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "id_too_long"));
    }

    #[test]
    fn id_numeric_is_error() {
        let r = make_record(
            &[("ID","42"),("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","0")],
            Some(""),
            &[],
        );
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "id_numeric"), "{:?}", errors(&diags));
    }

    #[test]
    fn id_alphanumeric_not_flagged() {
        let r = make_record(
            &[("ID","L1HS"),("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","0")],
            Some(""),
            &[],
        );
        let diags = lint_record(&r, None);
        assert!(!has_check(&diags, "id_numeric"));
    }

    #[test]
    fn ac_bad_format() {
        let r = make_record(
            &[("AC","BADFORMAT"),("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","0")],
            Some(""),
            &[],
        );
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "ac_format"));
    }

    #[test]
    fn ac_good_format_7digit() {
        let r = make_record(
            &[("AC","DF0001234"),("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","0")],
            Some(""),
            &[],
        );
        let diags = lint_record(&r, None);
        assert!(!has_check(&diags, "ac_format"));
    }

    #[test]
    fn ac_good_format_9digit_versioned() {
        let r = make_record(
            &[("AC","DR000123456.2"),("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","0")],
            Some(""),
            &[],
        );
        let diags = lint_record(&r, None);
        assert!(!has_check(&diags, "ac_format"));
    }

    #[test]
    fn sq_mismatch() {
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","3")],
            Some("ACGT"),
            &[("s1","ACGT")],
        );
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "sq_mismatch"));
    }

    #[test]
    fn seq_dash_is_error() {
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","1")],
            Some("ACGT"),
            &[("s1","AC-T")],
        );
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "seq_invalid_chars"));
    }

    #[test]
    fn seq_tilde_is_error() {
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","1")],
            Some("ACGT"),
            &[("s1","AC~T")],
        );
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "seq_invalid_chars"));
    }

    #[test]
    fn rf_invalid_char() {
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","1")],
            Some("AC-T"),   // '-' not allowed in RF
            &[("s1","ACGT")],
        );
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "rf_invalid_chars"));
    }

    #[test]
    fn rf_x_is_valid() {
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","1")],
            Some("xxxx"),
            &[("s1","ACGT")],
        );
        let diags = lint_record(&r, None);
        assert!(!has_check(&diags, "rf_invalid_chars"));
    }

    #[test]
    fn valid_header_no_error() {
        let mut r = RawDfamRecord::default();
        r.header = "# STOCKHOLM 1.0".to_string();
        r.record_num = 1;
        // check_header fires before required-fields; just check no header error
        let diags = lint_record(&r, None);
        assert!(!has_check(&diags, "invalid_header"), "{:?}", diags);
    }

    #[test]
    fn bad_header_is_error() {
        use std::io::Cursor;
        use crate::dfam::record::iter_records;
        // Parser matches "# STOCKHOLM" prefix loosely; lint validates the version
        let stk = "# STOCKHOLM foo\n#=GF DE x\n#=GF AU x\n#=GF TP x\n#=GF OC x\n#=GF SQ 0\n#=GC RF \ns1 \n//\n";
        let records: Vec<_> = iter_records(Cursor::new(stk))
            .collect::<Result<_, _>>()
            .unwrap();
        let diags = lint_record(&records[0], None);
        assert!(has_check(&diags, "invalid_header"), "{:?}", diags);
    }

    #[test]
    fn unknown_annotation_prefix_is_error() {
        // #=AU is not a valid Stockholm prefix; should suggest #=GF AU
        use std::io::Cursor;
        use crate::dfam::record::iter_records;
        let stk = "\
# STOCKHOLM 1.0\n\
#=GF DE    Test family\n\
#=GF AU    Smith J\n\
#=GF TP    X\n\
#=GF OC    Mus musculus\n\
#=GF SQ    1\n\
#=GC RF    ACGT\n\
#=AU Robert Hubley, John Smith\n\
s1          ACGT\n\
//\n";
        let records: Vec<_> = iter_records(Cursor::new(stk))
            .collect::<Result<_, _>>()
            .unwrap();
        let diags = lint_record(&records[0], None);
        assert!(has_check(&diags, "unknown_annotation"),
            "expected unknown_annotation, got: {:?}",
            diags.iter().map(|d| d.check).collect::<Vec<_>>());
        let d = diags.iter().find(|d| d.check == "unknown_annotation").unwrap();
        assert_eq!(d.severity, Severity::Error);
        assert!(d.message.contains("#=AU"), "{}", d.message);
    }

    #[test]
    fn bare_hash_annotation_no_equals_is_error() {
        // "#AU ..." is missing both the "=" and the "GF" — should also be caught
        use std::io::Cursor;
        use crate::dfam::record::iter_records;
        let stk = "\
# STOCKHOLM 1.0\n\
#=GF DE    Test family\n\
#=GF AU    Smith J\n\
#=GF TP    X\n\
#=GF OC    Mus musculus\n\
#=GF SQ    1\n\
#=GC RF    ACGT\n\
#AU Robert Hubley, John Smith\n\
s1          ACGT\n\
//\n";
        let records: Vec<_> = iter_records(Cursor::new(stk))
            .collect::<Result<_, _>>()
            .unwrap();
        let diags = lint_record(&records[0], None);
        assert!(has_check(&diags, "unknown_annotation"),
            "expected unknown_annotation, got: {:?}",
            diags.iter().map(|d| d.check).collect::<Vec<_>>());
    }

    #[test]
    fn rf_consensus_match_no_warn() {
        // Four identical ACGT sequences → consensus is ACGT; RF matches.
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","4")],
            Some("ACGT"),
            &[("s1","ACGT"),("s2","ACGT"),("s3","ACGT"),("s4","ACGT")],
        );
        let diags = lint_record(&r, None);
        assert!(!has_check(&diags, "rf_consensus_mismatch"), "{:?}", diags);
    }

    #[test]
    fn rf_consensus_mismatch_warns() {
        // Sequences are all ACGT but RF claims TTTT.
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","4")],
            Some("TTTT"),
            &[("s1","ACGT"),("s2","ACGT"),("s3","ACGT"),("s4","ACGT")],
        );
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "rf_consensus_mismatch"), "{:?}", diags);
    }

    #[test]
    fn rf_consensus_dot_gap_mismatch_warns() {
        // Consensus column 2 is all-gap ('.') but RF marks it as non-gap ('x').
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","4")],
            Some("AxGT"),
            &[("s1","A.GT"),("s2","A.GT"),("s3","A.GT"),("s4","A.GT")],
        );
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "rf_consensus_mismatch"), "{:?}", diags);
    }

    #[test]
    fn missing_terminator_is_error() {
        let mut r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","0")],
            Some(""),
            &[],
        );
        r.terminated = false;
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "missing_terminator"), "{:?}", errors(&diags));
    }

    #[test]
    fn oc_common_name_suggestion() {
        let mut tax = HashSet::new();
        tax.insert("Rattus norvegicus".to_string());
        let mut common = HashMap::new();
        common.insert("rat".to_string(), "Rattus norvegicus".to_string());

        let cache = Cache {
            classification: None,
            taxonomy: Some(tax),
            taxonomy_common: Some(common),
            dfam_names: None,
        };

        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","rat"),("SQ","0")],
            Some(""),
            &[],
        );
        let diags = lint_record(&r, Some(&cache));
        let oc_diag = diags.iter().find(|d| d.check == "oc_unknown").unwrap();
        assert!(oc_diag.message.contains("Rattus norvegicus"), "{}", oc_diag.message);
        assert!(oc_diag.message.contains("common name"), "{}", oc_diag.message);
    }

    #[test]
    fn oc_fuzzy_typo_suggestion() {
        let mut tax = HashSet::new();
        tax.insert("Homo sapiens".to_string());

        let cache = Cache {
            classification: None,
            taxonomy: Some(tax),
            taxonomy_common: Some(HashMap::new()),
            dfam_names: None,
        };

        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","Homo sapien"),("SQ","0")],
            Some(""),
            &[],
        );
        let diags = lint_record(&r, Some(&cache));
        let oc_diag = diags.iter().find(|d| d.check == "oc_unknown").unwrap();
        assert!(oc_diag.message.contains("Homo sapiens"), "{}", oc_diag.message);
    }

    #[test]
    fn duplicate_ids_no_ac_is_error() {
        let mut r1 = RawDfamRecord::default();
        r1.record_num = 1;
        r1.gf.push(("ID".to_string(), "Fam1".to_string()));

        let mut r2 = RawDfamRecord::default();
        r2.record_num = 2;
        r2.gf.push(("ID".to_string(), "Fam1".to_string()));

        let mut r3 = RawDfamRecord::default();
        r3.record_num = 3;
        r3.gf.push(("ID".to_string(), "Fam2".to_string()));

        let diags = check_duplicate_ids(&[r1, r2, r3]);
        // Two records share "Fam1" with no AC → two errors
        assert_eq!(diags.len(), 2);
        assert!(diags.iter().all(|d| d.severity == Severity::Error));
        assert!(diags.iter().all(|d| d.check == "duplicate_id"));
        assert!(diags.iter().any(|d| d.message.contains("Fam1")));
    }

    #[test]
    fn duplicate_id_with_ac_is_info_update() {
        // r1: new record, no AC
        let mut r1 = RawDfamRecord::default();
        r1.record_num = 1;
        r1.gf.push(("ID".to_string(), "Fam1".to_string()));

        // r2: update record, has AC
        let mut r2 = RawDfamRecord::default();
        r2.record_num = 2;
        r2.gf.push(("ID".to_string(), "Fam1".to_string()));
        r2.gf.push(("AC".to_string(), "DF0001234".to_string()));

        let diags = check_duplicate_ids(&[r1, r2]);
        // r1 (no AC) → error; r2 (has AC) → info
        let errors: Vec<_> = diags.iter().filter(|d| d.severity == Severity::Error).collect();
        let infos: Vec<_> = diags.iter().filter(|d| d.severity == Severity::Info).collect();
        assert_eq!(errors.len(), 1, "expected one error for record without AC");
        assert_eq!(errors[0].check, "duplicate_id");
        assert_eq!(infos.len(), 1, "expected one info for update record with AC");
        assert_eq!(infos[0].check, "duplicate_id_update");
    }

    #[test]
    fn duplicate_ids_both_have_ac_are_both_info() {
        let mut r1 = RawDfamRecord::default();
        r1.record_num = 1;
        r1.gf.push(("ID".to_string(), "Fam1".to_string()));
        r1.gf.push(("AC".to_string(), "DF0001234".to_string()));

        let mut r2 = RawDfamRecord::default();
        r2.record_num = 2;
        r2.gf.push(("ID".to_string(), "Fam1".to_string()));
        r2.gf.push(("AC".to_string(), "DF0001234".to_string()));

        let diags = check_duplicate_ids(&[r1, r2]);
        assert!(diags.iter().all(|d| d.severity == Severity::Info));
        assert!(diags.iter().all(|d| d.check == "duplicate_id_update"));
    }
}
