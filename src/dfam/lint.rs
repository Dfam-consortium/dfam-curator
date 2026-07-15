/// Tier-1, tier-2, and network lint checks for Dfam Stockholm records.
use std::collections::{HashMap, HashSet};

use crate::consensus::{build_consensus_from_sequences, ConsensusParams};
use crate::dfam::cache::Cache;
use crate::dfam::record::RawDfamRecord;
use dfam_stk_io::{is_gap, IDVersion, DFAM_GAP};

// ── Character sets ────────────────────────────────────────────────────────────

/// IUPAC/IUB nucleotide codes (upper- and lower-case).
const IUB: &[u8] = b"ACGTRYSWKMBDHVNacgtrymkswhbvdn";

/// Valid non-gap characters in sequence rows: IUB codes.  Gaps are handled separately,
/// since Stockholm allows four gap characters but Dfam standardizes on `.` (see `is_gap`).
const SEQ_VALID: &[u8] = b"ACGTRYSWKMBDHVNacgtrymkswhbvdn";

/// Valid characters in `#=GC RF`: IUB codes, `.`, and `X`/`x`.
const RF_VALID: &[u8] = b"ACGTRYSWKMBDHVNacgtrymkswhbvdnXx.";

/// Valid characters in `#=GC MM` (HMMER model mask): `m` for a masked column, `.` otherwise.
const MM_VALID: &[u8] = b"m.";

/// All recognised `#=GC` annotations.
const KNOWN_GC_TAGS: &[&str] = &["RF", "MM"];

/// All recognised `#=GF` tags (for unknown-tag detection).
const KNOWN_GF_TAGS: &[&str] = &[
    "AC", "ID", "DE", "AU", "SE", "TP", "OC", "SQ",
    "TD", "CT", "RN", "RT", "RA", "RM", "RL", "RD", "DR", "CC", "**", "KD", "BM",
];

// ── Consensus type (`#=GF CT`) ────────────────────────────────────────────────

/// A reserved word for the optional `#=GF CT` (consensus type) field.
struct ConsensusType {
    /// The reserved word as it appears in the file (matched case-insensitively).
    word: &'static str,
    /// Whether `#=GC RF` is expected to be reproducible by calling the consensus
    /// from the alignment.  `false` for consensus types that a curator authored
    /// by hand, where a mismatch with the called consensus is the normal state.
    rf_is_called: bool,
}

/// The reserved words accepted in `#=GF CT`.  Extend this list to add new types.
const CONSENSUS_TYPES: &[ConsensusType] = &[
    ConsensusType { word: "handbuilt", rf_is_called: false },
];

fn consensus_type(word: &str) -> Option<&'static ConsensusType> {
    CONSENSUS_TYPES.iter().find(|ct| ct.word.eq_ignore_ascii_case(word.trim()))
}

/// `true` if the record's `#=GF CT` marks its `#=GC RF` line as one that was not
/// called from the alignment (e.g. `handbuilt`), so consensus checks should be skipped.
pub fn rf_is_handcurated(r: &RawDfamRecord) -> bool {
    r.gf_first("CT")
        .and_then(consensus_type)
        .is_some_and(|ct| !ct.rf_is_called)
}

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
    check_ct(record, &mut d);
    check_kd(record, &mut d);
    check_ref_blocks(record, &mut d);
    check_block_format(record, &mut d);
    check_mm(record, &mut d);
    check_unknown_tags(record, &mut d);
    check_unknown_annotations(record, &mut d);

    // A block-format record parses into one row per sequence *per block*, so every
    // alignment-derived check would report a downstream artefact of the misparse
    // (SQ=2 vs "4 rows", a nonsense consensus, ragged lengths) rather than a real
    // problem.  Report the block format alone and let the curator unwrap the file first.
    if record.block_separator_line.is_none() {
        check_sq(record, &mut d);
        check_rf(record, &mut d);
        check_rf_consensus(record, &mut d);
        check_sequences(record, &mut d);
    }

    if let Some(cache) = cache {
        tier2_tp(record, cache, &mut d);
        tier2_oc(record, cache, &mut d);
        tier2_id(record, cache, &mut d);
    }

    d
}

/// Cross-record check: summarise how many records carry an `AC` field.
///
/// An `AC` is not something a submitter assigns — it marks the record as an update
/// to a family that already exists in Dfam.  Curators sometimes add one because they
/// assume it is a required field, which would silently turn a new submission into a
/// replacement of an unrelated family.  Report the count so that is visible.
///
/// Returns file-level diagnostics (the caller prints them with label `FILE`).
pub fn check_ac_summary(records: &[RawDfamRecord]) -> Vec<Diagnostic> {
    /// At most this many record→accession pairs are named before eliding the rest.
    const MAX_LISTED: usize = 5;

    let with_ac: Vec<String> = records
        .iter()
        .filter_map(|r| {
            let ac = r.gf_all("AC").into_iter().find(|ac| !ac.trim().is_empty())?;
            Some(format!("{} = {}", r.label(), ac.trim()))
        })
        .collect();

    if with_ac.is_empty() {
        return Vec::new();
    }

    let mut listed = with_ac.iter().take(MAX_LISTED).cloned().collect::<Vec<_>>().join(", ");
    if with_ac.len() > MAX_LISTED {
        listed.push_str(&format!(", and {} more", with_ac.len() - MAX_LISTED));
    }

    let clause = if with_ac.len() == 1 {
        "1 record carries an AC, so it will replace the released Dfam family it names \
         rather than create a new one"
            .to_string()
    } else {
        format!(
            "{} records carry an AC, so they will replace the released Dfam families they \
             name rather than create new ones",
            with_ac.len()
        )
    };

    vec![info(
        "ac_update_records",
        format!(
            "{}: {}.  Dfam assigns AC on release — remove it from new submissions.",
            clause, listed,
        ),
    )]
}

/// Cross-record check: encourage ORCIDs by counting records whose `AU` field names
/// at least one author without one.
///
/// An ORCID disambiguates a curator from everyone who shares their name, so Dfam can
/// credit them reliably.  Supplying one is optional, hence INFO: this reports a count
/// rather than flagging individual records.
///
/// Returns file-level diagnostics (the caller prints them with label `FILE`).
pub fn check_orcid_summary(records: &[RawDfamRecord]) -> Vec<Diagnostic> {
    // A record is "incomplete" if any author token lacks an ORCID: prefix.  Records with
    // no AU at all are skipped — missing_required_field already covers those.
    let mut incomplete = 0usize;
    let mut with_au = 0usize;

    for r in records {
        let authors: Vec<&str> = r
            .gf_all("AU")
            .iter()
            .flat_map(|au| au.split(';'))
            .map(str::trim)
            .filter(|t| !t.is_empty())
            .collect();

        if authors.is_empty() {
            continue;
        }
        with_au += 1;
        if authors.iter().any(|t| !t.starts_with("ORCID:")) {
            incomplete += 1;
        }
    }

    if incomplete == 0 {
        return Vec::new();
    }

    vec![info(
        "orcid_missing",
        format!(
            "{} of {} records with an AU field name at least one author without an ORCID.  \
             ORCIDs are optional, but they let Dfam credit curators unambiguously — \
             prefix the name with the identifier, e.g. \
             '#=GF AU    ORCID:0000-0001-2345-6789 Barbara McClintock'.",
            incomplete, with_au,
        ),
    )]
}

/// Cross-record check: report how many records carry sequence identifiers in a
/// legacy Smitten format (V0/V1, or a mix) rather than the current standard (V2).
///
/// Older identifiers such as `chr1:100-200` (V1) or `chr1_100_200` (V0) are still
/// parseable, and DisCoord can verify their coordinates against a reference — so
/// this is not an error.  But Dfam stores identifiers in the canonical V2 form
/// (`assembly:sequence:start-end_orient`, e.g. `hg38:chr1:1000-2000_+`), so report
/// the count as INFO to nudge curators to upgrade before submission.
///
/// A record is counted once regardless of how many of its rows are legacy; rows
/// whose identifiers are already V2, or that carry no parseable coordinates at all
/// (bare consensus labels), do not contribute.
///
/// Returns file-level diagnostics (the caller prints them with label `FILE`).
pub fn check_seqid_format_summary(records: &[RawDfamRecord]) -> Vec<Diagnostic> {
    /// At most this many record labels are named before eliding the rest.
    const MAX_LISTED: usize = 5;

    let legacy: Vec<String> = records
        .iter()
        .filter(|r| {
            r.sequences.iter().any(|s| {
                matches!(
                    s.inferred_version,
                    Some(IDVersion::V0 | IDVersion::V1 | IDVersion::Mixed)
                )
            })
        })
        .map(|r| r.label())
        .collect();

    if legacy.is_empty() {
        return Vec::new();
    }

    let mut listed = legacy.iter().take(MAX_LISTED).cloned().collect::<Vec<_>>().join(", ");
    if legacy.len() > MAX_LISTED {
        listed.push_str(&format!(", and {} more", legacy.len() - MAX_LISTED));
    }

    let clause = if legacy.len() == 1 {
        "1 record has".to_string()
    } else {
        format!("{} records have", legacy.len())
    };

    vec![info(
        "seqid_legacy_format",
        format!(
            "{} sequence identifiers in a legacy Smitten format (V0/V1): {}.  They \
             are parseable and their coordinates can be verified, but Dfam stores \
             identifiers in the standard V2 form (e.g. 'hg38:chr1:1000-2000_+').  \
             Run the records through DisCoord to rewrite them before submission.",
            clause, listed,
        ),
    )]
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

// ── Cross-record checks ───────────────────────────────────────────────────────

/// File-level check: warn once if any record contains RT/RA/RL fields.
///
/// These are standard Stockholm fields that Dfam does not import; curators
/// who fill them in may believe Dfam will display that metadata, but it won't.
pub fn check_unused_citation_fields(records: &[RawDfamRecord]) -> Vec<Diagnostic> {
    let found: Vec<&str> = ["RT", "RA", "RL"]
        .iter()
        .copied()
        .filter(|tag| records.iter().any(|r| r.gf_has(tag)))
        .collect();
    if found.is_empty() {
        return vec![];
    }
    vec![warn(
        "citation_fields_unused",
        format!(
            "{} field(s) found in one or more records; these standard Stockholm fields \
             are not imported by Dfam and will be silently ignored",
            found.join(", ")
        ),
    )]
}

/// Normalize a raw RD value to a bare DOI string (strips doi.org URL prefix).
pub fn normalize_doi(raw: &str) -> &str {
    let raw = raw.trim();
    raw.strip_prefix("https://doi.org/")
        .or_else(|| raw.strip_prefix("http://doi.org/"))
        .unwrap_or(raw)
}

/// Network check: validate every unique RM (PMID) and RD (DOI) in the records.
///
/// Deduplicates identifiers across records and performs one HTTP HEAD request
/// per unique value.  On network errors the check is skipped conservatively
/// (the identifier is not flagged).  Returns one `(label, Diagnostic)` pair per
/// record that references a bad identifier, so callers can print per-record lines.
pub fn check_citations_network(records: &[RawDfamRecord]) -> Vec<(String, Diagnostic)> {
    // Collect unique identifiers → list of record labels that reference each.
    let mut pmids: HashMap<String, Vec<String>> = HashMap::new();
    let mut dois:  HashMap<String, Vec<String>> = HashMap::new();

    for r in records {
        let label = r.label();
        for v in r.gf_all("RM") {
            let pmid = v.trim().to_string();
            if !pmid.is_empty() {
                pmids.entry(pmid).or_default().push(label.clone());
            }
        }
        for v in r.gf_all("RD") {
            let doi = normalize_doi(v).to_string();
            if !doi.is_empty() {
                dois.entry(doi).or_default().push(label.clone());
            }
        }
    }

    let mut out: Vec<(String, Diagnostic)> = Vec::new();

    for (pmid, labels) in &pmids {
        if !pmid_exists(pmid) {
            let diag = err(
                "pmid_unknown",
                format!("PubMed ID {:?} was not found at pubmed.ncbi.nlm.nih.gov", pmid),
            );
            for label in labels {
                out.push((label.clone(), diag.clone()));
            }
        }
    }

    for (doi, labels) in &dois {
        if !doi_exists(doi) {
            let diag = err(
                "doi_unknown",
                format!("DOI {:?} could not be resolved via doi.org", doi),
            );
            for label in labels {
                out.push((label.clone(), diag.clone()));
            }
        }
    }

    out
}

fn pmid_exists(pmid: &str) -> bool {
    let url = format!("https://pubmed.ncbi.nlm.nih.gov/{}/", pmid);
    match ureq::head(&url)
        .set("User-Agent", "dfam-curator/stk-lint")
        .call()
    {
        Ok(_) => true,
        Err(ureq::Error::Status(404, _)) => false,
        Err(_) => true, // network error — be conservative, don't flag
    }
}

fn doi_exists(doi: &str) -> bool {
    // bioRxiv/medRxiv append a version suffix (e.g. "v2") that doi.org does not register.
    // Strip it before resolving so "10.1101/2024.01.27.577580v2" → "10.1101/2024.01.27.577580".
    let resolved = if doi.starts_with("10.1101/") {
        if let Some(v_pos) = doi.rfind('v') {
            let suffix = &doi[v_pos + 1..];
            if !suffix.is_empty() && suffix.chars().all(|c| c.is_ascii_digit()) {
                &doi[..v_pos]
            } else {
                doi
            }
        } else {
            doi
        }
    } else {
        doi
    };
    let url = format!("https://doi.org/{}", resolved);
    match ureq::head(&url)
        .set("User-Agent", "dfam-curator/stk-lint")
        .call()
    {
        Ok(_) => true,
        Err(ureq::Error::Status(404, _)) => false,
        Err(_) => true,
    }
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
                    "AC {:?} does not match DF/DR + 7 or 9 digits (e.g. DF0000001 or \
                     DF000000001).  AC is assigned by Dfam when a family is released — \
                     do not provide one for new submissions; remove the field.",
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

/// Validate every `AU` line on the record.  A record may carry several `AU` lines, and
/// each line several semicolon-separated authors; ORCID uniqueness is enforced across all
/// of them, since the same curator credited on two lines is still a duplicate.
fn check_au(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    let au_lines = r.gf_all("AU");
    if au_lines.is_empty() {
        return;
    }

    let mut seen_orcids: Vec<String> = Vec::new();

    for au in au_lines {
        let au = au.trim();
        if au.is_empty() {
            d.push(err("empty_field", "AU field is present but empty"));
            continue;
        }
        check_au_line(au, &mut seen_orcids, d);
    }
}

/// Validate one `AU` line's worth of semicolon-separated authors.
///
/// `seen_orcids` accumulates across the record's `AU` lines so duplicates are caught
/// wherever they appear.
fn check_au_line(au: &str, seen_orcids: &mut Vec<String>, d: &mut Vec<Diagnostic>) {
    for raw_token in au.split(';') {
        let token = raw_token.trim();
        if token.is_empty() {
            continue;
        }

        // Commas are not valid in author tokens — catches pure comma lists and mixed separators.
        if token.contains(',') {
            d.push(warn(
                "au_format",
                format!("AU token {:?} contains a comma; use semicolons to separate authors", token),
            ));
            continue;
        }

        // Strip optional ORCID prefix: ORCID:xxxx-xxxx-xxxx-xxxx
        let name = if let Some(rest) = token.strip_prefix("ORCID:") {
            let (orcid, name_part) = rest.split_once(' ').unwrap_or((rest, ""));
            if !is_valid_orcid(orcid) {
                d.push(warn(
                    "au_format",
                    format!(
                        "AU ORCID {:?} does not match the xxxx-xxxx-xxxx-xxxx pattern \
                         (last group may end in X)",
                        orcid
                    ),
                ));
            } else {
                let orcid_upper = orcid.to_uppercase();
                if seen_orcids.contains(&orcid_upper) {
                    d.push(err(
                        "au_format",
                        format!("AU ORCID {:?} appears more than once in this record", orcid),
                    ));
                } else {
                    seen_orcids.push(orcid_upper);
                }
            }
            name_part.trim()
        } else {
            token
        };

        if name.is_empty() {
            d.push(warn("au_format", format!("AU {:?}: ORCID prefix present but no name follows it", token)));
            continue;
        }

        // Colons are not valid in names (the ORCID: prefix was already consumed above).
        if name.contains(':') {
            d.push(warn(
                "au_format",
                format!("AU token {:?} contains a colon; colons are only valid as part of an ORCID: prefix", token),
            ));
            continue;
        }

        // Must have at least two words.
        if !name.contains(' ') {
            d.push(warn(
                "au_format",
                format!("AU token {:?} has no space; expected 'First Last' format", token),
            ));
            continue;
        }

        // Flag any word that contains a '.': catches "R.", "J.", "A.F.A.", "F.", etc.
        if name.split_whitespace().any(|w| w.contains('.')) {
            d.push(warn(
                "au_format",
                format!(
                    "AU token {:?} appears to use abbreviated initials; \
                     use full 'First Last' format",
                    token
                ),
            ));
            continue;
        }

        let words: Vec<&str> = name.split_whitespace().collect();

        // Flag single-letter first word: catches "B McClintock".
        // Middle initials without periods (e.g. "Arian F Smit") are still accepted.
        if words.first().map_or(false, |w| {
            w.len() == 1 && w.chars().next().map_or(false, |c| c.is_ascii_alphabetic())
        }) {
            d.push(err(
                "au_format",
                format!(
                    "AU token {:?} appears to use an abbreviated first name; \
                     use full 'First Last' format",
                    token
                ),
            ));
            continue;
        }

        // Flag old 'Last Initial' style: last word is a single uppercase letter.
        // Catches "Smith J", "Watts E".
        if words.last().map_or(false, |w| {
            w.len() == 1 && w.chars().next().map_or(false, |c| c.is_ascii_uppercase())
        }) {
            d.push(err(
                "au_format",
                format!(
                    "AU token {:?} appears to use 'Last Initial' format; \
                     use 'First Last' format",
                    token
                ),
            ));
        }
    }
}

/// Validate an ORCID identifier string (the bare `xxxx-xxxx-xxxx-xxxx` part,
/// without the `ORCID:` prefix).  The last group may end with `X` (check digit).
fn is_valid_orcid(s: &str) -> bool {
    let parts: Vec<&str> = s.split('-').collect();
    if parts.len() != 4 {
        return false;
    }
    parts[..3].iter().all(|p| p.len() == 4 && p.chars().all(|c| c.is_ascii_digit()))
        && parts[3].len() == 4
        && {
            let (digits, check) = parts[3].split_at(3);
            digits.chars().all(|c| c.is_ascii_digit())
                && (check == "X" || check.chars().all(|c| c.is_ascii_digit()))
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

fn check_ct(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    let Some(ct) = r.gf_first("CT") else { return };
    let value = ct.trim();

    if value.is_empty() {
        d.push(err("empty_field", "CT field is present but empty"));
        return;
    }
    if consensus_type(value).is_none() {
        let known: Vec<&str> = CONSENSUS_TYPES.iter().map(|c| c.word).collect();
        d.push(err(
            "ct_unknown",
            format!(
                "CT value {:?} is not a recognised consensus type (allowed: {})",
                value,
                known.join(", ")
            ),
        ));
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
    let has_rd = r.gf_has("RD");

    if has_rn && !has_rm && !has_rd {
        d.push(warn(
            "ref_block_incomplete",
            "RN is present but neither RM (PubMed ID) nor RD (DOI) found; at least one is required",
        ));
    }
    if (has_rm || has_rd) && !has_rn {
        d.push(err("ref_block_incomplete", "RM/RD is present but no RN (reference number) found"));
    }

    // RM/RD/RT/RA/RL must follow their RN line — only meaningful when RN exists.
    // If RN is absent entirely, ref_block_incomplete already covers it.
    if has_rn {
        const NEEDS_RN: &[&str] = &["RM", "RD", "RT", "RA", "RL"];
        let mut seen_rn = false;
        for (tag, _) in &r.gf {
            if tag == "RN" {
                seen_rn = true;
            } else if NEEDS_RN.contains(&tag.as_str()) && !seen_rn {
                d.push(err(
                    "ref_block_order",
                    format!("{} appears before RN; RN must precede all publication fields", tag),
                ));
                break;
            }
        }
    }

    // Validate RD format: must be a bare DOI (10.xxx/...) or a doi.org URL.
    for rd in r.gf_all("RD") {
        let raw = rd.trim();
        let bare = raw
            .strip_prefix("https://doi.org/")
            .or_else(|| raw.strip_prefix("http://doi.org/"))
            .unwrap_or(raw);
        if !bare.starts_with("10.") || bare.len() < 8 {
            d.push(err(
                "rd_format",
                format!(
                    "RD {:?} does not look like a valid DOI; expected a bare DOI \
                     (e.g. '10.1093/nar/gkl1049') or a doi.org URL \
                     (e.g. 'https://doi.org/10.1093/nar/gkl1049')",
                    raw
                ),
            ));
        }
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

/// Validate the optional `#=GC MM` model mask line (HMMER).
///
/// Each column is `m` (the column lies in a masked range, so hmmbuild emits background
/// frequencies for the corresponding match state) or `.` (unmasked).  The line must be
/// the same width as the rest of the alignment.
fn check_mm(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    let Some(mm) = r.gc.get("MM") else { return };

    if let Some(bad) = first_invalid(mm, MM_VALID) {
        d.push(err(
            "mm_invalid_chars",
            format!(
                "#=GC MM contains invalid character {:?}; only 'm' (masked column) \
                 and '.' (unmasked) are allowed",
                bad
            ),
        ));
    }

    // Width must agree with the alignment.  Prefer RF as the reference width when present,
    // since check_rf has already compared it against every sequence row.
    let (expected, what) = match r.gc.get("RF") {
        Some(rf) => (rf.len(), "#=GC RF"),
        None => match r.sequences.first() {
            Some(row) => (row.aligned_seq.len(), "the alignment"),
            None => return,
        },
    };

    if mm.len() != expected {
        d.push(err(
            "mm_length_mismatch",
            format!(
                "#=GC MM length {} does not match {} length {}",
                mm.len(), what, expected
            ),
        ));
    }
}

fn check_rf_consensus(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    // A hand-curated consensus is not expected to reproduce the called one.
    if rf_is_handcurated(r) {
        return;
    }

    let rf = match r.gc.get("RF") {
        Some(rf) if !rf.is_empty() => rf,
        _ => return,
    };
    if r.sequences.is_empty() {
        return;
    }

    // STK files may use any Stockholm gap character; the consensus builder uses '-'.
    let converted: Vec<Vec<u8>> = r.sequences.iter()
        .map(|row| row.aligned_seq.bytes().map(|b| if is_gap(b) { b'-' } else { b }).collect())
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

/// Reject HMMER's interleaved *block* format.
///
/// Stockholm permits an alignment to be split into blocks separated by blank lines, with
/// every sequence appearing once per block.  Dfam does not accept this: each sequence must
/// appear complete on a single line.  The parser records the blank line that separated two
/// groups of sequence rows.
fn check_block_format(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    if let Some(line) = r.block_separator_line {
        d.push(err(
            "block_format",
            format!(
                "blank line at line {} splits the sequence section into blocks; \
                 Stockholm interleaved block format is not supported by Dfam — \
                 each sequence must appear complete on its own single line",
                line
            ),
        ));
    }
}

fn check_sequences(r: &RawDfamRecord, d: &mut Vec<Diagnostic>) {
    let mut first_len: Option<usize> = None;

    for row in &r.sequences {
        // Character validation: report the first bad character per sequence, and
        // separately the first non-'.' gap character (legal Stockholm, but not the
        // Dfam convention).
        let mut nonstandard_gap: Option<char> = None;

        for (pos, &b) in row.aligned_seq.as_bytes().iter().enumerate() {
            if is_gap(b) {
                if b != DFAM_GAP && nonstandard_gap.is_none() {
                    nonstandard_gap = Some(b as char);
                }
                continue;
            }
            if !SEQ_VALID.contains(&b) {
                d.push(err(
                    "seq_invalid_chars",
                    format!(
                        "sequence {:?} contains invalid character {:?} at position {}",
                        row.original_id,
                        b as char,
                        pos + 1,
                    ),
                ));
                break;
            }
        }

        if let Some(gap) = nonstandard_gap {
            d.push(warn(
                "seq_nonstandard_gap",
                format!(
                    "sequence {:?} uses {:?} as a gap character; Stockholm permits \
                     '-', '.', '_' and '~', but Dfam has standardized on '.'",
                    row.original_id, gap,
                ),
            ));
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

    // Dfam consumes only RF and MM; other Stockholm #=GC annotations (SS_cons, PP_cons, …)
    // are carried through the file but ignored on import.
    let mut unknown: Vec<&str> = r.gc.keys()
        .map(String::as_str)
        .filter(|tag| !KNOWN_GC_TAGS.contains(tag))
        .collect();
    unknown.sort_unstable();
    for tag in unknown {
        d.push(info(
            "unknown_gc_tag",
            format!(
                "unrecognised #=GC annotation {:?}; Dfam uses RF (consensus) and \
                 MM (model mask) — others are not imported",
                tag
            ),
        ));
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
            &[("DE","A test family"),("AU","John Smith"),("TP","Interspersed_Repeat;Unknown"),
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
    fn seq_dash_gap_warns() {
        // '-' is a legal Stockholm gap, but Dfam standardizes on '.': WARN, not ERROR.
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","1")],
            Some("ACGT"),
            &[("s1","AC-T")],
        );
        let diags = lint_record(&r, None);
        assert!(!has_check(&diags, "seq_invalid_chars"), "{:?}", diags);
        let g = diags.iter().find(|d| d.check == "seq_nonstandard_gap").expect("expected warn");
        assert_eq!(g.severity, Severity::Warn);
        assert!(g.message.contains("'-'"), "{}", g.message);
    }

    #[test]
    fn seq_tilde_and_underscore_gaps_warn() {
        for gap in ["AC~T", "AC_T"] {
            let r = make_record(
                &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","1")],
                Some("ACGT"),
                &[("s1", gap)],
            );
            let diags = lint_record(&r, None);
            assert!(!has_check(&diags, "seq_invalid_chars"), "{}: {:?}", gap, diags);
            assert!(has_check(&diags, "seq_nonstandard_gap"), "{}: {:?}", gap, diags);
        }
    }

    #[test]
    fn seq_dot_gap_is_clean() {
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","1")],
            Some("AC.T"),
            &[("s1","AC.T")],
        );
        let diags = lint_record(&r, None);
        assert!(!has_check(&diags, "seq_nonstandard_gap"), "{:?}", diags);
        assert!(!has_check(&diags, "seq_invalid_chars"), "{:?}", diags);
    }

    #[test]
    fn seq_truly_invalid_char_is_still_error() {
        // '*' is not a residue and not a gap.
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","1")],
            Some("ACGT"),
            &[("s1","AC*T")],
        );
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "seq_invalid_chars"), "{:?}", diags);
    }

    #[test]
    fn nonstandard_gap_reported_once_per_sequence() {
        // Many gaps in one row produce a single warning, not one per column.
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","1")],
            Some("ACGTACGT"),
            &[("s1","A--T--GT")],
        );
        let diags = lint_record(&r, None);
        let n = diags.iter().filter(|d| d.check == "seq_nonstandard_gap").count();
        assert_eq!(n, 1, "{:?}", diags);
    }

    #[test]
    fn dash_gapped_alignment_still_consensus_checked() {
        // Gaps written as '-' must be understood as gaps by the consensus caller,
        // not treated as residues.
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","4")],
            Some("AC.T"),
            &[("s1","AC-T"),("s2","AC-T"),("s3","AC-T"),("s4","AC-T")],
        );
        let diags = lint_record(&r, None);
        assert!(!has_check(&diags, "rf_consensus_mismatch"), "{:?}", diags);
    }

    #[test]
    fn block_format_is_error() {
        let stk = "\
# STOCKHOLM 1.0\n\
#=GF DE    Test family\n\
#=GF AU    Barbara McClintock\n\
#=GF TP    X\n\
#=GF OC    Mus musculus\n\
#=GF SQ    2\n\
s1         ACGT\n\
s2         ACGT\n\
\n\
s1         TTTT\n\
s2         TTTT\n\
//\n";
        let diags = lint_stk(stk);
        assert!(has_check(&diags, "block_format"), "{:?}", diags);
        let d = diags.iter().find(|d| d.check == "block_format").unwrap();
        assert_eq!(d.severity, Severity::Error);
        assert!(d.message.contains("line 9"), "{}", d.message);
    }

    #[test]
    fn trailing_blank_line_is_not_block_format() {
        // A blank line before '//' separates nothing and must not be flagged.
        let stk = "\
# STOCKHOLM 1.0\n\
#=GF DE    Test family\n\
#=GF AU    Barbara McClintock\n\
#=GF TP    X\n\
#=GF OC    Mus musculus\n\
#=GF SQ    1\n\
#=GC RF    ACGT\n\
s1         ACGT\n\
\n\
//\n";
        let diags = lint_stk(stk);
        assert!(!has_check(&diags, "block_format"), "{:?}", diags);
    }

    #[test]
    fn blank_line_among_header_fields_is_not_block_format() {
        let stk = "\
# STOCKHOLM 1.0\n\
#=GF DE    Test family\n\
\n\
#=GF AU    Barbara McClintock\n\
#=GF TP    X\n\
#=GF OC    Mus musculus\n\
#=GF SQ    1\n\
#=GC RF    ACGT\n\
s1         ACGT\n\
//\n";
        let diags = lint_stk(stk);
        assert!(!has_check(&diags, "block_format"), "{:?}", diags);
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
#=AU Barbara McClintock, Roy Britten\n\
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
#AU Barbara McClintock, Roy Britten\n\
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

    /// Parse a one-record Stockholm string and lint it.
    fn lint_stk(stk: &str) -> Vec<Diagnostic> {
        use std::io::Cursor;
        use crate::dfam::record::iter_records;
        let records: Vec<_> = iter_records(Cursor::new(stk))
            .collect::<Result<_, _>>()
            .unwrap();
        lint_record(&records[0], None)
    }

    /// A minimal valid record with the given extra annotation lines spliced in.
    fn stk_with(extra: &str) -> String {
        format!(
            "# STOCKHOLM 1.0\n\
             #=GF DE    Test family\n\
             #=GF AU    Barbara McClintock\n\
             #=GF TP    X\n\
             #=GF OC    Mus musculus\n\
             #=GF SQ    1\n\
             #=GC RF    ACGTACGT\n\
             {}\
             s1         ACGTACGT\n\
             //\n",
            extra
        )
    }

    #[test]
    fn mm_valid_is_accepted() {
        let diags = lint_stk(&stk_with("#=GC MM    ..mmmm..\n"));
        assert!(!has_check(&diags, "mm_invalid_chars"), "{:?}", diags);
        assert!(!has_check(&diags, "mm_length_mismatch"), "{:?}", diags);
        assert!(!has_check(&diags, "unknown_gc_tag"), "{:?}", diags);
    }

    #[test]
    fn mm_absent_is_fine() {
        let diags = lint_stk(&stk_with(""));
        assert!(!has_check(&diags, "mm_invalid_chars"), "{:?}", diags);
        assert!(!has_check(&diags, "mm_length_mismatch"), "{:?}", diags);
    }

    #[test]
    fn mm_invalid_char_is_error() {
        // 'M' uppercase and 'x' are not valid mask characters; only 'm' and '.'.
        let diags = lint_stk(&stk_with("#=GC MM    ..MMMM..\n"));
        assert!(has_check(&diags, "mm_invalid_chars"), "{:?}", diags);
        let d = diags.iter().find(|d| d.check == "mm_invalid_chars").unwrap();
        assert_eq!(d.severity, Severity::Error);
    }

    #[test]
    fn mm_length_mismatch_is_error() {
        let diags = lint_stk(&stk_with("#=GC MM    ..mm\n"));
        assert!(has_check(&diags, "mm_length_mismatch"), "{:?}", diags);
    }

    #[test]
    fn unknown_gc_tag_is_info() {
        let diags = lint_stk(&stk_with("#=GC SS_cons    ........\n"));
        assert!(has_check(&diags, "unknown_gc_tag"), "{:?}", diags);
        let d = diags.iter().find(|d| d.check == "unknown_gc_tag").unwrap();
        assert_eq!(d.severity, Severity::Info);
    }

    #[test]
    fn au_second_line_is_validated() {
        // The first AU line is well-formed; the second uses the old 'Last Initial' style.
        let r = make_record(
            &[("DE","x"),("AU","Barbara McClintock"),("AU","Smith J"),
              ("TP","x"),("OC","x"),("SQ","0")],
            Some(""), &[],
        );
        let diags = lint_record(&r, None);
        let au: Vec<&Diagnostic> = diags.iter().filter(|d| d.check == "au_format").collect();
        assert_eq!(au.len(), 1, "{:?}", diags);
        assert_eq!(au[0].severity, Severity::Error);
        assert!(au[0].message.contains("Smith J"), "{}", au[0].message);
    }

    #[test]
    fn au_duplicate_orcid_across_lines_is_error() {
        let r = make_record(
            &[("DE","x"),
              ("AU","ORCID:0000-0001-2345-6789 Barbara McClintock"),
              ("AU","ORCID:0000-0001-2345-6789 Roy Britten"),
              ("TP","x"),("OC","x"),("SQ","0")],
            Some(""), &[],
        );
        let diags = lint_record(&r, None);
        let dup = diags.iter().find(|d| d.check == "au_format").expect("expected au_format");
        assert_eq!(dup.severity, Severity::Error);
        assert!(dup.message.contains("more than once"), "{}", dup.message);
    }

    #[test]
    fn au_multiple_valid_lines_are_clean() {
        let r = make_record(
            &[("DE","x"),
              ("AU","ORCID:0000-0001-2345-6789 Barbara McClintock; Roy Britten"),
              ("AU","ORCID:0000-0002-1825-0097 Josiah Carberry"),
              ("TP","x"),("OC","x"),("SQ","0")],
            Some(""), &[],
        );
        let diags = lint_record(&r, None);
        assert!(!has_check(&diags, "au_format"), "{:?}", diags);
    }

    #[test]
    fn orcid_summary_counts_records_missing_orcids() {
        let no_orcid = make_record(
            &[("DE","x"),("AU","Barbara McClintock"),("TP","x"),("OC","x"),("SQ","0")],
            Some(""), &[],
        );
        let with_orcid = make_record(
            &[("DE","x"),("AU","ORCID:0000-0001-2345-6789 Barbara McClintock"),
              ("TP","x"),("OC","x"),("SQ","0")],
            Some(""), &[],
        );
        let diags = check_orcid_summary(&[no_orcid, with_orcid]);
        assert_eq!(diags.len(), 1);
        assert_eq!(diags[0].check, "orcid_missing");
        assert_eq!(diags[0].severity, Severity::Info);
        assert!(diags[0].message.starts_with("1 of 2 records"), "{}", diags[0].message);
    }

    #[test]
    fn orcid_summary_flags_record_where_only_some_authors_have_orcids() {
        let r = make_record(
            &[("DE","x"),("AU","ORCID:0000-0001-2345-6789 Barbara McClintock; Roy Britten"),
              ("TP","x"),("OC","x"),("SQ","0")],
            Some(""), &[],
        );
        let diags = check_orcid_summary(&[r]);
        assert_eq!(diags.len(), 1);
        assert!(diags[0].message.starts_with("1 of 1 records"), "{}", diags[0].message);
    }

    #[test]
    fn orcid_summary_silent_when_all_authors_credited() {
        // Multiple AU lines, every author carrying an ORCID.
        let r = make_record(
            &[("DE","x"),
              ("AU","ORCID:0000-0001-2345-6789 Barbara McClintock"),
              ("AU","ORCID:0000-0002-1825-0097 Josiah Carberry"),
              ("TP","x"),("OC","x"),("SQ","0")],
            Some(""), &[],
        );
        assert!(check_orcid_summary(&[r]).is_empty());
    }

    #[test]
    fn orcid_summary_ignores_records_without_au() {
        // Missing AU is already an error; it should not inflate the ORCID denominator.
        let no_au = make_record(
            &[("DE","x"),("TP","x"),("OC","x"),("SQ","0")],
            Some(""), &[],
        );
        let no_orcid = make_record(
            &[("DE","x"),("AU","Barbara McClintock"),("TP","x"),("OC","x"),("SQ","0")],
            Some(""), &[],
        );
        let diags = check_orcid_summary(&[no_au, no_orcid]);
        assert!(diags[0].message.starts_with("1 of 1 records"), "{}", diags[0].message);
    }

    #[test]
    fn ac_summary_counts_update_records() {
        let with_ac = make_record(
            &[("AC","DR000000001"),("ID","famA"),("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","0")],
            Some(""), &[],
        );
        let without_ac = make_record(
            &[("ID","famB"),("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","0")],
            Some(""), &[],
        );
        let diags = check_ac_summary(&[with_ac, without_ac]);
        assert_eq!(diags.len(), 1);
        assert_eq!(diags[0].check, "ac_update_records");
        assert_eq!(diags[0].severity, Severity::Info);
        assert!(diags[0].message.starts_with("1 record carries an AC"), "{}", diags[0].message);
        // The accession itself is named, not just the record label.
        assert!(diags[0].message.contains("famA = DR000000001")
                || diags[0].message.contains("DR000000001"), "{}", diags[0].message);
    }

    #[test]
    fn ac_summary_silent_when_no_ac() {
        let r = make_record(
            &[("ID","famA"),("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","0")],
            Some(""), &[],
        );
        assert!(check_ac_summary(&[r]).is_empty());
    }

    #[test]
    fn ac_summary_ignores_empty_ac() {
        // An empty AC is already reported per-record as empty_field; it is not an update.
        let r = make_record(
            &[("AC","  "),("ID","famA"),("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","0")],
            Some(""), &[],
        );
        assert!(check_ac_summary(&[r]).is_empty());
    }

    /// Build a record whose sequence rows are parsed by Smitten (so `inferred_version`
    /// is populated), unlike `make_record` which uses `new_raw`.
    fn record_with_seq_ids(ids: &[&str]) -> RawDfamRecord {
        let mut r = RawDfamRecord::default();
        r.record_num = 1;
        r.header = "# STOCKHOLM 1.0".to_string();
        r.terminated = true;
        r.gf.push(("ID".to_string(), "famA".to_string()));
        for id in ids {
            r.sequences.push(SeqRow::from_name_seq(id, "ACGT"));
        }
        r
    }

    #[test]
    fn seqid_summary_flags_legacy_v1_identifier() {
        // `chr1:100-200` (no strand suffix) is V1: parseable, but not the V2 standard.
        let r = record_with_seq_ids(&["chr1:100-200"]);
        assert_eq!(r.sequences[0].inferred_version, Some(IDVersion::V1));
        let diags = check_seqid_format_summary(&[r]);
        assert_eq!(diags.len(), 1);
        assert_eq!(diags[0].check, "seqid_legacy_format");
        assert_eq!(diags[0].severity, Severity::Info);
        assert!(diags[0].message.starts_with("1 record has"), "{}", diags[0].message);
    }

    #[test]
    fn seqid_summary_flags_legacy_v0_identifier() {
        // `chr1_100_200` uses underscore separators: V0.
        let r = record_with_seq_ids(&["chr1_100_200"]);
        assert_eq!(r.sequences[0].inferred_version, Some(IDVersion::V0));
        let diags = check_seqid_format_summary(&[r]);
        assert_eq!(diags.len(), 1);
        assert_eq!(diags[0].check, "seqid_legacy_format");
    }

    #[test]
    fn seqid_summary_silent_for_v2_identifiers() {
        // `hg38:chr1:1000-2000_+` is the standard V2 form — no diagnostic.
        let r = record_with_seq_ids(&["hg38:chr1:1000-2000_+"]);
        assert_eq!(r.sequences[0].inferred_version, Some(IDVersion::V2));
        assert!(check_seqid_format_summary(&[r]).is_empty());
    }

    #[test]
    fn seqid_summary_counts_each_record_once() {
        // A record with both a legacy and a V2 row is counted a single time.
        let mixed_rows = record_with_seq_ids(&["chr1:100-200", "hg38:chr1:1000-2000_+"]);
        let clean = record_with_seq_ids(&["hg38:chr2:5-9_-"]);
        let diags = check_seqid_format_summary(&[mixed_rows, clean]);
        assert_eq!(diags.len(), 1);
        assert!(diags[0].message.starts_with("1 record has"), "{}", diags[0].message);
    }

    #[test]
    fn ct_handbuilt_skips_consensus_check() {
        // RF disagrees with the called consensus, but CT declares it hand-built.
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","4"),("CT","handbuilt")],
            Some("TTTT"),
            &[("s1","ACGT"),("s2","ACGT"),("s3","ACGT"),("s4","ACGT")],
        );
        let diags = lint_record(&r, None);
        assert!(!has_check(&diags, "rf_consensus_mismatch"), "{:?}", diags);
        assert!(!has_check(&diags, "unknown_gf_tag"), "{:?}", diags);
    }

    #[test]
    fn ct_unknown_word_is_error() {
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","1"),("CT","artisanal")],
            Some("ACGT"),
            &[("s1","ACGT")],
        );
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "ct_unknown"), "{:?}", diags);
        // An unrecognised CT does not suppress the consensus check.
        let d = diags.iter().find(|d| d.check == "ct_unknown").unwrap();
        assert_eq!(d.severity, Severity::Error);
        assert!(d.message.contains("handbuilt"), "{}", d.message);
    }

    #[test]
    fn ct_empty_is_error() {
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","1"),("CT","")],
            Some("ACGT"),
            &[("s1","ACGT")],
        );
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "empty_field"), "{:?}", diags);
    }

    #[test]
    fn ct_absent_still_checks_consensus() {
        let r = make_record(
            &[("DE","x"),("AU","x"),("TP","x"),("OC","x"),("SQ","4")],
            Some("TTTT"),
            &[("s1","ACGT"),("s2","ACGT"),("s3","ACGT"),("s4","ACGT")],
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

    // ── AU format tests ───────────────────────────────────────────────────────

    fn au_warns(au: &str) -> bool {
        let r = make_record(
            &[("DE","x"),("AU", au),("TP","x"),("OC","x"),("SQ","0")],
            Some(""),
            &[],
        );
        lint_record(&r, None).iter().any(|d| d.check == "au_format")
    }

    #[test]
    fn au_full_name_is_clean() {
        assert!(!au_warns("Barbara McClintock"));
    }

    #[test]
    fn au_multiple_full_names_semicolon_is_clean() {
        assert!(!au_warns("Barbara McClintock; Pita Enriquez-Lopez"));
    }

    #[test]
    fn au_orcid_full_name_is_clean() {
        assert!(!au_warns("ORCID:0000-0001-2345-6789 Barbara McClintock"));
    }

    #[test]
    fn au_orcid_with_x_check_digit_is_clean() {
        assert!(!au_warns("ORCID:0000-0001-2345-678X Barbara McClintock"));
    }

    #[test]
    fn au_last_initial_warns() {
        assert!(au_warns("Smith J"));
    }

    #[test]
    fn au_abbreviated_initial_dot_warns() {
        assert!(au_warns("E. Watts"));
    }

    #[test]
    fn au_middle_initial_dot_warns() {
        assert!(au_warns("Elena J. Watts"));
    }

    #[test]
    fn au_chained_initials_warns() {
        assert!(au_warns("A.F.A. Smit"));
    }

    #[test]
    fn au_single_letter_first_name_warns() {
        assert!(au_warns("B McClintock"));
    }

    #[test]
    fn au_middle_initial_no_dot_is_clean() {
        // Single-letter middle initial without period is acceptable.
        assert!(!au_warns("Arian F Smit"));
    }

    #[test]
    fn au_no_space_warns() {
        assert!(au_warns("Watts"));
    }

    #[test]
    fn au_comma_separator_warns() {
        assert!(au_warns("Barbara McClintock, Roy Britten"));
    }

    #[test]
    fn au_mixed_separator_warns() {
        assert!(au_warns("Barbara McClintock; Roy Britten, John Smith"));
    }

    #[test]
    fn au_colon_in_name_warns() {
        assert!(au_warns("Barbara: McClintock"));
    }

    #[test]
    fn au_duplicate_orcid_warns() {
        assert!(au_warns(
            "ORCID:0000-0001-2345-6789 Barbara McClintock; \
             ORCID:0000-0001-2345-6789 Roy Britten"
        ));
    }

    #[test]
    fn au_distinct_orcids_clean() {
        assert!(!au_warns(
            "ORCID:0000-0001-2345-6789 Barbara McClintock; \
             ORCID:0000-0009-8765-4321 Roy Britten"
        ));
    }

    #[test]
    fn au_bad_orcid_format_warns() {
        assert!(au_warns("ORCID:0000-0001-2345 Barbara McClintock"));
    }

    #[test]
    fn au_orcid_no_name_warns() {
        assert!(au_warns("ORCID:0000-0001-2345-6789"));
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

    #[test]
    fn rm_before_rn_is_error() {
        let r = make_record(
            &[("AU","x"),("RM","12345"),("RN","[1]")],
            None,
            &[],
        );
        let diags = lint_record(&r, None);
        assert!(has_check(&diags, "ref_block_order"), "{:?}", diags);
    }

    #[test]
    fn rn_before_rm_is_ok() {
        let r = make_record(
            &[("AU","x"),("RN","[1]"),("RM","12345")],
            None,
            &[],
        );
        let diags = lint_record(&r, None);
        assert!(!has_check(&diags, "ref_block_order"), "{:?}", diags);
    }
}
