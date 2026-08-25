/// Read Repbase's IG-derived family record (metadata + consensus).
///
/// The format carries EMBL-style content, but every metadata line is prefixed with
/// `;` (the IG comment marker).  A two-letter tag follows the `;` (`ID`, `DE`, `KW`,
/// `OS`, `OC`, `RN`/`RA`/`RT`/`RL`, `CC`, `SQ`, …); `XX` lines are separators.  After
/// the `;SQ` line the consensus begins: a bare identifier line, then the wrapped
/// (ungapped) consensus sequence to end of file.
///
/// ```text
/// ;ID   Mariner-N5_CyaStr DNA   ; PLN   ; 6225 BP
/// ;XX
/// ;DE   DNA transposon from the Cyathus striatus genome, consensus.
/// ;KW   Mariner/Tc1; DNA transposon; Transposable Element; nonautonomous;
/// ;KW   Mariner-N5_CyaStr.
/// ;OS   Cyathus striatus
/// ;RN   [1]  ()
/// ;RA   Bao,W.
/// ;RT   DNA transposons from the Cyathus striatus genome.
/// ;RL   Direct Submission to RR (8-Jul-2026)
/// ;CC   ~96% identical to consensus.
/// ;SQ   Sequence 6225 BP; 1923 A; ...
/// Mariner-N5_CyaStr
/// CTGGATAATTTCGACC....
/// ```
///
/// This module only *parses* the record into [`IgFamilyRecord`] — it does not
/// translate any field to Stockholm (that is a later phase).
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// One reference block (`RN`/`RA`/`RT`/`RL`) from a family record.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct IgReference {
    /// `RN` value, e.g. `[1]  ()`.
    pub number: Option<String>,
    /// `RA` — authors.
    pub authors: Option<String>,
    /// `RT` — title.
    pub title: Option<String>,
    /// `RL` — location / journal.
    pub location: Option<String>,
}

/// One `FT` feature-table entry (e.g. a `CDS`) from a family record.
///
/// Content is held exactly as Repbase wrote it — the location string is *not*
/// parsed or validated here, and the translation is *not* checked against the
/// consensus.  Both happen in the translation phase.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct IgFeature {
    /// Feature key, e.g. `CDS`.
    pub key: String,
    /// Location as written, with any line-wrapping rejoined — e.g. `490..4569`
    /// or `join(312..602,1346..2053,2317..2817)`.  Note that Repbase encodes the
    /// minus strand with *descending* coordinates (`join(5133..5095,4587..4357)`)
    /// rather than EMBL's `complement()`.
    pub location: String,
    /// Qualifiers in document order, as `(name, value)` with the surrounding
    /// quotes removed — e.g. `("product", "Gypsy-13_AnMou-I_1p")`.
    pub qualifiers: Vec<(String, String)>,
}

/// A parsed IG/Repbase family record.
///
/// Fields hold the raw parsed content; interpretation (KW→class, OS→OC, DE→CC, …)
/// happens in the translation phase, not here.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct IgFamilyRecord {
    /// Identifier — the first whitespace token of the `ID` line
    /// (e.g. `Mariner-N5_CyaStr`).
    pub id: String,
    /// The full `ID` line content after the tag (molecule type, division, length).
    pub id_line: String,
    /// `DE` description, wrapped lines joined with a space.
    pub description: Option<String>,
    /// `AC` accession as written (Repbase often uses `.` to mean "none").
    pub accession: Option<String>,
    /// `DT` date line(s).
    pub date: Option<String>,
    /// `KW` keywords: all lines joined, split on `;`, trimmed, trailing `.` removed.
    pub keywords: Vec<String>,
    /// `OS` organism (species name).
    pub organism: Option<String>,
    /// `OC` lineage tokens (joined, split on `;`, trimmed, trailing `.` removed).
    pub oc_lineage: Vec<String>,
    /// Reference blocks in document order.
    pub references: Vec<IgReference>,
    /// `CC` comment lines (one entry per source line).
    pub comments: Vec<String>,
    /// `FT` feature-table entries in document order (empty for the many Repbase
    /// families with no coding capacity).
    pub features: Vec<IgFeature>,
    /// `SQ` summary line content, if present.
    pub sq_summary: Option<String>,
    /// The bare identifier line that introduces the consensus (after `SQ`).
    pub consensus_name: Option<String>,
    /// Concatenated consensus sequence (ungapped).
    pub consensus: Vec<u8>,
}

/// Column at which an `FT` location/qualifier begins, for records that carry no
/// `FH` header to name it.  Repbase writes 19 throughout RepBase31.06.
const DEFAULT_LOC_COL: usize = 19;

/// Split an `FT` body into `(key, content)` at the location column.
///
/// `body` has had the leading `;` removed, so it carries EMBL's layout verbatim:
/// `FT   CDS           490..4569` — tag at `0..2`, key field at `2..loc_col`,
/// location or qualifier from `loc_col`.  An empty key field means the line
/// continues the feature above it.
fn split_ft_line(body: &str, loc_col: usize) -> (String, String) {
    let key_end = loc_col.min(body.len());
    let key = body.get(2..key_end).unwrap_or("").trim().to_string();
    let content = body.get(loc_col..).unwrap_or("").trim().to_string();
    (key, content)
}

/// Fold a wrapped `FT` continuation line into the feature under construction.
///
/// A line beginning `/` opens a new qualifier.  Otherwise it continues the
/// location (while no qualifier has appeared yet) or the value of the most recent
/// qualifier.  EMBL wraps `/translation` mid-token, so that one rejoins with no
/// separator; free-text qualifiers rejoin with a space.
fn push_ft_continuation(f: &mut IgFeature, content: &str) {
    if let Some(q) = content.strip_prefix('/') {
        let (name, value) = match q.split_once('=') {
            Some((n, v)) => (
                n.trim().to_string(),
                v.trim_start().trim_start_matches('"').to_string(),
            ),
            None => (q.trim().to_string(), String::new()),
        };
        f.qualifiers.push((name, value));
    } else if let Some((name, value)) = f.qualifiers.last_mut() {
        if name != "translation" && !value.is_empty() {
            value.push(' ');
        }
        value.push_str(content);
    } else {
        f.location.push_str(content);
    }
}

/// Drop the closing quote from each qualifier value (the opening one is removed as
/// the qualifier is parsed).
fn finish_feature(mut f: IgFeature) -> IgFeature {
    for (_, v) in f.qualifiers.iter_mut() {
        if v.ends_with('"') {
            v.pop();
        }
    }
    f
}

/// Parse an IG/Repbase family record file into an [`IgFamilyRecord`].
pub fn read(path: &Path) -> io::Result<IgFamilyRecord> {
    let f = BufReader::new(std::fs::File::open(path)?);
    let mut rec = IgFamilyRecord::default();

    // Wrapped multi-line fields are accumulated raw, then post-processed.
    let mut de_parts: Vec<String> = Vec::new();
    let mut kw_parts: Vec<String> = Vec::new();
    let mut oc_parts: Vec<String> = Vec::new();
    let mut cur_ref: Option<IgReference> = None;
    let mut cur_feat: Option<IgFeature> = None;

    // Where an FT location/qualifier begins.  The `FH` header is self-describing
    // ("FH   Key           Location/Qualifiers"), so the column is read from it
    // rather than assumed; `DEFAULT_LOC_COL` only applies to records with no `FH`.
    let mut loc_col = DEFAULT_LOC_COL;

    // `after_sq` is set by the SQ line; the next bare line is the consensus name.
    let mut after_sq = false;
    let mut in_sequence = false;

    for line in f.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        if in_sequence {
            // Everything past the consensus name is sequence data.
            rec.consensus.extend(trimmed.bytes());
            continue;
        }

        if let Some(body) = trimmed.strip_prefix(';') {
            if body.len() < 2 {
                continue; // bare ';' separator
            }
            let tag = &body[..2];
            let content = body[2..].trim();

            match tag {
                "XX" => {} // separator
                "FH" => {
                    // "FH   Key           Location/Qualifiers" names its own columns.
                    if let Some(i) = body.find("Location") {
                        loc_col = i;
                    }
                }
                "FT" => {
                    let (key, content) = split_ft_line(body, loc_col);
                    if key.is_empty() {
                        // Continuation of the feature already under construction.
                        if let Some(f) = cur_feat.as_mut() {
                            push_ft_continuation(f, &content);
                        }
                    } else {
                        if let Some(f) = cur_feat.take() {
                            rec.features.push(finish_feature(f));
                        }
                        cur_feat = Some(IgFeature {
                            key,
                            location: content,
                            qualifiers: Vec::new(),
                        });
                    }
                }
                "ID" => {
                    rec.id_line = content.to_string();
                    rec.id = content
                        .split_whitespace()
                        .next()
                        .unwrap_or("")
                        .to_string();
                }
                "DE" => de_parts.push(content.to_string()),
                "AC" => rec.accession = Some(content.to_string()),
                "DT" => rec.date = Some(content.to_string()),
                "KW" => kw_parts.push(content.to_string()),
                "OS" => rec.organism = Some(content.to_string()),
                "OC" => oc_parts.push(content.to_string()),
                "RN" => {
                    if let Some(r) = cur_ref.take() {
                        rec.references.push(r);
                    }
                    cur_ref = Some(IgReference {
                        number: Some(content.to_string()),
                        ..Default::default()
                    });
                }
                "RA" => cur_ref.get_or_insert_with(Default::default).authors = Some(content.to_string()),
                "RT" => cur_ref.get_or_insert_with(Default::default).title = Some(content.to_string()),
                "RL" => cur_ref.get_or_insert_with(Default::default).location = Some(content.to_string()),
                "CC" => rec.comments.push(content.to_string()),
                "SQ" => {
                    rec.sq_summary = Some(content.to_string());
                    after_sq = true;
                }
                _ => {} // unknown tag (e.g. DR) — ignored for now
            }
            continue;
        }

        // A bare (non-';') line: only expected as the consensus name after SQ.
        if after_sq {
            rec.consensus_name = Some(trimmed.to_string());
            in_sequence = true;
        }
        // Otherwise it's unexpected content before SQ; ignore it.
    }

    if let Some(r) = cur_ref.take() {
        rec.references.push(r);
    }
    if let Some(f) = cur_feat.take() {
        rec.features.push(finish_feature(f));
    }

    rec.description = join_nonempty(&de_parts, " ");
    rec.keywords = split_semicolon_tokens(&kw_parts);
    rec.oc_lineage = split_semicolon_tokens(&oc_parts);

    Ok(rec)
}

/// Join wrapped field lines with `sep`, returning `None` if all were empty.
fn join_nonempty(parts: &[String], sep: &str) -> Option<String> {
    let joined = parts
        .iter()
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect::<Vec<_>>()
        .join(sep);
    if joined.is_empty() { None } else { Some(joined) }
}

/// Split each wrapped line on `;`, trim each token, drop a trailing `.`, and
/// discard empties.  Used for `KW` and `OC`.
///
/// Lines are split independently rather than joined first: EMBL wraps these
/// `;`-separated lists only at token boundaries, and Repbase sometimes writes a
/// standalone continuation line with no trailing `;` (e.g. the species name as the
/// first `OC` entry), so joining with a space would fuse it onto the next token.
fn split_semicolon_tokens(parts: &[String]) -> Vec<String> {
    parts
        .iter()
        .flat_map(|line| line.split(';'))
        .map(|t| t.trim().trim_end_matches('.').trim())
        .filter(|t| !t.is_empty())
        .map(|t| t.to_string())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    const SAMPLE: &str = "\
;ID   Mariner-N5_CyaStr DNA   ; PLN   ; 6225 BP
;XX
;DE   DNA transposon from the Cyathus striatus genome, consensus.
;XX
;AC   .
;XX
;DT   18-MAR-2026 (Created)
;XX
;KW   Mariner/Tc1; DNA transposon; Transposable Element; nonautonomous;
;KW   Mariner-N5_CyaStr.
;XX
;OS   Cyathus striatus
;XX
;OC   Cyathus striatus
;OC   Eukaryota; Fungi; Dikarya; Basidiomycota; Agaricomycotina;
;OC   Nidulariaceae; Cyathus.
;XX
;RN   [1]  ()
;RA   Bao,W.
;RT   DNA transposons from the Cyathus striatus genome.
;RL   Direct Submission to RR (8-Jul-2026)
;XX
;CC   ~96% identical to consensus.
;XX
;SQ   Sequence 6225 BP; 1923 A; 1173 C; 1331 G; 1797 T; 1 other;
Mariner-N5_CyaStr
CTGGATAATTTCGACCAAGTGT
TCGGACCCTTTCAGGTGTGACG
";

    fn parse(body: &str, name: &str) -> IgFamilyRecord {
        let path = std::env::temp_dir().join(name);
        std::fs::File::create(&path).unwrap().write_all(body.as_bytes()).unwrap();
        read(&path).unwrap()
    }

    #[test]
    fn id_is_first_token() {
        let r = parse(SAMPLE, "ig_fam_id.ig");
        assert_eq!(r.id, "Mariner-N5_CyaStr");
        assert_eq!(r.id_line, "Mariner-N5_CyaStr DNA   ; PLN   ; 6225 BP");
    }

    #[test]
    fn description_and_accession() {
        let r = parse(SAMPLE, "ig_fam_de.ig");
        assert_eq!(
            r.description.as_deref(),
            Some("DNA transposon from the Cyathus striatus genome, consensus.")
        );
        assert_eq!(r.accession.as_deref(), Some("."));
    }

    #[test]
    fn keywords_split_across_lines_and_dot_stripped() {
        let r = parse(SAMPLE, "ig_fam_kw.ig");
        assert_eq!(
            r.keywords,
            vec![
                "Mariner/Tc1",
                "DNA transposon",
                "Transposable Element",
                "nonautonomous",
                "Mariner-N5_CyaStr",
            ]
        );
    }

    #[test]
    fn organism_and_lineage() {
        let r = parse(SAMPLE, "ig_fam_os.ig");
        assert_eq!(r.organism.as_deref(), Some("Cyathus striatus"));
        assert_eq!(r.oc_lineage.first().map(String::as_str), Some("Cyathus striatus"));
        assert!(r.oc_lineage.iter().any(|t| t == "Fungi"));
        assert_eq!(r.oc_lineage.last().map(String::as_str), Some("Cyathus"));
    }

    #[test]
    fn single_reference_block() {
        let r = parse(SAMPLE, "ig_fam_ref.ig");
        assert_eq!(r.references.len(), 1);
        let rf = &r.references[0];
        assert_eq!(rf.number.as_deref(), Some("[1]  ()"));
        assert_eq!(rf.authors.as_deref(), Some("Bao,W."));
        assert_eq!(rf.title.as_deref(), Some("DNA transposons from the Cyathus striatus genome."));
        assert_eq!(rf.location.as_deref(), Some("Direct Submission to RR (8-Jul-2026)"));
    }

    #[test]
    fn consensus_name_and_sequence() {
        let r = parse(SAMPLE, "ig_fam_seq.ig");
        assert_eq!(r.consensus_name.as_deref(), Some("Mariner-N5_CyaStr"));
        assert_eq!(r.consensus, b"CTGGATAATTTCGACCAAGTGTTCGGACCCTTTCAGGTGTGACG");
        assert_eq!(r.comments, vec!["~96% identical to consensus."]);
    }

    #[test]
    fn multiple_reference_blocks() {
        let body = "\
;ID   X DNA
;RN   [1]
;RA   Author A
;RN   [2]
;RA   Author B
;RT   Second title
;SQ   Sequence 4 BP;
X
ACGT
";
        let r = parse(body, "ig_fam_multiref.ig");
        assert_eq!(r.references.len(), 2);
        assert_eq!(r.references[0].number.as_deref(), Some("[1]"));
        assert_eq!(r.references[0].authors.as_deref(), Some("Author A"));
        assert_eq!(r.references[1].number.as_deref(), Some("[2]"));
        assert_eq!(r.references[1].title.as_deref(), Some("Second title"));
    }

    #[test]
    fn no_feature_table_yields_no_features() {
        // Most Repbase families are nonautonomous: an FH header, no FT lines.
        let body = "\
;ID   X DNA
;FH   Key           Location/Qualifiers
;SQ   Sequence 4 BP;
X
ACGT
";
        assert!(parse(body, "ig_fam_nofeat.ig").features.is_empty());
    }

    #[test]
    fn feature_table_parses_location_and_qualifiers() {
        // Columns matter: key at 5, location/qualifier at 19 (after the ';').
        let body = "\
;ID   X DNA
;FH   Key           Location/Qualifiers
;FT   CDS           490..4569
;FT                 /product=\"X_1p\"
;FT                 /pseudo
;SQ   Sequence 4 BP;
X
ACGT
";
        let r = parse(body, "ig_fam_feat.ig");
        assert_eq!(r.features.len(), 1);
        let f = &r.features[0];
        assert_eq!(f.key, "CDS");
        assert_eq!(f.location, "490..4569");
        assert_eq!(
            f.qualifiers,
            vec![
                ("product".to_string(), "X_1p".to_string()),
                ("pseudo".to_string(), String::new()),
            ]
        );
    }

    #[test]
    fn wrapped_location_and_qualifiers_rejoin_correctly() {
        // A wrapped location concatenates; free text rejoins with a space;
        // /translation rejoins with none (EMBL wraps it mid-token).
        let body = "\
;ID   X DNA
;FH   Key           Location/Qualifiers
;FT   CDS           join(312..602,1346..2053,
;FT                 2317..2817)
;FT                 /note=\"SAP domain, zinc finger,
;FT                 reverse transcriptase.\"
;FT                 /translation=\"MEVTDKVAELVESFTRTGLVKKCEAKNLSTSGTKEEL
;FT                 AARLANLSESEERGAEQ\"
;SQ   Sequence 4 BP;
X
ACGT
";
        let f = &parse(body, "ig_fam_wrap.ig").features[0];
        assert_eq!(f.location, "join(312..602,1346..2053,2317..2817)");
        assert_eq!(f.qualifiers[0].1, "SAP domain, zinc finger, reverse transcriptase.");
        assert_eq!(
            f.qualifiers[1].1,
            "MEVTDKVAELVESFTRTGLVKKCEAKNLSTSGTKEELAARLANLSESEERGAEQ"
        );
    }

    #[test]
    fn multiple_features_are_kept_separate() {
        let body = "\
;ID   X DNA
;FH   Key           Location/Qualifiers
;FT   CDS           265..1518
;FT                 /product=\"orf1\"
;FT   CDS           2314..5043
;FT                 /product=\"orf2\"
;SQ   Sequence 4 BP;
X
ACGT
";
        let r = parse(body, "ig_fam_two.ig");
        assert_eq!(r.features.len(), 2);
        assert_eq!(r.features[0].location, "265..1518");
        assert_eq!(r.features[0].qualifiers[0].1, "orf1");
        assert_eq!(r.features[1].location, "2314..5043");
        assert_eq!(r.features[1].qualifiers[0].1, "orf2");
    }

    #[test]
    fn location_column_is_taken_from_the_fh_header() {
        // If Repbase ever shifts the layout, FH names the new column and the
        // parser follows it rather than mis-slicing at the default of 19.
        let body = "\
;ID   X DNA
;FH   Key   Location/Qualifiers
;FT   CDS   490..4569
;FT         /product=\"X_1p\"
;SQ   Sequence 4 BP;
X
ACGT
";
        let f = &parse(body, "ig_fam_shift.ig").features[0];
        assert_eq!(f.key, "CDS");
        assert_eq!(f.location, "490..4569");
        assert_eq!(f.qualifiers[0].1, "X_1p");
    }
}
