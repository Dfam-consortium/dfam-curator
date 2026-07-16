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
    /// `SQ` summary line content, if present.
    pub sq_summary: Option<String>,
    /// The bare identifier line that introduces the consensus (after `SQ`).
    pub consensus_name: Option<String>,
    /// Concatenated consensus sequence (ungapped).
    pub consensus: Vec<u8>,
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
                "XX" | "FH" => {} // separators / fixed header
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
}
