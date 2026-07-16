/// Translate a Repbase IG family record (+ its IG MSA) into a Dfam Stockholm record.
///
/// This is the glue behind `stk repbase-import`.  The family record supplies the
/// metadata (`#=GF` fields); the MSA supplies the consensus (`#=GC RF`) and the
/// aligned instance rows.  Field mapping:
///
/// | Repbase (IG) | Dfam (Stockholm) |
/// |---|---|
/// | `ID` identifier token | `#=GF ID` |
/// | *(boilerplate)* | `#=GF DE` = "Repbase TE family" |
/// | `KW` (up to "Transposable Element") | `#=GF TP` via class lookup (stub) |
/// | `OS` species name | `#=GF OC` (name only, not the `OC` lineage) |
/// | `RN`/`RA`/`RT`/`RL` | Stockholm reference block |
/// | `DE` + `CC` | `#=GF CC` |
/// | MSA consensus (reference row) | `#=GC RF` |
/// | MSA instances | sequence rows (id `name:start-end_orient`) |
use crate::alignment::{MultiAlign, Orientation};
use crate::dfam::lint::rf_consensus_status;
use crate::dfam::record::{RawDfamRecord, SeqRow};
use crate::io::ig_family::{IgFamilyRecord, IgReference};
use dfam_stk_io::is_gap;

/// Result of translating an IG family record + MSA: the Stockholm record plus any
/// non-fatal warnings (e.g. the classification could not be resolved).
pub struct RepbaseImport {
    pub record: RawDfamRecord,
    pub warnings: Vec<String>,
}

/// Build the classification lookup key from the keywords: the `KW` tokens up to
/// **and including** the marker `"Transposable Element"`, joined with `"; "`.
///
/// e.g. `["Mariner/Tc1", "DNA transposon", "Transposable Element", "nonautonomous", …]`
/// → `"Mariner/Tc1; DNA transposon; Transposable Element"`.  The tokens after the
/// marker (autonomy, family name) vary per family and are dropped.  When the marker
/// is absent, all keywords are used.
pub fn class_lookup_key(keywords: &[String]) -> String {
    let mut key: Vec<String> = Vec::new();
    for k in keywords {
        key.push(k.clone());
        if k.eq_ignore_ascii_case("Transposable Element") {
            break;
        }
    }
    key.join("; ")
}

/// Repbase keyword classification → Dfam `#=GF TP` lineage.
///
/// Seed table; extend as more Repbase superfamilies are encountered.  Keys are the
/// [`class_lookup_key`] form (matched case-insensitively).
const KW_CLASS_TABLE: &[(&str, &str)] = &[
    ("Gypsy; LTR Retrotransposon; Transposable Element",
     "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Gypsy"),
    ("Copia; LTR Retrotransposon; Transposable Element",
     "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Ty1-Copia"),
    ("BEL; LTR Retrotransposon; Transposable Element",
     "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Bel-Pao"),
    ("hAT; DNA transposon; Transposable Element",
     "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT"),
    ("Mariner/Tc1; DNA transposon; Transposable Element",
     "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner"),
    ("L1; Non-LTR Retrotransposon; Transposable Element",
     "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;L1-like;L1-group;L1"),
    ("ERV1; Endogenous Retrovirus; Transposable Element",
     "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV1"),
    ("MuDR; DNA transposon; Transposable Element",
     "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Mutator-like;MuDR"),
    ("ERV2; Endogenous Retrovirus; Transposable Element",
     "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2-group;ERV2"),
    ("Harbinger; DNA transposon; Transposable Element",
     "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PIF-Harbinger;Harbinger"),
    ("SINE2/tRNA; SINE; Non-LTR Retrotransposon; Transposable Element",
     "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;Unknown_LINE-dependent"),
    ("DNA transposon; Transposable Element",
     "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase"),
];

/// Map Repbase keywords to a Dfam classification string (`#=GF TP` value), or
/// `None` when the keyword combination is not in the seed table.
pub fn lookup_class(keywords: &[String]) -> Option<String> {
    let key = class_lookup_key(keywords);
    KW_CLASS_TABLE
        .iter()
        .find(|(k, _)| k.eq_ignore_ascii_case(&key))
        .map(|(_, v)| v.to_string())
}

/// Translate a family record and its MSA into a single Dfam Stockholm record.
pub fn to_stk_record(family: &IgFamilyRecord, msa: &MultiAlign) -> RepbaseImport {
    let mut warnings = Vec::new();
    let mut rec = RawDfamRecord::default();
    rec.record_num = 1;
    rec.terminated = true;

    // Partition references.  A "Direct Submission to RR / Repbase Reports" is not a
    // real citation: its author becomes the Dfam curator (AU) and its date a curator
    // note (**).  Every other reference is kept as a normal Stockholm ref block.
    let mut au_authors: Vec<String> = Vec::new();
    let mut submission_notes: Vec<String> = Vec::new();
    let mut normal_refs: Vec<&IgReference> = Vec::new();
    for r in &family.references {
        let is_rr = r.location.as_deref().is_some_and(is_repbase_reports_submission);
        if is_rr {
            if let Some(a) = &r.authors {
                au_authors.push(a.clone());
            }
            // Prefer the date in the RL's parentheses; fall back to the record's DT.
            let date = r
                .location
                .as_deref()
                .and_then(extract_submission_date)
                .or_else(|| family.date.clone());
            if let Some(d) = date {
                submission_notes.push(format!("Repbase Reports Submission Date: {}", d));
            }
        } else {
            normal_refs.push(r);
        }
    }

    // ── #=GF fields, in Dfam order ──────────────────────────────────────────
    rec.gf.push(("ID".into(), family.id.clone()));
    rec.gf.push(("DE".into(), "Repbase TE family".into()));

    if !au_authors.is_empty() {
        rec.gf.push(("AU".into(), au_authors.join("; ")));
    }

    match lookup_class(&family.keywords) {
        Some(tp) => rec.gf.push(("TP".into(), tp)),
        None => warnings.push(format!(
            "could not map Repbase keywords to a Dfam classification \
             (lookup key: {:?}); #=GF TP left unset",
            class_lookup_key(&family.keywords)
        )),
    }

    // OS → OC, taxon name only (the ;OC lineage is intentionally dropped).
    if let Some(os) = &family.organism {
        rec.gf.push(("OC".into(), os.clone()));
    }

    rec.gf.push(("SQ".into(), msa.num_instances().to_string()));

    // Real reference blocks, renumbered sequentially: RN then RT/RA/RL.
    for (i, r) in normal_refs.iter().enumerate() {
        rec.gf.push(("RN".into(), format!("[{}]", i + 1)));
        if let Some(t) = &r.title {
            rec.gf.push(("RT".into(), t.clone()));
        }
        if let Some(a) = &r.authors {
            rec.gf.push(("RA".into(), a.clone()));
        }
        if let Some(l) = &r.location {
            rec.gf.push(("RL".into(), l.clone()));
        }
    }

    // Repbase DE → Dfam CC, followed by any Repbase CC lines.
    if let Some(de) = &family.description {
        rec.gf.push(("CC".into(), de.clone()));
    }
    for cc in &family.comments {
        rec.gf.push(("CC".into(), cc.clone()));
    }

    // Curator notes (**) carrying the Repbase Reports submission date(s).
    for note in submission_notes {
        rec.gf.push(("**".into(), note));
    }

    // ── #=GC RF (consensus) and sequence rows (instances) ───────────────────
    if let Some(reference) = msa.sequences.first() {
        rec.gc.insert("RF".into(), seq_to_string(&reference.seq));
    }
    for row in msa.sequences.iter().skip(1) {
        rec.sequences
            .push(SeqRow::new_raw(instance_id(row), seq_to_string(&row.seq)));
    }

    // ── Consensus validation (advisory warnings) ────────────────────────────
    // These do not change the output; they surface consensus discrepancies at
    // import time so real divergent examples can be collected before deciding how
    // to handle them.
    warnings.extend(validate_consensus(family, msa, &rec));

    RepbaseImport { record: rec, warnings }
}

/// Advisory consensus checks run during import.
///
/// (A) The `--msa` consensus (reference row) and the `--record` file's own
/// consensus should describe the same ungapped sequence — a check only the import
/// can make, since `stk lint` never sees the two source files.
///
/// (B) The consensus (`#=GC RF`) should reproduce the consensus called from the
/// aligned instances — the same test `stk lint` reports as `rf_consensus_mismatch`,
/// run here (via [`rf_consensus_status`]) so the import alerts without a separate
/// lint pass.
fn validate_consensus(
    family: &IgFamilyRecord,
    msa: &MultiAlign,
    rec: &RawDfamRecord,
) -> Vec<String> {
    let mut w = Vec::new();

    // (A) MSA consensus vs record consensus.
    if let Some(reference) = msa.sequences.first() {
        let msa_cons = ungapped_upper(&reference.seq);
        let rec_cons = ungapped_upper(&family.consensus);
        if !msa_cons.is_empty() && !rec_cons.is_empty() && msa_cons != rec_cons {
            w.push(format!(
                "consensus mismatch: the --msa consensus ({} bp) and the --record \
                 consensus ({} bp) are not the same sequence",
                msa_cons.len(),
                rec_cons.len()
            ));
        }
    }

    // (B) RF vs the consensus called from the aligned instances (mirrors lint).
    if rf_consensus_status(rec) == Some(false) {
        w.push(
            "consensus mismatch: the MSA consensus does not match the consensus \
             called from the aligned instances (stk lint reports this as \
             rf_consensus_mismatch)"
                .to_string(),
        );
    }

    w
}

/// Strip gap/space characters and upper-case, for consensus comparison.
fn ungapped_upper(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .filter(|&&b| !is_gap(b) && b != b' ')
        .map(|b| b.to_ascii_uppercase())
        .collect()
}

/// `true` if a reference location denotes a Repbase Reports direct submission
/// (e.g. `"Direct Submission to RR (8-Jul-2026)"`), which Dfam does not treat as a
/// citation.  Matches the words "Repbase Reports" or a standalone `RR` token.
fn is_repbase_reports_submission(location: &str) -> bool {
    location.to_ascii_lowercase().contains("repbase reports")
        || location
            .split(|c: char| !c.is_ascii_alphanumeric())
            .any(|tok| tok == "RR")
}

/// Extract the date from the last parenthetical group of a reference location,
/// e.g. `"Direct Submission to RR (8-Jul-2026)"` → `"8-Jul-2026"`.
fn extract_submission_date(location: &str) -> Option<String> {
    let open = location.rfind('(')?;
    let close = location[open..].find(')')? + open;
    let inner = location[open + 1..close].trim();
    (!inner.is_empty()).then(|| inner.to_string())
}

/// Convert a gapped alignment row to a Stockholm-safe string: flanking space
/// padding becomes `.` (a literal space would break column parsing).  Interior
/// `-` gaps are left as-is for clean-on-write to normalize.
fn seq_to_string(seq: &[u8]) -> String {
    seq.iter()
        .map(|&b| if b == b' ' { '.' } else { b as char })
        .collect()
}

/// Format an instance identifier as `name:start-end_orient`, or bare `name` when
/// the row carries no coordinates.
fn instance_id(row: &crate::alignment::SequenceRow) -> String {
    if row.seq_start == 0 && row.seq_end == 0 {
        return row.name.clone();
    }
    let orient = match row.orient {
        Orientation::Forward => '+',
        Orientation::Reverse => '-',
    };
    format!("{}:{}-{}_{}", row.name, row.seq_start, row.seq_end, orient)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::{MultiAlign, SequenceRow};
    use crate::io::ig_family::{IgFamilyRecord, IgReference};

    fn sample_family() -> IgFamilyRecord {
        IgFamilyRecord {
            id: "Mariner-N5_CyaStr".into(),
            id_line: "Mariner-N5_CyaStr DNA".into(),
            description: Some("DNA transposon from Cyathus striatus, consensus.".into()),
            accession: Some(".".into()),
            date: Some("18-MAR-2026 (Created)".into()),
            keywords: vec![
                "Mariner/Tc1".into(),
                "DNA transposon".into(),
                "Transposable Element".into(),
                "nonautonomous".into(),
            ],
            organism: Some("Cyathus striatus".into()),
            oc_lineage: vec!["Eukaryota".into(), "Fungi".into()],
            references: vec![IgReference {
                number: Some("[1]  ()".into()),
                authors: Some("Bao,W.".into()),
                title: Some("DNA transposons from Cyathus striatus.".into()),
                location: Some("Direct Submission to RR (8-Jul-2026)".into()),
            }],
            comments: vec!["~96% identical to consensus.".into()],
            sq_summary: Some("Sequence 8 BP;".into()),
            consensus_name: Some("Mariner-N5_CyaStr".into()),
            // Ungapped form of the sample MSA consensus row (ACGT-ACGT → ACGTACGT).
            consensus: b"ACGTACGT".to_vec(),
        }
    }

    fn sample_msa() -> MultiAlign {
        let cons = SequenceRow::new("Mariner-N5_CyaStr", b"ACGT-ACGT".to_vec());
        let mut inst = SequenceRow::new("JAOVFP01_1", b"ACGT-ACGT".to_vec());
        inst.seq_start = 1;
        inst.seq_end = 9;
        MultiAlign::from_sequences(cons, vec![inst])
    }

    fn gf<'a>(r: &'a RawDfamRecord, tag: &str) -> Vec<&'a str> {
        r.gf.iter().filter(|(t, _)| t == tag).map(|(_, v)| v.as_str()).collect()
    }

    #[test]
    fn maps_core_metadata_fields() {
        let out = to_stk_record(&sample_family(), &sample_msa());
        let r = &out.record;
        assert_eq!(gf(r, "ID"), ["Mariner-N5_CyaStr"]);
        assert_eq!(gf(r, "DE"), ["Repbase TE family"]);
        assert_eq!(gf(r, "OC"), ["Cyathus striatus"]); // OS name, not lineage
        assert_eq!(gf(r, "SQ"), ["1"]);
    }

    #[test]
    fn description_and_comment_become_cc() {
        let out = to_stk_record(&sample_family(), &sample_msa());
        let cc = gf(&out.record, "CC");
        assert_eq!(
            cc,
            [
                "DNA transposon from Cyathus striatus, consensus.",
                "~96% identical to consensus.",
            ]
        );
    }

    #[test]
    fn repbase_reports_ref_becomes_au_and_note() {
        // The sample's only reference is a Direct Submission to RR: author → AU,
        // date → ** note, and it is NOT emitted as a Stockholm reference block.
        let out = to_stk_record(&sample_family(), &sample_msa());
        let r = &out.record;
        assert_eq!(gf(r, "AU"), ["Bao,W."]);
        assert!(gf(r, "RN").is_empty(), "RR submission must not be a ref block");
        assert!(gf(r, "RT").is_empty());
        assert_eq!(gf(r, "**"), ["Repbase Reports Submission Date: 8-Jul-2026"]);
    }

    #[test]
    fn real_reference_kept_as_stockholm_block() {
        let mut fam = sample_family();
        fam.references = vec![
            IgReference {
                number: Some("[1]".into()),
                authors: Some("Bao,W.".into()),
                title: None,
                location: Some("Direct Submission to RR (8-Jul-2026)".into()),
            },
            IgReference {
                number: Some("[2]".into()),
                authors: Some("Smith,J.; Doe,A.".into()),
                title: Some("A real paper".into()),
                location: Some("J Mol Biol 2020".into()),
            },
        ];
        let out = to_stk_record(&fam, &sample_msa());
        let r = &out.record;
        // RR ref → AU + ** ; real ref → ref block, renumbered to [1].
        assert_eq!(gf(r, "AU"), ["Bao,W."]);
        assert_eq!(gf(r, "RN"), ["[1]"]);
        assert_eq!(gf(r, "RA"), ["Smith,J.; Doe,A."]);
        assert_eq!(gf(r, "RT"), ["A real paper"]);
        assert_eq!(gf(r, "RL"), ["J Mol Biol 2020"]);
        assert_eq!(gf(r, "**"), ["Repbase Reports Submission Date: 8-Jul-2026"]);
    }

    #[test]
    fn known_class_maps_to_tp() {
        // sample_family()'s Mariner keywords are in the seed table.
        let out = to_stk_record(&sample_family(), &sample_msa());
        assert_eq!(
            gf(&out.record, "TP"),
            ["Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner"]
        );
        assert!(out.warnings.is_empty(), "unexpected warnings: {:?}", out.warnings);
    }

    #[test]
    fn unmapped_class_warns_and_leaves_tp_unset() {
        let mut fam = sample_family();
        fam.keywords = vec![
            "Bogus/Unknown".into(),
            "Mystery transposon".into(),
            "Transposable Element".into(),
        ];
        let out = to_stk_record(&fam, &sample_msa());
        assert!(gf(&out.record, "TP").is_empty());
        assert_eq!(out.warnings.len(), 1);
        assert!(
            out.warnings[0].contains("Bogus/Unknown; Mystery transposon; Transposable Element"),
            "{}",
            out.warnings[0]
        );
    }

    #[test]
    fn seed_table_lookups() {
        let kw = |s: &[&str]| s.iter().map(|x| x.to_string()).collect::<Vec<_>>();
        assert!(lookup_class(&kw(&["Gypsy", "LTR Retrotransposon", "Transposable Element"]))
            .unwrap()
            .ends_with(";Gypsy-ERV;Gypsy"));
        assert_eq!(
            lookup_class(&kw(&["DNA transposon", "Transposable Element"])).as_deref(),
            Some("Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase")
        );
        assert!(lookup_class(&kw(&[
            "SINE2/tRNA",
            "SINE",
            "Non-LTR Retrotransposon",
            "Transposable Element"
        ]))
        .unwrap()
        .contains("Unknown_LINE-dependent"));
    }

    #[test]
    fn consensus_becomes_rf_and_instances_become_rows() {
        let out = to_stk_record(&sample_family(), &sample_msa());
        let r = &out.record;
        assert_eq!(r.gc.get("RF").map(String::as_str), Some("ACGT-ACGT"));
        assert_eq!(r.sequences.len(), 1);
        assert_eq!(r.sequences[0].original_id, "JAOVFP01_1:1-9_+");
        assert_eq!(r.sequences[0].aligned_seq, "ACGT-ACGT");
    }

    #[test]
    fn matching_consensus_produces_no_consensus_warning() {
        let out = to_stk_record(&sample_family(), &sample_msa());
        assert!(
            !out.warnings.iter().any(|w| w.contains("consensus mismatch")),
            "unexpected consensus warning: {:?}",
            out.warnings
        );
    }

    #[test]
    fn msa_vs_record_consensus_mismatch_warns() {
        let mut fam = sample_family();
        fam.consensus = b"TTTTTTTT".to_vec(); // differs from the MSA consensus (ACGTACGT)
        let out = to_stk_record(&fam, &sample_msa());
        assert!(
            out.warnings
                .iter()
                .any(|w| w.contains("--msa consensus") && w.contains("--record consensus")),
            "{:?}",
            out.warnings
        );
    }

    #[test]
    fn rf_vs_called_consensus_mismatch_warns() {
        // RF (reference row) says AAAA, but every instance calls CCCC.
        let cons = SequenceRow::new("cons", b"AAAA".to_vec());
        let i1 = SequenceRow::new("i1", b"CCCC".to_vec());
        let i2 = SequenceRow::new("i2", b"CCCC".to_vec());
        let msa = MultiAlign::from_sequences(cons, vec![i1, i2]);
        let mut fam = sample_family();
        fam.consensus = b"AAAA".to_vec(); // match the MSA consensus so only check (B) fires
        let out = to_stk_record(&fam, &msa);
        assert!(
            out.warnings.iter().any(|w| w.contains("rf_consensus_mismatch")),
            "{:?}",
            out.warnings
        );
    }

    #[test]
    fn class_lookup_key_includes_transposable_element() {
        let kws = vec![
            "Mariner/Tc1".to_string(),
            "DNA transposon".to_string(),
            "Transposable Element".to_string(),
            "nonautonomous".to_string(),
        ];
        assert_eq!(
            class_lookup_key(&kws),
            "Mariner/Tc1; DNA transposon; Transposable Element"
        );
    }
}
