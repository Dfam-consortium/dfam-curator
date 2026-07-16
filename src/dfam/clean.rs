/// Normalization applied when writing Dfam Stockholm output.
///
/// Every `stk` subcommand that emits a Stockholm file cleans each record on the
/// way out — unless `--no-clean` is given — so the result conforms to Dfam
/// conventions.  This is the write-side counterpart to the lint checks that only
/// *report* the same issues (`seq_nonstandard_gap`, `ac_format`): lint warns,
/// writing cleans.
///
/// Current transforms:
/// - Gap characters `-`, `_`, `~` in sequence rows and the `#=GC RF` line → `.`.
/// - A 7-digit `#=GF AC` accession widened to the 9-digit standard.
use crate::dfam::record::RawDfamRecord;
use dfam_stk_io::{is_gap, DFAM_GAP};

/// Tally of the changes [`clean_record`] made, for reporting back to the user.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct CleanReport {
    /// Sequence rows whose gap characters were rewritten to `.`.
    pub gap_rows: usize,
    /// `#=GC RF` lines whose gap characters were rewritten (0 or 1 per record).
    pub gap_rf: usize,
    /// `#=GF AC` fields widened from 7 to 9 digits.
    pub ac_padded: usize,
}

impl CleanReport {
    /// `true` if any change was made.
    pub fn changed(&self) -> bool {
        self.gap_rows > 0 || self.gap_rf > 0 || self.ac_padded > 0
    }

    /// Accumulate another report's counts into this one.
    pub fn merge(&mut self, other: CleanReport) {
        self.gap_rows += other.gap_rows;
        self.gap_rf += other.gap_rf;
        self.ac_padded += other.ac_padded;
    }

    /// A one-line human summary of what changed, or `None` when nothing did.
    pub fn summary(&self) -> Option<String> {
        if !self.changed() {
            return None;
        }
        let mut parts = Vec::new();
        if self.gap_rows > 0 {
            parts.push(format!(
                "{} sequence row{} regapped to '.'",
                self.gap_rows,
                if self.gap_rows == 1 { "" } else { "s" }
            ));
        }
        if self.gap_rf > 0 {
            parts.push(format!(
                "{} RF line{} regapped to '.'",
                self.gap_rf,
                if self.gap_rf == 1 { "" } else { "s" }
            ));
        }
        if self.ac_padded > 0 {
            parts.push(format!(
                "{} AC widened to 9 digits",
                self.ac_padded
            ));
        }
        Some(parts.join(", "))
    }
}

/// Rewrite every non-`.` Stockholm gap character (`-`, `_`, `~`) in `s` to `.`.
///
/// Returns `true` if anything changed.  Only ASCII gap bytes are swapped for the
/// ASCII `.`, so UTF-8 validity is preserved for any (already-ASCII) sequence data.
fn normalize_gaps(s: &mut String) -> bool {
    let mut changed = false;
    let out: Vec<u8> = s
        .bytes()
        .map(|b| {
            if is_gap(b) && b != DFAM_GAP {
                changed = true;
                DFAM_GAP
            } else {
                b
            }
        })
        .collect();
    if changed {
        *s = String::from_utf8(out).expect("gap normalization only swaps ASCII bytes");
    }
    changed
}

/// Widen a 7-digit `DF`/`DR` accession to the 9-digit standard, preserving the
/// prefix and any `.N` version suffix.
///
/// Returns the padded accession only when `ac` is a widenable 7-digit accession.
/// Already-9-digit, malformed, or empty accessions return `None` and are left
/// untouched (they are lint's concern, not the cleaner's).
fn pad_ac(ac: &str) -> Option<String> {
    let ac = ac.trim();

    // Split off an optional ".N" version suffix (the '.' is kept with `version`).
    let (base, version) = match ac.rfind('.') {
        Some(pos)
            if pos + 1 < ac.len() && ac[pos + 1..].chars().all(|c| c.is_ascii_digit()) =>
        {
            (&ac[..pos], &ac[pos..])
        }
        _ => (ac, ""),
    };

    if base.len() < 2 {
        return None;
    }
    let (prefix, digits) = base.split_at(2);
    if !matches!(prefix, "DF" | "DR") {
        return None;
    }
    if digits.len() != 7 || !digits.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }

    // 7 → 9 digits: prepend two zeros to the numeric block (value-preserving).
    Some(format!("{}00{}{}", prefix, digits, version))
}

/// Normalize `record` in place to Dfam-standard output conventions.
///
/// See the module docs for the transforms applied.  Returns a [`CleanReport`]
/// tallying what changed so callers can summarise it for the user.
pub fn clean_record(record: &mut RawDfamRecord) -> CleanReport {
    let mut report = CleanReport::default();

    for row in record.sequences.iter_mut() {
        if normalize_gaps(&mut row.aligned_seq) {
            report.gap_rows += 1;
        }
    }

    if let Some(rf) = record.gc.get_mut("RF") {
        if normalize_gaps(rf) {
            report.gap_rf += 1;
        }
    }

    for (tag, value) in record.gf.iter_mut() {
        if tag == "AC" {
            if let Some(padded) = pad_ac(value) {
                *value = padded;
                report.ac_padded += 1;
            }
        }
    }

    report
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dfam::record::{RawDfamRecord, SeqRow};

    #[test]
    fn normalize_gaps_rewrites_all_gap_flavors() {
        let mut s = "AC-GT_A~T.G".to_string();
        assert!(normalize_gaps(&mut s));
        assert_eq!(s, "AC.GT.A.T.G");
    }

    #[test]
    fn normalize_gaps_noop_when_already_dotted() {
        let mut s = "AC.GT.".to_string();
        assert!(!normalize_gaps(&mut s));
        assert_eq!(s, "AC.GT.");
    }

    #[test]
    fn pad_ac_widens_seven_to_nine() {
        assert_eq!(pad_ac("DF0001234").as_deref(), Some("DF000001234"));
        assert_eq!(pad_ac("DR0000001").as_deref(), Some("DR000000001"));
    }

    #[test]
    fn pad_ac_preserves_version_suffix() {
        assert_eq!(pad_ac("DF0001234.2").as_deref(), Some("DF000001234.2"));
    }

    #[test]
    fn pad_ac_leaves_nine_digit_untouched() {
        assert_eq!(pad_ac("DF000001234"), None);
        assert_eq!(pad_ac("DR000000001.3"), None);
    }

    #[test]
    fn pad_ac_ignores_malformed() {
        assert_eq!(pad_ac("BADFORMAT"), None);
        assert_eq!(pad_ac("DF12345"), None); // 5 digits
        assert_eq!(pad_ac("XY0001234"), None); // wrong prefix
        assert_eq!(pad_ac(""), None);
    }

    #[test]
    fn clean_record_normalizes_sequences_rf_and_ac() {
        let mut r = RawDfamRecord::default();
        r.gf.push(("AC".to_string(), "DF0001234".to_string()));
        r.gc.insert("RF".to_string(), "AC-GT".to_string());
        r.sequences.push(SeqRow::new_raw("s1", "AC-GT"));
        r.sequences.push(SeqRow::new_raw("s2", "AC_GT"));

        let report = clean_record(&mut r);

        assert_eq!(report.gap_rows, 2);
        assert_eq!(report.gap_rf, 1);
        assert_eq!(report.ac_padded, 1);
        assert!(report.changed());
        assert_eq!(r.gf_first("AC"), Some("DF000001234"));
        assert_eq!(r.gc.get("RF").map(String::as_str), Some("AC.GT"));
        assert_eq!(r.sequences[0].aligned_seq, "AC.GT");
        assert_eq!(r.sequences[1].aligned_seq, "AC.GT");
    }

    #[test]
    fn clean_record_clean_input_reports_no_change() {
        let mut r = RawDfamRecord::default();
        r.gf.push(("AC".to_string(), "DF000001234".to_string()));
        r.gc.insert("RF".to_string(), "AC.GT".to_string());
        r.sequences.push(SeqRow::new_raw("s1", "AC.GT"));

        let report = clean_record(&mut r);
        assert!(!report.changed());
        assert_eq!(report.summary(), None);
    }
}
