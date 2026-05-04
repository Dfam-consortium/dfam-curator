/// Edit operations for `#=GF` annotations in a `RawDfamRecord`.
use crate::dfam::record::RawDfamRecord;

/// A single edit operation applied to a record's GF fields.
#[derive(Debug, Clone)]
pub enum Op {
    /// Replace all occurrences of `tag` with `value`.
    /// If the tag is absent the entry is appended at the end.
    Set { tag: String, value: String },

    /// Remove every occurrence of `tag`.
    Delete { tag: String },

    /// Append a new `(tag, value)` pair.
    /// Intended for legitimately multi-valued fields such as OC and CC.
    Append { tag: String, value: String },
}

/// Apply a slice of operations to `record` in order.
///
/// Operations are applied in the following fixed sequence regardless of
/// command-line position: all `Delete` ops first, then all `Set` ops, then
/// all `Append` ops.  This ensures a `--delete` + `--set` pair for the same
/// tag behaves predictably.
pub fn apply_ops(record: &mut RawDfamRecord, ops: &[Op]) {
    for op in ops.iter().filter(|o| matches!(o, Op::Delete { .. })) {
        let Op::Delete { tag } = op else { unreachable!() };
        record.gf.retain(|(t, _)| t != tag);
    }

    for op in ops.iter().filter(|o| matches!(o, Op::Set { .. })) {
        let Op::Set { tag, value } = op else { unreachable!() };
        let mut placed = false;
        record.gf.retain_mut(|(t, v)| {
            if t == tag {
                if !placed {
                    *v = value.clone();
                    placed = true;
                    true
                } else {
                    false // drop extra occurrences
                }
            } else {
                true
            }
        });
        if !placed {
            record.gf.push((tag.clone(), value.clone()));
        }
    }

    for op in ops.iter().filter(|o| matches!(o, Op::Append { .. })) {
        let Op::Append { tag, value } = op else { unreachable!() };
        record.gf.push((tag.clone(), value.clone()));
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dfam::record::RawDfamRecord;

    fn make(gf: &[(&str, &str)]) -> RawDfamRecord {
        let mut r = RawDfamRecord::default();
        for (t, v) in gf {
            r.gf.push((t.to_string(), v.to_string()));
        }
        r
    }

    fn gf(r: &RawDfamRecord, tag: &str) -> Vec<String> {
        r.gf_all(tag).iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn set_replaces_existing() {
        let mut r = make(&[("AU", "Old Author"), ("DE", "desc")]);
        apply_ops(&mut r, &[Op::Set { tag: "AU".into(), value: "New Author".into() }]);
        assert_eq!(gf(&r, "AU"), ["New Author"]);
        assert_eq!(gf(&r, "DE"), ["desc"]); // unchanged
    }

    #[test]
    fn set_adds_when_absent() {
        let mut r = make(&[("DE", "desc")]);
        apply_ops(&mut r, &[Op::Set { tag: "AU".into(), value: "Smith J".into() }]);
        assert_eq!(gf(&r, "AU"), ["Smith J"]);
    }

    #[test]
    fn set_collapses_duplicate_tags() {
        let mut r = make(&[("OC", "root"), ("OC", "Mammalia"), ("OC", "Mus musculus")]);
        apply_ops(&mut r, &[Op::Set { tag: "OC".into(), value: "Homo sapiens".into() }]);
        assert_eq!(gf(&r, "OC"), ["Homo sapiens"]);
    }

    #[test]
    fn delete_removes_tag() {
        let mut r = make(&[("AU", "Smith"), ("SE", "old source"), ("DE", "desc")]);
        apply_ops(&mut r, &[Op::Delete { tag: "SE".into() }]);
        assert!(gf(&r, "SE").is_empty());
        assert_eq!(gf(&r, "AU"), ["Smith"]);
    }

    #[test]
    fn append_adds_new_entry() {
        let mut r = make(&[("OC", "Mammalia")]);
        apply_ops(&mut r, &[Op::Append { tag: "OC".into(), value: "Mus musculus".into() }]);
        assert_eq!(gf(&r, "OC"), ["Mammalia", "Mus musculus"]);
    }

    #[test]
    fn delete_then_set_same_tag() {
        let mut r = make(&[("AU", "Wrong Author")]);
        apply_ops(&mut r, &[
            Op::Delete { tag: "AU".into() },
            Op::Set    { tag: "AU".into(), value: "Right Author".into() },
        ]);
        assert_eq!(gf(&r, "AU"), ["Right Author"]);
    }
}
