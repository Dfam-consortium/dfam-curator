/// Edit operations for `#=GF` annotations in a `RawDfamRecord`.
use crate::dfam::record::RawDfamRecord;
use regex::Regex;

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

    /// Apply a regex substitution to every existing value for `tag`.
    /// If `all` is true, replaces every match within each value; otherwise
    /// replaces only the first match.
    Sub { tag: String, pattern: Regex, replacement: String, all: bool },
}

/// Apply a slice of operations to `record` in order.
///
/// Operations are applied in the following fixed sequence regardless of
/// command-line position: Delete → Set → Append → Sub.
///
/// Sub runs last so it transforms all values present after Set and Append
/// have completed, regardless of their origin.  It applies once per value
/// string (or globally within each string when `all` is true).
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

    for op in ops.iter().filter(|o| matches!(o, Op::Sub { .. })) {
        let Op::Sub { tag, pattern, replacement, all } = op else { unreachable!() };
        for (t, v) in record.gf.iter_mut() {
            if t == tag {
                *v = if *all {
                    pattern.replace_all(v, replacement.as_str()).into_owned()
                } else {
                    pattern.replace(v, replacement.as_str()).into_owned()
                };
            }
        }
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

    #[test]
    fn sub_first_match_only() {
        let mut r = make(&[("ID", "MyFam-")]);
        apply_ops(&mut r, &[Op::Sub {
            tag: "ID".into(),
            pattern: Regex::new(r"^(.*)-$").unwrap(),
            replacement: "$1".into(),
            all: false,
        }]);
        assert_eq!(gf(&r, "ID"), ["MyFam"]);
    }

    #[test]
    fn sub_global_replaces_all_matches() {
        let mut r = make(&[("DE", "foo bar foo")]);
        apply_ops(&mut r, &[Op::Sub {
            tag: "DE".into(),
            pattern: Regex::new("foo").unwrap(),
            replacement: "baz".into(),
            all: true,
        }]);
        assert_eq!(gf(&r, "DE"), ["baz bar baz"]);
    }

    #[test]
    fn sub_without_global_replaces_first_only() {
        let mut r = make(&[("DE", "foo bar foo")]);
        apply_ops(&mut r, &[Op::Sub {
            tag: "DE".into(),
            pattern: Regex::new("foo").unwrap(),
            replacement: "baz".into(),
            all: false,
        }]);
        assert_eq!(gf(&r, "DE"), ["baz bar foo"]);
    }

    #[test]
    fn sub_applies_to_all_tag_lines() {
        let mut r = make(&[("OC", "Mus-musculus"), ("OC", "Homo-sapiens")]);
        apply_ops(&mut r, &[Op::Sub {
            tag: "OC".into(),
            pattern: Regex::new("-").unwrap(),
            replacement: " ".into(),
            all: false,
        }]);
        assert_eq!(gf(&r, "OC"), ["Mus musculus", "Homo sapiens"]);
    }

    #[test]
    fn sub_after_set_transforms_new_value() {
        let mut r = make(&[("ID", "old")]);
        apply_ops(&mut r, &[
            Op::Set { tag: "ID".into(), value: "new-value-".into() },
            Op::Sub {
                tag: "ID".into(),
                pattern: Regex::new(r"^(.*)-$").unwrap(),
                replacement: "$1".into(),
                all: false,
            },
        ]);
        assert_eq!(gf(&r, "ID"), ["new-value"]);
    }

    #[test]
    fn sub_after_append_transforms_appended_values() {
        let mut r = make(&[("OC", "existing")]);
        apply_ops(&mut r, &[
            Op::Append { tag: "OC".into(), value: "appended-".into() },
            Op::Sub {
                tag: "OC".into(),
                pattern: Regex::new(r"-$").unwrap(),
                replacement: "".into(),
                all: false,
            },
        ]);
        assert_eq!(gf(&r, "OC"), ["existing", "appended"]);
    }

    #[test]
    fn sub_negated_character_class() {
        let mut r = make(&[("DE", "abc123def")]);
        apply_ops(&mut r, &[Op::Sub {
            tag: "DE".into(),
            pattern: Regex::new(r"[^a-z]+").unwrap(),
            replacement: "_".into(),
            all: true,
        }]);
        assert_eq!(gf(&r, "DE"), ["abc_def"]);
    }

    #[test]
    fn sub_no_match_is_noop() {
        let mut r = make(&[("ID", "MyFam")]);
        apply_ops(&mut r, &[Op::Sub {
            tag: "ID".into(),
            pattern: Regex::new("NOMATCH").unwrap(),
            replacement: "x".into(),
            all: false,
        }]);
        assert_eq!(gf(&r, "ID"), ["MyFam"]);
    }

    #[test]
    fn sub_absent_tag_is_noop() {
        let mut r = make(&[("DE", "desc")]);
        apply_ops(&mut r, &[Op::Sub {
            tag: "ID".into(),
            pattern: Regex::new("x").unwrap(),
            replacement: "y".into(),
            all: false,
        }]);
        assert!(gf(&r, "ID").is_empty());
        assert_eq!(gf(&r, "DE"), ["desc"]);
    }
}
