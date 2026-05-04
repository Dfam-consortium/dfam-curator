/// Cache management for tier-2 validation data.
///
/// Four optional data files live under a single cache directory:
///   - `taxonomy.tsv`        — one NCBI scientific taxon name per line
///   - `taxonomy-common.tsv` — common/alternate name TAB scientific name (for suggestions)
///   - `classification.txt`  — one valid Dfam TP string per line
///   - `dfam-names.txt`      — one Dfam family ID per line
///
/// The cache directory defaults to `$STK_CACHE_DIR` → `~/.cache/stk`.
/// Each file is optional; absent files cause the corresponding tier-2 check
/// to be skipped (with an informational message from the caller).
use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::path::{Path, PathBuf};

/// Loaded tier-2 validation data.  Each field is `None` when the
/// corresponding cache file is absent.
pub struct Cache {
    pub classification:  Option<HashSet<String>>,
    pub taxonomy:        Option<HashSet<String>>,
    /// Alternate name → scientific name (keys are lowercased for
    /// case-insensitive lookup).  Populated from `taxonomy-common.tsv`.
    pub taxonomy_common: Option<HashMap<String, String>>,
    pub dfam_names:      Option<HashSet<String>>,
}

/// Resolve the cache directory: `--cache-dir` arg → `$STK_CACHE_DIR` → `~/.cache/stk`.
pub fn cache_dir() -> PathBuf {
    if let Ok(d) = std::env::var("STK_CACHE_DIR") {
        return PathBuf::from(d);
    }
    dirs::cache_dir()
        .unwrap_or_else(|| PathBuf::from("."))
        .join("stk")
}

/// Load all available cache files from `dir`.
pub fn load_cache(dir: &Path) -> Cache {
    Cache {
        classification:  load_line_set(dir.join("classification.txt")),
        taxonomy:        load_line_set(dir.join("taxonomy.tsv")),
        taxonomy_common: load_common_map(dir.join("taxonomy-common.tsv")),
        dfam_names:      load_line_set(dir.join("dfam-names.txt")),
    }
}

/// Report which tier-2 checks are unavailable due to missing cache files.
/// Returns a list of `(filename, description)` pairs.
pub fn missing_cache_files(cache: &Cache, dir: &Path) -> Vec<(PathBuf, &'static str)> {
    let mut missing = Vec::new();
    if cache.classification.is_none() {
        missing.push((dir.join("classification.txt"), "TP classification validation"));
    }
    if cache.taxonomy.is_none() {
        missing.push((dir.join("taxonomy.tsv"), "OC NCBI taxonomy validation"));
    }
    if cache.taxonomy_common.is_none() {
        missing.push((dir.join("taxonomy-common.tsv"), "OC common-name and synonym suggestions"));
    }
    if cache.dfam_names.is_none() {
        missing.push((dir.join("dfam-names.txt"), "ID collision check against Dfam"));
    }
    missing
}

fn load_line_set(path: PathBuf) -> Option<HashSet<String>> {
    let f = std::fs::File::open(&path).ok()?;
    let reader = std::io::BufReader::new(f);
    let mut set = HashSet::new();
    for line in reader.lines() {
        if let Ok(line) = line {
            let trimmed = line.trim().to_string();
            if !trimmed.is_empty() && !trimmed.starts_with('#') {
                set.insert(trimmed);
            }
        }
    }
    Some(set)
}

/// Load `taxonomy-common.tsv` as a lowercase-key map → scientific name.
///
/// File format: `<alternate_name>\t<scientific_name>` per line.
/// Keys are lowercased so callers can use `query.to_lowercase()` for lookup.
fn load_common_map(path: PathBuf) -> Option<HashMap<String, String>> {
    let f = std::fs::File::open(&path).ok()?;
    let reader = std::io::BufReader::new(f);
    let mut map = HashMap::new();
    for line in reader.lines() {
        if let Ok(line) = line {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            if let Some((alt, sci)) = trimmed.split_once('\t') {
                map.insert(alt.trim().to_lowercase(), sci.trim().to_string());
            }
        }
    }
    Some(map)
}
