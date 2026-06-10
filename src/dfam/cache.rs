/// Cache management for tier-2 validation data.
///
/// Four optional data files live under a single cache directory:
///   - `taxonomy.tsv`        — one NCBI scientific taxon name per line
///   - `taxonomy-common.tsv` — common/alternate name TAB scientific name (for suggestions)
///   - `classification.tsv`  — Dfam TP classification TSV: short_form TAB long_form
///                             Both forms are accepted in the TP field.
///   - `dfam-names.txt`      — one Dfam family name per line (case-insensitive)
///
/// The cache directory defaults to `$STK_CACHE_DIR` → `~/.cache/stk`.
/// Each file is optional; absent files cause the corresponding tier-2 check
/// to be skipped (with an informational message from the caller).
use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::path::{Path, PathBuf};

use crate::dfam::fetch::{
    fetch_dfam_file, fetch_taxonomy, is_stale,
    CACHE_MAX_AGE_DAYS, DFAM_CLASS_URL, DFAM_NAMES_URL,
};

/// Controls how `refresh_cache` decides whether to fetch each file.
pub enum RefreshMode {
    /// Check mtime; if older than `CACHE_MAX_AGE_DAYS`, use `curl -z`.
    Auto,
    /// Always force-download unconditionally, ignoring age and `curl -z`.
    Force,
}

/// Attempt to refresh cache files according to `mode`.
///
/// Creates the cache directory if it does not exist. Network and fetch
/// failures are reported to stderr and do not propagate — callers continue
/// with whatever cached data is available.
pub fn refresh_cache(dir: &Path, mode: RefreshMode) {
    if let Err(e) = std::fs::create_dir_all(dir) {
        eprintln!(
            "WARN  [cache] cannot create cache directory {}: {e}",
            dir.display()
        );
        return;
    }

    let force = matches!(mode, RefreshMode::Force);

    refresh_one(
        "classification.tsv",
        || fetch_dfam_file(DFAM_CLASS_URL, &dir.join("classification.tsv"), force),
        force || is_stale(&dir.join("classification.tsv"), CACHE_MAX_AGE_DAYS),
        force,
    );

    refresh_one(
        "dfam-names.txt",
        || fetch_dfam_file(DFAM_NAMES_URL, &dir.join("dfam-names.txt"), force),
        force || is_stale(&dir.join("dfam-names.txt"), CACHE_MAX_AGE_DAYS),
        force,
    );

    // taxonomy.tsv + taxonomy-common.tsv share one download
    let tax_stale = force
        || is_stale(&dir.join("taxonomy.tsv"), CACHE_MAX_AGE_DAYS)
        || is_stale(&dir.join("taxonomy-common.tsv"), CACHE_MAX_AGE_DAYS);
    if tax_stale {
        if force {
            eprintln!("Cache: taxonomy files — force-downloading (NCBI taxonomy is ~60 MB) ...");
        } else {
            eprintln!(
                "Cache: taxonomy files are stale or missing — checking for update \
                 (may take a moment, NCBI taxonomy is ~60 MB) ..."
            );
        }
        match fetch_taxonomy(dir, force) {
            Ok(true)  => eprintln!("Cache: taxonomy files updated."),
            Ok(false) => eprintln!("Cache: taxonomy files are already up to date."),
            Err(e)    => eprintln!("WARN  [cache] failed to refresh taxonomy: {e}"),
        }
    }
}

fn refresh_one<F>(name: &str, fetch: F, should: bool, force: bool)
where
    F: FnOnce() -> anyhow::Result<bool>,
{
    if !should { return; }
    if force {
        eprintln!("Cache: {name} — force-downloading ...");
    } else {
        eprintln!("Cache: {name} is stale or missing — checking for update ...");
    }
    match fetch() {
        Ok(true)  => eprintln!("Cache: {name} updated."),
        Ok(false) => eprintln!("Cache: {name} is already up to date."),
        Err(e)    => eprintln!("WARN  [cache] failed to refresh {name}: {e}"),
    }
}

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
        classification:  load_classification_tsv(dir.join("classification.tsv")),
        taxonomy:        load_line_set(dir.join("taxonomy.tsv")),
        taxonomy_common: load_common_map(dir.join("taxonomy-common.tsv")),
        dfam_names:      load_names_lowercase(dir.join("dfam-names.txt")),
    }
}

/// Report which tier-2 checks are unavailable due to missing cache files.
/// Returns a list of `(filename, description)` pairs.
pub fn missing_cache_files(cache: &Cache, dir: &Path) -> Vec<(PathBuf, &'static str)> {
    let mut missing = Vec::new();
    if cache.classification.is_none() {
        missing.push((dir.join("classification.tsv"), "TP classification validation"));
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

/// Load `classification.tsv`: two-column TSV where col1 is the short form
/// (e.g. `DNA/TIR`) and col2 is the long semicolon-delimited form
/// (e.g. `Interspersed_Repeat;Transposable_Element;DNA_Transposon;TIR`).
/// Both forms are inserted into the set so either is accepted in the TP field.
/// Lines with only one column are accepted as-is (backwards compat).
fn load_classification_tsv(path: PathBuf) -> Option<HashSet<String>> {
    let f = std::fs::File::open(&path).ok()?;
    let reader = std::io::BufReader::new(f);
    let mut set = HashSet::new();
    for line in reader.lines() {
        if let Ok(line) = line {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            if let Some((short, long)) = trimmed.split_once('\t') {
                let short = short.trim();
                let long = long.trim();
                if !short.is_empty() { set.insert(short.to_string()); }
                if !long.is_empty()  { set.insert(long.to_string()); }
            } else {
                set.insert(trimmed.to_string());
            }
        }
    }
    Some(set)
}

/// Load `dfam-names.txt` with all names lowercased for case-insensitive lookup.
fn load_names_lowercase(path: PathBuf) -> Option<HashSet<String>> {
    let f = std::fs::File::open(&path).ok()?;
    let reader = std::io::BufReader::new(f);
    let mut set = HashSet::new();
    for line in reader.lines() {
        if let Ok(line) = line {
            let trimmed = line.trim();
            if !trimmed.is_empty() && !trimmed.starts_with('#') {
                set.insert(trimmed.to_lowercase());
            }
        }
    }
    Some(set)
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
