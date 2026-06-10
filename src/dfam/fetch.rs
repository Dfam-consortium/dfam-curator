//! Low-level fetch helpers shared by `update-cache` and the auto-refresh
//! logic in `refresh_cache`.
//!
//! Requires `curl(1)` (and `unzip(1)` for taxonomy) on PATH.

use anyhow::{bail, Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::Command;
use std::time::{Duration, SystemTime};

pub const CACHE_MAX_AGE_DAYS: u64 = 60;

pub const TAXONOMY_URL: &str = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip";
pub const DFAM_CLASS_URL: &str =
    "https://www.dfam.org/releases/current/infrastructure/class_ns.tsv";
pub const DFAM_NAMES_URL: &str =
    "https://www.dfam.org/releases/current/infrastructure/name_ns.tsv";

const ALTERNATE_CLASSES: &[&str] = &[
    "common name",
    "genbank common name",
    "blast name",
    "equivalent name",
];

/// Returns true if `path` is missing or its mtime exceeds `max_age_days`.
pub fn is_stale(path: &Path, max_age_days: u64) -> bool {
    match path.metadata().and_then(|m| m.modified()) {
        Ok(mtime) => {
            SystemTime::now()
                .duration_since(mtime)
                .unwrap_or(Duration::MAX)
                > Duration::from_secs(max_age_days * 86_400)
        }
        Err(_) => true,
    }
}

fn touch(path: &Path) {
    let _ = Command::new("touch").arg(path).status();
}

/// Fetch a single Dfam TSV file from `url` to `dest` atomically.
///
/// If `force` is false and `dest` already exists, passes `-z dest` to curl so
/// the download is skipped when the server file is unchanged (HTTP 304).
/// On HTTP 304, `dest`'s mtime is updated so the 60-day clock resets.
///
/// Returns `true` if new content was written, `false` on HTTP 304.
pub fn fetch_dfam_file(url: &str, dest: &Path, force: bool) -> Result<bool> {
    let dir = dest.parent().unwrap_or(Path::new("."));
    let tmp = tempfile::Builder::new()
        .tempfile_in(dir)
        .context("cannot create temp file in cache dir")?;
    let tmp_path = tmp.path().to_path_buf();

    let mut cmd = Command::new("curl");
    cmd.args(["-fSL", "--progress-bar", url, "-o"]).arg(&tmp_path);
    if !force && dest.exists() {
        cmd.arg("-z").arg(dest);
    }

    let status = cmd.status().context("cannot run curl — is it on your PATH?")?;
    if !status.success() {
        bail!("curl failed with exit status {:?}", status.code());
    }

    if tmp_path.metadata().map(|m| m.len()).unwrap_or(0) > 0 {
        std::fs::rename(&tmp_path, dest)
            .with_context(|| format!("cannot write {}", dest.display()))?;
        let _ = tmp.keep();
        Ok(true)
    } else {
        // HTTP 304 Not Modified — bump mtime so the 60-day clock resets
        if dest.exists() {
            touch(dest);
        }
        Ok(false)
    }
}

/// Download the NCBI taxonomy zip and write `taxonomy.tsv` + `taxonomy-common.tsv`
/// under `cdir`.
///
/// If `force` is false and `taxonomy.tsv` exists, its mtime is used as
/// `If-Modified-Since` so the ~60 MB download is skipped on HTTP 304.
///
/// Returns `true` if the files were (re)written, `false` on HTTP 304.
pub fn fetch_taxonomy(cdir: &Path, force: bool) -> Result<bool> {
    let sci_path    = cdir.join("taxonomy.tsv");
    let common_path = cdir.join("taxonomy-common.tsv");

    let tmp = tempfile::Builder::new()
        .suffix(".zip")
        .tempfile_in(cdir)
        .context("cannot create temp file in cache dir")?;
    let tmp_path = tmp.path().to_path_buf();

    let mut cmd = Command::new("curl");
    cmd.args(["-fSL", "--progress-bar", TAXONOMY_URL, "-o"]).arg(&tmp_path);
    if !force && sci_path.exists() {
        cmd.arg("-z").arg(&sci_path);
    }

    let status = cmd.status().context("cannot run curl — is it on your PATH?")?;
    if !status.success() {
        bail!("curl failed with exit status {:?}", status.code());
    }

    if tmp_path.metadata().map(|m| m.len()).unwrap_or(0) == 0 {
        touch(&sci_path);
        touch(&common_path);
        return Ok(false);
    }

    eprintln!("Extracting names.dmp ...");
    let out = Command::new("unzip")
        .args(["-p"])
        .arg(&tmp_path)
        .arg("names.dmp")
        .output()
        .context("cannot run unzip — is it on your PATH?")?;

    if !out.status.success() {
        bail!("unzip failed: {}", String::from_utf8_lossy(&out.stderr).trim());
    }

    eprintln!("Parsing names.dmp ...");

    let mut sci_map: HashMap<u32, String> = HashMap::new();
    let mut alt_list: Vec<(u32, String)>  = Vec::new();

    for line in BufReader::new(out.stdout.as_slice()).lines() {
        let line = line.context("I/O error reading unzip output")?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 7 { continue; }
        let Ok(id) = f[0].trim().parse::<u32>() else { continue; };
        let name  = f[2].trim();
        let class = f[6].trim();
        if class == "scientific name" {
            sci_map.insert(id, name.to_string());
        } else if ALTERNATE_CLASSES.contains(&class) {
            alt_list.push((id, name.to_string()));
        }
    }

    {
        let mut w = File::create(&sci_path)
            .with_context(|| format!("cannot write {}", sci_path.display()))?;
        for name in sci_map.values() {
            writeln!(w, "{name}")?;
        }
        eprintln!("Wrote {} scientific names → {}", sci_map.len(), sci_path.display());
    }

    {
        let mut w = File::create(&common_path)
            .with_context(|| format!("cannot write {}", common_path.display()))?;
        let mut n = 0usize;
        for (id, alt) in &alt_list {
            if let Some(sci) = sci_map.get(id) {
                writeln!(w, "{alt}\t{sci}")?;
                n += 1;
            }
        }
        eprintln!("Wrote {n} alternate name mappings → {}", common_path.display());
    }

    Ok(true)
}

/// Count non-empty, non-comment lines in a file (used by `update-cache info`).
pub fn count_data_lines(path: &Path) -> usize {
    std::fs::File::open(path)
        .ok()
        .map(|f| {
            BufReader::new(f)
                .lines()
                .filter_map(|l| l.ok())
                .filter(|l| {
                    let t = l.trim();
                    !t.is_empty() && !t.starts_with('#')
                })
                .count()
        })
        .unwrap_or(0)
}
