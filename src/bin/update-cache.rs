/// update-cache — download and populate the stk-lint tier-2 validation cache.
///
/// Populates `~/.cache/stk/` (or `--cache-dir`) with:
///   taxonomy.tsv        — NCBI scientific taxon names
///   taxonomy-common.tsv — common names / synonyms → scientific name
///   classification.tsv  — Dfam TP classification strings (short and long forms)
///   dfam-names.txt      — Dfam family names
///
/// Requires: curl(1) and unzip(1) on the PATH.
use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use dfam_curator::dfam::cache::cache_dir;

const TAXONOMY_URL: &str = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip";
const DFAM_CLASS_URL: &str =
    "https://www.dfam.org/releases/current/infrastructure/class_ns.tsv";
const DFAM_NAMES_URL: &str =
    "https://www.dfam.org/releases/current/infrastructure/name_ns.tsv";

/// Name classes from names.dmp that we index for common-name lookup.
/// "scientific name" is handled separately; synonyms and misnomers are excluded.
const ALTERNATE_CLASSES: &[&str] = &[
    "common name",
    "genbank common name",
    "blast name",
    "equivalent name",
];

#[derive(Parser, Debug)]
#[command(
    name = "update-cache",
    about = "Populate the stk-lint tier-2 validation cache",
    version
)]
struct Args {
    /// Override the cache directory.
    /// Defaults to $STK_CACHE_DIR, then ~/.cache/stk.
    #[arg(long)]
    cache_dir: Option<PathBuf>,

    /// Re-download / regenerate files even if they already exist.
    #[arg(long)]
    force: bool,

    #[command(subcommand)]
    command: Option<Cmd>,
}

#[derive(Subcommand, Debug)]
enum Cmd {
    /// Download NCBI taxonomy and write taxonomy.tsv + taxonomy-common.tsv.
    Taxonomy,
    /// Download Dfam TP classifications and write classification.tsv.
    Classifications,
    /// Download Dfam family names and write dfam-names.txt.
    Names,
    /// Show cache status.
    Info,
    /// Download all data, then show cache status. [default]
    All,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let cdir = args.cache_dir.unwrap_or_else(cache_dir);

    std::fs::create_dir_all(&cdir)
        .with_context(|| format!("cannot create cache directory {}", cdir.display()))?;

    eprintln!("Cache directory: {}", cdir.display());

    match args.command {
        None | Some(Cmd::All) => {
            fetch_taxonomy(&cdir, args.force)?;
            fetch_dfam_classifications(&cdir, args.force)?;
            fetch_dfam_names(&cdir, args.force)?;
            print_info(&cdir);
        }
        Some(Cmd::Taxonomy)        => fetch_taxonomy(&cdir, args.force)?,
        Some(Cmd::Classifications) => fetch_dfam_classifications(&cdir, args.force)?,
        Some(Cmd::Names)           => fetch_dfam_names(&cdir, args.force)?,
        Some(Cmd::Info)            => print_info(&cdir),
    }

    Ok(())
}

// ── Taxonomy ──────────────────────────────────────────────────────────────────

fn fetch_taxonomy(cdir: &Path, force: bool) -> Result<()> {
    let sci_path    = cdir.join("taxonomy.tsv");
    let common_path = cdir.join("taxonomy-common.tsv");

    if sci_path.exists() && common_path.exists() && !force {
        eprintln!("taxonomy files already exist; skipping (use --force to re-download)");
        return Ok(());
    }

    eprintln!("Downloading NCBI taxonomy from {TAXONOMY_URL}");

    let tmp = tempfile::Builder::new()
        .suffix(".zip")
        .tempfile_in(cdir)
        .context("cannot create temp file in cache dir")?;
    let tmp_path = tmp.path().to_path_buf();

    let status = Command::new("curl")
        .args(["-fSL", "--progress-bar", TAXONOMY_URL, "-o"])
        .arg(&tmp_path)
        .status()
        .context("cannot run curl — is curl installed and on your PATH?")?;

    if !status.success() {
        bail!("curl failed with exit status {:?}", status.code());
    }

    eprintln!("Extracting names.dmp …");

    let out = Command::new("unzip")
        .args(["-p"])
        .arg(&tmp_path)
        .arg("names.dmp")
        .output()
        .context("cannot run unzip — is unzip installed and on your PATH?")?;

    if !out.status.success() {
        bail!("unzip failed: {}", String::from_utf8_lossy(&out.stderr).trim());
    }

    eprintln!("Parsing names.dmp (single pass) …");

    // names.dmp: fields separated by <TAB>|<TAB>, trailing <TAB>|
    //   tax_id | name_txt | unique_name | name_class |
    // Splitting on '\t' gives meaningful fields at indices 0, 2, 4, 6.
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

    // ── Write taxonomy.tsv ────────────────────────────────────────────────────
    {
        let mut w = File::create(&sci_path)
            .with_context(|| format!("cannot write {}", sci_path.display()))?;
        for name in sci_map.values() {
            writeln!(w, "{name}")?;
        }
        eprintln!("Wrote {} scientific names → {}", sci_map.len(), sci_path.display());
    }

    // ── Write taxonomy-common.tsv ─────────────────────────────────────────────
    // Format: <alternate_name>\t<scientific_name>   (original casing preserved)
    // Keys are lowercased when loaded into Cache for case-insensitive lookup.
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

    Ok(())
}

// ── Dfam data ─────────────────────────────────────────────────────────────────

/// Download `url` to `dest` atomically via a temp file.
fn fetch_url_to_file(url: &str, dest: &Path) -> Result<()> {
    let dir = dest.parent().unwrap_or(Path::new("."));
    let tmp = tempfile::Builder::new()
        .tempfile_in(dir)
        .context("cannot create temp file in cache dir")?;
    let tmp_path = tmp.path().to_path_buf();

    let status = Command::new("curl")
        .args(["-fSL", "--progress-bar", url, "-o"])
        .arg(&tmp_path)
        .status()
        .context("cannot run curl — is curl installed and on your PATH?")?;

    if !status.success() {
        bail!("curl failed with exit status {:?}", status.code());
    }

    std::fs::rename(&tmp_path, dest)
        .with_context(|| format!("cannot write {}", dest.display()))?;
    let _ = tmp.keep();

    Ok(())
}

/// Count non-empty, non-comment lines in a file.
fn count_data_lines(path: &Path) -> usize {
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

fn fetch_dfam_classifications(cdir: &Path, force: bool) -> Result<()> {
    let dest = cdir.join("classification.tsv");
    if dest.exists() && !force {
        eprintln!("classification.tsv already exists; skipping (use --force to re-download)");
        return Ok(());
    }
    eprintln!("Downloading Dfam classifications from {DFAM_CLASS_URL}");
    fetch_url_to_file(DFAM_CLASS_URL, &dest)?;
    let n = count_data_lines(&dest);
    eprintln!("Wrote {n} classification entries → {}", dest.display());
    Ok(())
}

fn fetch_dfam_names(cdir: &Path, force: bool) -> Result<()> {
    let dest = cdir.join("dfam-names.txt");
    if dest.exists() && !force {
        eprintln!("dfam-names.txt already exists; skipping (use --force to re-download)");
        return Ok(());
    }
    eprintln!("Downloading Dfam names from {DFAM_NAMES_URL}");
    fetch_url_to_file(DFAM_NAMES_URL, &dest)?;
    let n = count_data_lines(&dest);
    eprintln!("Wrote {n} Dfam family names → {}", dest.display());
    Ok(())
}

// ── Info ──────────────────────────────────────────────────────────────────────

fn print_info(cdir: &Path) {
    const FILES: &[(&str, &str)] = &[
        ("taxonomy.tsv",        "NCBI scientific taxon names"),
        ("taxonomy-common.tsv", "NCBI alternate name mappings"),
        ("classification.tsv",  "Dfam TP classifications"),
        ("dfam-names.txt",      "Dfam family names"),
    ];
    println!();
    println!("Cache directory: {}", cdir.display());
    for (file, desc) in FILES {
        let status = if cdir.join(file).exists() { "ok     " } else { "missing" };
        println!("  [{status}]  {file}  ({desc})");
    }
    println!();
    println!("Run 'update-cache' or 'update-cache all' to refresh all files.");
}
