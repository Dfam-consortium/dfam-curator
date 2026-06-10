/// update-cache — download and populate the stk-lint tier-2 validation cache.
///
/// Populates `~/.cache/stk/` (or `--cache-dir`) with:
///   taxonomy.tsv        — NCBI scientific taxon names
///   taxonomy-common.tsv — common names / synonyms → scientific name
///   classification.tsv  — Dfam TP classification strings (short and long forms)
///   dfam-names.txt      — Dfam family names
///
/// Without --force, each file is fetched with curl -z so only changed files
/// are downloaded.  With --force, all files are downloaded unconditionally.
///
/// Requires: curl(1) and unzip(1) on the PATH.
use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use std::path::{Path, PathBuf};

use dfam_curator::dfam::{
    cache::cache_dir,
    fetch::{
        count_data_lines, fetch_dfam_file, fetch_taxonomy,
        DFAM_CLASS_URL, DFAM_NAMES_URL, TAXONOMY_URL,
    },
};

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

    /// Re-download all files unconditionally, ignoring curl's conditional
    /// GET (-z).
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
            run_taxonomy(&cdir, args.force)?;
            run_classifications(&cdir, args.force)?;
            run_names(&cdir, args.force)?;
            print_info(&cdir);
        }
        Some(Cmd::Taxonomy)        => run_taxonomy(&cdir, args.force)?,
        Some(Cmd::Classifications) => run_classifications(&cdir, args.force)?,
        Some(Cmd::Names)           => run_names(&cdir, args.force)?,
        Some(Cmd::Info)            => print_info(&cdir),
    }

    Ok(())
}

fn run_taxonomy(cdir: &Path, force: bool) -> Result<()> {
    eprintln!("Checking NCBI taxonomy from {TAXONOMY_URL}");
    if !force {
        eprintln!("(~60 MB — only downloaded if the server file has changed)");
    }
    match fetch_taxonomy(cdir, force)? {
        true  => eprintln!("Taxonomy files updated."),
        false => eprintln!("Taxonomy files are already up to date."),
    }
    Ok(())
}

fn run_classifications(cdir: &Path, force: bool) -> Result<()> {
    let dest = cdir.join("classification.tsv");
    eprintln!("Checking Dfam classifications from {DFAM_CLASS_URL}");
    match fetch_dfam_file(DFAM_CLASS_URL, &dest, force)? {
        true  => eprintln!(
            "Wrote {} classification entries → {}",
            count_data_lines(&dest),
            dest.display()
        ),
        false => eprintln!("classification.tsv is already up to date."),
    }
    Ok(())
}

fn run_names(cdir: &Path, force: bool) -> Result<()> {
    let dest = cdir.join("dfam-names.txt");
    eprintln!("Checking Dfam names from {DFAM_NAMES_URL}");
    match fetch_dfam_file(DFAM_NAMES_URL, &dest, force)? {
        true  => eprintln!(
            "Wrote {} Dfam family names → {}",
            count_data_lines(&dest),
            dest.display()
        ),
        false => eprintln!("dfam-names.txt is already up to date."),
    }
    Ok(())
}

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
