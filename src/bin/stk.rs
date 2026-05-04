/// stk — Dfam Stockholm file toolkit.
///
/// Subcommands:
///   lint   Validate one or more STK files and report diagnostics.
///   edit   Modify #=GF annotation fields across records.
///
/// Exit status for `lint`: 0 = clean, 1 = at least one ERROR, 2 = I/O failure.
/// Exit status for `edit`: 0 = success, 2 = I/O failure.
use anyhow::Context;
use clap::{Args, Parser, Subcommand, ValueEnum};
use std::collections::HashMap;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

use dfam_curator::dfam::{
    cache::{cache_dir, load_cache, missing_cache_files},
    edit::{apply_ops, Op},
    lint::{check_duplicate_ids, lint_record, Diagnostic, Severity},
    record::{iter_records, RawDfamRecord},
};

// ── Top-level CLI ─────────────────────────────────────────────────────────────

#[derive(Parser)]
#[command(
    name = "stk",
    about = "Dfam Stockholm file toolkit",
    version,
    subcommand_required = true,
    arg_required_else_help = true,
)]
struct Cli {
    #[command(subcommand)]
    command: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Validate STK files and report structural and semantic diagnostics.
    Lint(LintArgs),
    /// Edit #=GF annotation fields in STK files.
    Edit(EditArgs),
    /// Extract one or more records from a multi-record STK file.
    Extract(ExtractArgs),
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Cmd::Lint(args) => run_lint(args),
        Cmd::Edit(args) => run_edit(args),
        Cmd::Extract(args) => run_extract(args),
    }
}

// ── lint subcommand ───────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, ValueEnum)]
enum MinSev {
    Error,
    Warn,
    Info,
}

impl From<MinSev> for Severity {
    fn from(m: MinSev) -> Self {
        match m {
            MinSev::Error => Severity::Error,
            MinSev::Warn  => Severity::Warn,
            MinSev::Info  => Severity::Info,
        }
    }
}

#[derive(Args, Debug)]
struct LintArgs {
    /// One or more Stockholm files to check.
    #[arg(required = true)]
    input: Vec<PathBuf>,

    /// Override the cache directory for tier-2 validation data.
    /// Defaults to $STK_CACHE_DIR, then ~/.cache/stk.
    #[arg(long)]
    cache_dir: Option<PathBuf>,

    /// Minimum severity to report.
    #[arg(long, value_enum, default_value = "info")]
    min_severity: MinSev,

    /// Suppress informational notices about missing tier-2 cache files.
    #[arg(long)]
    no_cache_warn: bool,

    /// Reference genome (FASTA or .2bit) for coordinate validation.
    /// When provided, each sequence row's coordinates are checked against the
    /// reference.  Fixable issues (half-open intervals, small shifts, wrong
    /// strand) are reported as WARN; coordinates that cannot be located within
    /// ±3 bp are also reported as WARN.
    #[arg(long, value_name = "FILE")]
    genome: Option<PathBuf>,
}

fn run_lint(args: LintArgs) -> anyhow::Result<()> {
    let min_sev: Severity = args.min_severity.into();
    let cdir = args.cache_dir.unwrap_or_else(cache_dir);
    let cache = load_cache(&cdir);

    if !args.no_cache_warn {
        for (path, desc) in missing_cache_files(&cache, &cdir) {
            eprintln!(
                "INFO  [cache] {:?} not found — {} will be skipped \
                 (run update-cache to populate)",
                path, desc
            );
        }
    }

    // Load reference genome once if --genome was supplied.
    let genome_map: Option<HashMap<String, Vec<u8>>> = match &args.genome {
        None => None,
        Some(genome_path) => {
            let s = genome_path
                .to_str()
                .ok_or_else(|| anyhow::anyhow!("genome path is not valid UTF-8"))?;
            eprintln!("Loading reference genome: {}", s);
            Some(dfam_coord::load_reference(s)
                .with_context(|| format!("cannot load genome {:?}", genome_path))?)
        }
    };

    let cache_ref = Some(&cache);
    let mut any_error = false;

    for path in &args.input {
        let filename = path.display().to_string();
        let f = std::fs::File::open(path)
            .with_context(|| format!("cannot open {}", path.display()))?;
        let mut records: Vec<RawDfamRecord> = Vec::new();

        for result in iter_records(BufReader::new(f)) {
            let record = result
                .with_context(|| format!("parse error in {}", path.display()))?;
            let label = record.label();

            let mut all_diags = lint_record(&record, cache_ref);
            if let Some(ref gmap) = genome_map {
                all_diags.extend(coord_lint(&record, gmap));
            }

            for d in &all_diags {
                if d.severity >= min_sev {
                    if d.severity == Severity::Error {
                        any_error = true;
                    }
                    println!(
                        "{}\t{}\t{}\t{}\t{}",
                        filename, label, d.severity, d.check, d.message
                    );
                }
            }
            records.push(record);
        }

        for d in check_duplicate_ids(&records) {
            if d.severity >= min_sev {
                println!(
                    "{}\tFILE\t{}\t{}\t{}",
                    filename, d.severity, d.check, d.message
                );
            }
        }
    }

    std::process::exit(if any_error { 1 } else { 0 });
}

/// Run DisCoord phase-1 coordinate validation on the sequences in one record.
///
/// Reports WARN diagnostics for:
/// - `seq_coord_fixed`: coordinates were fixable (half-open, small shift, or
///   wrong strand).  The status string describes the repair applied.
/// - `seq_coord_invalid`: the sequence_id is present in the reference but the
///   sequence could not be located within ±3 bp in either orientation.
///
/// Sequences whose `sequence_id` is absent from the reference are silently
/// skipped (they likely belong to a different assembly).
fn coord_lint(
    record: &RawDfamRecord,
    genome_map: &HashMap<String, Vec<u8>>,
) -> Vec<Diagnostic> {
    let mut seq_records = dfam_coord::records_from_rows(&record.sequences, "");
    dfam_coord::validate_sequences(&mut seq_records, genome_map, false);

    let mut diags = Vec::new();
    for sr in &seq_records {
        let orig = sr.original_id.as_deref().unwrap_or("?");
        match sr.validated.as_deref() {
            Some("valid") => {}
            Some(status) if status.starts_with("fixed") => {
                diags.push(Diagnostic {
                    severity: Severity::Warn,
                    check: "seq_coord_fixed",
                    message: format!("sequence {:?} coordinates repaired: {}", orig, status),
                });
            }
            None if genome_map.contains_key(&sr.sequence_id) => {
                diags.push(Diagnostic {
                    severity: Severity::Warn,
                    check: "seq_coord_invalid",
                    message: format!(
                        "sequence {:?} coordinates do not match reference within ±3 bp",
                        orig
                    ),
                });
            }
            _ => {}
        }
    }
    diags
}

// ── extract subcommand ────────────────────────────────────────────────────────

#[derive(Args, Debug)]
#[command(
    after_help = "Examples:\n  \
                  stk extract --select MyFam families.stk\n  \
                  stk extract --select 3 families.stk\n  \
                  stk extract --select MyFam -o MyFam.stk families.stk"
)]
struct ExtractArgs {
    /// Record to extract.
    /// A purely numeric value selects by 1-based record number (e.g. --select 3).
    /// Any non-numeric value selects by exact #=GF ID match (e.g. --select MyFam).
    #[arg(long, value_name = "SELECT", required = true)]
    select: String,

    /// Write output to FILE instead of stdout.
    #[arg(long, short = 'o', value_name = "FILE")]
    output: Option<PathBuf>,

    /// One or more Stockholm files to search.
    #[arg(required = true)]
    input: Vec<PathBuf>,
}

fn run_extract(args: ExtractArgs) -> anyhow::Result<()> {
    let mut out: Box<dyn Write> = match &args.output {
        None => Box::new(BufWriter::new(std::io::stdout())),
        Some(p) => Box::new(BufWriter::new(
            std::fs::File::create(p)
                .with_context(|| format!("cannot create {}", p.display()))?,
        )),
    };

    let mut found = false;
    for path in &args.input {
        let f = std::fs::File::open(path)
            .with_context(|| format!("cannot open {}", path.display()))?;

        for result in iter_records(BufReader::new(f)) {
            let record = result
                .with_context(|| format!("parse error in {}", path.display()))?;

            if record_selected(&record, &args.select) {
                record.write_to(&mut out)?;
                found = true;
            }
        }
    }

    out.flush()?;

    if !found {
        anyhow::bail!("no record matching {:?} found", args.select);
    }
    Ok(())
}

// ── edit subcommand ───────────────────────────────────────────────────────────

#[derive(Args, Debug)]
#[command(
    after_help = "Operations are applied in a fixed sequence: --delete first, \
                  then --set, then --append.  This ensures that \
                  '--delete AU --set AU \"New\"' always produces the expected result.\n\n\
                  Examples:\n  \
                  stk edit --set AU \"Hubley R\" families.stk\n  \
                  stk edit --delete SE --set DE \"Updated desc\" families.stk\n  \
                  stk edit --select MyFam --append OC \"Mus musculus\" families.stk\n  \
                  stk edit --select 3 --set AU \"Hubley R\" families.stk\n  \
                  stk edit --set AU \"Hubley R\" -o fixed.stk families.stk"
)]
struct EditArgs {
    /// Set (or add) a GF field, replacing any existing occurrences.
    /// May be specified multiple times: --set AU "Smith J" --set DE "New desc"
    #[arg(long = "set", value_names = ["TAG", "VALUE"], num_args = 2,
          action = clap::ArgAction::Append)]
    set: Vec<String>,

    /// Remove all occurrences of a GF tag.
    /// May be specified multiple times: --delete SE --delete TD
    #[arg(long = "delete", value_name = "TAG",
          action = clap::ArgAction::Append)]
    delete: Vec<String>,

    /// Append a new GF field (for multi-valued fields such as OC or CC).
    /// May be specified multiple times: --append OC "Mus musculus"
    #[arg(long = "append", value_names = ["TAG", "VALUE"], num_args = 2,
          action = clap::ArgAction::Append)]
    append: Vec<String>,

    /// Only edit records matching SELECT.
    /// A purely numeric value selects by 1-based record number (e.g. --select 3).
    /// Any non-numeric value selects by exact #=GF ID match (e.g. --select MyFam).
    /// Records that do not match are passed through unchanged.
    ///
    /// NOTE: because numeric strings are always interpreted as record numbers,
    /// IDs that consist solely of digits cannot be targeted by name with this flag.
    /// stk lint will report such IDs as an error.
    #[arg(long, value_name = "SELECT")]
    select: Option<String>,

    /// Write output to FILE instead of stdout.
    #[arg(long, short = 'o', value_name = "FILE")]
    output: Option<PathBuf>,

    /// One or more Stockholm files to edit.
    #[arg(required = true)]
    input: Vec<PathBuf>,
}

/// Return `true` if `record` matches the `--select` value.
///
/// A purely numeric string is interpreted as a 1-based record number within
/// the file.  Any non-numeric string is matched against the record's `#=GF ID`
/// field.
///
/// This means IDs that consist solely of digits can never be targeted by name;
/// `stk lint` reports such IDs as an error (`id_numeric`) to prevent ambiguity.
fn record_selected(record: &RawDfamRecord, select: &str) -> bool {
    if let Ok(n) = select.parse::<usize>() {
        record.record_num == n
    } else {
        record.gf_first("ID").map(str::trim) == Some(select)
    }
}

fn run_edit(args: EditArgs) -> anyhow::Result<()> {
    // Build the operation list.  --set and --append consume 2 values each.
    let mut ops: Vec<Op> = Vec::new();
    for pair in args.delete.iter() {
        ops.push(Op::Delete { tag: pair.clone() });
    }
    for pair in args.set.chunks_exact(2) {
        ops.push(Op::Set { tag: pair[0].clone(), value: pair[1].clone() });
    }
    for pair in args.append.chunks_exact(2) {
        ops.push(Op::Append { tag: pair[0].clone(), value: pair[1].clone() });
    }

    if ops.is_empty() {
        anyhow::bail!(
            "no edit operations given; \
             specify at least one of --set, --delete, or --append"
        );
    }

    let mut out: Box<dyn Write> = match &args.output {
        None => Box::new(BufWriter::new(std::io::stdout())),
        Some(p) => Box::new(BufWriter::new(
            std::fs::File::create(p)
                .with_context(|| format!("cannot create {}", p.display()))?,
        )),
    };

    for path in &args.input {
        let f = std::fs::File::open(path)
            .with_context(|| format!("cannot open {}", path.display()))?;

        for result in iter_records(BufReader::new(f)) {
            let mut record = result
                .with_context(|| format!("parse error in {}", path.display()))?;

            let selected = args
                .select
                .as_deref()
                .map(|sel| record_selected(&record, sel))
                .unwrap_or(true);

            if selected {
                apply_ops(&mut record, &ops);
            }

            record.write_to(&mut out)?;
        }
    }

    out.flush()?;
    Ok(())
}
