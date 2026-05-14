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
use regex::Regex;
use std::collections::HashMap;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

use dfam_curator::{
    consensus::{build_consensus_from_sequences, ConsensusParams},
    dfam::{
        cache::{cache_dir, load_cache, missing_cache_files},
        edit::{apply_ops, Op},
        lint::{check_duplicate_ids, lint_record, Diagnostic, Severity},
        record::{iter_records, iter_records_raw, RawDfamRecord},
    },
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

    if args.genome.is_none() {
        eprintln!(
            "WARN  [genome] Without a reference assembly, sequence coordinates cannot be \
             verified.  Provide a FASTA or 2bit file with --genome for full validation."
        );
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
                        "sequence {:?} coordinates do not match reference within +/-3 bp",
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
    after_help = "Operations are applied in a fixed sequence: --delete, then --set, \
                  then --append, then --sub.  Running --sub last means it transforms \
                  values however they arrived — pre-existing, set, or appended.\n\n\
                  --sub EXPR format: /PATTERN/REPLACEMENT/[g]\n  \
                  The first character is the delimiter (any character works, e.g. |…|…|).\n  \
                  Omit the trailing /g to replace only the first match within each value;\n  \
                  add /g to replace every match.  For multi-valued fields (OC, CC) the \
                  substitution is applied to each line independently.\n  \
                  Capture groups use $1, $2, … in the replacement.\n\n\
                  Examples:\n  \
                  stk edit --set AU \"Hubley R\" families.stk\n  \
                  stk edit --delete SE --set DE \"Updated desc\" families.stk\n  \
                  stk edit --select MyFam --append OC \"Mus musculus\" families.stk\n  \
                  stk edit --select 3 --set AU \"Hubley R\" families.stk\n  \
                  stk edit --set AU \"Hubley R\" -o fixed.stk families.stk\n  \
                  stk edit --sub ID \"/^(.*)-$/$1/\" families.stk\n  \
                  stk edit --sub DE \"/foo/bar/g\" families.stk\n  \
                  stk edit --set ID \"new-\" --sub ID \"/^(.*)-$/$1/\" families.stk\n  \
                  stk edit --update-consensus families.stk\n  \
                  stk edit --select MyFam --update-consensus families.stk"
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

    /// Apply a regex substitution to all values of a GF tag.
    /// Format: /PATTERN/REPLACEMENT/[g]  (first char is delimiter; /g = replace all matches).
    /// Capture groups use $1, $2, … in the replacement.
    /// May be specified multiple times: --sub ID "/^(.*)-$/$1/" --sub DE "/old/new/"
    #[arg(long = "sub", value_names = ["TAG", "EXPR"], num_args = 2,
          action = clap::ArgAction::Append)]
    sub: Vec<String>,

    /// Recompute the #=GC RF consensus from the aligned sequences.
    /// Because this loads all sequences into MSA data structures it is more
    /// expensive than plain field edits, so it must be requested explicitly.
    /// The new consensus is written as the #=GC RF line; any existing RF value
    /// is replaced.
    #[arg(long)]
    update_consensus: bool,

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

/// Parse a substitution expression of the form `DELIM PATTERN DELIM REPLACEMENT DELIM [g]`.
/// The delimiter is the first character of the string (typically `/` or `|`).
/// Returns `(pattern, replacement, all)`.
fn parse_sub_expr(expr: &str) -> anyhow::Result<(Regex, String, bool)> {
    let mut chars = expr.chars();
    let delim = chars
        .next()
        .ok_or_else(|| anyhow::anyhow!("substitution expression is empty"))?;

    // Split the remainder on the delimiter.  We need exactly three parts:
    // pattern, replacement, flags (flags may be empty).
    let rest = &expr[delim.len_utf8()..];
    let parts: Vec<&str> = rest.splitn(3, delim).collect();
    if parts.len() != 3 {
        anyhow::bail!(
            "substitution expression {:?} must have the form {}PAT{}REPL{}[g]",
            expr, delim, delim, delim
        );
    }
    let (pattern_str, replacement, flags) = (parts[0], parts[1], parts[2]);

    for c in flags.chars() {
        if c != 'g' {
            anyhow::bail!("unknown flag {:?} in substitution {:?}", c, expr);
        }
    }
    let all = flags.contains('g');

    let pattern = Regex::new(pattern_str)
        .with_context(|| format!("invalid regex {:?} in substitution {:?}", pattern_str, expr))?;

    Ok((pattern, replacement.to_string(), all))
}

/// Rebuild the `#=GC RF` line from the aligned sequences in `record`.
///
/// All sequence rows are treated as instances and fed to the CpG-aware
/// consensus caller.  Does nothing when the record has no sequence rows.
fn update_rf_consensus(record: &mut RawDfamRecord) {
    let sequences = &record.sequences[..];
    if sequences.is_empty() {
        return;
    }

    // Normalize '.' -> '-' so the consensus caller treats them as gaps,
    // matching what the MultiAlign reader does when loading STK files.
    let owned: Vec<Vec<u8>> = sequences
        .iter()
        .map(|row| {
            row.aligned_seq
                .bytes()
                .map(|b| if b == b'.' { b'-' } else { b })
                .collect()
        })
        .collect();
    let seq_refs: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();

    let params = ConsensusParams::default();
    let consensus = build_consensus_from_sequences(&seq_refs, &params);

    // Use '.' for gap positions — Stockholm convention for the RF line.
    let cons_str: String = consensus.into_iter()
        .map(|b| if b == b'-' { '.' } else { b as char })
        .collect();
    record.gc.insert("RF".to_string(), cons_str);
}

fn run_edit(args: EditArgs) -> anyhow::Result<()> {
    // Build the operation list.  --set, --append, and --sub consume 2 values each.
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
    for pair in args.sub.chunks_exact(2) {
        let (pattern, replacement, all) = parse_sub_expr(&pair[1])?;
        ops.push(Op::Sub { tag: pair[0].clone(), pattern, replacement, all });
    }

    if ops.is_empty() && !args.update_consensus {
        anyhow::bail!(
            "no edit operations given; \
             specify at least one of --set, --delete, --append, --sub, or --update-consensus"
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

        for result in iter_records_raw(BufReader::new(f)) {
            let mut record = result
                .with_context(|| format!("parse error in {}", path.display()))?;

            let selected = args
                .select
                .as_deref()
                .map(|sel| record_selected(&record, sel))
                .unwrap_or(true);

            if selected {
                apply_ops(&mut record, &ops);
                if args.update_consensus {
                    update_rf_consensus(&mut record);
                }
            }

            record.write_to(&mut out)?;
        }
    }

    out.flush()?;
    Ok(())
}
