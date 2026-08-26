//! Calibrate the repetitive-input probe across a corpus of families.
//!
//! Reports the distribution of self-HSP counts so the threshold can be set from
//! real families rather than from a synthetic control. Prints one line per
//! family over the threshold, then a summary.
//!
//! ```sh
//! cargo run -p cons-core --release --example probe_corpus -- matrix.txt f1.fa f2.fa ...
//! ```

use aln_core::{io, SubstMatrix};
use aln_engine::{AlignMode, AlignParams, SearchParams};
use aln_rmblast::{RmblastEngine, RmblastOptions};
use cons_core::probe::{probe_repetitive, ProbeConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let matrix_path = args.next().ok_or("usage: probe_corpus <matrix> <family.fa>...")?;
    let files: Vec<String> = args.collect();
    if files.is_empty() {
        return Err("no family files given".into());
    }
    let matrix = SubstMatrix::from_file(&matrix_path)?;
    let align = AlignParams {
        mode: AlignMode::Local,
        min_score: 150,
        traceback: true,
        ..AlignParams::from_matrix(&matrix)
    };
    let search = SearchParams {
        matrix: Some(matrix),
        gap_init: -(align.gap_open as i32),
        ins_gap_ext: -(align.gap_extend as i32),
        del_gap_ext: -(align.gap_extend as i32),
        min_match: 7,
        min_score: align.min_score.max(0),
        mask_level: 101,
        cores: Some(1),
        score_mode: if std::env::var("PROBE_CA").is_ok() {
            aln_engine::engine::ScoreMode::ComplexityAdjusted
        } else {
            aln_engine::engine::ScoreMode::Basic
        },
        ..Default::default()
    };
    // DUST **off**, deliberately. DUST masks low-complexity query sequence, so
    // with it on a pure tandem repeat self-aligns to nothing and the probe
    // reports 0 — it erases exactly the signal being measured. Measured: the
    // `(GGGAGG)n` synthetic scores 0 with DUST on and 133 with it off.
    let engine = RmblastEngine::new(search, RmblastOptions { dust: false, ..Default::default() })?;
    let mut cfg = ProbeConfig::default();
    if let Ok(v) = std::env::var("PROBE_MAX_HSPS") {
        cfg.max_hsps = v.parse().unwrap_or(cfg.max_hsps);
    }

    let t0 = std::time::Instant::now();
    let mut worsts: Vec<(usize, String)> = Vec::new();
    let mut flagged_families = 0usize;
    let mut total_flagged_instances = 0usize;
    let mut skipped = 0usize;

    for f in &files {
        let seqs = match io::read_fasta_file(f) {
            Ok(s) if !s.is_empty() => s,
            _ => {
                skipped += 1;
                continue;
            }
        };
        let rep = probe_repetitive(&engine, &seqs, &cfg)?;
        let name = std::path::Path::new(f)
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| f.clone());
        if !rep.flagged.is_empty() {
            flagged_families += 1;
            total_flagged_instances += rep.flagged.len();
            println!(
                "FLAG {name}: {} of {} probed over {} (worst {} self-HSPs, median {})",
                rep.flagged.len(),
                rep.probed.len(),
                cfg.max_hsps,
                rep.worst(),
                rep.median()
            );
        }
        worsts.push((rep.worst(), name));
    }

    worsts.sort_by_key(|w| std::cmp::Reverse(w.0));
    let n = worsts.len();
    let counts: Vec<usize> = worsts.iter().map(|w| w.0).collect();
    let mut asc = counts.clone();
    asc.sort_unstable();
    let pct = |p: f64| asc[((n as f64 - 1.0) * p) as usize];

    println!("\n=== {n} families probed in {:.1}s ({skipped} unreadable) ===", t0.elapsed().as_secs_f64());
    println!("  window {} bp, sample {}, threshold >{}", cfg.window, cfg.sample, cfg.max_hsps);
    println!("  worst-instance self-HSPs per family:");
    println!("    median {}  p90 {}  p99 {}  max {}", pct(0.50), pct(0.90), pct(0.99), asc[n - 1]);
    println!("  families flagged: {flagged_families} ({:.1}%), instances flagged: {total_flagged_instances}",
             100.0 * flagged_families as f64 / n as f64);
    println!("  top 10 families by worst instance:");
    for (c, name) in worsts.iter().take(10) {
        println!("    {name:<44} {c}");
    }
    Ok(())
}
