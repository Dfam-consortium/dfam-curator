/// Pretty-print an MSA in the Perl Linup–compatible block format.
///
/// Output matches `MultAln::printAlignments(blockSize => 100, showCons => 1)` followed
/// by the Linup footer (Avg Kimura, INFO, Cons length, FASTA consensus).
use std::io::{self, Write};

use crate::alignment::{MultiAlign, Orientation};
use crate::kimura;

/// Write the alignment to `out` in Linup pretty-print format.
///
/// `consensus` is the full-width gapped consensus (same length as `msa.width()`).
/// `block_size` matches Perl's `blockSize` parameter (Linup uses 100).
pub fn write(
    msa: &MultiAlign,
    consensus: &[u8],
    out: &mut dyn Write,
    block_size: usize,
) -> io::Result<()> {
    let width = msa.width();
    if width == 0 {
        return Ok(());
    }
    let ref_row = match msa.reference() {
        Some(r) => r,
        None => return Ok(()),
    };

    // ── Layout constants ──────────────────────────────────────────────────────

    let ref_display = format!("ref:{}", ref_row.name);
    // Match Perl printAlignments: maxIDLen is computed from the raw name length
    // (without the "ref:" prefix).  The consensus row sets the floor at 9.
    // This means ref_display may be longer than max_id_len and is truncated
    // when printed, exactly as Perl does with its printf("%-*s", ...) format.
    let mut max_id_len = "consensus".len().max(ref_row.name.len());
    for inst in &msa.sequences[1..] {
        max_id_len = max_id_len.max(inst.name.len());
    }

    // The reference start coordinate.  `seq_start` is 1-based, so a value of 0
    // means "unset" (a consensus-style reference with no source coordinate); such
    // a reference is numbered from 1, matching the consensus row.
    let ref_start = ref_row.seq_start.max(1);

    // Widest coordinate across reference, consensus, and all instance start/end.
    let ref_ungapped = count_alpha(&ref_row.seq) as u64;
    let ref_last = ref_start.saturating_add(ref_ungapped).saturating_sub(1);
    let cons_ungapped = count_alpha(consensus) as u64;
    let mut max_coord: u64 = ref_start.max(ref_last).max(cons_ungapped);
    for inst in &msa.sequences[1..] {
        max_coord = max_coord.max(inst.seq_start).max(inst.seq_end);
    }
    let max_coord_len = max_coord.to_string().len().max(1);

    // ── Per-instance precomputed values ───────────────────────────────────────

    let n = msa.sequences.len().saturating_sub(1);

    // Perl's getAlignedStart/End: first and last ALPHABETIC (non-gap, non-space)
    // character column.  Leading/trailing gap characters are treated like spaces.
    let base_starts: Vec<usize> = (0..n)
        .map(|i| {
            msa.sequences[i + 1]
                .seq
                .iter()
                .position(|&b| b.is_ascii_alphabetic())
                .unwrap_or(width)
        })
        .collect();
    let base_ends: Vec<usize> = (0..n)
        .map(|i| {
            msa.sequences[i + 1]
                .seq
                .iter()
                .rposition(|&b| b.is_ascii_alphabetic())
                .unwrap_or(0)
        })
        .collect();

    // Sort by gapped alignment start (Perl: sorted by getAlignedStart).
    let mut sorted: Vec<usize> = (0..n).collect();
    sorted.sort_by_key(|&i| base_starts[i]);

    let mut line_id = vec![0usize; n];
    for (rank, &i) in sorted.iter().enumerate() {
        line_id[i] = rank + 1;
    }

    // Per-instance coordinate trackers.
    // Reverse-strand sequences count DOWN from seq_end.
    let mut coord: Vec<u64> = (0..n)
        .map(|i| {
            let inst = &msa.sequences[i + 1];
            if inst.orient == Orientation::Reverse {
                inst.seq_end
            } else {
                inst.seq_start
            }
        })
        .collect();

    // ── Emit blocks ───────────────────────────────────────────────────────────

    let mut ref_pos: u64 = ref_start;
    let mut cons_pos: u64 = 1;
    let mut line_start = 0usize;

    // Perl separator width formula.
    let sep_width = max_id_len + 2 + max_coord_len + 1 + block_size + 1 + max_coord_len + 1 + 7;
    let sep_line: String = "~".repeat(sep_width);

    while line_start < width {
        let line_end = (line_start + block_size - 1).min(width - 1);

        let cslice = &consensus[line_start..=line_end];
        let rslice = &ref_row.seq[line_start..=line_end];

        // ── Consensus row ─────────────────────────────────────────────────────
        // Perl: end = startPos + numLetters - 1  (gives startPos-1 when empty)
        let cl = count_alpha(cslice) as u64;
        let cons_end = if cl > 0 { cons_pos + cl - 1 } else { cons_pos.saturating_sub(1) };
        write_row(out, "consensus", max_id_len, cons_pos, max_coord_len, cslice, block_size, cons_end, None)?;
        if cl > 0 { cons_pos += cl; }

        // ── Diff indicator ────────────────────────────────────────────────────
        write!(out, "{}", " ".repeat(max_id_len + max_coord_len + 2))?;
        for (&c, &r) in cslice.iter().zip(rslice.iter()) {
            out.write_all(&[diff_char(c, r)])?;
        }
        writeln!(out)?;

        // ── Reference row ─────────────────────────────────────────────────────
        // Perl truncates the reference display name to max_id_len characters
        // (it uses %-*s which right-pads but does not widen; our name column
        // was sized from the raw name without "ref:" so the prefix may overflow).
        let ref_label = if ref_display.len() > max_id_len {
            &ref_display[..max_id_len]
        } else {
            ref_display.as_str()
        };
        let rl = count_alpha(rslice) as u64;
        let ref_end = if rl > 0 { ref_pos + rl - 1 } else { ref_pos.saturating_sub(1) };
        write_row(out, ref_label, max_id_len, ref_pos, max_coord_len, rslice, block_size, ref_end, None)?;
        if rl > 0 { ref_pos += rl; }

        // ── Separator ─────────────────────────────────────────────────────────
        writeln!(out, "{}", sep_line)?;

        // ── Instance rows ─────────────────────────────────────────────────────
        for &i in &sorted {
            let bs = base_starts[i];
            let be = base_ends[i];

            // Perl skip: next if start >= lineEnd; next if end <= lineStart
            // where start/end are the first/last alphabetic positions.
            if bs >= line_end || be <= line_start {
                continue;
            }

            // Build the block sequence using only the alphabetic region.
            // Positions outside [bs, be] are shown as spaces (Perl behaviour:
            // leading/trailing gaps are treated as padding, not alignment chars).
            let bslice = build_block(
                &msa.sequences[i + 1].seq,
                line_start,
                line_end,
                bs,
                be,
            );

            let inst = &msa.sequences[i + 1];
            let nl = count_upper(&bslice) as u64; // uppercase only (Perl: tr/A-Z/)

            // Displayed start: when empty (nl==0) Perl shows coord-1.
            let disp_start = if nl == 0 {
                coord[i].saturating_sub(1)
            } else {
                coord[i]
            };
            let end_adj = nl.saturating_sub(1);
            let disp_end = if inst.orient == Orientation::Reverse {
                disp_start.saturating_sub(end_adj)
            } else {
                disp_start + end_adj
            };

            write_row(out, &inst.name, max_id_len, disp_start, max_coord_len, &bslice, block_size, disp_end, Some(line_id[i]))?;

            if inst.orient == Orientation::Reverse {
                coord[i] = coord[i].saturating_sub(nl);
            } else {
                coord[i] += nl;
            }
        }

        writeln!(out)?;
        writeln!(out)?;

        line_start = line_end + 1;
    }

    // ── Footer ────────────────────────────────────────────────────────────────

    let instances: Vec<&[u8]> = msa.sequences[1..].iter().map(|s| s.seq.as_slice()).collect();
    // mean_kimura returns percent (0–100); Linup prints the 0–1 value.
    let avg_k    = kimura::mean_kimura(consensus, &instances, false) / 100.0;
    let avg_kadj = kimura::mean_kimura(consensus, &instances, true)  / 100.0;
    writeln!(out, "Avg Kimura Div: {:.2}", avg_k)?;
    writeln!(out, "Avg Kimura Div (CpG adjusted): {:.2}", avg_kadj)?;

    // If the reference sequence is in old Dfam occupancy format (only x/X/-/./space),
    // print warnings matching Perl's Linup.  Otherwise print a plain INFO line.
    let is_occupancy_rf = ref_row.seq.iter().all(|&b| {
        matches!(b, b'x' | b'X' | b'-' | b'.' | b' ')
    });
    if is_occupancy_rf {
        writeln!(out, "WARNING: RF is in the old Dfam occupancy format, needs updating.")?;
        // Convert consensus to occupancy: bases → 'x', gaps/dots/spaces unchanged.
        let cons_occ: Vec<u8> = consensus.iter().map(|&b| {
            if b == b'-' || b == b'.' || b == b' ' { b } else { b'x' }
        }).collect();
        if ref_row.seq != cons_occ.as_slice() {
            write!(out, "WARNING: RF differs from consensus in occupancy. First difference at:")?;
            let mut bp_pos: u64 = 0;
            for (&r, &c) in ref_row.seq.iter().zip(cons_occ.iter()) {
                if c == b'x' { bp_pos += 1; }
                if r != c { write!(out, "{}", bp_pos)?; break; }
            }
            writeln!(out)?;
        }
    } else if ref_row.seq != consensus {
        writeln!(out, "INFO: Consensus differs from reference sequence.")?;
    }

    let cons_seq: Vec<u8> = consensus.iter().copied().filter(|&b| b != b'-').collect();
    writeln!(out, "Cons length: {}", cons_seq.len())?;

    // FASTA block: Perl emits "\n\n>name\nseq\n\n".
    write!(out, "\n\n>{}\n", ref_row.name)?;
    out.write_all(&cons_seq)?;
    write!(out, "\n\n")?;

    Ok(())
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Build the displayed block bytes for one instance over `[line_start, line_end]`.
///
/// Positions outside the alphabetic region `[bs, be]` are returned as spaces,
/// matching Perl's treatment of leading/trailing gap characters as padding.
fn build_block(
    seq: &[u8],
    line_start: usize,
    line_end: usize,
    base_start: usize,
    base_end: usize,
) -> Vec<u8> {
    let block_len = line_end - line_start + 1;
    let mut out = vec![b' '; block_len];

    let content_start = base_start.max(line_start);
    let content_end   = base_end.min(line_end);
    if content_start <= content_end && content_start < seq.len() {
        let src_end = content_end.min(seq.len() - 1);
        let dst_off = content_start - line_start;
        let src_len = src_end - content_start + 1;
        out[dst_off..dst_off + src_len]
            .copy_from_slice(&seq[content_start..=src_end]);
    }
    out
}

/// Write one alignment row in Linup format.
fn write_row(
    out: &mut dyn Write,
    name: &str,
    max_id_len: usize,
    start: u64,
    max_coord_len: usize,
    seq: &[u8],
    block_size: usize,
    end: u64,
    index: Option<usize>,
) -> io::Result<()> {
    let start_s = start.to_string();
    let end_s   = end.to_string();

    write!(out, "{}{} ", name, " ".repeat(max_id_len.saturating_sub(name.len())))?;
    write!(out, "{}{} ", " ".repeat(max_coord_len.saturating_sub(start_s.len())), start_s)?;
    out.write_all(seq)?;
    if seq.len() < block_size {
        write!(out, "{}", " ".repeat(block_size - seq.len()))?;
    }
    write!(out, "    {}", end_s)?;
    if let Some(idx) = index {
        write!(out, " [{}]", idx)?;
    }
    writeln!(out)
}

/// Diff character between consensus `c` and reference `r`.
fn diff_char(c: u8, r: u8) -> u8 {
    let cu = c.to_ascii_uppercase();
    let ru = r.to_ascii_uppercase();
    if cu == ru { return b' '; }
    if c == b'-' || r == b'-' { return b'-'; }
    match (cu, ru) {
        (b'C', b'T') | (b'T', b'C') | (b'A', b'G') | (b'G', b'A') => b'i',
        (b'G', b'T') | (b'T', b'G') | (b'G', b'C') | (b'C', b'G')
        | (b'C', b'A') | (b'A', b'C') | (b'A', b'T') | (b'T', b'A') => b'v',
        _ => b'?',
    }
}

/// Count alphabetic bytes (both cases) — for consensus/reference position tracking.
fn count_alpha(bytes: &[u8]) -> usize {
    bytes.iter().filter(|&&b| b.is_ascii_alphabetic()).count()
}

/// Count uppercase alphabetic bytes — for instance position tracking (Perl: tr/A-Z/).
fn count_upper(bytes: &[u8]) -> usize {
    bytes.iter().filter(|&&b| b.is_ascii_uppercase()).count()
}
