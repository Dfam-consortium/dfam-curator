//! Call a consensus from an MSA supplied in the `msabuild` `ROW` format.
//!
//! Feeding this the *C++'s* MSA dump isolates the consensus caller from
//! everything upstream: same columns in, so any difference is the caller.
//!
//! ```sh
//! cargo run -p cons-core --release --example callcons -- m.cpp.txt [min]
//! ```

use aln_core::giri;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let path = args.next().ok_or("usage: callcons <msa.txt> [min_non_gap]")?;
    let min: usize = args.next().map(|s| s.parse()).transpose()?.unwrap_or(1);

    let text = std::fs::read_to_string(&path)?;
    let rows: Vec<Vec<u8>> = text
        .lines()
        .filter_map(|l| l.strip_prefix("ROW\t"))
        .map(|rest| {
            // ROW <TAB> name <TAB> sequence
            let seq = rest.split_once('\t').map(|x| x.1).unwrap_or("");
            // GIRI writes padding as '<' / '>'; canonicalise to a space.
            seq.bytes()
                .map(|b| if b == b'<' || b == b'>' { b' ' } else { b })
                .collect()
        })
        .collect();

    if rows.is_empty() {
        return Err("no ROW records found".into());
    }
    let refs: Vec<&[u8]> = rows.iter().map(|r| r.as_slice()).collect();

    let gapped = giri::get_consensus(&refs, min);
    let ungapped = giri::strip_gaps(&gapped);
    eprintln!(
        "rows={} width={} gapped_len={} ungapped_len={}",
        rows.len(),
        rows[0].len(),
        gapped.len(),
        ungapped.len()
    );
    println!("{}", String::from_utf8_lossy(&gapped));
    Ok(())
}
