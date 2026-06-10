/// Random-access reader for UCSC 2bit genome files.
///
/// Builds a lightweight in-memory index (sequence names → file offsets and
/// N-block tables) on open, then extracts individual base ranges on demand
/// via `pread`-style I/O.  The packed DNA data is never fully decoded;
/// only the bytes covering the requested range are read from disk.
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::os::unix::fs::FileExt; // read_exact_at — Linux/macOS only
use std::path::Path;

struct SeqEntry {
    dna_size: u64,
    n_starts: Vec<u32>,
    n_sizes: Vec<u32>,
    /// Absolute byte offset in the file where the packed 2-bit DNA begins.
    packed_dna_offset: u64,
}

/// Random-access reader for a 2bit file.
///
/// `fetch` uses `pread` so multiple calls can be issued without seeking;
/// this also means the struct is `Send + Sync` and `fetch` takes `&self`.
pub struct TwoBitReader {
    file: File,
    seqs: HashMap<String, SeqEntry>,
}

impl TwoBitReader {
    /// Open a 2bit file and build the index.  Does not read any packed DNA.
    pub fn open(path: &Path) -> io::Result<Self> {
        let mut f = File::open(path)?;

        let mut buf4 = [0u8; 4];
        f.read_exact(&mut buf4)?;
        let is_le = match u32::from_be_bytes(buf4) {
            0x1A412743 => false,
            0x4327411A => true,
            _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid 2bit signature")),
        };

        macro_rules! ru32 {
            () => {{
                let mut b = [0u8; 4];
                f.read_exact(&mut b)?;
                if is_le { u32::from_le_bytes(b) } else { u32::from_be_bytes(b) }
            }};
        }
        macro_rules! ru64 {
            () => {{
                let mut b = [0u8; 8];
                f.read_exact(&mut b)?;
                if is_le { u64::from_le_bytes(b) } else { u64::from_be_bytes(b) }
            }};
        }

        let version = ru32!();
        if version > 1 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "unsupported 2bit version"));
        }
        let seq_count = ru32!() as usize;
        ru32!(); // reserved

        // Read the index table: (name, byte offset to sequence record).
        let mut index: Vec<(String, u64)> = Vec::with_capacity(seq_count);
        for _ in 0..seq_count {
            let mut len_buf = [0u8; 1];
            f.read_exact(&mut len_buf)?;
            let mut name_buf = vec![0u8; len_buf[0] as usize];
            f.read_exact(&mut name_buf)?;
            let name = String::from_utf8(name_buf)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let offset = if version == 0 { ru32!() as u64 } else { ru64!() };
            index.push((name, offset));
        }

        // Sort by file offset so we read the sequence headers sequentially.
        index.sort_unstable_by_key(|(_, off)| *off);

        // Parse each sequence record header to capture N-blocks and the
        // absolute offset of the packed DNA payload.
        let mut seqs = HashMap::with_capacity(seq_count);
        for (name, seq_off) in index {
            f.seek(SeekFrom::Start(seq_off))?;

            let dna_size = ru32!() as u64;
            let n_count  = ru32!() as usize;

            let mut n_starts = Vec::with_capacity(n_count);
            let mut n_sizes  = Vec::with_capacity(n_count);
            for _ in 0..n_count { n_starts.push(ru32!()); }
            for _ in 0..n_count { n_sizes.push(ru32!()); }

            let mask_count = ru32!() as usize;
            // Skip mask-block starts, mask-block sizes, and the reserved word.
            f.seek(SeekFrom::Current((mask_count as i64) * 8 + 4))?;

            let packed_dna_offset = f.stream_position()?;
            seqs.insert(name, SeqEntry { dna_size, n_starts, n_sizes, packed_dna_offset });
        }

        Ok(TwoBitReader { file: f, seqs })
    }

    /// Fetch bases in the 0-based half-open range `[start, end)` from `chrom`.
    ///
    /// Returns uppercase ASCII `ACGTN`.  Only the bytes that cover the
    /// requested positions are read from disk.
    pub fn fetch(&self, chrom: &str, start: u64, end: u64) -> io::Result<Vec<u8>> {
        let entry = self.seqs.get(chrom).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("chromosome {:?} not found in 2bit file", chrom),
            )
        })?;

        if end > entry.dna_size {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "range {}..{} exceeds {} length {}",
                    start, end, chrom, entry.dna_size
                ),
            ));
        }

        let length = (end - start) as usize;

        // 2-bit packing: T=0, C=1, A=2, G=3, MSB first within each byte.
        const BASES: [u8; 4] = [b'T', b'C', b'A', b'G'];

        let first_byte = start / 4;
        let last_byte  = (end - 1) / 4;
        let byte_count = (last_byte - first_byte + 1) as usize;

        let mut packed = vec![0u8; byte_count];
        self.file
            .read_exact_at(&mut packed, entry.packed_dna_offset + first_byte)?;

        let mut result = Vec::with_capacity(length);
        for pos in start..end {
            let byte_idx = ((pos / 4) - first_byte) as usize;
            let shift    = (3 - (pos % 4)) * 2;
            result.push(BASES[((packed[byte_idx] >> shift) & 0x03) as usize]);
        }

        // Overlay N-blocks (genomic gaps / unknown bases).
        for (&ns, &nz) in entry.n_starts.iter().zip(entry.n_sizes.iter()) {
            let ns = ns as u64;
            let ne = ns + nz as u64;
            let ov_s = ns.max(start);
            let ov_e = ne.min(end);
            if ov_s < ov_e {
                for i in ov_s..ov_e {
                    result[(i - start) as usize] = b'N';
                }
            }
        }

        Ok(result)
    }

    /// Returns `true` if the named chromosome is present.
    pub fn contains(&self, chrom: &str) -> bool {
        self.seqs.contains_key(chrom)
    }

    /// Returns the total length of the named chromosome, or `None` if absent.
    pub fn seq_len(&self, chrom: &str) -> Option<u64> {
        self.seqs.get(chrom).map(|e| e.dna_size)
    }
}
