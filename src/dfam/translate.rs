/// Translation of nucleotide sequence under the standard genetic code (NCBI table 1),
/// used to check an `#=GF FT` `/translation` against the consensus.
///
/// IUB ambiguity codes are expanded, and a codon yields a residue whenever all of its
/// expansions agree (`GTN` -> `V`, `AGR` -> `R`), otherwise `X`.  This is what Repbase
/// does, and it matters: without it only 80.8% of RepBase31.06 translations reproduce,
/// against 99.0% with it.

/// Standard code residues, indexed by codon with bases ordered `T`, `C`, `A`, `G`.
const STANDARD_CODE: &[u8; 64] =
    b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/// Unambiguous base indices (into `T`, `C`, `A`, `G`) that an IUB code stands for.
/// Anything unrecognised, including the `X` allowed in `#=GC RF`, is treated as `N`.
fn expansions(base: u8) -> &'static [usize] {
    match base.to_ascii_uppercase() {
        b'T' | b'U' => &[0],
        b'C' => &[1],
        b'A' => &[2],
        b'G' => &[3],
        b'Y' => &[0, 1],
        b'R' => &[2, 3],
        b'W' => &[0, 2],
        b'S' => &[1, 3],
        b'K' => &[0, 3],
        b'M' => &[1, 2],
        b'B' => &[0, 1, 3],
        b'D' => &[0, 2, 3],
        b'H' => &[0, 1, 2],
        b'V' => &[1, 2, 3],
        _ => &[0, 1, 2, 3],
    }
}

/// Translate one codon; `X` when its expansions do not all give the same residue.
pub fn translate_codon(codon: &[u8]) -> u8 {
    let mut residue = None;
    for &a in expansions(codon[0]) {
        for &b in expansions(codon[1]) {
            for &c in expansions(codon[2]) {
                let r = STANDARD_CODE[16 * a + 4 * b + c];
                match residue {
                    None => residue = Some(r),
                    Some(seen) if seen != r => return b'X',
                    Some(_) => {}
                }
            }
        }
    }
    residue.unwrap_or(b'X')
}

/// Translate whole codons; a trailing partial codon is ignored.
pub fn translate(nt: &[u8]) -> Vec<u8> {
    nt.chunks_exact(3).map(translate_codon).collect()
}

/// IUB complement.  Unrecognised characters are returned unchanged, and then read as
/// `N` by `translate_codon`.
pub fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other, // S, W and N are their own complements
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_codons() {
        assert_eq!(translate(b"ATGTGGTAATAGTGA"), b"MW***");
        assert_eq!(translate(b"atgtgg"), b"MW");
        assert_eq!(translate(b"AUG"), b"M");
    }

    #[test]
    fn ambiguity_resolves_when_all_expansions_agree() {
        assert_eq!(translate(b"GTNACNAGR"), b"VTR");
        // AAN is N or K; NNN is anything.
        assert_eq!(translate(b"AANNNN"), b"XX");
    }

    #[test]
    fn partial_trailing_codon_is_ignored() {
        assert_eq!(translate(b"ATGTG"), b"M");
    }

    #[test]
    fn complement_is_an_involution_over_iub() {
        for &b in b"ACGTRYSWKMBDHVN" {
            assert_eq!(complement(complement(b)), b, "{}", b as char);
        }
    }
}
