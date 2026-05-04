// Compact index: map ASCII byte -> matrix row/column index (255 = invalid)

static IDX: [u8; 256] = {
    let mut t = [255u8; 256];
    // A=0 R=1 G=2 C=3 Y=4 T=5 K=6 M=7 S=8 W=9 N=10 X=11 Z=12 V=13 H=14 D=15 B=16 -=17
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'R' as usize] = 1;
    t[b'r' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'C' as usize] = 3;
    t[b'c' as usize] = 3;
    t[b'Y' as usize] = 4;
    t[b'y' as usize] = 4;
    t[b'T' as usize] = 5;
    t[b't' as usize] = 5;
    t[b'K' as usize] = 6;
    t[b'k' as usize] = 6;
    t[b'M' as usize] = 7;
    t[b'm' as usize] = 7;
    t[b'S' as usize] = 8;
    t[b's' as usize] = 8;
    t[b'W' as usize] = 9;
    t[b'w' as usize] = 9;
    t[b'N' as usize] = 10;
    t[b'n' as usize] = 10;
    t[b'X' as usize] = 11;
    t[b'x' as usize] = 11;
    t[b'Z' as usize] = 12;
    t[b'z' as usize] = 12;
    t[b'V' as usize] = 13;
    t[b'v' as usize] = 13;
    t[b'H' as usize] = 14;
    t[b'h' as usize] = 14;
    t[b'D' as usize] = 15;
    t[b'd' as usize] = 15;
    t[b'B' as usize] = 16;
    t[b'b' as usize] = 16;
    t[b'-' as usize] = 17;
    t
};

pub const ALPHA_LEN: usize = 18; // 17 IUPAC + gap

/// Lookup a matrix index for a byte (A/R/G/C/Y/T/K/M/S/W/N/X/Z/V/H/D/B/-).
/// Returns None for unrecognised bytes.
#[inline]
pub fn alpha_idx(b: u8) -> Option<usize> {
    let i = IDX[b as usize];
    if i == 255 { None } else { Some(i as usize) }
}

/// Return the canonical uppercase byte for a matrix index.
pub fn alpha_byte(idx: usize) -> u8 {
    const LUT: [u8; ALPHA_LEN] =
        [b'A', b'R', b'G', b'C', b'Y', b'T', b'K', b'M',
         b'S', b'W', b'N', b'X', b'Z', b'V', b'H', b'D', b'B', b'-'];
    LUT[idx]
}

/// The 18×18 scoring matrix (rows = consensus candidate, cols = observed base).
///
/// Derived directly from the Perl buildConsensusFromArray matrix, with the gap
/// row/column added (gap vs. any letter = -6, gap vs. gap = +3).
///
/// Layout: MATRIX[row][col], both indexed by alpha_idx().
#[rustfmt::skip]
pub static MATRIX: [[i32; ALPHA_LEN]; ALPHA_LEN] = [
    //       A    R    G    C    Y    T    K    M    S    W    N    X    Z    V    H    D    B   [-]
    /* A */ [ 9,   0,  -8, -15, -16, -17, -13,  -3, -11,  -4,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* R */ [ 2,   1,   1, -15, -15, -16,  -7,  -6,  -6,  -7,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* G */ [-4,   3,  10, -14, -14, -15,  -2,  -9,  -2,  -9,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* C */ [-15, -14, -14,  10,   3,  -4,  -9,  -2,  -2,  -9,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* Y */ [-16, -15, -15,   1,   1,   2,  -6,  -7,  -6,  -7,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* T */ [-17, -16, -15,  -8,   0,   9,  -3, -13, -11,  -4,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* K */ [-11,  -6,  -2, -11,  -7,  -3,  -2, -11,  -6,  -7,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* M */ [ -3,  -7, -11,  -2,  -6, -11, -11,  -2,  -6,  -7,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* S */ [ -9,  -5,  -2,  -2,  -5,  -9,  -5,  -5,  -2,  -9,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* W */ [ -4,  -8, -11, -11,  -8,  -4,  -8,  -8, -11,  -4,  -2,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* N */ [ -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -1,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* X */ [ -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -7,  -3,  -3,  -3,  -3,  -3,  -6],
    /* Z */ [ -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -6],
    /* V */ [ -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -6],
    /* H */ [ -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -6],
    /* D */ [ -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -6],
    /* B */ [ -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -3,  -6],
    /* - */ [ -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,   3],
];

/// Score of aligning consensus candidate `row_byte` against observed base `col_byte`.
/// Returns 0 for any unrecognised character (treated as N-like).
#[inline]
pub fn score(row_byte: u8, col_byte: u8) -> i32 {
    match (alpha_idx(row_byte), alpha_idx(col_byte)) {
        (Some(r), Some(c)) => MATRIX[r][c],
        _ => 0,
    }
}
