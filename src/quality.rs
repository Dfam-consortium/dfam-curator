/// Low-quality alignment column detection using the Ruzzo-Tompa algorithm.
///
/// This mirrors the Perl `getLowScoringAlignmentColumns()` /
/// `_ruzzoTompaFindAllMaximalScoringSubsequences()` in MultAln.pm.
///
/// Workflow:
/// 1. Compute a per-column quality score using the substitution matrix.
/// 2. Negate the scores so low-quality regions become high-scoring targets.
/// 3. Apply Ruzzo-Tompa to find all maximal-scoring contiguous windows.
/// 4. Return the resulting ranges as low-quality blocks.
use crate::alignment::MultiAlign;
use crate::matrix::{ALPHA_LEN, MATRIX};

/// A contiguous range of alignment columns identified as low-quality.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LowQualityBlock {
    /// First column of the block (0-based, inclusive).
    pub start: usize,
    /// Last column of the block (0-based, inclusive).
    pub end: usize,
    /// Cumulative negated-quality score of the block (higher = worse quality).
    pub score: i64,
}

/// Compute per-column alignment quality scores.
///
/// For each column the quality score is the best-fitting IUPAC consensus score,
/// summed over all observed bases in that column.  Higher = better aligned.
pub fn column_scores(msa: &MultiAlign) -> Vec<i64> {
    msa.build_profile(false)
        .iter()
        .map(|col| best_column_score(col))
        .collect()
}

fn best_column_score(col: &[u32; ALPHA_LEN]) -> i64 {
    let mut best: i64 = i64::MIN;
    for cand in 0..ALPHA_LEN {
        let mut score: i64 = 0;
        for obs in 0..ALPHA_LEN {
            let cnt = col[obs] as i64;
            if cnt > 0 {
                score += cnt * MATRIX[cand][obs] as i64;
            }
        }
        if score > best {
            best = score;
        }
    }
    if best == i64::MIN { 0 } else { best }
}

/// Find all low-quality alignment blocks using the Ruzzo-Tompa algorithm.
///
/// The column quality scores are negated before passing to the algorithm so
/// that low-quality (negative-score) regions become the high-scoring targets.
///
/// Only blocks with total negated score > `threshold` are returned.
/// The Perl default threshold is 1.
pub fn find_low_quality_blocks(msa: &MultiAlign, threshold: i64) -> Vec<LowQualityBlock> {
    let scores = column_scores(msa);
    let negated: Vec<i64> = scores.iter().map(|&s| -s).collect();
    ruzzo_tompa_max_segments(&negated, threshold)
}

/// Ruzzo-Tompa algorithm for finding all maximal positive-sum contiguous
/// subsequences.
///
/// Direct port of `_ruzzoTompaFindAllMaximalScoringSubsequences` from
/// MultAln.pm.  Returns only segments with total score > `threshold`.
///
/// A segment [l, r] is maximal if:
/// - Its score is positive (or > threshold)
/// - No proper extension of the segment has a higher score
/// - No proper sub-segment has the same or higher score
///
/// The algorithm runs in O(n) time using a running prefix-sum and a stack of
/// candidate left endpoints.
pub fn ruzzo_tompa_max_segments(scores: &[i64], threshold: i64) -> Vec<LowQualityBlock> {
    let n = scores.len();

    // Per-interval data (mirrors the Perl @I, @L, @R, @Lidx arrays).
    struct Interval {
        start: usize, // left endpoint (inclusive, 0-based) — Perl Lidx[k]
        end: usize,   // right endpoint (exclusive, 0-based) — Perl I[k][1]
        l: i64,       // prefix sum before `start` — Perl L[k]
        r: i64,       // prefix sum at `end` (inclusive) — Perl R[k]
    }

    let mut ivals: Vec<Interval> = Vec::new();
    let mut total: i64 = 0;
    let mut k: usize = 0; // number of current valid intervals

    for i in 0..n {
        total += scores[i];

        // Only consider strictly positive input values (mirrors the Perl `if $b[$i] > 0`).
        if scores[i] <= 0 {
            continue;
        }

        // Initialise a new interval covering just position i.
        if k < ivals.len() {
            ivals[k].start = i;
            ivals[k].end = i + 1;
            ivals[k].l = total - scores[i];
            ivals[k].r = total;
        } else {
            ivals.push(Interval {
                start: i,
                end: i + 1,
                l: total - scores[i],
                r: total,
            });
        }

        // Merge with preceding intervals as long as the current right-endpoint
        // prefix sum is higher and there is a prior interval with lower left-
        // endpoint prefix sum.
        loop {
            // Scan right-to-left for the largest j < k with L[j] < L[k].
            let maxj = (0..k).rev().find(|&j| ivals[j].l < ivals[k].l);

            match maxj {
                Some(j) if ivals[j].r < ivals[k].r => {
                    // Merge: extend interval j to include position i.
                    ivals[j].end = i + 1;
                    ivals[j].r = total;
                    k = j; // roll back and check if j can merge further
                }
                _ => {
                    k += 1;
                    break;
                }
            }
        }
    }

    // Collect valid intervals (indices 0..k) that exceed the threshold.
    ivals[..k]
        .iter()
        .filter(|iv| iv.r - iv.l > threshold)
        .map(|iv| LowQualityBlock {
            start: iv.start,
            end: iv.end - 1, // convert half-open to inclusive
            score: iv.r - iv.l,
        })
        .collect()
}

/// Annotate a `MultiAlign` with its low-quality blocks (stored in
/// `msa.low_quality_blocks` for later use).
pub fn annotate_low_quality(msa: &mut MultiAlign, threshold: i64) {
    let blocks = find_low_quality_blocks(msa, threshold);
    msa.low_quality_blocks = blocks.into_iter().map(|b| (b.start, b.end)).collect();
}

#[cfg(test)]
mod tests {
    use super::*;

    /// All-positive input (columns with negated scores all negative) — no blocks.
    #[test]
    fn all_good_no_blocks() {
        let negated = vec![-100i64, -100, -100, -100];
        let blocks = ruzzo_tompa_max_segments(&negated, 0);
        assert!(blocks.is_empty());
    }

    /// One contiguous bad region flanked by good columns.
    #[test]
    fn one_bad_block() {
        // quality scores: [10, -5, -5, -5, 10], negated: [-10, 5, 5, 5, -10]
        let negated = vec![-10i64, 5, 5, 5, -10];
        let blocks = ruzzo_tompa_max_segments(&negated, 0);
        assert_eq!(blocks.len(), 1, "expected 1 block, got {:?}", blocks);
        assert_eq!(blocks[0].start, 1);
        assert_eq!(blocks[0].end, 3);
        assert_eq!(blocks[0].score, 15);
    }

    /// Two separated bad regions split by a very good column.
    #[test]
    fn two_separated_bad_blocks() {
        // quality scores: [-5, -5, 100, -5, -5], negated: [5, 5, -100, 5, 5]
        let negated = vec![5i64, 5, -100, 5, 5];
        let blocks = ruzzo_tompa_max_segments(&negated, 0);
        assert_eq!(blocks.len(), 2, "expected 2 blocks, got {:?}", blocks);
        assert_eq!(blocks[0].start, 0);
        assert_eq!(blocks[0].end, 1);
        assert_eq!(blocks[1].start, 3);
        assert_eq!(blocks[1].end, 4);
    }

    /// Perl documentation example: scores 4,-5,3,-3,1,2,-2,2,-2,1,5
    /// Expected intervals (0-based half-open): [0,1], [2,3], [4,11]
    #[test]
    fn perl_doc_example() {
        let scores = vec![4i64, -5, 3, -3, 1, 2, -2, 2, -2, 1, 5];
        let blocks = ruzzo_tompa_max_segments(&scores, 0);
        assert_eq!(blocks.len(), 3, "expected 3 blocks, got {:?}", blocks);
        assert_eq!(blocks[0].start, 0); assert_eq!(blocks[0].end, 0); // [0,1) → [0,0]
        assert_eq!(blocks[1].start, 2); assert_eq!(blocks[1].end, 2); // [2,3) → [2,2]
        assert_eq!(blocks[2].start, 4); assert_eq!(blocks[2].end, 10); // [4,11) → [4,10]
        assert_eq!(blocks[2].score, 7); // 1+2-2+2-2+1+5 = 7
    }
}
