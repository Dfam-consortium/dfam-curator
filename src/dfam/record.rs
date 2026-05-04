/// Re-exports from the shared `stk-io` crate.
pub use dfam_stk_io::{iter_records, SeqRow, StkRecord, StkRecordIter};

/// Backward-compatible alias so existing code referencing `RawDfamRecord` continues to compile.
pub type RawDfamRecord = StkRecord;
