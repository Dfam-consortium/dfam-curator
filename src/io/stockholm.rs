//! Stockholm 1.0, re-exported from `dfam-stk-io`.
//!
//! Both the parser and the `MultiAlign` bridge now live in dfam-lib, so
//! RepeatAfterMe and anything else can read Stockholm without depending on
//! dfam-curator. This module stays as a re-export so call sites here are
//! unaffected.

pub use dfam_stk_io::msa::{multialign_from_record, read, read_select, write};
