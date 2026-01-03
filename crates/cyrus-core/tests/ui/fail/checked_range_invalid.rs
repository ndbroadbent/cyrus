//! Test: checked_range! rejects invalid ranges at compile time.
//!
//! This should fail to compile because start >= end.

use cyrus_core::checked_range;
use cyrus_core::types::range::CheckedRange;

fn main() {
    // ERROR: checked_range! requires start < end
    const BAD: CheckedRange<usize> = checked_range!(10, 5);
    let _ = BAD;
}
