//! Compile tests using trybuild.
//!
//! - `tests/ui/pass/*.rs` - must compile successfully
//! - `tests/ui/fail/*.rs` - must fail to compile
//!
//! To update expected error output after intentional changes:
//! ```bash
//! TRYBUILD=overwrite cargo test --test trybuild
//! ```

#[test]
fn ui() {
    let t = trybuild::TestCases::new();
    t.pass("tests/ui/pass/*.rs");
    t.compile_fail("tests/ui/fail/*.rs");
}
