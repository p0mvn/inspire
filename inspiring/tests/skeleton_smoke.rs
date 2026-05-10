//! Phase 4 placeholder. Confirms the public API surface compiles; every
//! assertion here is a no-op against the stubs that Phase 5+ will replace.
//!
//! The full Phase 9 test suite (`tests/lemma1_trace.rs`,
//! `tests/transform_correctness.rs`, …, `tests/inspiring_vs_cdks_recursion.rs`)
//! supersedes this file.

use inspiring::{GadgetParams, InspiringError, RlweParams};

#[test]
fn public_api_surface_is_wired_up() {
    let err = RlweParams::new(
        8,
        12289,
        4,
        3.2,
        GadgetParams {
            bits_per: 3,
            ell: 5,
        },
    )
    .expect_err("Phase 4 stub must return an error");

    assert!(
        matches!(err, InspiringError::Internal(_)),
        "expected Internal stub variant, got {err:?}",
    );
}
