//! Confirms the public API surface compiles as implementation phases replace
//! the original skeleton stubs.
//!
//! The full Phase 9 test suite (`tests/lemma1_trace.rs`,
//! `tests/transform_correctness.rs`, …, `tests/inspiring_vs_cdks_recursion.rs`)
//! supersedes this file.

use inspiring::{GadgetParams, RlweParams};

#[test]
fn public_api_surface_is_wired_up() {
    let params = RlweParams::new(
        8,
        12289,
        4,
        3.2,
        GadgetParams {
            bits_per: 3,
            ell: 5,
        },
    )
    .expect("Stage 1 parameters should construct");

    assert_eq!(params.d, 8);
    assert_eq!(params.q, 12289);
    assert_eq!(params.p, 4);
    assert_eq!(params.delta, 3072);
    assert_eq!((params.d as u64 * params.d_inv) % params.q, 1);
}
