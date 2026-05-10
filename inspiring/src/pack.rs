//! Top-level `pack` entry point — Algorithm 1 of [eprint 2025/1352].
//!
//! [eprint 2025/1352]: https://eprint.iacr.org/2025/1352
//!
//! See SPEC.md §8 for the offline/online split, §10 for the symbol table,
//! and §9 for the structural comparison with CDKS.
//!
//! Phase 4 status: type alias and signature only.

use spiral_rs::poly::PolyMatrixNTT;

use crate::error::InspiringError;
use crate::lwe::LweBatch;
use crate::preprocess::PackPreprocessed;

/// An RLWE ciphertext under the base secret `s̃`.
///
/// Internally a `[2, 1]` `PolyMatrixNTT` (the spiral-rs convention) wrapped
/// in a newtype so callers can't accidentally mix it up with intermediate
/// pseudo-ciphertexts.
pub struct RlweCiphertext<'a> {
    /// `(c1, c2)` stacked vertically. `inner.rows == 2`, `inner.cols == 1`.
    pub inner: PolyMatrixNTT<'a>,
}

/// `InspiRING.Pack(b, pre) -> RlweCiphertext` — Algorithm 1.
///
/// **Online** entry point. Takes the `d` `b_k` scalars (via [`LweBatch`])
/// and a [`PackPreprocessed`] cache; returns a single RLWE ciphertext under
/// `s̃` that decrypts to `Σ_{k=0}^{d-1} m_k · X^k`.
///
/// API invariants (SPEC.md §10):
///
/// 1. Deterministic: no fresh randomness is sampled in this function.
/// 2. Calls `KS.Switch` exactly `d − 1` times. Asserted by
///    `tests/inspiring_vs_cdks_recursion.rs`.
/// 3. Touches `pre.kg` and `pre.kh` only via their precomputed
///    automorphic images cached in `pre.kg_images_left` and
///    `pre.kg_images_right`.
///
/// Phase 4 status: stub.
pub fn pack<'a>(
    _b: &LweBatch,
    _pre: &'a PackPreprocessed<'a>,
) -> Result<RlweCiphertext<'a>, InspiringError> {
    Err(InspiringError::Internal(
        "pack::pack not yet implemented (Phase 8)",
    ))
}
