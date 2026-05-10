//! LWE ciphertext type and the LWE→RLWE embedding.
//!
//! The embedding is **bit-identical** to the one used in CDKS \[18\] —
//! we keep this primitive verbatim. See SPEC.md §9.f for the rationale.
//!
//! Phase 4 status: type declarations only. The embedding bodies
//! (`a_tilde`, `b_tilde`) live in Phase 5.

use spiral_rs::poly::PolyMatrixRaw;

use crate::error::InspiringError;
use crate::params::RlweParams;

/// A single LWE ciphertext under an LWE secret `s ∈ Z_q^d`:
///
/// `b = -⟨a, s⟩ + e + Δ·m  mod q`.
///
/// `a` is `d` coefficients in `Z_q`, `b` is one coefficient in `Z_q`.
/// SPEC.md §1 / §4.
#[derive(Debug, Clone)]
pub struct LweCiphertext {
    /// LWE random component, length `d`.
    pub a: Vec<u64>,
    /// LWE pseudorandom component (single scalar in `Z_q`).
    pub b: u64,
}

/// A batch of `d` LWE ciphertexts to be packed by [`crate::pack::pack`] into
/// one RLWE ciphertext.
///
/// Per SPEC.md §8 (offline/online split), only the `b_k` scalars are needed
/// in the **online** call to [`crate::pack::pack`]; the `a_k` vectors are
/// consumed during preprocessing. We still store both here so a single
/// `LweBatch` value is round-trippable through the all-online execution
/// (used in `tests/offline_online_equivalence.rs`).
#[derive(Debug, Clone)]
pub struct LweBatch {
    /// All `d` LWE ciphertexts. `inner.len() == params.d`.
    pub inner: Vec<LweCiphertext>,
}

impl LweBatch {
    /// Validate the batch shape against the parameter set.
    ///
    /// Phase 4 status: stub.
    pub fn validate(&self, _params: &RlweParams) -> Result<(), InspiringError> {
        Err(InspiringError::Internal(
            "LweBatch::validate not yet implemented (Phase 5)",
        ))
    }
}

/// `ã := Σ_{i=0}^{d-1} a[i] · X^{-i}` (paper Equation 1, the LWE-to-RLWE
/// embedding for the random component). Matches CDKS \[18\] verbatim.
///
/// Phase 4 status: stub.
pub fn a_tilde<'a>(_params: &'a RlweParams, _a: &[u64]) -> PolyMatrixRaw<'a> {
    unimplemented!("lwe::a_tilde is implemented in Phase 5")
}

/// `b̃ := b · X^0` — a constant polynomial. Trivial half of Equation 1.
///
/// Phase 4 status: stub.
pub fn b_tilde<'a>(_params: &'a RlweParams, _b: u64) -> PolyMatrixRaw<'a> {
    unimplemented!("lwe::b_tilde is implemented in Phase 5")
}
