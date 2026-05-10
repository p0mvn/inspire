//! The intermediate ciphertext `IRCtx = (â, b̃)` and the algorithms that
//! produce it.
//!
//! Two stages:
//!
//! - **Stage 1, [`transform`]** — converts a single LWE ciphertext into an
//!   `IRCtx` whose secret is the structured vector
//!   `ŝ[j] = τ_g^j(s̃)`, `ŝ[j+d/2] = τ_h(τ_g^j(s̃))`. SPEC.md §4
//!   (paper Appendix B).
//! - **Stage 2, [`aggregate`]** — sums `Σ_{k=0}^{d-1} IRCtx(m_k) · X^k`
//!   homomorphically. SPEC.md §5.
//!
//! Phase 4 status: type and signature declarations only.

use spiral_rs::poly::{PolyMatrixNTT, PolyMatrixRaw};

use crate::lwe::LweCiphertext;
use crate::params::RlweParams;

/// Intermediate ciphertext `(â, b̃)` of paper §3.2 / SPEC.md §4.
///
/// - `a_hat` is `d` polynomials in NTT form, satisfying
///   `a_hat[j] = d^{-1} · τ_g^j(ã)` for `j ∈ [0, d/2)` and
///   `a_hat[j + d/2] = d^{-1} · τ_h(τ_g^j(ã))` for `j ∈ [0, d/2)`.
/// - `b_tilde` is one polynomial, the constant polynomial
///   `b · X^0` for Stage 1, and the running `Σ_k b_k · X^k` for Stage 2.
///
/// Note: `Debug` / `Clone` are not derived because [`PolyMatrixNTT`] and
/// [`PolyMatrixRaw`] do not implement them upstream; Phase 5 adds
/// hand-written impls if/when tests need them.
pub struct IRCtx<'a> {
    /// `d` ring elements in NTT form. SPEC.md §4.
    pub a_hat: Vec<PolyMatrixNTT<'a>>,
    /// Single ring element, in **coefficient** form so Stage 2's `X^k`
    /// shifts are cheap. SPEC.md §5.
    pub b_tilde: PolyMatrixRaw<'a>,
}

/// Stage 1 (`TRANSFORM` of paper Algorithm 1; SPEC.md §4):
///
/// `(a, b) ↦ IRCtx` such that decryption under `ŝ` yields
/// the constant polynomial `m̂(X) = m`.
///
/// All output is **CRS-side** (preprocessable) except `b_tilde`, which
/// depends on the LWE `b` scalar. Phase 8's [`crate::preprocess::PackPreprocessed`]
/// caches the `a_hat` half across many packs.
///
/// Phase 4 status: stub.
pub fn transform<'a>(_params: &'a RlweParams, _ct: &LweCiphertext) -> IRCtx<'a> {
    unimplemented!("intermediate::transform is implemented in Phase 5")
}

/// Stage 2 (paper Algorithm 1; SPEC.md §5): compute
///
/// `(â_agg, b̃_agg) = Σ_{k=0}^{d-1} IRCtx(m_k) · X^k`,
///
/// using `X^k` as a coefficient-form monomial shift on `b_tilde` and as a
/// (cached) NTT-form monomial multiply on each `a_hat[j]` slot. The choice
/// of where to absorb the `X^k` factor minimises NTT round-trips; see
/// SPEC.md §5.
///
/// Phase 4 status: stub.
pub fn aggregate<'a>(_params: &'a RlweParams, _per_index: &[IRCtx<'a>]) -> IRCtx<'a> {
    unimplemented!("intermediate::aggregate is implemented in Phase 6")
}
