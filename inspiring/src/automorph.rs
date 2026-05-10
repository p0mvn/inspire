//! Galois automorphisms `τ_g`, `τ_h`, and iterated `τ_g^j`.
//!
//! See SPEC.md §2 (Galois group) and §3 (Lemma 1, the trace operator).
//!
//! - `τ_g(p)(X) = p(X^5)` generates the `Z_{d/2}` factor of `Gal(R)`.
//! - `τ_h(p)(X) = p(X^{2d-1})` generates the `Z_2` factor.
//!
//! Both are realised by [`spiral_rs::poly::automorph_alloc`] which is
//! generic in the exponent. We add helpers for the iterated `τ_g^j`
//! (we cache the precomputed exponents `5^j mod 2d`) and a
//! NTT-form wrapper that round-trips through coefficient form (Phase 11
//! will replace it with an in-place NTT-slot permutation).

use spiral_rs::poly::{PolyMatrixNTT, PolyMatrixRaw};

/// The fixed generator of the `Z_{d/2}` factor of `Gal(R)`, per SPEC.md §2.
pub const G: u64 = 5;

/// `h = 2d − 1`, the generator of the `Z_2` factor of `Gal(R)`,
/// per SPEC.md §2.
#[must_use]
pub const fn h(d: usize) -> u64 {
    (2 * d as u64) - 1
}

/// `5^j mod 2d`, the exponent passed to `spiral_rs::poly::automorph` to
/// realise `τ_g^j`. SPEC.md §2.
///
/// Phase 4 status: stub.
#[must_use]
pub fn tau_g_pow(_j: usize, _d: usize) -> u64 {
    unimplemented!("automorph::tau_g_pow is implemented in Phase 5")
}

/// In-place application of `τ_t` to a coefficient-form polynomial matrix.
/// Trivial passthrough to [`spiral_rs::poly::automorph`]; declared here so
/// callers don't need to import spiral-rs directly.
///
/// Phase 4 status: stub.
pub fn tau_raw<'a>(_a: &PolyMatrixRaw<'a>, _t: u64) -> PolyMatrixRaw<'a> {
    unimplemented!("automorph::tau_raw is implemented in Phase 5")
}

/// `τ_t` for an NTT-form polynomial matrix. Phase 5 implements this as a
/// round-trip through coefficient form (see `docs/spiral-rs-mapping.md`
/// §3). Phase 11 hardening replaces the body with an in-place NTT-slot
/// permutation; the public signature is stable.
///
/// Phase 4 status: stub.
pub fn tau_ntt<'a>(_a: &PolyMatrixNTT<'a>, _t: u64) -> PolyMatrixNTT<'a> {
    unimplemented!("automorph::tau_ntt is implemented in Phase 5")
}

/// Lemma 1's trace `Tr(p) = Σ_{j=0}^{d/2-1} τ_g^j(p) + τ_h ∘ τ_g^j(p)`
/// (SPEC.md §3). Used by `tests/lemma1_trace.rs`.
///
/// Phase 4 status: stub. Delegates to `tau_raw`/`tau_ntt`/spiral-rs add.
pub fn trace<'a>(_p: &PolyMatrixRaw<'a>) -> PolyMatrixRaw<'a> {
    unimplemented!("automorph::trace is implemented in Phase 5")
}
