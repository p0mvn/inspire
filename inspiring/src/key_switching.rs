//! Key-switching primitives `KS.Setup` and `KS.Switch`, plus helpers to
//! compute automorphic images `τ_g^{k-1}(K_g)` of a base matrix locally
//! (without extra key material). See SPEC.md §6 (Stage 3) and §9.b
//! (the structural reason InspiRING needs only two base KS matrices vs.
//! CDKS's `lg d`).
//!
//! Phase 4 status: type and signature declarations only. The Phase 7
//! implementation patterns its body on the inline KS body of
//! `spiral_rs::server::coefficient_expansion` (lines 80–103 of
//! `spiral-rs/src/server.rs` at the pinned revision); we cannot call
//! `coefficient_expansion` directly because it is fused with
//! Spiral-PIR's expansion loop. See `docs/spiral-rs-mapping.md` §3.

use rand_chacha::ChaCha20Rng;
use spiral_rs::params::Params as SpiralParams;
use spiral_rs::poly::PolyMatrixNTT;

use crate::params::RlweParams;

/// A single key-switching matrix `K`. Internally a `[2, ℓ]` `PolyMatrixNTT`
/// (the row-2-by-cols-ℓ shape used by spiral-rs's gadget machinery).
///
/// `K = KS.Setup(s', s)` lets one transform a ciphertext under `s'`
/// (one of `τ_g(s̃)`, `τ_h(s̃)`, …) into one under `s = s̃`. SPEC.md §6.
///
/// Note: `Debug` / `Clone` are not derived because [`PolyMatrixNTT`] does
/// not implement them upstream; Phase 7 adds hand-written impls if needed.
pub struct KeySwitchingMatrix<'a> {
    /// The encrypted gadget-scaled secret. Shape `[2, ℓ]`.
    pub mat: PolyMatrixNTT<'a>,
}

/// `KS.Setup(s_from, s_to)` — encrypt the gadget-scaled `s_from` under
/// `s_to` to produce a key-switching matrix, per SPEC.md §6 / paper §2.
///
/// Phase 4 status: stub. See [`docs/spiral-rs-mapping.md` §3](../../docs/spiral-rs-mapping.md)
/// for the implementation plan.
pub fn ks_setup<'a>(
    _params: &'a RlweParams,
    _spiral: &'a SpiralParams,
    _s_from_ntt: &PolyMatrixNTT<'a>,
    _s_to_ntt: &PolyMatrixNTT<'a>,
    _rng: &mut ChaCha20Rng,
) -> KeySwitchingMatrix<'a> {
    unimplemented!("key_switching::ks_setup is implemented in Phase 7")
}

/// `KS.Switch(K, c)` — apply a key-switching matrix to an RLWE
/// ciphertext `c = (c1, c2)`. Returns a new ciphertext under `s_to`.
///
/// The body mirrors the inline KS pattern in `spiral-rs/src/server.rs`
/// lines 80–103: gadget-invert `c1` (raw), NTT-forward, multiply by
/// `K.mat`, add `(0, c2)`. See SPEC.md §6.
///
/// **Test-only instrumentation**: in `cfg(test)` builds a thread-local
/// counter is incremented on every call. `tests/inspiring_vs_cdks_recursion.rs`
/// asserts the counter equals exactly `d − 1` per call to
/// [`crate::pack::pack`]. Tampering with this is a production-blocker.
///
/// Phase 4 status: stub.
pub fn ks_switch<'a>(
    _k: &KeySwitchingMatrix<'a>,
    _c1: &PolyMatrixNTT<'a>,
    _c2: &PolyMatrixNTT<'a>,
) -> (PolyMatrixNTT<'a>, PolyMatrixNTT<'a>) {
    unimplemented!("key_switching::ks_switch is implemented in Phase 7")
}

/// Compute `τ_g^{k-1}(K_g)` from `K_g` without any extra key material.
/// The image is just `K_g` with `τ_g^{k-1}` applied component-wise to
/// each polynomial of the matrix. SPEC.md §6 / Appendix C.
///
/// Phase 4 status: stub.
#[must_use]
pub fn automorphic_image<'a>(_k: &KeySwitchingMatrix<'a>, _t: u64) -> KeySwitchingMatrix<'a> {
    unimplemented!("key_switching::automorphic_image is implemented in Phase 7")
}

/// Test-only thread-local counter for `KS.Switch` calls. Used by
/// `tests/inspiring_vs_cdks_recursion.rs` to assert the linear-cascade
/// `KS.Switch` count of exactly `d − 1` per pack — the runtime structural
/// guard against accidental CDKS-style implementation drift (SPEC.md §9.h).
#[cfg(test)]
pub mod ks_call_count {
    use std::cell::Cell;

    thread_local! {
        static COUNTER: Cell<u64> = const { Cell::new(0) };
    }

    /// Reset to 0. Call before a measured `pack`.
    pub fn reset() {
        COUNTER.with(|c| c.set(0));
    }

    /// Increment by one. Called from inside `ks_switch`.
    pub fn inc() {
        COUNTER.with(|c| c.set(c.get() + 1));
    }

    /// Read the current count.
    #[must_use]
    pub fn get() -> u64 {
        COUNTER.with(Cell::get)
    }
}
