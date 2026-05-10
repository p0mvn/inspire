//! # `inspiring` — InspiRING.Pack ring-packing crate
//!
//! A standalone Rust implementation of **Algorithm 1 ([`InspiRING.Pack`])
//! from the InsPIRe paper** (eprint 2025/1352, Mahdavi–Patel–Seo–Yeo, 2025).
//! The crate exposes a single primitive:
//!
//! ```ignore
//! pub fn pack(lwe_b: LweBatch, pre: &PackPreprocessed) -> RlweCiphertext
//! ```
//!
//! which compresses `d` LWE ciphertexts (each of LWE dimension `d`) into a
//! single RLWE ciphertext of degree `d`, using exactly **two** key-switching
//! matrices `K_g` and `K_h`. See [`SPEC.md`] for the mathematical contract
//! and [`docs/spiral-rs-mapping.md`] for the spiral-rs primitive audit.
//!
//! ## Locked-in scope
//!
//! - Algorithm 1 only. No `PartialPack`, no PIR layers.
//! - Built on [`spiral-rs`](https://github.com/menonsamir/spiral-rs) pinned to
//!   `rev = 6929441` (matching the reference Google implementation).
//! - Production posture: offline/online split (CRS model), full unit and
//!   integration tests, statistical noise validation against Theorem 2,
//!   benchmarks reproducing paper Table 5, CI, rustdoc.
//!
//! ## Crate map
//!
//! | Module | Concept (paper §) | Phase |
//! |---|---|---|
//! | [`params`] | `RlweParams`, `GadgetParams`, validators | Phase 4 |
//! | [`lwe`] | `LweCiphertext`, batch type, embedding (Eq. 1) | Phase 4 / 5 |
//! | [`automorph`] | `τ_g`, `τ_h`, `τ_g^j` (§2 + Lemma 1) | Phase 4 / 5 |
//! | [`intermediate`] | `IRCtx`, Stage 1 `transform`, Stage 2 `aggregate` | Phase 5 / 6 |
//! | [`collapse`] | `collapse_one`, `collapse_half`, `collapse` (Stage 3) | Phase 7 |
//! | [`key_switching`] | `KS.Setup`, `KS.Switch`, automorphic images | Phase 7 |
//! | [`preprocess`] | `PackPreprocessed` (CRS-model offline cache) | Phase 8 |
//! | [`mod@pack`] | top-level `pack` (Algorithm 1) | Phase 8 |
//! | [`error`] | `InspiringError` | Phase 4 |
//!
//! [`SPEC.md`]: https://github.com/<TBD>/inspiring/blob/main/SPEC.md
//! [`docs/spiral-rs-mapping.md`]: https://github.com/<TBD>/inspiring/blob/main/docs/spiral-rs-mapping.md
//! [`InspiRING.Pack`]: https://eprint.iacr.org/2025/1352
//!
//! ## Public API invariants
//!
//! These are also asserted by tests (`tests/inspiring_vs_cdks_recursion.rs`):
//!
//! 1. [`PackPreprocessed::build`](preprocess::PackPreprocessed::build) accepts
//!    **exactly two** key-switching matrices, `kg` and `kh`.
//! 2. A single call to [`pack::pack`] invokes `KS.Switch` exactly `d − 1` times.
//! 3. [`pack::pack`] is a deterministic function of `(lwe_b, pre)` —
//!    no fresh randomness is sampled in the online path.
//!
//! ## Toolchain & platform
//!
//! - **Nightly Rust** (pinned by `rust-toolchain.toml`): inherited from
//!   `spiral-rs`'s `#![feature(stdarch_x86_avx512)]`.
//! - **AVX-512** target feature: required by `spiral-rs`'s NTT inner loops
//!   *and* — more subtly — by a correctness bug in `spiral-rs`'s scalar
//!   `multiply` fallback (see the comment block above the `compile_error!`
//!   gate in this file, plus `docs/spiral-rs-mapping.md` §1). The crate
//!   refuses to build without `target_feature = "avx512f"` to make this a
//!   compile-time error rather than a silent run-time miscomputation. CI
//!   runs on `x86_64-unknown-linux-gnu`; the crate is not portable to
//!   `aarch64-*` without a spiral-rs port.
//!
//! See [`docs/spiral-rs-mapping.md`] for the full audit of inherited
//! constraints.

#![cfg_attr(docsrs, feature(doc_cfg))]
// Phase 11 hardening will flip these `warn`s to `deny` and re-enable
// `clippy::pedantic`. During Phase 4 (skeleton) and Phase 5–8
// (implementation) we keep them as `warn` so half-built modules don't break
// local rustdoc runs, and we deliberately do NOT opt into pedantic — it
// flags every stub signature here as `#[must_use]`-missing, which is the
// wrong call for skeleton fns that always `unimplemented!()`.
#![warn(missing_docs)]
#![warn(rust_2018_idioms)]

// =============================================================================
//  Compile-time AVX-512 gate (CORRECTNESS, NOT PERFORMANCE)
// =============================================================================
//
// `inspiring` *requires* the `avx512f` target feature at build time. This is
// not a performance optimisation — it is a correctness requirement, because
// the only scalar (non-SIMD) `multiply` path in our pinned `spiral-rs`
// revision is **silently buggy** for our parameter regime.
//
// Pin-point of the upstream bug
// -----------------------------
//   `spiral-rs` rev `6929441c6551769b7d099d3af3df347cde3bae7b`
//   `src/arith.rs:28-33`, function `multiply_add_modular`:
//
//       pub fn multiply_add_modular(params: &Params,
//                                   a: u64, b: u64, x: u64, c: usize) -> u64 {
//           if params.crt_count == 1 {
//               return multiply_uint_mod(a, b, params.moduli[c]);  // BUG
//           }
//           barrett_coeff_u64(params, a * b + x, c)
//       }
//
// The `crt_count == 1` branch returns `a * b mod q` and **drops the
// accumulator `x`**. The function is called from `multiply_add_poly`
// (`src/poly.rs:404`), which is in turn the inner loop of the scalar
// `multiply(res, a, b)` (`src/poly.rs:543`, gated `cfg(not(target_feature =
// "avx2"))`). The AVX2 sibling at `src/poly.rs:566` accumulates products in
// 64-bit lanes and reduces *after* the loop, so it is correct.
//
// Net effect on InspiRING: every `KS.Switch` (which multiplies a `[2, ℓ]`
// key-switching matrix by a `[ℓ, 1]` digit column) silently keeps **only the
// last gadget term** instead of the gadget sum. For our default gadget the
// last digit is the high-order base-`z` digit, which is zero for any
// coefficient `< z^{ℓ-1}` — i.e. every legitimately-bounded ciphertext
// coefficient. So `KS.Switch` returns the zero polynomial, the cascade
// `Collapse` step destroys the plaintext, and the only test that exercises a
// gadget-sum (the full `transform → aggregate → collapse` round-trip) fails
// while every smaller test passes.
//
// `inspiring` runs with `crt_count == 1` (a single `q`, c.f.
// `RlweParams::new` in `params.rs`), so we hit this branch on every
// non-AVX2 build. The fix is to ensure the AVX2 / AVX-512 codegen path is
// always taken, which is what the gate below enforces.
//
// Why `avx512f` and not `avx2`
// ----------------------------
// `spiral-rs`'s NTT inner loops use `_mm512_*` intrinsics gated on
// `cfg(target_feature = "avx2")` (the gate is a misnomer — see
// `docs/spiral-rs-mapping.md` §1). On a host that has AVX2 but not AVX-512,
// the build succeeds but execution traps with `SIGILL`. Gating on
// `avx512f` keeps the compile-time check honest about the runtime hardware
// requirement. `target-cpu=skylake-avx512` (the default in
// `.cargo/config.toml`) and `target-cpu=native` on any AVX-512 host both
// satisfy this gate.
//
// `docsrs` is exempted because docs.rs's builders don't expose AVX-512 and
// rustdoc never executes the code anyway.
#[cfg(not(any(target_feature = "avx512f", docsrs)))]
compile_error!(
    "inspiring requires an AVX-512 build (e.g. RUSTFLAGS='-C target-cpu=skylake-avx512' or \
     '-C target-cpu=native' on an AVX-512 host). This is a correctness requirement, not a \
     performance one: the scalar fallback in spiral-rs rev 6929441 has a bug in \
     `arith::multiply_add_modular` for `crt_count == 1` that silently zeroes out KS.Switch \
     results. See the comment in src/lib.rs above this `compile_error!` and \
     docs/spiral-rs-mapping.md §1 for the full pin-point."
);

pub mod automorph;
pub mod collapse;
pub mod error;
pub mod intermediate;
pub mod key_switching;
pub mod lwe;
pub mod pack;
pub mod params;
pub mod preprocess;

pub use error::InspiringError;
pub use lwe::{LweBatch, LweCiphertext};
pub use pack::{pack, RlweCiphertext};
pub use params::{GadgetParams, RlweParams};
pub use preprocess::PackPreprocessed;
