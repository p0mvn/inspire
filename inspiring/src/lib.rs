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
//! - **AVX-512** target feature: required by `spiral-rs`'s NTT inner loops.
//!   The crate as a whole therefore targets `x86_64-unknown-linux-gnu` (CI)
//!   and is not portable to `aarch64-*` without a spiral-rs port.
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
