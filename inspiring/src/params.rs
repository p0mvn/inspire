//! Parameter sets for `inspiring`.
//!
//! See [SPEC.md §1](../../SPEC.md) for the symbol table this module mirrors.
//! Every public field is paper-named, every constraint is paper-justified.
//!
//! Phase 4 status: declarations and validators only. Phase 5+ phases attach
//! NTT tables and convert this struct into a [`spiral_rs::params::Params`]
//! for the underlying ring arithmetic.

use spiral_rs::params::Params as SpiralParams;

use crate::error::InspiringError;

/// Gadget-decomposition parameters `(z, ℓ)` used by `KS.Setup` /
/// `KS.Switch`.
///
/// Per SPEC.md §1:
///
/// > `g_z = [1, z, z², …, z^{ℓ-1}]^⊤ ∈ Z_q^ℓ`,
/// > `ℓ = ⌈log q / log z⌉`,
/// > `g_z^{-1}: Z_q → Z^{1×ℓ}` returns digit decomposition with each digit
/// > in `[-z/2, z/2)`.
///
/// We require `z = 2^bits_per` (a power of two) so that the underlying
/// `spiral-rs` gadget — whose digit width is fixed at *bit-decomposition*
/// granularity — agrees with InspiRING's specification. See
/// [`docs/spiral-rs-mapping.md` §2](../../docs/spiral-rs-mapping.md).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GadgetParams {
    /// `bits_per = log₂(z)`. Must equal what
    /// [`spiral_rs::gadget::get_bits_per`] would return for `(modulus, ℓ)`.
    pub bits_per: u32,
    /// Number of digits `ℓ`. Must satisfy `z^ℓ ≥ q`.
    pub ell: usize,
}

impl GadgetParams {
    /// `z = 2^bits_per`.
    #[must_use]
    pub const fn z(self) -> u64 {
        1u64 << self.bits_per
    }
}

/// Public parameters for the `inspiring` ring-packing scheme.
///
/// Naming follows SPEC.md §1 exactly. Where this differs from
/// [`spiral_rs::params::Params`] (which is tailored to Spiral-PIR), the
/// helper [`RlweParams::to_spiral_params`] (added in Phase 5) does the
/// translation.
#[derive(Debug, Clone)]
pub struct RlweParams {
    /// Ring degree `d`. Power of two. Both the LWE dim and the RLWE degree.
    pub d: usize,
    /// Modulus `q`. Must be **odd** (so that `d^{-1} mod q` exists; see
    /// SPEC.md §1 and Lemma 1 in §3).
    pub q: u64,
    /// Plaintext modulus `p`. Messages live in `Z_p`.
    pub p: u64,
    /// Subgaussian parameter `σ_χ` of the noise distribution.
    pub sigma_chi: f64,
    /// Gadget-decomposition `(z, ℓ)`.
    pub gadget: GadgetParams,
    /// Cached `Δ = ⌊q / p⌋`.
    pub delta: u64,
    /// Cached `d^{-1} mod q`.
    pub d_inv: u64,
}

impl RlweParams {
    /// Construct an [`RlweParams`] and validate the constraints from
    /// [SPEC.md §1](../../SPEC.md):
    ///
    /// - `d` is a power of two.
    /// - `q` is odd (required for `d^{-1} mod q`).
    /// - `p ≥ 2` and `p ≤ q`.
    /// - `gadget.z()^ell ≥ q` (gadget covers the modulus).
    /// - `gadget.bits_per` matches what `spiral_rs::gadget::get_bits_per`
    ///   would return — see `docs/spiral-rs-mapping.md` §2.
    ///
    /// Phase 4 status: declarations only. Body is a stub; Phase 5
    /// implements the validator and computes `delta`, `d_inv`.
    pub fn new(
        _d: usize,
        _q: u64,
        _p: u64,
        _sigma_chi: f64,
        _gadget: GadgetParams,
    ) -> Result<Self, InspiringError> {
        Err(InspiringError::Internal(
            "RlweParams::new not yet implemented (Phase 5)",
        ))
    }

    /// Convert to a [`spiral_rs::params::Params`] suitable for passing into
    /// the underlying NTT/poly stack. Spiral-PIR-specific fields
    /// (`t_conv`, `t_gsw`, `db_dim_*`, …) are filled with safe no-op
    /// defaults documented in [`docs/spiral-rs-mapping.md` §4](../../docs/spiral-rs-mapping.md).
    ///
    /// Phase 4 status: stub.
    #[must_use = "the spiral-rs Params must be plumbed into PolyMatrix* allocators"]
    pub fn to_spiral_params(&self) -> SpiralParams {
        unimplemented!("RlweParams::to_spiral_params is implemented in Phase 5")
    }
}
