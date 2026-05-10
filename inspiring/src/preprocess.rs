//! `PackPreprocessed`: the CRS-model offline cache.
//!
//! See SPEC.md §8 (offline / online split). Every quantity in Algorithm 1
//! that depends only on `(A, K_g, K_h)` (and not on the LWE `b` scalars)
//! is materialised here, in NTT form, so the online [`crate::pack::pack`]
//! call is a pure function of `(b_0, …, b_{d-1}, &PackPreprocessed)`.
//!
//! Phase 4 status: type and signature declarations only.

use spiral_rs::poly::PolyMatrixNTT;

use crate::error::InspiringError;
use crate::key_switching::KeySwitchingMatrix;
use crate::params::RlweParams;

/// All preprocessable data for a single CRS `A` and a single pair of
/// key-switching matrices `(K_g, K_h)`.
///
/// **API invariant (SPEC.md §10)**: this struct holds **exactly two**
/// key-switching matrices. Any reviewer asked to add a third should
/// instead read SPEC.md §9.h and the test `tests/inspiring_vs_cdks_recursion.rs`.
///
/// Phase 4 status: declarations only. Phase 8 fills in the body.
pub struct PackPreprocessed<'a> {
    /// Underlying parameter set.
    pub params: &'a RlweParams,

    /// Per-LWE-slot Stage-1 result: `a_hat[k][j]` is the `j`-th
    /// component of `IRCtx`'s `â` for input slot `k`. `a_hat.len() == d`,
    /// `a_hat[k].len() == d`. SPEC.md §4.
    ///
    /// All NTT-form. CRS-side, fully preprocessable.
    pub a_hat: Vec<Vec<PolyMatrixNTT<'a>>>,

    /// Stage-2 aggregated `â_agg = Σ_k a_hat[k] · X^k`. SPEC.md §5.
    pub a_agg: Vec<PolyMatrixNTT<'a>>,

    /// `K_g`: the base key-switching matrix for the `τ_g`-cycle.
    pub kg: KeySwitchingMatrix<'a>,

    /// `K_h`: the final-step key-switching matrix that folds the
    /// `τ_h(s̃)` share into `s̃`.
    pub kh: KeySwitchingMatrix<'a>,

    /// Cache of `τ_g^{k-1}(K_g)` for `k ∈ [1, d/2)`, plus
    /// `τ_h(τ_g^{k-1}(K_g))` for the second half. Computed once per
    /// CRS so the online path never invokes an automorphism on `K_g`.
    /// SPEC.md §6.
    pub kg_images_left: Vec<KeySwitchingMatrix<'a>>,
    /// Same as `kg_images_left` but pre-composed with `τ_h` for the
    /// right-half collapse.
    pub kg_images_right: Vec<KeySwitchingMatrix<'a>>,
}

impl<'a> PackPreprocessed<'a> {
    /// Build all CRS-side data from `(A, K_g, K_h)`. Online callers then
    /// call [`crate::pack::pack`] with just the `b_k` scalars.
    ///
    /// API invariant: this signature accepts exactly two key-switching
    /// matrices. Adding a third is a breaking change and a CDKS-drift
    /// red flag (SPEC.md §9.h).
    ///
    /// Phase 4 status: stub.
    pub fn build(
        _params: &'a RlweParams,
        _crs: &PolyMatrixNTT<'a>,
        _kg: KeySwitchingMatrix<'a>,
        _kh: KeySwitchingMatrix<'a>,
    ) -> Result<Self, InspiringError> {
        Err(InspiringError::Internal(
            "PackPreprocessed::build not yet implemented (Phase 8)",
        ))
    }
}
