//! Stage 3 of `InspiRING.Pack`: the collapse from `IRCtx` (a wide
//! `(d+1)`-element ciphertext) down to a 2-element RLWE ciphertext under
//! the base secret `s̃`.
//!
//! Three layered subroutines, exactly as in Algorithm 1 / Appendix C
//! (see SPEC.md §6):
//!
//! - [`collapse_one`] — one key-switch step.
//! - [`collapse_half`] — `d/2 - 1` `collapse_one` calls applied to one
//!   half of `a_agg` using automorphic images of `K_g`, optionally
//!   pre-composed with `τ_h`.
//! - [`collapse`] — runs `collapse_half` twice, then a final `KS.Switch`
//!   with `K_h` to fold the `τ_h(s̃)` share into `s̃`.
//!
//! **Linear-cascade invariant** (SPEC.md §6 + §9):
//!
//! `# KS.Switch calls per pack = (d/2 - 1) + (d/2 - 1) + 1 = d - 1`.
//!
//! A CDKS-style implementation would have `(d - 1) · log₂ d` calls.
//! `tests/inspiring_vs_cdks_recursion.rs` asserts this empirically by
//! reading `key_switching::ks_call_count::get()` after each pack.
//!
//! Phase 4 status: signatures only.

use crate::intermediate::IRCtx;
use crate::key_switching::KeySwitchingMatrix;
use crate::pack::RlweCiphertext;
use crate::params::RlweParams;

/// `CollapseOne` — one cascade step. Takes the running collapse state
/// (a `(2 + remaining)`-element pseudo-ciphertext) plus the appropriate
/// automorphic image of the base KS matrix, and produces a state that
/// has one fewer element. SPEC.md §6 / paper Appendix C.
///
/// Phase 4 status: stub.
pub fn collapse_one<'a, 'b>(_state: &mut CollapseState<'a>, _k_image: &KeySwitchingMatrix<'b>) {
    unimplemented!("collapse::collapse_one is implemented in Phase 7")
}

/// `CollapseHalf` — runs `d/2 - 1` `collapse_one` calls over one half
/// (either the `τ_g^j` half or the `τ_h ∘ τ_g^j` half) of `a_agg`.
///
/// SPEC.md §6.
///
/// Phase 4 status: stub.
pub fn collapse_half<'a, 'b>(
    _state: &mut CollapseState<'a>,
    _kg_images: &[KeySwitchingMatrix<'b>],
) {
    unimplemented!("collapse::collapse_half is implemented in Phase 7")
}

/// `Collapse` — full Stage 3. Runs `collapse_half` twice, then a final
/// `KS.Switch` with `K_h`. Output is an RLWE ciphertext under `s̃`.
///
/// SPEC.md §6.
///
/// Phase 4 status: stub.
pub fn collapse<'a, 'b>(
    _params: &'a RlweParams,
    _agg: IRCtx<'a>,
    _kg_images_left: &[KeySwitchingMatrix<'b>],
    _kg_images_right: &[KeySwitchingMatrix<'b>],
    _kh: &KeySwitchingMatrix<'b>,
) -> RlweCiphertext<'a> {
    unimplemented!("collapse::collapse is implemented in Phase 7")
}

/// Running state of the collapse cascade. At each step it carries
/// `(c1, c2, …)` where the head two slots are the proto-RLWE pair
/// being assembled and the tail slots are the as-yet-untouched part
/// of `a_agg`. SPEC.md §6 / Appendix C.
///
/// Phase 4 status: opaque placeholder.
pub struct CollapseState<'a> {
    /// Reserved for the Phase 7 implementation. Kept private so the
    /// invariant set is small enough to audit.
    _phantom: std::marker::PhantomData<&'a ()>,
}
