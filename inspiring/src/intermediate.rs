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
//! Stage 1 status: [`transform`] is implemented; [`aggregate`] remains a
//! Stage 2 stub.

use spiral_rs::poly::{to_ntt_alloc, PolyMatrix, PolyMatrixNTT, PolyMatrixRaw};

use crate::automorph::{h, tau_g_pow, tau_raw};
use crate::lwe::LweCiphertext;
use crate::lwe::{a_tilde, b_tilde};
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
pub fn transform<'a>(params: &'a RlweParams, ct: &LweCiphertext) -> IRCtx<'a> {
    assert_eq!(
        ct.a.len(),
        params.d,
        "intermediate::transform expects an LWE vector of length d"
    );

    let a_tilde = a_tilde(params, &ct.a);
    let b_tilde = b_tilde(params, ct.b);
    let mut a_hat = Vec::with_capacity(params.d);
    a_hat.resize_with(params.d, || PolyMatrixNTT::zero(&params.spiral, 1, 1));

    let two_d = 2 * params.d as u64;
    let h_d = h(params.d);
    for j in 0..(params.d / 2) {
        let gj = tau_g_pow(j, params.d);
        a_hat[j] = scaled_tau_ntt(params, &a_tilde, gj);
        a_hat[j + params.d / 2] = scaled_tau_ntt(params, &a_tilde, (gj * h_d) % two_d);
    }

    IRCtx { a_hat, b_tilde }
}

fn scaled_tau_ntt<'a>(
    params: &'a RlweParams,
    a_tilde: &PolyMatrixRaw<'a>,
    exponent: u64,
) -> PolyMatrixNTT<'a> {
    let mut raw = tau_raw(a_tilde, exponent);
    for coeff in raw.get_poly_mut(0, 0) {
        *coeff = ((*coeff as u128 * params.d_inv as u128) % params.q as u128) as u64;
    }
    to_ntt_alloc(&raw)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::automorph::h;
    use crate::params::GadgetParams;
    use spiral_rs::poly::from_ntt_alloc;

    fn params() -> RlweParams {
        RlweParams::new(
            8,
            12289,
            4,
            3.2,
            GadgetParams {
                bits_per: 3,
                ell: 5,
            },
        )
        .expect("valid params")
    }

    fn raw_coeffs(poly: &PolyMatrixRaw<'_>) -> Vec<u64> {
        poly.get_poly(0, 0).to_vec()
    }

    fn ntt_coeffs(poly: &PolyMatrixNTT<'_>) -> Vec<u64> {
        raw_coeffs(&from_ntt_alloc(poly))
    }

    fn scalar_mul(poly: &[u64], scalar: u64, q: u64) -> Vec<u64> {
        poly.iter()
            .map(|coeff| ((u128::from(*coeff) * u128::from(scalar)) % u128::from(q)) as u64)
            .collect()
    }

    #[test]
    fn scaled_tau_ntt_applies_automorphism_and_d_inverse_scaling() {
        let params = params();
        let a_tilde = a_tilde(&params, &[3, 8, 13, 21, 34, 55, 89, 144]);
        let exponent = tau_g_pow(2, params.d);
        let expected = scalar_mul(
            raw_coeffs(&tau_raw(&a_tilde, exponent)).as_slice(),
            params.d_inv,
            params.q,
        );

        assert_eq!(
            ntt_coeffs(&scaled_tau_ntt(&params, &a_tilde, exponent)),
            expected
        );
    }

    #[test]
    fn transform_returns_constant_b_tilde_and_d_a_hat_slots() {
        let params = params();
        let ct = LweCiphertext {
            a: vec![3, 8, 13, 21, 34, 55, 89, 144],
            b: params.q + 7,
        };

        let ictx = transform(&params, &ct);

        assert_eq!(ictx.a_hat.len(), params.d);
        assert_eq!(raw_coeffs(&ictx.b_tilde), vec![7, 0, 0, 0, 0, 0, 0, 0]);
        for slot in &ictx.a_hat {
            assert_eq!(slot.rows, 1);
            assert_eq!(slot.cols, 1);
        }
    }

    #[test]
    fn transform_a_hat_matches_algorithm_1_slots() {
        let params = params();
        let ct = LweCiphertext {
            a: vec![3, 8, 13, 21, 34, 55, 89, 144],
            b: 9876,
        };
        let a_tilde = a_tilde(&params, &ct.a);
        let ictx = transform(&params, &ct);
        let h_d = h(params.d);
        let two_d = 2 * params.d as u64;

        for j in 0..(params.d / 2) {
            let gj = tau_g_pow(j, params.d);
            assert_eq!(
                ntt_coeffs(&ictx.a_hat[j]),
                ntt_coeffs(&scaled_tau_ntt(&params, &a_tilde, gj))
            );
            assert_eq!(
                ntt_coeffs(&ictx.a_hat[j + params.d / 2]),
                ntt_coeffs(&scaled_tau_ntt(&params, &a_tilde, (gj * h_d) % two_d))
            );
        }
    }

    #[test]
    fn transform_a_hat_is_independent_of_b() {
        let params = params();
        let a = vec![1, 1, 2, 3, 5, 8, 13, 21];
        let left = transform(&params, &LweCiphertext { a: a.clone(), b: 1 });
        let right = transform(&params, &LweCiphertext { a, b: 12288 });

        for (left_slot, right_slot) in left.a_hat.iter().zip(&right.a_hat) {
            assert_eq!(ntt_coeffs(left_slot), ntt_coeffs(right_slot));
        }
        assert_ne!(raw_coeffs(&left.b_tilde), raw_coeffs(&right.b_tilde));
    }

    #[test]
    #[should_panic(expected = "intermediate::transform expects an LWE vector of length d")]
    fn transform_panics_on_wrong_lwe_shape() {
        let params = params();
        let _ = transform(
            &params,
            &LweCiphertext {
                a: vec![1, 2, 3],
                b: 0,
            },
        );
    }
}
