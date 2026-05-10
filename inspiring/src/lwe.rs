//! LWE ciphertext type and the LWE→RLWE embedding.
//!
//! The embedding is **bit-identical** to the one used in CDKS \[18\] —
//! we keep this primitive verbatim. See SPEC.md §9.f for the rationale.
//!
//! Stage 1 status: LWE shape validation and the Equation 1 embedding are
//! implemented.

use spiral_rs::poly::{PolyMatrix, PolyMatrixRaw};

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
    pub fn validate(&self, params: &RlweParams) -> Result<(), InspiringError> {
        if self.inner.len() != params.d {
            return Err(InspiringError::LweShape(format!(
                "expected {} LWE ciphertexts, got {}",
                params.d,
                self.inner.len()
            )));
        }

        for (idx, ct) in self.inner.iter().enumerate() {
            if ct.a.len() != params.d {
                return Err(InspiringError::LweShape(format!(
                    "ciphertext {idx} has a length {}, expected {}",
                    ct.a.len(),
                    params.d
                )));
            }
        }

        Ok(())
    }
}

/// `ã := Σ_{i=0}^{d-1} a[i] · X^{-i}` (paper Equation 1, the LWE-to-RLWE
/// embedding for the random component). Matches CDKS \[18\] verbatim.
///
pub fn a_tilde<'a>(params: &'a RlweParams, a: &[u64]) -> PolyMatrixRaw<'a> {
    assert_eq!(
        a.len(),
        params.d,
        "lwe::a_tilde expects an LWE vector of length d"
    );

    let mut out = PolyMatrixRaw::zero(&params.spiral, 1, 1);
    let poly = out.get_poly_mut(0, 0);
    poly[0] = a[0] % params.q;
    for (i, coeff) in a.iter().enumerate().skip(1) {
        let reduced = coeff % params.q;
        poly[params.d - i] = if reduced == 0 { 0 } else { params.q - reduced };
    }
    out
}

/// `b̃ := b · X^0` — a constant polynomial. Trivial half of Equation 1.
///
pub fn b_tilde<'a>(params: &'a RlweParams, b: u64) -> PolyMatrixRaw<'a> {
    PolyMatrixRaw::single_value(&params.spiral, b % params.q)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::GadgetParams;

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

    fn coeffs(poly: &PolyMatrixRaw<'_>) -> Vec<u64> {
        poly.get_poly(0, 0).to_vec()
    }

    #[test]
    fn validate_accepts_exactly_d_ciphertexts_with_d_coefficients_each() {
        let params = params();
        let batch = LweBatch {
            inner: (0..params.d)
                .map(|_| LweCiphertext {
                    a: vec![0; params.d],
                    b: 0,
                })
                .collect(),
        };

        assert!(batch.validate(&params).is_ok());
    }

    #[test]
    fn validate_rejects_wrong_batch_length() {
        let params = params();
        let batch = LweBatch {
            inner: vec![LweCiphertext {
                a: vec![0; params.d],
                b: 0,
            }],
        };

        assert!(matches!(
            batch.validate(&params),
            Err(InspiringError::LweShape(_))
        ));
    }

    #[test]
    fn validate_rejects_wrong_ciphertext_length() {
        let params = params();
        let mut batch = LweBatch {
            inner: (0..params.d)
                .map(|_| LweCiphertext {
                    a: vec![0; params.d],
                    b: 0,
                })
                .collect(),
        };
        batch.inner[3].a.pop();

        assert!(matches!(
            batch.validate(&params),
            Err(InspiringError::LweShape(_))
        ));
    }

    #[test]
    fn a_tilde_embeds_lwe_a_with_negative_exponents() {
        let params = params();
        let a = vec![1, 2, 0, 4, 5, 6, 7, 8];

        assert_eq!(
            coeffs(&a_tilde(&params, &a)),
            vec![1, 12281, 12282, 12283, 12284, 12285, 0, 12287]
        );
    }

    #[test]
    fn a_tilde_reduces_inputs_mod_q() {
        let params = params();
        let a = vec![12290, 12291, 12289, 12293, 12294, 12295, 12296, 12297];

        assert_eq!(
            coeffs(&a_tilde(&params, &a)),
            vec![1, 12281, 12282, 12283, 12284, 12285, 0, 12287]
        );
    }

    #[test]
    #[should_panic(expected = "lwe::a_tilde expects an LWE vector of length d")]
    fn a_tilde_panics_on_wrong_shape() {
        let params = params();
        let _ = a_tilde(&params, &[1, 2, 3]);
    }

    #[test]
    fn b_tilde_returns_constant_polynomial_reduced_mod_q() {
        let params = params();

        assert_eq!(
            coeffs(&b_tilde(&params, params.q + 17)),
            vec![17, 0, 0, 0, 0, 0, 0, 0]
        );
    }
}
