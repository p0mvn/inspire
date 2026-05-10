use inspiring::automorph::{h, tau_g_pow, trace};
use inspiring::intermediate::transform;
use inspiring::lwe::{a_tilde, b_tilde, LweBatch};
use inspiring::{GadgetParams, LweCiphertext, RlweParams};
use spiral_rs::poly::{from_ntt_alloc, PolyMatrix};

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
    .expect("valid tiny test parameters")
}

fn raw_coeffs(poly: &spiral_rs::poly::PolyMatrixRaw<'_>) -> Vec<u64> {
    poly.get_poly(0, 0).to_vec()
}

fn ntt_coeffs(poly: &spiral_rs::poly::PolyMatrixNTT<'_>) -> Vec<u64> {
    raw_coeffs(&from_ntt_alloc(poly))
}

fn add_assign(acc: &mut [u64], rhs: &[u64], q: u64) {
    for (out, value) in acc.iter_mut().zip(rhs) {
        *out = (*out + value) % q;
    }
}

fn scalar_mul(poly: &[u64], scalar: u64, q: u64) -> Vec<u64> {
    poly.iter()
        .map(|coeff| ((u128::from(*coeff) * u128::from(scalar)) % u128::from(q)) as u64)
        .collect()
}

fn tau_coeffs(poly: &[u64], exponent: u64, q: u64) -> Vec<u64> {
    let d = poly.len();
    let mut out = vec![0; d];

    for (i, coeff) in poly.iter().enumerate() {
        let exp = (i as u64 * exponent) % (2 * d as u64);
        let (idx, value) = if exp < d as u64 {
            (exp as usize, *coeff % q)
        } else {
            let reduced = coeff % q;
            (
                (exp - d as u64) as usize,
                if reduced == 0 { 0 } else { q - reduced },
            )
        };
        out[idx] = (out[idx] + value) % q;
    }

    out
}

fn negacyclic_mul(lhs: &[u64], rhs: &[u64], q: u64) -> Vec<u64> {
    let d = lhs.len();
    let mut out = vec![0; d];

    for (i, lhs_coeff) in lhs.iter().enumerate() {
        for (j, rhs_coeff) in rhs.iter().enumerate() {
            let product = (u128::from(*lhs_coeff) * u128::from(*rhs_coeff) % u128::from(q)) as u64;
            let degree = i + j;
            if degree < d {
                out[degree] = (out[degree] + product) % q;
            } else if product != 0 {
                out[degree - d] = (out[degree - d] + q - product) % q;
            }
        }
    }

    out
}

fn embedded_a(a: &[u64], q: u64) -> Vec<u64> {
    let d = a.len();
    let mut out = vec![0; d];
    out[0] = a[0] % q;
    for (i, coeff) in a.iter().enumerate().skip(1) {
        let reduced = coeff % q;
        out[d - i] = if reduced == 0 { 0 } else { q - reduced };
    }
    out
}

fn s_hat_from_s(s: &[u64], q: u64) -> Vec<Vec<u64>> {
    let d = s.len();
    let two_d = 2 * d as u64;
    let h_d = h(d);
    let mut out = Vec::with_capacity(d);

    for j in 0..(d / 2) {
        out.push(tau_coeffs(s, tau_g_pow(j, d), q));
    }
    for j in 0..(d / 2) {
        let gj = tau_g_pow(j, d);
        out.push(tau_coeffs(s, (gj * h_d) % two_d, q));
    }

    out
}

#[test]
fn transform_matches_algorithm_1_formula() {
    let params = params();
    let ct = LweCiphertext {
        a: vec![3, 8, 13, 21, 34, 55, 89, 144],
        b: 9876,
    };

    let ictx = transform(&params, &ct);
    let a_tilde = embedded_a(&ct.a, params.q);

    assert_eq!(raw_coeffs(&ictx.b_tilde), vec![ct.b, 0, 0, 0, 0, 0, 0, 0]);
    assert_eq!(ictx.a_hat.len(), params.d);
    for slot in &ictx.a_hat {
        assert_eq!(slot.rows, 1);
        assert_eq!(slot.cols, 1);
    }

    let two_d = 2 * params.d as u64;
    let h_d = h(params.d);
    for j in 0..(params.d / 2) {
        let gj = tau_g_pow(j, params.d);
        let expected_left = scalar_mul(&tau_coeffs(&a_tilde, gj, params.q), params.d_inv, params.q);
        let expected_right = scalar_mul(
            &tau_coeffs(&a_tilde, (gj * h_d) % two_d, params.q),
            params.d_inv,
            params.q,
        );

        assert_eq!(ntt_coeffs(&ictx.a_hat[j]), expected_left);
        assert_eq!(ntt_coeffs(&ictx.a_hat[j + params.d / 2]), expected_right);
    }
}

#[test]
fn transform_decrypts_under_widened_secret_to_constant_message() {
    let params = params();
    let a = vec![7, 2, 19, 4, 11, 6, 5, 17];
    let s = vec![3, 1, 4, 1, 5, 9, 2, 6];
    let message = 37;
    let inner_product = a.iter().zip(&s).fold(0_u64, |acc, (ai, si)| {
        (acc + (u128::from(*ai) * u128::from(*si) % u128::from(params.q)) as u64) % params.q
    });
    let b = (params.q + message - inner_product) % params.q;
    let ct = LweCiphertext { a, b };

    let ictx = transform(&params, &ct);
    let s_hat = s_hat_from_s(&s, params.q);
    let mut decrypted = raw_coeffs(&ictx.b_tilde);
    for (a_hat_slot, s_hat_slot) in ictx.a_hat.iter().zip(&s_hat) {
        add_assign(
            &mut decrypted,
            &negacyclic_mul(&ntt_coeffs(a_hat_slot), s_hat_slot, params.q),
            params.q,
        );
    }

    let mut expected = vec![0; params.d];
    expected[0] = message;
    assert_eq!(decrypted, expected);
}

#[test]
fn a_hat_is_independent_of_b() {
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
fn embedding_helpers_match_stage_1_contract() {
    let params = params();
    let ct = LweCiphertext {
        a: vec![1, 2, 0, 4, 5, 6, 7, 8],
        b: params.q + 9,
    };

    assert_eq!(
        raw_coeffs(&a_tilde(&params, &ct.a)),
        embedded_a(&ct.a, params.q)
    );
    assert_eq!(
        raw_coeffs(&b_tilde(&params, ct.b)),
        vec![9, 0, 0, 0, 0, 0, 0, 0]
    );
}

#[test]
fn trace_keeps_only_scaled_constant_coefficient() {
    let params = params();
    let mut poly = a_tilde(&params, &[4, 1, 9, 2, 6, 5, 3, 8]);
    poly.get_poly_mut(0, 0)[0] = 42;

    let traced = raw_coeffs(&trace(&poly));
    let mut expected = vec![0; params.d];
    expected[0] = (params.d as u64 * 42) % params.q;
    assert_eq!(traced, expected);
}

#[test]
fn params_and_batch_validation_are_enforced() {
    assert!(RlweParams::new(
        7,
        12289,
        4,
        3.2,
        GadgetParams {
            bits_per: 3,
            ell: 5,
        },
    )
    .is_err());
    assert!(RlweParams::new(
        8,
        12288,
        4,
        3.2,
        GadgetParams {
            bits_per: 3,
            ell: 5,
        },
    )
    .is_err());

    let params = params();
    let valid = LweBatch {
        inner: (0..params.d)
            .map(|_| LweCiphertext {
                a: vec![0; params.d],
                b: 0,
            })
            .collect(),
    };
    assert!(valid.validate(&params).is_ok());

    let invalid = LweBatch {
        inner: vec![LweCiphertext {
            a: vec![0; params.d - 1],
            b: 0,
        }],
    };
    assert!(invalid.validate(&params).is_err());
}
