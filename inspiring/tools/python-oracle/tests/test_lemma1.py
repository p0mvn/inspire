"""Stage 3: tests for the Lemma 1 trace operator.

This is the **conceptual heart of InspiRING**: an operator built purely
from Galois automorphisms that extracts the constant coefficient of a
polynomial. SPEC.md section 3 has the full proof; the tests here verify
both the headline claim (``Tr(p) = d * c_0``) and the proof's intermediate
step (``half_trace`` follows Lemma 5).

Test groups:

1. ``TestHalfTraceLemma5`` -- the "half-trace" intermediate. By Lemma 5
   (with gamma = d/2), summing ``tau_g^j(p)`` for ``j in [0, d/2)`` gives
   ``(d/2) * (c_0 + c_{d/2} * X^{d/2})`` -- everything else cancels. This
   is the firewall for bugs in ``tau`` and the exponent-folding logic.

2. ``TestTraceLemma1`` -- the headline result: ``Tr(p) = d * c_0``.
   Random samples at d=8 and d=16; KAT at d=8.

3. ``TestTraceProofStructure`` -- ``trace(p) == half_trace(p) + tau_h(half_trace(p))``.
   This is the form the SPEC.md proof actually uses; verifying both forms
   produce the same answer pins down the implementation.

4. ``TestTraceLinearity`` -- additive linearity (since each ``tau_g`` is
   linear). Used implicitly throughout Stage 8.

5. ``TestTraceWithDInv`` -- ``d_inv * Tr(p)`` recovers a clean ``c_0 * X^0``.
   This is the actual usage pattern in Stage 8 (TRANSFORM).

6. ``TestTraceCancellations`` -- specific inputs that exercise the
   ``X^{d/2}`` cancellation: a polynomial supported only on ``X^{d/2}``
   has trace zero (Lemma 1's "cancellation step" hits hardest here).
"""

from __future__ import annotations

import random

import pytest

from inspiring_oracle.automorph import (
    half_trace,
    tau_g_pow,
    tau_h,
    trace,
)
from inspiring_oracle.params import ORACLE_SMALL, ORACLE_TINY
from inspiring_oracle.ring import add, scalar_mul


def rand_poly(rng: random.Random, d: int, q: int) -> list[int]:
    return [rng.randrange(q) for _ in range(d)]


@pytest.fixture
def rng():
    return random.Random(0xCABBA6E)


# --------------------------------------------------------------------------
# Lemma 5: the half-trace intermediate
# --------------------------------------------------------------------------


class TestHalfTraceLemma5:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_half_trace_keeps_only_c_0_and_c_d_over_2(self, params, rng) -> None:
        """``pi_{d/2}(p)`` has ``(d/2)*c_0`` at position 0, ``(d/2)*c_{d/2}``
        at position ``d/2``, and zeros everywhere else.
        """
        d, q = params.d, params.q
        for _ in range(50):
            p = rand_poly(rng, d, q)
            ht = half_trace(p, q)
            expected = [0] * d
            expected[0] = (d // 2 * p[0]) % q
            expected[d // 2] = (d // 2 * p[d // 2]) % q
            assert ht == expected, (
                f"half_trace mismatch at d={d}: got {ht}, expected {expected}"
            )

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_half_trace_of_zero_is_zero(self, params) -> None:
        d, q = params.d, params.q
        zero = [0] * d
        assert half_trace(zero, q) == zero

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_half_trace_is_additive(self, params, rng) -> None:
        """``pi_{d/2}`` is linear (sum of linear maps). Used in Stage 8."""
        d, q = params.d, params.q
        for _ in range(20):
            a = rand_poly(rng, d, q)
            b = rand_poly(rng, d, q)
            assert half_trace(add(a, b, q), q) == add(
                half_trace(a, q),
                half_trace(b, q),
                q,
            )

    def test_half_trace_kat_d8(self) -> None:
        """KAT at d=8: pick simple p so the result is checkable by eye.

        p = [1, 0, 0, 0, 5, 0, 0, 0]    (only c_0 = 1 and c_{d/2} = 5)
        Expected:
          half_trace(p) = (d/2)*(c_0 + c_{d/2}*X^{d/2})
                        = 4*(1 + 5*X^4)
                        = [4, 0, 0, 0, 20, 0, 0, 0]
        """
        q = 12289
        p = [1, 0, 0, 0, 5, 0, 0, 0]
        assert half_trace(p, q) == [4, 0, 0, 0, 20, 0, 0, 0]

    def test_half_trace_kat_kills_odd_coefficients_at_d8(self) -> None:
        """KAT at d=8: a polynomial with only X^1 should half-trace to zero.

        At d=8, half_trace's nonzero positions are 0 and 4. X^1 contributes
        to neither, so the entire result is zero.
        """
        d, q = 8, 12289
        p = [0, 7, 0, 0, 0, 0, 0, 0]  # 7 * X
        assert half_trace(p, q) == [0] * d


# --------------------------------------------------------------------------
# Lemma 1: the trace
# --------------------------------------------------------------------------


class TestTraceLemma1:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_trace_equals_d_times_c_0_random(self, params, rng) -> None:
        """The headline claim: ``Tr(p) = d * c_0``."""
        d, q = params.d, params.q
        for _ in range(100):
            p = rand_poly(rng, d, q)
            result = trace(p, q)
            expected = [(d * p[0]) % q] + [0] * (d - 1)
            assert result == expected, (
                f"trace failed at d={d}: got {result}, expected {expected}"
            )

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_trace_of_zero_is_zero(self, params) -> None:
        d, q = params.d, params.q
        assert trace([0] * d, q) == [0] * d

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_trace_of_constant_polynomial(self, params, rng) -> None:
        """A constant polynomial ``c * X^0`` has trace ``d * c * X^0``."""
        d, q = params.d, params.q
        for _ in range(20):
            c = rng.randrange(q)
            const = [c] + [0] * (d - 1)
            assert trace(const, q) == [(d * c) % q] + [0] * (d - 1)

    def test_trace_kat_d8(self) -> None:
        """KAT at d=8, q=12289: every coefficient except c_0 is killed.

        p = [5, 7, 3, 11, 9, 13, 6, 15]
        Expected: trace(p) = [d * c_0, 0, ..., 0] = [40, 0, 0, 0, 0, 0, 0, 0]
        """
        q = 12289
        p = [5, 7, 3, 11, 9, 13, 6, 15]
        assert trace(p, q) == [40, 0, 0, 0, 0, 0, 0, 0]

    def test_trace_kat_d8_with_high_c_0(self) -> None:
        """KAT exercising mod-q wraparound at d=8: c_0 close to q.

        p = [12000, 1, 2, 3, 4, 5, 6, 7], d * c_0 = 8 * 12000 = 96000.
        96000 mod 12289 = 96000 - 7*12289 = 96000 - 86023 = 9977.
        """
        d, q = 8, 12289
        p = [12000, 1, 2, 3, 4, 5, 6, 7]
        expected = [(d * 12000) % q] + [0] * (d - 1)
        assert expected[0] == 9977  # double-check the arithmetic in the test
        assert trace(p, q) == expected


# --------------------------------------------------------------------------
# The proof structure: Tr = pi_{d/2} + tau_h o pi_{d/2}
# --------------------------------------------------------------------------


class TestTraceProofStructure:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_trace_equals_half_trace_plus_tau_h_of_half_trace(
        self, params, rng
    ) -> None:
        """The form SPEC.md's proof actually uses (and our implementation)."""
        d, q = params.d, params.q
        for _ in range(20):
            p = rand_poly(rng, d, q)
            ht = half_trace(p, q)
            assert trace(p, q) == add(ht, tau_h(ht, q), q)

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_trace_equals_full_galois_sum_form(self, params, rng) -> None:
        """The form SPEC.md *states*: sum over j of (tau_g^j(p) + tau_h(tau_g^j(p))).

        This is the literal Lemma 1 formula. Should be byte-identical to
        our half_trace-based implementation by the linearity of tau_h.
        """
        d, q = params.d, params.q
        for _ in range(10):
            p = rand_poly(rng, d, q)
            spec_form = [0] * d
            for j in range(d // 2):
                tj = tau_g_pow(p, j, q)
                spec_form = add(spec_form, tj, q)
                spec_form = add(spec_form, tau_h(tj, q), q)
            assert trace(p, q) == spec_form


# --------------------------------------------------------------------------
# Linearity (used implicitly in Stage 8)
# --------------------------------------------------------------------------


class TestTraceLinearity:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_trace_is_additive(self, params, rng) -> None:
        d, q = params.d, params.q
        for _ in range(20):
            a = rand_poly(rng, d, q)
            b = rand_poly(rng, d, q)
            assert trace(add(a, b, q), q) == add(trace(a, q), trace(b, q), q)

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_trace_commutes_with_scalar_mul(self, params, rng) -> None:
        d, q = params.d, params.q
        for _ in range(20):
            a = rand_poly(rng, d, q)
            k = rng.randrange(q)
            assert trace(scalar_mul(a, k, q), q) == scalar_mul(trace(a, q), k, q)


# --------------------------------------------------------------------------
# The Stage 8 usage pattern: scale by d^{-1} to extract c_0 cleanly
# --------------------------------------------------------------------------


class TestTraceWithDInv:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_d_inv_times_trace_extracts_c_0(self, params, rng) -> None:
        """The actual Stage 8 idiom: ``d_inv * Tr(p) = c_0 * X^0``.

        This is why ``RlweParams`` insists on odd ``q`` -- it needs
        ``d^{-1} mod q`` to exist. Without this, the trace gives ``d*c_0``
        but we can't divide back down to ``c_0`` cleanly.
        """
        d, q, d_inv = params.d, params.q, params.d_inv
        for _ in range(50):
            p = rand_poly(rng, d, q)
            scaled = scalar_mul(trace(p, q), d_inv, q)
            assert scaled == [p[0]] + [0] * (d - 1)


# --------------------------------------------------------------------------
# The X^{d/2} cancellation: where the second half of the trace earns its keep
# --------------------------------------------------------------------------


class TestTraceCancellations:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_trace_kills_pure_x_pow_d_over_2(self, params) -> None:
        """A polynomial supported only on ``X^{d/2}`` has trace = 0.

        ``half_trace`` puts ``(d/2)*c_{d/2}`` at position ``d/2``;
        ``tau_h`` flips its sign (since ``tau_h(X^{d/2}) = -X^{d/2}``);
        their sum cancels. This is the **second half of the Lemma 1 proof**.

        If this test fails, ``tau_h`` is buggy or the trace is not folding
        in the second half correctly.
        """
        d, q = params.d, params.q
        p = [0] * d
        p[d // 2] = 7
        assert trace(p, q) == [0] * d

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_trace_kills_every_non_constant_monomial(self, params) -> None:
        """For each k in [1, d), trace of ``c * X^k`` is zero."""
        d, q = params.d, params.q
        for k in range(1, d):
            p = [0] * d
            p[k] = 5
            assert trace(p, q) == [0] * d, (
                f"trace failed to kill X^{k} at d={d}"
            )

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_trace_preserves_only_constant_coefficient(self, params, rng) -> None:
        """Equivalent restatement of Lemma 1: trace ignores everything but c_0."""
        d, q = params.d, params.q
        for _ in range(20):
            c0 = rng.randrange(q)
            # Vary the non-constant coefficients arbitrarily; trace should
            # depend only on c_0.
            p1 = [c0] + [rng.randrange(q) for _ in range(d - 1)]
            p2 = [c0] + [rng.randrange(q) for _ in range(d - 1)]
            assert trace(p1, q) == trace(p2, q)
