"""Stage 2: tests for Galois automorphisms ``tau_g``, ``tau_h``.

Test groups:

1. ``TestMonomialKAT_tau_g_d8`` -- hand-checked table for ``tau_5`` acting
   on every monomial in R_q at d=8. This is the firewall against any bug
   in the exponent-folding logic of ``tau``.
2. ``TestMonomialKAT_tau_h_d8`` -- hand-checked table for ``tau_{15}`` at
   d=8 (the ``X^k -> -X^{d-k}`` flip).
3. ``TestOrders`` -- ``ord(tau_g) = d/2`` and ``ord(tau_h) = 2``.
4. ``TestComposition`` -- ``tau_a o tau_b == tau_{a*b mod 2d}``.
5. ``TestCommutativity`` -- ``tau_g`` and ``tau_h`` commute (since
   ``Z*_{2d}`` is abelian).
6. ``TestRingHomomorphism`` -- ``tau_g`` (and ``tau_h``) preserve ``+``
   and ``*`` (the defining property of an automorphism).
7. ``TestSpecialValues`` -- ``tau_g(0) = 0``, ``tau_g(c) = c`` for any
   constant polynomial.
8. ``TestGaloisGroupStructure`` -- the d distinct elements
   ``{5^j * h^b mod 2d}`` exhaust ``Z*_{2d}``.
9. ``TestArbitraryG`` -- ``tau`` works for any odd ``g``, not just ``G``
   or ``h``.
"""

from __future__ import annotations

import random

import pytest

from inspiring_oracle.automorph import G, h, tau, tau_g, tau_g_pow, tau_h
from inspiring_oracle.params import ORACLE_SMALL, ORACLE_TINY
from inspiring_oracle.ring import add, mul


def rand_poly(rng: random.Random, d: int, q: int) -> list[int]:
    return [rng.randrange(q) for _ in range(d)]


def x_pow(k: int, d: int) -> list[int]:
    out = [0] * d
    out[k] = 1
    return out


@pytest.fixture
def rng():
    return random.Random(0xBEEF)


# --------------------------------------------------------------------------
# Hand-checked monomial tables
# --------------------------------------------------------------------------


class TestMonomialKAT_tau_g_d8:
    """Table for ``tau_5`` at d=8.

    Each row says ``tau_g(X^i) = sign * X^position``. Computed by hand:
    ``(i * 5) mod 16``, then folded into ``[0, 8)`` with sign flip when
    the reduced exponent is ``>= 8``.

        i  i*5  i*5 mod 16  fold   tau_g(X^i)
        -  ---  ----------  -----  ----------
        0   0      0          +0    +X^0
        1   5      5          +5    +X^5
        2  10     10        -(10-8)  -X^2
        3  15     15        -(15-8)  -X^7
        4  20      4          +4    +X^4
        5  25      9        -(9-8)   -X^1
        6  30     14        -(14-8)  -X^6
        7  35      3          +3    +X^3
    """

    KAT = [
        (0, 0, +1),
        (1, 5, +1),
        (2, 2, -1),
        (3, 7, -1),
        (4, 4, +1),
        (5, 1, -1),
        (6, 6, -1),
        (7, 3, +1),
    ]

    @pytest.mark.parametrize("i,pos,sign", KAT)
    def test_tau_g_on_monomial(self, i: int, pos: int, sign: int) -> None:
        d, q = 8, 12289
        result = tau_g(x_pow(i, d), q)
        expected = [0] * d
        expected[pos] = sign % q
        assert result == expected, (
            f"tau_g(X^{i}): got {result}, expected {expected}"
        )


class TestMonomialKAT_tau_h_d8:
    """Table for ``tau_h`` at d=8 (``h = 15``, so ``i*h mod 16 = -i mod 16``).

    ``tau_h(X^k) = X^{-k}`` interpreted in R_q. For ``k = 0`` this is
    ``X^0 = 1``; for ``k > 0`` we have ``X^{-k} = X^{2d - k} = X^{d-k} * X^d
    = -X^{d-k}``.
    """

    @pytest.mark.parametrize("k", list(range(8)))
    def test_tau_h_on_monomial(self, k: int) -> None:
        d, q = 8, 12289
        result = tau_h(x_pow(k, d), q)
        expected = [0] * d
        if k == 0:
            expected[0] = 1
        else:
            expected[d - k] = q - 1  # -1 mod q
        assert result == expected, (
            f"tau_h(X^{k}): got {result}, expected {expected}"
        )


# --------------------------------------------------------------------------
# Orders
# --------------------------------------------------------------------------


class TestOrders:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_g_iterated_d_over_2_times_is_identity(self, params, rng) -> None:
        d, q = params.d, params.q
        a = rand_poly(rng, d, q)
        b = a
        for _ in range(d // 2):
            b = tau_g(b, q)
        assert b == a

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_g_pow_d_over_2_is_identity(self, params, rng) -> None:
        d, q = params.d, params.q
        a = rand_poly(rng, d, q)
        assert tau_g_pow(a, d // 2, q) == a

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_g_pow_smaller_exponents_are_not_identity(self, params, rng) -> None:
        """For 1 <= j < d/2, ``tau_g^j != identity`` (so the order is exactly d/2).

        The order of 5 in ``Z*_{2d}`` is exactly ``d/2``; together with
        ``tau`` being a faithful action, no smaller power is the identity.
        """
        d, q = params.d, params.q
        a = rand_poly(rng, d, q)
        for j in range(1, d // 2):
            assert tau_g_pow(a, j, q) != a, (
                f"tau_g^{j} should not be identity but it equals input"
            )

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_h_squared_is_identity(self, params, rng) -> None:
        d, q = params.d, params.q
        a = rand_poly(rng, d, q)
        assert tau_h(tau_h(a, q), q) == a


# --------------------------------------------------------------------------
# Composition: tau_a o tau_b = tau_{a*b mod 2d}
# --------------------------------------------------------------------------


class TestComposition:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_g_pow_zero_is_identity(self, params, rng) -> None:
        d, q = params.d, params.q
        a = rand_poly(rng, d, q)
        assert tau_g_pow(a, 0, q) == a

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_g_pow_one_equals_tau_g(self, params, rng) -> None:
        d, q = params.d, params.q
        a = rand_poly(rng, d, q)
        assert tau_g_pow(a, 1, q) == tau_g(a, q)

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_g_pow_composes_via_addition_of_exponents(self, params, rng) -> None:
        d, q = params.d, params.q
        a = rand_poly(rng, d, q)
        for i in range(d // 2):
            for j in range(d // 2):
                lhs = tau_g_pow(tau_g_pow(a, i, q), j, q)
                rhs = tau_g_pow(a, (i + j) % (d // 2), q)
                assert lhs == rhs, f"tau_g^{i} o tau_g^{j} != tau_g^{(i+j)}"

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_composition_via_index_multiplication(self, params, rng) -> None:
        """``tau_{g1} o tau_{g2} == tau_{g1 * g2 mod 2d}`` for any odd indices."""
        d, q = params.d, params.q
        a = rand_poly(rng, d, q)
        odd_indices = [g for g in range(1, 2 * d, 2)]
        # Test a representative subset to keep it fast
        for g1 in odd_indices[:6]:
            for g2 in odd_indices[:6]:
                lhs = tau(tau(a, g2, q), g1, q)
                rhs = tau(a, (g1 * g2) % (2 * d), q)
                assert lhs == rhs, (
                    f"tau_{g1} o tau_{g2} != tau_{(g1*g2) % (2*d)}"
                )

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_with_g_one_is_identity(self, params, rng) -> None:
        """``tau_1(p)(X) = p(X^1) = p`` -- the identity automorphism."""
        d, q = params.d, params.q
        a = rand_poly(rng, d, q)
        assert tau(a, 1, q) == a


# --------------------------------------------------------------------------
# Commutativity (Z*_{2d} is abelian)
# --------------------------------------------------------------------------


class TestCommutativity:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_g_and_tau_h_commute(self, params, rng) -> None:
        d, q = params.d, params.q
        for _ in range(20):
            a = rand_poly(rng, d, q)
            assert tau_g(tau_h(a, q), q) == tau_h(tau_g(a, q), q)


# --------------------------------------------------------------------------
# Ring-homomorphism property: tau_g(a + b) = tau_g(a) + tau_g(b)
#                             tau_g(a * b) = tau_g(a) * tau_g(b)
# --------------------------------------------------------------------------


class TestRingHomomorphism:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_g_preserves_addition(self, params, rng) -> None:
        d, q = params.d, params.q
        for _ in range(20):
            a = rand_poly(rng, d, q)
            b = rand_poly(rng, d, q)
            assert tau_g(add(a, b, q), q) == add(tau_g(a, q), tau_g(b, q), q)

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_g_preserves_multiplication(self, params, rng) -> None:
        d, q = params.d, params.q
        for _ in range(20):
            a = rand_poly(rng, d, q)
            b = rand_poly(rng, d, q)
            assert tau_g(mul(a, b, q), q) == mul(tau_g(a, q), tau_g(b, q), q)

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_h_preserves_addition(self, params, rng) -> None:
        d, q = params.d, params.q
        for _ in range(20):
            a = rand_poly(rng, d, q)
            b = rand_poly(rng, d, q)
            assert tau_h(add(a, b, q), q) == add(tau_h(a, q), tau_h(b, q), q)

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_h_preserves_multiplication(self, params, rng) -> None:
        d, q = params.d, params.q
        for _ in range(20):
            a = rand_poly(rng, d, q)
            b = rand_poly(rng, d, q)
            assert tau_h(mul(a, b, q), q) == mul(tau_h(a, q), tau_h(b, q), q)


# --------------------------------------------------------------------------
# Special values
# --------------------------------------------------------------------------


class TestSpecialValues:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_g_fixes_zero(self, params) -> None:
        d, q = params.d, params.q
        zero = [0] * d
        assert tau_g(zero, q) == zero
        assert tau_h(zero, q) == zero

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_fixes_constant_polynomial(self, params, rng) -> None:
        """``tau_g(c) = c`` for any constant polynomial ``c`` in Z_q.

        Constants live in the fixed field of every Galois automorphism.
        """
        d, q = params.d, params.q
        c = rng.randrange(q)
        const = [c] + [0] * (d - 1)
        assert tau_g(const, q) == const
        assert tau_h(const, q) == const
        for j in range(d // 2):
            assert tau_g_pow(const, j, q) == const


# --------------------------------------------------------------------------
# Galois group structure
# --------------------------------------------------------------------------


class TestGaloisGroupStructure:
    """The d elements ``{5^j * h^b mod 2d}`` exhaust ``Z*_{2d}``."""

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_galois_group_has_d_distinct_elements(self, params) -> None:
        d = params.d
        elements = set()
        for j in range(d // 2):
            for b in range(2):
                g = (pow(G, j, 2 * d) * pow(h(d), b, 2 * d)) % (2 * d)
                elements.add(g)
        assert len(elements) == d, f"expected d={d} distinct elements, got {len(elements)}"

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_all_galois_group_elements_are_odd(self, params) -> None:
        """``Z*_{2d}`` for ``2d`` a power of 2 is exactly the odd residues."""
        d = params.d
        for j in range(d // 2):
            for b in range(2):
                g = (pow(G, j, 2 * d) * pow(h(d), b, 2 * d)) % (2 * d)
                assert g % 2 == 1, f"element {g} is not odd"

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_galois_elements_match_odd_residues_in_2d(self, params) -> None:
        """The d generated elements equal exactly the d odd residues in ``[1, 2d)``."""
        d = params.d
        generated = set()
        for j in range(d // 2):
            for b in range(2):
                g = (pow(G, j, 2 * d) * pow(h(d), b, 2 * d)) % (2 * d)
                generated.add(g)
        odd_residues = set(range(1, 2 * d, 2))
        assert generated == odd_residues


# --------------------------------------------------------------------------
# Inverse: tau_g^{-1} = tau_{g^{-1} mod 2d}
# --------------------------------------------------------------------------


class TestInverse:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_g_pow_j_then_d_over_2_minus_j_is_identity(self, params, rng) -> None:
        """``tau_g^j composed with tau_g^{(d/2) - j}`` is identity in the cyclic group."""
        d, q = params.d, params.q
        a = rand_poly(rng, d, q)
        for j in range(d // 2):
            forward = tau_g_pow(a, j, q)
            recovered = tau_g_pow(forward, (d // 2) - j, q)
            assert recovered == a, (
                f"tau_g^{j} composed with tau_g^{(d//2) - j} != identity"
            )

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_tau_h_is_self_inverse(self, params, rng) -> None:
        d, q = params.d, params.q
        a = rand_poly(rng, d, q)
        assert tau_h(tau_h(a, q), q) == a
