"""Stage 6: tests for the gadget vector ``g_z`` and its inverse ``g_z^{-1}``.

Six test groups, layered from primitive to composite:

1. ``TestGadgetVector`` -- ``g_z = [1, z, z^2, ...]`` shape and contents.

2. ``TestSignedDigitDecompositionScalar`` -- the core ``g_z^{-1}`` operator
   on scalars: shape, bounds, reconstruction, KAT for boundary cases.

3. ``TestSignedDigitDecompositionPoly`` -- coefficient-wise extension to
   polynomials. The shape is ``ell x d`` (rows are gadget levels) which
   matches how Stage 7's ``KS.Switch`` consumes them.

4. ``TestRingLevelReconstruction`` -- the **gadget identity**::

        p = sum_i (z^i * digits[i])    in R_q

   This is the equation that makes ``KS.Switch`` work: when ``KS.Switch``
   multiplies ``digits`` by the gadget-encoded secret-shares, the gadget
   identity is what makes the original ``s_in * a`` term reconstitute
   itself for cancellation against the secret-key term.

5. ``TestNoiseBoundOnDigits`` -- the Theorem 2 input::

        max |digit| <= z / 2

   Pinned down explicitly because the ``z^2 / 4`` term in the per-step
   noise bound flows directly from this inequality.

6. ``TestEdgeAndBoundary`` -- defensive checks for ``x = 0``, ``x = q-1``,
   negative ``x``, and ``x > q``. The ``% q`` normalization makes these
   degenerate cases of the random-x test, but exercising them by name
   catches "I refactored away the % q" regressions.
"""

from __future__ import annotations

import random

import pytest

from inspiring_oracle.gadget import gz, gz_inv_poly, gz_inv_scalar
from inspiring_oracle.params import ORACLE_SMALL, ORACLE_TINY
from inspiring_oracle.ring import add, scalar_mul


@pytest.fixture
def rng() -> random.Random:
    return random.Random(0xCABBA6E)


# --------------------------------------------------------------------------
# 1. Gadget vector
# --------------------------------------------------------------------------


class TestGadgetVector:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_length_equals_ell(self, params) -> None:
        assert len(gz(params)) == params.ell

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_first_entry_is_one(self, params) -> None:
        assert gz(params)[0] == 1

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_entries_are_powers_of_z_mod_q(self, params) -> None:
        g = gz(params)
        for i in range(params.ell):
            assert g[i] == pow(params.z, i, params.q), (
                f"g[{i}] = {g[i]} != z^{i} = {pow(params.z, i, params.q)}"
            )


# --------------------------------------------------------------------------
# 2. Signed digit decomposition (scalar)
# --------------------------------------------------------------------------


class TestSignedDigitDecompositionScalar:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_length_equals_ell(self, params, rng) -> None:
        for _ in range(20):
            x = rng.randrange(params.q)
            assert len(gz_inv_scalar(x, params)) == params.ell

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_digits_in_signed_range(self, params, rng) -> None:
        """Every digit satisfies ``-z/2 <= d < z/2``."""
        z = params.z
        for _ in range(2000):
            x = rng.randrange(params.q)
            for d in gz_inv_scalar(x, params):
                assert -(z // 2) <= d < (z // 2), (
                    f"digit {d} out of [-{z//2}, {z//2}) for x={x}"
                )

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_reconstruction_random(self, params, rng) -> None:
        """1000 random ``x``: ``sum_i z^i * digits[i] == x mod q``."""
        z, q = params.z, params.q
        for _ in range(1000):
            x = rng.randrange(q)
            digits = gz_inv_scalar(x, params)
            recon = sum(pow(z, i, q) * d for i, d in enumerate(digits)) % q
            assert recon == x, (
                f"reconstruction failed: x={x}, digits={digits}, got {recon}"
            )

    def test_kat_zero(self) -> None:
        """``x = 0`` decomposes to all zeros."""
        assert gz_inv_scalar(0, ORACLE_TINY) == [0, 0, 0, 0, 0]

    def test_kat_one(self) -> None:
        """``x = 1`` decomposes to ``[1, 0, 0, 0, 0]``."""
        assert gz_inv_scalar(1, ORACLE_TINY) == [1, 0, 0, 0, 0]

    def test_kat_below_z_half_kept_positive(self) -> None:
        """``x in [1, z/2)`` keeps its positive digit at level 0."""
        for x in (1, 2, 3):  # z/2 = 4 at ORACLE_TINY
            assert gz_inv_scalar(x, ORACLE_TINY) == [x, 0, 0, 0, 0]

    def test_kat_z_half_rebalances_negative(self) -> None:
        """``x = z/2`` is exactly the boundary: digit becomes ``-z/2``,
        carry ``+1`` to next level. For ``z = 8``: ``digit_0 = -4``,
        ``digit_1 = 1``, reconstructs to ``-4 + 1*8 = 4``.
        """
        assert gz_inv_scalar(4, ORACLE_TINY) == [-4, 1, 0, 0, 0]

    def test_kat_above_z_half_rebalances(self) -> None:
        """For ``z = 8``: ``x = 5`` -> ``[-3, 1, 0, 0, 0]``, since
        ``-3 + 1*8 = 5``.
        """
        assert gz_inv_scalar(5, ORACLE_TINY) == [-3, 1, 0, 0, 0]

    def test_kat_z_minus_one(self) -> None:
        """For ``z = 8``: ``x = 7`` (= z-1) -> ``[-1, 1, 0, 0, 0]``."""
        assert gz_inv_scalar(7, ORACLE_TINY) == [-1, 1, 0, 0, 0]

    def test_kat_exact_z(self) -> None:
        """For ``z = 8``: ``x = 8`` -> ``[0, 1, 0, 0, 0]``."""
        assert gz_inv_scalar(8, ORACLE_TINY) == [0, 1, 0, 0, 0]

    def test_kat_double_rebalance(self) -> None:
        """``x = 12`` at ``z = 8``: digit_0 = 12%8 = 4, rebalance to -4
        with carry. New x = 12 + 8 = 20, /8 = 2 (note: NOT 1, the carry
        bumped it). digit_1 = 2, no rebalance. Reconstruct: -4 + 2*8 = 12.
        """
        assert gz_inv_scalar(12, ORACLE_TINY) == [-4, 2, 0, 0, 0]


# --------------------------------------------------------------------------
# 3. Signed digit decomposition (polynomial)
# --------------------------------------------------------------------------


class TestSignedDigitDecompositionPoly:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_shape_is_ell_by_d(self, params, rng) -> None:
        """Returned shape is ``[ell rows][d cols]``."""
        p = [rng.randrange(params.q) for _ in range(params.d)]
        digits = gz_inv_poly(p, params)
        assert len(digits) == params.ell
        for row in digits:
            assert len(row) == params.d

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_coefficient_wise_reconstruction(self, params, rng) -> None:
        """For each coefficient ``k`` of ``p``, ``sum_i z^i * digits[i][k] == p[k] mod q``."""
        z, q = params.z, params.q
        for _ in range(50):
            p = [rng.randrange(q) for _ in range(params.d)]
            digits = gz_inv_poly(p, params)
            for k in range(params.d):
                recon = (
                    sum(pow(z, i, q) * digits[i][k] for i in range(params.ell))
                    % q
                )
                assert recon == p[k], (
                    f"coef {k} reconstruction failed: p[{k}] = {p[k]}, "
                    f"got {recon}"
                )

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_all_digits_in_signed_range(self, params, rng) -> None:
        z = params.z
        for _ in range(50):
            p = [rng.randrange(params.q) for _ in range(params.d)]
            for row in gz_inv_poly(p, params):
                assert all(-(z // 2) <= d < (z // 2) for d in row)


# --------------------------------------------------------------------------
# 4. The ring-level gadget identity (this is what makes KS.Switch work)
# --------------------------------------------------------------------------


class TestRingLevelReconstruction:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_p_equals_sum_of_z_i_times_digits_i(self, params, rng) -> None:
        """The ring identity used in Stage 7's ``KS.Switch`` derivation::

            p == sum_i (z^i * digits[i])   in R_q

        Computed using the same ring primitives KS.Switch uses (``add``,
        ``scalar_mul``), so this test pins down compatibility with the
        downstream consumer.
        """
        q = params.q
        for _ in range(50):
            p = [rng.randrange(q) for _ in range(params.d)]
            digits = gz_inv_poly(p, params)
            recon = [0] * params.d
            for i in range(params.ell):
                recon = add(recon, scalar_mul(digits[i], pow(params.z, i, q), q), q)
            assert recon == p, f"ring reconstruction failed: p={p}, recon={recon}"


# --------------------------------------------------------------------------
# 5. The size bound that drives Theorem 2's z^2/4 term
# --------------------------------------------------------------------------


class TestNoiseBoundOnDigits:
    """Pin down the size bound ``|digit| <= z/2`` that drives the
    Theorem 2 per-step noise factor of ``z^2 / 4`` (vs. unsigned's
    ``z^2``). This is the entire reason we use signed digits.
    """

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_max_abs_digit_at_most_z_half_scalar(self, params, rng) -> None:
        z = params.z
        max_abs = 0
        for _ in range(5000):
            x = rng.randrange(params.q)
            for d in gz_inv_scalar(x, params):
                max_abs = max(max_abs, abs(d))
        assert max_abs <= z // 2, f"max |digit| = {max_abs} > z/2 = {z//2}"

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_max_abs_digit_actually_reached(self, params) -> None:
        """The bound is **tight**: there exists ``x`` with a digit of
        absolute value exactly ``z/2``. Specifically, ``x = z/2`` decomposes
        to ``[-z/2, 1, 0, ...]``.
        """
        z = params.z
        digits = gz_inv_scalar(z // 2, params)
        assert digits[0] == -(z // 2)


# --------------------------------------------------------------------------
# 6. Edge / boundary inputs
# --------------------------------------------------------------------------


class TestEdgeAndBoundary:
    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_zero_decomposes_to_all_zeros(self, params) -> None:
        assert gz_inv_scalar(0, params) == [0] * params.ell

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_q_minus_one_reconstructs(self, params) -> None:
        x = params.q - 1
        digits = gz_inv_scalar(x, params)
        recon = (
            sum(pow(params.z, i, params.q) * d for i, d in enumerate(digits))
            % params.q
        )
        assert recon == x

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_x_above_q_reduces_first(self, params) -> None:
        """``gz_inv_scalar(q + 5)`` decomposes the same as ``gz_inv_scalar(5)``."""
        assert gz_inv_scalar(params.q + 5, params) == gz_inv_scalar(5, params)

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_negative_x_reduces_first(self, params) -> None:
        """``gz_inv_scalar(-1)`` decomposes the same as ``gz_inv_scalar(q - 1)``."""
        assert gz_inv_scalar(-1, params) == gz_inv_scalar(params.q - 1, params)

    @pytest.mark.parametrize("params", [ORACLE_TINY, ORACLE_SMALL])
    def test_empty_polynomial_returns_ell_empty_rows(self, params) -> None:
        """``gz_inv_poly([])`` returns ``ell`` empty rows -- defensive but cheap."""
        result = gz_inv_poly([], params)
        assert len(result) == params.ell
        assert all(row == [] for row in result)
