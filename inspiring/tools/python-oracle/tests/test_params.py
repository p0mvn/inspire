"""Stage 0: smoke tests for RlweParams.

Validates the dataclass invariants and confirms both presets satisfy the
Theorem 2 correctness bound. These tests are intentionally small; the real
correctness work begins at Stage 1.
"""

from __future__ import annotations

import math

import pytest

from inspiring_oracle.params import ORACLE_SMALL, ORACLE_TINY, RlweParams


class TestModularInverse:
    def test_d_inv_is_correct_for_oracle_tiny(self) -> None:
        p = ORACLE_TINY
        assert (p.d * p.d_inv) % p.q == 1

    def test_d_inv_is_correct_for_oracle_small(self) -> None:
        p = ORACLE_SMALL
        assert (p.d * p.d_inv) % p.q == 1


class TestValidation:
    def test_even_q_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="q must be odd"):
            RlweParams(d=8, q=12290, p=4, sigma=3.2, z=8, ell=5)

    def test_non_power_of_two_d_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="power of 2"):
            RlweParams(d=12, q=12289, p=4, sigma=3.2, z=8, ell=5)

    def test_d_zero_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="power of 2"):
            RlweParams(d=0, q=12289, p=4, sigma=3.2, z=8, ell=5)

    def test_undersized_gadget_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="gadget too short"):
            # 8^4 = 4096 < q = 12289
            RlweParams(d=8, q=12289, p=4, sigma=3.2, z=8, ell=4)

    def test_negative_sigma_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="sigma"):
            RlweParams(d=8, q=12289, p=4, sigma=-1.0, z=8, ell=5)


class TestNoiseBudget:
    def test_oracle_tiny_satisfies_correctness_bound(self) -> None:
        assert ORACLE_TINY.correctness_ok(num_sigmas=6.0)

    def test_oracle_small_satisfies_correctness_bound(self) -> None:
        assert ORACLE_SMALL.correctness_ok(num_sigmas=6.0)

    def test_noise_budget_formula_is_sqrt_of_theorem2_bound(self) -> None:
        p = ORACLE_TINY
        expected = math.sqrt(p.ell * p.d ** 2 * p.z ** 2 * p.sigma ** 2 / 4.0)
        assert p.noise_budget_sigma == pytest.approx(expected)


class TestDelta:
    def test_delta_is_floor_q_over_p(self) -> None:
        assert ORACLE_TINY.delta == 12289 // 4
        assert ORACLE_SMALL.delta == 65537 // 4
