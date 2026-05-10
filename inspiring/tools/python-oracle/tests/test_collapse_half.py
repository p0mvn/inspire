"""Tests for ``collapse.collapse_half`` -- the InspiRING crown jewel.

``CollapseHalf`` is where InspiRING's headline efficiency win shows up:
all ``d/2 - 1`` per-step key-switching matrices are derived from the
**single** base matrix ``K_g`` by entry-wise application of Galois
automorphisms (``apply_automorph``, Stage 7). CDKS by contrast needs
``log d`` distinct base matrices.

Test groups
-----------

1. **TestInputValidation** -- rejects wrong input sizes; rejects invalid
   ``rho`` values.

2. **TestOutputShape** -- output is a single polynomial pair; both have
   length ``d`` and coefficients in ``[0, q)``.

3. **TestCorrectnessLeftHalf** -- the firewall for ``rho = "identity"``:
   build a wider ciphertext under ``s_hat_left[j] = tau_g^j(s_tilde)``,
   collapse, decrypt under ``s_tilde``, recover the plaintext.

4. **TestCorrectnessRightHalf** -- the firewall for ``rho = "tau_h"``:
   build under ``s_hat_right[j] = tau_h(tau_g^j(s_tilde))``, collapse,
   decrypt under ``tau_h(s_tilde)``, recover the plaintext.

5. **TestSwitchCallCount** -- each ``collapse_half`` does **exactly**
   ``d/2 - 1`` ``KS.Switch`` calls. Stage 12 needs both halves to
   contribute exactly this count to hit the SPEC's ``d - 1`` total.

6. **TestRandomComponentInvariant** -- the output ``a_out`` depends only
   on ``(a, K_g, rho)``, not on ``b``. The basis for offline-precomputing
   the entire ``a``-trace of the collapse pipeline.

7. **TestNoiseGrowthBounded** -- after the full half-collapse, residual
   noise per coefficient is bounded by ``Delta / 2`` (decryption budget).
   The variance bound is ``(d/2 - 1) * sigma_one_ks^2 ~= sigma_pack^2 / 2``.

8. **TestDeterminism** -- pure function of inputs.

9. **TestUsesOnlyOneBaseMatrix** (structural) -- ``collapse_half`` takes
   exactly **one** ``K_g`` matrix as input; this is encoded in the API
   and is the property that distinguishes InspiRING from CDKS at the
   type level. Verified by the type signature (no runtime test needed),
   but documented here for future Rust-port preservation.

10. **TestEdgeCases** -- all-zero message; all-max message.
"""

from __future__ import annotations

import math
import random

import pytest

from inspiring_oracle import key_switching, rlwe
from inspiring_oracle.automorph import G, h, tau
from inspiring_oracle.collapse import collapse_half
from inspiring_oracle.params import ORACLE_SMALL, ORACLE_TINY, RlweParams
from inspiring_oracle.ring import mul as ring_mul
from inspiring_oracle.ring import sub as ring_sub
from inspiring_oracle.wide_helpers import (
    build_wide_ciphertext,
    decrypt_wide,
    extract_wide_noise,
)

PARAMS = [ORACLE_TINY, ORACLE_SMALL]


# ---------------------------------------------------------------------------
# Wider-secret builders that match the s_hat structure SPEC.md section 4
# ---------------------------------------------------------------------------


def _s_hat_left(s_tilde: list[int], params: RlweParams) -> list[list[int]]:
    """Build ``s_hat_left[j] = tau_g^j(s_tilde)`` for j in [0, d/2)."""
    d, q = params.d, params.q
    two_d = 2 * d
    return [tau(s_tilde, pow(G, j, two_d), q) for j in range(d // 2)]


def _s_hat_right(s_tilde: list[int], params: RlweParams) -> list[list[int]]:
    """Build ``s_hat_right[j] = tau_h(tau_g^j(s_tilde))`` for j in [0, d/2)."""
    d, q = params.d, params.q
    two_d = 2 * d
    h_d = h(d)
    return [
        tau(s_tilde, (pow(G, j, two_d) * h_d) % two_d, q)
        for j in range(d // 2)
    ]


def _per_step_ks_noise_sigma_bound(params: RlweParams) -> float:
    """Theorem 2 single-step bound: sigma_one_ks <= sqrt(ell * d / 4) * z * sigma."""
    return math.sqrt(params.ell * params.d / 4.0) * params.z * params.sigma


def _half_noise_sigma_bound(params: RlweParams) -> float:
    """Theorem 2 half-collapse bound: sqrt(d/2 - 1) * sigma_one_ks (variance adds)."""
    return math.sqrt(params.d // 2 - 1) * _per_step_ks_noise_sigma_bound(params)


# ---------------------------------------------------------------------------
# 1. Input validation
# ---------------------------------------------------------------------------


class TestInputValidation:
    @pytest.mark.parametrize("params", PARAMS)
    def test_rejects_wrong_a_length(self, params: RlweParams) -> None:
        rng = random.Random(0xC11_01)
        s_tilde = rlwe.keygen(params, rng)
        K_g = key_switching.setup(
            tau(s_tilde, G, params.q), s_tilde, params, rng
        )
        a = [[0] * params.d for _ in range(params.d // 2 + 1)]  # one too long
        b = [0] * params.d
        with pytest.raises(ValueError, match="d/2"):
            collapse_half(a, b, K_g, "identity", params)

    @pytest.mark.parametrize("params", PARAMS)
    def test_rejects_invalid_rho(self, params: RlweParams) -> None:
        rng = random.Random(0xC11_02)
        s_tilde = rlwe.keygen(params, rng)
        K_g = key_switching.setup(
            tau(s_tilde, G, params.q), s_tilde, params, rng
        )
        a = [[0] * params.d for _ in range(params.d // 2)]
        b = [0] * params.d
        with pytest.raises(ValueError, match="rho"):
            collapse_half(a, b, K_g, "invalid", params)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# 2. Output shape
# ---------------------------------------------------------------------------


class TestOutputShape:
    @pytest.mark.parametrize("params", PARAMS)
    @pytest.mark.parametrize("rho", ["identity", "tau_h"])
    def test_output_is_single_polynomial_pair(
        self, params: RlweParams, rho: str
    ) -> None:
        rng = random.Random(0xC11_10)
        s_tilde = rlwe.keygen(params, rng)
        K_g = key_switching.setup(
            tau(s_tilde, G, params.q), s_tilde, params, rng
        )
        s_wide = (
            _s_hat_left(s_tilde, params)
            if rho == "identity"
            else _s_hat_right(s_tilde, params)
        )
        m_bar = [0] * params.d
        a, b, _ = build_wide_ciphertext(s_wide, m_bar, params, rng)
        a_out, b_out = collapse_half(a, b, K_g, rho, params)  # type: ignore[arg-type]
        assert isinstance(a_out, list)
        assert len(a_out) == params.d
        assert all(0 <= c < params.q for c in a_out)
        assert isinstance(b_out, list)
        assert len(b_out) == params.d
        assert all(0 <= c < params.q for c in b_out)


# ---------------------------------------------------------------------------
# 3. Correctness firewall: LEFT half (rho = identity)
# ---------------------------------------------------------------------------


class TestCorrectnessLeftHalf:
    """Build under s_hat_left, collapse with rho=identity, decrypt under s_tilde."""

    @pytest.mark.parametrize("params", PARAMS)
    def test_random_messages(self, params: RlweParams) -> None:
        rng = random.Random(0xC11_20)
        s_tilde = rlwe.keygen(params, rng)
        s_hat_left = _s_hat_left(s_tilde, params)
        K_g = key_switching.setup(
            tau(s_tilde, G, params.q), s_tilde, params, rng
        )
        for _ in range(20):
            m_bar = [rng.randrange(params.p) for _ in range(params.d)]
            a, b, _ = build_wide_ciphertext(s_hat_left, m_bar, params, rng)
            a_out, b_out = collapse_half(a, b, K_g, "identity", params)
            recovered = decrypt_wide([a_out], b_out, [s_tilde], params)
            assert recovered == m_bar


# ---------------------------------------------------------------------------
# 4. Correctness firewall: RIGHT half (rho = tau_h)
# ---------------------------------------------------------------------------


class TestCorrectnessRightHalf:
    """Build under s_hat_right, collapse with rho=tau_h, decrypt under tau_h(s_tilde)."""

    @pytest.mark.parametrize("params", PARAMS)
    def test_random_messages(self, params: RlweParams) -> None:
        rng = random.Random(0xC11_30)
        s_tilde = rlwe.keygen(params, rng)
        s_hat_right = _s_hat_right(s_tilde, params)
        K_g = key_switching.setup(
            tau(s_tilde, G, params.q), s_tilde, params, rng
        )
        s_target = tau(s_tilde, h(params.d), params.q)
        for _ in range(20):
            m_bar = [rng.randrange(params.p) for _ in range(params.d)]
            a, b, _ = build_wide_ciphertext(s_hat_right, m_bar, params, rng)
            a_out, b_out = collapse_half(a, b, K_g, "tau_h", params)
            recovered = decrypt_wide([a_out], b_out, [s_target], params)
            assert recovered == m_bar


# ---------------------------------------------------------------------------
# 5. Switch call count
# ---------------------------------------------------------------------------


class TestSwitchCallCount:
    """Each collapse_half does exactly d/2 - 1 KS.Switch calls."""

    @pytest.mark.parametrize("params", PARAMS)
    @pytest.mark.parametrize("rho", ["identity", "tau_h"])
    def test_exact_count(self, params: RlweParams, rho: str) -> None:
        rng = random.Random(0xC11_40)
        s_tilde = rlwe.keygen(params, rng)
        K_g = key_switching.setup(
            tau(s_tilde, G, params.q), s_tilde, params, rng
        )
        s_wide = (
            _s_hat_left(s_tilde, params)
            if rho == "identity"
            else _s_hat_right(s_tilde, params)
        )
        a, b, _ = build_wide_ciphertext(s_wide, [0] * params.d, params, rng)
        key_switching.reset_switch_counter()
        collapse_half(a, b, K_g, rho, params)  # type: ignore[arg-type]
        assert key_switching.switch_call_count() == params.d // 2 - 1


# ---------------------------------------------------------------------------
# 6. Random-component invariant
# ---------------------------------------------------------------------------


class TestRandomComponentInvariant:
    """a_out depends only on (a, K_g, rho), not on b."""

    @pytest.mark.parametrize("params", PARAMS)
    @pytest.mark.parametrize("rho", ["identity", "tau_h"])
    def test_varying_b_does_not_change_a_out(
        self, params: RlweParams, rho: str
    ) -> None:
        rng = random.Random(0xC11_50)
        s_tilde = rlwe.keygen(params, rng)
        K_g = key_switching.setup(
            tau(s_tilde, G, params.q), s_tilde, params, rng
        )
        # Synthetic a (any uniform values are fine -- we don't decrypt here).
        a = [
            [rng.randrange(params.q) for _ in range(params.d)]
            for _ in range(params.d // 2)
        ]
        b1 = [rng.randrange(params.q) for _ in range(params.d)]
        b2 = [rng.randrange(params.q) for _ in range(params.d)]
        a_out_1, _ = collapse_half(a, b1, K_g, rho, params)  # type: ignore[arg-type]
        a_out_2, _ = collapse_half(a, b2, K_g, rho, params)  # type: ignore[arg-type]
        assert a_out_1 == a_out_2


# ---------------------------------------------------------------------------
# 7. Noise growth bounded
# ---------------------------------------------------------------------------


class TestNoiseGrowthBounded:
    """After the full half-collapse, residual noise stays under Delta / 2."""

    @pytest.mark.parametrize("params", PARAMS)
    @pytest.mark.parametrize("rho", ["identity", "tau_h"])
    def test_noise_under_decryption_budget(
        self, params: RlweParams, rho: str
    ) -> None:
        rng = random.Random(0xC11_60)
        s_tilde = rlwe.keygen(params, rng)
        K_g = key_switching.setup(
            tau(s_tilde, G, params.q), s_tilde, params, rng
        )
        s_wide = (
            _s_hat_left(s_tilde, params)
            if rho == "identity"
            else _s_hat_right(s_tilde, params)
        )
        s_target = (
            s_tilde if rho == "identity"
            else tau(s_tilde, h(params.d), params.q)
        )
        budget = params.delta // 2
        for _ in range(20):
            m_bar = [rng.randrange(params.p) for _ in range(params.d)]
            a, b, _ = build_wide_ciphertext(s_wide, m_bar, params, rng)
            a_out, b_out = collapse_half(a, b, K_g, rho, params)  # type: ignore[arg-type]
            e = extract_wide_noise([a_out], b_out, [s_target], m_bar, params)
            max_abs = max(abs(c) for c in e)
            assert max_abs < budget, (
                f"noise overflowed Delta/2 = {budget}: max |e| = {max_abs}"
            )

    @pytest.mark.parametrize("params", PARAMS)
    def test_ks_noise_below_subgaussian_bound(self, params: RlweParams) -> None:
        """With zero input noise the residual IS the half-collapse KS noise.

        Bounded by 6 * sigma_half (subgaussian, variance adds across
        d/2 - 1 independent KS steps).
        """
        rng = random.Random(0xC11_61)
        s_tilde = rlwe.keygen(params, rng)
        s_hat_left = _s_hat_left(s_tilde, params)
        K_g = key_switching.setup(
            tau(s_tilde, G, params.q), s_tilde, params, rng
        )
        # Build a zero-noise wider ciphertext by hand.
        d, q = params.d, params.q
        a = [
            [rng.randrange(q) for _ in range(d)]
            for _ in range(d // 2)
        ]
        b = [0] * d
        for i in range(d // 2):
            b = ring_sub(b, ring_mul(a[i], s_hat_left[i], q), q)
        a_out, b_out = collapse_half(a, b, K_g, "identity", params)
        e = extract_wide_noise([a_out], b_out, [s_tilde], [0] * d, params)
        max_abs = max(abs(c) for c in e)
        bound = 6 * _half_noise_sigma_bound(params)
        assert max_abs < bound, (
            f"max |e_half| = {max_abs} > 6 * sigma_half = {bound:.0f}"
        )


# ---------------------------------------------------------------------------
# 8. Determinism
# ---------------------------------------------------------------------------


class TestDeterminism:
    @pytest.mark.parametrize("params", PARAMS)
    @pytest.mark.parametrize("rho", ["identity", "tau_h"])
    def test_same_inputs_same_output(
        self, params: RlweParams, rho: str
    ) -> None:
        rng = random.Random(0xC11_70)
        s_tilde = rlwe.keygen(params, rng)
        K_g = key_switching.setup(
            tau(s_tilde, G, params.q), s_tilde, params, rng
        )
        s_wide = (
            _s_hat_left(s_tilde, params)
            if rho == "identity"
            else _s_hat_right(s_tilde, params)
        )
        m_bar = [rng.randrange(params.p) for _ in range(params.d)]
        a, b, _ = build_wide_ciphertext(s_wide, m_bar, params, rng)
        out1 = collapse_half(a, b, K_g, rho, params)  # type: ignore[arg-type]
        out2 = collapse_half(a, b, K_g, rho, params)  # type: ignore[arg-type]
        assert out1 == out2


# ---------------------------------------------------------------------------
# 10. Edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    @pytest.mark.parametrize("params", PARAMS)
    @pytest.mark.parametrize("rho", ["identity", "tau_h"])
    def test_all_zero_message(self, params: RlweParams, rho: str) -> None:
        rng = random.Random(0xC11_80)
        s_tilde = rlwe.keygen(params, rng)
        K_g = key_switching.setup(
            tau(s_tilde, G, params.q), s_tilde, params, rng
        )
        s_wide = (
            _s_hat_left(s_tilde, params)
            if rho == "identity"
            else _s_hat_right(s_tilde, params)
        )
        s_target = (
            s_tilde if rho == "identity"
            else tau(s_tilde, h(params.d), params.q)
        )
        m_bar = [0] * params.d
        a, b, _ = build_wide_ciphertext(s_wide, m_bar, params, rng)
        a_out, b_out = collapse_half(a, b, K_g, rho, params)  # type: ignore[arg-type]
        assert decrypt_wide([a_out], b_out, [s_target], params) == m_bar

    @pytest.mark.parametrize("params", PARAMS)
    @pytest.mark.parametrize("rho", ["identity", "tau_h"])
    def test_all_max_message(self, params: RlweParams, rho: str) -> None:
        rng = random.Random(0xC11_81)
        s_tilde = rlwe.keygen(params, rng)
        K_g = key_switching.setup(
            tau(s_tilde, G, params.q), s_tilde, params, rng
        )
        s_wide = (
            _s_hat_left(s_tilde, params)
            if rho == "identity"
            else _s_hat_right(s_tilde, params)
        )
        s_target = (
            s_tilde if rho == "identity"
            else tau(s_tilde, h(params.d), params.q)
        )
        m_bar = [params.p - 1] * params.d
        a, b, _ = build_wide_ciphertext(s_wide, m_bar, params, rng)
        a_out, b_out = collapse_half(a, b, K_g, rho, params)  # type: ignore[arg-type]
        assert decrypt_wide([a_out], b_out, [s_target], params) == m_bar
