"""Tests for ``intermediate.aggregate`` -- Stage 2 of Algorithm 1.

Aggregate fuses ``d`` IRCtxs (each carrying one LWE message in slot 0)
into a single IRCtx whose ``m_hat`` polynomial has one message per
coefficient slot. Like ``transform``, it is a noise-free algebraic
rearrangement; the ``X^k`` shifts route each input's noise to a unique
output slot, so noise is moved but never amplified or summed.

Test groups
-----------

1. **TestInputValidation** -- reject any input list whose length is not
   exactly ``d``.

2. **TestOutputShape** -- output IRCtx has the right shape and all
   coefficients are in canonical ``[0, q)`` form.

3. **TestBTildeIsLiteralBValues** -- the SPEC.md section 5 algebraic
   simplification: for fresh ``transform`` outputs,
   ``b_tilde_agg = [b_0, b_1, ..., b_{d-1}]``. The only piece of stage 2
   that touches per-query data, and a key invariant for the offline /
   online split.

4. **TestDecryptionRoundtrip** -- the firewall: encrypt ``d`` LWEs of
   ``(m_0, ..., m_{d-1})``, transform each, aggregate, decrypt under
   ``s_hat`` -> recover ``[m_0, ..., m_{d-1}]``. Many random samples plus
   targeted edge cases (all-zero messages, single non-zero slot, all
   max-value).

5. **TestNoiseRoutingPerSlot** -- the noise extracted from slot ``k`` of
   the aggregated IRCtx equals the original LWE noise of input ``k``.
   No mixing across slots, no amplification.

6. **TestAHatAggIndependentOfBs** -- for fixed ``a_k``-vectors,
   ``a_hat_agg`` does not change as the ``b_k``-values change. This is
   the foundation of the offline/online split: the server precomputes
   ``a_hat_agg`` once per CRS-fixed ``a``-vectors and only assembles the
   trivial ``b_tilde_agg`` per query.

7. **TestZeroIRCtxsAggregateToZero** -- aggregating ``d`` IRCtxs whose
   inputs are all-zero LWE ciphertexts (``a = 0``, ``b = 0``) yields the
   zero IRCtx.

8. **TestSingleNonzeroSlot** -- aggregating ``d - 1`` zero IRCtxs and one
   non-zero IRCtx at position ``k`` yields an aggregate whose only
   non-zero plaintext slot is ``k``.

9. **TestDeterminism** -- aggregate is a pure function of its inputs.

10. **TestRecoveredPolynomialMatchesPerSlotBetweenSampleAndAgg** --
    cross-check that aggregating ``d`` random IRCtxs produces a
    decryption polynomial whose slot ``k`` coefficient (centered) equals
    ``Delta * m_k + e_k_lwe`` for each ``k``. Pins down both the message
    routing and the noise routing in one pass.
"""

from __future__ import annotations

import random

import pytest

from inspiring_oracle import lwe, rlwe
from inspiring_oracle.decrypt_under_s_hat import (
    decrypt_polynomial_under_s_hat,
    decrypt_under_s_hat,
    extract_noise_under_s_hat,
    s_hat_from_s_tilde,
)
from inspiring_oracle.intermediate import IRCtx, aggregate, transform
from inspiring_oracle.params import ORACLE_SMALL, ORACLE_TINY, RlweParams

PARAMS = [ORACLE_TINY, ORACLE_SMALL]


def _build_zero_irctx(params: RlweParams) -> IRCtx:
    """An IRCtx with all-zero a_hat and all-zero b_tilde."""
    return IRCtx(
        a_hat=[[0] * params.d for _ in range(params.d)],
        b_tilde=[0] * params.d,
    )


def _encrypt_d_messages(
    s: list[int],
    messages: list[int],
    params: RlweParams,
    rng: random.Random,
) -> tuple[list[lwe.LweCiphertext], list[int]]:
    """Encrypt ``messages`` under ``s``, returning ciphertexts and per-msg noise."""
    cts = [lwe.encrypt(s, m, params, rng) for m in messages]
    noises = [
        lwe.extract_noise(s, ct, m, params)
        for ct, m in zip(cts, messages, strict=True)
    ]
    return cts, noises


# ---------------------------------------------------------------------------
# 1. Input validation
# ---------------------------------------------------------------------------


class TestInputValidation:
    @pytest.mark.parametrize("params", PARAMS)
    def test_rejects_too_few(self, params: RlweParams) -> None:
        irctxs = [_build_zero_irctx(params) for _ in range(params.d - 1)]
        with pytest.raises(ValueError, match="exactly d"):
            aggregate(irctxs, params)

    @pytest.mark.parametrize("params", PARAMS)
    def test_rejects_too_many(self, params: RlweParams) -> None:
        irctxs = [_build_zero_irctx(params) for _ in range(params.d + 1)]
        with pytest.raises(ValueError, match="exactly d"):
            aggregate(irctxs, params)

    @pytest.mark.parametrize("params", PARAMS)
    def test_empty_list_raises(self, params: RlweParams) -> None:
        with pytest.raises(ValueError, match="exactly d"):
            aggregate([], params)


# ---------------------------------------------------------------------------
# 2. Output shape
# ---------------------------------------------------------------------------


class TestOutputShape:
    @pytest.mark.parametrize("params", PARAMS)
    def test_shapes_and_ranges(self, params: RlweParams) -> None:
        rng = random.Random(0xA66_01)
        s = lwe.keygen(params, rng)
        irctxs = [
            transform(lwe.encrypt(s, rng.randrange(params.p), params, rng), params)
            for _ in range(params.d)
        ]
        agg = aggregate(irctxs, params)
        assert len(agg.a_hat) == params.d
        for poly in agg.a_hat:
            assert len(poly) == params.d
            assert all(0 <= c < params.q for c in poly)
        assert len(agg.b_tilde) == params.d
        assert all(0 <= c < params.q for c in agg.b_tilde)


# ---------------------------------------------------------------------------
# 3. b_tilde_agg = literal [b_0, ..., b_{d-1}] (SPEC simplification)
# ---------------------------------------------------------------------------


class TestBTildeIsLiteralBValues:
    """For fresh transform inputs, b_tilde_agg[k] == b_k mod q.

    SPEC.md section 5 calls this out explicitly: since each b_tilde_k is
    constant, the X^k shift just places b_k at slot k with no wrap-around
    and no sign flip. The d coefficients of b_tilde_agg are literally the
    LWE b-values -- the only piece of stage 2 that depends on per-query
    data, and the algebraic basis for the cheap online phase.
    """

    @pytest.mark.parametrize("params", PARAMS)
    def test_b_tilde_agg_equals_b_vector(self, params: RlweParams) -> None:
        rng = random.Random(0xA66_02)
        s = lwe.keygen(params, rng)
        cts = [
            lwe.encrypt(s, rng.randrange(params.p), params, rng)
            for _ in range(params.d)
        ]
        irctxs = [transform(ct, params) for ct in cts]
        agg = aggregate(irctxs, params)
        assert agg.b_tilde == [ct.b % params.q for ct in cts]


# ---------------------------------------------------------------------------
# 4. Decryption round-trip (the firewall)
# ---------------------------------------------------------------------------


class TestDecryptionRoundtrip:
    """Encrypt d LWEs of m_0, ..., m_{d-1}; transform each; aggregate;
    decrypt under s_hat; recover [m_0, ..., m_{d-1}].
    """

    @pytest.mark.parametrize("params", PARAMS)
    def test_random_messages(self, params: RlweParams) -> None:
        rng = random.Random(0xA66_10)
        s = lwe.keygen(params, rng)
        s_tilde = rlwe.s_tilde_from_s(s, params)
        s_hat = s_hat_from_s_tilde(s_tilde, params)
        for _ in range(20):
            messages = [rng.randrange(params.p) for _ in range(params.d)]
            cts, _ = _encrypt_d_messages(s, messages, params, rng)
            irctxs = [transform(ct, params) for ct in cts]
            agg = aggregate(irctxs, params)
            recovered = decrypt_under_s_hat(agg, s_hat, params)
            assert recovered == messages

    @pytest.mark.parametrize("params", PARAMS)
    def test_all_zero_messages(self, params: RlweParams) -> None:
        rng = random.Random(0xA66_11)
        s = lwe.keygen(params, rng)
        s_tilde = rlwe.s_tilde_from_s(s, params)
        s_hat = s_hat_from_s_tilde(s_tilde, params)
        messages = [0] * params.d
        cts, _ = _encrypt_d_messages(s, messages, params, rng)
        irctxs = [transform(ct, params) for ct in cts]
        agg = aggregate(irctxs, params)
        assert decrypt_under_s_hat(agg, s_hat, params) == messages

    @pytest.mark.parametrize("params", PARAMS)
    def test_all_max_messages(self, params: RlweParams) -> None:
        rng = random.Random(0xA66_12)
        s = lwe.keygen(params, rng)
        s_tilde = rlwe.s_tilde_from_s(s, params)
        s_hat = s_hat_from_s_tilde(s_tilde, params)
        messages = [params.p - 1] * params.d
        cts, _ = _encrypt_d_messages(s, messages, params, rng)
        irctxs = [transform(ct, params) for ct in cts]
        agg = aggregate(irctxs, params)
        assert decrypt_under_s_hat(agg, s_hat, params) == messages

    @pytest.mark.parametrize("params", PARAMS)
    def test_alternating_messages(self, params: RlweParams) -> None:
        """Force every other slot to be max value -- catches off-by-one
        sign errors in the X^k shift."""
        rng = random.Random(0xA66_13)
        s = lwe.keygen(params, rng)
        s_tilde = rlwe.s_tilde_from_s(s, params)
        s_hat = s_hat_from_s_tilde(s_tilde, params)
        messages = [params.p - 1 if k % 2 == 0 else 0 for k in range(params.d)]
        cts, _ = _encrypt_d_messages(s, messages, params, rng)
        irctxs = [transform(ct, params) for ct in cts]
        agg = aggregate(irctxs, params)
        assert decrypt_under_s_hat(agg, s_hat, params) == messages


# ---------------------------------------------------------------------------
# 5. Noise routing per slot
# ---------------------------------------------------------------------------


class TestNoiseRoutingPerSlot:
    """Slot k of the aggregated noise polynomial equals the LWE noise of
    input k. No cross-slot mixing, no amplification; the X^k shift simply
    routes each input's slot-0 noise to its dedicated output slot.
    """

    @pytest.mark.parametrize("params", PARAMS)
    def test_per_slot_noise_matches_per_input_lwe_noise(
        self, params: RlweParams
    ) -> None:
        rng = random.Random(0xA66_20)
        s = lwe.keygen(params, rng)
        s_tilde = rlwe.s_tilde_from_s(s, params)
        s_hat = s_hat_from_s_tilde(s_tilde, params)
        messages = [rng.randrange(params.p) for _ in range(params.d)]
        cts, lwe_noises = _encrypt_d_messages(s, messages, params, rng)
        irctxs = [transform(ct, params) for ct in cts]
        agg = aggregate(irctxs, params)
        agg_noise = extract_noise_under_s_hat(agg, s_hat, messages, params)
        assert agg_noise == lwe_noises


# ---------------------------------------------------------------------------
# 6. a_hat_agg independent of b's (offline/online split)
# ---------------------------------------------------------------------------


class TestAHatAggIndependentOfBs:
    """For fixed a_k-vectors, varying the b_k-values does not change
    a_hat_agg. The online phase only needs to assemble the trivial
    b_tilde_agg = [b_0, ..., b_{d-1}].
    """

    @pytest.mark.parametrize("params", PARAMS)
    def test_random_b_vectors_do_not_change_a_hat(
        self, params: RlweParams
    ) -> None:
        rng = random.Random(0xA66_30)
        a_vecs = [
            [rng.randrange(params.q) for _ in range(params.d)]
            for _ in range(params.d)
        ]
        b1_vec = [rng.randrange(params.q) for _ in range(params.d)]
        b2_vec = [rng.randrange(params.q) for _ in range(params.d)]
        irctxs1 = [
            transform(lwe.LweCiphertext(a=a, b=b), params)
            for a, b in zip(a_vecs, b1_vec, strict=True)
        ]
        irctxs2 = [
            transform(lwe.LweCiphertext(a=a, b=b), params)
            for a, b in zip(a_vecs, b2_vec, strict=True)
        ]
        agg1 = aggregate(irctxs1, params)
        agg2 = aggregate(irctxs2, params)
        assert agg1.a_hat == agg2.a_hat
        # Sanity: b_tilde_agg DOES change with b's.
        assert agg1.b_tilde == [b % params.q for b in b1_vec]
        assert agg2.b_tilde == [b % params.q for b in b2_vec]


# ---------------------------------------------------------------------------
# 7. Zero IRCtxs aggregate to zero IRCtx
# ---------------------------------------------------------------------------


class TestZeroIRCtxsAggregateToZero:
    @pytest.mark.parametrize("params", PARAMS)
    def test_d_zero_irctxs(self, params: RlweParams) -> None:
        irctxs = [_build_zero_irctx(params) for _ in range(params.d)]
        agg = aggregate(irctxs, params)
        zero_a_hat = [[0] * params.d for _ in range(params.d)]
        assert agg.a_hat == zero_a_hat
        assert agg.b_tilde == [0] * params.d


# ---------------------------------------------------------------------------
# 8. Single non-zero slot
# ---------------------------------------------------------------------------


class TestSingleNonzeroSlot:
    """If only input k has a non-zero LWE message, the aggregate's
    decryption polynomial has a non-zero plaintext only at slot k.
    """

    @pytest.mark.parametrize("params", PARAMS)
    @pytest.mark.parametrize("nonzero_slot", [0, 1, 3])
    def test_only_slot_k_is_nonzero(
        self, params: RlweParams, nonzero_slot: int
    ) -> None:
        nonzero_slot = nonzero_slot % params.d
        rng = random.Random(0xA66_40 + nonzero_slot)
        s = lwe.keygen(params, rng)
        s_tilde = rlwe.s_tilde_from_s(s, params)
        s_hat = s_hat_from_s_tilde(s_tilde, params)
        messages = [0] * params.d
        messages[nonzero_slot] = params.p - 1
        cts, _ = _encrypt_d_messages(s, messages, params, rng)
        irctxs = [transform(ct, params) for ct in cts]
        agg = aggregate(irctxs, params)
        recovered = decrypt_under_s_hat(agg, s_hat, params)
        assert recovered[nonzero_slot] == params.p - 1
        for k in range(params.d):
            if k != nonzero_slot:
                assert recovered[k] == 0


# ---------------------------------------------------------------------------
# 9. Determinism
# ---------------------------------------------------------------------------


class TestDeterminism:
    @pytest.mark.parametrize("params", PARAMS)
    def test_same_inputs_same_output(self, params: RlweParams) -> None:
        rng = random.Random(0xA66_50)
        s = lwe.keygen(params, rng)
        cts = [
            lwe.encrypt(s, rng.randrange(params.p), params, rng)
            for _ in range(params.d)
        ]
        irctxs = [transform(ct, params) for ct in cts]
        a = aggregate(irctxs, params)
        b = aggregate(irctxs, params)
        assert a == b


# ---------------------------------------------------------------------------
# 10. Per-slot raw decryption coefficients match Delta * m_k + e_k
# ---------------------------------------------------------------------------


class TestRecoveredPolynomialMatchesPerSlotMessageAndNoise:
    """The recovered raw polynomial slot k (centered) equals
    ``Delta * m_k + e_k_lwe`` exactly for every k.

    Pins down both the message routing (Delta * m_k at slot k) and the
    noise routing (e_k_lwe at slot k) in a single assertion. This is the
    sharpest available statement about Stage 2's algebraic structure.
    """

    @pytest.mark.parametrize("params", PARAMS)
    def test_centered_per_slot_decomposition(self, params: RlweParams) -> None:
        rng = random.Random(0xA66_60)
        s = lwe.keygen(params, rng)
        s_tilde = rlwe.s_tilde_from_s(s, params)
        s_hat = s_hat_from_s_tilde(s_tilde, params)
        messages = [rng.randrange(params.p) for _ in range(params.d)]
        cts, lwe_noises = _encrypt_d_messages(s, messages, params, rng)
        irctxs = [transform(ct, params) for ct in cts]
        agg = aggregate(irctxs, params)
        raw = decrypt_polynomial_under_s_hat(agg, s_hat, params)
        for k in range(params.d):
            # Compute residual = (raw[k] - Delta * m_k) mod q, then center.
            residual_modq = (raw[k] - params.delta * messages[k]) % params.q
            residual = (
                residual_modq if residual_modq <= params.q // 2
                else residual_modq - params.q
            )
            assert residual == lwe_noises[k], (
                f"slot {k}: residual {residual} != lwe noise {lwe_noises[k]}"
            )
