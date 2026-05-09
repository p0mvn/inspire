"""RLWE parameters for the InspiRING.Pack oracle.

Two presets are provided:

* ``ORACLE_TINY`` (``d = 8``) - fast unit tests; finishes in milliseconds.
* ``ORACLE_SMALL`` (``d = 16``) - small-but-non-trivial fixtures; seconds.

Both satisfy Theorem 2's correctness condition (``6 * sigma_pack < delta / 2``)
with substantial slack so the oracle's stage-by-stage tests do not flake on
noise-budget edge cases.

Symbol mapping (SPEC.md section 10):

    d       ring dimension; power of 2
    q       ciphertext modulus; must be odd so ``d^{-1} mod q`` exists
    p       plaintext modulus
    sigma   chi (discrete-Gaussian) standard deviation
    z       gadget base
    ell     gadget length; chosen so ``z**ell >= q``
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class RlweParams:
    """Parameters for the cyclotomic ring R_q = Z_q[X] / (X^d + 1)."""

    d: int
    q: int
    p: int
    sigma: float
    z: int
    ell: int

    def __post_init__(self) -> None:
        if self.d <= 0 or (self.d & (self.d - 1)) != 0:
            raise ValueError(f"d must be a positive power of 2, got d = {self.d}")
        if self.q <= 1:
            raise ValueError(f"q must be > 1, got q = {self.q}")
        if self.q % 2 == 0:
            raise ValueError(
                f"q must be odd so d^-1 mod q exists, got q = {self.q}"
            )
        if self.p <= 1:
            raise ValueError(f"p must be > 1, got p = {self.p}")
        if self.z <= 1:
            raise ValueError(f"z must be > 1, got z = {self.z}")
        if self.ell <= 0:
            raise ValueError(f"ell must be positive, got ell = {self.ell}")
        if self.z ** self.ell < self.q:
            raise ValueError(
                f"gadget too short: z**ell = {self.z**self.ell} < q = {self.q}"
            )
        if self.sigma <= 0:
            raise ValueError(f"sigma must be positive, got sigma = {self.sigma}")

    @property
    def d_inv(self) -> int:
        """The multiplicative inverse of d mod q. Used by Stage 1 (TRANSFORM)."""
        return pow(self.d, -1, self.q)

    @property
    def delta(self) -> int:
        """Plaintext scaling factor floor(q / p). LWE encrypts m as Delta * m."""
        return self.q // self.p

    @property
    def noise_budget_sigma(self) -> float:
        """SPEC.md section 7 (Theorem 2) upper bound on sigma_pack.

        Final noise variance after a full pack is bounded by
        ``ell * d^2 * z^2 * sigma^2 / 4``.
        """
        return math.sqrt(self.ell * self.d ** 2 * self.z ** 2 * self.sigma ** 2 / 4.0)

    def correctness_ok(self, num_sigmas: float = 6.0) -> bool:
        """True if Theorem 2's bound leaves enough headroom for decryption.

        With sub-Gaussian noise, ``num_sigmas = 6`` corresponds to a
        per-coefficient decryption-failure probability around 2e-9.
        """
        return num_sigmas * self.noise_budget_sigma < self.delta / 2.0


ORACLE_TINY = RlweParams(d=8, q=12289, p=4, sigma=3.2, z=8, ell=5)
ORACLE_SMALL = RlweParams(d=16, q=65537, p=4, sigma=3.2, z=16, ell=5)
