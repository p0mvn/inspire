"""Gadget vector ``g_z`` and signed decomposition ``g_z^{-1}`` (SPEC.md section 1).

A **gadget** is a pair ``(g_z, g_z^{-1})`` that lets you take any element
``x in Z_q`` and split it into ``ell`` "small" digits ``d_0, ..., d_{ell-1}``
in base ``z`` such that

    sum_i z^i * d_i  ==  x  (mod q)            <-- reconstruction identity
    |d_i|  <=  z / 2  for every i              <-- size bound

The gadget vector is just the powers of ``z``::

    g_z = [1, z, z^2, ..., z^{ell-1}]   in Z_q^ell

so the reconstruction is literally an inner product: ``<g_z^{-1}(x), g_z> == x``.

Why this matters in InspiRING: every key-switching step (Stage 7's
``KS.Switch``, used ``d - 1`` times in Stage 12's ``collapse``) computes

    (a', b') = (0, b) + g_z^{-1}(a) . K

where ``K = [w_i, y_i]`` is a length-``ell`` vector of RLWE encryptions of
``s_in * z^i``. The gadget identity makes the message survive
(``g_z^{-1}(a) . g_z == a``, so the ``s_in * a`` term reappears and
cancels), while the noise picks up a factor of ``g_z^{-1}(a) . e`` --
which is **small** precisely because every digit is bounded by ``z/2``.

**Signed vs. unsigned digits**: The decomposition returns digits in
``[-z/2, z/2)`` (signed), not ``[0, z)`` (unsigned). This is the source
of the ``z^2 / 4`` term in Theorem 2's noise bound (SPEC.md section 7):

* Signed: ``max |digit|^2 = (z/2)^2 = z^2 / 4``.
* Unsigned: ``max |digit|^2 = (z-1)^2 ~ z^2``.

So switching from unsigned to signed digits cuts the per-step KS noise
worst-case bound by ``4x``. Free win, applied universally in modern HE.

The decomposition algorithm is "balanced base-z": at each step,

* Compute the unsigned digit ``d = x mod z``.
* If ``d >= z/2``, rebalance: ``d -= z`` (now in ``[-z/2, 0)``) and
  carry ``+1`` to the next step (equivalent to ``x += z`` before the
  ``//= z``).
* Advance ``x //= z``.

The carry is what makes the rebalancing work: ``(x + z) // z == x // z + 1``.
"""

from __future__ import annotations

from inspiring_oracle.params import RlweParams


def gz(params: RlweParams) -> list[int]:
    """The gadget vector ``g_z = [1, z, z^2, ..., z^{ell-1}] mod q``.

    Length ``params.ell``; entries in ``[0, q)``.
    """
    return [pow(params.z, i, params.q) for i in range(params.ell)]


def gz_inv_scalar(x: int, params: RlweParams) -> list[int]:
    """Signed base-``z`` decomposition of one scalar ``x`` mod ``q``.

    Returns a length-``ell`` list of digits, each in ``[-z/2, z/2)``,
    such that ``sum_i z^i * digits[i] == x mod q``.

    Negative or ``>= q`` inputs are reduced first via ``x % q``.
    """
    z, ell, q = params.z, params.ell, params.q
    half = z // 2
    x = x % q
    digits: list[int] = []
    for _ in range(ell):
        d = x % z
        if d >= half:
            d -= z
            x += z
        digits.append(d)
        x //= z
    return digits


def gz_inv_poly(p: list[int], params: RlweParams) -> list[list[int]]:
    """Signed base-``z`` decomposition extended **coefficient-wise** to ``R_q``.

    Returns ``ell`` polynomials of length ``len(p)``, where ``result[i][k]``
    is the ``i``-th digit of ``p[k]``. Reconstruction:

        p == sum_i (z^i * result[i])    (in R_q, ring addition + scalar mul)

    Or coefficient-wise: ``p[k] == sum_i z^i * result[i][k] mod q``.

    The "rotated" indexing -- ``ell`` outer, ``len(p)`` inner -- matches
    the way ``KS.Switch`` consumes the digits: it does ``ell`` ring
    multiplications, one per gadget level, against the ``ell`` rows of
    the key-switching matrix ``K``.
    """
    return [[gz_inv_scalar(c, params)[i] for c in p] for i in range(params.ell)]
