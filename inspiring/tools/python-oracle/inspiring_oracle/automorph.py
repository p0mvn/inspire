"""Galois automorphisms of R_q = Z_q[X] / (X^d + 1) and the Lemma 1 trace.

The Galois group ``Gal(R_q / Z_q)`` is isomorphic to ``(Z / 2dZ)*`` and, for
``d`` a power of 2, has structure ``Z_{d/2} x Z_2``. The two canonical
generators (SPEC.md section 2) are:

* ``G = 5`` -- generator of the order-``d/2`` cyclic subgroup. Verified by
  ``5^{d/2} = 1 mod 2d`` and no smaller positive power of 5 reaches 1.
* ``h(d) = 2d - 1`` (which is ``-1 mod 2d``) -- the order-2 generator.

Every element ``g`` of ``Z*_{2d}`` -- equivalently, every odd integer in
``[1, 2d)`` -- can be written uniquely as ``5^j * h^b mod 2d`` for some
``0 <= j < d/2`` and ``b in {0, 1}``. There are therefore ``d`` distinct
Galois automorphisms; the trace operator (Stage 3) sums precisely these.

How ``tau_g`` acts on a polynomial: every monomial ``X^i`` is sent to
``X^{i*g mod 2d}``. Because ``X^{2d} = 1`` in R_q (a consequence of
``X^d = -1``), exponents reduce mod ``2d``. When the reduced exponent lands
in ``[d, 2d)``, it folds back to ``[0, d)`` with a sign flip via the
negacyclic rule. This single line of reasoning is the algorithmic content
of every function below.

The trace operator ``Tr`` (Stage 3) is the only "non-obvious" function in
this module; everything else is just primitive Galois machinery. ``Tr`` is
the conceptual heart of InspiRING -- a constant-term extractor built from
sums of automorphisms (SPEC.md section 3, Lemma 1).
"""

from __future__ import annotations

from inspiring_oracle.ring import add as _ring_add

G: int = 5
"""SPEC.md section 2: canonical generator of the order-``d/2`` subgroup of Z*_{2d}."""


def h(d: int) -> int:
    """Canonical generator of the order-2 subgroup: ``2d - 1`` (== -1 mod 2d)."""
    return 2 * d - 1


def tau(p: list[int], g: int, q: int) -> list[int]:
    """Apply the Galois automorphism ``tau_g(p)(X) = p(X^g)``.

    Algorithm: for each input coefficient ``p[i]`` (the coefficient of
    ``X^i``), compute the new exponent ``e = (i * g) mod 2d``; if ``e < d``
    it lands at position ``e`` with positive sign, otherwise at position
    ``e - d`` with negative sign (the negacyclic fold).

    ``g`` must be in ``Z*_{2d}``. Because ``2d`` is a power of 2, this is
    equivalent to ``g`` being odd. The function does not enforce this;
    passing an even ``g`` produces a non-bijective output.
    """
    d = len(p)
    two_d = 2 * d
    out = [0] * d
    for i, c in enumerate(p):
        if c == 0:
            continue
        e = (i * g) % two_d
        if e < d:
            out[e] = (out[e] + c) % q
        else:
            out[e - d] = (out[e - d] - c) % q
    return out


def tau_g(p: list[int], q: int) -> list[int]:
    """Apply ``tau_5``, the canonical order-``d/2`` generator."""
    return tau(p, G, q)


def tau_g_pow(p: list[int], j: int, q: int) -> list[int]:
    """Apply ``tau_g^j == tau_{5^j mod 2d}`` as a single ``tau`` call.

    This is mathematically the same as iterating ``tau_g`` ``j`` times but
    cheaper -- one ``tau`` call instead of ``j`` of them.
    """
    d = len(p)
    return tau(p, pow(G, j, 2 * d), q)


def tau_h(p: list[int], q: int) -> list[int]:
    """Apply ``tau_h`` where ``h = 2d - 1 = -1 mod 2d`` (the order-2 generator).

    Concretely: ``tau_h(p)(X) = p(X^{-1})`` interpreted in R_q. This sends
    ``X^k`` to ``-X^{d-k}`` for ``k > 0`` and fixes the constant term.
    """
    d = len(p)
    return tau(p, h(d), q)


# --------------------------------------------------------------------------
# The trace operator (SPEC.md section 3, Lemma 1)
# --------------------------------------------------------------------------


def half_trace(p: list[int], q: int) -> list[int]:
    """Compute ``pi_{d/2}(p) := sum_{j=0}^{d/2 - 1} tau_g^j(p)``.

    By SPEC.md Lemma 5 (with ``gamma = d/2``), the result is

        pi_{d/2}(p) = (d/2) * (c_0 + c_{d/2} * X^{d/2})

    -- i.e. positions ``0`` and ``d/2`` of the output equal ``(d/2)*c_0``
    and ``(d/2)*c_{d/2}`` respectively, and every other position is zero.

    This is the "first half" of the Lemma 1 trace; ``trace`` then folds in
    a single ``tau_h`` application to cancel the ``X^{d/2}`` coefficient.

    Exposed as a public function because the Lemma 5 intermediate is a
    powerful firewall for bugs in ``tau``: if ``half_trace`` produces
    nonzero values at any position other than ``0`` and ``d/2``, the
    monomial-folding logic in ``tau`` is wrong.
    """
    d = len(p)
    two_d = 2 * d
    out = [0] * d
    for j in range(d // 2):
        gj = pow(G, j, two_d)
        out = _ring_add(out, tau(p, gj, q), q)
    return out


def trace(p: list[int], q: int) -> list[int]:
    """The Lemma 1 trace operator: extracts the constant coefficient.

    Defined in SPEC.md section 3 as

        Tr(p) := sum_{j=0}^{d/2 - 1} (tau_g^j(p) + tau_h o tau_g^j(p))

    By Lemma 1 this equals ``d * c_0`` (returned as a constant polynomial:
    position 0 holds ``d * p[0] mod q``, every other position is zero).

    Implementation note: equivalent to ``half_trace(p) + tau_h(half_trace(p))``
    by linearity of ``tau_h`` (it factors out of the inner sum). This form
    uses ``d/2 + 1`` ``tau`` calls instead of ``d`` and matches the proof
    structure (Lemma 5 + the X^{d/2} cancellation step).

    Stage 8 (TRANSFORM) immediately scales by ``d^{-1} mod q`` to recover
    just ``c_0`` -- this is why ``RlweParams`` requires ``q`` to be odd.
    """
    ht = half_trace(p, q)
    return _ring_add(ht, tau_h(ht, q), q)
