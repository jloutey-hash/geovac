"""Backing test for Paper 56 `thm:endo_rigidity` (the PS-4 endomorphism-rigidity
theorem) and its closed-form dimension.

Why this file exists (2026-08-22 FULL certifying run, code dimension). The
theorem's closed form was WRONG and had survived because PS-4 had no test at
all -- 872 of the paper's headline 5,864 residuals rested on a `debug/` driver
that `tests/` never exercised. The stated polynomial was

    n(n+1)(n+2)(3n+5)/24   ->   2, 11, 35, 85, 175      (WRONG)

while the sum it claims to equal, and the values quoted alongside it, are

    sum_{1<=j<=i<=n} (i+1)(j+1)  ->  4, 19, 55, 125, 245  (CORRECT)

The correct factor is (3n+13), not (3n+5). A three-line test would have caught
it on day one, which is exactly the point of the claim->artifact rule.
"""
from __future__ import annotations

import pytest


def endo_dim_sum(n_max: int) -> int:
    """dim End_compat(O_{n_max}) as the paper's block-lower-triangular sum.

    Sector (n, l) with 0 <= l <= n contributes n+1 basis vectors, and a
    compatible endomorphism may map sector i into sector j only for j <= i,
    giving one (i+1)x(j+1) block per ordered pair.
    """
    return sum((i + 1) * (j + 1)
               for i in range(1, n_max + 1)
               for j in range(1, i + 1))


def endo_dim_closed_form(n_max: int) -> int:
    """The corrected closed form printed in `thm:endo_rigidity`."""
    return n_max * (n_max + 1) * (n_max + 2) * (3 * n_max + 13) // 24


@pytest.mark.parametrize("n_max", range(1, 13))
def test_closed_form_equals_the_sum(n_max):
    """The published closed form must reproduce the sum at every cutoff."""
    assert endo_dim_closed_form(n_max) == endo_dim_sum(n_max)


def test_published_values():
    """The values quoted in `thm:endo_rigidity` for n_max = 2, 3, 4, 5."""
    assert [endo_dim_sum(n) for n in (2, 3, 4, 5)] == [19, 55, 125, 245]
    assert [endo_dim_closed_form(n) for n in (2, 3, 4, 5)] == [19, 55, 125, 245]


def test_the_retired_polynomial_is_rejected():
    """Tripwire: the pre-2026-08-22 factor (3n+5) must NOT reproduce the sum.

    Guards against a silent revert. If this ever passes, the closed form has
    been changed back to the wrong one.
    """
    def retired(n):
        return n * (n + 1) * (n + 2) * (3 * n + 5) // 24

    assert retired(2) != endo_dim_sum(2)
    assert [retired(n) for n in (1, 2, 3, 4, 5)] == [2, 11, 35, 85, 175]


def test_block_lower_triangular_structure():
    """The sum is the block-lower-triangular count, not the full endomorphism
    algebra: it must be strictly smaller than dim End_Q(O_{n_max}) for n>=2."""
    for n_max in range(2, 8):
        total_dim = sum(i + 1 for i in range(1, n_max + 1))
        assert endo_dim_sum(n_max) < total_dim ** 2
        # and strictly larger than the block-DIAGONAL subalgebra
        block_diag = sum((i + 1) ** 2 for i in range(1, n_max + 1))
        assert endo_dim_sum(n_max) > block_diag


def test_closed_form_is_a_symbolic_identity():
    """The closed form is provable, not merely verified at sampled cutoffs.

    Added 2026-08-22 (DELTA run) when `rem:tc2c_closed_form` dropped its
    "Empirically across n_max in {1,2,3,4}" framing: the paper now asserts the
    identity, so the backing must prove it rather than sample it.

        sum_{1<=j<=i<=n} (i+1)(j+1)
          = (1/2) * sum_{i=1..n} i(i+1)(i+3)
          = n(n+1)(3n^2 + 19n + 26)/24
          = n(n+1)(n+2)(3n+13)/24,     since 3n^2+19n+26 = (n+2)(3n+13).
    """
    sp = pytest.importorskip("sympy")
    i, j, n = sp.symbols("i j n", positive=True, integer=True)

    # Inner sum over j = 1..i, then outer over i = 1..n, in closed form.
    inner = sp.summation((i + 1) * (j + 1), (j, 1, i))
    total = sp.simplify(sp.summation(inner, (i, 1, n)))

    closed = n * (n + 1) * (n + 2) * (3 * n + 13) / 24
    assert sp.simplify(sp.expand(total - closed)) == 0

    # and the factorisation the derivation turns on
    assert sp.factor(3 * n ** 2 + 19 * n + 26) == (n + 2) * (3 * n + 13)
