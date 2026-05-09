"""
Regression tests for ``geovac.hypergeometric_slater``.
=======================================================

Validates the Coulomb Slater R^k integral evaluators
(``compute_rk_algebraic`` exact-Fraction ground truth and
``compute_rk_float`` fast-float dispatch path) across:

1. Small n (n <= 5): exact match between float and Fraction paths.
2. Large n (n >= 6): fix for the catastrophic-cancellation bug —
   ``compute_rk_float`` was returning wrong answers (e.g. -89 Ha instead
   of +0.04 Ha for R^0(12s,12s,12s,12s)) and crashing for k>=2 same-l-zero
   quartets.  The fix dispatches large-n calls to ``compute_rk_algebraic``
   and casts to float.
3. Gaunt-violating quartets raise an informative ``ValueError`` instead of
   crashing with ``factorial() not defined for negative values``.

The fix is documented in
``debug/precision_catalogue_he_2s_singlet_triplet_memo.md`` §2.
"""

from __future__ import annotations

from fractions import Fraction

import pytest

from geovac.hypergeometric_slater import (
    _FLOAT_PATH_MAX_N,
    compute_rk_algebraic,
    compute_rk_float,
    get_rk_float,
)


# Tolerance for the float vs exact-Fraction comparison.  The exact path is
# converted to float at the comparison point; the only floating-point error
# enters in that single cast.
_REL_TOL = 1e-10


def _rel_err(vf: float, ve_frac: Fraction) -> float:
    ve = float(ve_frac)
    if ve == 0.0 and vf == 0.0:
        return 0.0
    denom = max(abs(ve), 1e-30)
    return abs(vf - ve) / denom


# ----------------------------------------------------------------------
# Small-n regression: float and exact paths must agree, and the fast
# pure-float branch must still be used (no perf regression).
# ----------------------------------------------------------------------


class TestSmallN:
    """At n <= 4 the fast pure-float path is taken bit-identical."""

    @pytest.mark.parametrize("args", [
        (1, 0, 1, 0, 1, 0, 1, 0, 0),      # 1s,1s,1s,1s k=0
        (2, 0, 2, 0, 2, 0, 2, 0, 0),      # 2s,2s,2s,2s k=0
        (1, 0, 2, 0, 1, 0, 2, 0, 0),      # 1s,2s,1s,2s k=0
        (2, 0, 2, 1, 2, 0, 2, 1, 1),      # s,p Coulomb k=1
        (3, 1, 3, 2, 3, 1, 3, 2, 1),      # p,d k=1
        (3, 0, 3, 1, 3, 0, 3, 1, 1),      # 3s,3p k=1
        (4, 0, 4, 0, 4, 0, 4, 0, 0),      # 4s,4s k=0  (boundary)
        (4, 1, 4, 1, 4, 1, 4, 1, 0),      # 4p,4p k=0  (boundary)
        (4, 1, 4, 1, 4, 1, 4, 1, 2),      # 4p,4p k=2  (boundary)
    ])
    def test_float_matches_exact(self, args):
        vf = compute_rk_float(*args)
        ve = compute_rk_algebraic(*args)
        assert _rel_err(vf, ve) < _REL_TOL, (
            f"args={args}: float={vf!r}, exact={float(ve)!r}, "
            f"rel={_rel_err(vf, ve):.2e}"
        )

    def test_threshold_value(self):
        """The dispatch threshold _FLOAT_PATH_MAX_N is 4."""
        # If this changes intentionally, update the regression tests below.
        assert _FLOAT_PATH_MAX_N == 4


# ----------------------------------------------------------------------
# Large-n regression: this is the bug fix.  Pre-fix, compute_rk_float
# returned non-physical negative values (e.g. -89 Ha for n=12 ssss k=0)
# and crashed for k >= 2 same-shell s,s.  After the fix, compute_rk_float
# delegates to the exact Fraction path for max(n) > 5.
# ----------------------------------------------------------------------


class TestLargeN:
    """Same-shell same-l quartets at n >= 6 must now be correct.

    The fast cases (n <= 12) are run unmarked; the slow same-shell n=13
    and n=15 cases are marked ``@pytest.mark.slow`` because each call
    traverses the exact Fraction path which is ~10-30 seconds at those
    cutoffs.
    """

    @pytest.mark.parametrize("n", [6, 7, 8, 10, 12])
    def test_ssss_k0(self, n):
        args = (n, 0, n, 0, n, 0, n, 0, 0)
        vf = compute_rk_float(*args)
        ve = compute_rk_algebraic(*args)
        assert _rel_err(vf, ve) < _REL_TOL, (
            f"n={n} ssss k=0: float={vf!r}, exact={float(ve)!r}"
        )
        # Sanity: the value should be positive and small (R^0(ns,ns,ns,ns)
        # decays as ~1/n with n).  Pre-fix this returned huge negative
        # numbers (-89 at n=12, -5e7 at n=15).
        assert 0.0 < vf < 1.0, (
            f"n={n} ssss k=0: value {vf!r} outside physical range"
        )

    @pytest.mark.slow
    @pytest.mark.parametrize("n", [13, 15])
    def test_ssss_k0_slow(self, n):
        args = (n, 0, n, 0, n, 0, n, 0, 0)
        vf = compute_rk_float(*args)
        ve = compute_rk_algebraic(*args)
        assert _rel_err(vf, ve) < _REL_TOL
        assert 0.0 < vf < 1.0

    @pytest.mark.parametrize("n", [6, 8, 12])
    def test_pppp_k0(self, n):
        args = (n, 1, n, 1, n, 1, n, 1, 0)
        vf = compute_rk_float(*args)
        ve = compute_rk_algebraic(*args)
        assert _rel_err(vf, ve) < _REL_TOL

    @pytest.mark.parametrize("n", [6, 8, 12])
    def test_pppp_k2(self, n):
        """Same-shell p,p k=2 — previously crashed for some patterns."""
        args = (n, 1, n, 1, n, 1, n, 1, 2)
        vf = compute_rk_float(*args)
        ve = compute_rk_algebraic(*args)
        assert _rel_err(vf, ve) < _REL_TOL

    @pytest.mark.slow
    @pytest.mark.parametrize("n", [15])
    def test_pppp_slow(self, n):
        for k in [0, 2]:
            args = (n, 1, n, 1, n, 1, n, 1, k)
            vf = compute_rk_float(*args)
            ve = compute_rk_algebraic(*args)
            assert _rel_err(vf, ve) < _REL_TOL, f"n={n} k={k}"

    @pytest.mark.parametrize("n", [6, 8, 12])
    def test_spsp_k1(self, n):
        """Mixed-l same-shell, k=1."""
        args = (n, 0, n, 1, n, 0, n, 1, 1)
        vf = compute_rk_float(*args)
        ve = compute_rk_algebraic(*args)
        assert _rel_err(vf, ve) < _REL_TOL

    @pytest.mark.parametrize("args", [
        (10, 0, 12, 0, 10, 0, 12, 0, 0),    # mixed-n s
        (8, 1, 12, 1, 8, 1, 12, 1, 0),      # mixed-n p
        (6, 0, 8, 1, 6, 0, 8, 1, 1),        # mixed-n s,p k=1
    ])
    def test_mixed_n(self, args):
        vf = compute_rk_float(*args)
        ve = compute_rk_algebraic(*args)
        assert _rel_err(vf, ve) < _REL_TOL, (
            f"args={args}: float={vf!r}, exact={float(ve)!r}"
        )

    @pytest.mark.slow
    def test_mixed_n_slow(self):
        args = (15, 1, 8, 2, 15, 1, 8, 2, 1)
        vf = compute_rk_float(*args)
        ve = compute_rk_algebraic(*args)
        assert _rel_err(vf, ve) < _REL_TOL


# ----------------------------------------------------------------------
# Gaunt-violating quartets: the formula's outer-integral power m1 = b-k-1
# (or m2 = a-k-1) goes negative when k > l_a + l_b.  This corresponds to a
# request that the angular Gaunt coefficient would zero out anyway.  Pre-
# fix this crashed with ``factorial() not defined for negative values``;
# the fix raises an informative ValueError instead.
# ----------------------------------------------------------------------


class TestGauntViolation:
    """k > l_a + l_b should raise ValueError, not crash."""

    def test_ssss_k2_raises(self):
        with pytest.raises(ValueError) as exc_info:
            compute_rk_float(12, 0, 12, 0, 12, 0, 12, 0, 2)
        # The error message should be informative — mention Gaunt.
        msg = str(exc_info.value).lower()
        assert "gaunt" in msg or "negative" in msg

    def test_ssss_k3_raises(self):
        with pytest.raises(ValueError):
            compute_rk_float(12, 0, 12, 0, 12, 0, 12, 0, 3)


# ----------------------------------------------------------------------
# Cached-API regression: get_rk_float must produce the same answers as
# compute_rk_float (it's just a memoisation wrapper).
# ----------------------------------------------------------------------


class TestCacheConsistency:

    @pytest.mark.parametrize("args", [
        (2, 0, 2, 0, 2, 0, 2, 0, 0),     # small-n fast path
        (8, 1, 8, 1, 8, 1, 8, 1, 0),     # large-n exact path
        (12, 0, 12, 0, 12, 0, 12, 0, 0), # large-n exact path
    ])
    def test_get_matches_compute(self, args):
        v1 = compute_rk_float(*args)
        v2 = get_rk_float(*args)
        # Cached version is identical (same float).
        assert v1 == v2 or abs(v1 - v2) < _REL_TOL * max(abs(v1), 1e-30)

    @pytest.mark.slow
    def test_get_matches_compute_n15(self):
        args = (15, 0, 15, 0, 15, 0, 15, 0, 0)
        v1 = compute_rk_float(*args)
        v2 = get_rk_float(*args)
        assert v1 == v2 or abs(v1 - v2) < _REL_TOL * max(abs(v1), 1e-30)


# ----------------------------------------------------------------------
# Slow regression at n=20 (large quartet, exact path is several seconds).
# Marked slow so the standard fast suite stays quick.
# ----------------------------------------------------------------------


@pytest.mark.slow
class TestSlow:

    def test_n20_ssss(self):
        args = (20, 0, 20, 0, 20, 0, 20, 0, 0)
        vf = compute_rk_float(*args)
        ve = compute_rk_algebraic(*args)
        assert _rel_err(vf, ve) < _REL_TOL
