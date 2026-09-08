"""Guards for the grid-free (exact) Slater/secular evaluator in
``debug/sturmian_exact_slater.py``.

Why this file exists
--------------------
``geovac/sturmian_secular.py`` evaluates Paper 60's secular matrix by trapezoid
quadrature on a FIXED box ``R_MAX = 60`` and L2-normalises each hydrogenic
orbital ON that truncated box.  Configurations carry weighted charge
``Q_nu = 1/R_nu`` and ``a = Q_nu/n``, so a 10s-14s orbital at ``Q ~ 1`` reaches
~100-200 bohr: the truncation error is systematic IN K, which is the axis
Paper 60 fits its 1-norm exponents along.  These guards pin the exact
(box-free, quadrature-free) evaluator that removes it.

Each test states the WRONG ANSWER it rejects, per the guard-writing rule.
"""
from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("mpmath")
import mpmath as mp  # noqa: E402

# The debug module raises mpmath's GLOBAL working precision at import (it has to:
# see test_module_applies_its_precision_at_import).  Do not let that leak into the
# rest of a full-suite run -- restore the session default immediately and re-raise
# it only for the duration of each test in this file.
_DPS_BEFORE_IMPORT = mp.mp.dps
import debug.sturmian_exact_slater as ex  # noqa: E402

mp.mp.dps = _DPS_BEFORE_IMPORT


@pytest.fixture(autouse=True)
def _high_precision():
    """Run every test in this file at dps=60 and restore the session default after."""
    saved = mp.mp.dps
    ex.set_dps(60)
    yield
    mp.mp.dps = saved


# ======================================================================================
# 1. Equal-exponent case: must reproduce the INDEPENDENT high-precision route in
#    geovac/hypergeometric_slater.py (T-kernel double loop) exactly.
#    Rejects: an algebra error in the suffix-sum (A, W1, W2) reformulation, which
#    would give a wrong-but-smooth answer that no self-consistency check would see.
# ======================================================================================
_UNIT_CASES = [
    (1, 0, 1, 0, 1, 0, 1, 0, 0),
    (2, 1, 2, 1, 2, 1, 2, 1, 2),
    (3, 0, 2, 0, 3, 1, 2, 1, 1),
    (4, 2, 3, 1, 5, 2, 4, 3, 1),
    (6, 1, 4, 3, 5, 2, 6, 0, 2),
    (7, 3, 7, 3, 7, 3, 7, 3, 4),
    (8, 0, 7, 1, 6, 2, 5, 3, 1),
    (9, 4, 8, 2, 9, 3, 8, 3, 5),
]


@pytest.mark.parametrize("case", _UNIT_CASES)
def test_unit_exponent_matches_hypergeometric_slater(case):
    """R^k at a = 1/n reproduces hypergeometric_slater's mpmath T-kernel route."""
    n1, l1, n3, l3, n2, l2, n4, l4, k = case
    from geovac.hypergeometric_slater import _compute_rk_mpmath

    ref = _compute_rk_mpmath(*case, dps=80)          # returns a Python float
    mine = float(ex.slater_rk(ex.orb_unit(n1, l1), ex.orb_unit(n3, l3),
                              ex.orb_unit(n2, l2), ex.orb_unit(n4, l4), k))
    assert ref != 0.0
    assert abs(mine - ref) / abs(ref) < 1e-13, (case, mine, ref)


def test_unit_exponent_F0_1s1s_is_five_eighths():
    """The one closed form the paper quotes: F^0(1s,1s;1s,1s) = 5/8, bit-exact."""
    v = ex.slater_rk(ex.orb_unit(1, 0), ex.orb_unit(1, 0),
                     ex.orb_unit(1, 0), ex.orb_unit(1, 0), 0)
    assert abs(float(v) - 0.625) < 1e-15


def test_mixed_exponent_is_not_the_equal_exponent_answer():
    """Sanity: the mixed-charge integral genuinely differs from the unit-exponent one.

    Rejects a silent fallback in which orb_config() collapsed to orb_unit() and the
    'mixed-exponent generalisation' was really evaluating the a = 1/n case.
    """
    o_u = ex.orb_unit(2, 0)
    o_c = ex.orb_config(2, 0, 2, 3)     # a = 3/sqrt(13), not 1/2
    assert float(ex.rate(o_c)) != pytest.approx(float(ex.rate(o_u)), rel=1e-6)
    a = float(ex.slater_rk(o_u, o_u, o_u, o_u, 0))
    b = float(ex.slater_rk(o_c, o_c, o_c, o_c, 0))
    assert abs(a - b) / abs(a) > 0.1


# ======================================================================================
# 2. Precision.  The Laguerre expansion cancels ~17 digits at n = 14; mpmath's
#    DEFAULT dps is 15, which is float64-equivalent.
#    Rejects: a module that records _DPS but never applies it (the bug that made the
#    K = 244..340 ladder rungs blow up by 9 orders before assert_precision existed).
# ======================================================================================
def test_module_applies_its_precision_at_import():
    """Importing the module must RAISE mpmath's working precision, not just record it.

    Rejects exactly the bug found while building this: the module held _DPS = 60 but
    never executed ``mp.mp.dps = _DPS``, so it silently ran at mpmath's default
    dps=15 (float64-equivalent) and the K = 244..340 ladder rungs came out 9 orders
    too large.
    """
    import importlib

    saved = mp.mp.dps
    try:
        mp.mp.dps = 15
        importlib.reload(ex)
        assert mp.mp.dps >= 50, "module import did not raise mpmath precision"
    finally:
        mp.mp.dps = saved
        ex.set_dps(60)


def test_assert_precision_actually_fires():
    saved = mp.mp.dps
    try:
        mp.mp.dps = 15
        with pytest.raises(RuntimeError):
            ex.assert_precision(50)
    finally:
        mp.mp.dps = saved
        ex.set_dps(60)


def test_result_is_dps_independent_at_n14():
    """A high-n mixed-charge quartet is stable from dps=50 to dps=150.

    Rejects: a working precision too low for the cancellation, which would make
    every exponent in this file precision-dependent (the Levin/log failure mode).
    """
    o1 = ex.orb_config(11, 3, 11, 14)
    o3 = ex.orb_config(14, 3, 11, 14)
    o2 = ex.orb_config(13, 2, 13, 13)
    o4 = ex.orb_config(13, 2, 13, 13)
    try:
        ex.set_dps(150)
        ref = ex.slater_rk(o1, o3, o2, o4, 2)
        ex.set_dps(50)
        lo = ex.slater_rk(o1, o3, o2, o4, 2)
        assert abs(lo - ref) / abs(ref) < mp.mpf(10) ** -25
        ex.set_dps(15)
        bad = ex.slater_rk(o1, o3, o2, o4, 2)
        # and confirm the cancellation is real: dps=15 is NOT good enough
        assert abs(bad - ref) / abs(ref) > mp.mpf(10) ** -13
    finally:
        ex.set_dps(60)


# ======================================================================================
# 3. The grid engine converges TO the exact evaluator: O(dr^2) at a converged box.
#    Rejects: an exact evaluator that is merely self-consistent but disagrees with
#    the (independently written, quadrature-based) production builder.
# ======================================================================================
def _regrid(rmax, ngrid):
    import geovac.sturmian_secular as ss
    ss.R_MAX, ss.N_GRID = rmax, ngrid
    ss.r = np.linspace(1e-7, rmax, ngrid)
    ss.dr = ss.r[1] - ss.r[0]
    ss.r2 = ss.r * ss.r
    ss.reset_caches()
    return ss


def test_grid_converges_to_exact_at_second_order():
    """||M||_1(grid) -> ||M||_1(exact) as dr^2 at a box large enough for n=3."""
    import geovac.sturmian_secular as ss

    tuples = ex.gen_configs(1, {0: 3, 1: 3})            # K = 9
    ref = float(np.abs(ex.build_exact_M(ex.build_exact_configs(tuples))).sum())
    errs = []
    saved = (ss.R_MAX, ss.N_GRID)
    try:
        for ngrid in (36000, 72000, 144000):
            s = _regrid(120.0, ngrid)
            M = s.build_M(s.build_configs(tuples))
            errs.append(abs(float(np.abs(M).sum()) - ref) / ref)
    finally:
        _regrid(*saved)
    assert errs[0] > errs[1] > errs[2]
    for a, b in zip(errs, errs[1:]):
        assert 3.5 < a / b < 4.5, f"not O(dr^2): {errs}"
    assert errs[-1] < 1e-8


def test_production_box_error_grows_with_K():
    """The R_MAX = 60 box error is SYSTEMATIC IN K, not a fixed offset.

    This is the defect the exact evaluator exists to remove: at K = 52 the
    production box misses ||M||_1 by ~1e-3 relative and the miss keeps growing,
    so any exponent fitted along K inherits a K-dependent bias.

    Rejects: 'the grid is converged because E(1s^2) is stable' -- E(1s^2) is a
    K = 1 quantity and is bit-identical at every box size.
    """
    import geovac.sturmian_secular as ss

    saved = (ss.R_MAX, ss.N_GRID)
    rels = []
    try:
        for (lmax, n) in ((1, 3), (2, 5), (3, 6)):
            tuples = ex.gen_configs(lmax, {l: n for l in range(lmax + 1)})
            ref = float(np.abs(ex.build_exact_M(ex.build_exact_configs(tuples))).sum())
            s = _regrid(60.0, 18000)
            M = s.build_M(s.build_configs(tuples))
            rels.append(abs(float(np.abs(M).sum()) - ref) / ref)
    finally:
        _regrid(*saved)
    assert rels[0] < rels[1] < rels[2], rels
    assert rels[0] < 1e-5 and rels[-1] > 1e-4, rels


# ======================================================================================
# 4. ||T0||_1 = Z sum_nu R_nu is a CLOSED-FORM combinatorial sum (no ERI, no radial
#    integral at all).  Derived in debug/sturmian_T0_closed_form.py:
#
#        sum_{m<=a<=b<=N} sqrt(a^-2 + b^-2) = N ln N + (c1 - H_{m-1}) N + O(log N),
#        c1 = gamma + sqrt2 + ln2 - 2 - ln(1+sqrt2) = -0.1967971792...
#        K   = 2 N^2 - 4 N + 4          (lmax = 3, exactly)
#    ==> ||T0||_1 ~ Z sqrt(2K) ln(K/2):  K^{1/2} TIMES A LOG, not a power law.
#
#    Rejects: reading the fitted "K^0.70" as an exponent.  A genuine power law has a
#    flat local slope; this one falls monotonically toward 1/2.
# ======================================================================================
_C1 = (0.5772156649015328606 + np.sqrt(2) + np.log(2) - 2 - np.log(1 + np.sqrt(2)))


def _A_exact(m: int, N: int) -> float:
    a = np.arange(m, N + 1, dtype=np.float64)
    inv2 = 1.0 / a ** 2
    return float(sum(np.sqrt(inv2[i] + inv2[i:]).sum() for i in range(len(a))))


def _sum_Rnu(lmax: int, N: int) -> float:
    return sum(_A_exact(l + 1, N) for l in range(lmax + 1))


def test_T0_matches_the_configuration_sum_exactly():
    """||T0||_1 from the assembled matrix == Z * the pure combinatorial sum."""
    tuples = ex.gen_configs(3, {l: 6 for l in range(4)})
    cfgs = ex.build_exact_configs(tuples)
    assert abs(2.0 * sum(c.Rnu for c in cfgs) - 2.0 * _sum_Rnu(3, 6)) < 1e-10


def test_config_count_closed_form():
    """K = 2N^2 - 4N + 4 exactly for lmax = 3 (used to invert N(K))."""
    for N in range(4, 40):
        assert len(ex.gen_configs(3, {l: N for l in range(4)})) == 2 * N * N - 4 * N + 4


def test_T0_leading_constant_is_the_closed_form():
    """(A_1(N) - N ln N)/N -> c1 = gamma + sqrt2 + ln2 - 2 - ln(1+sqrt2).

    Rejects a mis-derived constant: the residual must SHRINK toward c1, and the
    nearby wrong candidate gamma + phi - 1 = +0.0444 is excluded outright.
    """
    res = [( _A_exact(1, N) - N * np.log(N)) / N for N in (2000, 8000, 20000)]
    assert all(abs(r - _C1) < 0.02 for r in res), res
    assert abs(res[-1] - _C1) < abs(res[0] - _C1)
    wrong = 0.5772156649 + (1 - np.sqrt(2) + np.log(1 + np.sqrt(2))) - 1
    assert abs(res[-1] - wrong) > 0.2


def test_T0_local_slope_falls_toward_one_half():
    """The 'K^0.70' is 1/2 + O(1/log K), not a converged exponent.

    Rejects: a flat exponent.  Over N = 10 -> 3000 the local slope must fall
    monotonically and get below 0.60, which no power law K^0.70 can do.
    """
    Ns = [10, 14, 30, 80, 200, 500, 1300, 3000]
    Ks = [2 * N * N - 4 * N + 4 for N in Ns]
    S = [_sum_Rnu(3, N) for N in Ns]
    sl = [np.log(S[i + 1] / S[i]) / np.log(Ks[i + 1] / Ks[i]) for i in range(len(Ns) - 1)]
    assert all(a > b for a, b in zip(sl, sl[1:])), sl
    assert sl[0] > 0.69 and sl[-1] < 0.60, sl
    assert sl[-1] > 0.50, sl
