"""Regression tests for the isoenergetic generalized-Sturmian He secular builder.

Backs four load-bearing Paper 60 (group2) claims, ported from the
``debug/sturmian_he_lmax.py`` driver into ``geovac/sturmian_secular.py``:

  (i)   single-config 1s^2  ->  E = -2.847 Ha  (textbook variational He);
  (ii)  the entrywise 1-norm ``||M||_1`` grows SUBLINEARLY with the config count K,
        approaching ``~K^0.84`` for the full s+p+d+f basis (Paper 60 ``eq:sublinear``);
  (iii) restoring the L2 overlap metric S (generalized eigenproblem ``M B = p S B``)
        is ill-conditioned -- ``cond(S)`` grows from ~4 into the thousands
        (the paper's "4 -> 3673");
  (iv)  ``[eq:secular]`` the interelectron matrix T' is a matrix of PURE NUMBERS,
        independent of the nuclear charge Z.

The K-ladder (ii) and metric-divergence (iii) tests are marked ``slow`` (build_M at
K~100 is ~15 s; the ladder/divergence sweeps aggregate several such builds).
"""
import numpy as np
import pytest

from geovac.sturmian_secular import (
    build_M,
    build_Tprime,
    build_configs,
    gen_configs,
    solve,
    solve_with_metric,
)


# --------------------------------------------------------------------------------------
# Claim (i): single-config 1s^2 = -2.847 Ha (textbook variational He).
# --------------------------------------------------------------------------------------
def test_single_config_1s2() -> None:
    """1s^2 single configuration recovers the textbook variational He energy -2.847 Ha."""
    E, one_norm, K, M = solve([(0, 1, 1)])
    assert K == 1
    assert M.shape == (1, 1)
    assert abs(E - (-2.84766)) < 3e-3, f"E(1s^2) = {E:.5f}, expected -2.84766"


# --------------------------------------------------------------------------------------
# Claim (iv): eq:secular -- the interelectron matrix T' is Z-independent (pure numbers).
# Fast: a small s+p basis (K=9).  T' is built at the reference weighted charge
# Q_nu = 1/R_nu, so no nuclear charge Z enters it; only diag(Z*R_nu) carries Z.
# --------------------------------------------------------------------------------------
def test_Tprime_Z_independence() -> None:
    """T' (interelectron block of M) is bit-identical across two nuclear charges Z."""
    cfgs = build_configs(gen_configs(1, {0: 3, 1: 3}))
    Rnu = np.array([c.Rnu for c in cfgs])

    # Recover T' from build_M at two charges by subtracting the Z-dependent diagonal.
    M_z2 = build_M(cfgs, Z=2.0)
    M_z3 = build_M(cfgs, Z=3.0)
    Tp_from_z2 = M_z2 - np.diag(2.0 * Rnu)
    Tp_from_z3 = M_z3 - np.diag(3.0 * Rnu)

    residual = float(np.max(np.abs(Tp_from_z2 - Tp_from_z3)))
    assert residual < 1e-12, f"T' varies with Z at level {residual:.3e} (should be ~0)"

    # And it equals the standalone build_Tprime (which never sees Z at all).
    Tp_direct = build_Tprime(cfgs)
    assert np.max(np.abs(Tp_from_z2 - Tp_direct)) < 1e-12


# --------------------------------------------------------------------------------------
# Claim (ii): eq:sublinear -- ||M||_1 ~ K^p with p < 1 when l>0 configs are present.
#
# MEASURED: at the largest tractable scale used here (nested s+p+d+f up to K=100,
# build ~35 s aggregate) the fit gives p ~ 0.84.  The exponent is SCALE-DEPENDENT: it
# rises with basis size from ~0.78-0.80 at small K toward the full-K asymptote ~0.84
# (K~164, s8p8d8f8 + N=10 in the driver).  We assert a band bracketing the measured
# value, not the asymptote itself.
# --------------------------------------------------------------------------------------
@pytest.mark.slow
def test_onenorm_sublinear() -> None:
    """Entrywise 1-norm ||M||_1 grows sublinearly (p<1) with l>0 configs present."""
    scale_sets = [
        gen_configs(1, {0: 3, 1: 3}),            # K=9
        gen_configs(1, {0: 4, 1: 4}),            # K=16
        gen_configs(2, {0: 5, 1: 5, 2: 5}),      # K=31
        gen_configs(2, {0: 6, 1: 6, 2: 6}),      # K=46
        gen_configs(3, {0: 7, 1: 7, 2: 7, 3: 7}),  # K=74
        gen_configs(3, {0: 8, 1: 8, 2: 8, 3: 8}),  # K=100
    ]
    Ks, Ls = [], []
    for cts in scale_sets:
        _E, one_norm, K, _M = solve(cts)
        Ks.append(K)
        Ls.append(one_norm)

    assert len(Ks) >= 4
    p = float(np.polyfit(np.log(Ks), np.log(Ls), 1)[0])

    # Core claim: strictly sublinear.
    assert p < 1.0, f"||M||_1 exponent {p:.3f} is not sublinear"
    # Defensible band bracketing the measured value (~0.84 at K=100; asymptote ~0.84).
    assert 0.6 < p < 0.95, f"||M||_1 exponent {p:.3f} outside expected sublinear band"


# --------------------------------------------------------------------------------------
# Claim (iii): the L2 overlap metric S is ill-conditioned -- cond(S) grows from ~4
# into the thousands (paper "4 -> 3673") as the spdf basis grows.  The metric-free
# standard eigenproblem is the well-conditioned resolution.
# --------------------------------------------------------------------------------------
@pytest.mark.slow
def test_L2_metric_divergence() -> None:
    """cond(S) grows monotonically from ~O(1) into the thousands as the basis grows."""
    conds = []
    for N in (3, 4, 6, 8):
        lmax = min(3, N - 1)
        nmp = {l: N for l in range(lmax + 1)}
        _E_free, _E_S, cond_S, K = solve_with_metric(gen_configs(lmax, nmp))
        conds.append(cond_S)

    # Monotone growth across the whole sweep.
    for a, b in zip(conds, conds[1:]):
        assert b > a, f"cond(S) not monotone increasing: {conds}"

    # The smallest basis is well-conditioned (~4); the largest blows up to ~3673 (paper's "4→3673").
    assert conds[0] < 20.0, f"cond(S) at N=3 = {conds[0]:.1f} unexpectedly large"
    assert 2000.0 < conds[-1] < 6000.0, f"cond(S) at N=8 = {conds[-1]:.1f} not ~3673 (paper 4→3673)"
