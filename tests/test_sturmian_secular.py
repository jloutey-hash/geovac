"""Regression tests for the isoenergetic generalized-Sturmian He secular builder.

Backs four load-bearing Paper 60 (group2) claims, ported from the
``debug/sturmian_he_lmax.py`` driver into ``geovac/sturmian_secular.py``:

  (i)   single-config 1s^2  ->  E = -2.847 Ha  (textbook variational He);
  (ii)  the entrywise 1-norm ``||M||_1`` grows SUBLINEARLY with the config count K,
        approaching ``~K^0.82`` over the window K = 74..164 on a CONVERGED radial
        box (Paper 60 ``eq:sublinear``).  ``K^0.84`` is RETIRED -- a 60-bohr-box
        artifact.  This file measures 0.842 on the module default box;  that is
        box-limited, which is why the assertion below is a BAND, not a value;
  (iii) restoring the L2 overlap metric S (generalized eigenproblem ``M B = p S B``)
        is NOT ill-conditioned -- the "~4 into the thousands" growth is a
        radial-box artifact (converged cond(S) = 16.0, growing ~0.12*K).  The
        paper's "4 -> 3673" is RETIRED [retracted 2026-09-07: p60-l2-metric-diverges];  the real reason to keep the
        metric-free form is eq:scale_lock, not conditioning;
  (iv)  ``[eq:secular]`` the interelectron matrix T' is a matrix of PURE NUMBERS,
        independent of the nuclear charge Z.

The K-ladder (ii) and metric-conditioning (iii) tests are marked ``slow`` (build_M at
K~100 is ~15 s; the ladder/divergence sweeps aggregate several such builds).
"""
import numpy as np
import pytest

from geovac.sturmian_secular import (  # noqa: F401
    build_S, build_configs,
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
# is box-dependent.  There is NO asymptote: on a converged box the local slope FALLS
# monotonically (0.827 -> 0.766 across K=100..514), so no window fit is stable.  We
# assert a BAND bracketing the measured value, never an asymptote.
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
    # Defensible band bracketing the measured value (0.842 on this box; converged
    # window fit 0.82).  Deliberately a band: there is no asymptote to pin.
    assert 0.6 < p < 0.95, f"||M||_1 exponent {p:.3f} outside expected sublinear band"


# --------------------------------------------------------------------------------------
# Claim (iii), AS CORRECTED 2026-09-07: the L2 overlap metric S is NOT
# ill-conditioned.  The "4 -> 3673" growth is a radial-box artifact [retracted 2026-09-07: p60-l2-metric-diverges]; converged,
# cond(S) = 16.0 and grows ~0.12*K -- an ordinary Gram matrix.  This test pins the
# artifact AS an artifact.  (The metric-free form is preferable for a different
# reason entirely: it exists only at the locked scale, Paper 60 eq:scale_lock.)
# --------------------------------------------------------------------------------------
@pytest.mark.slow
def test_L2_metric_conditioning_is_box_dependent() -> None:
    """cond(S) is a RADIAL-DOMAIN artifact above n_max^2 = R_MAX, not a property
    of the L2 metric.

    Replaces test_L2_metric_divergence (2026-09-07), which asserted
    `2000 < cond(S)[N=8] < 6000` -- i.e. it PINNED the artifact and would have
    failed if anyone repaired the engine. Paper 60's "cond(S) climbs 4 -> 3673" [retracted 2026-09-07: p60-l2-metric-diverges]
    was withdrawn on the strength of this measurement.

    What is actually true, verified at R_MAX = 60/120/240/480:

        N=3  K=10   n_max^2=9    4.07  4.07  4.07  4.07   <- box-independent
        N=4  K=20   n_max^2=16   5.69  5.69  5.69  5.69   <- box-independent
        N=6  K=52   n_max^2=36   9.79  9.76  9.76  9.76   <- box-independent
        N=8  K=100  n_max^2=64   3673  16.01 16.01 16.01  <- ARTIFACT

    The switch-on is exactly where n_max^2 first exceeds R_MAX: orbitals that
    outgrow the domain are truncated and then renormalised, which manufactures
    near-degenerate rows in the Gram matrix. Converged, cond(S) grows mildly
    (~0.12*K).

    The wrong answers this rejects:
      (a) cond(S) blowing up at a basis whose orbitals FIT the domain -- that
          would mean the L2 metric really is ill-conditioned;
      (b) cond(S) being box-INdependent where the orbitals overflow -- that
          would mean the 3673 was real after all;
    NOT claimed: this does not fire when the module's default R_MAX changes.
    It sets the domain itself, so a repair of the engine default leaves it green
    -- deliberately. The OLD assertion failed on repair, which is exactly the
    defect being removed here; a guard that punishes the fix is worse than none.
    (Fire-tested 2026-09-07: R_MAX 60 -> 400 correctly does NOT fire; removing
    the grid renormalisation in hyd_radial DOES.)
    """
    import numpy as np

    import geovac.sturmian_secular as _S

    def _cond_at(N: int, rmax: float) -> float:
        r_o, dr_o, r2_o, rm_o = _S.r, _S.dr, _S.r2, _S.R_MAX
        try:
            npts = int(_S.N_GRID * rmax / 60.0)
            _S.R_MAX = rmax
            _S.r = np.linspace(1e-7, rmax, npts)
            _S.dr = _S.r[1] - _S.r[0]
            _S.r2 = _S.r * _S.r
            lmax = min(3, N - 1)
            cfgs = build_configs(gen_configs(lmax, {l: N for l in range(lmax + 1)}))
            return float(np.linalg.cond(build_S(cfgs)))
        finally:
            _S.r, _S.dr, _S.r2, _S.R_MAX = r_o, dr_o, r2_o, rm_o

    # (1) Where the orbitals FIT the domain, cond(S) does not depend on it.
    for N in (3, 4, 6):
        c60, c240 = _cond_at(N, 60.0), _cond_at(N, 240.0)
        assert N * N < 60, f"N={N} was supposed to fit inside the 60-bohr box"
        assert abs(c60 - c240) / c240 < 5e-3, (
            f"N={N}: cond(S) moved with the box ({c60:.3f} vs {c240:.3f}) "
            f"even though n_max^2={N*N} fits inside it")
        assert c60 < 20.0, f"N={N}: cond(S)={c60:.1f} unexpectedly large"

    # (2) Where they do NOT fit, the small box inflates it by orders of magnitude.
    c60_8, c240_8 = _cond_at(8, 60.0), _cond_at(8, 240.0)
    assert 8 * 8 > 60, "N=8 was supposed to overflow the 60-bohr box"
    assert c60_8 > 100.0 * c240_8, (
        f"N=8: the 60-bohr box no longer inflates cond(S) "
        f"({c60_8:.1f} vs {c240_8:.2f}) -- has the engine been repaired?")

    # (3) Converged, it is an ordinary Gram matrix, not an ill-conditioned one.
    assert c240_8 < 40.0, (
        f"converged cond(S) at N=8 is {c240_8:.1f}; Paper 60's withdrawn "
        f"'4 -> 3673' claim would need this to be in the thousands")  # [retracted 2026-09-07: p60-l2-metric-diverges]
