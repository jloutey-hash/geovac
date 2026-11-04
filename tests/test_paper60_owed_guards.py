"""Paper 60 -- owed guards from the group2 /qa FULL run (2026-11-03/04, docs/qa/group2.debt.md).

(1) eq:sublinear and its s-only twin, frozen on the CONVERGED box.  The paper prints
    ||M||_1 ~ K^0.82 (full s+p+d+f, helium, K = 74..164) and K^0.73 (s-only, the same
    n_max = 7..10 window); the registry provenance gives the 4-dp values 0.8193 and 0.7261,
    measured on R_MAX >= 3 n_max^2.  The existing tests only assert sublinear BANDS
    (0.6 < p < 0.95; p < 1), which would accept the retired 0.842 / 0.78 -- the figures a
    60-bohr box produced.  Measured before writing this guard (2026-11-04, box
    max(80, 3 n^2), uniform 18000 points; the graded 24000-point mesh agrees to 1e-5):

        s-only  n = 7..10   K = 28, 36, 45, 55      ||M||_1 = 33.73 .. 55.07    p = 0.7261
        spdf    n = 7..10   K = 74, 100, 130, 164   ||M||_1 = 85.06 .. 163.25   p = 0.8193

    Pinned at +-0.005, which rejects both retired values (0.842 is 0.023 away; 0.78 is 0.054).

(2) sec:resource conditioning numbers (tab:resource + the ratio-dependence caveat), from the
    tracked `_gaussian_metric` / `_sw_metric_mom` helpers of tests/test_paper60_sturmian.py,
    H2+ at R = 2 bohr, N = 2 n_max = 4..16.  Measured 2026-11-04:

        ratio 1.6   cond(N=16) = 1.467e6   fitted exponent 6.26   (prose: ~1.5e6, N^6)
        ratio 2     cond(N=16) = 1.991e4                          (tab:resource: 1.99e4)
        ratio 3     fitted exponent 1.571; cond(N=16)/SW = 4.23    (prose: N^1.6, ~4x)
        SW full     cond(N=16) = 91.84                            (tab:resource: 91.8)

    NOTE: the prose calls the ratio-3 factor "the SW prefactor"; what is ~4.2x is the ratio of
    condition numbers at N = 16.  The ratio of fitted power-law prefactors is ~9.2.  This guard
    pins the N = 16 ratio, which is the number the prose's 4x matches.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import geovac.sturmian_secular as SS                                       # noqa: E402
import geovac.sturmian_variational as SV                                   # noqa: E402
from test_paper60_sturmian import _gaussian_metric, _sw_metric_mom         # noqa: E402


@pytest.fixture(scope="module", autouse=True)
def _restore_secular_grid():
    """`SV.set_grid` mutates module globals of `geovac.sturmian_secular`; restore them."""
    saved = dict(r=SS.r, dr=SS.dr, r2=SS.r2, R_MAX=SS.R_MAX, N_GRID=SS.N_GRID,
                 fwd=SS._ctrap_fwd, rev=SS._ctrap_rev)
    yield
    SS.r, SS.dr, SS.r2 = saved['r'], saved['dr'], saved['r2']
    SS.R_MAX, SS.N_GRID = saved['R_MAX'], saved['N_GRID']
    SS._ctrap_fwd, SS._ctrap_rev = saved['fwd'], saved['rev']
    SS.reset_caches()


def _converged_onenorm_exponent(lmax: int):
    K, L = [], []
    for n in (7, 8, 9, 10):
        SV.set_grid(max(80, 3 * n * n), 18000)        # converged box, R_MAX >= 3 n_max^2
        cfgs = SS.build_configs(SS.gen_configs(lmax, {l: n for l in range(lmax + 1)}))
        M = SS.build_M(cfgs, Z=2.0)
        K.append(len(cfgs))
        L.append(float(np.abs(M).sum()))
    return K, float(np.polyfit(np.log(K), np.log(L), 1)[0])


@pytest.mark.slow
def test_paper60_sonly_onenorm_exponent_converged_box():
    K, p = _converged_onenorm_exponent(0)
    assert K == [28, 36, 45, 55], K
    assert abs(p - 0.7261) < 0.005, f"s-only converged-box exponent moved: {p:.4f}"


@pytest.mark.slow
def test_paper60_spdf_onenorm_exponent_converged_box():
    K, p = _converged_onenorm_exponent(3)
    assert K == [74, 100, 130, 164], K
    assert abs(p - 0.8193) < 0.005, f"spdf converged-box exponent moved: {p:.4f}"


# ---------------------------------------------------------------- (2) conditioning numbers
R_BOND = 2.0
NMAX = (2, 4, 6, 8)
NS = np.array([2 * n for n in NMAX], float)


def _conds(ratio=None):
    if ratio is None:
        return np.array([np.linalg.cond(_sw_metric_mom(n, R_BOND)) for n in NMAX])
    return np.array([np.linalg.cond(_gaussian_metric(n, R_BOND, ratio=ratio)) for n in NMAX])


def _exp(c):
    return float(np.polyfit(np.log(NS), np.log(c), 1)[0])


def test_paper60_gaussian_ratio16_blows_up():
    c = _conds(1.6)
    assert abs(c[-1] / 1.467e6 - 1) < 0.02, f"ratio-1.6 cond at N=16: {c[-1]:.4g}"
    assert 5.9 < _exp(c) < 6.6, f"ratio-1.6 exponent (prose N^6): {_exp(c):.3f}"


def test_paper60_gaussian_ratio2_table_value():
    c = _conds(2.0)
    assert abs(c[-1] / 1.99e4 - 1) < 0.01, f"ratio-2 cond at N=16 (tab:resource 1.99e4): {c[-1]:.4g}"


def test_paper60_sw_full_table_value():
    c = _conds()
    assert abs(c[-1] - 91.8) < 0.5, f"SW full cond at N=16 (tab:resource 91.8): {c[-1]:.4g}"


def test_paper60_gaussian_ratio3_exponent_and_factor():
    g3, sw = _conds(3.0), _conds()
    assert abs(_exp(g3) - 1.57) < 0.03, f"ratio-3 exponent (prose N^1.6): {_exp(g3):.3f}"
    assert abs(g3[-1] / sw[-1] - 4.2) < 0.15, f"ratio-3 / SW at N=16 (prose ~4x): {g3[-1] / sw[-1]:.3f}"
