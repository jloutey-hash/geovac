"""Guard: the general-m triangle class of <F^2 V_ee> (debug/lih_pi_m16.tri_wcomb, used by the dense
get_tri_matrix_m and the sparse lih_gvee_sparse.tri_entries) respects the electron-relabeling
symmetry of |Psi|^2.

What wrong answer this rejects: the pre-2026-09-26 middle-mode weight, which gave the (A>=1, B=0)
combos pi/2 instead of pi (debug/sprint_gvee_defect_audit_memo.md).  Those combos exist only when the
middle (shared-vertex) density carries circulation, and the bug is L/R-asymmetric, so the 16
alpha<->beta-equivalent triangles of the orbit  f_12 f_13 / r_23  split into two values (alpha-middle
vs beta-middle).  It shifted the linexp g_Vee by up to 0.007 Ha (5.9 sigma vs VMC) while every
dense==sparse guard stayed green (both used the same wrong weight).  Fire test built in: with
M._WCOMB_LEGACY = True the orbit spread must be large.

Needs the untracked caches debug/data/lih_pi_m16_prim.npz + lih_pi_m16_fci_pi_{canon,no}.npz; SKIPs
otherwise.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEBUG = os.path.join(ROOT, "debug")
DATA = os.path.join(DEBUG, "data")
_NEED = ["lih_pi_m16_prim.npz", "lih_pi_m16_fci_pi_canon.npz", "lih_pi_m16_fci_pi_no.npz"]

pytestmark = [pytest.mark.slow,
              pytest.mark.skipif(not all(os.path.exists(os.path.join(DATA, f)) for f in _NEED),
                                 reason="pi M=16 caches not present (run debug/lih_pi_m16.py g0)")]

TRUNC = 5e-2
# the orbit of f_12 f_13 / r_23 under (12), (34), (13)(24): alpha-middle and beta-middle members
ORBIT = [((1, 2), (1, 3), (2, 3)), ((1, 2), (1, 4), (2, 4)), ((1, 2), (2, 3), (1, 3)),
         ((1, 2), (2, 4), (1, 4)), ((3, 4), (1, 3), (1, 4)), ((3, 4), (1, 4), (1, 3)),
         ((3, 4), (2, 3), (2, 4)), ((3, 4), (2, 4), (2, 3))]


@pytest.fixture(scope="module")
def setup():
    import geovac.lih_r12ci.kernels as KG
    import geovac.lih_r12ci.pi_m16 as M
    import geovac.lih_r12ci.gvee_sparse as S
    old = KG.USE_EXACT_NEUMANN
    KG.USE_EXACT_NEUMANN = True
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    refm, _ = M.build_ref_from_prims(ref, "guard pi tiny tri-sym", cfg_tag='pi', trunc=TRUNC, verbose=False)
    gem = M.GeminalM('linexp')
    yield dict(M=M, S=S, gem=gem, refm=refm)
    KG.USE_EXACT_NEUMANN = old


def _orbit_values(st):
    S, gem, refm = st['S'], st['gem'], st['refm']
    vals = []
    for p, q, r in ORBIT:
        ctx = S.Ctx(gem, refm, 1.0)
        col = S.TriCollector()
        rem_f, kk, gk = S._triple_args(p, q, r, gem)
        assert S.vee_reduce_sparse(rem_f, r, kk, gk, ctx, col) == 0.0 and col.recs, "not a triangle"
        vals.append(S.eval_triangles(ctx, col))
    return np.array(vals)


def test_weight_is_left_right_symmetric(setup):
    M = setup['M']
    assert not M._WCOMB_LEGACY
    for A in range(1, 5):
        for s in (+1, -1):
            assert M.tri_wcomb(A, s, 0, 0) == M.tri_wcomb(0, 0, A, s) == pytest.approx(np.pi)
    assert M.tri_wcomb(0, 0, 0, 0) == pytest.approx(2 * np.pi)
    assert M.tri_wcomb(2, 1, 1, -1) == pytest.approx(np.pi / 2)


def test_reference_has_circulating_middle_densities(setup):
    assert np.any(setup['refm'].mA != 0)


def test_alpha_beta_orbit_is_degenerate(setup):
    v = _orbit_values(setup)
    assert abs(v).max() > 1e-3
    assert v.max() - v.min() < 1e-6 * abs(v).max(), v


def test_firetest_legacy_weight_breaks_the_orbit(setup):
    M = setup['M']
    M._WCOMB_LEGACY = True
    try:
        v = _orbit_values(setup)
    finally:
        M._WCOMB_LEGACY = False
    assert v.max() - v.min() > 1e-5 * abs(v).max(), ("legacy half-weight did not split the orbit -> "
                                                    "guard asleep", v)
