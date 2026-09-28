"""Guard: the O(nnz) g_Vee reducer (debug/lih_gvee_sparse.py) equals the dense reducer
(debug/lih_pi_m16.hg_Vee_m, the oracle) on a small pi-containing reference.

What wrong answer this rejects: any sparse restriction that drops or mis-weights index tuples
(a skipped topology class, a wrong circulation / L_z mask, a dropped triangle mode combo, a
mis-indexed compound field).  Both reducers are run from IDENTICAL random streams (the triangle's
randomized eigh/SVD), so equality is to roundoff, not to the ~1e-7 randomized-SVD noise.
Built-in fire test: with TREE4-type (nested-chain) kept-pair terms dropped, equality must FAIL.

Needs the untracked caches debug/data/lih_pi_m16_prim.npz + lih_pi_m16_fci_pi_{canon,no}.npz
(written by debug/lih_pi_m16.py g0); SKIPs otherwise.  Fire tests: debug/firetest_lih_gvee_sparse.py.
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

TRUNC = 5e-2          # tiny pi reference: 5 dets, nact 36, 18 m!=0 pairs, all topology classes live
TOL = 1e-10


@pytest.fixture(scope="module")
def setup():
    sys.path.insert(0, DEBUG)
    sys.path.insert(0, ROOT)
    import geovac.lih_r12ci.kernels as KG
    import geovac.lih_r12ci.triangle as TRI
    import lih_pi_m16 as M
    import lih_gvee_sparse as S
    saved_flag = KG.USE_EXACT_NEUMANN
    KG.USE_EXACT_NEUMANN = True
    ref = M.build_primitives(M.CFG, M.PRIM_CACHE_PI, verbose=False)
    refm, _info = M.build_ref_from_prims(ref, "guard pi tiny", cfg_tag='pi', trunc=TRUNC, verbose=False)
    gem = M.GeminalM('linexp')
    Fbar = M.fbar(gem, refm)

    def reset_rng():
        TRI.rng = np.random.default_rng(20260921)
        M._RSVD_RNG = np.random.default_rng(20260923)
        TRI._MODE_CACHE.clear()
        M._TRIM_CACHE.clear()

    reset_rng()
    dense = M.hg_Vee_m(gem, refm, Fbar, force_general=True)
    # restore the module-global flag on teardown: leaking it True made the
    # legacy-anchored tests/test_lih_r12ci*.py fail when run in the same session
    yield dict(M=M, S=S, gem=gem, refm=refm, Fbar=Fbar, dense=dense, reset_rng=reset_rng)
    KG.USE_EXACT_NEUMANN = saved_flag


def _sparse(st, mem_cap_gb=2.0):
    st['reset_rng']()
    return st['S'].hg_Vee_sparse(st['gem'], st['refm'], st['Fbar'], mem_cap_gb=mem_cap_gb)


def test_reference_exercises_pi(setup):
    """the guard reference must carry m!=0 pairs (else the m-machinery is untested)."""
    assert np.any(setup['refm'].mA != 0)


def test_sparse_equals_dense(setup):
    names = ['h_Vee', 'g_Vee', 'Vee', 'FVee', 'F2Vee']
    sp = _sparse(setup)
    for nm, d, s in zip(names, setup['dense'], sp):
        assert abs(d - s) < TOL * max(1.0, abs(d)), (nm, d, s, d - s)


def test_sparse_equals_dense_tiny_chunks(setup):
    """streaming path: a memory cap that forces many chunks gives the same answer."""
    sp = _sparse(setup, mem_cap_gb=0.02)
    assert abs(sp[1] - setup['dense'][1]) < TOL * max(1.0, abs(setup['dense'][1]))


def test_firetest_dropping_tree4_breaks_equality(setup):
    """built-in fire test: skipping the nested-chain (TREE4) kept-pair terms must break equality."""
    S = setup['S']
    S._DROP_NESTED = True
    try:
        sp = _sparse(setup)
    finally:
        S._DROP_NESTED = False
    assert abs(sp[4] - setup['dense'][4]) > 1e-4, "dropping TREE4 did not change F2Vee -> guard asleep"
