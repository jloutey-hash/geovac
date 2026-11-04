"""Paper 18 thm:ac_factorization -- both directions, on production data.

Guard (group3 guard pass, 2026-11-03); previously backed only by a debug toy.

Claim (Paper 18, scoped 2026-11-03): for D = D_GV (x) 1_F + gamma_GV (x) D_F
with gamma_GV an outer grading ANTICOMMUTING with D_GV, the cross term of D^2
vanishes and Tr e^{-tD^2} = Tr e^{-tD_GV^2} Tr e^{-tD_F^2}.  The production
combined Dirac (geovac.almost_commutative, Paper 32 Sprint H1) uses
gamma_GV = sign(D_GV), which COMMUTES with D_GV, so the factorisation fails at
non-zero Yukawa; it holds for the graded product D_GV (x) gamma_F + 1 (x) D_F.

Data: production CH Dirac D_GV at n_max = 2 and the production electroweak
D_F at Yukawa (e, nu) = (0.3, 0.2), so D_F != 0 (non-vacuous).
Holds-direction runs on the chirality-doubled model D' = sigma_z (x) D_GV with
gamma' = sigma_x (x) 1 (anticommuting).
Wrong answers excluded: "the sign(D_GV) twist also factorises" (fails-test);
"a commuting grading also factorises" (holds-tests fire if gamma' is swapped
for a commuting one).
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.linalg import expm

from geovac.almost_commutative import minimal_electroweak_triple

TS = (0.1, 0.5, 1.0)


@pytest.fixture(scope="module")
def tri():
    t = minimal_electroweak_triple(2, yukawa_e=0.3, yukawa_nu=0.2)
    assert np.linalg.norm(t.D_F()) > 0.1
    return t


def _tr_heat(D: np.ndarray, t: float) -> float:
    return float(np.trace(expm(-t * D @ D)).real)


def _rel_dev(D, DA, DB) -> float:
    return max(abs(_tr_heat(D, t) / (_tr_heat(DA, t) * _tr_heat(DB, t)) - 1) for t in TS)


def test_holds_for_anticommuting_outer_grading(tri):
    DG, DF = tri.D_GV(), tri.D_F()
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sz = np.diag([1.0, -1.0]).astype(complex)
    Dp = np.kron(sz, DG)                       # doubled outer Dirac
    gp = np.kron(sx, np.eye(DG.shape[0]))      # anticommuting outer grading
    assert np.linalg.norm(gp @ Dp + Dp @ gp) < 1e-12
    nF = DF.shape[0]
    D = np.kron(Dp, np.eye(nF)) + np.kron(gp, DF)
    # D^2 = D'^2 (x) 1 + 1 (x) D_F^2 (no cross term)
    D2_pred = np.kron(Dp @ Dp, np.eye(nF)) + np.kron(np.eye(Dp.shape[0]), DF @ DF)
    assert np.linalg.norm(D @ D - D2_pred) < 1e-10
    assert _rel_dev(D, Dp, DF) < 1e-12


def test_holds_for_graded_product(tri):
    DG, DF, gF = tri.D_GV(), tri.D_F(), tri.gamma_F()
    assert np.linalg.norm(gF @ DF + DF @ gF) < 1e-12
    D = np.kron(DG, gF) + np.kron(np.eye(DG.shape[0]), DF)
    assert _rel_dev(D, DG, DF) < 1e-12


def test_fails_for_production_sign_twist(tri):
    DG = tri.D_GV()
    g = tri.gamma_GV()
    assert np.linalg.norm(g @ DG - DG @ g) < 1e-12          # commutes: hypothesis fails
    assert np.linalg.norm(g @ DG + DG @ g) > 1.0
    D = tri.dirac_combined()
    assert _rel_dev(D, DG, tri.D_F()) > 0.05
