"""Paper 32, G3 / G4a: which gradings anticommute with which combined SM Dirac.

geovac/standard_model_triple.py (and almost_commutative.py) build
    D_mod = D_GV (x) 1_F + chi (x) D_F,
with chi = diag(chirality label) of the full-Dirac basis.  chi COMMUTES with
D_GV, so D_mod is not the grading-twisted product.  The production grading is
gamma_GV = sigma_x (x) 1 (chirality_grading.build_gamma_GV), which
ANTICOMMUTES with D_GV.  Measured facts pinned here, n_max = 1, 2, 3:

  * D_F = 0 at zero Yukawa, so every statement "at Y = 0" is about D_GV alone;
  * with D_mod, the product grading gamma_GV (x) gamma_F fails linearly in Y
    (16 sqrt2, 32 sqrt2, 32 sqrt5 at uniform Y = 1);
  * with the untwisted sum D_GV (x) 1 + 1 (x) D_F, or the grading-twisted
    D_GV (x) 1 + gamma_GV (x) D_F, the product grading anticommutes with D at
    every Yukawa AND the reality condition J D = +D J holds;
  * s (x) gamma_F with s = chi: ||{s (x) gamma_F, D}|| = 2 ||D_GV|| sqrt(dim H_F)
    at Y = 0 -- an identity ({chi, D_GV} = 2 chi D_GV with chi D_GV = |D_GV|),
    not a measured obstruction.

Wrong answers rejected: "the product grading grades the combined Dirac at zero
Yukawa only" read as structural (it holds at every Yukawa for both standard
products); "the module's D is the twisted product"; and a reality condition
lost by the grading-compatible products.
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.chirality_grading import build_gamma_GV
from geovac.standard_model_triple import standard_model_triple

N_MAX = (1, 2, 3)
# ||{gamma_GV (x) gamma_F, D_mod}||_F at uniform Yukawa Y = 1.
MOD_AT_Y1 = {1: 16 * np.sqrt(2), 2: 32 * np.sqrt(2), 3: 32 * np.sqrt(5)}


def _triple(n_max: int, y: float):
    return standard_model_triple(n_max, yukawa_nu=y, yukawa_e=y, yukawa_u=y, yukawa_d=y)


def _anti(g, D) -> float:
    return float(np.linalg.norm(g @ D + D @ g))


def _diracs(T, sx):
    I_F, I_G = np.eye(T.dim_F), np.eye(T.dim_GV)
    base = np.kron(T.D_GV(), I_F)
    return {"module": T.dirac_combined(),
            "untwisted": base + np.kron(I_G, T.D_F()),
            "grading-twisted": base + np.kron(sx, T.D_F())}


@pytest.mark.parametrize("n_max", N_MAX)
def test_module_twist_commutes_with_dgv(n_max):
    T = _triple(n_max, 1.0)
    chi, D = T.gamma_GV(), T.D_GV()
    sx = build_gamma_GV(n_max).matrix
    assert np.linalg.norm(chi @ D - D @ chi) < 1e-10
    assert np.linalg.norm(sx @ D + D @ sx) < 1e-10
    assert np.linalg.norm(_triple(n_max, 0.0).D_F()) == 0.0


@pytest.mark.parametrize("n_max", N_MAX)
def test_product_grading_by_construction(n_max):
    sx = build_gamma_GV(n_max).matrix
    for y in (0.0, 0.5, 1.0):
        T = _triple(n_max, y)
        g = np.kron(sx, T.gamma_F())
        U = np.asarray(T.real_structure_combined())
        Uinv = np.linalg.inv(U)
        for name, D in _diracs(T, sx).items():
            assert np.linalg.norm(U @ D.conj() @ Uinv - D) < 1e-9, (n_max, y, name)
            expected = MOD_AT_Y1[n_max] * y if name == "module" else 0.0
            assert abs(_anti(g, D) - expected) < 1e-8, (n_max, y, name)


@pytest.mark.parametrize("n_max", N_MAX)
def test_sign_grading_norm_is_an_identity(n_max):
    T = _triple(n_max, 0.0)
    val = _anti(T.chirality_combined(), T.dirac_combined())
    assert abs(val - 2 * np.linalg.norm(T.D_GV()) * np.sqrt(T.dim_F)) < 1e-9
