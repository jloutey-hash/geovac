"""He / He-like fixed-k=Z hydrogenic CI on the PHYSICAL Hamiltonian (2026-09-25).

Backs the re-measured claims in Papers 13, FCI-atoms, 26 and 27 after the
retirement of the "graph-native CI 0.216%" headline.  In the k=Z hydrogenic
(n,l,m) basis T + V_ne is exactly diagonal (-Z^2/2n^2), so the physical
one-body operator has no off-diagonal elements; the legacy default of
``casimir_ci.build_graph_native_fci`` (h1_method='hybrid') adds a (1/16)*A
graph hop that is not part of the atom's Hamiltonian.

Driver with the full numbers (n_max = 2..7, both variants):
debug/he_graph_regression.py -> debug/data/he_graph_regression*.log.

Speed: V_ee is exactly linear in k, so H(Z) = Z^2 D1 + Z V1 is recomposed
from two cheap builds; the recomposition is itself checked below.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from geovac.casimir_ci import (  # noqa: E402
    build_graph_native_fci, build_fci_matrix, _build_graph_h1,
)
from geovac.lattice import GeometricLattice  # noqa: E402
from debug.archive.misc.entanglement_geometry import (  # noqa: E402
    build_1rdm_from_singlet_ci, compute_entanglement_measures,
)

E_HE = -2.903724377034
E_HM = -0.527751016544
E_HF_LIMIT_HE = -2.86168


def _pct(E, Eex=E_HE):
    return 100.0 * (E - Eex) / abs(Eex)


def _e0(H):
    return float(np.linalg.eigvalsh(H)[0])


_CACHE: dict = {}


def _blocks(n_max):
    """(D1, V1, configs, n_spatial): physical H(Z) = Z^2 D1 + Z V1."""
    if n_max not in _CACHE:
        H1 = build_graph_native_fci(Z=1, n_max=n_max, h1_method='exact')
        H2 = build_graph_native_fci(Z=2, n_max=n_max, h1_method='exact')
        D1 = 0.5 * (H2 - 2.0 * H1)
        V1 = H1 - D1
        lat = GeometricLattice(max_n=n_max)
        orbs = list(lat.states)
        cfg = [(i, j) for i in range(len(orbs)) for j in range(i, len(orbs))
               if orbs[i][2] + orbs[j][2] == 0]
        _CACHE[n_max] = (D1, V1, cfg, len(orbs))
    return _CACHE[n_max]


def _entropy(ci, cfg, ns):
    rho = build_1rdm_from_singlet_ci(ci, cfg, ns)
    return float(compute_entanglement_measures(rho)['von_neumann_entropy'])


# ---------------------------------------------------------------------------
# Paper 13 / FCI-atoms: the physical ladder
# ---------------------------------------------------------------------------

def test_physical_h1_is_diagonal_hydrogenic():
    h1, orbs = _build_graph_h1(2, 4, include_hop=False)
    assert np.count_nonzero(h1 - np.diag(np.diag(h1))) == 0
    for i, (n, _l, _m) in enumerate(orbs):
        assert h1[i, i] == -4.0 / (2.0 * n * n)


@pytest.mark.parametrize("n_max,E_ref,pct_ref", [
    (2, -2.833405, 2.422), (3, -2.842908, 2.094), (4, -2.846189, 1.981),
])
def test_he_physical_ladder(n_max, E_ref, pct_ref):
    """Paper 13 / FCI-atoms: fixed-k=Z CI with He's physical Hamiltonian."""
    E = _e0(build_graph_native_fci(Z=2, n_max=n_max, h1_method='exact'))
    assert E == pytest.approx(E_ref, abs=2e-6)
    assert _pct(E) == pytest.approx(pct_ref, abs=2e-3)
    assert E > E_HE                      # variational
    assert E > E_HF_LIMIT_HE             # the fixed-k bound-state basis does
                                         # not even reach the HF limit


def test_he_physical_equals_independent_k_route():
    """Same Hamiltonian via build_fci_matrix(k=Z) (get_h1_element route)."""
    for n in (2, 3):
        E_a = _e0(build_graph_native_fci(Z=2, n_max=n, h1_method='exact'))
        E_b = _e0(build_fci_matrix(Z=2, n_max=n, k_orb=2.0))
        assert abs(E_a - E_b) < 1e-10


def test_recomposition_is_exact():
    D1, V1, _, _ = _blocks(3)
    H3 = build_graph_native_fci(Z=3, n_max=3, h1_method='exact')
    assert np.max(np.abs(H3 - (9 * D1 + 3 * V1))) < 1e-12


# ---------------------------------------------------------------------------
# Paper 13 / 26: H^-, no "Z_c" boundary without the hop
# ---------------------------------------------------------------------------

def test_hminus_physical_bound_and_variational():
    D1, V1, _, _ = _blocks(4)
    E = _e0(D1 + V1)
    assert E < -0.5                      # bound relative to H + e-
    assert E > E_HM                      # variational (no over-binding)
    assert _pct(E, E_HM) == pytest.approx(2.22, abs=0.02)


def test_no_sub_exact_regime_without_hop():
    """E_phys(Z) stays above the 1/Z-series exact energy on Z in [1, 3]."""
    coef = [-1.0, 0.625, -0.157666405, 0.008698991, -0.000888587,
            -0.001036372, -0.000612917, -0.000372187, -0.000242004]
    D1, V1, _, _ = _blocks(3)
    for Z in np.linspace(1.0, 3.0, 21):
        e_ex = sum(c * Z ** (2 - k) for k, c in enumerate(coef))
        assert _e0(Z * Z * D1 + Z * V1) > e_ex + 5e-3


# ---------------------------------------------------------------------------
# Paper 26 Table I (physical columns) and S ~ Z^-2.52
# ---------------------------------------------------------------------------

def test_paper26_physical_table_and_scaling():
    D1, V1, cfg, ns = _blocks(4)
    ref = {2: (1.723, 0.0392), 3: (1.880, 0.0113), 10: (1.981, 0.0007)}
    Zs, Ss = [2, 3, 4, 5, 6, 8, 10], []
    for Z in Zs:
        H = Z * Z * D1 + Z * V1
        w, U = np.linalg.eigh(H)
        ci = U[:, 0]
        Ecorr = w[0] - H[0, 0]
        V = Z * V1
        f_vee = float(ci @ (V - np.diag(np.diag(V))) @ ci) / Ecorr
        S = _entropy(ci, cfg, ns)
        Ss.append(S)
        if Z in ref:
            assert f_vee == pytest.approx(ref[Z][0], abs=2e-3)
            assert S == pytest.approx(ref[Z][1], abs=6e-5)
    slope = np.polyfit(np.log(Zs), np.log(Ss), 1)[0]
    assert slope == pytest.approx(-2.522, abs=5e-3)


# ---------------------------------------------------------------------------
# Paper 27: Table ep1 physical rows, EP-2c physical fit, sec:cusp numbers
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_max,E_ref,S_ref", [(2, -2.8334, 0.03352), (3, -2.8429, 0.03772)])
def test_paper27_ep1_physical_rows(n_max, E_ref, S_ref):
    D1, V1, cfg, ns = _blocks(n_max)
    w, U = np.linalg.eigh(4 * D1 + 2 * V1)
    assert w[0] == pytest.approx(E_ref, abs=5e-5)
    assert _entropy(U[:, 0], cfg, ns) == pytest.approx(S_ref, abs=5e-6)
    # S_kin is EXACTLY zero: the one-body ground state is the 1s^2 determinant
    _, Uk = np.linalg.eigh(4 * D1)
    assert _entropy(Uk[:, 0], cfg, ns) == 0.0


def test_paper27_ep2c_physical_fit():
    D1, V1, cfg, ns = _blocks(3)
    Zs, S, wt = [2, 3, 4, 6, 8, 10], [], []
    for Z in Zs:
        Hk, V = Z * Z * D1, Z * V1
        _, U = np.linalg.eigh(Hk + V)
        S.append(_entropy(U[:, 0], cfg, ns))
        off = V - np.diag(np.diag(V))          # H1 eigenbasis = config basis
        wt.append(np.linalg.norm(off) / np.linalg.norm(Hk))
    g, b = np.polyfit(np.log(wt), np.log(S), 1)
    assert g == pytest.approx(2.565, abs=5e-3)
    assert np.exp(b) == pytest.approx(16.5, rel=0.02)
    assert S[0] == pytest.approx(3.772e-2, rel=1e-3)
    assert wt[0] == pytest.approx(0.0877, abs=2e-4)
    # w~ is exactly proportional to 1/Z for the physical H1
    assert np.allclose(np.array(wt) * np.array(Zs), wt[0] * 2, rtol=1e-12)


def test_paper27_cusp_numbers_physical():
    D1, V1, cfg, ns = _blocks(3)
    h1, V = 4 * D1, 2 * V1
    diag_frac = np.linalg.norm(np.diag(np.diag(V))) / np.linalg.norm(V)
    comm = np.linalg.norm(h1 @ V - V @ h1) / np.linalg.norm(h1)
    assert diag_frac == pytest.approx(0.945, abs=1e-3)
    assert comm == pytest.approx(0.106, abs=1e-3)


def test_paper27_gamma_surface_physical_nmax3():
    """tab:gamma_surface physical column at n_max=3: global (Z in
    {2,3,4,6,10,15,30}) = 2.321, local Z:15->30 = 1.979 (< 2)."""
    D1, V1, cfg, ns = _blocks(3)
    d = {}
    for Z in (2, 3, 4, 6, 10, 15, 30):
        Hk, V = Z * Z * D1, Z * V1
        _, U = np.linalg.eigh(Hk + V)
        off = V - np.diag(np.diag(V))
        d[Z] = (_entropy(U[:, 0], cfg, ns), np.linalg.norm(off) / np.linalg.norm(Hk))
    Zw = sorted(d)
    g = np.polyfit(np.log([d[z][1] for z in Zw]), np.log([d[z][0] for z in Zw]), 1)[0]
    loc = np.log(d[30][0] / d[15][0]) / np.log(d[30][1] / d[15][1])
    assert g == pytest.approx(2.321, abs=2e-3)
    assert loc == pytest.approx(1.979, abs=1e-3)
    assert loc < 2.0
