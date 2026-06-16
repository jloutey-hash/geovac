"""
Verification test for Paper 24 (Bargmann-Segal lattice), Section
sec:entanglement-rigidity:  the two-fermion harmonic-oscillator
entanglement-rigidity corollary.

Paper 24 claims (corollary to the HO rigidity theorem, Thm. rigidity):

    The spatial 1-RDM von Neumann entropy of the closed-shell
    two-fermion ground state on the Bargmann-Segal graph is
    identically zero for any central two-body interaction V(r_12),
    and the kinetic HO Hamiltonian commutes with the interaction:

        S_full = 0,
        || [H_HO, V] ||_F / || H_HO ||_F < 1e-15.

    Verified at N_max in {2, 3} with the Minnesota NN singlet
    potential at hbar*omega = 10 MeV.

Mechanism (Moshinsky-Talmi lab-to-relative transformation): any
central V(r_12) on the 3D HO basis preserves the total HO quantum
number N_tot = N_rel + N_CM, so V is block-diagonal in N_tot, hence
commutes with the (also N_tot-diagonal) one-body H_HO; the closed-shell
ground state lies entirely in the lowest N_tot=0 sector, which is
one-dimensional --- a single Slater determinant |(0s)^2>, occupation
(2, 0, ..., 0), zero entanglement entropy.

This is the structural dual of the *nonzero* universal Coulomb-S^3
two-variable scaling of Paper 27.  Paper 24 cites the empirical
verification to "Paper 27, Sec. VII.A", but the in-suite Paper 27
test file (tests/test_paper27_entropy.py) imports its 1-RDM / entropy
helpers from an archived debug module; this file is the self-contained
backing test pinned to the Paper 24 corollary (closes the /qa group3
C1 coverage gap).  The ground-state machinery is the production module
geovac.nuclear.ho_two_fermion; the 1-RDM and von Neumann entropy are
computed by the small self-contained helpers below (no debug imports).

Reference values (frozen, hbar*omega = 10 MeV):
    N_max = 2:  n_spatial=10, n_configs=15, E0=22.185 MeV, S=0, occ=(2,0,0,0)
    N_max = 3:  n_spatial=20, n_configs=42, E0=22.185 MeV, S=0, occ=(2,0,0,0)
    rel commutator norm ~ 2.7e-16 (N_max=2), ~ 2.4e-16 (N_max=3)
"""

from __future__ import annotations

import os
import sys
from typing import List, Tuple

import numpy as np
import pytest

# Ensure project root on path (mirrors conftest).
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from geovac.nuclear.ho_two_fermion import (  # noqa: E402
    build_decomposed_ho_hamiltonians,
)


# ---------------------------------------------------------------------------
# Self-contained 1-RDM / entropy helpers (no archived-debug dependency).
# These reproduce debug/archive/misc/entanglement_geometry.py for the
# spatial-singlet two-fermion CI, so the Paper 24 corollary is backed
# entirely by in-suite production code.
# ---------------------------------------------------------------------------

def _build_singlet_1rdm(
    ci_coeffs: np.ndarray,
    configs: List[Tuple[int, int]],
    n_spatial: int,
) -> np.ndarray:
    """Spatial 1-RDM (Tr = 2) from spin-adapted singlet CI coefficients.

    Each config I = (i, j) is the spatial singlet
    [phi_i(1) phi_j(2) + phi_j(1) phi_i(2)] / N_IJ, with N_IJ = sqrt(2)
    for i != j and 1 for i == j.  Tracing out electron 2 gives

        rho_{a,c} = 2 * sum_{I,J} c_I c_J / (N_I N_J)
                      * sum_{(a,b) in I, (c,d) in J} delta(b, d).
    """
    rho = np.zeros((n_spatial, n_spatial))
    for I, (i, j) in enumerate(configs):
        for J, (p, q) in enumerate(configs):
            coeff = ci_coeffs[I] * ci_coeffs[J]
            if abs(coeff) < 1e-16:
                continue
            N_I = np.sqrt(2.0) if i != j else 1.0
            N_J = np.sqrt(2.0) if p != q else 1.0
            I_perms = [(i, j)] + ([(j, i)] if i != j else [])
            J_perms = [(p, q)] + ([(q, p)] if p != q else [])
            for a, b in I_perms:
                for c, d in J_perms:
                    if b == d:
                        rho[a, c] += coeff / (N_I * N_J)
    rho *= 2.0  # spin trace: closed-shell singlet, Tr = N_electrons = 2
    return rho


def _von_neumann_entropy(rho: np.ndarray) -> Tuple[float, np.ndarray]:
    """Return (S, occupations) for a 2-electron spatial 1-RDM (Tr = 2).

    Natural-orbital occupations n_i sum to 2; the normalized 1-RDM has
    eigenvalues lambda_i = n_i / 2 and S = -sum lambda_i log lambda_i.
    """
    occ = np.linalg.eigvalsh(rho)
    occ = np.sort(occ)[::-1]
    occ = np.maximum(occ, 0.0)
    lam = occ / 2.0
    S = 0.0
    for x in lam:
        if x > 1e-15:
            S -= x * np.log(x)
    return float(S), occ


# Frozen reference values (hbar*omega = 10 MeV, Minnesota NN singlet).
_REFERENCE = {
    2: {'n_spatial': 10, 'n_configs': 15, 'E0_MeV': 22.185},
    3: {'n_spatial': 20, 'n_configs': 42, 'E0_MeV': 22.185},
}


def _build(N_max: int):
    return build_decomposed_ho_hamiltonians(N_max=N_max, hw=10.0)


@pytest.mark.parametrize('N_max', [2, 3])
def test_paper24_ho_zero_entanglement_entropy(N_max):
    """S_full = 0 and the GS is a single Slater determinant |(0s)^2>.

    Paper 24 corollary, entanglement-rigidity leg: the 2-fermion HO
    ground state has zero von Neumann entropy with occupation
    (2, 0, ..., 0).
    """
    data = _build(N_max)
    ref = _REFERENCE[N_max]

    # Build shape sanity.
    assert data['n_spatial'] == ref['n_spatial']
    assert data['n_configs'] == ref['n_configs']

    H_full = data['H_full']
    configs = data['configs']
    n_spatial = data['n_spatial']

    eigs, vecs = np.linalg.eigh(H_full)
    ci = vecs[:, 0]

    rho = _build_singlet_1rdm(ci, configs, n_spatial)
    # 1-RDM is a valid 2-electron density: trace = 2.
    assert np.trace(rho) == pytest.approx(2.0, abs=1e-10), \
        f'Tr(rho) = {np.trace(rho):.6f} at N_max={N_max}'

    S, occ = _von_neumann_entropy(rho)

    # Single Slater determinant: occupation (2, 0, 0, ...).
    assert occ[0] == pytest.approx(2.0, abs=1e-10), \
        f'top occupation = {occ[0]:.6f} at N_max={N_max}'
    for k in range(1, min(4, len(occ))):
        assert abs(occ[k]) < 1e-10, \
            f'occupation[{k}] = {occ[k]:.3e} nonzero at N_max={N_max}'

    # Von Neumann entropy identically zero.
    assert S < 1e-10, f'S_full = {S:.3e} (expected 0) at N_max={N_max}'

    # Ground-state energy is the N_tot=0-sector |(0s)^2> determinant,
    # basis-independent (same at N_max=2 and N_max=3).
    assert eigs[0] == pytest.approx(ref['E0_MeV'], abs=5e-2), \
        f'GS energy = {eigs[0]:.4f} MeV (expected {ref["E0_MeV"]}) ' \
        f'at N_max={N_max}'


@pytest.mark.parametrize('N_max', [2, 3])
def test_paper24_ho_kinetic_interaction_commute(N_max):
    """|| [H_HO, V] ||_F / || H_HO ||_F < 1e-12 (Moshinsky-Talmi N_tot).

    Paper 24 corollary, operator-commutativity leg: the HO kinetic
    Hamiltonian commutes with any central two-body interaction because
    V is block-diagonal in the total HO quantum number N_tot.
    """
    data = _build(N_max)
    H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
    H_vee = data['H_vee_full']

    # The one-body off-diagonal block is identically zero in the
    # Bargmann eigenbasis (pure HO is diagonal).
    assert np.linalg.norm(data['H_h1_offdiag']) < 1e-12, \
        'H_h1_offdiag should be identically zero for the pure HO'

    C = H_kin @ H_vee - H_vee @ H_kin
    rel_norm = np.linalg.norm(C) / np.linalg.norm(H_kin)
    # Paper 24 states < 1e-15; observed ~2.4e-16 to 2.7e-16. Assert the
    # paper-quoted ceiling-with-headroom 1e-12.
    assert rel_norm < 1e-12, \
        f'|| [H_HO, V] ||_F / || H_HO ||_F = {rel_norm:.3e} ' \
        f'exceeds noise floor at N_max={N_max}'


def test_paper24_ho_gs_energy_basis_independent():
    """The GS energy is identical at N_max=2 and N_max=3.

    Direct consequence of the GS living entirely in the
    one-dimensional N_tot=0 sector: enlarging the basis adds only
    higher-N_tot configurations that do not couple to the ground state.
    """
    e2 = np.linalg.eigvalsh(_build(2)['H_full'])[0]
    e3 = np.linalg.eigvalsh(_build(3)['H_full'])[0]
    assert e2 == pytest.approx(e3, abs=1e-9), \
        f'GS energy basis-dependent: N_max=2 -> {e2:.6f}, ' \
        f'N_max=3 -> {e3:.6f} MeV'
