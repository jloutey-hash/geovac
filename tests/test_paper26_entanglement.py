"""Verification tests for Paper 26 (Entanglement Structure of the Angular
Momentum Eigenbasis).

Backs the four load-bearing headlines of Paper 26 against the reference
driver debug/archive/misc/energy_entanglement_decoupling.py (the same module
Paper 27's tests import).  Prior to this file the driver COMPUTED these
quantities but no test ASSERTED them (a /qa group6 first-cert coverage gap):

1. Energy-entanglement decoupling (abstract + Sec. II, tab:decoupling):
   the off-diagonal one-body Hamiltonian (graph Laplacian, kappa=-1/16)
   contributes 10-39% of the correlation energy but < 0.2% of the
   entanglement entropy; the off-diagonal V_ee contributes ~100%.
2. Entanglement entropy scaling S ~ Z^{-2.56} across He-like ions.
3. Basis-intrinsic sparsity (Sec. III): the angular-momentum eigenbasis
   gives ERI density 42.4% (265/625) at n_max=2; any rotation fills to
   99.2%; the nonzero ERI count is Z-independent.

The entropy reported is the von Neumann entropy of the normalized one-body
reduced density matrix (natural-orbital occupation entropy), matching the
paper's corrected Sec. II definition.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from debug.archive.misc.energy_entanglement_decoupling import (  # noqa: E402
    part1_isoelectronic_decoupling,
    build_orbital_integrals,
    rotate_integrals,
    count_significant,
    scipy_matrix_exp,
)

_Z_VALUES = [2, 3, 4, 5, 6, 8, 10]


@pytest.fixture(scope="module")
def decoupling():
    """Run the n_max=4 isoelectronic decoupling once (~2 s) and share it."""
    return part1_isoelectronic_decoupling(n_max=4)


# ---------------------------------------------------------------------------
# Headline 1 -- energy-entanglement decoupling (the lead result)
# ---------------------------------------------------------------------------

def test_paper26_energy_entanglement_decoupling(decoupling):
    """h1_offdiag: 10-39% of E_corr but <0.2% of S; V_ee_offdiag: ~100% of S."""
    frac_E_h1 = [decoupling[str(z)]['frac_E_h1_offdiag'] for z in _Z_VALUES]
    frac_S_h1 = [decoupling[str(z)]['frac_S_h1_offdiag'] for z in _Z_VALUES]
    frac_S_vee = [decoupling[str(z)]['frac_S_vee_offdiag'] for z in _Z_VALUES]

    # Energy: the graph off-diagonal carries a 10-39% band of the correlation
    # energy (paper abstract).  Genuine test -- the band is the claim.
    assert min(frac_E_h1) == pytest.approx(0.10, abs=0.02), \
        f'min h1 E-fraction {min(frac_E_h1):.3f} not ~10%'
    assert max(frac_E_h1) == pytest.approx(0.39, abs=0.02), \
        f'max h1 E-fraction {max(frac_E_h1):.3f} not ~39%'
    assert all(0.09 <= f <= 0.40 for f in frac_E_h1)

    # Entanglement: the SAME graph off-diagonal is entanglement-inert (<0.2%).
    assert max(abs(f) for f in frac_S_h1) < 2.0e-3, \
        f'h1 S-fraction {max(abs(f) for f in frac_S_h1):.2e} not <0.2%'

    # V_ee off-diagonal carries essentially all the entanglement (~100%).
    assert min(frac_S_vee) > 0.99, \
        f'V_ee S-fraction {min(frac_S_vee):.4f} not ~100%'

    # The decoupling is the >100x separation between the two channels:
    # h1 moves energy but not entanglement; V_ee moves entanglement.
    for z in _Z_VALUES:
        f_E = abs(decoupling[str(z)]['frac_E_h1_offdiag'])
        f_S = abs(decoupling[str(z)]['frac_S_h1_offdiag'])
        assert f_E / max(f_S, 1e-12) > 50.0, f'no decoupling at Z={z}'


# ---------------------------------------------------------------------------
# Headline 2 -- entanglement entropy scaling S ~ Z^{-2.56}
# ---------------------------------------------------------------------------

def test_paper26_entanglement_entropy_z_scaling(decoupling):
    """von Neumann entropy of the normalized 1-RDM scales as S ~ Z^{-2.56}."""
    S = np.array([decoupling[str(z)]['S_full'] for z in _Z_VALUES])
    logZ = np.log(np.array(_Z_VALUES, dtype=float))
    slope = np.polyfit(logZ, np.log(S), 1)[0]
    # Paper: S ~ Z^{-2.56}.  Reproduces -2.563.
    assert slope == pytest.approx(-2.56, abs=0.10), \
        f'entropy scaling exponent {slope:.3f} not ~ -2.56'
    # Monotone decrease (a genuine power-law, not noise).
    assert np.all(np.diff(S) < 0)


# ---------------------------------------------------------------------------
# Headline 3 -- basis-intrinsic sparsity: 42.4% -> 99.2% step, Z-independent
# ---------------------------------------------------------------------------

def test_paper26_basis_intrinsic_sparsity_step_function():
    """Graph-basis ERI density 42.4% (265/625); any rotation fills to >99%."""
    h1, eri, orbitals, n_spatial = build_orbital_integrals(2.0, n_max=2)
    total_eri = n_spatial ** 4
    n_graph = count_significant(eri)

    # Graph (angular-momentum) basis: 265/625 = 42.4% density (paper).
    assert n_graph == 265, f'graph ERI count {n_graph} != 265'
    assert n_graph / total_eri == pytest.approx(0.424, abs=0.01)

    # Any rotation off the identity fills the tensor to near-complete density.
    rng = np.random.RandomState(42)
    G = rng.randn(n_spatial, n_spatial)
    G = (G - G.T) / 2.0
    G /= np.linalg.norm(G)
    U = scipy_matrix_exp(0.1 * G)  # a small but nonzero rotation
    _, eri_rot = rotate_integrals(h1, eri, U)
    dens_rot = count_significant(eri_rot) / total_eri
    assert dens_rot > 0.99, f'rotated density {dens_rot:.3f} not >99%'
    # The transition is a step at the identity: sparse -> dense.
    assert dens_rot - n_graph / total_eri > 0.5


def test_paper26_eri_count_is_z_independent():
    """The nonzero ERI count is set by geometry, not by Z (Z-independent)."""
    counts = []
    for Z in (2.0, 6.0, 10.0):
        _, eri, _, _ = build_orbital_integrals(Z, n_max=2)
        counts.append(count_significant(eri))
    assert len(set(counts)) == 1, f'ERI count Z-dependent: {counts}'
    assert counts[0] == 265


# ---------------------------------------------------------------------------
# Headline 4 -- core-valence decoupling (the composed-factorization justification)
# ---------------------------------------------------------------------------

def test_paper26_core_valence_decoupling_thresholds():
    """Core-valence mutual information: ~2e-3 by Z=4 (Be), <1e-3 for Z>=5 (B),
    and O(1) for He/Li (no core-valence separation) -- the quantitative
    justification for the composed fiber-bundle factorization (Paper 26 Sec V,
    Table II). Recomputed from the FCI ground state via run_atom."""
    from debug.archive.misc.entanglement_first_row import run_atom
    # run_atom returns (results, MI, s_single, orbitals); results[0] carries I_cv.
    # n_max=2 matches the paper's Table II source (debug/data/entanglement_first_row.json).
    def icv(Z, ne):
        return run_atom(Z, ne, 2, 'X')[0]['I_core_valence']
    # Be (Z=4): the threshold datum, I_cv ~ 2e-3 (small but not yet <1e-3).
    be = icv(4, 4)
    assert be == pytest.approx(0.0022, abs=5e-4), f'Be I_cv = {be:.4e}'
    assert be < 1e-2  # already decoupled by Z=4
    # B (Z=5): first row below 1e-3.
    b = icv(5, 5)
    assert abs(b) < 1e-3, f'B I_cv = {b:.4e} not < 1e-3'
    # Li (Z=3): still strongly coupled (no clean core-valence split yet).
    li = icv(3, 3)
    assert li > 0.1, f'Li I_cv = {li:.4e} should be O(0.1)'
