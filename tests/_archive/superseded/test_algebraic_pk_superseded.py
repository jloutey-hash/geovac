"""
Tests for the algebraic Phillips-Kleinman projector.

Validates:
  1. Core eigenvector persistence in CoreScreening
  2. Algebraic projector shape and content from AbInitioPK
  3. Backward compatibility of existing Gaussian PK modes
  4. Algebraic PK smoke test (bound PES)
  5. Projector rank (must be 1)
  6. Channel coupling (off-diagonal elements in (0,0) block)
"""

import numpy as np
import pytest

from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK
from geovac.composed_diatomic import ComposedDiatomicSolver
from geovac.level4_multichannel import (
    build_angular_hamiltonian,
    _channel_list,
)


# ---------- Fixture: solved Li core ----------

@pytest.fixture(scope="module")
def li_core() -> CoreScreening:
    """Solve Li core once for all tests in this module."""
    cs = CoreScreening(Z=3, l_max=2, n_alpha=100)
    cs.solve(verbose=False)
    return cs


@pytest.fixture(scope="module")
def li_pk(li_core) -> AbInitioPK:
    """AbInitioPK from solved Li core."""
    return AbInitioPK(li_core, n_core=2)


# ---------- Test 1: Core eigenvector persistence ----------

def test_core_eigenvector_persisted(li_core):
    """CoreScreening.solve() must store core_l0_wavefunction."""
    assert li_core.core_l0_wavefunction is not None
    assert li_core.core_eigenvector is not None
    assert li_core.core_n_alpha is not None
    assert li_core.core_h_alpha is not None
    assert li_core.core_R_representative is not None


def test_core_l0_normalized(li_core):
    """The l=0 core wavefunction must be normalized to 1."""
    l0 = li_core.core_l0_wavefunction
    h = li_core.core_h_alpha
    norm_sq = np.sum(l0**2) * h
    assert abs(norm_sq - 1.0) < 1e-6, f"norm^2 = {norm_sq}, expected 1.0"


def test_core_l0_shape(li_core):
    """l=0 slice has shape (n_alpha,)."""
    assert li_core.core_l0_wavefunction.shape == (li_core.core_n_alpha,)


def test_core_eigenvector_shape(li_core):
    """Full eigenvector has shape ((l_max+1) * n_alpha,)."""
    n_l = li_core.l_max + 1
    expected = n_l * li_core.core_n_alpha
    assert li_core.core_eigenvector.shape == (expected,)


# ---------- Test 2: Projector shape ----------

def test_algebraic_projector_keys(li_pk):
    """algebraic_projector() returns dict with all required keys."""
    proj = li_pk.algebraic_projector(atom='A')
    required_keys = {
        'mode', 'core_l0_wavefunction', 'core_n_alpha',
        'core_h_alpha', 'energy_shift', 'E_core_per_electron',
        'E_val', 'atom',
    }
    assert required_keys.issubset(proj.keys())


def test_algebraic_projector_mode(li_pk):
    """Mode must be 'algebraic'."""
    proj = li_pk.algebraic_projector()
    assert proj['mode'] == 'algebraic'


def test_algebraic_projector_energy_shift_positive(li_pk):
    """Energy shift must be positive (repulsive barrier)."""
    proj = li_pk.algebraic_projector()
    assert proj['energy_shift'] > 0


def test_algebraic_projector_wavefunction_shape(li_pk):
    """core_l0_wavefunction shape matches core_n_alpha."""
    proj = li_pk.algebraic_projector()
    assert proj['core_l0_wavefunction'].shape == (proj['core_n_alpha'],)


# ---------- Test 3: Backward compatibility ----------

def test_gaussian_pk_unchanged():
    """Existing Gaussian PK (ab_initio) must still produce reasonable R_eq."""
    solver = ComposedDiatomicSolver.LiH_ab_initio(l_max=2, verbose=False)
    solver.solve_core()
    # Just verify the Gaussian PK parameters are set (full PES is too slow
    # for a unit test — the existing test_composed_diatomic.py covers that).
    assert solver.pk_A > 0
    assert solver.pk_B > 0
    assert solver.pk_potentials is not None
    assert len(solver.pk_potentials) == 1
    assert 'C_core' in solver.pk_potentials[0]
    assert 'beta_core' in solver.pk_potentials[0]


# ---------- Test 4: Algebraic PK smoke test ----------

def test_algebraic_pk_solver_runs():
    """LiH_algebraic_pk() must solve core without crashing."""
    solver = ComposedDiatomicSolver.LiH_algebraic_pk(l_max=2, verbose=False)
    solver.solve_core()

    assert solver.pk_projector is not None
    assert solver.pk_projector['mode'] == 'algebraic'
    assert solver.pk_projector['energy_shift'] > 0
    # No Gaussian pk_potentials in algebraic mode
    assert solver.pk_potentials is None


def test_algebraic_pk_single_R():
    """Algebraic PK must produce a finite electronic energy at a single R."""
    solver = ComposedDiatomicSolver.LiH_algebraic_pk(
        l_max=2, n_alpha=60, verbose=False,
    )
    solver.solve_core()
    E_elec = solver._solve_valence_at_R(R=3.0, n_Re=100)
    assert np.isfinite(E_elec), f"E_elec = {E_elec} is not finite"


# ---------- Test 5: Projector rank ----------

def test_projector_rank_one(li_pk):
    """The algebraic PK contribution to H must have rank 1."""
    proj = li_pk.algebraic_projector()

    # Build a small Level 4 angular Hamiltonian with only algebraic PK
    l_max = 2
    n_alpha = 30
    h = (np.pi / 2) / (n_alpha + 1)
    alpha_grid = (np.arange(n_alpha) + 1) * h

    # Heteronuclear channel list (LiH: Z_A=1 eff, Z_B=1)
    channels = _channel_list(l_max, homonuclear=False)
    n_ch = len(channels)
    N = n_ch * n_alpha

    # Build H with algebraic PK only (no other potentials)
    # We need the full H to isolate the PK contribution, so build with
    # and without the projector and take the difference.
    R_e = 2.0
    rho = 3.0 / (2.0 * R_e)

    H_with = build_angular_hamiltonian(
        alpha_grid, rho, R_e, l_max=l_max, Z_A=1.0, Z_B=1.0,
        pk_projector=proj,
    )
    H_without = build_angular_hamiltonian(
        alpha_grid, rho, R_e, l_max=l_max, Z_A=1.0, Z_B=1.0,
    )
    H_pk = H_with - H_without

    # Rank-1 matrix has exactly 1 nonzero eigenvalue
    evals = np.linalg.eigvalsh(H_pk)
    nonzero = np.abs(evals) > 1e-10 * np.max(np.abs(evals))
    rank = int(np.sum(nonzero))
    assert rank == 1, f"PK projector rank = {rank}, expected 1"


# ---------- Test 6: Channel coupling ----------

def test_projector_has_offdiag_in_00_block(li_pk):
    """The algebraic PK matrix must have off-diagonal elements within
    the (0,0) channel block (unlike diagonal-only Gaussian PK)."""
    proj = li_pk.algebraic_projector()

    l_max = 2
    n_alpha = 30
    h = (np.pi / 2) / (n_alpha + 1)
    alpha_grid = (np.arange(n_alpha) + 1) * h

    R_e = 2.0
    rho = 3.0 / (2.0 * R_e)

    H_with = build_angular_hamiltonian(
        alpha_grid, rho, R_e, l_max=l_max, Z_A=1.0, Z_B=1.0,
        pk_projector=proj,
    )
    H_without = build_angular_hamiltonian(
        alpha_grid, rho, R_e, l_max=l_max, Z_A=1.0, Z_B=1.0,
    )
    H_pk = H_with - H_without

    # (0,0) channel is always channel 0 in sorted order
    block_00 = H_pk[0:n_alpha, 0:n_alpha]

    # Extract off-diagonal elements
    offdiag = block_00 - np.diag(np.diag(block_00))
    max_offdiag = np.max(np.abs(offdiag))
    assert max_offdiag > 1e-10, (
        f"PK (0,0) block has no off-diagonal elements "
        f"(max = {max_offdiag})"
    )
