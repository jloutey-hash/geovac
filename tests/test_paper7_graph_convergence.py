"""Paper 7, "New contributions" item 1 -- discrete graph convergence [MEASURED].

Paper 7 states that the GeoVac paraboloid lattice (nodes (n, l, m), edges
from L+- and T+-) converges numerically, through n_max = 30 and with no
rate proven, to the Laplace-Beltrami operator on unit S^3.  Until
2026-09-02 no test in tests/ built the n_max = 30 spectrum for that claim
(claim_test_matrix row "NO-TEST", trunk DELTA #1 CODE-A).  The s/p-splitting
half of the convergence story is covered by tests/test_trunk_qa_splitting.py
(Paper 1); this file covers the other half -- the energy the production
Hamiltonian actually converges to.

What is measured (production objects only: AtomicSolver builds
H = kappa Z^2 (D - A) on GeometricLattice with kappa = -1/16):

    n_max :   5        8        10       15       20       30
    E_0   : -0.41363 -0.46375 -0.47571 -0.48883 -0.49364 -0.49713
    error :  17.27%    7.25%    4.86%    2.23%    1.27%    0.57%

against the S^3 n = 1 harmonic's Fock image E = -1/2 (the hydrogen ground
state).  Measured 2026-09-02; the sequence is monotone, one-sided (never
below -1/2), and sub-percent by n_max = 30.  NO rate is asserted, matching
the paper's wording.

Honest scope note (recorded, not asserted): the bottom of the graph
spectrum is DENSE -- at n_max = 30 the six lowest eigenvalues lie within
0.06% of each other -- so the graph spectrum is not the Rydberg ladder
level by level (the -(n^2 - 1) ladder is the continuum S^3 object, not the
graph's; see Paper 1 and the kappa memory note).  The convergence this
file pins is the ground-state energy, which is what "converges to the
Laplace-Beltrami operator" is measured by at the level of the spectrum's
edge.  Tiers: [MEASURED] against the exact -1/2.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse.linalg import eigsh

from geovac.atomic_solver import AtomicSolver

E_EXACT = -0.5
GRID = [5, 8, 10, 15, 20, 30]

# Measured 2026-09-02 (this file's docstring); used as pins with margin.
MEASURED_E0 = {5: -0.41363, 8: -0.46375, 10: -0.47571,
               15: -0.48883, 20: -0.49364, 30: -0.49713}


def _ground_energy(max_n: int) -> float:
    solver = AtomicSolver(max_n, Z=1)
    E, _ = eigsh(solver.H, k=1, which="SA")
    return float(E[0])


@pytest.fixture(scope="module")
def ground_energies() -> dict:
    return {n: _ground_energy(n) for n in GRID}


def test_ground_energy_reproduces_measured_values(ground_energies):
    """Pins the measured sequence (5e-4 absolute) so a lattice or kappa
    change is visible here, not only in the endpoint."""
    for n in GRID:
        assert abs(ground_energies[n] - MEASURED_E0[n]) < 5e-4, (
            f"n_max={n}: E0={ground_energies[n]:.5f}, "
            f"measured 2026-09-02 {MEASURED_E0[n]:.5f}")


def test_ground_energy_error_decreases_monotonically(ground_energies):
    """The error against -1/2 shrinks at every step of the grid (no rate
    asserted -- Paper 7 proves none)."""
    errs = [abs(ground_energies[n] - E_EXACT) for n in GRID]
    for (n_a, e_a), (n_b, e_b) in zip(zip(GRID, errs), zip(GRID[1:], errs[1:])):
        assert e_b < e_a, (
            f"error grew from n_max={n_a} ({e_a:.5f}) to n_max={n_b} ({e_b:.5f})")


def test_ground_energy_never_overshoots(ground_energies):
    """One-sided approach: E0 > -1/2 at every n_max on the grid."""
    for n in GRID:
        assert ground_energies[n] > E_EXACT, (
            f"n_max={n}: E0={ground_energies[n]:.5f} below -1/2")


def test_ground_energy_subpercent_at_nmax_30(ground_energies):
    """The n_max = 30 endpoint Paper 7 names: sub-percent, and a two-sided
    window around the measured 0.57% so the assertion can fail in both
    directions (a one-sided '< 1%' could not detect a silent improvement
    that would make the paper's wording stale)."""
    err30 = abs(ground_energies[30] - E_EXACT) / abs(E_EXACT)
    assert 0.004 < err30 < 0.0075, f"n_max=30 relative error {err30*100:.3f}%"


def test_bottom_of_spectrum_is_dense_at_nmax_30():
    """Scope guard for the docstring's honest note: the six lowest
    eigenvalues at n_max = 30 lie within 0.1% of each other, so no test in
    this file may be read as a level-by-level Rydberg match."""
    solver = AtomicSolver(30, Z=1)
    E, _ = eigsh(solver.H, k=6, which="SA")
    E = np.sort(E)
    spread = (E[-1] - E[0]) / abs(E[0])
    assert spread < 1e-3, f"lowest-six spread {spread*100:.4f}%"
