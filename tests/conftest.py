"""
Shared pytest fixtures for GeoVac test suite.
"""

import sys
import os
import pytest

# Ensure project root is on path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


@pytest.fixture(scope="session")
def hydrogen_solver():
    """Shared hydrogen AtomicSolver (max_n=10) for reuse across tests."""
    from geovac import AtomicSolver, UNIVERSAL_KINETIC_SCALE
    return AtomicSolver(max_n=10, Z=1, kinetic_scale=UNIVERSAL_KINETIC_SCALE)


@pytest.fixture(scope="session")
def helium_solver():
    """Shared helium MoleculeHamiltonian (max_n=5) for reuse across tests."""
    from geovac import MoleculeHamiltonian
    mol = MoleculeHamiltonian(
        nuclei=[(0.0, 0.0, 0.0)],
        nuclear_charges=[2],
        max_n=5,
    )
    return mol
