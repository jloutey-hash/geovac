"""
ARCHIVED 2026-05-23 (Cleanup Track B): TestMolecularEncoding class extracted from
tests/test_qubit_encoding.py. MolecularLatticeIndex was the LCAO molecular
encoding (CLAUDE.md §3 documented dead-end: 'Graph Laplacian kinetic energy is
R-independent; need natural geometry where separation occurs'). It was removed
from geovac/lattice_index.py in v2.7.0 (commit 8d692a0) as part of the LCAO →
natural geometry transition documented in Paper FCI-M (29-version diagnostic arc,
v0.9.8-0.9.37).

The molecular qubit encoding equivalent in current code is geovac/composed_qubit.py
(separately tested in tests/test_composed_qubit.py and friends), which uses the
composed natural-geometry architecture (Paper 17) instead of the dead-end LCAO
approach.
"""

import warnings
import pytest
import numpy as np

from openfermion import jordan_wigner, get_sparse_operator

from geovac.lattice_index import LatticeIndex, MolecularLatticeIndex
from geovac.qubit_encoding import (
    JordanWignerEncoder,
    PauliAnalysis,
    build_fermion_op_from_integrals,
    fit_pauli_scaling,
)

warnings.filterwarnings('ignore', category=UserWarning)


class TestMolecularEncoding:
    """Tests for MolecularLatticeIndex encoding."""

    @pytest.fixture(scope="class")
    def h2_mol(self) -> "MolecularLatticeIndex":
        """H2 molecule at R=1.4 bohr, nmax=2."""
        return MolecularLatticeIndex(
            Z_A=1, Z_B=1, nmax_A=2, nmax_B=2,
            R=1.4, n_electrons=2, n_bridges=10,
            vee_method='slater_full',
        )

    def test_h2_encoding(self, h2_mol: "MolecularLatticeIndex") -> None:
        """H2 molecule produces valid qubit encoding."""
        enc = JordanWignerEncoder(h2_mol)
        analysis = enc.analyze()
        assert analysis.n_pauli_terms > 0
        assert analysis.n_qubits == h2_mol.n_sp

    def test_h2_nuclear_repulsion(self, h2_mol: "MolecularLatticeIndex") -> None:
        """Nuclear repulsion is automatically extracted."""
        enc = JordanWignerEncoder(h2_mol)
        assert enc._v_nn == pytest.approx(1.0 / 1.4, abs=1e-10)
