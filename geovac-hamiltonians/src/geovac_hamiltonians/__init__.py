"""
geovac-hamiltonians: Structurally sparse qubit Hamiltonians from geometric quantum chemistry
=============================================================================================

Provides pre-built qubit Hamiltonians for molecular systems using the GeoVac
composed natural geometry framework.  Hamiltonians are returned as OpenFermion
``QubitOperator`` objects wrapped in ``GeoVacHamiltonian``, which additionally
exports to Qiskit and PennyLane.

Quick start::

    from geovac_hamiltonians import hamiltonian

    H = hamiltonian('LiH')
    print(H.n_qubits, H.n_terms)   # 30, 334

    # Export to your framework of choice
    of_op  = H.to_openfermion()
    qk_op  = H.to_qiskit()       # requires qiskit
    pl_op  = H.to_pennylane()     # requires pennylane

Supported systems: H2, He, LiH, BeH2, H2O.
"""

__version__ = "0.2.0"

from geovac_hamiltonians._ecosystem_export import (
    GeoVacHamiltonian,
    hamiltonian,
)
from geovac_hamiltonians._pk_partitioning import pk_classical_energy

__all__ = [
    "GeoVacHamiltonian",
    "hamiltonian",
    "pk_classical_energy",
    "__version__",
]
