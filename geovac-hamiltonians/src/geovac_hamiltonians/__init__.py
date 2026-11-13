"""
geovac-hamiltonians: qubit Hamiltonians from the GeoVac library
===============================================================

A thin front door to the main ``geovac`` package.  ``hamiltonian`` and
``GeoVacHamiltonian`` ARE ``geovac.ecosystem_export``'s objects, re-exported,
so this package cannot drift from the library the papers measure.

(Until 0.3.0 this package carried its own copy of the builder modules, last
synced at v5.20.0 and still on the retired pair-diagonal ERI rule: it returned
LiH 334 terms where the library gives 838.  CHANGELOG v6.7.24.)

Quick start::

    from geovac_hamiltonians import hamiltonian

    H = hamiltonian('LiH')
    print(H.n_qubits, H.n_terms)

    of_op  = H.to_openfermion()
    qk_op  = H.to_qiskit()       # requires qiskit
    pl_op  = H.to_pennylane()    # requires pennylane

Supported systems: the 37-system library (``geovac.ecosystem_export``).
"""

__version__ = "0.3.0"

from geovac.ecosystem_export import GeoVacHamiltonian, hamiltonian
from geovac.pk_partitioning import pk_classical_energy

__all__ = [
    "GeoVacHamiltonian",
    "hamiltonian",
    "pk_classical_energy",
    "__version__",
]
