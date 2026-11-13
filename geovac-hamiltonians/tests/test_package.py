"""Smoke tests for the geovac-hamiltonians front door (0.3.0+).

The package re-exports geovac.ecosystem_export; counts are pinned in the main repository's tests
(tests/test_geovac_hamiltonians_shim.py and the library tests), not duplicated here.
"""
import warnings


def test_import_and_version():
    import geovac_hamiltonians
    from geovac_hamiltonians import GeoVacHamiltonian, hamiltonian
    assert callable(hamiltonian) and GeoVacHamiltonian is not None
    assert geovac_hamiltonians.__version__ == "0.3.0"


def test_reexports_library():
    import geovac.ecosystem_export as ee
    from geovac_hamiltonians import hamiltonian
    assert hamiltonian is ee.hamiltonian


def test_lih_builds_and_exports():
    from openfermion import QubitOperator
    from geovac_hamiltonians import hamiltonian
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        h = hamiltonian("LiH")
    assert h.n_qubits == 30 and h.n_terms == 838
    op = h.to_openfermion()
    assert isinstance(op, QubitOperator) and len(op.terms) == 838
