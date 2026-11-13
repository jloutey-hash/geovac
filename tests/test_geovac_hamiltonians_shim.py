"""Guard: the geovac-hamiltonians package is a front door to geovac, not a copy of it.

Until 0.3.0 the package vendored its own builder modules, last synced at v5.20.0 and still on the retired
pair-diagonal ERI rule (LiH 334 terms vs the library's 838; CHANGELOG v6.7.24). Each test names the wrong
answer it rejects.
"""
import os
import sys
import warnings

import pytest

_SRC = os.path.join(os.path.dirname(__file__), "..", "geovac-hamiltonians", "src")


@pytest.fixture(scope="module")
def pkg():
    sys.path.insert(0, os.path.abspath(_SRC))
    try:
        import geovac_hamiltonians
        yield geovac_hamiltonians
    finally:
        sys.path.remove(os.path.abspath(_SRC))


def test_hamiltonian_is_the_library_function(pkg):
    """Rejects a re-vendored builder: any copy, however fresh, is a different function object."""
    import geovac.ecosystem_export as ee
    assert pkg.hamiltonian is ee.hamiltonian
    assert pkg.GeoVacHamiltonian is ee.GeoVacHamiltonian


def test_no_vendored_modules(pkg):
    """Rejects the 0.2.0 layout: the package directory holds only __init__.py."""
    d = os.path.dirname(pkg.__file__)
    assert sorted(f for f in os.listdir(d) if f.endswith(".py")) == ["__init__.py"]


def test_lih_exact_rule_count(pkg):
    """LiH through the package: 838 terms incl. identity (exact global-M_L rule). Rejects the retired 334."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        h = pkg.hamiltonian("LiH")
    assert (h.n_qubits, h.n_terms) == (30, 838)
