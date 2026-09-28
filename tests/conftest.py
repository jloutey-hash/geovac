"""
Shared pytest fixtures for GeoVac test suite.
"""

import sys
import os
import pytest

# Ensure project root is on path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def pytest_addoption(parser):
    parser.addoption("--slow", action="store_true", default=False,
                     help="run slow tests")


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--slow"):
        skip_slow = pytest.mark.skip(reason="need --slow option to run")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)


@pytest.fixture(scope="session")
def hamiltonian_cache():
    """Session-scoped cache for ``ecosystem_export.hamiltonian()`` builds.

    Each unique ``(system, **kwargs)`` is built once and reused across tests
    in the same pytest session. The build is the expensive step (full
    balanced cross-V_ne multipole + JW transform); reusing it cuts the
    `test_ecosystem_export.py` default sweep from ~9 min to ~30 s when most
    tests touch the same handful of systems (LiH, BeH2, H2O, etc.).

    Usage from a test that previously called ``hamiltonian('LiH', ...)``:

        def test_foo(hamiltonian_cache):
            H = hamiltonian_cache('LiH')
            assert H.n_qubits > 0

    Kwargs are part of the cache key (e.g. ``hamiltonian_cache('LiH', max_n=3)``
    and ``hamiltonian_cache('LiH', max_n=4)`` are independent entries).
    ``verbose`` is dropped from the key because it never affects the build.

    Sprint Test-Slim 2026-06-07.
    """
    from geovac.ecosystem_export import hamiltonian as _hamiltonian
    cache = {}

    def get(system, **kwargs):
        kwargs.pop('verbose', None)
        key = (system, tuple(sorted(kwargs.items())))
        if key not in cache:
            cache[key] = _hamiltonian(system, verbose=False, **kwargs)
        return cache[key]

    return get


@pytest.fixture(scope="session")
def hydrogen_solver():
    """Shared hydrogen AtomicSolver (max_n=10) for reuse across tests."""
    from geovac import AtomicSolver, UNIVERSAL_KINETIC_SCALE
    return AtomicSolver(max_n=10, Z=1, kinetic_scale=UNIVERSAL_KINETIC_SCALE)


# (2026-09-26) The `helium_solver` fixture was removed: it had no consumers and
# built the archived MoleculeHamiltonian, whose atomic default kinetic scale
# -0.10298808 was FITTED to experiment. For helium use
# geovac.hyperspherical_radial.solve_helium (see tests/test_hyperspherical_he.py).
