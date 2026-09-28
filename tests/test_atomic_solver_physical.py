"""AtomicSolver produces the physical hydrogen-like spectrum by default.

Added 2026-09-28 (PI-approved).  Before this date ``AtomicSolver(max_n, Z)``
built only H = kappa * Z^2 * (D - A), whose low spectrum is a dense cluster
near -Z^2/2 rather than the Rydberg ladder.  The default is now
``hamiltonian='exact'``: H = diag(-Z^2/(2 n^2)), the one-electron Coulomb
Hamiltonian T + V_ne in the hydrogenic (k = Z) basis, which is exactly diagonal.
The graph model survives as ``hamiltonian='graph'``.

Wrong answers these tests reject:
  * the default reverting to the graph model (spectrum -0.4936, -0.4936, ...
    instead of -1/2, -1/8 x4, ...);
  * a wrong Z scaling (e.g. -Z/(2n^2)) or wrong degeneracy (2l+1 or n instead
    of n^2) in the exact model;
  * the legacy graph model silently changing (regression pin at max_n = 20).
"""
from collections import Counter
from fractions import Fraction

import numpy as np
import pytest

from geovac.atomic_solver import AtomicSolver, solve_atom


def _expected_multiset(Z: int, max_n: int) -> Counter:
    """Exact rational levels -Z^2/(2n^2), each with multiplicity n^2."""
    out: Counter = Counter()
    for n in range(1, max_n + 1):
        out[Fraction(-Z * Z, 2 * n * n)] += n * n
    return out


def _as_multiset(values: np.ndarray, max_n: int) -> Counter:
    """Snap floats to the nearest candidate level (tolerance 1e-12)."""
    out: Counter = Counter()
    for v in values:
        f = Fraction(float(v)).limit_denominator(2 * max_n * max_n)
        assert abs(float(f) - float(v)) < 1e-12, (v, f)
        out[f] += 1
    return out


@pytest.mark.parametrize("Z", [1, 2, 3])
@pytest.mark.parametrize("max_n", [1, 2, 3, 4, 5, 6])
def test_default_is_hydrogenic_spectrum(Z: int, max_n: int) -> None:
    solver = AtomicSolver(max_n, Z)               # no kwarg: the default
    assert solver.hamiltonian == 'exact'
    N = solver.n_states
    assert N == max_n * (max_n + 1) * (2 * max_n + 1) // 6

    expected = _expected_multiset(Z, max_n)
    # (a) full dense spectrum of the matrix actually built
    dense = np.linalg.eigvalsh(solver.H.toarray())
    assert _as_multiset(dense, max_n) == expected
    # (b) the solver's own eigen-API returns the same, ascending
    e, psi = solver.compute_ground_state(n_states=N)
    assert np.all(np.diff(e) >= 0)
    assert _as_multiset(e, max_n) == expected
    np.testing.assert_allclose(solver.H @ psi, psi * e, atol=1e-14)
    np.testing.assert_allclose(psi.T @ psi, np.eye(N), atol=1e-14)
    # (c) each eigenvector lives on a single shell n with E = -Z^2/(2n^2)
    for k in range(N):
        n = solver.lattice.states[int(np.argmax(np.abs(psi[:, k])))][0]
        assert e[k] == pytest.approx(-Z * Z / (2.0 * n * n), abs=1e-14)


@pytest.mark.parametrize("Z", [1, 2, 3])
def test_default_ground_state_and_first_excited(Z: int) -> None:
    e, _ = AtomicSolver(6, Z).compute_ground_state(n_states=5)
    np.testing.assert_allclose(
        e, [-Z * Z / 2.0] + [-Z * Z / 8.0] * 4, rtol=0, atol=1e-14)
    E0, _ = solve_atom(Z=Z, max_n=6)
    assert E0 == pytest.approx(-Z * Z / 2.0, abs=1e-14)


def test_graph_model_legacy_cluster_pinned() -> None:
    """Regression pin for the LABELLED legacy model (not physical)."""
    solver = AtomicSolver(20, 1, hamiltonian='graph')
    e = np.sort(solver.compute_ground_state(n_states=4)[0])
    np.testing.assert_allclose(
        e, [-0.49363618, -0.49361237, -0.49323372, -0.49323178],
        rtol=0, atol=5e-8)
    # the physical first excited level -1/8 is nowhere near
    assert abs(e[1] - (-0.125)) > 0.3


def test_invalid_arguments_raise() -> None:
    with pytest.raises(ValueError):
        AtomicSolver(3, 1, hamiltonian='laplacian')
    with pytest.raises(ValueError):
        AtomicSolver(3, 1, kinetic_scale=-0.1)    # scale is graph-only
    # graph mode still accepts a custom scale
    AtomicSolver(3, 1, kinetic_scale=-0.1, hamiltonian='graph')
