"""
Tests for Numba-accelerated FCI Hamiltonian assembly.

Validates that the Numba kernel reproduces the pure-Python DirectCISolver
to machine precision across multiple systems.
"""

import time
import warnings

import numpy as np
import pytest
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

from geovac.lattice_index import LatticeIndex
from geovac.direct_ci import DirectCISolver

try:
    from geovac.direct_ci_numba import (
        NUMBA_AVAILABLE,
        assemble_hamiltonian_numba,
        build_colex_to_lex,
        sd_basis_to_array,
        warmup_jit,
    )
except ImportError:
    NUMBA_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not NUMBA_AVAILABLE, reason="Numba not installed"
)


def _build_numba_H(idx: LatticeIndex, solver: DirectCISolver) -> csr_matrix:
    """Assemble Hamiltonian via Numba and return as CSR."""
    diag_idx, diag_val, off_row, off_col, off_val = \
        assemble_hamiltonian_numba(
            idx.sd_basis, solver._H1_dense, solver._eri_4d,
            solver._h1_diag_arr, solver.n_sp, idx.threshold,
            solver._spatial_targets, solver.n_spatial,
        )
    n_sd = solver.n_sd
    H_diag = csr_matrix((diag_val, (diag_idx, diag_idx)), shape=(n_sd, n_sd))
    H_upper = csr_matrix((off_val, (off_row, off_col)), shape=(n_sd, n_sd))
    return (H_upper + H_upper.T + H_diag).tocsr()


def _solve_and_compare(n_el: int, max_n: int, Z: int, name: str) -> None:
    """Compare Python vs Numba assembly and assert exact match."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        idx = LatticeIndex(
            n_electrons=n_el, max_n=max_n, nuclear_charge=Z,
            vee_method='slater_full', fci_method='direct',
        )

    solver = DirectCISolver(idx)
    H_py = solver._assemble_python()
    H_nb = _build_numba_H(idx, solver)

    max_diff = abs(H_py - H_nb).max()
    assert max_diff < 1e-12, (
        f"{name}: max matrix element diff = {max_diff:.2e}"
    )

    rng = np.random.RandomState(42)
    v0 = rng.randn(solver.n_sd)
    E_py, _ = eigsh(H_py, k=1, which="SA", v0=v0)
    E_nb, _ = eigsh(H_nb, k=1, which="SA", v0=v0)

    assert abs(E_py[0] - E_nb[0]) < 1e-10, (
        f"{name}: E_py={E_py[0]:.10f}, E_nb={E_nb[0]:.10f}, "
        f"diff={abs(E_py[0]-E_nb[0]):.2e}"
    )


class TestCombinatoricalIndex:
    """Verify SD→index mapping via combinatorial number system."""

    def test_colex_mapping_he(self) -> None:
        """Colex→lex mapping is valid permutation for He nmax=3."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            idx = LatticeIndex(
                n_electrons=2, max_n=3, nuclear_charge=2,
                vee_method='slater_full', fci_method='direct',
            )
        c2l = build_colex_to_lex(idx.sd_basis, idx.n_sp)
        assert len(c2l) == idx.n_sd
        assert np.all(c2l >= 0)
        assert len(set(c2l)) == idx.n_sd  # all unique

    def test_colex_mapping_li(self) -> None:
        """Colex→lex mapping for Li nmax=3 (3-electron)."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            idx = LatticeIndex(
                n_electrons=3, max_n=3, nuclear_charge=3,
                vee_method='slater_full', fci_method='direct',
            )
        c2l = build_colex_to_lex(idx.sd_basis, idx.n_sp)
        assert len(c2l) == idx.n_sd
        assert np.all(c2l >= 0)


class TestNumbaVsPython:
    """Numba assembly must match Python to machine precision."""

    def test_he_nmax2(self) -> None:
        """He nmax=2: exact matrix match."""
        _solve_and_compare(2, 2, 2, "He nmax=2")

    def test_he_nmax3(self) -> None:
        """He nmax=3: exact matrix match."""
        _solve_and_compare(2, 3, 2, "He nmax=3")

    def test_li_nmax3(self) -> None:
        """Li nmax=3 (3-electron): exact matrix match."""
        _solve_and_compare(3, 3, 3, "Li nmax=3")

    def test_be_nmax3(self) -> None:
        """Be nmax=3 (4-electron, 20k SDs): exact matrix match."""
        _solve_and_compare(4, 3, 4, "Be nmax=3")


class TestNumbaSpeedup:
    """Verify Numba provides meaningful speedup."""

    def test_li_speedup(self) -> None:
        """Li nmax=3: Numba at least 5x faster than Python."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            idx = LatticeIndex(
                n_electrons=3, max_n=3, nuclear_charge=3,
                vee_method='slater_full', fci_method='direct',
            )

        solver = DirectCISolver(idx)

        # Python timing
        t0 = time.perf_counter()
        solver._assemble_python()
        dt_py = time.perf_counter() - t0

        # Numba timing (kernel already warm)
        t0 = time.perf_counter()
        _build_numba_H(idx, solver)
        dt_nb = time.perf_counter() - t0

        speedup = dt_py / dt_nb
        print(f"  Li nmax=3 speedup: {speedup:.1f}x "
              f"(py={dt_py:.3f}s, nb={dt_nb:.3f}s)")
        assert speedup > 5.0, (
            f"Speedup only {speedup:.1f}x, expected > 5x"
        )


class TestIntegration:
    """Verify Numba is used automatically via DirectCISolver."""

    def test_auto_numba(self) -> None:
        """DirectCISolver.assemble_hamiltonian() uses Numba when available."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            idx = LatticeIndex(
                n_electrons=2, max_n=3, nuclear_charge=2,
                vee_method='slater_full', fci_method='direct',
            )
        solver = DirectCISolver(idx)

        # Should use Numba automatically
        H = solver.assemble_hamiltonian()
        rng = np.random.RandomState(42)
        v0 = rng.randn(solver.n_sd)
        E, _ = eigsh(H, k=1, which="SA", v0=v0)

        # Compare against Python path
        H_py = solver._assemble_python()
        E_py, _ = eigsh(H_py, k=1, which="SA", v0=v0)

        assert abs(E[0] - E_py[0]) < 1e-10
