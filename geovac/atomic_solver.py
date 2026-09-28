"""
Atomic Solver - single-electron hydrogen-like atoms on the (n, l, m) lattice
============================================================================

Two Hamiltonians are available, selected by the ``hamiltonian`` keyword:

``'exact'`` (default, since 2026-09-28)
    The physical one-electron Coulomb Hamiltonian T + V_ne expressed in the
    hydrogenic basis |n, l, m> (k = Z): it is exactly diagonal,

        H = diag(-Z^2 / (2 n^2)),

    so the spectrum is the Bohr/Rydberg ladder with degeneracy n^2 at every
    cutoff max_n.  This mirrors the He fix (``casimir_ci`` ``h1_method='exact'``).

``'graph'`` (legacy, opt-in)
    H = kappa * Z^2 * (D - A) with kappa = -1/16, the graph Laplacian of the
    ``GeometricLattice`` edge set.  Kept for reproducibility of graph-level
    results.  Its excited states are NOT physical: the lattice has no
    l-changing edges, the Laplacian splits into max_n disconnected l-blocks,
    and the low end of its spectrum is a dense cluster near -Z^2/2 (at
    max_n = 20: -0.4936, -0.4936, -0.4932, ... for Z = 1), not a Rydberg
    ladder.  The ground state approaches -Z^2/2 only by construction,
    E_0 = kappa * Z^2 * lambda_max with lambda_max -> 8 = 2 d_max
    (an Observation, not a derivation; Paper 7 / Paper 0).

Author: GeoVac Development Team
Date: February 2026 (exact default added 2026-09-28)
"""

import numpy as np
from scipy.sparse import diags, csr_matrix
from scipy.sparse.linalg import eigsh
from typing import Tuple, Optional

try:
    from .lattice import GeometricLattice
except ImportError:
    from lattice import GeometricLattice


KAPPA_GRAPH: float = -1 / 16
_HAMILTONIANS = ('exact', 'graph')


class AtomicSolver:
    """
    Single-electron hydrogen-like atomic solver on the (n, l, m) lattice.

    Parameters
    ----------
    max_n : int
        Maximum principal quantum number.
    Z : int, optional
        Nuclear charge (default 1).
    kinetic_scale : float, optional
        Base graph scale kappa (default -1/16).  Only used by
        ``hamiltonian='graph'``; the stored ``self.kinetic_scale`` is
        kappa * Z^2 in either mode.  Passing a non-default value together
        with ``hamiltonian='exact'`` raises ``ValueError`` (it would be
        silently ignored otherwise).
    hamiltonian : {'exact', 'graph'}, optional
        ``'exact'`` (default): H = diag(-Z^2/(2 n^2)), the physical one-electron
        Coulomb Hamiltonian in the hydrogenic basis.
        ``'graph'``: the legacy H = kappa * Z^2 * (D - A).  Its excited states
        are NOT physical, and E_0 -> -Z^2/2 only as kappa * lambda_max by
        construction (lambda_max -> 8).  See the module docstring.

    The lattice (``self.lattice``: states, adjacency) is the same object in
    both modes; only ``self.H`` differs.
    """

    def __init__(self, max_n: int, Z: int = 1, kinetic_scale: float = KAPPA_GRAPH,
                 hamiltonian: str = 'exact'):
        if hamiltonian not in _HAMILTONIANS:
            raise ValueError(
                f"hamiltonian must be one of {_HAMILTONIANS}, got {hamiltonian!r}")
        if hamiltonian == 'exact' and kinetic_scale != KAPPA_GRAPH:
            raise ValueError(
                "kinetic_scale only applies to hamiltonian='graph'; the exact "
                "Hamiltonian diag(-Z^2/(2n^2)) has no scale parameter")
        self.max_n = max_n
        self.Z = Z
        self.hamiltonian = hamiltonian
        # Graph scale: E ~ Z^2, so kappa is multiplied by Z^2 (graph mode only)
        self.kinetic_scale = kinetic_scale * (Z ** 2)

        # Build lattice
        self.lattice = GeometricLattice(max_n=max_n)
        self.n_states = self.lattice.num_states

        self.H = self._build_hamiltonian()

    def _build_hamiltonian(self) -> csr_matrix:
        """Build H for the selected ``self.hamiltonian`` mode."""
        if self.hamiltonian == 'exact':
            return self._build_exact_hamiltonian()
        return self._build_graph_hamiltonian()

    def _build_exact_hamiltonian(self) -> csr_matrix:
        """
        Physical one-electron Hamiltonian in the hydrogenic basis (k = Z).

        T + V_ne is exactly diagonal on |n, l, m>: H_ii = -Z^2 / (2 n_i^2).
        """
        n_vals = np.array([s[0] for s in self.lattice.states], dtype=float)
        diag = -(float(self.Z) ** 2) / (2.0 * n_vals ** 2)
        return diags(diag, 0, shape=(self.n_states, self.n_states), format='csr')

    def _build_graph_hamiltonian(self) -> csr_matrix:
        """
        Legacy graph Hamiltonian H = kinetic_scale * (D - A).

        NOT the physical spectrum beyond the ground state (see class docstring).
        """
        adjacency = self.lattice.adjacency

        # Degree matrix
        degree = np.array(adjacency.sum(axis=1)).flatten()
        D = diags(degree, 0, shape=(self.n_states, self.n_states), format='csr')

        # Graph Laplacian
        laplacian = D - adjacency

        # Hamiltonian
        H = self.kinetic_scale * laplacian

        return H

    def compute_ground_state(self, n_states: int = 1) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the lowest eigenvalues and eigenvectors.

        Parameters
        ----------
        n_states : int, optional
            Number of states to compute (default: 1)

        Returns
        -------
        energies : np.ndarray, shape (n_states,)
            Eigenvalues (energies) in ascending order
        wavefunctions : np.ndarray, shape (n_basis, n_states)
            Eigenvectors (wavefunctions), normalized.  In ``'exact'`` mode
            these are the basis vectors |n, l, m> (stable order within a
            degenerate shell).
        """
        if self.hamiltonian == 'exact':
            k = min(n_states, self.n_states)
            diag = self.H.diagonal()
            order = np.argsort(diag, kind='stable')[:k]
            energies = diag[order]
            wavefunctions = np.zeros((self.n_states, k))
            wavefunctions[order, np.arange(k)] = 1.0
            return energies, wavefunctions

        # Graph mode: sparse eigensolver (unchanged legacy behaviour)
        k = min(n_states, self.n_states - 1)
        energies, wavefunctions = eigsh(self.H, k=k, which='SA')

        return energies, wavefunctions

    def compute_expectation_value(self, operator: csr_matrix, state_index: int = 0) -> float:
        """
        Compute expectation value of an operator.

        Parameters
        ----------
        operator : csr_matrix
            Operator matrix
        state_index : int, optional
            Index of state to use (default: 0 = ground state)

        Returns
        -------
        expectation : float
            <ψ|O|ψ>
        """
        _, wavefunctions = self.compute_ground_state(n_states=state_index+1)
        psi = wavefunctions[:, state_index]

        expectation = psi.conj() @ (operator @ psi)

        return np.real(expectation)

    def __repr__(self) -> str:
        return (f"AtomicSolver(max_n={self.max_n}, Z={self.Z}, "
                f"n_states={self.n_states}, "
                f"hamiltonian={self.hamiltonian!r}, "
                f"kinetic_scale={self.kinetic_scale:.6f})")


def solve_hydrogen(max_n: int = 30, kinetic_scale: float = KAPPA_GRAPH,
                   hamiltonian: str = 'exact') -> Tuple[float, np.ndarray]:
    """
    Convenience function: hydrogen ground state.

    Parameters
    ----------
    max_n : int, optional
        Basis size (default 30)
    kinetic_scale : float, optional
        Graph scale kappa (default -1/16); ``hamiltonian='graph'`` only.
    hamiltonian : {'exact', 'graph'}, optional
        See ``AtomicSolver``.  Default ``'exact'``.

    Returns
    -------
    energy : float
        Ground state energy in Hartree
    wavefunction : np.ndarray
        Ground state wavefunction (normalized)

    Examples
    --------
    >>> E, psi = solve_hydrogen(max_n=30)
    >>> print(f"H ground state: {E:.6f} Ha")
    H ground state: -0.500000 Ha
    >>> E, _ = solve_hydrogen(max_n=30, hamiltonian='graph')
    >>> print(f"{E:.6f}")   # kappa * lambda_max, by construction
    -0.497131
    """
    solver = AtomicSolver(max_n=max_n, Z=1, kinetic_scale=kinetic_scale,
                          hamiltonian=hamiltonian)
    energies, wavefunctions = solver.compute_ground_state(n_states=1)

    return energies[0], wavefunctions[:, 0]


def solve_atom(Z: int, max_n: int = 30, kinetic_scale: float = KAPPA_GRAPH,
               hamiltonian: str = 'exact') -> Tuple[float, np.ndarray]:
    """
    Ground state of a single-electron hydrogen-like ion.

    With the default ``hamiltonian='exact'`` the energy is exactly -Z^2/2.
    With ``hamiltonian='graph'`` it is kappa * Z^2 * lambda_max (legacy).

    Parameters
    ----------
    Z : int
        Nuclear charge (Z=1 for H, Z=2 for He+, etc.)
    max_n : int, optional
        Basis size
    kinetic_scale : float, optional
        Graph scale kappa (default -1/16); ``hamiltonian='graph'`` only.
    hamiltonian : {'exact', 'graph'}, optional
        See ``AtomicSolver``.  Default ``'exact'``.

    Returns
    -------
    energy : float
        Ground state energy
    wavefunction : np.ndarray
        Ground state wavefunction

    Examples
    --------
    >>> E, psi = solve_atom(Z=2, max_n=30)   # He+: -2.0 Ha
    >>> E, psi = solve_atom(Z=3, max_n=30)   # Li2+: -4.5 Ha
    """
    solver = AtomicSolver(max_n=max_n, Z=Z, kinetic_scale=kinetic_scale,
                          hamiltonian=hamiltonian)
    energies, wavefunctions = solver.compute_ground_state(n_states=1)

    return energies[0], wavefunctions[:, 0]


if __name__ == "__main__":
    import sys
    import io

    # Set UTF-8 encoding for Windows
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

    print("=" * 70)
    print("ATOMIC SOLVER: exact (default) vs legacy graph Hamiltonian, Z = 1")
    print("=" * 70)
    for max_n in [10, 20, 30]:
        ex = AtomicSolver(max_n=max_n, Z=1)
        gr = AtomicSolver(max_n=max_n, Z=1, hamiltonian='graph')
        e_ex, _ = ex.compute_ground_state(n_states=5)
        e_gr, _ = gr.compute_ground_state(n_states=5)
        print(f"max_n={max_n:>3d}  exact: "
              + ", ".join(f"{e:.4f}" for e in np.sort(e_ex)))
        print(f"           graph: "
              + ", ".join(f"{e:.4f}" for e in np.sort(e_gr))
              + "   (NOT physical beyond E_0)")

