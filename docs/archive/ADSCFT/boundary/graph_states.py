"""
Graph Boundary - Interface to GeoVac

Provides a clean interface to the boundary (CFT/graph) theory
implemented in `geovac/` for use in AdS/CFT translation.

This wraps `geovac.GeometricLattice` and `geovac.AtomicSolver`
to provide the graph-theoretic data needed for boundary-to-bulk
translation.

Author: GeoVac Development Team
Date: February 14, 2026
Version: 0.1.0-alpha
Status: Experimental
"""

import numpy as np
from typing import List, Tuple, Optional
import sys
import os

# Add parent directory to path to import geovac
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from geovac.lattice import GeometricLattice
from geovac.atomic_solver import AtomicSolver


class GraphBoundary:
    """
    Interface to boundary (CFT/graph) theory from geovac.

    Wraps `GeometricLattice` and `AtomicSolver` to provide:
    - List of quantum states (graph nodes)
    - Adjacency matrix (graph edges)
    - Laplacian matrix (graph operator)
    - Eigenvalues (energies)

    This data represents the "boundary theory" that will be
    translated to bulk (geometric) theory via AdS/CFT.

    Parameters
    ----------
    Z : int
        Nuclear charge (Z=1 for hydrogen)
    max_n : int
        Maximum principal quantum number
    lepton_mass_factor : float, default=1.0
        Mass scaling (1.0 for electron, 206.768 for muon)

    Attributes
    ----------
    lattice : GeometricLattice
        Graph structure
    solver : AtomicSolver
        Quantum solver (if needed)
    states : List[Tuple[int,int,int]]
        List of (n,ℓ,m) quantum states
    dim : int
        Number of states

    Examples
    --------
    >>> boundary = GraphBoundary(Z=1, max_n=5)
    >>> print(f"Boundary has {boundary.dim} states")
    Boundary has 35 states
    >>> states = boundary.get_states()
    >>> print(f"First 3 states: {states[:3]}")
    First 3 states: [(1,0,0), (2,0,0), (2,1,-1)]
    """

    def __init__(
        self,
        Z: int = 1,
        max_n: int = 5,
        lepton_mass_factor: float = 1.0
    ):
        """
        Initialize graph boundary interface.

        Parameters
        ----------
        Z : int
            Nuclear charge
        max_n : int
            Maximum principal quantum number
        lepton_mass_factor : float
            Mass scaling factor (1.0 for electron, 206.768 for muon)
        """
        self.Z = Z
        self.max_n = max_n
        self.lepton_mass_factor = lepton_mass_factor

        # Create lattice (boundary structure)
        self.lattice = GeometricLattice(max_n=max_n)

        # Create solver (optional, for energies)
        self.solver: Optional[AtomicSolver] = None

        # Extract states
        self.states = self.lattice.states  # Direct attribute access
        self.dim = len(self.states)

    def get_states(self) -> List[Tuple[int, int, int]]:
        """
        Get list of quantum states (graph nodes).

        Returns
        -------
        states : List[Tuple[int,int,int]]
            List of (n,ℓ,m) tuples
        """
        return self.states.copy()

    def get_adjacency_matrix(self):
        """
        Get graph adjacency matrix.

        Returns
        -------
        A : scipy.sparse matrix
            Adjacency matrix (0/1 edges between states)
        """
        return self.lattice.adjacency

    def get_laplacian_matrix(self):
        """
        Get graph Laplacian matrix L = D - A.

        Returns
        -------
        L : scipy.sparse matrix
            Graph Laplacian operator
        """
        import scipy.sparse as sp
        A = self.lattice.adjacency
        degrees = np.array(A.sum(axis=1)).flatten()
        D = sp.diags(degrees, format='csr')
        return D - A

    def get_eigenvalues(self, k: int = 10) -> np.ndarray:
        """
        Get eigenvalues of graph Laplacian.

        Note: This requires initializing the solver.

        Parameters
        ----------
        k : int
            Number of eigenvalues to compute

        Returns
        -------
        eigenvalues : np.ndarray
            Smallest k eigenvalues of Laplacian
        """
        if self.solver is None:
            # Initialize solver
            self.solver = AtomicSolver(
                Z=self.Z,
                max_n=self.max_n,
                lepton_mass_factor=self.lepton_mass_factor
            )

        # Solve for energies
        energies, _ = self.solver.solve(n_states=k)
        return energies

    def get_state_index(self, n: int, l: int, m: int) -> int:
        """
        Get index of state (n,ℓ,m) in state list.

        Parameters
        ----------
        n, l, m : int
            Quantum numbers

        Returns
        -------
        index : int
            Index in state list (0 to dim-1)

        Raises
        ------
        KeyError
            If state not in lattice
        """
        return self.lattice._state_index[(n, l, m)]

    def has_state(self, n: int, l: int, m: int) -> bool:
        """Check if state exists in boundary"""
        return (n, l, m) in self.lattice._state_index

    def __repr__(self) -> str:
        return f"GraphBoundary(Z={self.Z}, max_n={self.max_n}, dim={self.dim})"

    def __len__(self) -> int:
        return self.dim


if __name__ == "__main__":
    # Demo: Create boundary and inspect
    print("\nGraph Boundary Demo")
    print("="*70)

    boundary = GraphBoundary(Z=1, max_n=5)

    print(f"\nCreated: {boundary}")
    print(f"Total states: {boundary.dim}")

    # Show states
    print("\nFirst 10 States (n,ℓ,m):")
    print("  " + "-"*40)
    states = boundary.get_states()
    for i, (n, l, m) in enumerate(states[:10]):
        print(f"  {i:2d}. ({n},{l},{m:+d})")

    # Graph structure
    A = boundary.get_adjacency_matrix()
    L = boundary.get_laplacian_matrix()

    print(f"\nGraph Structure:")
    print(f"  Adjacency matrix: {A.shape}, {A.nnz} non-zero")
    print(f"  Laplacian matrix: {L.shape}, {L.nnz} non-zero")

    # Eigenvalues (energies)
    print(f"\nComputing eigenvalues...")
    energies = boundary.get_eigenvalues(k=5)

    print(f"\nFirst 5 Energies (graph eigenvalues):")
    print("  " + "-"*40)
    for i, E in enumerate(energies):
        print(f"  E_{i} = {E:.6f} Ry")

    print("\n" + "="*70)
