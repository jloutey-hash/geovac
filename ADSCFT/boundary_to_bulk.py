"""
Boundary-to-Bulk Translation - AdS/CFT Correspondence

Implements the dictionary between:
- Boundary Theory (CFT/Graph): Quantum states as abstract nodes
- Bulk Theory (AdS/Geometric): Quantum states embedded in 3D space

This is the core of the AdS/CFT correspondence for GeoVac.

Translation Map:
---------------
Boundary (CFT)         ↔  Bulk (AdS)
--------------------------------------------------------------------------------
Graph nodes (n,ℓ,m)    ↔  3D coordinates (x,y,z) on paraboloid
Adjacency matrix       ↔  Geometric proximity
Laplacian eigenvalues  ↔  Energies from geometric Hamiltonian
Graph connectivity     ↔  Symplectic structure (plaquette areas)
Spectral dimension     ↔  Geometric dimension
Holographic entropy    ↔  Geometric area

Usage:
------
The translator allows moving from the well-validated boundary theory
(graph methods in `geovac/`) to the bulk theory (geometric methods)
needed for fine structure and hyperfine calculations.

Examples:
---------
>>> from ADSCFT.boundary_to_bulk import BoundaryBulkTranslator
>>> from geovac import AtomicSolver
>>>
>>> # Create boundary solver
>>> solver = AtomicSolver(Z=1, max_n=10)
>>>
>>> # Translate to bulk
>>> translator = BoundaryBulkTranslator(max_n=10)
>>> bulk_lattice = translator.embed_to_bulk()
>>>
>>> # Compute geometric quantities not available on boundary
>>> from ADSCFT.bulk.symplectic import compute_shell_capacity
>>> S_1 = compute_shell_capacity(bulk_lattice, n=1)
>>> print(f"Symplectic capacity S_1 = {S_1:.6f}")
Symplectic capacity S_1 = 0.433013

Author: GeoVac Development Team
Date: February 14, 2026
Version: 0.1.0-alpha
Status: Experimental
"""

import numpy as np
from typing import Dict, List, Tuple, Optional

from ADSCFT.boundary.graph_states import GraphBoundary
from ADSCFT.bulk.paraboloid_lattice import ParaboloidLattice
from ADSCFT.bulk.symplectic import (
    compute_shell_capacity,
    compute_capacity_series,
    compute_plaquette_area
)


class BoundaryBulkTranslator:
    """
    Translates between boundary (graph) and bulk (geometric) theories.

    This class implements the AdS/CFT correspondence by providing:
    1. Embedding: Graph states → 3D geometric coordinates
    2. Extraction: Geometric data → Graph properties
    3. Validation: Consistency checks between boundary and bulk

    Parameters
    ----------
    max_n : int
        Maximum principal quantum number
    Z : int, default=1
        Nuclear charge (for consistency checks)
    lepton_mass_factor : float, default=1.0
        Mass scaling factor

    Attributes
    ----------
    boundary : GraphBoundary
        Boundary (CFT/graph) representation
    bulk : ParaboloidLattice
        Bulk (AdS/geometric) representation

    Examples
    --------
    >>> translator = BoundaryBulkTranslator(max_n=10)
    >>> print(f"Boundary: {translator.boundary}")
    Boundary: GraphBoundary(Z=1, max_n=10, dim=385)
    >>> print(f"Bulk: {translator.bulk}")
    Bulk: ParaboloidLattice(max_n=10, dim=385)
    """

    def __init__(
        self,
        max_n: int,
        Z: int = 1,
        lepton_mass_factor: float = 1.0
    ):
        """
        Initialize boundary-to-bulk translator.

        Parameters
        ----------
        max_n : int
            Maximum principal quantum number
        Z : int
            Nuclear charge
        lepton_mass_factor : float
            Mass scaling factor (1.0 for electron, 206.768 for muon)
        """
        self.max_n = max_n
        self.Z = Z
        self.lepton_mass_factor = lepton_mass_factor

        # Create boundary (graph) representation
        self.boundary = GraphBoundary(
            Z=Z,
            max_n=max_n,
            lepton_mass_factor=lepton_mass_factor
        )

        # Create bulk (geometric) representation
        self.bulk = ParaboloidLattice(max_n=max_n)

        # Validate consistency
        self._validate_consistency()

    def _validate_consistency(self):
        """
        Validate that boundary and bulk have same state space.

        Checks:
        1. Same number of states
        2. Same quantum number lists
        """
        if self.boundary.dim != self.bulk.dim:
            raise ValueError(
                f"Dimension mismatch: boundary has {self.boundary.dim} states, "
                f"bulk has {self.bulk.dim} states"
            )

        # Check state lists match
        boundary_states = set(self.boundary.get_states())
        bulk_states = set(self.bulk.states)

        if boundary_states != bulk_states:
            missing_in_bulk = boundary_states - bulk_states
            missing_in_boundary = bulk_states - boundary_states
            raise ValueError(
                f"State mismatch:\n"
                f"  Missing in bulk: {missing_in_bulk}\n"
                f"  Missing in boundary: {missing_in_boundary}"
            )

    def embed_to_bulk(self) -> ParaboloidLattice:
        """
        Embed boundary (graph) states into bulk (geometric) lattice.

        Returns the ParaboloidLattice with 3D coordinates.

        Returns
        -------
        bulk_lattice : ParaboloidLattice
            Geometric lattice with 3D coordinates

        Examples
        --------
        >>> translator = BoundaryBulkTranslator(max_n=5)
        >>> bulk = translator.embed_to_bulk()
        >>> coord = bulk.get_coordinate(n=2, l=1, m=0)
        >>> print(f"State |2,1,0⟩ → (x,y,z) = ({coord[0]:.3f}, {coord[1]:.3f}, {coord[2]:.3f})")
        State |2,1,0⟩ → (x,y,z) = (2.828, 0.000, -0.250)
        """
        return self.bulk

    def extract_from_bulk(self) -> GraphBoundary:
        """
        Extract boundary (graph) properties from bulk (geometric) lattice.

        Currently returns the stored boundary representation.
        Future: Could reconstruct graph from geometric proximity.

        Returns
        -------
        boundary : GraphBoundary
            Graph representation
        """
        return self.boundary

    def get_coordinate_for_state(self, n: int, l: int, m: int) -> np.ndarray:
        """
        Get 3D bulk coordinates for boundary state.

        Parameters
        ----------
        n, l, m : int
            Quantum numbers

        Returns
        -------
        coord : np.ndarray, shape (3,)
            Cartesian coordinates [x, y, z]

        Examples
        --------
        >>> translator = BoundaryBulkTranslator(max_n=5)
        >>> coord = translator.get_coordinate_for_state(2, 1, 0)
        >>> print(f"(x,y,z) = ({coord[0]:.3f}, {coord[1]:.3f}, {coord[2]:.3f})")
        (x,y,z) = (2.828, 0.000, -0.250)
        """
        return self.bulk.get_coordinate(n, l, m)

    def compute_symplectic_capacity(self, n: int) -> float:
        """
        Compute symplectic capacity S_n for shell n.

        This is a purely geometric quantity not available on the boundary!

        Parameters
        ----------
        n : int
            Principal quantum number (energy shell)

        Returns
        -------
        S_n : float
            Total symplectic capacity (sum of plaquette areas)

        Examples
        --------
        >>> translator = BoundaryBulkTranslator(max_n=10)
        >>> S_1 = translator.compute_symplectic_capacity(n=1)
        >>> print(f"Ground state capacity: S_1 = {S_1:.6f}")
        Ground state capacity: S_1 = 0.433013
        """
        return compute_shell_capacity(self.bulk, n)

    def compute_capacity_series(self, n_max: Optional[int] = None) -> np.ndarray:
        """
        Compute symplectic capacities for shells n=1 to n_max.

        Parameters
        ----------
        n_max : int, optional
            Maximum shell (default: self.max_n - 1)

        Returns
        -------
        capacities : np.ndarray
            Array of S_n values

        Examples
        --------
        >>> translator = BoundaryBulkTranslator(max_n=10)
        >>> S_series = translator.compute_capacity_series(n_max=5)
        >>> for n, S_n in enumerate(S_series, start=1):
        ...     print(f"S_{n} = {S_n:.6f}")
        S_1 = 0.433013
        S_2 = 3.141593
        ...
        """
        if n_max is None:
            n_max = self.max_n - 1

        return compute_capacity_series(self.bulk, n_max=n_max)

    def validate_correspondence(self, verbose: bool = False) -> Dict[str, bool]:
        """
        Validate the AdS/CFT correspondence.

        Checks:
        1. State space consistency (boundary ↔ bulk)
        2. Dimension match
        3. Quantum number preservation

        Parameters
        ----------
        verbose : bool
            Print validation results

        Returns
        -------
        results : Dict[str, bool]
            Validation results for each check
        """
        results = {}

        # Check 1: Dimension match
        dim_match = (self.boundary.dim == self.bulk.dim)
        results['dimension_match'] = dim_match

        # Check 2: State space match
        boundary_states = set(self.boundary.get_states())
        bulk_states = set(self.bulk.states)
        state_match = (boundary_states == bulk_states)
        results['state_space_match'] = state_match

        # Check 3: All boundary states have bulk coordinates
        all_embedded = all(
            self.bulk.has_state(n, l, m)
            for n, l, m in self.boundary.get_states()
        )
        results['all_states_embedded'] = all_embedded

        if verbose:
            print("="*70)
            print("AdS/CFT CORRESPONDENCE VALIDATION")
            print("="*70)
            print(f"\nBoundary: {self.boundary}")
            print(f"Bulk:     {self.bulk}")
            print("\nChecks:")
            print(f"  Dimension match:       {dim_match} {'✓' if dim_match else '✗'}")
            print(f"  State space match:     {state_match} {'✓' if state_match else '✗'}")
            print(f"  All states embedded:   {all_embedded} {'✓' if all_embedded else '✗'}")

            all_pass = all(results.values())
            print(f"\nStatus: {'✓ PASS - Correspondence verified' if all_pass else '✗ FAIL - Inconsistency detected'}")
            print("="*70)

        return results

    def __repr__(self) -> str:
        return (
            f"BoundaryBulkTranslator("
            f"max_n={self.max_n}, "
            f"boundary_dim={self.boundary.dim}, "
            f"bulk_dim={self.bulk.dim})"
        )


if __name__ == "__main__":
    # Demo: Boundary-to-bulk translation
    print("\nBoundary-to-Bulk Translation Demo")
    print("="*70)

    # Create translator
    translator = BoundaryBulkTranslator(max_n=10, Z=1)

    print(f"\nTranslator: {translator}")

    # Validate correspondence
    print()
    results = translator.validate_correspondence(verbose=True)

    # Show some translations
    print("\nBoundary States → Bulk Coordinates:")
    print("  " + "-"*60)
    print("  State (n,ℓ,m)  →  Position (x, y, z)")
    print("  " + "-"*60)

    for n, l, m in [(1,0,0), (2,0,0), (2,1,0), (3,2,-1)]:
        if translator.boundary.has_state(n, l, m):
            coord = translator.get_coordinate_for_state(n, l, m)
            print(f"  ({n},{l},{m:+d})        →  ({coord[0]:8.3f}, {coord[1]:8.3f}, {coord[2]:8.3f})")

    # Compute symplectic capacities (bulk-only quantity)
    print("\nSymplectic Capacities (Bulk-Only Quantity):")
    print("  " + "-"*60)
    print("  Shell n  →  S_n (sum of plaquette areas)")
    print("  " + "-"*60)

    capacities = translator.compute_capacity_series(n_max=5)
    for n, S_n in enumerate(capacities, start=1):
        print(f"  n = {n}    →  S_{n} = {S_n:.6f}")

    print("\nNote: S_n values are NOT available in pure boundary (graph) theory!")
    print("      This demonstrates the power of the bulk (geometric) embedding.")

    print("\n" + "="*70)
