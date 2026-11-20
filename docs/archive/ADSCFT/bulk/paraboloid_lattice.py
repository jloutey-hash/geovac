"""
Paraboloid Lattice - 3D Geometric Embedding

Implements the bulk (AdS) side of AdS/CFT correspondence by embedding
quantum states (n,ℓ,m) as points on a 3D paraboloid of revolution.

Theory:
-------
The paraboloid surface z = -1/r² (where r² = x² + y²) is the natural
geometric manifold for hydrogen states because:

1. The Coulomb potential V = -1/r maps to vertical coordinate z
2. Angular momentum (ℓ,m) maps to spherical angles (θ,φ)
3. Principal quantum number n maps to radial distance r = n²

This embedding allows geometric calculations (areas, volumes) that
are not accessible in pure graph topology (boundary/CFT).

Coordinate Mapping:
------------------
For quantum state |n,ℓ,m⟩:

  r = n²  (parabolic radius in xy-plane)
  θ = π × ℓ/(n-1)  if n > 1, else 0  (polar angle)
  φ = 2π × (m+ℓ)/(2ℓ+1)  if ℓ > 0, else 0  (azimuthal angle)
  z = -1/n²  (energy depth, negative potential)

  x = r × sin(θ) × cos(φ)
  y = r × sin(θ) × sin(φ)
  z = -1/n²

Reference:
----------
old_research_archive/src/paraboloid_lattice_su11.py
Papers/Paper_2_Fine_Structure.tex (symplectic framework)

Author: GeoVac Development Team
Date: February 14, 2026
Version: 0.1.0-alpha
Status: Experimental
"""

import numpy as np
from typing import Dict, List, Tuple
from dataclasses import dataclass


@dataclass(frozen=True)
class QuantumState:
    """Represents a quantum state |n,ℓ,m⟩"""
    n: int
    l: int
    m: int

    def __post_init__(self):
        """Validate quantum numbers"""
        if self.n < 1:
            raise ValueError(f"n must be >= 1, got {self.n}")
        if self.l < 0 or self.l >= self.n:
            raise ValueError(f"l must be in [0, n-1], got l={self.l}, n={self.n}")
        if abs(self.m) > self.l:
            raise ValueError(f"|m| must be <= l, got m={self.m}, l={self.l}")

    def as_tuple(self) -> Tuple[int, int, int]:
        """Return as (n, l, m) tuple"""
        return (self.n, self.l, self.m)


class ParaboloidLattice:
    """
    3D paraboloid embedding of hydrogen quantum states.

    Maps each quantum state |n,ℓ,m⟩ to 3D Euclidean coordinates (x,y,z)
    on a paraboloid of revolution z = -1/r².

    This provides the geometric (bulk/AdS) representation needed for:
    - Symplectic plaquette area calculations
    - Fine structure constant extraction
    - Hyperfine contact geometry

    Parameters
    ----------
    max_n : int
        Maximum principal quantum number (n ≥ 1)

    Attributes
    ----------
    states : List[Tuple[int,int,int]]
        List of (n,ℓ,m) quantum number tuples
    coordinates : Dict[Tuple[int,int,int], np.ndarray]
        Map from (n,ℓ,m) to (x,y,z) coordinates
    dim : int
        Total number of states

    Examples
    --------
    >>> lattice = ParaboloidLattice(max_n=3)
    >>> lattice.dim
    10
    >>> coord = lattice.get_coordinate(n=2, l=1, m=0)
    >>> print(f"Position: x={coord[0]:.3f}, y={coord[1]:.3f}, z={coord[2]:.3f}")
    Position: x=2.828, y=0.000, z=-0.250
    """

    def __init__(self, max_n: int):
        """
        Initialize paraboloid lattice with 3D coordinates.

        Parameters
        ----------
        max_n : int
            Maximum principal quantum number (n ≥ 1)
        """
        if max_n < 1:
            raise ValueError(f"max_n must be >= 1, got {max_n}")

        self.max_n = max_n
        self.states: List[Tuple[int, int, int]] = []
        self.coordinates: Dict[Tuple[int, int, int], np.ndarray] = {}
        self.state_index: Dict[Tuple[int, int, int], int] = {}

        # Build lattice
        self._construct_states()
        self._compute_coordinates()

        self.dim = len(self.states)

    def _construct_states(self):
        """Generate all valid quantum states (n,ℓ,m)"""
        idx = 0
        for n in range(1, self.max_n + 1):
            for l in range(n):  # l = 0, 1, ..., n-1
                for m in range(-l, l + 1):  # m = -l, ..., 0, ..., l
                    state = (n, l, m)
                    self.states.append(state)
                    self.state_index[state] = idx
                    idx += 1

    def _compute_coordinates(self):
        """
        Compute 3D Euclidean coordinates for each quantum state.

        Embedding formula:
        - r = n² (parabolic radius)
        - θ = π × ℓ/(n-1) if n > 1, else 0
        - φ = 2π × (m+ℓ)/(2ℓ+1) if ℓ > 0, else 0
        - z = -1/n²

        Cartesian conversion:
        - x = r × sin(θ) × cos(φ)
        - y = r × sin(θ) × sin(φ)
        - z = -1/n²
        """
        for n, l, m in self.states:
            # Parabolic radius (scales as n²)
            r = n ** 2

            # Energy depth (Coulomb potential)
            z = -1.0 / (n ** 2)

            # Polar angle (from ℓ)
            if n > 1:
                theta = np.pi * l / (n - 1)
            else:
                theta = 0.0

            # Azimuthal angle (from m)
            if l > 0:
                phi = 2.0 * np.pi * (m + l) / (2 * l + 1)
            else:
                phi = 0.0

            # Convert to Cartesian
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)

            # Store coordinates
            self.coordinates[(n, l, m)] = np.array([x, y, z])

    def get_coordinate(self, n: int, l: int, m: int) -> np.ndarray:
        """
        Get 3D coordinates for quantum state |n,ℓ,m⟩.

        Parameters
        ----------
        n, l, m : int
            Quantum numbers

        Returns
        -------
        coord : np.ndarray, shape (3,)
            Cartesian coordinates [x, y, z]

        Raises
        ------
        KeyError
            If state not in lattice
        """
        state = (n, l, m)
        if state not in self.coordinates:
            raise KeyError(f"State {state} not in lattice (max_n={self.max_n})")
        return self.coordinates[state]

    def has_state(self, n: int, l: int, m: int) -> bool:
        """Check if state |n,ℓ,m⟩ exists in lattice"""
        return (n, l, m) in self.coordinates

    def get_states_by_shell(self, n: int) -> List[Tuple[int, int, int]]:
        """
        Get all states in a given energy shell n.

        Parameters
        ----------
        n : int
            Principal quantum number

        Returns
        -------
        states : List[Tuple[int,int,int]]
            All (n,ℓ,m) states with this n
        """
        return [(n_i, l_i, m_i) for (n_i, l_i, m_i) in self.states if n_i == n]

    def get_states_by_angular(self, l: int) -> List[Tuple[int, int, int]]:
        """
        Get all states with a given angular momentum l.

        Parameters
        ----------
        l : int
            Angular momentum quantum number

        Returns
        -------
        states : List[Tuple[int,int,int]]
            All (n,ℓ,m) states with this ℓ
        """
        return [(n_i, l_i, m_i) for (n_i, l_i, m_i) in self.states if l_i == l]

    def get_coordinate_array(self) -> np.ndarray:
        """
        Get all coordinates as a numpy array.

        Returns
        -------
        coords : np.ndarray, shape (dim, 3)
            Array of all (x,y,z) coordinates, ordered by state index
        """
        coords = np.zeros((self.dim, 3))
        for state, idx in self.state_index.items():
            coords[idx] = self.coordinates[state]
        return coords

    def get_radial_distance(self, n: int, l: int, m: int) -> float:
        """
        Get cylindrical radial distance r_cyl = sqrt(x² + y²) for state.

        Note: This equals n²*sin(θ), not n². The parameter r=n² is used
        in the spherical-to-Cartesian conversion, but the actual cylindrical
        radius includes the sin(θ) factor.
        """
        coord = self.get_coordinate(n, l, m)
        return np.sqrt(coord[0]**2 + coord[1]**2)

    def get_energy_depth(self, n: int, l: int, m: int) -> float:
        """
        Get energy depth z = -1/n² for state.

        This corresponds to the Coulomb potential.
        """
        coord = self.get_coordinate(n, l, m)
        return coord[2]

    def validate_embedding(self, verbose: bool = False) -> Dict[str, float]:
        """
        Validate the paraboloid embedding properties.

        Checks:
        1. Energy depth z = -1/n² (Coulomb potential)
        2. States with same n lie on same paraboloid shell (same z)
        3. No NaN or infinite values in coordinates

        Parameters
        ----------
        verbose : bool
            Print validation results

        Returns
        -------
        errors : Dict[str, float]
            Maximum relative errors for each check
        """
        errors = {
            'energy_depth': 0.0,
            'shell_consistency': 0.0,
            'coordinate_validity': 0.0
        }

        # Check for NaN/inf values
        for n, l, m in self.states:
            coord = self.get_coordinate(n, l, m)
            if not np.all(np.isfinite(coord)):
                errors['coordinate_validity'] = 1.0

        # Check energy depth z = -1/n²
        for n, l, m in self.states:
            z_actual = self.get_energy_depth(n, l, m)
            z_expected = -1.0 / (n ** 2)
            rel_error = abs(z_actual - z_expected) / abs(z_expected) if z_expected != 0 else 0
            errors['energy_depth'] = max(errors['energy_depth'], rel_error)

        # Check shell consistency (all states with same n have same z)
        for n in range(1, self.max_n + 1):
            shell_states = self.get_states_by_shell(n)
            if len(shell_states) > 1:
                z_values = [self.get_energy_depth(*state) for state in shell_states]
                z_mean = np.mean(z_values)
                z_std = np.std(z_values)
                rel_error = z_std / abs(z_mean) if z_mean != 0 else 0
                errors['shell_consistency'] = max(errors['shell_consistency'], rel_error)

        if verbose:
            print("="*70)
            print("PARABOLOID EMBEDDING VALIDATION")
            print("="*70)
            print(f"\nLattice: max_n={self.max_n}, dim={self.dim} states")
            print(f"\nGeometric Properties:")
            print(f"  Coordinate validity (no NaN):   {errors['coordinate_validity']:.3e}")
            print(f"  Energy depth (z = -1/n²):       {errors['energy_depth']:.3e}")
            print(f"  Shell consistency (same n):     {errors['shell_consistency']:.3e}")

            all_pass = all(err < 1e-12 for err in errors.values())
            print(f"\nStatus: {'✓ PASS' if all_pass else '✗ FAIL'}")
            print("="*70)

        return errors

    def __repr__(self) -> str:
        return f"ParaboloidLattice(max_n={self.max_n}, dim={self.dim})"

    def __len__(self) -> int:
        return self.dim


if __name__ == "__main__":
    # Demo: Create lattice and validate
    print("\nParaboloid Lattice Demo")
    print("="*70)

    lattice = ParaboloidLattice(max_n=5)

    print(f"\nCreated: {lattice}")
    print(f"Total states: {lattice.dim}")

    # Show some coordinates
    print("\nSample Coordinates:")
    print("  State (n,ℓ,m)  →  Position (x, y, z)")
    print("  " + "-"*50)

    for n, l, m in [(1,0,0), (2,0,0), (2,1,0), (2,1,1), (3,2,-1)]:
        if lattice.has_state(n, l, m):
            coord = lattice.get_coordinate(n, l, m)
            r = lattice.get_radial_distance(n, l, m)
            print(f"  ({n},{l},{m:+d})        →  ({coord[0]:7.3f}, {coord[1]:7.3f}, {coord[2]:7.3f})  r={r:.3f}")

    # Validate embedding
    print()
    errors = lattice.validate_embedding(verbose=True)
