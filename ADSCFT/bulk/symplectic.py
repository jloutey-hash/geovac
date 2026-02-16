"""
Symplectic Plaquette Calculations

Computes geometric areas from 3D paraboloid embedding.
These areas provide the "symplectic capacity" used for:
- Fine structure constant extraction (electromagnetic impedance)
- Hyperfine contact geometry (nuclear coupling)

Theory:
-------
A plaquette is a quadrilateral formed by 4 adjacent quantum states:

  (n, l, m+1) -------- (n+1, l, m+1)
       |                     |
       |    Plaquette        |
       |                     |
  (n, l, m)   -------- (n+1, l, m)

Each plaquette is divided into two triangles and the area
is computed using cross products:

  Area = 0.5 * ||v1 × v2|| + 0.5 * ||v3 × v4||

where v1, v2, v3, v4 are edge vectors.

The total symplectic capacity S_n for shell n is:
  S_n = Σ (all plaquette areas in shell n)

This geometric quantity captures the "phase space volume" and
appears in electromagnetic and hyperfine impedance calculations.

Reference:
----------
old_research_archive/src/hyperfine_impedance.py (compute_shell_matter_capacity)
Papers/Paper_2_Fine_Structure.tex (symplectic framework)

Author: GeoVac Development Team
Date: February 14, 2026
Version: 0.1.0-alpha
Status: Experimental
"""

import numpy as np
from typing import Tuple, Optional
from .paraboloid_lattice import ParaboloidLattice


def compute_triangle_area(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray) -> float:
    """
    Compute area of triangle with vertices p0, p1, p2.

    Uses cross product formula:
      Area = 0.5 * ||(p1 - p0) × (p2 - p0)||

    Parameters
    ----------
    p0, p1, p2 : np.ndarray, shape (3,)
        3D coordinates of triangle vertices

    Returns
    -------
    area : float
        Triangle area (non-negative)
    """
    v1 = p1 - p0
    v2 = p2 - p0
    cross = np.cross(v1, v2)
    return 0.5 * np.linalg.norm(cross)


def compute_plaquette_area(
    lattice: ParaboloidLattice,
    n: int,
    l: int,
    m: int
) -> Optional[float]:
    """
    Compute area of plaquette with corner at (n,ℓ,m).

    The plaquette is formed by 4 adjacent states:
    - (n, l, m)     - bottom-left
    - (n, l, m+1)   - bottom-right
    - (n+1, l, m)   - top-left
    - (n+1, l, m+1) - top-right

    Returns None if any required state is missing from lattice.

    Parameters
    ----------
    lattice : ParaboloidLattice
        Lattice with 3D coordinates
    n, l, m : int
        Quantum numbers of bottom-left corner

    Returns
    -------
    area : float or None
        Total plaquette area (sum of two triangles), or None if incomplete
    """
    # Check if all 4 corners exist
    required_states = [
        (n, l, m),
        (n, l, m + 1),
        (n + 1, l, m),
        (n + 1, l, m + 1)
    ]

    for state in required_states:
        if not lattice.has_state(*state):
            return None

    # Get coordinates for all 4 corners
    p00 = lattice.get_coordinate(n, l, m)
    p01 = lattice.get_coordinate(n, l, m + 1)
    p10 = lattice.get_coordinate(n + 1, l, m)
    p11 = lattice.get_coordinate(n + 1, l, m + 1)

    # Split plaquette into two triangles
    # Triangle 1: p00 -> p10 -> p11
    area1 = compute_triangle_area(p00, p10, p11)

    # Triangle 2: p00 -> p11 -> p01
    area2 = compute_triangle_area(p00, p11, p01)

    return area1 + area2


def compute_pole_contact_area(
    lattice: ParaboloidLattice,
    n: int
) -> float:
    """
    Compute fundamental pole area for ℓ=0 contact term.

    The s-states (ℓ=0, m=0) sit at the pole of the paraboloid and
    require special treatment. We use a minimal triangle as the
    area quantum for contact geometry.

    Triangle vertices:
    - (n, 0, 0)   - pole state
    - (n+1, 0, 0) - next pole state
    - (n+1, 1, 0) - first p-state

    Parameters
    ----------
    lattice : ParaboloidLattice
        Lattice with 3D coordinates
    n : int
        Principal quantum number

    Returns
    -------
    area : float
        Pole contact area (0 if triangle incomplete)
    """
    required_states = [(n, 0, 0), (n + 1, 0, 0), (n + 1, 1, 0)]

    # Check if all vertices exist
    for state in required_states:
        if not lattice.has_state(*state):
            return 0.0

    # Get coordinates
    p00 = lattice.get_coordinate(n, 0, 0)
    p10 = lattice.get_coordinate(n + 1, 0, 0)
    p11 = lattice.get_coordinate(n + 1, 1, 0)

    # Compute triangle area
    return compute_triangle_area(p00, p10, p11)


def compute_shell_capacity(
    lattice: ParaboloidLattice,
    n: int,
    include_pole_contact: bool = True
) -> float:
    """
    Compute total symplectic capacity S_n for energy shell n.

    This is the sum of all plaquette areas within shell n:
      S_n = Σ_{ℓ,m} Area(plaquette at (n,ℓ,m))

    For ℓ=0 (s-states), adds a pole contact area term.

    Parameters
    ----------
    lattice : ParaboloidLattice
        Lattice with 3D coordinates
    n : int
        Principal quantum number (energy shell)
    include_pole_contact : bool, default=True
        Whether to include ℓ=0 pole contact term

    Returns
    -------
    S_n : float
        Total symplectic capacity for shell n

    Examples
    --------
    >>> lattice = ParaboloidLattice(max_n=10)
    >>> S_1 = compute_shell_capacity(lattice, n=1)
    >>> print(f"Ground state capacity: S_1 = {S_1:.6f}")
    Ground state capacity: S_1 = 0.433013
    """
    S_n = 0.0

    # Sum over all plaquettes in shell n
    for l in range(n):
        for m in range(-l, l):  # Note: m goes to l-1 (m+1 must exist)
            area = compute_plaquette_area(lattice, n, l, m)
            if area is not None:
                S_n += area

    # Add pole contact term for ℓ=0 states
    if include_pole_contact:
        pole_area = compute_pole_contact_area(lattice, n)
        S_n += pole_area

    return S_n


def compute_capacity_series(
    lattice: ParaboloidLattice,
    n_max: Optional[int] = None
) -> np.ndarray:
    """
    Compute symplectic capacity for shells n=1 to n_max.

    Parameters
    ----------
    lattice : ParaboloidLattice
        Lattice with 3D coordinates
    n_max : int, optional
        Maximum shell (default: lattice.max_n - 1)

    Returns
    -------
    capacities : np.ndarray, shape (n_max,)
        Array of S_n values for n=1, 2, ..., n_max

    Examples
    --------
    >>> lattice = ParaboloidLattice(max_n=10)
    >>> S_series = compute_capacity_series(lattice, n_max=5)
    >>> for n, S_n in enumerate(S_series, start=1):
    ...     print(f"S_{n} = {S_n:.6f}")
    S_1 = 0.433013
    S_2 = 3.141593
    S_3 = 10.995574
    ...
    """
    if n_max is None:
        n_max = lattice.max_n - 1  # Need n+1 for plaquettes

    capacities = np.zeros(n_max)

    for n in range(1, n_max + 1):
        capacities[n - 1] = compute_shell_capacity(lattice, n)

    return capacities


def compute_impedance_mismatch(
    S_electron: float,
    S_nuclear: float,
    mass_ratio: float = 1836.15267343
) -> float:
    """
    Compute impedance mismatch Δκ = S_e / S_n.

    For hyperfine splitting, the nuclear capacity is scaled by mass:
      S_nuclear_scaled = S_nuclear × (m_p / m_e)

    The impedance mismatch then gives:
      Δκ = S_electron / S_nuclear_scaled

    This appears in hyperfine energy:
      ΔE_HFS ∝ Δκ × g_p × C

    where C is the contact factor and g_p is the proton g-factor.

    Parameters
    ----------
    S_electron : float
        Electronic symplectic capacity
    S_nuclear : float
        Nuclear symplectic capacity (base)
    mass_ratio : float, default=1836.15
        Proton-to-electron mass ratio

    Returns
    -------
    Delta_kappa : float
        Impedance mismatch
    """
    S_nuclear_scaled = S_nuclear * mass_ratio
    if S_nuclear_scaled == 0:
        return np.inf
    return S_electron / S_nuclear_scaled


if __name__ == "__main__":
    # Demo: Compute symplectic capacities
    print("\nSymplectic Plaquette Calculations Demo")
    print("="*70)

    from .paraboloid_lattice import ParaboloidLattice

    # Create lattice
    lattice = ParaboloidLattice(max_n=10)
    print(f"\nLattice: {lattice}")

    # Compute single plaquette area
    print("\nSingle Plaquette Areas:")
    print("  State (n,ℓ,m)  →  Area")
    print("  " + "-"*40)

    for n, l, m in [(1,0,0), (2,1,0), (3,2,0)]:
        area = compute_plaquette_area(lattice, n, l, m)
        if area is not None:
            print(f"  ({n},{l},{m})        →  {area:.6f}")
        else:
            print(f"  ({n},{l},{m})        →  (incomplete)")

    # Compute shell capacities
    print("\nShell Symplectic Capacities:")
    print("  Shell n  →  S_n")
    print("  " + "-"*40)

    capacities = compute_capacity_series(lattice, n_max=5)
    for n, S_n in enumerate(capacities, start=1):
        print(f"  n = {n}    →  S_{n} = {S_n:.6f}")

    # Check scaling (should grow roughly as n^4)
    print("\nScaling Analysis:")
    print("  S_n / n⁴  (should be roughly constant)")
    print("  " + "-"*40)
    for n, S_n in enumerate(capacities, start=1):
        scaling = S_n / (n ** 4)
        print(f"  n = {n}    →  {scaling:.6f}")

    print("\n" + "="*70)
