"""
Bulk Theory (AdS) - Geometric Embedding

This module implements the bulk (Anti-de Sitter space) side of the
AdS/CFT correspondence for the GeoVac framework.

Key components:
- 3D paraboloid embedding of quantum states
- Symplectic plaquette area calculations
- Geometric impedance computations

Author: GeoVac Development Team
Date: February 14, 2026
Status: Experimental/Theoretical
"""

from .paraboloid_lattice import ParaboloidLattice
from .symplectic import compute_plaquette_area, compute_shell_capacity

__all__ = [
    'ParaboloidLattice',
    'compute_plaquette_area',
    'compute_shell_capacity',
]
