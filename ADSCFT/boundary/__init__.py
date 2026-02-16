"""
Boundary Theory (CFT) - Graph Representation

This module implements the boundary (Conformal Field Theory) side of the
AdS/CFT correspondence for the GeoVac framework.

The boundary theory is what we currently have in `geovac/`:
- Quantum states as graph nodes
- Hamiltonian as graph Laplacian
- Spectral methods for eigenvalues

This module provides an interface layer to translate `geovac.GeometricLattice`
and `geovac.AtomicSolver` into a form suitable for AdS/CFT translation.

Author: GeoVac Development Team
Date: February 14, 2026
Status: Experimental/Theoretical
"""

from .graph_states import GraphBoundary

__all__ = ['GraphBoundary']
