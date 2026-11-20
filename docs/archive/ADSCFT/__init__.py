"""
AdS/CFT Correspondence for GeoVac

This package implements the holographic duality (AdS/CFT correspondence)
between boundary theory (graph/CFT) and bulk theory (geometric/AdS).

**Status:** Experimental/Theoretical - May not be in production releases
**Precision:** Less precise than core graph methods
**Purpose:** Provides geometric embedding needed for fine structure and
            detailed contact geometry calculations

Modules:
--------
- boundary: CFT/graph representation (current geovac implementation)
- bulk: AdS/geometric representation (3D paraboloid embedding)
- boundary_to_bulk: Translation layer between theories

Quick Start:
-----------
>>> from ADSCFT import BoundaryBulkTranslator
>>>
>>> # Create translator
>>> translator = BoundaryBulkTranslator(max_n=10)
>>>
>>> # Get bulk lattice with 3D coordinates
>>> bulk = translator.embed_to_bulk()
>>>
>>> # Compute geometric quantities
>>> S_1 = translator.compute_symplectic_capacity(n=1)
>>> print(f"Symplectic capacity: S_1 = {S_1:.6f}")
Symplectic capacity: S_1 = 0.433013

Reference:
----------
See ADSCFT/README.md for detailed documentation.

Author: GeoVac Development Team
Date: February 14, 2026
Version: 0.1.0-alpha
"""

from ADSCFT.boundary_to_bulk import BoundaryBulkTranslator
from ADSCFT.boundary.graph_states import GraphBoundary
from ADSCFT.bulk.paraboloid_lattice import ParaboloidLattice
from ADSCFT.bulk.symplectic import (
    compute_plaquette_area,
    compute_shell_capacity,
    compute_capacity_series,
    compute_impedance_mismatch
)
from ADSCFT.fine_structure import FineStructureCalculator
from ADSCFT.proton_radius import ProtonRadiusCalculator, solve_proton_radius_puzzle
from ADSCFT.hyperfine_impedance import HyperfineImpedanceCalculator
from ADSCFT.muonic_hydrogen import (
    MuonicHydrogenSolver,
    MUON_ELECTRON_MASS_RATIO,
    MUON_REDUCED_MASS_RATIO,
    solve_muonic_hydrogen,
)
from ADSCFT.holographic_analysis import (
    compute_spectral_dimension,
    compute_holographic_entropy,
    extract_central_charge,
    compare_holographic_properties,
)

__version__ = '0.1.0-alpha'

__all__ = [
    'BoundaryBulkTranslator',
    'GraphBoundary',
    'ParaboloidLattice',
    'compute_plaquette_area',
    'compute_shell_capacity',
    'compute_capacity_series',
    'compute_impedance_mismatch',
    'FineStructureCalculator',
    'ProtonRadiusCalculator',
    'solve_proton_radius_puzzle',
    'HyperfineImpedanceCalculator',
    'MuonicHydrogenSolver',
    'MUON_ELECTRON_MASS_RATIO',
    'MUON_REDUCED_MASS_RATIO',
    'solve_muonic_hydrogen',
    'compute_spectral_dimension',
    'compute_holographic_entropy',
    'extract_central_charge',
    'compare_holographic_properties',
]
