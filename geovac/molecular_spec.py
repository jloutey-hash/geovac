"""
Molecular Specification for Composed Geometry (Track BH, v2.0.30)
=================================================================

Dataclass-based molecular specification that encodes the block structure
of composed-geometry Hamiltonians.  Each molecule is described as a list
of OrbitalBlock objects plus a nuclear repulsion constant.  The general
builder in ``composed_qubit.build_composed_hamiltonian`` consumes a
MolecularSpec and produces the second-quantized Hamiltonian.

Author: GeoVac Development Team
Date: April 2026
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class OrbitalBlock:
    """One orbital sub-block in a composed-geometry Hamiltonian.

    Each block maps to either 1 or 2 sub-blocks in the full orbital array:

    * ``core`` / ``lone_pair`` / ``bond_pair``: 1 sub-block (center orbitals only)
    * ``bond`` with ``has_h_partner=True``: 2 sub-blocks (center + partner)

    Attributes
    ----------
    label : str
        Human-readable label, e.g. ``"Li_core"``, ``"LiH_bond"``.
    block_type : str
        One of ``"core"``, ``"bond"``, ``"lone_pair"``, ``"bond_pair"``.
    Z_center : float
        Nuclear charge for the center-side orbital wavefunctions.
    n_electrons : int
        Number of electrons assigned to this block.
    max_n : int
        Maximum principal quantum number for center orbitals.
    has_h_partner : bool
        If True, this bond block also has partner (H-side) orbitals.
    Z_partner : float
        Nuclear charge for the partner-side orbitals (when ``has_h_partner``).
    max_n_partner : int
        Maximum n for partner orbitals.  0 means same as ``max_n``.
    pk_A : float
        PK barrier height (Ha*bohr^2).  0 means no PK on this block.
    pk_B : float
        PK barrier width exponent (bohr^-2).  0 means no PK.
    """

    label: str
    block_type: str
    Z_center: float
    n_electrons: int
    max_n: int
    has_h_partner: bool = False
    Z_partner: float = 1.0
    max_n_partner: int = 0
    pk_A: float = 0.0
    pk_B: float = 0.0

    def __post_init__(self) -> None:
        if self.max_n_partner == 0 and self.has_h_partner:
            self.max_n_partner = self.max_n


@dataclass
class MolecularSpec:
    """Complete specification of a composed-geometry molecule.

    Attributes
    ----------
    name : str
        Molecule name (e.g. ``"LiH"``).
    blocks : list of OrbitalBlock
        Ordered list of orbital blocks.
    nuclear_repulsion_constant : float
        Combined V_NN + V_cross + E_core constant (Ha).
    description : str
        Optional human-readable description.
    """

    name: str
    blocks: List[OrbitalBlock]
    nuclear_repulsion_constant: float
    description: str = ''
