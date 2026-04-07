"""
Atomic Classifier for Composed Geometry (Track BG/CH/CR, v2.2.0)
=================================================================

Classifies atoms by nuclear charge Z, returning the composed-geometry
block specification needed by the general builder.  Structure types
follow Paper 16 (Chemical Periodicity as Representation Theory):

  Type A: single electron (H)            -- Level 1 system
  Type B: closed shell pair (He, Ne, Ar) -- Level 3 system
  Type C: single valence over closed core -- simplest composed (Li, Na, K)
  Type D: s² valence over closed core     -- paired valence (Be, Mg, Ca)
  Type E: open p-shell valence            -- exchange coupling important (B-F, Al-Cl, Ga-Br)

Coverage:
  First row (Z=1-10):   He-like cores, PK from Level 3 solver
  Second row (Z=11-18): [Ne] frozen cores (10e)
  Third row s-block (Z=19-20): [Ar] frozen cores (18e)
  Third row d-block (Z=21-30): OUT OF SCOPE (transition metals)
  Third row p-block (Z=31-36): [Ar]3d¹⁰ frozen cores (28e)

Author: GeoVac Development Team
Date: April 2026
"""

from dataclasses import dataclass, field
from typing import Dict, Optional


# ---------------------------------------------------------------------------
# PK anchor values from Paper 17, Table 1 (Li²⁺ core, inv2 method)
# ---------------------------------------------------------------------------
_PK_COMPUTED = {
    3: {'A': 6.93, 'B': 7.00},     # Li²⁺ core (Paper 17 Table 1, inv2)
    4: {'A': 13.01, 'B': 12.53},   # Be³⁺ core (Paper 17 Table 1, inv2)
    5: {'A': 21.40, 'B': 18.46},   # B⁴⁺ core (Track BI, Level 3 solver, inv2)
    6: {'A': 31.37, 'B': 25.54},   # C⁵⁺ core (Track BI, Level 3 solver, inv2)
    7: {'A': 43.09, 'B': 33.05},   # N⁶⁺ core (Track BI, Level 3 solver, inv2)
    9: {'A': 71.80, 'B': 48.61},   # F⁸⁺ core (Track BI, Level 3 solver, inv2)
}

# Li²⁺ anchor for Z² scaling
_PK_LI_A = 6.93
_PK_LI_B = 7.00
_PK_LI_Z = 3


# ---------------------------------------------------------------------------
# Atom database (Z=1-10, first row)
# ---------------------------------------------------------------------------
_ATOM_DATA = {
    1: {
        'structure_type': 'A',
        'n_core_electrons': 0,
        'n_valence_electrons': 1,
        'Z_eff_valence': 1.0,
        'core_config': 'none',
        'valence_config': '1s1',
        'period': 1,
        'group_type': 'hydrogen',
        'pk_source': 'none',
    },
    2: {
        'structure_type': 'B',
        'n_core_electrons': 0,
        'n_valence_electrons': 2,
        'Z_eff_valence': 2.0,
        'core_config': 'none',
        'valence_config': '1s2',
        'period': 1,
        'group_type': 'noble_gas',
        'pk_source': 'none',
    },
    3: {
        'structure_type': 'C',
        'n_core_electrons': 2,
        'n_valence_electrons': 1,
        'Z_eff_valence': 1.0,
        'core_config': '1s2',
        'valence_config': '2s1',
        'period': 2,
        'group_type': 'alkali_metal',
        'pk_source': 'computed',
    },
    4: {
        'structure_type': 'D',
        'n_core_electrons': 2,
        'n_valence_electrons': 2,
        'Z_eff_valence': 2.0,
        'core_config': '1s2',
        'valence_config': '2s2',
        'period': 2,
        'group_type': 'alkaline_earth',
        'pk_source': 'computed',
    },
    5: {
        'structure_type': 'E',
        'n_core_electrons': 2,
        'n_valence_electrons': 3,
        'Z_eff_valence': 3.0,
        'core_config': '1s2',
        'valence_config': '2s2 2p1',
        'period': 2,
        'group_type': 'p_block',
        'pk_source': 'computed',
    },
    6: {
        'structure_type': 'E',
        'n_core_electrons': 2,
        'n_valence_electrons': 4,
        'Z_eff_valence': 4.0,
        'core_config': '1s2',
        'valence_config': '2s2 2p2',
        'period': 2,
        'group_type': 'p_block',
        'pk_source': 'computed',
    },
    7: {
        'structure_type': 'E',
        'n_core_electrons': 2,
        'n_valence_electrons': 5,
        'Z_eff_valence': 5.0,
        'core_config': '1s2',
        'valence_config': '2s2 2p3',
        'period': 2,
        'group_type': 'p_block',
        'pk_source': 'computed',
    },
    8: {
        'structure_type': 'E',
        'n_core_electrons': 2,
        'n_valence_electrons': 6,
        'Z_eff_valence': 6.0,
        'core_config': '1s2',
        'valence_config': '2s2 2p4',
        'period': 2,
        'group_type': 'p_block',
        'pk_source': 'z2_scaled',
    },
    9: {
        'structure_type': 'E',
        'n_core_electrons': 2,
        'n_valence_electrons': 7,
        'Z_eff_valence': 7.0,
        'core_config': '1s2',
        'valence_config': '2s2 2p5',
        'period': 2,
        'group_type': 'p_block',
        'pk_source': 'computed',
    },
    10: {
        'structure_type': 'B',
        'n_core_electrons': 2,
        'n_valence_electrons': 8,
        'Z_eff_valence': 8.0,
        'core_config': '1s2',
        'valence_config': '2s2 2p6',
        'period': 2,
        'group_type': 'noble_gas',
        'pk_source': 'z2_scaled',
    },
    # ----- Second row (Z=11-18) -----
    11: {
        'structure_type': 'C',
        'n_core_electrons': 10,
        'n_valence_electrons': 1,
        'Z_eff_valence': 1.0,
        'core_config': '1s2 2s2 2p6',
        'valence_config': '3s1',
        'period': 3,
        'group_type': 'alkali_metal',
        'pk_source': 'frozen_core',
    },
    12: {
        'structure_type': 'D',
        'n_core_electrons': 10,
        'n_valence_electrons': 2,
        'Z_eff_valence': 2.0,
        'core_config': '1s2 2s2 2p6',
        'valence_config': '3s2',
        'period': 3,
        'group_type': 'alkaline_earth',
        'pk_source': 'frozen_core',
    },
    13: {
        'structure_type': 'E',
        'n_core_electrons': 10,
        'n_valence_electrons': 3,
        'Z_eff_valence': 3.0,
        'core_config': '1s2 2s2 2p6',
        'valence_config': '3s2 3p1',
        'period': 3,
        'group_type': 'p_block',
        'pk_source': 'frozen_core',
    },
    14: {
        'structure_type': 'E',
        'n_core_electrons': 10,
        'n_valence_electrons': 4,
        'Z_eff_valence': 4.0,
        'core_config': '1s2 2s2 2p6',
        'valence_config': '3s2 3p2',
        'period': 3,
        'group_type': 'p_block',
        'pk_source': 'frozen_core',
    },
    15: {
        'structure_type': 'E',
        'n_core_electrons': 10,
        'n_valence_electrons': 5,
        'Z_eff_valence': 5.0,
        'core_config': '1s2 2s2 2p6',
        'valence_config': '3s2 3p3',
        'period': 3,
        'group_type': 'p_block',
        'pk_source': 'frozen_core',
    },
    16: {
        'structure_type': 'E',
        'n_core_electrons': 10,
        'n_valence_electrons': 6,
        'Z_eff_valence': 6.0,
        'core_config': '1s2 2s2 2p6',
        'valence_config': '3s2 3p4',
        'period': 3,
        'group_type': 'p_block',
        'pk_source': 'frozen_core',
    },
    17: {
        'structure_type': 'E',
        'n_core_electrons': 10,
        'n_valence_electrons': 7,
        'Z_eff_valence': 7.0,
        'core_config': '1s2 2s2 2p6',
        'valence_config': '3s2 3p5',
        'period': 3,
        'group_type': 'p_block',
        'pk_source': 'frozen_core',
    },
    18: {
        'structure_type': 'B',
        'n_core_electrons': 10,
        'n_valence_electrons': 8,
        'Z_eff_valence': 8.0,
        'core_config': '1s2 2s2 2p6',
        'valence_config': '3s2 3p6',
        'period': 3,
        'group_type': 'noble_gas',
        'pk_source': 'frozen_core',
    },
    # ----- Third row s-block (Z=19-20), [Ar] frozen core -----
    19: {
        'structure_type': 'C',
        'n_core_electrons': 18,
        'n_valence_electrons': 1,
        'Z_eff_valence': 1.0,
        'core_config': '1s2 2s2 2p6 3s2 3p6',
        'valence_config': '4s1',
        'period': 4,
        'group_type': 'alkali_metal',
        'pk_source': 'frozen_core',
    },
    20: {
        'structure_type': 'D',
        'n_core_electrons': 18,
        'n_valence_electrons': 2,
        'Z_eff_valence': 2.0,
        'core_config': '1s2 2s2 2p6 3s2 3p6',
        'valence_config': '4s2',
        'period': 4,
        'group_type': 'alkaline_earth',
        'pk_source': 'frozen_core',
    },
    # ----- Z=21-30 (transition metals): OUT OF SCOPE -----
    # Handled by explicit check in classify_atom()
    # ----- Third row p-block (Z=31-36), [Ar]3d10 frozen core -----
    31: {
        'structure_type': 'E',
        'n_core_electrons': 28,
        'n_valence_electrons': 3,
        'Z_eff_valence': 3.0,
        'core_config': '1s2 2s2 2p6 3s2 3p6 3d10',
        'valence_config': '4s2 4p1',
        'period': 4,
        'group_type': 'p_block',
        'pk_source': 'frozen_core',
    },
    32: {
        'structure_type': 'E',
        'n_core_electrons': 28,
        'n_valence_electrons': 4,
        'Z_eff_valence': 4.0,
        'core_config': '1s2 2s2 2p6 3s2 3p6 3d10',
        'valence_config': '4s2 4p2',
        'period': 4,
        'group_type': 'p_block',
        'pk_source': 'frozen_core',
    },
    33: {
        'structure_type': 'E',
        'n_core_electrons': 28,
        'n_valence_electrons': 5,
        'Z_eff_valence': 5.0,
        'core_config': '1s2 2s2 2p6 3s2 3p6 3d10',
        'valence_config': '4s2 4p3',
        'period': 4,
        'group_type': 'p_block',
        'pk_source': 'frozen_core',
    },
    34: {
        'structure_type': 'E',
        'n_core_electrons': 28,
        'n_valence_electrons': 6,
        'Z_eff_valence': 6.0,
        'core_config': '1s2 2s2 2p6 3s2 3p6 3d10',
        'valence_config': '4s2 4p4',
        'period': 4,
        'group_type': 'p_block',
        'pk_source': 'frozen_core',
    },
    35: {
        'structure_type': 'E',
        'n_core_electrons': 28,
        'n_valence_electrons': 7,
        'Z_eff_valence': 7.0,
        'core_config': '1s2 2s2 2p6 3s2 3p6 3d10',
        'valence_config': '4s2 4p5',
        'period': 4,
        'group_type': 'p_block',
        'pk_source': 'frozen_core',
    },
    36: {
        'structure_type': 'B',
        'n_core_electrons': 28,
        'n_valence_electrons': 8,
        'Z_eff_valence': 8.0,
        'core_config': '1s2 2s2 2p6 3s2 3p6 3d10',
        'valence_config': '4s2 4p6',
        'period': 4,
        'group_type': 'noble_gas',
        'pk_source': 'frozen_core',
    },
}


@dataclass
class AtomClassification:
    """Classification of an atom for composed geometry construction.

    Attributes
    ----------
    Z : int
        Nuclear charge.
    N_electrons : int
        Total electron count (neutral atom).
    structure_type : str
        Paper 16 structure type: 'A', 'B', 'C', 'D', or 'E'.
    n_core_electrons : int
        Number of core electrons (0 for H, He).
    n_valence_electrons : int
        Number of valence electrons.
    Z_eff_valence : float
        Effective nuclear charge seen by valence electrons.
    nu : int
        Angular quantum number N - 2 (Paper 16).
    mu_free : float
        Pauli centrifugal cost 2 * nu^2.
    pk_params : Optional[Dict[str, float]]
        Phillips-Kleinman parameters {'A': float, 'B': float}, or None.
    pk_source : str
        Origin of PK params: 'computed', 'z2_scaled', or 'none'.
    core_config : str
        Core electron configuration, e.g. '1s2' or 'none'.
    valence_config : str
        Valence electron configuration, e.g. '2s2 2p3'.
    period : int
        Period in the periodic table.
    group_type : str
        Descriptive group: 'hydrogen', 'noble_gas', 'alkali_metal',
        'alkaline_earth', 'p_block'.
    supported : bool
        Whether this atom is supported for composed geometry.
    support_note : str
        Explanation if unsupported.
    """
    Z: int
    N_electrons: int
    structure_type: str
    n_core_electrons: int
    n_valence_electrons: int
    Z_eff_valence: float
    nu: int
    mu_free: float
    pk_params: Optional[Dict[str, float]]
    pk_source: str
    core_config: str
    valence_config: str
    period: int
    group_type: str
    supported: bool = True
    support_note: str = ''


def pk_params_z2_scaled(Z: int) -> Dict[str, float]:
    """Z^2 scaling of PK parameters from Li^2+ anchor.

    Parameters
    ----------
    Z : int
        Nuclear charge of the atom with the He-like 1s^2 core.

    Returns
    -------
    Dict[str, float]
        {'A': float, 'B': float} scaled PK parameters.

    Notes
    -----
    A scales with core-valence energy gap ~ Z^2.
    B scales with <1/r^2>_core ~ Z^2.
    Anchor: Li (Z=3), A=6.93, B=7.00 from Paper 17 Table 1.
    """
    scale = (Z / _PK_LI_Z) ** 2
    return {
        'A': _PK_LI_A * scale,
        'B': _PK_LI_B * scale,
    }


def pk_params_from_formulas(
    Z: int, E_core: float, r_eff: float
) -> Dict[str, float]:
    """Ab initio PK parameters from Paper 17 formulas.

    Parameters
    ----------
    Z : int
        Nuclear charge.
    E_core : float
        Core orbital energy (Ha).
    r_eff : float
        Effective core radius (bohr).

    Returns
    -------
    Dict[str, float]
        {'A': float, 'B': float} ab initio PK parameters.

    Notes
    -----
    From Paper 17 Sec IV:
        A = -2 * E_core  (energy gap)
        B = 1 / r_eff^2  (spatial confinement)
    """
    A = -2.0 * E_core
    B = 1.0 / (r_eff ** 2)
    return {'A': A, 'B': B}


def classify_atom(Z: int) -> AtomClassification:
    """Classify atom by nuclear charge Z for composed geometry.

    Parameters
    ----------
    Z : int
        Nuclear charge (must be >= 1).

    Returns
    -------
    AtomClassification
        Full classification including structure type, PK parameters,
        electron configuration, and support status.

    Raises
    ------
    ValueError
        If Z < 1.
    """
    if Z < 1:
        raise ValueError(f"Nuclear charge Z must be >= 1, got {Z}")

    # Transition metals (Z=21-30): explicitly out of scope
    if 21 <= Z <= 30:
        raise NotImplementedError(
            f"Transition metals (Z={Z}, Z=21-30) are out of scope. "
            "Multi-reference d-electron correlations require treatment "
            "beyond the composed geometry framework. "
            "See SCOPE_BOUNDARY.md for details."
        )

    # Unsupported atoms (Z > 36)
    if Z > 36:
        N = Z
        nu = max(N - 2, 0)
        mu_free = 2.0 * nu ** 2 if nu > 0 else 0.0
        return AtomClassification(
            Z=Z,
            N_electrons=N,
            structure_type='unknown',
            n_core_electrons=0,
            n_valence_electrons=N,
            Z_eff_valence=0.0,
            nu=nu,
            mu_free=mu_free,
            pk_params=None,
            pk_source='none',
            core_config='unknown',
            valence_config='unknown',
            period=0,
            group_type='unknown',
            supported=False,
            support_note=(
                f'Z={Z} is beyond supported range (Z=1-36). '
                'Extension requires frozen-core data for heavier '
                'noble gas cores.'
            ),
        )

    data = _ATOM_DATA[Z]
    N = Z  # neutral atom
    nu = max(N - 2, 0)
    mu_free = 2.0 * nu ** 2 if nu > 0 else 0.0

    # Determine PK parameters
    pk_source = data['pk_source']
    if pk_source == 'none':
        pk_params = None
    elif pk_source == 'computed':
        pk_params = dict(_PK_COMPUTED[Z])
    elif pk_source == 'frozen_core':
        # Frozen-core PK not yet computed; params pending
        pk_params = None
    else:
        # z2_scaled
        pk_params = pk_params_z2_scaled(Z)

    return AtomClassification(
        Z=Z,
        N_electrons=N,
        structure_type=data['structure_type'],
        n_core_electrons=data['n_core_electrons'],
        n_valence_electrons=data['n_valence_electrons'],
        Z_eff_valence=data['Z_eff_valence'],
        nu=nu,
        mu_free=mu_free,
        pk_params=pk_params,
        pk_source=pk_source,
        core_config=data['core_config'],
        valence_config=data['valence_config'],
        period=data['period'],
        group_type=data['group_type'],
    )
