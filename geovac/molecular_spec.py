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
    l_min: int = 0
    center_nucleus_idx: int = -1
    partner_nucleus_idx: int = -1

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
    core_method: str = 'pk'  # 'pk' or 'downfolded'
    nuclei: Optional[List[Dict]] = None  # [{'Z': float, 'position': (x,y,z), 'label': str}]


# ---------------------------------------------------------------------------
# Transition metal hydride spec factories (Track CZ/DA extension, v2.8.0)
# ---------------------------------------------------------------------------

# Element symbols for Z=21-30
_TM_ELEMENTS = {
    21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn',
    26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
}

# Experimental equilibrium bond lengths in bohr (CRC / NIST)
_TM_HYDRIDE_REQ = {
    21: 3.49, 22: 3.32, 23: 3.23, 24: 3.17, 25: 3.23,
    26: 3.12, 27: 3.02, 28: 2.88, 29: 2.83, 30: 3.08,
}


def transition_metal_hydride_spec(
    Z: int,
    R: Optional[float] = None,
    max_n_bond: int = 2,
    max_n_d: int = 3,
) -> MolecularSpec:
    """Create a MolecularSpec for a first-row transition metal hydride (XH).

    Block structure:
      - Bond block (sigma): screened TM center (Z_eff = Z - 18) + H partner.
        Encodes the 4s electrons that participate in bonding.
      - d-block (lone pair): d-only orbitals (l_min=2) on the TM center.
        Encodes the 3d electrons.

    The [Ar] core (18 electrons) is treated as a frozen core via
    ``neon_core.FrozenCore``.

    Parameters
    ----------
    Z : int
        Nuclear charge of the transition metal (21 <= Z <= 30).
    R : float, optional
        Internuclear distance in bohr. Defaults to experimental R_eq.
    max_n_bond : int
        Maximum n for bond-block orbitals (default 2).
    max_n_d : int
        Maximum n for d-block orbitals (default 3, minimum for l=2).

    Returns
    -------
    MolecularSpec
        Complete specification for the XH hydride.
    """
    if Z < 21 or Z > 30:
        raise ValueError(f"Z={Z} is not a first-row transition metal (21-30).")

    from geovac.neon_core import FrozenCore

    n_s = _TM_N_S[Z]
    n_d_val = _TM_N_D[Z]
    element = _TM_ELEMENTS[Z]

    if R is None:
        R = _TM_HYDRIDE_REQ[Z]

    if max_n_d < 3:
        raise ValueError(f"max_n_d must be >= 3 for d-orbitals, got {max_n_d}.")

    Z_eff_val = float(Z - 18)
    Z_H = 1.0

    # Frozen core
    fc = FrozenCore(Z)
    fc.solve()
    E_core = fc.energy

    # Nuclear repulsion
    V_NN = float(Z) * Z_H / R

    # Core-to-H attraction
    from geovac.composed_qubit import _v_cross_nuc_frozen_core
    V_cross = _v_cross_nuc_frozen_core(Z, Z_H, R)

    nuclear_repulsion = V_NN + V_cross + E_core

    blocks = [
        OrbitalBlock(
            label=f'{element}H_sigma',
            block_type='bond',
            Z_center=Z_eff_val,
            n_electrons=n_s,
            max_n=max_n_bond,
            has_h_partner=True,
            Z_partner=Z_H,
            max_n_partner=max_n_bond,
        ),
        OrbitalBlock(
            label=f'{element}_3d',
            block_type='lone_pair',
            Z_center=Z_eff_val,
            n_electrons=n_d_val,
            max_n=max_n_d,
            l_min=2,
        ),
    ]

    return MolecularSpec(
        name=f'{element}H',
        blocks=blocks,
        nuclear_repulsion_constant=nuclear_repulsion,
        description=(
            f'{element}H: [Ar] frozen core (18e), '
            f'{n_s} bond e- + {n_d_val} d e- = {n_s + n_d_val} encoded'
        ),
    )


# d-electron and s-electron counts for Z=21-30
_TM_N_D = {
    21: 1, 22: 2, 23: 3, 24: 5, 25: 5,
    26: 6, 27: 7, 28: 8, 29: 10, 30: 10,
}
_TM_N_S = {
    21: 2, 22: 2, 23: 2, 24: 1, 25: 2,
    26: 2, 27: 2, 28: 2, 29: 1, 30: 2,
}
# Legacy alias used in transition_metal_hydride_spec
_ATOM_DATA_N_S = _TM_N_S


# Convenience aliases
def sch_spec(**kw) -> MolecularSpec:
    """ScH (scandium hydride) spec."""
    return transition_metal_hydride_spec(21, **kw)

def tih_spec(**kw) -> MolecularSpec:
    """TiH (titanium hydride) spec."""
    return transition_metal_hydride_spec(22, **kw)

def vh_spec(**kw) -> MolecularSpec:
    """VH (vanadium hydride) spec."""
    return transition_metal_hydride_spec(23, **kw)

def crh_spec(**kw) -> MolecularSpec:
    """CrH (chromium hydride) spec."""
    return transition_metal_hydride_spec(24, **kw)

def mnh_spec(**kw) -> MolecularSpec:
    """MnH (manganese hydride) spec."""
    return transition_metal_hydride_spec(25, **kw)

def feh_spec(**kw) -> MolecularSpec:
    """FeH (iron hydride) spec."""
    return transition_metal_hydride_spec(26, **kw)

def coh_spec(**kw) -> MolecularSpec:
    """CoH (cobalt hydride) spec."""
    return transition_metal_hydride_spec(27, **kw)

def nih_spec(**kw) -> MolecularSpec:
    """NiH (nickel hydride) spec."""
    return transition_metal_hydride_spec(28, **kw)

def cuh_spec(**kw) -> MolecularSpec:
    """CuH (copper hydride) spec."""
    return transition_metal_hydride_spec(29, **kw)

def znh_spec(**kw) -> MolecularSpec:
    """ZnH (zinc hydride) spec."""
    return transition_metal_hydride_spec(30, **kw)


# ---------------------------------------------------------------------------
# General main-group hydride spec factory
# ---------------------------------------------------------------------------

# Element data: symbol, period, n_core_electrons, core_type
_ELEMENT_DATA = {
    # First row (explicit 1s² core in Hamiltonian)
    3:  ('Li', 2, 2,  'explicit'),
    4:  ('Be', 2, 2,  'explicit'),
    5:  ('B',  2, 2,  'explicit'),
    6:  ('C',  2, 2,  'explicit'),
    7:  ('N',  2, 2,  'explicit'),
    8:  ('O',  2, 2,  'explicit'),
    9:  ('F',  2, 2,  'explicit'),
    # Second row ([Ne] frozen core)
    11: ('Na', 3, 10, 'frozen'),
    12: ('Mg', 3, 10, 'frozen'),
    13: ('Al', 3, 10, 'frozen'),
    14: ('Si', 3, 10, 'frozen'),
    15: ('P',  3, 10, 'frozen'),
    16: ('S',  3, 10, 'frozen'),
    17: ('Cl', 3, 10, 'frozen'),
    # Third row s-block ([Ar] frozen core)
    19: ('K',  4, 18, 'frozen'),
    20: ('Ca', 4, 18, 'frozen'),
    # Third row p-block ([Ar]3d10 frozen core)
    31: ('Ga', 4, 28, 'frozen'),
    32: ('Ge', 4, 28, 'frozen'),
    33: ('As', 4, 28, 'frozen'),
    34: ('Se', 4, 28, 'frozen'),
    35: ('Br', 4, 28, 'frozen'),
    36: ('Kr', 4, 28, 'frozen'),
}

# Experimental R_eq in bohr for XHn hydrides (CRC / NIST)
_HYDRIDE_REQ = {
    'LiH':  3.015, 'BeH2': 2.502,
    'CH4':  2.050, 'NH3':  1.912, 'H2O':  1.809, 'HF':   1.733,
    'NaH':  3.566, 'MgH2': 3.261,
    'SiH4': 2.800, 'PH3':  2.680, 'H2S':  2.534, 'HCl':  2.409,
    'KH':   4.243, 'CaH2': 3.807,
    'GeH4': 2.870, 'AsH3': 2.820, 'H2Se': 2.760, 'HBr':  2.670,
}

# He-like core energies for first-row atoms (2-electron core, Z=Z_heavy)
# Uses hyperspherical solver values where available, else variational -(Z-5/16)²
_FIRST_ROW_CORE_ENERGY = {
    3: -7.2799,    # Li²⁺, hyperspherical (Paper 13)
    4: -13.65,     # Be²⁺, variational estimate
    5: -21.97,     # B³⁺
    6: -32.35,     # C⁴⁺
    7: -44.73,     # N⁵⁺
    8: -59.10,     # O⁶⁺
    9: -75.47,     # F⁷⁺
}

# Ab initio PK parameters from atomic classifier (Paper 17, Track BI)
_PK_PARAMS = {
    3:  (6.93,  7.00),    # Li, Paper 17 Table 1
    4:  (13.01, 12.53),   # Be
    5:  (21.4,  18.46),   # B
    6:  (31.37, 25.54),   # C
    7:  (43.09, 33.05),   # N
    8:  (49.28, 49.78),   # O
    9:  (71.8,  48.61),   # F
}


def _hydride_formula(symbol: str, n_H: int) -> str:
    """Generate chemical formula for a hydride: XHn.

    Uses conventional naming: H₂O not OH₂, HF not FH, HCl not ClH, HBr not BrH.
    """
    # Halogen and chalcogen hydrides: H comes first by convention
    _H_FIRST = {'O', 'S', 'Se', 'Te', 'F', 'Cl', 'Br', 'I'}
    if symbol in _H_FIRST:
        if n_H == 0:
            return symbol
        elif n_H == 1:
            return f'H{symbol}'
        else:
            return f'H{n_H}{symbol}'
    else:
        if n_H == 0:
            return symbol
        elif n_H == 1:
            return f'{symbol}H' if symbol != 'H' else 'H2'
        else:
            return f'{symbol}H{n_H}'


def _n_hydrogens(n_valence: int, structure_type: str) -> int:
    """Number of hydrogen atoms in a standard hydride.

    Alkali (1 val) → 1H, alkaline earth (2 val) → 2H,
    p-block: C/Si/Ge → 4H, N/P/As → 3H, O/S/Se → 2H, F/Cl/Br → 1H.
    """
    if structure_type in ('C',):  # alkali
        return 1
    elif structure_type in ('D',):  # alkaline earth
        return 2
    else:  # E-type p-block: 8 - n_valence lone-pair electrons → n_H = 8 - n_val
        # n_val=4 → 4H, n_val=5 → 3H, n_val=6 → 2H, n_val=7 → 1H
        return 8 - n_valence


def hydride_spec(
    Z: int,
    R: Optional[float] = None,
    max_n: int = 2,
    include_pk: bool = True,
    core_method: str = 'pk',
) -> MolecularSpec:
    """Create a MolecularSpec for a main-group hydride XHn.

    Handles first-row (Z=3-9, explicit 1s² core), second-row (Z=11-17,
    [Ne] frozen core), third-row s-block (Z=19-20, [Ar] frozen core),
    and third-row p-block (Z=31-36, [Ar]3d¹⁰ frozen core).

    Block structure (automatic from atomic classifier):
      - Core block (first-row only): Z_heavy, 2 electrons
      - Bond blocks: Z_eff center + H partner, 2 electrons each
      - Lone pair blocks: Z_eff center only, 2 electrons each
      - PK on all valence center-side orbitals (first-row only)

    Parameters
    ----------
    Z : int
        Nuclear charge of the heavy atom.
    R : float, optional
        Bond distance in bohr. Defaults to experimental equilibrium.
    max_n : int
        Maximum principal quantum number for all blocks (default 2).
    include_pk : bool
        If True (default), include PK pseudopotential on first-row
        valence blocks. Frozen-core molecules never have PK.

    Returns
    -------
    MolecularSpec
    """
    if Z not in _ELEMENT_DATA:
        raise ValueError(
            f"Z={Z} not supported. Supported: {sorted(_ELEMENT_DATA.keys())}"
        )

    from geovac.atomic_classifier import classify_atom

    symbol, period, n_core_e, core_type = _ELEMENT_DATA[Z]
    cls = classify_atom(Z)
    n_valence = cls.n_valence_electrons
    Z_eff = cls.Z_eff_valence
    n_H = _n_hydrogens(n_valence, cls.structure_type)
    n_bonds = n_H
    n_lone_pairs = (n_valence - n_bonds) // 2

    formula = _hydride_formula(symbol, n_H)

    if R is None:
        R = _HYDRIDE_REQ.get(formula)
        if R is None:
            raise ValueError(
                f"No default R_eq for {formula}. Provide R explicitly."
            )

    # PK parameters (first-row only)
    pk_A, pk_B = 0.0, 0.0
    if include_pk and core_type == 'explicit' and Z in _PK_PARAMS:
        pk_A, pk_B = _PK_PARAMS[Z]

    # --- Nuclear repulsion constant ---
    Z_H = 1.0
    if core_type == 'explicit':
        # First-row: He-like core energy
        E_core = _FIRST_ROW_CORE_ENERGY.get(Z, -(Z - 5.0 / 16) ** 2)
        V_NN = float(Z) * Z_H * n_H / R
        nuclear_repulsion = V_NN + E_core
    else:
        # Frozen core: FrozenCore energy + core-H electrostatic
        from geovac.neon_core import FrozenCore
        from geovac.composed_qubit import _v_cross_nuc_frozen_core

        fc = FrozenCore(Z)
        fc.solve()
        E_core = fc.energy
        V_NN = float(Z) * Z_H * n_H / R
        V_cross = _v_cross_nuc_frozen_core(Z, Z_H, R) * n_H
        nuclear_repulsion = V_NN + V_cross + E_core

    # --- Build blocks ---
    blocks: List[OrbitalBlock] = []

    # Core block (first-row only)
    if core_type == 'explicit':
        blocks.append(OrbitalBlock(
            label=f'{symbol}_core',
            block_type='core',
            Z_center=float(Z),
            n_electrons=2,
            max_n=max_n,
        ))

    # Bond blocks
    for i in range(n_bonds):
        suffix = f'_{i+1}' if n_bonds > 1 else ''
        blocks.append(OrbitalBlock(
            label=f'{formula}_bond{suffix}',
            block_type='bond',
            Z_center=float(Z_eff),
            n_electrons=2,
            max_n=max_n,
            has_h_partner=True,
            Z_partner=Z_H,
            max_n_partner=max_n,
            pk_A=pk_A,
            pk_B=pk_B,
        ))

    # Lone pair blocks
    for i in range(n_lone_pairs):
        suffix = f'_{i+1}' if n_lone_pairs > 1 else ''
        blocks.append(OrbitalBlock(
            label=f'{symbol}_lp{suffix}',
            block_type='lone_pair',
            Z_center=float(Z_eff),
            n_electrons=2,
            max_n=max_n,
            pk_A=pk_A,
            pk_B=pk_B,
        ))

    # Description
    n_encoded = 2 * n_bonds + 2 * n_lone_pairs + (2 if core_type == 'explicit' else 0)
    core_desc = f'explicit 1s² core' if core_type == 'explicit' else f'frozen core ({n_core_e}e)'
    desc = (
        f'{formula}: {core_desc}, '
        f'{n_bonds} bond(s) + {n_lone_pairs} LP(s), '
        f'{n_encoded} electrons encoded'
    )

    return MolecularSpec(
        name=formula,
        blocks=blocks,
        nuclear_repulsion_constant=nuclear_repulsion,
        description=desc,
        core_method=core_method,
    )


# ---------------------------------------------------------------------------
# Main-group hydride convenience aliases
# ---------------------------------------------------------------------------

# First row
def lih_spec(**kw) -> MolecularSpec:
    """LiH spec."""
    return hydride_spec(3, **kw)

def beh2_spec(**kw) -> MolecularSpec:
    """BeH₂ spec."""
    return hydride_spec(4, **kw)

def ch4_spec(**kw) -> MolecularSpec:
    """CH₄ spec."""
    return hydride_spec(6, **kw)

def nh3_spec(**kw) -> MolecularSpec:
    """NH₃ spec."""
    return hydride_spec(7, **kw)

def h2o_spec(**kw) -> MolecularSpec:
    """H₂O spec."""
    return hydride_spec(8, **kw)

def hf_spec(**kw) -> MolecularSpec:
    """HF spec."""
    return hydride_spec(9, **kw)


# Second row
def nah_spec(**kw) -> MolecularSpec:
    """NaH spec."""
    return hydride_spec(11, **kw)

def mgh2_spec(**kw) -> MolecularSpec:
    """MgH₂ spec."""
    return hydride_spec(12, **kw)

def sih4_spec(**kw) -> MolecularSpec:
    """SiH₄ spec."""
    return hydride_spec(14, **kw)

def ph3_spec(**kw) -> MolecularSpec:
    """PH₃ spec."""
    return hydride_spec(15, **kw)

def h2s_spec(**kw) -> MolecularSpec:
    """H₂S spec."""
    return hydride_spec(16, **kw)

def hcl_spec(**kw) -> MolecularSpec:
    """HCl spec."""
    return hydride_spec(17, **kw)


# Third row s-block
def kh_spec(**kw) -> MolecularSpec:
    """KH spec."""
    return hydride_spec(19, **kw)

def cah2_spec(**kw) -> MolecularSpec:
    """CaH₂ spec."""
    return hydride_spec(20, **kw)


# Third row p-block
def geh4_spec(**kw) -> MolecularSpec:
    """GeH₄ spec."""
    return hydride_spec(32, **kw)

def ash3_spec(**kw) -> MolecularSpec:
    """AsH₃ spec."""
    return hydride_spec(33, **kw)

def h2se_spec(**kw) -> MolecularSpec:
    """H₂Se spec."""
    return hydride_spec(34, **kw)

def hbr_spec(**kw) -> MolecularSpec:
    """HBr spec."""
    return hydride_spec(35, **kw)


# ---------------------------------------------------------------------------
# Multi-center molecule spec factories
# ---------------------------------------------------------------------------

# Experimental equilibrium bond lengths in bohr
_MULTI_CENTER_REQ = {
    'LiF': 2.955, 'CO': 2.132, 'N2': 2.074, 'F2': 2.668,
    'NaCl': 4.461, 'CH2O': 2.273, 'C2H2': 2.273, 'C2H6': 2.896,
}


def _diatomic_multi_center_spec(
    name: str,
    Z_A: int, Z_B: int,
    R: Optional[float] = None,
    max_n: int = 2,
) -> MolecularSpec:
    """Create a MolecularSpec for a homonuclear or heteronuclear diatomic
    with two heavy-atom centers (e.g. CO, N₂, F₂, LiF, NaCl).

    Each heavy atom contributes: core (if first-row) + bond + lone pairs.
    """
    from geovac.atomic_classifier import classify_atom

    if R is None:
        R = _MULTI_CENTER_REQ.get(name)
        if R is None:
            raise ValueError(f"No default R_eq for {name}. Provide R.")

    cls_A = classify_atom(Z_A)
    cls_B = classify_atom(Z_B)
    sym_A = _ELEMENT_DATA[Z_A][0] if Z_A in _ELEMENT_DATA else f'Z{Z_A}'
    sym_B = _ELEMENT_DATA[Z_B][0] if Z_B in _ELEMENT_DATA else f'Z{Z_B}'

    n_val_A = cls_A.n_valence_electrons
    n_val_B = cls_B.n_valence_electrons
    Z_eff_A = cls_A.Z_eff_valence
    Z_eff_B = cls_B.Z_eff_valence
    n_core_A = cls_A.n_core_electrons
    n_core_B = cls_B.n_core_electrons
    core_type_A = _ELEMENT_DATA.get(Z_A, ('?', 0, 0, 'frozen'))[3]
    core_type_B = _ELEMENT_DATA.get(Z_B, ('?', 0, 0, 'frozen'))[3]

    # Bond order and lone pair counts from standard chemistry
    # n_bonds shared bonds (2e each), remaining valence → lone pairs
    _BOND_DATA = {
        'LiF':  (1, 0, 3),   # 1 bond, 0 LP on Li, 3 LP on F
        'CO':   (3, 1, 1),   # triple bond, 1 LP on C, 1 LP on O
        'N2':   (3, 1, 1),   # triple bond, 1 LP per N
        'F2':   (1, 3, 3),   # single bond, 3 LP per F
        'NaCl': (1, 0, 3),   # 1 bond, 0 LP on Na, 3 LP on Cl
    }
    if name in _BOND_DATA:
        n_bonds, n_lp_A, n_lp_B = _BOND_DATA[name]
    else:
        # Default: single bond, remaining as lone pairs
        n_bonds = 1
        n_lp_A = (n_val_A - 1) // 2
        n_lp_B = (n_val_B - 1) // 2

    # PK parameters (first-row explicit core only)
    pk_A_a, pk_B_a = 0.0, 0.0
    pk_A_b, pk_B_b = 0.0, 0.0
    if core_type_A == 'explicit' and Z_A in _PK_PARAMS:
        pk_A_a, pk_B_a = _PK_PARAMS[Z_A]
    if core_type_B == 'explicit' and Z_B in _PK_PARAMS:
        pk_A_b, pk_B_b = _PK_PARAMS[Z_B]

    blocks: List[OrbitalBlock] = []

    # Atom A core
    if core_type_A == 'explicit':
        blocks.append(OrbitalBlock(
            label=f'{sym_A}_core', block_type='core',
            Z_center=float(Z_A), n_electrons=2, max_n=max_n,
            center_nucleus_idx=0,
        ))

    # Shared bond blocks (center on A, partner on B)
    for i in range(n_bonds):
        suffix = f'_{i+1}' if n_bonds > 1 else ''
        blocks.append(OrbitalBlock(
            label=f'{name}_bond{suffix}', block_type='bond',
            Z_center=float(Z_eff_A), n_electrons=2, max_n=max_n,
            has_h_partner=True, Z_partner=float(Z_eff_B),
            max_n_partner=max_n,
            pk_A=pk_A_a, pk_B=pk_B_a,
            center_nucleus_idx=0, partner_nucleus_idx=1,
        ))

    # Lone pairs on A
    for i in range(n_lp_A):
        suffix = f'_{i+1}' if n_lp_A > 1 else ''
        blocks.append(OrbitalBlock(
            label=f'{sym_A}_lp{suffix}', block_type='lone_pair',
            Z_center=float(Z_eff_A), n_electrons=2, max_n=max_n,
            pk_A=pk_A_a, pk_B=pk_B_a,
            center_nucleus_idx=0,
        ))

    # Atom B core
    if core_type_B == 'explicit':
        blocks.append(OrbitalBlock(
            label=f'{sym_B}_core', block_type='core',
            Z_center=float(Z_B), n_electrons=2, max_n=max_n,
            center_nucleus_idx=1,
        ))

    # Lone pairs on B
    for i in range(n_lp_B):
        suffix = f'_{i+1}' if n_lp_B > 1 else ''
        blocks.append(OrbitalBlock(
            label=f'{sym_B}_lp{suffix}', block_type='lone_pair',
            Z_center=float(Z_eff_B), n_electrons=2, max_n=max_n,
            pk_A=pk_A_b, pk_B=pk_B_b,
            center_nucleus_idx=1,
        ))

    # Nuclear repulsion
    V_NN = float(Z_A * Z_B) / R
    # Core energies
    E_core_A = _FIRST_ROW_CORE_ENERGY.get(Z_A, 0.0) if core_type_A == 'explicit' else 0.0
    E_core_B = _FIRST_ROW_CORE_ENERGY.get(Z_B, 0.0) if core_type_B == 'explicit' else 0.0

    # Frozen core contributions
    if core_type_A == 'frozen':
        from geovac.neon_core import FrozenCore
        from geovac.composed_qubit import _v_cross_nuc_frozen_core
        fc = FrozenCore(Z_A)
        fc.solve()
        E_core_A = fc.energy
        V_NN += _v_cross_nuc_frozen_core(Z_A, float(Z_B), R)
    if core_type_B == 'frozen':
        from geovac.neon_core import FrozenCore
        from geovac.composed_qubit import _v_cross_nuc_frozen_core
        fc = FrozenCore(Z_B)
        fc.solve()
        E_core_B = fc.energy
        V_NN += _v_cross_nuc_frozen_core(Z_B, float(Z_A), R)

    nuclear_repulsion = V_NN + E_core_A + E_core_B

    nuclei = [
        {'Z': float(Z_A), 'position': (0.0, 0.0, 0.0), 'label': sym_A},
        {'Z': float(Z_B), 'position': (0.0, 0.0, R), 'label': sym_B},
    ]

    return MolecularSpec(
        name=name,
        blocks=blocks,
        nuclear_repulsion_constant=nuclear_repulsion,
        description=f'{name}: {sym_A}(Z={Z_A}) + {sym_B}(Z={Z_B}), R={R:.3f} bohr',
        nuclei=nuclei,
    )


# Convenience factories
def lif_spec(R: Optional[float] = None, **kw) -> MolecularSpec:
    """LiF spec."""
    return _diatomic_multi_center_spec('LiF', 3, 9, R=R, **kw)

def co_spec(R: Optional[float] = None, **kw) -> MolecularSpec:
    """CO spec."""
    return _diatomic_multi_center_spec('CO', 6, 8, R=R, **kw)

def n2_spec(R: Optional[float] = None, **kw) -> MolecularSpec:
    """N₂ spec."""
    return _diatomic_multi_center_spec('N2', 7, 7, R=R, **kw)

def f2_spec(R: Optional[float] = None, **kw) -> MolecularSpec:
    """F₂ spec."""
    return _diatomic_multi_center_spec('F2', 9, 9, R=R, **kw)

def nacl_spec(R: Optional[float] = None, **kw) -> MolecularSpec:
    """NaCl spec."""
    return _diatomic_multi_center_spec('NaCl', 11, 17, R=R, **kw)
