"""
Composed Qubit Hamiltonians (LiH, BeH2)
========================================

Constructs second-quantized Hamiltonians for composed natural geometry
systems (Paper 17) and performs Jordan-Wigner transformation to count
Pauli terms.

**LiH** (core-valence diatomic):
    Core:    2 electrons on Li (Z=3), hydrogenic orbitals up to max_n_core
    Valence: screened Li orbitals (Z_eff=1) + H orbitals (Z=1) up to max_n_val

**BeH2** (linear triatomic, Approach A -- independent bond blocks):
    Core:   2 electrons on Be (Z=4), hydrogenic orbitals up to max_n_core
    Bond 1: screened Be orbitals (Z_eff=2) + H1 orbitals (Z=1) up to max_n_val
    Bond 2: screened Be orbitals (Z_eff=2) + H2 orbitals (Z=1) up to max_n_val

The pseudopotential approximation replaces explicit two-body core-valence
Coulomb with one-body Z_eff screening, so cross-block ERIs are zero.

Compatible with qubit_encoding.build_fermion_op_from_integrals() via
chemist-notation ERIs (pq|rs).

Author: GeoVac Development Team
Date: March 2026
"""

import json
import time
from math import factorial
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.special import genlaguerre

from openfermion import jordan_wigner

from geovac.qubit_encoding import build_fermion_op_from_integrals


# ---------------------------------------------------------------------------
# Default PK parameters from Paper 17, Table 1 (Li²⁺ core, inv2 method)
# These are used when AbInitioPK is not run (avoids expensive hyperspherical
# solve for quick qubit-encoding analysis).
# ---------------------------------------------------------------------------
_PK_DEFAULTS = {
    3: {'A': 6.93, 'B': 7.00, 'source': 'Paper 17 Table 1, Li2+ inv2'},
    4: {'A': 13.01, 'B': 12.53, 'source': 'Paper 17 Table 1, Be3+ inv2'},
}

# PK parameters for He-like cores (2 core electrons), Z²-scaled from Li²⁺.
# A scales with core-valence energy gap ∝ Z², B scales with ⟨1/r²⟩_core ∝ Z².
# Source: Li²⁺ (Z=3) values from Paper 17 Table 1, scaled by (Z/3)².
_PK_HELIKE_DEFAULTS = {
    3: {'A': 6.93, 'B': 7.00, 'source': 'Paper 17 Table 1, Li2+ inv2'},
    4: {'A': 12.32, 'B': 12.44,
        'source': 'Z²-scaled from Li2+ (Paper 17): A=6.93*(4/3)²=12.32, B=7.00*(4/3)²=12.44'},
    5: {'A': 21.40, 'B': 18.46,
        'source': 'Ab initio (Track BI, Level 3 inv2)'},
    6: {'A': 31.37, 'B': 25.54,
        'source': 'Ab initio (Track BI, Level 3 inv2)'},
    7: {'A': 43.09, 'B': 33.05,
        'source': 'Ab initio (Track BI, Level 3 inv2)'},
    8: {'A': 49.28, 'B': 49.78,
        'source': 'Z²-scaled from Li2+ (Paper 17): A=6.93*(8/3)²=49.28, B=7.00*(8/3)²=49.78'},
    9: {'A': 71.80, 'B': 48.61,
        'source': 'Ab initio (Track BI, Level 3 inv2)'},
}


# ---------------------------------------------------------------------------
# Published Gaussian Pauli term counts (Trenev et al., Quantum 2025)
# ---------------------------------------------------------------------------
# Source: D. Trenev, P. J. Ollitrault, S. M. Harwood, T. P. Gujarati,
#         S. Raman, A. Mezzacapo, S. Mostame,
#         "Refining resource estimation for the quantum computation of
#         molecular spectra through Trotter error analysis,"
#         Quantum (2025). arXiv:2311.03719.
#         Table 5, Jordan-Wigner encoding with 2-qubit reduction.
# ---------------------------------------------------------------------------

GAUSSIAN_LIH_PUBLISHED: Dict[str, Dict[str, Any]] = {
    'sto-3g': {
        'Q': 10, 'N_pauli': 276,
        'note': '2-qubit reduction from 12 -> 10 qubits',
    },
    '6-31g': {
        'Q': 20, 'N_pauli': 5851,
        'note': '2-qubit reduction',
    },
    'cc-pvdz': {
        'Q': 36, 'N_pauli': 63519,
        'note': '2-qubit reduction',
    },
}

GAUSSIAN_H2O_PUBLISHED: Dict[str, Dict[str, Any]] = {
    'sto-3g': {
        'Q': 12, 'N_pauli': 551,
        'note': '2-qubit reduction from 14 -> 12 qubits',
    },
    '6-31g': {
        'Q': 24, 'N_pauli': 8921,
        'note': '2-qubit reduction',
    },
    'cc-pvdz': {
        'Q': 46, 'N_pauli': 107382,
        'note': '2-qubit reduction',
    },
}

_TRENEV_REFERENCE = (
    "D. Trenev, P. J. Ollitrault, S. M. Harwood, T. P. Gujarati, "
    "S. Raman, A. Mezzacapo, S. Mostame, "
    "\"Refining resource estimation for the quantum computation of "
    "molecular spectra through Trotter error analysis,\" "
    "Quantum (2025). arXiv:2311.03719."
)


# ---------------------------------------------------------------------------
# Hydrogenic state enumeration
# ---------------------------------------------------------------------------

def _enumerate_states(max_n: int) -> List[Tuple[int, int, int]]:
    """Generate (n, l, m) states up to max_n in canonical order."""
    states: List[Tuple[int, int, int]] = []
    for n in range(1, max_n + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                states.append((n, l, m))
    return states


# ---------------------------------------------------------------------------
# Wigner 3j symbol (self-contained, matches LatticeIndex implementation)
# ---------------------------------------------------------------------------

def _wigner3j(j1: int, j2: int, j3: int,
              m1: int, m2: int, m3: int) -> float:
    """Wigner 3j symbol for integer arguments via Racah formula."""
    if m1 + m2 + m3 != 0:
        return 0.0
    if abs(j1 - j2) > j3 or j3 > j1 + j2:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0

    def _tri(a: int, b: int, c: int) -> float:
        return (factorial(a + b - c) * factorial(a - b + c)
                * factorial(-a + b + c)
                / factorial(a + b + c + 1))

    pre = ((-1) ** (j1 - j2 - m3)
           * np.sqrt(_tri(j1, j2, j3)
                     * factorial(j1 + m1) * factorial(j1 - m1)
                     * factorial(j2 + m2) * factorial(j2 - m2)
                     * factorial(j3 + m3) * factorial(j3 - m3)))

    t_min = max(0, j2 - j3 - m1, j1 - j3 + m2)
    t_max = min(j1 + j2 - j3, j1 - m1, j2 + m2)

    s = 0.0
    for t in range(t_min, t_max + 1):
        s += ((-1) ** t
              / (factorial(t)
                 * factorial(j3 - j2 + m1 + t)
                 * factorial(j3 - j1 - m2 + t)
                 * factorial(j1 + j2 - j3 - t)
                 * factorial(j1 - m1 - t)
                 * factorial(j2 + m2 - t)))
    return pre * s


# ---------------------------------------------------------------------------
# Gaunt angular coupling coefficient
# ---------------------------------------------------------------------------

def _ck_coefficient(la: int, ma: int, lc: int, mc: int, k: int) -> float:
    """
    c^k(l,m,l',m') = (-1)^m sqrt((2l+1)(2l'+1)) * (l k l'; 0 0 0) * (l k l'; -m q m')
    where q = mc - ma.
    """
    q = mc - ma
    pre = ((-1) ** ma * np.sqrt((2 * la + 1) * (2 * lc + 1)))
    w1 = _wigner3j(la, k, lc, 0, 0, 0)
    if abs(w1) < 1e-15:
        return 0.0
    w2 = _wigner3j(la, k, lc, -ma, q, mc)
    return pre * w1 * w2


# ---------------------------------------------------------------------------
# Radial wavefunctions and Slater R^k integrals
# ---------------------------------------------------------------------------

def _radial_wf_grid(
    Z: float, n: int, l: int, r_grid: np.ndarray,
) -> np.ndarray:
    """Normalized hydrogenic radial wavefunction R_{nl}(r) on grid."""
    rho = 2.0 * Z * r_grid / n
    L_poly = genlaguerre(n - l - 1, 2 * l + 1)(rho)
    wf = rho ** l * np.exp(-rho / 2.0) * L_poly
    norm_sq = np.trapezoid(wf ** 2 * r_grid ** 2, r_grid)
    if norm_sq < 1e-30:
        return np.zeros_like(r_grid)
    return wf / np.sqrt(norm_sq)


def _compute_rk_integrals_block(
    Z: float,
    states: List[Tuple[int, int, int]],
    n_grid: int = 2000,
) -> Dict[Tuple[int, ...], float]:
    """
    Compute R^k(n1l1, n2l2, n3l3, n4l4) Slater radial integrals for a
    single-center block at nuclear charge Z.

    Returns dict mapping (n1,l1,n2,l2,n3,l3,n4,l4,k) -> float.
    """
    unique_nl = sorted(set((n, l) for n, l, m in states))
    r_max = 80.0 / max(Z, 0.5)
    r_grid = np.linspace(0, r_max, n_grid + 1)[1:]
    dr = r_grid[1] - r_grid[0]

    # Pre-compute radial wavefunctions
    R_on_grid: Dict[Tuple[int, int], np.ndarray] = {}
    for n, l in unique_nl:
        R_on_grid[(n, l)] = _radial_wf_grid(Z, n, l, r_grid)

    # Determine needed integrals
    needed: List[Tuple[int, ...]] = []
    for n1, l1 in unique_nl:
        for n2, l2 in unique_nl:
            for n3, l3 in unique_nl:
                for n4, l4 in unique_nl:
                    k_max = min(l1 + l3, l2 + l4)
                    for k in range(0, k_max + 1):
                        if (l1 + l3 + k) % 2 != 0:
                            continue
                        if (l2 + l4 + k) % 2 != 0:
                            continue
                        needed.append((n1, l1, n2, l2, n3, l3, n4, l4, k))

    if not needed:
        return {}

    # Group by (n2,l2,n4,l4,k) to reuse Y^k potential
    from collections import defaultdict
    yk_groups: Dict[Tuple[int, ...], List[Tuple[int, ...]]] = defaultdict(list)
    for key in needed:
        n1, l1, n2, l2, n3, l3, n4, l4, k = key
        yk_key = (n2, l2, n4, l4, k)
        yk_groups[yk_key].append(key)

    # Compute Y^k potentials
    yk_cache: Dict[Tuple[int, ...], np.ndarray] = {}
    for yk_key in yk_groups:
        n2, l2, n4, l4, k = yk_key
        f_r = R_on_grid[(n2, l2)] * R_on_grid[(n4, l4)] * r_grid ** 2
        yk = np.zeros(n_grid)
        for i in range(n_grid):
            r1 = r_grid[i]
            inner = np.sum(f_r[:i + 1] * (r_grid[:i + 1] / r1) ** k) * dr
            if i + 1 < n_grid:
                outer = np.sum(
                    f_r[i + 1:]
                    * (r1 / r_grid[i + 1:]) ** k
                    / r_grid[i + 1:]
                ) * dr
            else:
                outer = 0.0
            yk[i] = inner / r1 + outer
        yk_cache[yk_key] = yk

    # Compute R^k integrals
    rk_cache: Dict[Tuple[int, ...], float] = {}
    for yk_key, keys in yk_groups.items():
        n2, l2, n4, l4, k = yk_key
        yk = yk_cache[yk_key]
        for key in keys:
            n1, l1, _, _, n3, l3, _, _, _ = key
            integrand = (R_on_grid[(n1, l1)] * R_on_grid[(n3, l3)]
                         * yk * r_grid ** 2)
            val = np.trapezoid(integrand, r_grid)
            rk_cache[key] = val

    return rk_cache


def _build_eri_block(
    Z: float,
    states: List[Tuple[int, int, int]],
    rk_cache: Dict[Tuple[int, ...], float],
) -> Dict[Tuple[int, int, int, int], float]:
    """
    Build ERI dict in PHYSICIST notation <ab|cd> for a single-center block.

    Uses Gaunt integral selection rules to enforce sparsity.
    Returns dict mapping (a, b, c, d) -> float (spatial orbital indices
    relative to the block).
    """
    n_sp = len(states)
    eri: Dict[Tuple[int, int, int, int], float] = {}

    # Pre-compute c^k table
    ck_table: Dict[Tuple[int, int, int], float] = {}
    for a in range(n_sp):
        la, ma = states[a][1], states[a][2]
        for c in range(n_sp):
            lc, mc = states[c][1], states[c][2]
            k_max = la + lc
            for k in range(0, k_max + 1):
                if (la + lc + k) % 2 != 0:
                    continue
                val = _ck_coefficient(la, ma, lc, mc, k)
                if abs(val) > 1e-15:
                    ck_table[(a, c, k)] = val

    # Group by (a, c) pair
    from collections import defaultdict
    ac_k_map: Dict[Tuple[int, int], List[Tuple[int, float]]] = defaultdict(list)
    for (a, c, k), val in ck_table.items():
        ac_k_map[(a, c)].append((k, val))

    for (a, c), ck_ac_list in ac_k_map.items():
        na, la, ma = states[a]
        nc, lc, mc = states[c]
        for (b, d), ck_bd_list in ac_k_map.items():
            nb, lb, mb = states[b]
            nd, ld, md = states[d]

            if ma + mb != mc + md:
                continue

            val = 0.0
            for k_ac, c_ac in ck_ac_list:
                for k_bd, c_bd in ck_bd_list:
                    if k_ac != k_bd:
                        continue
                    k = k_ac
                    rk_key = (na, la, nb, lb, nc, lc, nd, ld, k)
                    rk_val = rk_cache.get(rk_key)
                    if rk_val is None:
                        continue
                    val += c_ac * c_bd * rk_val

            if abs(val) > 1e-15:
                eri[(a, b, c, d)] = val

    return eri


def _physicist_to_chemist(
    eri_phys: Dict[Tuple[int, int, int, int], float],
    n_spatial: int,
    offset: int = 0,
) -> np.ndarray:
    """
    Convert physicist-notation ERI dict to chemist-notation dense tensor.

    Physicist: <ab|cd> = integral phi_a(1) phi_b(2) 1/r12 phi_c(1) phi_d(2)
    Chemist:  (pq|rs) = integral phi_p(1) phi_q(1) 1/r12 phi_r(2) phi_s(2)

    Relation: (pq|rs) = <pr|qs>  =>  chemist[p,q,r,s] = phys[p,r,q,s]
    Or equivalently: phys[a,b,c,d] -> chemist[a,c,b,d]
    """
    eri_chem = np.zeros((n_spatial, n_spatial, n_spatial, n_spatial))
    for (a, b, c, d), val in eri_phys.items():
        # phys <ab|cd> -> chemist (ac|bd)
        p, q, r, s = a + offset, c + offset, b + offset, d + offset
        eri_chem[p, q, r, s] = val
    return eri_chem


# ---------------------------------------------------------------------------
# PK pseudopotential matrix elements (Paper 17 Sec IV)
# ---------------------------------------------------------------------------

def _compute_pk_matrix_elements(
    Z_eff: float,
    states: List[Tuple[int, int, int]],
    A_pk: float,
    B_pk: float,
    n_grid: int = 4000,
) -> np.ndarray:
    """
    Compute <φ_p | V_PK(r) | φ_q> for same-center hydrogenic orbitals.

    V_PK(r) = A * exp(-B*r²) / r² is spherically symmetric, so the angular
    integral gives δ_{l_p,l_q} δ_{m_p,m_q}. Only the radial integral remains:

        <n_p l m | V_PK | n_q l m> = A ∫₀^∞ R_{n_p,l}(r;Z) R_{n_q,l}(r;Z)
                                      x exp(-B r²) dr

    Note: the r² from the volume element cancels the 1/r² in V_PK.

    Parameters
    ----------
    Z_eff : float
        Effective nuclear charge for the orbital wavefunctions.
    states : list of (n, l, m)
        Orbital quantum numbers.
    A_pk, B_pk : float
        PK barrier height (Ha·bohr²) and width exponent (bohr⁻²).
    n_grid : int
        Number of radial grid points.

    Returns
    -------
    h1_pk : ndarray of shape (M, M)
        PK matrix elements. Diagonal in (l, m), couples different n.
    """
    M = len(states)
    h1_pk = np.zeros((M, M))

    # Group states by (l, m) to exploit selection rule
    from collections import defaultdict
    lm_groups: Dict[Tuple[int, int], List[Tuple[int, int]]] = defaultdict(list)
    for idx, (n, l, m) in enumerate(states):
        lm_groups[(l, m)].append((idx, n))

    # Set up radial grid (extend far enough for diffuse orbitals)
    r_max = 80.0 / max(Z_eff, 0.3)
    r_grid = np.linspace(0, r_max, n_grid + 1)[1:]  # skip r=0

    # V_PK integrand weight: A * exp(-B*r²)  (the r²/r² cancels)
    vpk_weight = A_pk * np.exp(-B_pk * r_grid**2)

    # Pre-compute radial wavefunctions
    wf_cache: Dict[Tuple[int, int], np.ndarray] = {}
    for (l, m), group in lm_groups.items():
        for idx, n in group:
            if (n, l) not in wf_cache:
                wf_cache[(n, l)] = _radial_wf_grid(Z_eff, n, l, r_grid)

    # Compute matrix elements within each (l,m) block
    for (l, m), group in lm_groups.items():
        for i, (idx_p, n_p) in enumerate(group):
            R_p = wf_cache[(n_p, l)]
            for idx_q, n_q in group[i:]:  # exploit symmetry
                R_q = wf_cache[(n_q, l)]
                integrand = R_p * R_q * vpk_weight
                val = float(np.trapezoid(integrand, r_grid))
                h1_pk[idx_p, idx_q] = val
                h1_pk[idx_q, idx_p] = val  # Hermitian

    return h1_pk


# ---------------------------------------------------------------------------
# Cross-nuclear attraction (Paper 17 Eq. 4)
# ---------------------------------------------------------------------------

def _v_cross_nuc_1s(Z_core: float, n_core: int, Z_other: float,
                    R: float) -> float:
    """
    Analytical core-to-other-nucleus attraction for 1s^n_core core.

    V_cross = -n_core * Z_other * <1s_Z|1/r_other|1s_Z>
    where <1s_Z|1/r_B|1s_Z> = (1/R) * [1 - (1 + Z*R) * exp(-2*Z*R)]
    """
    zr = Z_core * R
    expectation = (1.0 / R) * (1.0 - (1.0 + zr) * np.exp(-2.0 * zr))
    return -n_core * Z_other * expectation


# ---------------------------------------------------------------------------
# General composed-geometry builder  (Track BH)
# ---------------------------------------------------------------------------

from geovac.molecular_spec import MolecularSpec, OrbitalBlock


def build_composed_hamiltonian(
    spec: MolecularSpec,
    pk_in_hamiltonian: bool = True,
    verbose: bool = False,
) -> Dict[str, Any]:
    """
    General builder for composed-geometry qubit Hamiltonians.

    Consumes a :class:`MolecularSpec` that describes the block structure and
    produces h1, ERI, fermion operator, and JW-transformed qubit operator.

    Parameters
    ----------
    spec : MolecularSpec
        Molecular specification (blocks, nuclear repulsion, etc.).
    pk_in_hamiltonian : bool
        If True (default), PK matrix elements are added to h1 for the
        quantum Hamiltonian.  If False, PK is computed and returned
        separately in ``h1_pk`` for classical 1-RDM post-processing.
    verbose : bool
        Print progress messages.

    Returns
    -------
    dict
        Keys include ``M``, ``Q``, ``N_pauli``, ``h1``, ``h1_pk``, ``eri``,
        ``nuclear_repulsion``, ``qubit_op``, ``fermion_op``, ``wall_time_s``,
        ``blocks``, and per-block metadata.
    """
    t0 = time.perf_counter()

    # ------------------------------------------------------------------
    # 1. Enumerate states per block and compute offsets
    # ------------------------------------------------------------------
    block_info: List[Dict[str, Any]] = []
    offset = 0
    for blk in spec.blocks:
        info: Dict[str, Any] = {
            'label': blk.label,
            'block_type': blk.block_type,
            'Z_center': blk.Z_center,
            'n_electrons': blk.n_electrons,
        }
        center_states = _enumerate_states(blk.max_n)
        info['center_states'] = center_states
        info['center_offset'] = offset
        info['center_M'] = len(center_states)
        offset += len(center_states)

        if blk.has_h_partner:
            partner_n = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            partner_states = _enumerate_states(partner_n)
            info['partner_states'] = partner_states
            info['partner_offset'] = offset
            info['partner_M'] = len(partner_states)
            info['Z_partner'] = blk.Z_partner
            offset += len(partner_states)
        else:
            info['partner_states'] = None
            info['partner_offset'] = None
            info['partner_M'] = 0
            info['Z_partner'] = None

        info['pk_A'] = blk.pk_A
        info['pk_B'] = blk.pk_B
        block_info.append(info)

    M = offset
    Q = 2 * M

    if verbose:
        print(f"[build_composed_hamiltonian] {spec.name}: M={M}, Q={Q}")

    # ------------------------------------------------------------------
    # 2. Build h1  (M x M diagonal)
    # ------------------------------------------------------------------
    h1 = np.zeros((M, M))
    for bi in block_info:
        # Center orbitals
        off_c = bi['center_offset']
        Z_c = bi['Z_center']
        for i, (n, l, m) in enumerate(bi['center_states']):
            h1[off_c + i, off_c + i] = -Z_c ** 2 / (2.0 * n ** 2)

        # Partner orbitals (if any)
        if bi['partner_states'] is not None:
            off_p = bi['partner_offset']
            Z_p = bi['Z_partner']
            for i, (n, l, m) in enumerate(bi['partner_states']):
                h1[off_p + i, off_p + i] = -Z_p ** 2 / (2.0 * n ** 2)

    # ------------------------------------------------------------------
    # 3. PK pseudopotential  (always on center orbitals of blocks with pk_A>0)
    # ------------------------------------------------------------------
    h1_pk_full = np.zeros((M, M))
    n_pk_nonzero = 0
    for bi in block_info:
        if bi['pk_A'] > 0.0:
            h1_pk_block = _compute_pk_matrix_elements(
                bi['Z_center'], bi['center_states'],
                bi['pk_A'], bi['pk_B'],
            )
            off_c = bi['center_offset']
            M_c = bi['center_M']
            for i in range(M_c):
                for j in range(M_c):
                    if abs(h1_pk_block[i, j]) > 1e-15:
                        h1_pk_full[off_c + i, off_c + j] = h1_pk_block[i, j]
                        n_pk_nonzero += 1

    if pk_in_hamiltonian:
        h1 = h1 + h1_pk_full

    # ------------------------------------------------------------------
    # 4. Build ERI  (M x M x M x M), block-diagonal in chemist notation
    # ------------------------------------------------------------------
    eri = np.zeros((M, M, M, M))

    # Cache: (Z, max_n) -> (rk_cache, eri_phys, states) to avoid recomputation
    rk_eri_cache: Dict[Tuple[float, int], Tuple[
        Dict[Tuple[int, ...], float],
        Dict[Tuple[int, int, int, int], float],
        List[Tuple[int, int, int]],
    ]] = {}

    def _get_rk_eri(Z: float, max_n: int) -> Tuple[
        Dict[Tuple[int, ...], float],
        Dict[Tuple[int, int, int, int], float],
        List[Tuple[int, int, int]],
    ]:
        key = (Z, max_n)
        if key not in rk_eri_cache:
            st = _enumerate_states(max_n)
            rk = _compute_rk_integrals_block(Z, st)
            ep = _build_eri_block(Z, st, rk)
            rk_eri_cache[key] = (rk, ep, st)
        return rk_eri_cache[key]

    for bi in block_info:
        # Center sub-block
        Z_c = bi['Z_center']
        off_c = bi['center_offset']
        max_n_c = len(set(n for n, l, m in bi['center_states']))  # infer from states
        # We need the max_n that was used to generate center_states.
        # Since _enumerate_states(max_n) generates all (n,l,m) up to max_n,
        # max_n = max(n for n,l,m in states).
        max_n_c = max(n for n, l, m in bi['center_states'])
        _, eri_phys_c, _ = _get_rk_eri(Z_c, max_n_c)

        for (a, b, c, d), val in eri_phys_c.items():
            # phys <ab|cd> -> chemist (ac|bd)
            eri[a + off_c, c + off_c, b + off_c, d + off_c] = val

        # Partner sub-block (if any)
        if bi['partner_states'] is not None:
            Z_p = bi['Z_partner']
            off_p = bi['partner_offset']
            max_n_p = max(n for n, l, m in bi['partner_states'])
            _, eri_phys_p, _ = _get_rk_eri(Z_p, max_n_p)

            for (a, b, c, d), val in eri_phys_p.items():
                eri[a + off_p, c + off_p, b + off_p, d + off_p] = val

    # Symmetrize: enforce (pq|rs) = (rs|pq) exactly
    eri = 0.5 * (eri + eri.transpose(2, 3, 0, 1))

    # ------------------------------------------------------------------
    # 5. JW transform
    # ------------------------------------------------------------------
    nuclear_repulsion = spec.nuclear_repulsion_constant

    if verbose:
        print(f"[build_composed_hamiltonian] Building fermion operator...")
    fermion_op = build_fermion_op_from_integrals(h1, eri, nuclear_repulsion)

    if verbose:
        print(f"[build_composed_hamiltonian] Jordan-Wigner transform...")
    qubit_op = jordan_wigner(fermion_op)
    N_pauli = len(qubit_op.terms)

    elapsed = time.perf_counter() - t0

    if verbose:
        print(f"[build_composed_hamiltonian] {spec.name}: Q={Q}, N_pauli={N_pauli}, "
              f"wall_time={elapsed:.1f}s")

    # ------------------------------------------------------------------
    # 6. Build result dict
    # ------------------------------------------------------------------
    n_eri_total = int(np.count_nonzero(np.abs(eri) > 1e-15))

    results: Dict[str, Any] = {
        'M': M,
        'Q': Q,
        'N_pauli': N_pauli,
        'nuclear_repulsion': nuclear_repulsion,
        'wall_time_s': elapsed,
        'n_eri_total': n_eri_total,
        'ERI_density_total': n_eri_total / max(1, M ** 4),
        'pk_in_hamiltonian': pk_in_hamiltonian,
        'n_pk_nonzero': n_pk_nonzero,
        'h1': h1,
        'h1_pk': h1_pk_full,
        'eri': eri,
        'qubit_op': qubit_op,
        'fermion_op': fermion_op,
        'blocks': [
            {
                'label': bi['label'],
                'block_type': bi['block_type'],
                'Z_center': bi['Z_center'],
                'center_offset': bi['center_offset'],
                'center_M': bi['center_M'],
                'partner_offset': bi['partner_offset'],
                'partner_M': bi['partner_M'],
                'Z_partner': bi['Z_partner'],
            }
            for bi in block_info
        ],
        'spec_name': spec.name,
    }

    return results


# ---------------------------------------------------------------------------
# Spec factory functions  (Track BH)
# ---------------------------------------------------------------------------

def lih_spec(
    max_n_core: int = 2,
    max_n_val: int = 2,
    R: float = 3.015,
    E_core: Optional[float] = None,
    include_pk: bool = True,
    A_pk: Optional[float] = None,
    B_pk: Optional[float] = None,
) -> MolecularSpec:
    """Create a MolecularSpec for LiH (core + bond pair).

    Parameters match :func:`build_composed_lih`.
    """
    Z_A = 3
    Z_B = 1
    n_core_electrons = 2
    Z_eff_val = Z_A - n_core_electrons  # 1

    if E_core is None:
        E_core = -7.2799

    # PK parameters
    if include_pk:
        if A_pk is None:
            A_pk = _PK_DEFAULTS[Z_A]['A']
        if B_pk is None:
            B_pk = _PK_DEFAULTS[Z_A]['B']
    else:
        A_pk = 0.0
        B_pk = 0.0

    # Nuclear repulsion constant
    V_NN = Z_A * Z_B / R
    V_cross = _v_cross_nuc_1s(Z_A, n_core_electrons, Z_B, R)
    nuclear_repulsion = V_NN + V_cross + E_core

    blocks = [
        OrbitalBlock(
            label='Li_core', block_type='core',
            Z_center=float(Z_A), n_electrons=2, max_n=max_n_core,
        ),
        OrbitalBlock(
            label='LiH_bond_Li', block_type='bond',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            has_h_partner=True, Z_partner=float(Z_B),
            max_n_partner=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
    ]
    return MolecularSpec('LiH', blocks, nuclear_repulsion)


def beh2_spec(
    max_n_core: int = 2,
    max_n_val: int = 2,
    R: float = 2.51,
    E_core: Optional[float] = None,
    include_pk: bool = True,
    A_pk: Optional[float] = None,
    B_pk: Optional[float] = None,
) -> MolecularSpec:
    """Create a MolecularSpec for BeH2 (core + 2 bond pairs).

    Parameters match :func:`build_composed_beh2`.
    """
    Z_A = 4
    Z_B = 1
    n_core_electrons = 2
    Z_eff_val = Z_A - n_core_electrons  # 2

    if E_core is None:
        E_core = -13.65

    # PK parameters -- uses _PK_HELIKE_DEFAULTS for Z=4
    if include_pk:
        if A_pk is None:
            A_pk = _PK_HELIKE_DEFAULTS[Z_A]['A']
        if B_pk is None:
            B_pk = _PK_HELIKE_DEFAULTS[Z_A]['B']
    else:
        A_pk = 0.0
        B_pk = 0.0

    # Nuclear repulsion: linear H-Be-H
    V_NN = 2.0 * Z_A * Z_B / R + Z_B * Z_B / (2.0 * R)
    V_cross = 2.0 * _v_cross_nuc_1s(Z_A, n_core_electrons, Z_B, R)
    nuclear_repulsion = V_NN + V_cross + E_core

    blocks = [
        OrbitalBlock(
            label='Be_core', block_type='core',
            Z_center=float(Z_A), n_electrons=2, max_n=max_n_core,
        ),
        OrbitalBlock(
            label='BeH2_bond1_Be', block_type='bond',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            has_h_partner=True, Z_partner=float(Z_B),
            max_n_partner=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
        OrbitalBlock(
            label='BeH2_bond2_Be', block_type='bond',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            has_h_partner=True, Z_partner=float(Z_B),
            max_n_partner=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
    ]
    return MolecularSpec('BeH2', blocks, nuclear_repulsion)


def h2o_spec(
    max_n_core: int = 2,
    max_n_val: int = 2,
    R_OH: float = 1.809,
    angle_HOH: float = 104.5,
    E_core: Optional[float] = None,
    include_pk: bool = True,
    A_pk: Optional[float] = None,
    B_pk: Optional[float] = None,
) -> MolecularSpec:
    """Create a MolecularSpec for H2O (core + 2 bonds + 2 lone pairs).

    Parameters match :func:`build_composed_h2o`.
    """
    Z_O = 8
    Z_H = 1
    n_core_electrons = 2
    Z_eff_val = Z_O - n_core_electrons  # 6

    if E_core is None:
        E_core = -(Z_O - 5.0 / 16.0) ** 2

    # PK parameters -- uses _PK_HELIKE_DEFAULTS for Z=8
    if include_pk:
        if A_pk is None:
            A_pk = _PK_HELIKE_DEFAULTS[Z_O]['A']
        if B_pk is None:
            B_pk = _PK_HELIKE_DEFAULTS[Z_O]['B']
    else:
        A_pk = 0.0
        B_pk = 0.0

    # Nuclear repulsion
    angle_rad = np.radians(angle_HOH)
    R_HH = 2.0 * R_OH * np.sin(angle_rad / 2.0)
    V_OH = 2.0 * Z_O * Z_H / R_OH
    V_HH = Z_H * Z_H / R_HH
    V_NN = V_OH + V_HH
    V_cross = 2.0 * _v_cross_nuc_1s(Z_O, n_core_electrons, Z_H, R_OH)
    nuclear_repulsion = V_NN + V_cross + E_core

    blocks = [
        OrbitalBlock(
            label='O_core', block_type='core',
            Z_center=float(Z_O), n_electrons=2, max_n=max_n_core,
        ),
        OrbitalBlock(
            label='H2O_bond1_O', block_type='bond',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            has_h_partner=True, Z_partner=float(Z_H),
            max_n_partner=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
        OrbitalBlock(
            label='H2O_bond2_O', block_type='bond',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            has_h_partner=True, Z_partner=float(Z_H),
            max_n_partner=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
        OrbitalBlock(
            label='H2O_lone1', block_type='lone_pair',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
        OrbitalBlock(
            label='H2O_lone2', block_type='lone_pair',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
    ]
    return MolecularSpec('H2O', blocks, nuclear_repulsion)


def h2_bond_pair_spec(
    max_n: int = 2,
    R: float = 1.4,
) -> MolecularSpec:
    """Create a MolecularSpec for H2 (single bond pair, no core).

    Parameters match :func:`build_h2_bond_pair`.
    """
    Z_eff = 1.0
    V_NN = 1.0 / R  # Z_A * Z_B / R = 1*1/R

    blocks = [
        OrbitalBlock(
            label='H2_bond_pair', block_type='bond_pair',
            Z_center=Z_eff, n_electrons=2, max_n=max_n,
        ),
    ]
    return MolecularSpec('H2', blocks, V_NN)


def he_spec(
    max_n: int = 2,
) -> MolecularSpec:
    """Create a MolecularSpec for He (single block, Z=2, no nuclear repulsion)."""
    blocks = [
        OrbitalBlock(
            label='He', block_type='bond_pair',
            Z_center=2.0, n_electrons=2, max_n=max_n,
        ),
    ]
    return MolecularSpec('He', blocks, 0.0)


def hf_spec(
    max_n_core: int = 2,
    max_n_val: int = 2,
    R: float = 1.733,
    E_core: Optional[float] = None,
    include_pk: bool = True,
    A_pk: Optional[float] = None,
    B_pk: Optional[float] = None,
) -> MolecularSpec:
    """Create a MolecularSpec for HF (core + 1 bond + 3 lone pairs).

    HF (hydrogen fluoride):
        Core:    F 1s² (Z=9), 2 electrons
        Bond:    F(Z_eff=7) + H(Z=1), 2 electrons
        Lone 1-3: F(Z_eff=7), 2 electrons each
        Total: 2 core + 2 bond + 6 lone = 10 electrons

    Parameters
    ----------
    max_n_core : int
        Maximum n for core orbitals on F (Z=9).
    max_n_val : int
        Maximum n for valence orbitals (screened-F and H, per group).
    R : float
        F-H bond distance in bohr (default: 1.733, expt 0.917 Å).
    E_core : float or None
        Core energy in Ha. If None, uses variational He-like estimate
        for F⁷⁺: -(Z-5/16)² = -(9-5/16)² Ha.
    include_pk : bool
        If True (default), include PK pseudopotential on F-side valence blocks.
    A_pk, B_pk : float or None
        PK parameters. If None, uses ab initio values for Z=9.
    """
    Z_F = 9
    Z_H = 1
    n_core_electrons = 2
    Z_eff_val = Z_F - n_core_electrons  # 7

    if E_core is None:
        E_core = -(Z_F - 5.0 / 16.0) ** 2

    # PK parameters
    if include_pk:
        if A_pk is None:
            A_pk = _PK_HELIKE_DEFAULTS[Z_F]['A']
        if B_pk is None:
            B_pk = _PK_HELIKE_DEFAULTS[Z_F]['B']
    else:
        A_pk = 0.0
        B_pk = 0.0

    # Nuclear repulsion
    V_FH = Z_F * Z_H / R
    V_NN = V_FH
    V_cross = _v_cross_nuc_1s(Z_F, n_core_electrons, Z_H, R)
    nuclear_repulsion = V_NN + V_cross + E_core

    blocks = [
        OrbitalBlock(
            label='F_core', block_type='core',
            Z_center=float(Z_F), n_electrons=2, max_n=max_n_core,
        ),
        OrbitalBlock(
            label='FH_bond', block_type='bond',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            has_h_partner=True, Z_partner=float(Z_H),
            max_n_partner=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
        OrbitalBlock(
            label='F_lone1', block_type='lone_pair',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
        OrbitalBlock(
            label='F_lone2', block_type='lone_pair',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
        OrbitalBlock(
            label='F_lone3', block_type='lone_pair',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
    ]
    return MolecularSpec('HF', blocks, nuclear_repulsion)


def nh3_spec(
    max_n_core: int = 2,
    max_n_val: int = 2,
    R_NH: float = 1.912,
    angle_HNH: float = 106.7,
    E_core: Optional[float] = None,
    include_pk: bool = True,
    A_pk: Optional[float] = None,
    B_pk: Optional[float] = None,
) -> MolecularSpec:
    """Create a MolecularSpec for NH3 (core + 3 bonds + 1 lone pair).

    NH3 (ammonia):
        Core:    N 1s² (Z=7), 2 electrons
        Bond 1-3: N(Z_eff=5) + H(Z=1), 2 electrons each
        Lone:    N(Z_eff=5), 2 electrons
        Total: 2 core + 6 bond + 2 lone = 10 electrons

    Parameters
    ----------
    max_n_core : int
        Maximum n for core orbitals on N (Z=7).
    max_n_val : int
        Maximum n for valence orbitals (screened-N and H, per group).
    R_NH : float
        N-H bond distance in bohr (default: 1.912, expt 1.012 Å).
    angle_HNH : float
        H-N-H bond angle in degrees (default: 106.7°, experimental).
    E_core : float or None
        Core energy in Ha. If None, uses variational He-like estimate
        for N⁵⁺: -(Z-5/16)² = -(7-5/16)² Ha.
    include_pk : bool
        If True (default), include PK pseudopotential on N-side valence blocks.
    A_pk, B_pk : float or None
        PK parameters. If None, uses ab initio values for Z=7.
    """
    Z_N = 7
    Z_H = 1
    n_core_electrons = 2
    Z_eff_val = Z_N - n_core_electrons  # 5

    if E_core is None:
        E_core = -(Z_N - 5.0 / 16.0) ** 2

    # PK parameters
    if include_pk:
        if A_pk is None:
            A_pk = _PK_HELIKE_DEFAULTS[Z_N]['A']
        if B_pk is None:
            B_pk = _PK_HELIKE_DEFAULTS[Z_N]['B']
    else:
        A_pk = 0.0
        B_pk = 0.0

    # Nuclear repulsion for pyramidal geometry
    angle_rad = np.radians(angle_HNH)
    R_HH = 2.0 * R_NH * np.sin(angle_rad / 2.0)
    V_NH = 3.0 * Z_N * Z_H / R_NH        # 3 N-H bonds
    V_HH = 3.0 * Z_H * Z_H / R_HH        # 3 H-H pairs
    V_NN = V_NH + V_HH
    V_cross = 3.0 * _v_cross_nuc_1s(Z_N, n_core_electrons, Z_H, R_NH)
    nuclear_repulsion = V_NN + V_cross + E_core

    blocks = [
        OrbitalBlock(
            label='N_core', block_type='core',
            Z_center=float(Z_N), n_electrons=2, max_n=max_n_core,
        ),
        OrbitalBlock(
            label='NH3_bond1_N', block_type='bond',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            has_h_partner=True, Z_partner=float(Z_H),
            max_n_partner=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
        OrbitalBlock(
            label='NH3_bond2_N', block_type='bond',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            has_h_partner=True, Z_partner=float(Z_H),
            max_n_partner=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
        OrbitalBlock(
            label='NH3_bond3_N', block_type='bond',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            has_h_partner=True, Z_partner=float(Z_H),
            max_n_partner=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
        OrbitalBlock(
            label='NH3_lone', block_type='lone_pair',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
    ]
    return MolecularSpec('NH3', blocks, nuclear_repulsion)


def ch4_spec(
    max_n_core: int = 2,
    max_n_val: int = 2,
    R_CH: float = 2.050,
    E_core: Optional[float] = None,
    include_pk: bool = True,
    A_pk: Optional[float] = None,
    B_pk: Optional[float] = None,
) -> MolecularSpec:
    """Create a MolecularSpec for CH4 (core + 4 bonds, no lone pairs).

    CH4 (methane):
        Core:    C 1s² (Z=6), 2 electrons
        Bond 1-4: C(Z_eff=4) + H(Z=1), 2 electrons each
        Total: 2 core + 8 bond = 10 electrons

    Parameters
    ----------
    max_n_core : int
        Maximum n for core orbitals on C (Z=6).
    max_n_val : int
        Maximum n for valence orbitals (screened-C and H, per group).
    R_CH : float
        C-H bond distance in bohr (default: 2.050, expt 1.085 Å).
    E_core : float or None
        Core energy in Ha. If None, uses variational He-like estimate
        for C⁴⁺: -(Z-5/16)² = -(6-5/16)² Ha.
    include_pk : bool
        If True (default), include PK pseudopotential on C-side valence blocks.
    A_pk, B_pk : float or None
        PK parameters. If None, uses ab initio values for Z=6.
    """
    Z_C = 6
    Z_H = 1
    n_core_electrons = 2
    Z_eff_val = Z_C - n_core_electrons  # 4

    if E_core is None:
        E_core = -(Z_C - 5.0 / 16.0) ** 2

    # PK parameters
    if include_pk:
        if A_pk is None:
            A_pk = _PK_HELIKE_DEFAULTS[Z_C]['A']
        if B_pk is None:
            B_pk = _PK_HELIKE_DEFAULTS[Z_C]['B']
    else:
        A_pk = 0.0
        B_pk = 0.0

    # Nuclear repulsion for tetrahedral geometry
    # In a perfect tetrahedron, H-H distance = R_CH * sqrt(8/3)
    R_HH = R_CH * np.sqrt(8.0 / 3.0)
    V_CH = 4.0 * Z_C * Z_H / R_CH        # 4 C-H bonds
    V_HH = 6.0 * Z_H * Z_H / R_HH        # 6 H-H pairs
    V_NN = V_CH + V_HH
    V_cross = 4.0 * _v_cross_nuc_1s(Z_C, n_core_electrons, Z_H, R_CH)
    nuclear_repulsion = V_NN + V_cross + E_core

    blocks = [
        OrbitalBlock(
            label='C_core', block_type='core',
            Z_center=float(Z_C), n_electrons=2, max_n=max_n_core,
        ),
        OrbitalBlock(
            label='CH4_bond1_C', block_type='bond',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            has_h_partner=True, Z_partner=float(Z_H),
            max_n_partner=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
        OrbitalBlock(
            label='CH4_bond2_C', block_type='bond',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            has_h_partner=True, Z_partner=float(Z_H),
            max_n_partner=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
        OrbitalBlock(
            label='CH4_bond3_C', block_type='bond',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            has_h_partner=True, Z_partner=float(Z_H),
            max_n_partner=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
        OrbitalBlock(
            label='CH4_bond4_C', block_type='bond',
            Z_center=float(Z_eff_val), n_electrons=2,
            max_n=max_n_val,
            has_h_partner=True, Z_partner=float(Z_H),
            max_n_partner=max_n_val,
            pk_A=A_pk, pk_B=B_pk,
        ),
    ]
    return MolecularSpec('CH4', blocks, nuclear_repulsion)


# ---------------------------------------------------------------------------
# Main builder (original, now delegates to general builder)
# ---------------------------------------------------------------------------

def build_composed_lih(
    max_n_core: int = 2,
    max_n_val: int = 2,
    R: float = 3.015,
    E_core: Optional[float] = None,
    include_pk: bool = True,
    pk_in_hamiltonian: Optional[bool] = None,
    A_pk: Optional[float] = None,
    B_pk: Optional[float] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Build the composed LiH qubit Hamiltonian and count Pauli terms.

    Orbital basis:
        Core block:    hydrogenic (n,l,m) on Li center, Z_A=3, up to max_n_core
        Valence block: hydrogenic on screened Li (Z_eff~1) up to max_n_val
                       + hydrogenic on H center (Z_B=1) up to max_n_val

    The pseudopotential approximation (Paper 17) means:
        - Cross-block ERIs are zero (core-valence Coulomb -> one-body Z_eff)
        - Cross-block h1 terms (PK, cross-center nuclear attraction) are set
          to zero in this first pass. See APPROXIMATIONS section below.

    APPROXIMATIONS (documented per CLAUDE.md theory check rule):
        1. PK pseudopotential on valence-Li orbitals is included when
           include_pk=True (default). PK on valence-H orbitals requires
           two-center integration (Li-centered potential on H-centered
           orbitals) and is not included -- see cross-center note below.
        2. Cross-block ERI = 0. This is physically justified by the
           pseudopotential approximation (Paper 17 Sec II).
        3. Valence Z_eff is approximated as Z_A - n_core = 1 (complete
           screening). The true Z_eff(r) is r-dependent; using a constant
           is exact at large r but approximate near the nucleus.
        4. E_core default uses Li²⁺ He-like value -7.2799 Ha (from
           hyperspherical solver, Paper 13).
        5. Cross-center nuclear attraction integrals (e.g., H orbital
           feeling Li nucleus) require two-center radial integration
           in prolate spheroidal coordinates. These are NOT included.
           Estimated impact: up to M_val_li * M_val_h additional nonzero
           h1 entries per center, but coupling is exponentially small at
           large R (screening length ~ 1/Z).

    Parameters
    ----------
    max_n_core : int
        Maximum n for core orbitals on Li (Z=3).
    max_n_val : int
        Maximum n for valence orbitals (both screened-Li and H).
    R : float
        Internuclear distance in bohr (default: 3.015, expt R_eq for LiH).
    E_core : float or None
        Core energy in Ha. If None, uses -7.2799 Ha (Li²⁺ ground state
        from hyperspherical solver at l_max=2).
    include_pk : bool
        If True (default), compute PK pseudopotential matrix elements on
        valence-Li orbitals and add to h1. Paper 17 Sec IV.
    A_pk : float or None
        PK barrier height (Ha·bohr²). If None, uses Paper 17 default for Z=3.
    B_pk : float or None
        PK barrier width exponent (bohr⁻²). If None, uses Paper 17 default.
    verbose : bool
        Print progress and results.

    Returns
    -------
    dict with keys:
        M, Q, N_pauli, M_core, M_val,
        ERI_density_total, ERI_density_core, ERI_density_val,
        h1, eri, nuclear_repulsion, qubit_op, fermion_op,
        states_core, states_val
    """
    t0 = time.perf_counter()

    # --- Physical parameters ---
    Z_A = 3      # Li nuclear charge
    Z_B = 1      # H nuclear charge
    n_core_electrons = 2
    Z_eff_val = Z_A - n_core_electrons  # = 1 (complete screening approx)

    if E_core is None:
        # Li²⁺ (He-like, Z=3) ground state from hyperspherical solver
        # Paper 13: He at Z=2 gives -2.9037 Ha; Z=3 gives -7.2799 Ha
        E_core = -7.2799

    # --- Define orbital basis ---
    states_core = _enumerate_states(max_n_core)
    M_core = len(states_core)

    # Valence: screened Li orbitals + H orbitals
    states_val_li = _enumerate_states(max_n_val)   # Z_eff orbitals on Li
    states_val_h = _enumerate_states(max_n_val)    # Z=1 orbitals on H
    M_val_li = len(states_val_li)
    M_val_h = len(states_val_h)
    M_val = M_val_li + M_val_h

    M = M_core + M_val
    Q = 2 * M  # spin-orbitals = qubits under JW

    if verbose:
        print(f"[composed_qubit] Orbital basis:")
        print(f"  Core (Z={Z_A}):         {M_core} orbitals, max_n={max_n_core}")
        print(f"  Val Li (Z_eff={Z_eff_val}): {M_val_li} orbitals, max_n={max_n_val}")
        print(f"  Val H  (Z={Z_B}):       {M_val_h} orbitals, max_n={max_n_val}")
        print(f"  Total: M={M} spatial, Q={Q} qubits")

    # --- Construct h1 (M x M) ---
    # Block-diagonal: exact hydrogenic eigenvalues on diagonal
    h1 = np.zeros((M, M))

    # Core block: -Z_A^2 / (2n^2) for each (n,l,m)
    for i, (n, l, m) in enumerate(states_core):
        h1[i, i] = -Z_A**2 / (2.0 * n**2)

    # Valence-Li block: -Z_eff^2 / (2n^2)
    offset_val_li = M_core
    for i, (n, l, m) in enumerate(states_val_li):
        idx = offset_val_li + i
        h1[idx, idx] = -Z_eff_val**2 / (2.0 * n**2)

    # Valence-H block: -Z_B^2 / (2n^2)
    offset_val_h = M_core + M_val_li
    for i, (n, l, m) in enumerate(states_val_h):
        idx = offset_val_h + i
        h1[idx, idx] = -Z_B**2 / (2.0 * n**2)

    # --- PK pseudopotential on valence-Li orbitals (Paper 17 Sec IV) ---
    # pk_in_hamiltonian controls whether PK is added to the quantum Hamiltonian.
    # When False (default when include_pk=True), PK is computed and returned
    # separately for classical post-processing via 1-RDM.
    if pk_in_hamiltonian is None:
        pk_in_hamiltonian = include_pk  # backward compat: old behavior

    # Full-size PK matrix (always computed when include_pk=True)
    h1_pk_full = np.zeros((M, M))
    n_pk_nonzero = 0
    if include_pk:
        # Resolve PK parameters
        if A_pk is None or B_pk is None:
            pk_defaults = _PK_DEFAULTS.get(Z_A)
            if pk_defaults is None:
                raise ValueError(
                    f"No default PK parameters for Z={Z_A}. "
                    f"Provide A_pk and B_pk explicitly."
                )
            if A_pk is None:
                A_pk = pk_defaults['A']
            if B_pk is None:
                B_pk = pk_defaults['B']

        if verbose:
            print(f"[composed_qubit] PK pseudopotential: A={A_pk:.2f}, B={B_pk:.2f}")
            print(f"  Source: {_PK_DEFAULTS.get(Z_A, {}).get('source', 'user-supplied')}")

        # PK on valence-Li orbitals (same-center, Z_eff)
        h1_pk_val_li = _compute_pk_matrix_elements(
            Z_eff_val, states_val_li, A_pk, B_pk,
        )
        # Store in full-size PK matrix
        for i in range(M_val_li):
            for j in range(M_val_li):
                if abs(h1_pk_val_li[i, j]) > 1e-15:
                    h1_pk_full[offset_val_li + i, offset_val_li + j] = h1_pk_val_li[i, j]
                    n_pk_nonzero += 1

        # Only add PK to h1 if requested
        if pk_in_hamiltonian:
            h1 += h1_pk_full

        if verbose:
            diag_pk = np.diag(h1_pk_val_li)
            print(f"  Val-Li PK nonzero entries: {n_pk_nonzero}")
            print(f"  PK diagonal range: [{diag_pk.min():.6f}, {diag_pk.max():.6f}] Ha")
            print(f"  PK in quantum Hamiltonian: {pk_in_hamiltonian}")

        # PK on valence-H orbitals: SKIPPED (two-center integration required)
        # The PK potential is centered on Li. For H-centered orbitals φ_H(r_H),
        # the integral <φ_H | V_PK(r_Li) | φ_H> requires expressing r_Li in
        # terms of r_H (r_Li² = r_H² + R² - 2*r_H*R*cos θ), leading to a
        # non-separable 3D integral that requires prolate spheroidal coordinates
        # or numerical quadrature. At R=3.015 bohr with B=7.0, the PK barrier
        # has decayed by exp(-7*3²) ~ 10⁻²⁸ at the H center, so the cross-center
        # PK contribution is negligible.

    # Cross-center nuclear attraction: SKIPPED (two-center integrals)
    # <φ_p^Li | -Z_B/r_B | φ_q^Li> and <φ_p^H | -Z_A_eff/r_A | φ_q^H>
    # require two-center integration. These would add at most:
    #   M_val_li² + M_val_h² additional entries in h1
    # For max_n_val=2: 5² + 5² = 50 entries (25 unique pairs per block).
    # Implementation requires prolate spheroidal auxiliary integrals --
    # deferred to future work. The impact is bounded by the cross-center
    # overlap which decays exponentially with R.
    n_cross_nuc_estimate = M_val_li**2 + M_val_h**2

    if pk_in_hamiltonian:
        h1_desc = "diagonal + PK"
    elif include_pk:
        h1_desc = "diagonal only (PK computed separately)"
    else:
        h1_desc = "diagonal only"
    if verbose:
        print(f"[composed_qubit] h1 constructed ({M}x{M}), {h1_desc}")
        if not include_pk:
            print(f"  NOTE: PK disabled (include_pk=False)")
        print(f"  Cross-center nuclear attraction: SKIPPED (two-center)")
        print(f"    Estimated missing entries: {n_cross_nuc_estimate}")

    # --- Construct ERI (M x M x M x M) in chemist notation ---
    # Block-diagonal by pseudopotential approximation
    eri = np.zeros((M, M, M, M))

    # Core-core block: Gaunt-integral ERIs at Z=3
    if verbose:
        print(f"[composed_qubit] Computing core R^k integrals (Z={Z_A})...")
    rk_core = _compute_rk_integrals_block(Z_A, states_core)
    eri_core_phys = _build_eri_block(Z_A, states_core, rk_core)
    n_eri_core = len(eri_core_phys)

    # Convert core block to chemist notation and place in full tensor
    for (a, b, c, d), val in eri_core_phys.items():
        # phys <ab|cd> -> chemist (ac|bd)
        eri[a, c, b, d] = val

    # Valence-Li block: Gaunt-integral ERIs at Z_eff
    if verbose:
        print(f"[composed_qubit] Computing val-Li R^k integrals (Z_eff={Z_eff_val})...")
    rk_val_li = _compute_rk_integrals_block(Z_eff_val, states_val_li)
    eri_val_li_phys = _build_eri_block(Z_eff_val, states_val_li, rk_val_li)
    n_eri_val_li = len(eri_val_li_phys)

    for (a, b, c, d), val in eri_val_li_phys.items():
        p, q, r, s = (a + offset_val_li, c + offset_val_li,
                       b + offset_val_li, d + offset_val_li)
        eri[p, q, r, s] = val

    # Valence-H block: Gaunt-integral ERIs at Z_B=1
    if verbose:
        print(f"[composed_qubit] Computing val-H R^k integrals (Z={Z_B})...")
    rk_val_h = _compute_rk_integrals_block(Z_B, states_val_h)
    eri_val_h_phys = _build_eri_block(Z_B, states_val_h, rk_val_h)
    n_eri_val_h = len(eri_val_h_phys)

    for (a, b, c, d), val in eri_val_h_phys.items():
        p, q, r, s = (a + offset_val_h, c + offset_val_h,
                       b + offset_val_h, d + offset_val_h)
        eri[p, q, r, s] = val

    # Cross blocks: zero by pseudopotential approximation (Paper 17 Sec II)
    # The core-valence electron repulsion is replaced by one-body Z_eff screening.

    # Symmetrize: enforce (pq|rs) = (rs|pq) exactly
    # (numerical integration introduces ~1e-4 asymmetry)
    eri = 0.5 * (eri + eri.transpose(2, 3, 0, 1))

    # --- ERI density statistics ---
    n_eri_total = int(np.count_nonzero(np.abs(eri) > 1e-15))
    eri_density_total = n_eri_total / max(1, M**4)
    eri_density_core = n_eri_core / max(1, M_core**4)
    # Valence density uses combined valence orbital count
    n_eri_val_total = n_eri_val_li + n_eri_val_h
    eri_density_val = n_eri_val_total / max(1, M_val**4)

    if verbose:
        print(f"[composed_qubit] ERI statistics:")
        print(f"  Core:  {n_eri_core} nonzero / {M_core**4} = {eri_density_core:.1%}")
        print(f"  Val:   {n_eri_val_total} nonzero / {M_val**4} = {eri_density_val:.1%}")
        print(f"  Total: {n_eri_total} nonzero / {M**4} = {eri_density_total:.1%}")

    # --- Nuclear repulsion ---
    # V_NN = Z_A * Z_B / R  (bare nuclear repulsion)
    V_NN = Z_A * Z_B / R

    # V_cross_nuc: core-to-H attraction (Paper 17 Eq. 4)
    V_cross = _v_cross_nuc_1s(Z_A, n_core_electrons, Z_B, R)

    # Total nuclear repulsion constant includes core energy
    nuclear_repulsion = V_NN + V_cross + E_core

    if verbose:
        print(f"[composed_qubit] Energy constants:")
        print(f"  V_NN        = {V_NN:.6f} Ha")
        print(f"  V_cross_nuc = {V_cross:.6f} Ha")
        print(f"  E_core      = {E_core:.6f} Ha")
        print(f"  nuc_repul   = {nuclear_repulsion:.6f} Ha")

    # --- JW transform and Pauli count ---
    if verbose:
        print(f"[composed_qubit] Building fermion operator...")
    fermion_op = build_fermion_op_from_integrals(h1, eri, nuclear_repulsion)

    if verbose:
        print(f"[composed_qubit] Jordan-Wigner transform...")
    qubit_op = jordan_wigner(fermion_op)
    N_pauli = len(qubit_op.terms)

    elapsed = time.perf_counter() - t0

    if verbose:
        print(f"\n{'='*60}")
        print(f"COMPOSED LiH QUBIT HAMILTONIAN -- RESULTS")
        print(f"{'='*60}")
        print(f"  Spatial orbitals M = {M}  (core={M_core}, val={M_val})")
        print(f"  Qubits           Q = {Q}")
        print(f"  Pauli terms      N = {N_pauli:,}")
        print(f"  ERI density (total)   = {eri_density_total:.2%}")
        print(f"  ERI density (core)    = {eri_density_core:.2%}")
        print(f"  ERI density (valence) = {eri_density_val:.2%}")
        print(f"  PK included          = {include_pk}"
              f"  (nonzero entries: {n_pk_nonzero})")
        print(f"  Cross-nuc skipped    = {n_cross_nuc_estimate} est. entries")
        print(f"  Gaussian LiH STO-3G   = 631 Pauli terms at 12 qubits")
        if Q == 12:
            ratio = 631 / max(1, N_pauli)
            print(f"  GeoVac advantage ratio = {ratio:.2f}x at Q=12")
        print(f"  Wall time = {elapsed:.1f}s")
        print(f"{'='*60}")

    # --- Save results ---
    results = {
        'M': M,
        'Q': Q,
        'N_pauli': N_pauli,
        'M_core': M_core,
        'M_val': M_val,
        'M_val_li': M_val_li,
        'M_val_h': M_val_h,
        'ERI_density_total': eri_density_total,
        'ERI_density_core': eri_density_core,
        'ERI_density_val': eri_density_val,
        'n_eri_core': n_eri_core,
        'n_eri_val_li': n_eri_val_li,
        'n_eri_val_h': n_eri_val_h,
        'n_eri_total': n_eri_total,
        'nuclear_repulsion': nuclear_repulsion,
        'V_NN': V_NN,
        'V_cross_nuc': V_cross,
        'E_core': E_core,
        'R_bohr': R,
        'max_n_core': max_n_core,
        'max_n_val': max_n_val,
        'Z_A': Z_A,
        'Z_B': Z_B,
        'Z_eff_val': Z_eff_val,
        'gaussian_lih_sto3g_pauli': 631,
        'gaussian_lih_sto3g_qubits': 12,
        'wall_time_s': elapsed,
        'include_pk': include_pk,
        'pk_in_hamiltonian': pk_in_hamiltonian,
        'A_pk': A_pk if include_pk else None,
        'B_pk': B_pk if include_pk else None,
        'n_pk_nonzero': n_pk_nonzero,
        'n_cross_nuc_estimate': n_cross_nuc_estimate,
        'approximations': [
            ('PK in quantum Hamiltonian' if pk_in_hamiltonian
             else 'PK computed separately (classical correction)')
            if include_pk else 'PK disabled',
            'PK on val-H skipped (two-center, negligible at R=3 bohr)',
            'Cross-center nuclear attraction skipped (two-center)',
            'Cross-block ERI = 0 (pseudopotential approximation)',
            'Valence Z_eff = Z_A - n_core (constant, not r-dependent)',
            'E_core from hyperspherical solver default (-7.2799 Ha)',
        ],
    }

    # Save JSON report
    out_path = Path(__file__).parent.parent / 'debug' / 'data' / 'composed_lih_pauli_analysis.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    if verbose:
        print(f"\nResults saved to {out_path}")

    # Attach non-serializable objects for programmatic use
    results['h1'] = h1
    results['h1_pk'] = h1_pk_full
    results['eri'] = eri
    results['qubit_op'] = qubit_op
    results['fermion_op'] = fermion_op
    results['states_core'] = states_core
    results['states_val_li'] = states_val_li
    results['states_val_h'] = states_val_h

    return results


# ---------------------------------------------------------------------------
# Scaling sweep (Part 1)
# ---------------------------------------------------------------------------

# He single-geometry reference data (Paper 14, validated)
_HE_SCALING_DATA = {
    2: {'Q': 10, 'N_pauli': 120},
    3: {'Q': 28, 'N_pauli': 2_659},
    4: {'Q': 60, 'N_pauli': 31_039},
    5: {'Q': 110, 'N_pauli': 227_338},
}


def composed_lih_scaling_sweep(
    max_n_values: Optional[List[int]] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Run composed LiH at multiple basis sizes and fit Pauli scaling exponent.

    Calls build_composed_lih(max_n_core=n, max_n_val=n) for each n,
    collects Pauli term counts, fits N_pauli = a x Q^α, and compares
    to He single-geometry (Q^3.15) and Gaussian molecular (Q^4.60) scaling.

    Parameters
    ----------
    max_n_values : list of int or None
        Basis sizes to sweep. Default [1, 2, 3].
    verbose : bool
        Print progress and comparison table.

    Returns
    -------
    dict with sweep data, fitted exponent, and comparison table.
    """
    if max_n_values is None:
        max_n_values = [1, 2, 3]

    sweep_data: List[Dict[str, Any]] = []

    for n in max_n_values:
        if verbose:
            print(f"\n{'='*60}")
            print(f"  SCALING SWEEP: max_n = {n}")
            print(f"{'='*60}")

        result = build_composed_lih(
            max_n_core=n, max_n_val=n, verbose=verbose,
        )

        # Count h1-only vs ERI Pauli contributions
        # Build h1-only version to measure h1 contribution
        fermion_h1_only = build_fermion_op_from_integrals(
            result['h1'], np.zeros_like(result['eri']),
            result['nuclear_repulsion'],
        )
        qubit_h1_only = jordan_wigner(fermion_h1_only)
        n_pauli_h1 = len(qubit_h1_only.terms)
        n_pauli_eri = result['N_pauli'] - n_pauli_h1

        entry = {
            'max_n': n,
            'M': result['M'],
            'Q': result['Q'],
            'N_pauli': result['N_pauli'],
            'N_pauli_h1': n_pauli_h1,
            'N_pauli_eri': n_pauli_eri,
            'ERI_density_total': result['ERI_density_total'],
            'n_pk_nonzero': result['n_pk_nonzero'],
            'wall_time_s': result['wall_time_s'],
        }
        sweep_data.append(entry)

    # --- Fit power law: N_pauli = a * Q^alpha ---
    Q_arr = np.array([d['Q'] for d in sweep_data], dtype=float)
    N_arr = np.array([d['N_pauli'] for d in sweep_data], dtype=float)

    alpha_fit = np.nan
    r_squared = np.nan
    a_fit = np.nan

    if len(Q_arr) >= 2 and np.all(Q_arr > 0) and np.all(N_arr > 0):
        log_Q = np.log(Q_arr)
        log_N = np.log(N_arr)
        # Linear regression in log-log space
        coeffs = np.polyfit(log_Q, log_N, 1)
        alpha_fit = float(coeffs[0])
        a_fit = float(np.exp(coeffs[1]))

        # R² calculation
        log_N_pred = np.polyval(coeffs, log_Q)
        ss_res = np.sum((log_N - log_N_pred)**2)
        ss_tot = np.sum((log_N - np.mean(log_N))**2)
        r_squared = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 1.0

    # --- Comparison table ---
    if verbose:
        print(f"\n{'='*70}")
        print(f"  COMPOSED LiH SCALING SWEEP -- RESULTS")
        print(f"{'='*70}")
        print(f"  {'max_n':>5} {'M':>4} {'Q':>4} {'N_pauli':>10}"
              f" {'N_h1':>8} {'N_eri':>8} {'ERI_dens':>10}")
        print(f"  {'-'*5} {'-'*4} {'-'*4} {'-'*10}"
              f" {'-'*8} {'-'*8} {'-'*10}")
        for d in sweep_data:
            print(f"  {d['max_n']:>5} {d['M']:>4} {d['Q']:>4}"
                  f" {d['N_pauli']:>10,}"
                  f" {d['N_pauli_h1']:>8,} {d['N_pauli_eri']:>8,}"
                  f" {d['ERI_density_total']:>10.2%}")

        print(f"\n  Fitted exponent: N_pauli = {a_fit:.2f} x Q^{alpha_fit:.2f}"
              f"  (R² = {r_squared:.4f})")

        print(f"\n  --- Scaling comparison ---")
        print(f"  {'System':>25} {'Exponent':>10} {'Notes':>30}")
        print(f"  {'-'*25} {'-'*10} {'-'*30}")
        print(f"  {'GeoVac composed LiH':>25} {alpha_fit:>10.2f}"
              f" {'this sweep':>30}")
        print(f"  {'GeoVac He (single-geom)':>25} {'3.15':>10}"
              f" {'Paper 14, nmax=2-5':>30}")
        print(f"  {'Gaussian molecular':>25} {'4.60':>10}"
              f" {'literature estimate':>30}")

        # He comparison at matching max_n
        print(f"\n  --- He single-geometry comparison at matching max_n ---")
        print(f"  {'max_n':>5} {'Q_comp':>7} {'N_comp':>10}"
              f" {'Q_He':>6} {'N_He':>10} {'ratio':>8}")
        for d in sweep_data:
            n = d['max_n']
            he = _HE_SCALING_DATA.get(n)
            if he:
                ratio = d['N_pauli'] / he['N_pauli'] if he['N_pauli'] > 0 else float('nan')
                print(f"  {n:>5} {d['Q']:>7} {d['N_pauli']:>10,}"
                      f" {he['Q']:>6} {he['N_pauli']:>10,}"
                      f" {ratio:>8.2f}x")
            else:
                print(f"  {n:>5} {d['Q']:>7} {d['N_pauli']:>10,}"
                      f" {'N/A':>6} {'N/A':>10} {'N/A':>8}")

        print(f"  Gaussian LiH STO-3G: Q=12, N=631 Pauli terms")
        print(f"{'='*70}")

    # --- Save results ---
    output = {
        'sweep_data': sweep_data,
        'fit': {
            'alpha': alpha_fit,
            'a': a_fit,
            'R_squared': r_squared,
            'formula': f'N_pauli = {a_fit:.2f} * Q^{alpha_fit:.2f}',
        },
        'he_single_geometry': {
            str(k): v for k, v in _HE_SCALING_DATA.items()
        },
        'he_exponent': 3.15,
        'gaussian_molecular_exponent': 4.60,
        'gaussian_lih_sto3g': {'Q': 12, 'N_pauli': 631},
    }

    out_path = (Path(__file__).parent.parent
                / 'debug' / 'data' / 'composed_lih_scaling_sweep.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    if verbose:
        print(f"\nResults saved to {out_path}")

    return output


# ---------------------------------------------------------------------------
# Full analysis with comparison table (Part 3)
# ---------------------------------------------------------------------------

def composed_lih_full_analysis(
    max_n_values: Optional[List[int]] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Run scaling sweep and produce comprehensive comparison analysis.

    Generates debug/data/composed_lih_scaling_analysis.json with:
    - Scaling sweep data (Q, N_pauli, ERI_density) at each max_n
    - Fitted power law exponent α with R²
    - Per-max_n breakdown: Pauli terms from h1 vs from ERI
    - Comparison: Gaussian LiH STO-3G (Q=12, 631 terms)
    - Comparison: He single-geometry at each matching max_n
    """
    if max_n_values is None:
        max_n_values = [1, 2, 3]

    # Run the sweep (includes PK by default)
    sweep_result = composed_lih_scaling_sweep(max_n_values, verbose=verbose)

    # Also run without PK at default basis to measure PK impact
    result_no_pk = build_composed_lih(
        max_n_core=2, max_n_val=2, include_pk=False, verbose=False,
    )
    result_with_pk = build_composed_lih(
        max_n_core=2, max_n_val=2, include_pk=True, verbose=False,
    )

    pk_impact = {
        'N_pauli_without_pk': result_no_pk['N_pauli'],
        'N_pauli_with_pk': result_with_pk['N_pauli'],
        'delta_pauli': result_with_pk['N_pauli'] - result_no_pk['N_pauli'],
        'n_pk_nonzero_h1_entries': result_with_pk['n_pk_nonzero'],
        'max_n_core': 2,
        'max_n_val': 2,
    }

    # Build He comparison rows
    he_comparison = []
    for d in sweep_result['sweep_data']:
        n = d['max_n']
        he = _HE_SCALING_DATA.get(n)
        row = {
            'max_n': n,
            'composed_Q': d['Q'],
            'composed_N_pauli': d['N_pauli'],
        }
        if he:
            row['he_Q'] = he['Q']
            row['he_N_pauli'] = he['N_pauli']
            row['ratio_composed_over_he'] = (
                d['N_pauli'] / he['N_pauli'] if he['N_pauli'] > 0 else None
            )
        he_comparison.append(row)

    analysis = {
        'scaling_sweep': sweep_result['sweep_data'],
        'power_law_fit': sweep_result['fit'],
        'pk_impact': pk_impact,
        'he_comparison': he_comparison,
        'gaussian_lih_sto3g': {'Q': 12, 'N_pauli': 631},
        'reference_exponents': {
            'geovac_he_single_geometry': 3.15,
            'gaussian_molecular': 4.60,
        },
        'cross_center_nuclear_attraction': {
            'status': 'skipped',
            'reason': 'Two-center integrals require prolate spheroidal auxiliary functions',
            'estimated_missing_entries_max_n2': 50,
        },
    }

    out_path = (Path(__file__).parent.parent
                / 'debug' / 'data' / 'composed_lih_scaling_analysis.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(analysis, f, indent=2)
    if verbose:
        print(f"\nFull analysis saved to {out_path}")

    return analysis


# ---------------------------------------------------------------------------
# Cross-center valence ERI bound (Deliverable 2)
# ---------------------------------------------------------------------------

def estimate_cross_center_eri_count(
    max_n: int,
) -> Dict[str, Any]:
    """
    Estimate the number of nonzero cross-center valence-valence ERIs
    that would exist if two-center integrals were included.

    For orbitals on center A with (l1,m1), (l3,m3) and center B with
    (l2,m2), (l4,m4), the Gaunt selection rules still apply:
      - Triangle inequality: |l1-l3| <= k <= l1+l3
      - Parity: l1+k+l3 = even (and l2+k+l4 = even)
      - m-conservation: m1+m2 = m3+m4

    This counts qualifying quadruplets (a on Li-val, b on H-val, c on
    Li-val, d on H-val) without computing actual two-center integrals.

    Parameters
    ----------
    max_n : int
        Maximum principal quantum number for valence orbitals.

    Returns
    -------
    dict with:
        n_cross_center : int
            Number of cross-center ERI quadruplets passing selection rules.
        n_same_center : int
            Number of same-center ERI entries (sum of both valence blocks).
        ratio : float
            n_cross_center / n_same_center.
        M_val : int
            Number of valence orbitals per center.
        upper_bound : int
            M_val^4 (absolute maximum without selection rules).
    """
    states = _enumerate_states(max_n)
    M_val = len(states)

    # Count same-center nonzero ERIs using Gaunt rules (one center)
    n_same_center = 0
    for a in range(M_val):
        la, ma = states[a][1], states[a][2]
        for c in range(M_val):
            lc, mc = states[c][1], states[c][2]
            # Find valid k values for (a,c) pair
            k_values_ac = []
            for k in range(abs(la - lc), la + lc + 1):
                if (la + lc + k) % 2 != 0:
                    continue
                # Check 3j m-selection: -ma + q + mc = 0 where q = mc - ma
                # This is always satisfied by definition of q
                # But need |q| <= k
                q_ac = mc - ma
                if abs(q_ac) <= k:
                    k_values_ac.append(k)

            if not k_values_ac:
                continue

            for b in range(M_val):
                lb, mb = states[b][1], states[b][2]
                for d in range(M_val):
                    ld, md = states[d][1], states[d][2]
                    # m-conservation: ma + mb = mc + md
                    if ma + mb != mc + md:
                        continue
                    # Check (b,d) Gaunt selection for at least one matching k
                    for k in k_values_ac:
                        if abs(lb - ld) <= k <= lb + ld:
                            if (lb + ld + k) % 2 == 0:
                                q_bd = md - mb
                                if abs(q_bd) <= k:
                                    n_same_center += 1
                                    break

    # Count cross-center ERIs: a,c on Li-val; b,d on H-val
    # Same Gaunt selection rules apply (angular parts factorize)
    n_cross_center = 0
    for a in range(M_val):
        la, ma = states[a][1], states[a][2]
        for c in range(M_val):
            lc, mc = states[c][1], states[c][2]
            k_values_ac = []
            for k in range(abs(la - lc), la + lc + 1):
                if (la + lc + k) % 2 != 0:
                    continue
                q_ac = mc - ma
                if abs(q_ac) <= k:
                    k_values_ac.append(k)

            if not k_values_ac:
                continue

            for b in range(M_val):
                lb, mb = states[b][1], states[b][2]
                for d in range(M_val):
                    ld, md = states[d][1], states[d][2]
                    if ma + mb != mc + md:
                        continue
                    for k in k_values_ac:
                        if abs(lb - ld) <= k <= lb + ld:
                            if (lb + ld + k) % 2 == 0:
                                q_bd = md - mb
                                if abs(q_bd) <= k:
                                    n_cross_center += 1
                                    break

    upper_bound = M_val ** 4
    ratio = n_cross_center / max(1, n_same_center)

    return {
        'max_n': max_n,
        'n_cross_center': n_cross_center,
        'n_same_center': n_same_center,
        'ratio': ratio,
        'M_val': M_val,
        'upper_bound': upper_bound,
        'cross_fraction_of_upper': n_cross_center / max(1, upper_bound),
        'same_fraction_of_upper': n_same_center / max(1, upper_bound),
    }


# ---------------------------------------------------------------------------
# Gaussian LiH comparison (Deliverable 3)
# ---------------------------------------------------------------------------

# Published / estimated Gaussian LiH Pauli term counts.
# STO-3G (6 spatial, Q=12): 631 Pauli terms from literature.
# Larger bases estimated using Gaussian molecular scaling Q^4.60.
_GAUSSIAN_LIH_DATA = {
    'sto-3g': {
        'M': 6, 'Q': 12, 'N_pauli': 631,
        'source': 'published (JW of STO-3G LiH)',
        'computed': False,
    },
    '6-31g': {
        'M': 11, 'Q': 22, 'N_pauli': None,  # filled by estimation
        'source': 'estimated from Q^4.60 scaling',
        'computed': False,
    },
    'cc-pvdz': {
        'M': 28, 'Q': 56, 'N_pauli': None,
        'source': 'estimated from Q^4.60 scaling',
        'computed': False,
    },
}

# Gaussian molecular scaling fit from Paper 14
_GAUSS_MOL_ALPHA = 4.60
_GAUSS_MOL_C = None  # calibrated from STO-3G anchor


def _estimate_gaussian_lih_pauli(Q: int) -> int:
    """Estimate Gaussian LiH Pauli count from STO-3G anchor + Q^4.60 scaling."""
    # Calibrate: 631 = C * 12^4.60  =>  C = 631 / 12^4.60
    Q_anchor = 12
    N_anchor = 631
    C = N_anchor / (Q_anchor ** _GAUSS_MOL_ALPHA)
    return int(round(C * Q ** _GAUSS_MOL_ALPHA))


def build_gaussian_lih_reference() -> Dict[str, Any]:
    """
    Build Gaussian LiH reference data for comparison.

    PySCF is not available on this platform, so we use the published
    STO-3G value (631 Pauli terms at Q=12) and estimate larger bases
    using the Gaussian molecular scaling Q^4.60 from Paper 14.

    Returns
    -------
    dict with basis-keyed data and equal-qubit comparison table.
    """
    data = {}
    for basis, info in _GAUSSIAN_LIH_DATA.items():
        entry = dict(info)
        if entry['N_pauli'] is None:
            entry['N_pauli'] = _estimate_gaussian_lih_pauli(entry['Q'])
        data[basis] = entry

    return {
        'gaussian_data': data,
        'scaling_exponent': _GAUSS_MOL_ALPHA,
        'anchor': {'basis': 'sto-3g', 'Q': 12, 'N_pauli': 631},
        'note': (
            'PySCF unavailable on Windows. STO-3G value is published; '
            '6-31G and cc-pVDZ values estimated from Q^4.60 molecular '
            'Gaussian scaling anchored at STO-3G.'
        ),
    }


def build_equal_qubit_comparison() -> List[Dict[str, Any]]:
    """
    For each GeoVac composed data point, find the nearest Gaussian basis
    and compare Pauli counts.

    Returns list of comparison rows.
    """
    geovac_points = [
        {'max_n': 1, 'Q': 6, 'N_pauli': 10},
        {'max_n': 2, 'Q': 30, 'N_pauli': 334},
        {'max_n': 3, 'Q': 84, 'N_pauli': 7879},
        {'max_n': 4, 'Q': 180, 'N_pauli': 92899},
    ]

    rows = []
    for gv in geovac_points:
        Q_gv = gv['Q']
        N_gauss = _estimate_gaussian_lih_pauli(Q_gv)
        ratio = N_gauss / max(1, gv['N_pauli'])
        rows.append({
            'max_n': gv['max_n'],
            'Q': Q_gv,
            'N_pauli_geovac': gv['N_pauli'],
            'N_pauli_gaussian_est': N_gauss,
            'gaussian_advantage_ratio': round(ratio, 1),
            'note': 'Gaussian estimated from Q^4.60 anchored at STO-3G',
        })

    return rows


# ---------------------------------------------------------------------------
# Published Gaussian LiH power-law fit (Trenev et al.)
# ---------------------------------------------------------------------------

def fit_gaussian_lih_published_exponent() -> Dict[str, Any]:
    """
    Fit a power law N_pauli = C * Q^alpha to the 3 published LiH data points
    from Trenev et al. (Quantum 2025, Table 5, JW with 2-qubit reduction).

    Returns
    -------
    dict with alpha, C, R_squared, and the fitted formula.
    """
    Q_arr = np.array([v['Q'] for v in GAUSSIAN_LIH_PUBLISHED.values()], dtype=float)
    N_arr = np.array([v['N_pauli'] for v in GAUSSIAN_LIH_PUBLISHED.values()], dtype=float)

    log_Q = np.log(Q_arr)
    log_N = np.log(N_arr)

    coeffs = np.polyfit(log_Q, log_N, 1)
    alpha = float(coeffs[0])
    C = float(np.exp(coeffs[1]))

    log_N_pred = np.polyval(coeffs, log_Q)
    ss_res = float(np.sum((log_N - log_N_pred)**2))
    ss_tot = float(np.sum((log_N - np.mean(log_N))**2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return {
        'alpha': alpha,
        'C': C,
        'R_squared': r_squared,
        'formula': f'N_pauli = {C:.4f} * Q^{alpha:.4f}',
        'n_points': len(Q_arr),
        'source': _TRENEV_REFERENCE,
        'note': 'JW encoding with 2-qubit reduction',
    }


def _interpolate_gaussian_published(Q: int) -> int:
    """Interpolate Gaussian LiH Pauli count using published power-law fit."""
    fit = fit_gaussian_lih_published_exponent()
    return int(round(fit['C'] * Q ** fit['alpha']))


def build_equal_qubit_comparison_v3() -> Dict[str, Any]:
    """
    Recompute equal-qubit comparison using interpolated Gaussian Pauli counts
    from the published Trenev et al. power law (not Q^4.60 extrapolation).

    Saves results to debug/data/composed_lih_scaling_sweep_v3.json.

    Returns
    -------
    dict with sweep_data, Gaussian fit, equal-qubit comparison, and references.
    """
    gauss_fit = fit_gaussian_lih_published_exponent()

    geovac_points = [
        {'max_n': 1, 'Q': 6, 'N_pauli': 10},
        {'max_n': 2, 'Q': 30, 'N_pauli': 334},
        {'max_n': 3, 'Q': 84, 'N_pauli': 7879},
        {'max_n': 4, 'Q': 180, 'N_pauli': 92899},
    ]

    comparison = []
    for gv in geovac_points:
        Q_gv = gv['Q']
        N_gauss = int(round(gauss_fit['C'] * Q_gv ** gauss_fit['alpha']))
        ratio = round(N_gauss / max(1, gv['N_pauli']), 1)
        comparison.append({
            'max_n': gv['max_n'],
            'Q': Q_gv,
            'N_pauli_geovac': gv['N_pauli'],
            'N_pauli_gaussian_published_interp': N_gauss,
            'advantage_ratio': ratio,
            'note': f'Gaussian interpolated from published power law Q^{gauss_fit["alpha"]:.2f}',
        })

    output = {
        'sweep_data': geovac_points,
        'geovac_fit_4point': {
            'alpha': 2.6773,
            'a': 0.0617,
            'R_squared': 0.990,
            'formula': 'N_pauli = 0.0617 * Q^2.6773',
            'n_points': 4,
        },
        'gaussian_published_fit': gauss_fit,
        'gaussian_lih_published': {
            k: v for k, v in GAUSSIAN_LIH_PUBLISHED.items()
        },
        'gaussian_h2o_published': {
            k: v for k, v in GAUSSIAN_H2O_PUBLISHED.items()
        },
        'equal_qubit_comparison_v3': comparison,
        'reference': _TRENEV_REFERENCE,
        'notes': [
            'v3: Gaussian comparison uses published Trenev et al. data with 2-qubit reduction',
            'v2 used Q^4.60 extrapolation from STO-3G anchor (overestimated Gaussian counts)',
            'Published Gaussian STO-3G has Q=10 (2-qubit reduction), not Q=12',
            'Advantage ratios are ~7x-22x range, not 128x-1746x as in v2',
        ],
    }

    out_path = (Path(__file__).parent.parent
                / 'debug' / 'data' / 'composed_lih_scaling_sweep_v3.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)

    return output


# ---------------------------------------------------------------------------
# Composed BeH2 qubit Hamiltonian (triatomic extension)
# ---------------------------------------------------------------------------
#
# BeH2 is linear H-Be-H with 6 electrons:
#   - 2 core electrons on Be (1s², He-like, Z=4)
#   - 4 valence electrons (2 per Be-H bond)
#
# Approach A (independent blocks): each bond pair gets its own copy of
# screened-Be valence orbitals plus its own H orbitals. No orbital sharing
# between bonds. Cross-bond ERIs are zero by construction.
#
# Orbital layout:
#   [Core: Be Z=4] [Bond1-Be: Z_eff=2] [Bond1-H1: Z=1] [Bond2-Be: Z_eff=2] [Bond2-H2: Z=1]
#
# ERI blocks (5 independent):
#   1. Core (Z=4)
#   2. Bond1-Be (Z_eff=2)
#   3. Bond1-H1 (Z=1)
#   4. Bond2-Be (Z_eff=2) -- identical integrals to block 2
#   5. Bond2-H2 (Z=1) -- identical integrals to block 3
# ---------------------------------------------------------------------------

def build_composed_beh2(
    max_n_core: int = 2,
    max_n_val: int = 2,
    R: float = 2.51,
    E_core: Optional[float] = None,
    include_pk: bool = True,
    pk_in_hamiltonian: Optional[bool] = None,
    A_pk: Optional[float] = None,
    B_pk: Optional[float] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Build the composed BeH2 qubit Hamiltonian and count Pauli terms.

    Uses Approach A: independent bond blocks with no cross-bond ERIs.

    Orbital basis:
        Core block:   hydrogenic (n,l,m) on Be center, Z=4, up to max_n_core
        Bond 1 block: screened Be (Z_eff=2) + H1 (Z=1), up to max_n_val
        Bond 2 block: screened Be (Z_eff=2) + H2 (Z=1), up to max_n_val

    APPROXIMATIONS:
        1. Approach A: each bond gets independent copies of Be valence orbitals.
           Cross-bond ERIs are zero by construction (no orbital sharing).
        2. PK pseudopotential on bond-Be orbitals uses Z²-scaled Li²⁺ values
           (A=12.32, B=12.44) since ab initio Be2+ PK is not yet computed.
        3. Cross-center nuclear attraction integrals not included.
        4. E_core from Z²-scaled He-like estimate: -13.65 Ha.
        5. Valence Z_eff = Z_A - n_core = 2 (constant screening).

    Parameters
    ----------
    max_n_core : int
        Maximum n for core orbitals on Be (Z=4).
    max_n_val : int
        Maximum n for valence orbitals (screened-Be and H, per bond).
    R : float
        Be-H bond distance in bohr (default: 2.51, expt equilibrium).
    E_core : float or None
        Core energy in Ha. If None, uses -13.65 Ha (Be2+ estimate).
    include_pk : bool
        If True (default), compute PK on bond-Be orbitals.
    A_pk, B_pk : float or None
        PK parameters. If None, uses Z²-scaled He-like defaults.
    verbose : bool
        Print progress and results.

    Returns
    -------
    dict with keys:
        M, Q, N_pauli, M_core, M_bond, M_bond_be, M_bond_h,
        ERI_density_total, h1, eri, nuclear_repulsion, qubit_op, etc.
    """
    t0 = time.perf_counter()

    # --- Physical parameters ---
    Z_A = 4      # Be nuclear charge
    Z_B = 1      # H nuclear charge (both H atoms)
    n_core_electrons = 2
    Z_eff_val = Z_A - n_core_electrons  # = 2

    if E_core is None:
        # Be2+ (He-like, Z=4) ground state estimate
        # Z²-scaling from He: E(He) = -2.9037, E(Be2+) ≈ -2.9037 x (4/2)² = -11.61
        # Better estimate using variational He-like: -(Z-5/16)² ≈ -13.60
        # Prompt suggests -13.65 Ha
        E_core = -13.65

    # --- Define orbital basis ---
    states_core = _enumerate_states(max_n_core)
    M_core = len(states_core)

    # Bond blocks: each bond has screened-Be + H orbitals
    states_bond_be = _enumerate_states(max_n_val)   # Z_eff orbitals
    states_bond_h = _enumerate_states(max_n_val)    # Z=1 orbitals
    M_bond_be = len(states_bond_be)
    M_bond_h = len(states_bond_h)
    M_bond = M_bond_be + M_bond_h  # per bond pair

    M = M_core + 2 * M_bond
    Q = 2 * M

    # Offsets into the full orbital array
    offset_core = 0
    offset_b1_be = M_core
    offset_b1_h = M_core + M_bond_be
    offset_b2_be = M_core + M_bond
    offset_b2_h = M_core + M_bond + M_bond_be

    if verbose:
        print(f"[composed_beh2] Orbital basis:")
        print(f"  Core (Z={Z_A}):           {M_core} orbitals, max_n={max_n_core}")
        print(f"  Bond-Be (Z_eff={Z_eff_val}):    {M_bond_be} orbitals x 2 bonds")
        print(f"  Bond-H  (Z={Z_B}):         {M_bond_h} orbitals x 2 bonds")
        print(f"  Per bond:              {M_bond} orbitals")
        print(f"  Total: M={M} spatial, Q={Q} qubits")

    # --- Construct h1 (M x M) ---
    h1 = np.zeros((M, M))

    # Core block: -Z_A^2 / (2n^2)
    for i, (n, l, m) in enumerate(states_core):
        h1[offset_core + i, offset_core + i] = -Z_A**2 / (2.0 * n**2)

    # Bond 1 -- Be valence: -Z_eff^2 / (2n^2)
    for i, (n, l, m) in enumerate(states_bond_be):
        h1[offset_b1_be + i, offset_b1_be + i] = -Z_eff_val**2 / (2.0 * n**2)

    # Bond 1 -- H1: -Z_B^2 / (2n^2)
    for i, (n, l, m) in enumerate(states_bond_h):
        h1[offset_b1_h + i, offset_b1_h + i] = -Z_B**2 / (2.0 * n**2)

    # Bond 2 -- Be valence: identical to bond 1 Be
    for i, (n, l, m) in enumerate(states_bond_be):
        h1[offset_b2_be + i, offset_b2_be + i] = -Z_eff_val**2 / (2.0 * n**2)

    # Bond 2 -- H2: identical to bond 1 H
    for i, (n, l, m) in enumerate(states_bond_h):
        h1[offset_b2_h + i, offset_b2_h + i] = -Z_B**2 / (2.0 * n**2)

    # --- PK pseudopotential on bond-Be orbitals ---
    if pk_in_hamiltonian is None:
        pk_in_hamiltonian = include_pk  # backward compat

    h1_pk_full = np.zeros((M, M))
    n_pk_nonzero = 0
    if include_pk:
        if A_pk is None or B_pk is None:
            pk_defaults = _PK_HELIKE_DEFAULTS.get(Z_A)
            if pk_defaults is None:
                raise ValueError(
                    f"No He-like PK parameters for Z={Z_A}. "
                    f"Provide A_pk and B_pk explicitly."
                )
            if A_pk is None:
                A_pk = pk_defaults['A']
            if B_pk is None:
                B_pk = pk_defaults['B']

        if verbose:
            print(f"[composed_beh2] PK pseudopotential: A={A_pk:.2f}, B={B_pk:.2f}")
            print(f"  Source: {_PK_HELIKE_DEFAULTS.get(Z_A, {}).get('source', 'user-supplied')}")

        # PK on bond-Be orbitals (same for both bonds)
        h1_pk_bond_be = _compute_pk_matrix_elements(
            Z_eff_val, states_bond_be, A_pk, B_pk,
        )

        # Store in full-size PK matrix (both bond blocks)
        for i in range(M_bond_be):
            for j in range(M_bond_be):
                if abs(h1_pk_bond_be[i, j]) > 1e-15:
                    h1_pk_full[offset_b1_be + i, offset_b1_be + j] = h1_pk_bond_be[i, j]
                    h1_pk_full[offset_b2_be + i, offset_b2_be + j] = h1_pk_bond_be[i, j]
                    n_pk_nonzero += 2

        # Only add PK to h1 if requested
        if pk_in_hamiltonian:
            h1 += h1_pk_full

        if verbose:
            diag_pk = np.diag(h1_pk_bond_be)
            print(f"  Bond-Be PK nonzero entries: {n_pk_nonzero} (across 2 bonds)")
            print(f"  PK diagonal range: [{diag_pk.min():.6f}, {diag_pk.max():.6f}] Ha")
            print(f"  PK in quantum Hamiltonian: {pk_in_hamiltonian}")

    if verbose:
        if pk_in_hamiltonian:
            h1_desc = "diagonal + PK"
        elif include_pk:
            h1_desc = "diagonal only (PK computed separately)"
        else:
            h1_desc = "diagonal only"
        print(f"[composed_beh2] h1 constructed ({M}x{M}), {h1_desc}")

    # --- Construct ERI (M x M x M x M) in chemist notation ---
    eri = np.zeros((M, M, M, M))

    # Block 1: Core (Z=4)
    if verbose:
        print(f"[composed_beh2] Computing core R^k integrals (Z={Z_A})...")
    rk_core = _compute_rk_integrals_block(Z_A, states_core)
    eri_core_phys = _build_eri_block(Z_A, states_core, rk_core)
    n_eri_core = len(eri_core_phys)

    for (a, b, c, d), val in eri_core_phys.items():
        eri[a, c, b, d] = val

    # Block 2: Bond1-Be (Z_eff=2)
    if verbose:
        print(f"[composed_beh2] Computing bond-Be R^k integrals (Z_eff={Z_eff_val})...")
    rk_bond_be = _compute_rk_integrals_block(Z_eff_val, states_bond_be)
    eri_bond_be_phys = _build_eri_block(Z_eff_val, states_bond_be, rk_bond_be)
    n_eri_bond_be = len(eri_bond_be_phys)

    for (a, b, c, d), val in eri_bond_be_phys.items():
        p = a + offset_b1_be
        q = c + offset_b1_be
        r = b + offset_b1_be
        s = d + offset_b1_be
        eri[p, q, r, s] = val

    # Block 3: Bond1-H1 (Z=1)
    if verbose:
        print(f"[composed_beh2] Computing bond-H R^k integrals (Z={Z_B})...")
    rk_bond_h = _compute_rk_integrals_block(Z_B, states_bond_h)
    eri_bond_h_phys = _build_eri_block(Z_B, states_bond_h, rk_bond_h)
    n_eri_bond_h = len(eri_bond_h_phys)

    for (a, b, c, d), val in eri_bond_h_phys.items():
        p = a + offset_b1_h
        q = c + offset_b1_h
        r = b + offset_b1_h
        s = d + offset_b1_h
        eri[p, q, r, s] = val

    # Block 4: Bond2-Be (Z_eff=2) -- same integrals as block 2
    for (a, b, c, d), val in eri_bond_be_phys.items():
        p = a + offset_b2_be
        q = c + offset_b2_be
        r = b + offset_b2_be
        s = d + offset_b2_be
        eri[p, q, r, s] = val

    # Block 5: Bond2-H2 (Z=1) -- same integrals as block 3
    for (a, b, c, d), val in eri_bond_h_phys.items():
        p = a + offset_b2_h
        q = c + offset_b2_h
        r = b + offset_b2_h
        s = d + offset_b2_h
        eri[p, q, r, s] = val

    # Symmetrize: enforce (pq|rs) = (rs|pq)
    eri = 0.5 * (eri + eri.transpose(2, 3, 0, 1))

    # --- ERI density statistics ---
    n_eri_total = int(np.count_nonzero(np.abs(eri) > 1e-15))
    eri_density_total = n_eri_total / max(1, M**4)
    eri_density_core = n_eri_core / max(1, M_core**4)
    n_eri_bond_total = 2 * (n_eri_bond_be + n_eri_bond_h)
    n_eri_per_bond = n_eri_bond_be + n_eri_bond_h

    if verbose:
        print(f"[composed_beh2] ERI statistics:")
        print(f"  Core:     {n_eri_core} nonzero / {M_core**4}")
        print(f"  Bond-Be:  {n_eri_bond_be} x 2 = {2*n_eri_bond_be} (Z_eff={Z_eff_val})")
        print(f"  Bond-H:   {n_eri_bond_h} x 2 = {2*n_eri_bond_h} (Z={Z_B})")
        print(f"  Total:    {n_eri_total} nonzero / {M**4} = {eri_density_total:.2%}")

    # --- Nuclear repulsion ---
    # Linear geometry: H1 at -R, Be at 0, H2 at +R
    # V_NN = Z_Be*Z_H/R (x2 for two Be-H bonds) + Z_H*Z_H/(2R) (H-H)
    V_NN = 2.0 * Z_A * Z_B / R + Z_B * Z_B / (2.0 * R)

    # Core-to-H attractions (2 symmetric contributions)
    V_cross_1 = _v_cross_nuc_1s(Z_A, n_core_electrons, Z_B, R)
    V_cross = 2.0 * V_cross_1  # symmetric: same R for both H atoms

    nuclear_repulsion = V_NN + V_cross + E_core

    if verbose:
        print(f"[composed_beh2] Energy constants:")
        print(f"  V_NN (Be-Hx2 + H-H)  = {V_NN:.6f} Ha")
        print(f"  V_cross_nuc (x2)      = {V_cross:.6f} Ha")
        print(f"  E_core (Be2+)         = {E_core:.6f} Ha")
        print(f"  nuc_repul             = {nuclear_repulsion:.6f} Ha")

    # --- JW transform and Pauli count ---
    if verbose:
        print(f"[composed_beh2] Building fermion operator...")
    fermion_op = build_fermion_op_from_integrals(h1, eri, nuclear_repulsion)

    if verbose:
        print(f"[composed_beh2] Jordan-Wigner transform...")
    qubit_op = jordan_wigner(fermion_op)
    N_pauli = len(qubit_op.terms)

    elapsed = time.perf_counter() - t0

    if verbose:
        print(f"\n{'='*60}")
        print(f"COMPOSED BeH2 QUBIT HAMILTONIAN -- RESULTS")
        print(f"{'='*60}")
        print(f"  Spatial orbitals M = {M}  (core={M_core}, bond={M_bond}x2)")
        print(f"  Qubits           Q = {Q}")
        print(f"  Pauli terms      N = {N_pauli:,}")
        print(f"  ERI density (total)   = {eri_density_total:.2%}")
        print(f"  ERI blocks: {n_eri_core} core + {n_eri_bond_be}x2 bond-Be + {n_eri_bond_h}x2 bond-H")
        print(f"  PK included          = {include_pk}"
              f"  (nonzero entries: {n_pk_nonzero})")
        print(f"  Wall time = {elapsed:.1f}s")
        print(f"{'='*60}")

    # --- Build results dict ---
    results = {
        'M': M,
        'Q': Q,
        'N_pauli': N_pauli,
        'M_core': M_core,
        'M_bond': M_bond,
        'M_bond_be': M_bond_be,
        'M_bond_h': M_bond_h,
        'ERI_density_total': eri_density_total,
        'ERI_density_core': eri_density_core,
        'n_eri_core': n_eri_core,
        'n_eri_bond_be': n_eri_bond_be,
        'n_eri_bond_h': n_eri_bond_h,
        'n_eri_per_bond': n_eri_per_bond,
        'n_eri_bond_total': n_eri_bond_total,
        'n_eri_total': n_eri_total,
        'nuclear_repulsion': nuclear_repulsion,
        'V_NN': V_NN,
        'V_cross_nuc': V_cross,
        'E_core': E_core,
        'R_bohr': R,
        'max_n_core': max_n_core,
        'max_n_val': max_n_val,
        'Z_A': Z_A,
        'Z_B': Z_B,
        'Z_eff_val': Z_eff_val,
        'wall_time_s': elapsed,
        'include_pk': include_pk,
        'pk_in_hamiltonian': pk_in_hamiltonian,
        'A_pk': A_pk if include_pk else None,
        'B_pk': B_pk if include_pk else None,
        'n_pk_nonzero': n_pk_nonzero,
        'approximations': [
            'Approach A: independent bond blocks, no cross-bond ERIs',
            ('PK in quantum Hamiltonian' if pk_in_hamiltonian
             else 'PK computed separately (classical correction)')
            if include_pk else 'PK disabled',
            'PK on bond-H skipped (two-center, negligible at R~2.5 bohr)',
            'Cross-center nuclear attraction skipped',
            'Cross-block ERI = 0 (pseudopotential + Approach A)',
            f'E_core = {E_core} Ha (Be2+ He-like estimate)',
            f'Valence Z_eff = {Z_eff_val} (constant screening)',
        ],
    }

    # Save JSON report
    out_path = (Path(__file__).parent.parent
                / 'debug' / 'data' / 'composed_beh2_pauli_analysis.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    json_results = {k: v for k, v in results.items()
                    if not isinstance(v, np.ndarray)}
    with open(out_path, 'w') as f:
        json.dump(json_results, f, indent=2)
    if verbose:
        print(f"\nResults saved to {out_path}")

    # Attach non-serializable objects
    results['h1'] = h1
    results['h1_pk'] = h1_pk_full
    results['eri'] = eri
    results['qubit_op'] = qubit_op
    results['fermion_op'] = fermion_op
    results['states_core'] = states_core
    results['states_bond_be'] = states_bond_be
    results['states_bond_h'] = states_bond_h

    return results


# ---------------------------------------------------------------------------
# Cross-bond ERI estimate for BeH2
# ---------------------------------------------------------------------------

def estimate_cross_bond_eri_count(
    max_n: int,
) -> Dict[str, Any]:
    """
    Estimate cross-bond ERIs that would exist if bond 1 and bond 2 shared
    Be valence orbitals (Approach B). Uses Gaunt selection rules.

    Parameters
    ----------
    max_n : int
        Maximum principal quantum number for valence orbitals.

    Returns
    -------
    dict with n_cross_bond, n_same_bond, ratio, M_bond.
    """
    states = _enumerate_states(max_n)
    M_per_center = len(states)
    M_bond = 2 * M_per_center  # Be-val + H per bond

    # Use the same-center ERI counting logic from estimate_cross_center_eri_count
    same_center_result = estimate_cross_center_eri_count(max_n)
    n_same_center = same_center_result['n_same_center']

    # Cross-bond ERIs: a,c from bond1 orbital set; b,d from bond2 orbital set.
    # Since both bond sets have identical angular quantum numbers,
    # the selection rule count is the same as cross-center.
    n_cross_bond = same_center_result['n_cross_center']

    return {
        'max_n': max_n,
        'n_cross_bond': n_cross_bond,
        'n_same_bond': n_same_center,
        'ratio': n_cross_bond / max(1, n_same_center),
        'M_bond': M_bond,
        'M_per_center': M_per_center,
        'upper_bound': M_bond ** 4,
        'note': ('Cross-bond ERIs would arise from shared Be valence orbitals '
                 'in Approach B. Approach A eliminates these by construction.'),
    }


# ---------------------------------------------------------------------------
# LiH vs BeH2 scaling comparison
# ---------------------------------------------------------------------------

def compare_lih_beh2_scaling(
    max_n_values: Optional[List[int]] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Run both LiH and BeH2 at multiple basis sizes, compare Pauli scaling.

    Parameters
    ----------
    max_n_values : list of int or None
        Basis sizes to sweep. Default [1, 2, 3].
    verbose : bool
        Print comparison table.

    Returns
    -------
    dict with sweep data, fitted exponents, and comparison.
    """
    if max_n_values is None:
        max_n_values = [1, 2, 3]

    lih_data: List[Dict[str, Any]] = []
    beh2_data: List[Dict[str, Any]] = []

    for n in max_n_values:
        if verbose:
            print(f"\n{'='*60}")
            print(f"  max_n = {n}: Building LiH...")
            print(f"{'='*60}")
        lih = build_composed_lih(max_n_core=n, max_n_val=n, verbose=verbose)
        lih_data.append({
            'max_n': n, 'M': lih['M'], 'Q': lih['Q'],
            'N_pauli': lih['N_pauli'],
            'ERI_density': lih['ERI_density_total'],
            'wall_time_s': lih['wall_time_s'],
        })

        if verbose:
            print(f"\n{'='*60}")
            print(f"  max_n = {n}: Building BeH2...")
            print(f"{'='*60}")
        beh2 = build_composed_beh2(max_n_core=n, max_n_val=n, verbose=verbose)
        beh2_data.append({
            'max_n': n, 'M': beh2['M'], 'Q': beh2['Q'],
            'N_pauli': beh2['N_pauli'],
            'ERI_density': beh2['ERI_density_total'],
            'wall_time_s': beh2['wall_time_s'],
            'n_eri_per_bond': beh2['n_eri_per_bond'],
        })

    # Fit power laws
    def _fit_power_law(data: List[Dict]) -> Dict[str, float]:
        Q_arr = np.array([d['Q'] for d in data], dtype=float)
        N_arr = np.array([d['N_pauli'] for d in data], dtype=float)
        if len(Q_arr) < 2 or not np.all(Q_arr > 0) or not np.all(N_arr > 0):
            return {'alpha': np.nan, 'a': np.nan, 'R_squared': np.nan}
        log_Q = np.log(Q_arr)
        log_N = np.log(N_arr)
        coeffs = np.polyfit(log_Q, log_N, 1)
        alpha = float(coeffs[0])
        a = float(np.exp(coeffs[1]))
        log_N_pred = np.polyval(coeffs, log_Q)
        ss_res = np.sum((log_N - log_N_pred)**2)
        ss_tot = np.sum((log_N - np.mean(log_N))**2)
        r_sq = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 1.0
        return {'alpha': alpha, 'a': a, 'R_squared': r_sq}

    lih_fit = _fit_power_law(lih_data)
    beh2_fit = _fit_power_law(beh2_data)

    # Cross-bond ERI estimates
    cross_bond_estimates = []
    for n in max_n_values:
        cross_bond_estimates.append(estimate_cross_bond_eri_count(n))

    if verbose:
        print(f"\n{'='*70}")
        print(f"  LiH vs BeH2 SCALING COMPARISON")
        print(f"{'='*70}")
        print(f"\n  {'System':>8} {'max_n':>5} {'M':>4} {'Q':>5} {'N_pauli':>10}"
              f" {'ERI_dens':>10} {'time':>8}")
        print(f"  {'-'*8} {'-'*5} {'-'*4} {'-'*5} {'-'*10}"
              f" {'-'*10} {'-'*8}")
        for ld, bd in zip(lih_data, beh2_data):
            print(f"  {'LiH':>8} {ld['max_n']:>5} {ld['M']:>4} {ld['Q']:>5}"
                  f" {ld['N_pauli']:>10,} {ld['ERI_density']:>10.2%}"
                  f" {ld['wall_time_s']:>7.1f}s")
            print(f"  {'BeH2':>8} {bd['max_n']:>5} {bd['M']:>4} {bd['Q']:>5}"
                  f" {bd['N_pauli']:>10,} {bd['ERI_density']:>10.2%}"
                  f" {bd['wall_time_s']:>7.1f}s")

        print(f"\n  Fitted exponents:")
        print(f"    LiH:  N_pauli = {lih_fit['a']:.4f} x Q^{lih_fit['alpha']:.4f}"
              f"  (R² = {lih_fit['R_squared']:.4f})")
        print(f"    BeH2: N_pauli = {beh2_fit['a']:.4f} x Q^{beh2_fit['alpha']:.4f}"
              f"  (R² = {beh2_fit['R_squared']:.4f})")
        print(f"    Gaussian molecular: Q^4.25 (Trenev et al.)")

        print(f"\n  Key question: does BeH2 fall on the same scaling law as LiH?")
        if abs(lih_fit['alpha'] - beh2_fit['alpha']) < 0.3:
            print(f"    -> YES, exponents within 0.3: d_alpha = "
                  f"{abs(lih_fit['alpha'] - beh2_fit['alpha']):.3f}")
        else:
            print(f"    -> SHIFT detected: d_alpha = "
                  f"{abs(lih_fit['alpha'] - beh2_fit['alpha']):.3f}")

        print(f"\n  Cross-bond ERI estimates (Approach B impact):")
        for cb in cross_bond_estimates:
            print(f"    max_n={cb['max_n']}: {cb['n_cross_bond']} cross-bond ERIs"
                  f" vs {cb['n_same_bond']} same-bond, ratio={cb['ratio']:.2f}")

        print(f"{'='*70}")

    # Save results
    output = {
        'lih_data': lih_data,
        'beh2_data': beh2_data,
        'lih_fit': lih_fit,
        'beh2_fit': beh2_fit,
        'cross_bond_estimates': cross_bond_estimates,
        'gaussian_exponent': 4.25,
    }

    out_path = (Path(__file__).parent.parent
                / 'debug' / 'data' / 'composed_beh2_pauli_analysis.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    if verbose:
        print(f"\nComparison saved to {out_path}")

    return output


# ---------------------------------------------------------------------------
# Composed H2O qubit Hamiltonian (bent triatomic, 7 blocks)
# ---------------------------------------------------------------------------
#
# H2O has 10 electrons:
#   - 2 core electrons on O (1s², He-like, Z=8)
#   - 8 valence electrons:
#       2 per O-H bond (x2 bonds)
#       2 per lone pair (x2 lone pairs)
#
# Approach A (independent blocks): each valence group gets its own copy
# of screened-O orbitals (Z_eff=6) plus H orbitals where applicable.
# Cross-block ERIs are zero by construction.
#
# Orbital layout (7 blocks):
#   [Core: O Z=8]
#   [Bond1-O: Z_eff=6] [Bond1-H1: Z=1]
#   [Bond2-O: Z_eff=6] [Bond2-H2: Z=1]
#   [Lone1-O: Z_eff=6]
#   [Lone2-O: Z_eff=6]
#
# ERI blocks (7 independent, but 4 O-side valence blocks share Z_eff=6):
#   1. Core (Z=8)
#   2. Bond1-O (Z_eff=6) -- same integrals as blocks 4, 6, 7
#   3. Bond1-H1 (Z=1) -- same integrals as block 5
#   4. Bond2-O (Z_eff=6) -- reuse block 2
#   5. Bond2-H2 (Z=1) -- reuse block 3
#   6. Lone1-O (Z_eff=6) -- reuse block 2
#   7. Lone2-O (Z_eff=6) -- reuse block 2
# ---------------------------------------------------------------------------

def build_composed_h2o(
    max_n_core: int = 2,
    max_n_val: int = 2,
    R_OH: float = 1.809,
    angle_HOH: float = 104.5,
    E_core: Optional[float] = None,
    include_pk: bool = True,
    pk_in_hamiltonian: Optional[bool] = None,
    A_pk: Optional[float] = None,
    B_pk: Optional[float] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Build the composed H2O qubit Hamiltonian and count Pauli terms.

    Uses Approach A: 7 independent orbital blocks with no cross-block ERIs.

    Orbital basis (7 blocks):
        Core:    O hydrogenic (Z=8) up to max_n_core
        Bond 1:  screened O (Z_eff=6) + H1 (Z=1), up to max_n_val
        Bond 2:  screened O (Z_eff=6) + H2 (Z=1), up to max_n_val
        Lone 1:  screened O (Z_eff=6), up to max_n_val
        Lone 2:  screened O (Z_eff=6), up to max_n_val

    APPROXIMATIONS:
        1. Approach A: each valence group gets independent copies of O valence
           orbitals. Cross-block ERIs are zero by construction.
        2. PK pseudopotential on all O-side valence blocks uses Z²-scaled Li²⁺
           values: A=49.28, B=49.78 (6.93*(8/3)², 7.00*(8/3)²).
        3. PK on bond-H orbitals skipped (two-center, negligible at R~1.8 bohr).
        4. Cross-center nuclear attraction integrals not included.
        5. E_core from first-order He-like estimate for O²⁺ (Z=8):
           E ≈ -Z² + (5/8)Z = -64 + 5 = -59.0 Ha. This uses the standard
           first-order perturbation theory result for the He-like isoelectronic
           sequence (variational: -(Z-5/16)² = -59.16 Ha). We use -59.16 Ha.
        6. Valence Z_eff = Z_O - 2 = 6 (constant screening).
        7. The 4 O-side valence blocks (bond1-O, bond2-O, lone1, lone2) all
           have identical R^k integrals (same Z_eff=6, same orbital set),
           computed once and reused 4 times.

    Parameters
    ----------
    max_n_core : int
        Maximum n for core orbitals on O (Z=8).
    max_n_val : int
        Maximum n for valence orbitals (screened-O and H, per group).
    R_OH : float
        O-H bond distance in bohr (default: 1.809, expt 0.957 Å).
    angle_HOH : float
        H-O-H bond angle in degrees (default: 104.5°, experimental).
    E_core : float or None
        Core energy in Ha. If None, uses -59.16 Ha (O²⁺ variational He-like
        estimate: -(Z-5/16)² = -(8-5/16)² = -59.16 Ha).
    include_pk : bool
        If True (default), compute PK on all O-side valence blocks.
    A_pk, B_pk : float or None
        PK parameters. If None, uses Z²-scaled He-like defaults for Z=8.
    verbose : bool
        Print progress and results.

    Returns
    -------
    dict with keys:
        M, Q, N_pauli, M_core, M_bond, M_bond_o, M_bond_h, M_lone,
        ERI_density_total, h1, eri, nuclear_repulsion, qubit_op, etc.
    """
    t0 = time.perf_counter()

    # --- Physical parameters ---
    Z_O = 8      # Oxygen nuclear charge
    Z_H = 1      # Hydrogen nuclear charge
    n_core_electrons = 2
    Z_eff_val = Z_O - n_core_electrons  # = 6

    if E_core is None:
        # O²⁺ (He-like, Z=8) variational estimate: -(Z-5/16)²
        E_core = -(Z_O - 5.0 / 16.0) ** 2

    # Geometry
    angle_rad = np.radians(angle_HOH)
    R_HH = 2.0 * R_OH * np.sin(angle_rad / 2.0)

    # --- Define orbital basis ---
    states_core = _enumerate_states(max_n_core)
    M_core = len(states_core)

    # Valence sub-blocks
    states_val_o = _enumerate_states(max_n_val)   # Z_eff=6 orbitals on O
    states_val_h = _enumerate_states(max_n_val)    # Z=1 orbitals on H
    M_val_o = len(states_val_o)
    M_val_h = len(states_val_h)
    M_bond = M_val_o + M_val_h  # per bond pair (O-side + H-side)
    M_lone = M_val_o            # lone pair (O-side only)

    M = M_core + 2 * M_bond + 2 * M_lone
    Q = 2 * M

    # Offsets into the full orbital array
    offset_core = 0
    offset_b1_o = M_core
    offset_b1_h = M_core + M_val_o
    offset_b2_o = M_core + M_bond
    offset_b2_h = M_core + M_bond + M_val_o
    offset_lp1 = M_core + 2 * M_bond
    offset_lp2 = M_core + 2 * M_bond + M_lone

    if verbose:
        print(f"[composed_h2o] Orbital basis:")
        print(f"  Core (Z={Z_O}):           {M_core} orbitals, max_n={max_n_core}")
        print(f"  Bond-O (Z_eff={Z_eff_val}):    {M_val_o} orbitals x 2 bonds")
        print(f"  Bond-H  (Z={Z_H}):         {M_val_h} orbitals x 2 bonds")
        print(f"  Lone-O (Z_eff={Z_eff_val}):    {M_val_o} orbitals x 2 lone pairs")
        print(f"  Per bond:              {M_bond} orbitals")
        print(f"  Per lone pair:         {M_lone} orbitals")
        print(f"  Total: M={M} spatial, Q={Q} qubits")

    # --- Construct h1 (M x M) ---
    h1 = np.zeros((M, M))

    # Core block: -Z_O^2 / (2n^2)
    for i, (n, l, m) in enumerate(states_core):
        h1[offset_core + i, offset_core + i] = -Z_O**2 / (2.0 * n**2)

    # Bond 1 O-side: -Z_eff^2 / (2n^2)
    for i, (n, l, m) in enumerate(states_val_o):
        h1[offset_b1_o + i, offset_b1_o + i] = -Z_eff_val**2 / (2.0 * n**2)

    # Bond 1 H1: -Z_H^2 / (2n^2)
    for i, (n, l, m) in enumerate(states_val_h):
        h1[offset_b1_h + i, offset_b1_h + i] = -Z_H**2 / (2.0 * n**2)

    # Bond 2 O-side
    for i, (n, l, m) in enumerate(states_val_o):
        h1[offset_b2_o + i, offset_b2_o + i] = -Z_eff_val**2 / (2.0 * n**2)

    # Bond 2 H2
    for i, (n, l, m) in enumerate(states_val_h):
        h1[offset_b2_h + i, offset_b2_h + i] = -Z_H**2 / (2.0 * n**2)

    # Lone pair 1
    for i, (n, l, m) in enumerate(states_val_o):
        h1[offset_lp1 + i, offset_lp1 + i] = -Z_eff_val**2 / (2.0 * n**2)

    # Lone pair 2
    for i, (n, l, m) in enumerate(states_val_o):
        h1[offset_lp2 + i, offset_lp2 + i] = -Z_eff_val**2 / (2.0 * n**2)

    # --- PK pseudopotential on all O-side valence blocks ---
    if pk_in_hamiltonian is None:
        pk_in_hamiltonian = include_pk  # backward compat

    h1_pk_full = np.zeros((M, M))
    n_pk_nonzero = 0
    if include_pk:
        if A_pk is None or B_pk is None:
            pk_defaults = _PK_HELIKE_DEFAULTS.get(Z_O)
            if pk_defaults is None:
                raise ValueError(
                    f"No He-like PK parameters for Z={Z_O}. "
                    f"Provide A_pk and B_pk explicitly."
                )
            if A_pk is None:
                A_pk = pk_defaults['A']
            if B_pk is None:
                B_pk = pk_defaults['B']

        if verbose:
            print(f"[composed_h2o] PK pseudopotential: A={A_pk:.2f}, B={B_pk:.2f}")
            print(f"  Source: {_PK_HELIKE_DEFAULTS.get(Z_O, {}).get('source', 'user-supplied')}")

        # PK on O-side valence orbitals (same for all 4 O-side blocks)
        h1_pk_val_o = _compute_pk_matrix_elements(
            Z_eff_val, states_val_o, A_pk, B_pk,
        )

        # Store in full-size PK matrix (all 4 O-side blocks)
        for offset in [offset_b1_o, offset_b2_o, offset_lp1, offset_lp2]:
            for i in range(M_val_o):
                for j in range(M_val_o):
                    if abs(h1_pk_val_o[i, j]) > 1e-15:
                        h1_pk_full[offset + i, offset + j] = h1_pk_val_o[i, j]
                        n_pk_nonzero += 1

        # Only add PK to h1 if requested
        if pk_in_hamiltonian:
            h1 += h1_pk_full

        if verbose:
            diag_pk = np.diag(h1_pk_val_o)
            print(f"  O-side PK nonzero entries: {n_pk_nonzero} (across 4 blocks)")
            print(f"  PK diagonal range: [{diag_pk.min():.6f}, {diag_pk.max():.6f}] Ha")
            print(f"  PK in quantum Hamiltonian: {pk_in_hamiltonian}")

    if verbose:
        if pk_in_hamiltonian:
            h1_desc = "diagonal + PK"
        elif include_pk:
            h1_desc = "diagonal only (PK computed separately)"
        else:
            h1_desc = "diagonal only"
        print(f"[composed_h2o] h1 constructed ({M}x{M}), {h1_desc}")

    # --- Construct ERI (M x M x M x M) in chemist notation ---
    eri = np.zeros((M, M, M, M))

    # Block 1: Core (Z=8)
    if verbose:
        print(f"[composed_h2o] Computing core R^k integrals (Z={Z_O})...")
    rk_core = _compute_rk_integrals_block(Z_O, states_core)
    eri_core_phys = _build_eri_block(Z_O, states_core, rk_core)
    n_eri_core = len(eri_core_phys)

    for (a, b, c, d), val in eri_core_phys.items():
        eri[a, c, b, d] = val

    # Blocks 2,4,6,7: O-side valence (Z_eff=6) -- compute ONCE, place 4 times
    if verbose:
        print(f"[composed_h2o] Computing O-side valence R^k integrals (Z_eff={Z_eff_val})...")
    rk_val_o = _compute_rk_integrals_block(Z_eff_val, states_val_o)
    eri_val_o_phys = _build_eri_block(Z_eff_val, states_val_o, rk_val_o)
    n_eri_val_o = len(eri_val_o_phys)

    for offset in [offset_b1_o, offset_b2_o, offset_lp1, offset_lp2]:
        for (a, b, c, d), val in eri_val_o_phys.items():
            p = a + offset
            q = c + offset
            r = b + offset
            s = d + offset
            eri[p, q, r, s] = val

    # Blocks 3,5: Bond-H (Z=1) -- compute ONCE, place 2 times
    if verbose:
        print(f"[composed_h2o] Computing bond-H R^k integrals (Z={Z_H})...")
    rk_val_h = _compute_rk_integrals_block(Z_H, states_val_h)
    eri_val_h_phys = _build_eri_block(Z_H, states_val_h, rk_val_h)
    n_eri_val_h = len(eri_val_h_phys)

    for offset in [offset_b1_h, offset_b2_h]:
        for (a, b, c, d), val in eri_val_h_phys.items():
            p = a + offset
            q = c + offset
            r = b + offset
            s = d + offset
            eri[p, q, r, s] = val

    # Symmetrize: enforce (pq|rs) = (rs|pq)
    eri = 0.5 * (eri + eri.transpose(2, 3, 0, 1))

    # --- ERI density statistics ---
    n_eri_total = int(np.count_nonzero(np.abs(eri) > 1e-15))
    eri_density_total = n_eri_total / max(1, M**4)
    eri_density_core = n_eri_core / max(1, M_core**4)
    n_eri_val_o_total = 4 * n_eri_val_o   # 4 O-side blocks
    n_eri_val_h_total = 2 * n_eri_val_h   # 2 H blocks
    n_eri_blocks_total = n_eri_core + n_eri_val_o_total + n_eri_val_h_total

    if verbose:
        print(f"[composed_h2o] ERI statistics:")
        print(f"  Core:    {n_eri_core} nonzero / {M_core**4}")
        print(f"  Val-O:   {n_eri_val_o} x 4 = {n_eri_val_o_total} (Z_eff={Z_eff_val})")
        print(f"  Val-H:   {n_eri_val_h} x 2 = {n_eri_val_h_total} (Z={Z_H})")
        print(f"  Total:   {n_eri_total} nonzero / {M**4} = {eri_density_total:.2%}")

    # --- Nuclear repulsion ---
    # V_OH: O-H repulsion (x2 for two O-H bonds)
    V_OH = 2.0 * Z_O * Z_H / R_OH
    # V_HH: H-H repulsion
    V_HH = Z_H * Z_H / R_HH

    V_NN = V_OH + V_HH

    # Core-to-H attractions: core 1s² electrons on O attract both H nuclei
    V_cross_1 = _v_cross_nuc_1s(Z_O, n_core_electrons, Z_H, R_OH)
    V_cross = 2.0 * V_cross_1  # same R_OH for both H atoms

    nuclear_repulsion = V_NN + V_cross + E_core

    if verbose:
        print(f"[composed_h2o] Energy constants:")
        print(f"  R_OH = {R_OH:.4f} bohr, angle = {angle_HOH:.1f} deg")
        print(f"  R_HH = {R_HH:.4f} bohr")
        print(f"  V_OH (x2)     = {V_OH:.6f} Ha")
        print(f"  V_HH          = {V_HH:.6f} Ha")
        print(f"  V_NN          = {V_NN:.6f} Ha")
        print(f"  V_cross_nuc (x2) = {V_cross:.6f} Ha")
        print(f"  E_core (O2+)  = {E_core:.6f} Ha")
        print(f"  nuc_repul     = {nuclear_repulsion:.6f} Ha")

    # --- JW transform and Pauli count ---
    if verbose:
        print(f"[composed_h2o] Building fermion operator...")
    fermion_op = build_fermion_op_from_integrals(h1, eri, nuclear_repulsion)

    if verbose:
        print(f"[composed_h2o] Jordan-Wigner transform...")
    qubit_op = jordan_wigner(fermion_op)
    N_pauli = len(qubit_op.terms)

    elapsed = time.perf_counter() - t0

    if verbose:
        print(f"\n{'='*60}")
        print(f"COMPOSED H2O QUBIT HAMILTONIAN -- RESULTS")
        print(f"{'='*60}")
        print(f"  Spatial orbitals M = {M}  (core={M_core},"
              f" bond={M_bond}x2, lone={M_lone}x2)")
        print(f"  Qubits           Q = {Q}")
        print(f"  Pauli terms      N = {N_pauli:,}")
        print(f"  ERI density (total)   = {eri_density_total:.2%}")
        print(f"  ERI blocks: {n_eri_core} core"
              f" + {n_eri_val_o}x4 val-O + {n_eri_val_h}x2 val-H")
        print(f"  PK included          = {include_pk}"
              f"  (nonzero entries: {n_pk_nonzero})")
        print(f"  Gaussian H2O STO-3G  = {GAUSSIAN_H2O_PUBLISHED['sto-3g']['N_pauli']}"
              f" Pauli terms at {GAUSSIAN_H2O_PUBLISHED['sto-3g']['Q']} qubits")
        print(f"  Wall time = {elapsed:.1f}s")
        print(f"{'='*60}")

    # --- Build results dict ---
    results = {
        'M': M,
        'Q': Q,
        'N_pauli': N_pauli,
        'M_core': M_core,
        'M_bond': M_bond,
        'M_bond_o': M_val_o,
        'M_bond_h': M_val_h,
        'M_lone': M_lone,
        'n_blocks': 7,
        'ERI_density_total': eri_density_total,
        'ERI_density_core': eri_density_core,
        'n_eri_core': n_eri_core,
        'n_eri_val_o': n_eri_val_o,
        'n_eri_val_h': n_eri_val_h,
        'n_eri_val_o_total': n_eri_val_o_total,
        'n_eri_val_h_total': n_eri_val_h_total,
        'n_eri_total': n_eri_total,
        'nuclear_repulsion': nuclear_repulsion,
        'V_OH': V_OH,
        'V_HH': V_HH,
        'V_NN': V_NN,
        'V_cross_nuc': V_cross,
        'E_core': E_core,
        'R_OH_bohr': R_OH,
        'R_HH_bohr': R_HH,
        'angle_HOH_deg': angle_HOH,
        'max_n_core': max_n_core,
        'max_n_val': max_n_val,
        'Z_O': Z_O,
        'Z_H': Z_H,
        'Z_eff_val': Z_eff_val,
        'wall_time_s': elapsed,
        'include_pk': include_pk,
        'pk_in_hamiltonian': pk_in_hamiltonian,
        'A_pk': A_pk if include_pk else None,
        'B_pk': B_pk if include_pk else None,
        'n_pk_nonzero': n_pk_nonzero,
        'approximations': [
            'Approach A: independent blocks (7), no cross-block ERIs',
            ('PK in quantum Hamiltonian' if pk_in_hamiltonian
             else 'PK computed separately (classical correction)')
            if include_pk else 'PK disabled',
            'PK on bond-H skipped (two-center, negligible at R~1.8 bohr)',
            'Cross-center nuclear attraction skipped',
            'Cross-block ERI = 0 (pseudopotential + Approach A)',
            f'E_core = {E_core:.4f} Ha (O2+ He-like variational estimate: -(Z-5/16)²)',
            f'Valence Z_eff = {Z_eff_val} (constant screening)',
            '4 O-side valence blocks share identical R^k integrals (computed once)',
        ],
    }

    # Save JSON report
    out_path = (Path(__file__).parent.parent
                / 'debug' / 'data' / 'composed_h2o_pauli_analysis.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    json_results = {k: v for k, v in results.items()
                    if not isinstance(v, np.ndarray)}
    with open(out_path, 'w') as f:
        json.dump(json_results, f, indent=2)
    if verbose:
        print(f"\nResults saved to {out_path}")

    # Attach non-serializable objects
    results['h1'] = h1
    results['h1_pk'] = h1_pk_full
    results['eri'] = eri
    results['qubit_op'] = qubit_op
    results['fermion_op'] = fermion_op
    results['states_core'] = states_core
    results['states_val_o'] = states_val_o
    results['states_val_h'] = states_val_h

    return results


# ---------------------------------------------------------------------------
# Multi-molecule scaling comparison (LiH / BeH2 / H2O / Gaussian)
# ---------------------------------------------------------------------------

def _fit_power_law_general(data: List[Dict]) -> Dict[str, float]:
    """Fit N_pauli = a * Q^alpha to a list of dicts with 'Q' and 'N_pauli'."""
    Q_arr = np.array([d['Q'] for d in data], dtype=float)
    N_arr = np.array([d['N_pauli'] for d in data], dtype=float)
    if len(Q_arr) < 2 or not np.all(Q_arr > 0) or not np.all(N_arr > 0):
        return {'alpha': np.nan, 'a': np.nan, 'R_squared': np.nan}
    log_Q = np.log(Q_arr)
    log_N = np.log(N_arr)
    coeffs = np.polyfit(log_Q, log_N, 1)
    alpha = float(coeffs[0])
    a = float(np.exp(coeffs[1]))
    log_N_pred = np.polyval(coeffs, log_Q)
    ss_res = np.sum((log_N - log_N_pred)**2)
    ss_tot = np.sum((log_N - np.mean(log_N))**2)
    r_sq = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 1.0
    return {'alpha': alpha, 'a': a, 'R_squared': r_sq}


def fit_gaussian_h2o_published_exponent() -> Dict[str, Any]:
    """
    Fit a power law N_pauli = C * Q^alpha to the 3 published H2O data points
    from Trenev et al. (Quantum 2025, Table 5, JW with 2-qubit reduction).
    """
    Q_arr = np.array([v['Q'] for v in GAUSSIAN_H2O_PUBLISHED.values()], dtype=float)
    N_arr = np.array([v['N_pauli'] for v in GAUSSIAN_H2O_PUBLISHED.values()], dtype=float)

    log_Q = np.log(Q_arr)
    log_N = np.log(N_arr)

    coeffs = np.polyfit(log_Q, log_N, 1)
    alpha = float(coeffs[0])
    C = float(np.exp(coeffs[1]))

    log_N_pred = np.polyval(coeffs, log_Q)
    ss_res = float(np.sum((log_N - log_N_pred)**2))
    ss_tot = float(np.sum((log_N - np.mean(log_N))**2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return {
        'alpha': alpha,
        'C': C,
        'R_squared': r_squared,
        'formula': f'N_pauli = {C:.4f} * Q^{alpha:.4f}',
        'n_points': len(Q_arr),
        'source': _TRENEV_REFERENCE,
        'note': 'JW encoding with 2-qubit reduction',
    }


def compare_all_molecules_scaling(
    max_n_values: Optional[List[int]] = None,
    verbose: bool = True,
    timeout_minutes: float = 20.0,
) -> Dict[str, Any]:
    """
    Run LiH, BeH2, and H2O at multiple basis sizes and compare Pauli scaling.

    Fits power law for each molecule and produces a combined comparison table
    including published Gaussian H2O data.

    Parameters
    ----------
    max_n_values : list of int or None
        Basis sizes to sweep. Default [1, 2, 3].
    verbose : bool
        Print comparison table.
    timeout_minutes : float
        If any single max_n takes longer than this, skip remaining points.

    Returns
    -------
    dict with sweep data, fitted exponents, and comparison.
    """
    if max_n_values is None:
        max_n_values = [1, 2, 3]

    lih_data: List[Dict[str, Any]] = []
    beh2_data: List[Dict[str, Any]] = []
    h2o_data: List[Dict[str, Any]] = []

    for n in max_n_values:
        t_start = time.perf_counter()

        if verbose:
            print(f"\n{'='*60}")
            print(f"  max_n = {n}: Building LiH...")
            print(f"{'='*60}")
        lih = build_composed_lih(max_n_core=n, max_n_val=n, verbose=verbose)
        lih_data.append({
            'max_n': n, 'M': lih['M'], 'Q': lih['Q'],
            'N_pauli': lih['N_pauli'],
            'ERI_density': lih['ERI_density_total'],
            'wall_time_s': lih['wall_time_s'],
        })

        if verbose:
            print(f"\n{'='*60}")
            print(f"  max_n = {n}: Building BeH2...")
            print(f"{'='*60}")
        beh2 = build_composed_beh2(max_n_core=n, max_n_val=n, verbose=verbose)
        beh2_data.append({
            'max_n': n, 'M': beh2['M'], 'Q': beh2['Q'],
            'N_pauli': beh2['N_pauli'],
            'ERI_density': beh2['ERI_density_total'],
            'wall_time_s': beh2['wall_time_s'],
        })

        if verbose:
            print(f"\n{'='*60}")
            print(f"  max_n = {n}: Building H2O...")
            print(f"{'='*60}")
        h2o = build_composed_h2o(max_n_core=n, max_n_val=n, verbose=verbose)
        h2o_data.append({
            'max_n': n, 'M': h2o['M'], 'Q': h2o['Q'],
            'N_pauli': h2o['N_pauli'],
            'ERI_density': h2o['ERI_density_total'],
            'wall_time_s': h2o['wall_time_s'],
        })

        elapsed_min = (time.perf_counter() - t_start) / 60.0
        if elapsed_min > timeout_minutes and n < max(max_n_values):
            if verbose:
                print(f"\n  WARNING: max_n={n} took {elapsed_min:.1f} min"
                      f" (limit={timeout_minutes:.0f} min). Stopping sweep.")
            break

    # Fit power laws
    lih_fit = _fit_power_law_general(lih_data)
    beh2_fit = _fit_power_law_general(beh2_data)
    h2o_fit = _fit_power_law_general(h2o_data)

    # Gaussian published fits
    gauss_lih_fit = fit_gaussian_lih_published_exponent()
    gauss_h2o_fit = fit_gaussian_h2o_published_exponent()

    if verbose:
        print(f"\n{'='*75}")
        print(f"  ALL-MOLECULE SCALING COMPARISON")
        print(f"{'='*75}")
        print(f"\n  {'System':>8} {'max_n':>5} {'M':>4} {'Q':>5} {'N_pauli':>10}"
              f" {'ERI_dens':>10} {'time':>8}")
        print(f"  {'-'*8} {'-'*5} {'-'*4} {'-'*5} {'-'*10}"
              f" {'-'*10} {'-'*8}")
        for i in range(len(lih_data)):
            for label, data in [('LiH', lih_data), ('BeH2', beh2_data), ('H2O', h2o_data)]:
                if i < len(data):
                    d = data[i]
                    print(f"  {label:>8} {d['max_n']:>5} {d['M']:>4} {d['Q']:>5}"
                          f" {d['N_pauli']:>10,} {d['ERI_density']:>10.2%}"
                          f" {d['wall_time_s']:>7.1f}s")
            print()

        print(f"  Fitted exponents (GeoVac composed):")
        print(f"    LiH:  N = {lih_fit['a']:.4f} x Q^{lih_fit['alpha']:.4f}"
              f"  (R² = {lih_fit['R_squared']:.4f})")
        print(f"    BeH2: N = {beh2_fit['a']:.4f} x Q^{beh2_fit['alpha']:.4f}"
              f"  (R² = {beh2_fit['R_squared']:.4f})")
        print(f"    H2O:  N = {h2o_fit['a']:.4f} x Q^{h2o_fit['alpha']:.4f}"
              f"  (R² = {h2o_fit['R_squared']:.4f})")

        print(f"\n  Gaussian published (Trenev et al.):")
        print(f"    LiH:  N = {gauss_lih_fit['C']:.4f} x Q^{gauss_lih_fit['alpha']:.4f}")
        print(f"    H2O:  N = {gauss_h2o_fit['C']:.4f} x Q^{gauss_h2o_fit['alpha']:.4f}")

        # Key question
        all_geovac_alphas = [lih_fit['alpha'], beh2_fit['alpha'], h2o_fit['alpha']]
        spread = max(all_geovac_alphas) - min(all_geovac_alphas)
        print(f"\n  Key question: do all molecules fall on the same ~Q^2.5 line?")
        print(f"    Exponent spread: {spread:.3f}")
        if spread < 0.5:
            print(f"    -> YES, consistent scaling across LiH/BeH2/H2O")
        else:
            print(f"    -> DIVERGENCE detected")

        # Published Gaussian H2O comparison
        print(f"\n  --- Published Gaussian H2O data ---")
        for basis, info in GAUSSIAN_H2O_PUBLISHED.items():
            print(f"    {basis:>10}: Q={info['Q']:>3}, N_pauli={info['N_pauli']:>7,}")

        print(f"{'='*75}")

    # Save results
    output = {
        'lih_data': lih_data,
        'beh2_data': beh2_data,
        'h2o_data': h2o_data,
        'lih_fit': lih_fit,
        'beh2_fit': beh2_fit,
        'h2o_fit': h2o_fit,
        'gaussian_lih_published': {k: v for k, v in GAUSSIAN_LIH_PUBLISHED.items()},
        'gaussian_h2o_published': {k: v for k, v in GAUSSIAN_H2O_PUBLISHED.items()},
        'gaussian_lih_fit': gauss_lih_fit,
        'gaussian_h2o_fit': gauss_h2o_fit,
        'reference': _TRENEV_REFERENCE,
    }

    out_path = (Path(__file__).parent.parent
                / 'debug' / 'data' / 'composed_h2o_pauli_analysis.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    if verbose:
        print(f"\nResults saved to {out_path}")

    return output


# ---------------------------------------------------------------------------
# H2 Bond-Pair Qubit Encoding (Track AZ)
# ---------------------------------------------------------------------------

def build_h2_bond_pair(
    max_n: int = 2,
    R: float = 1.4,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Build H2 qubit Hamiltonian using a single bond-pair block encoding.

    H2 has no core electrons — it is a single bond pair of 2 electrons in
    hydrogenic orbitals at Z_eff=1 (symmetric H-H). This is structurally
    identical to the He atomic encoding but with Z=1 instead of Z=2, plus
    nuclear repulsion V_NN = 1/R.

    The encoding uses hydrogenic (n,l,m) states on a single center with
    Z=1, which captures the correct angular momentum structure and Gaunt
    integral selection rules. Two-center effects (cross-nuclear attraction,
    prolate spheroidal splitting) are not included — this is the composed
    pipeline's single-block approximation, not the exact Level 2 solver.

    Parameters
    ----------
    max_n : int
        Maximum principal quantum number (default 2).
    R : float
        Internuclear distance in bohr (default 1.4, experimental H2 R_eq).
    verbose : bool
        Print progress and results.

    Returns
    -------
    dict with keys:
        M, Q, N_pauli, h1, eri, nuclear_repulsion, qubit_op, fermion_op,
        states, ERI_density, n_eri, wall_time_s
    """
    t0 = time.perf_counter()

    Z_eff = 1.0  # Each electron sees Z=1 in the bond pair
    Z_A = 1      # H nuclear charge
    Z_B = 1      # H nuclear charge

    # --- Define orbital basis ---
    states = _enumerate_states(max_n)
    M = len(states)
    Q = 2 * M  # spin-orbitals = qubits under JW

    if verbose:
        print(f"[build_h2_bond_pair] Orbital basis:")
        print(f"  Bond pair (Z_eff={Z_eff}): {M} orbitals, max_n={max_n}")
        print(f"  Qubits Q = {Q}")
        print(f"  R = {R:.4f} bohr")

    # --- Construct h1 (M x M) ---
    # Diagonal: exact hydrogenic eigenvalues at Z_eff=1
    h1 = np.zeros((M, M))
    for i, (n, l, m) in enumerate(states):
        h1[i, i] = -Z_eff**2 / (2.0 * n**2)

    if verbose:
        print(f"[build_h2_bond_pair] h1 constructed ({M}x{M}), diagonal only")

    # --- Construct ERI in physicist notation, then convert to chemist ---
    if verbose:
        print(f"[build_h2_bond_pair] Computing R^k integrals (Z={Z_eff})...")
    rk_cache = _compute_rk_integrals_block(Z_eff, states)
    eri_phys = _build_eri_block(Z_eff, states, rk_cache)
    n_eri = len(eri_phys)

    # Convert to chemist notation dense tensor
    eri = _physicist_to_chemist(eri_phys, M, offset=0)

    # Symmetrize: enforce (pq|rs) = (rs|pq)
    eri = 0.5 * (eri + eri.transpose(2, 3, 0, 1))

    n_eri_total = int(np.count_nonzero(np.abs(eri) > 1e-15))
    eri_density = n_eri_total / max(1, M**4)

    if verbose:
        print(f"[build_h2_bond_pair] ERI statistics:")
        print(f"  Physicist notation: {n_eri} nonzero")
        print(f"  Chemist tensor:    {n_eri_total} nonzero / {M**4} = {eri_density:.1%}")

    # --- Nuclear repulsion ---
    V_NN = Z_A * Z_B / R

    if verbose:
        print(f"[build_h2_bond_pair] V_NN = {V_NN:.6f} Ha")

    # --- JW transform ---
    if verbose:
        print(f"[build_h2_bond_pair] Building fermion operator...")
    fermion_op = build_fermion_op_from_integrals(h1, eri, V_NN)

    if verbose:
        print(f"[build_h2_bond_pair] Jordan-Wigner transform...")
    qubit_op = jordan_wigner(fermion_op)
    N_pauli = len(qubit_op.terms)

    elapsed = time.perf_counter() - t0

    if verbose:
        print(f"\n{'='*60}")
        print(f"H2 BOND-PAIR QUBIT HAMILTONIAN -- RESULTS")
        print(f"{'='*60}")
        print(f"  Spatial orbitals M = {M}")
        print(f"  Qubits           Q = {Q}")
        print(f"  Pauli terms      N = {N_pauli:,}")
        print(f"  ERI density          = {eri_density:.2%}")
        print(f"  R                    = {R:.4f} bohr")
        print(f"  V_NN                 = {V_NN:.6f} Ha")
        print(f"  Wall time = {elapsed:.1f}s")
        print(f"{'='*60}")

    results: Dict[str, Any] = {
        'M': M,
        'Q': Q,
        'N_pauli': N_pauli,
        'ERI_density': eri_density,
        'n_eri': n_eri,
        'n_eri_total': n_eri_total,
        'nuclear_repulsion': V_NN,
        'R_bohr': R,
        'max_n': max_n,
        'Z_eff': Z_eff,
        'wall_time_s': elapsed,
        'h1': h1,
        'eri': eri,
        'qubit_op': qubit_op,
        'fermion_op': fermion_op,
        'states': states,
    }

    return results


# ---------------------------------------------------------------------------
# New molecule wrapper functions (Track BJ)
# ---------------------------------------------------------------------------

def build_composed_hf(
    max_n_core: int = 2,
    max_n_val: int = 2,
    R: float = 1.733,
    E_core: Optional[float] = None,
    include_pk: bool = True,
    pk_in_hamiltonian: Optional[bool] = None,
    A_pk: Optional[float] = None,
    B_pk: Optional[float] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Build the composed HF qubit Hamiltonian and count Pauli terms.

    Orbital basis (5 blocks):
        Core:    F hydrogenic (Z=9) up to max_n_core
        Bond:    screened F (Z_eff=7) + H (Z=1), up to max_n_val
        Lone 1-3: screened F (Z_eff=7), up to max_n_val

    Parameters
    ----------
    max_n_core : int
        Maximum n for core orbitals on F (Z=9).
    max_n_val : int
        Maximum n for valence orbitals.
    R : float
        F-H bond distance in bohr (default: 1.733, expt 0.917 Å).
    E_core : float or None
        Core energy in Ha.
    include_pk : bool
        If True (default), compute PK on F-side valence blocks.
    pk_in_hamiltonian : bool or None
        If True, add PK to h1. If None, uses include_pk value.
    A_pk, B_pk : float or None
        PK parameters.
    verbose : bool
        Print progress and results.

    Returns
    -------
    dict
        Keys include M, Q, N_pauli, h1, h1_pk, eri, nuclear_repulsion,
        qubit_op, fermion_op, wall_time_s, blocks, spec_name.
    """
    spec = hf_spec(max_n_core, max_n_val, R, E_core, include_pk, A_pk, B_pk)
    if pk_in_hamiltonian is None:
        pk_in_hamiltonian = include_pk
    return build_composed_hamiltonian(
        spec, pk_in_hamiltonian=pk_in_hamiltonian, verbose=verbose,
    )


def build_composed_nh3(
    max_n_core: int = 2,
    max_n_val: int = 2,
    R_NH: float = 1.912,
    angle_HNH: float = 106.7,
    E_core: Optional[float] = None,
    include_pk: bool = True,
    pk_in_hamiltonian: Optional[bool] = None,
    A_pk: Optional[float] = None,
    B_pk: Optional[float] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Build the composed NH3 qubit Hamiltonian and count Pauli terms.

    Orbital basis (5 blocks):
        Core:    N hydrogenic (Z=7) up to max_n_core
        Bond 1-3: screened N (Z_eff=5) + H (Z=1), up to max_n_val
        Lone:    screened N (Z_eff=5), up to max_n_val

    Parameters
    ----------
    max_n_core : int
        Maximum n for core orbitals on N (Z=7).
    max_n_val : int
        Maximum n for valence orbitals.
    R_NH : float
        N-H bond distance in bohr (default: 1.912, expt 1.012 Å).
    angle_HNH : float
        H-N-H bond angle in degrees (default: 106.7°).
    E_core : float or None
        Core energy in Ha.
    include_pk : bool
        If True (default), compute PK on N-side valence blocks.
    pk_in_hamiltonian : bool or None
        If True, add PK to h1. If None, uses include_pk value.
    A_pk, B_pk : float or None
        PK parameters.
    verbose : bool
        Print progress and results.

    Returns
    -------
    dict
        Keys include M, Q, N_pauli, h1, h1_pk, eri, nuclear_repulsion,
        qubit_op, fermion_op, wall_time_s, blocks, spec_name.
    """
    spec = nh3_spec(
        max_n_core, max_n_val, R_NH, angle_HNH, E_core, include_pk, A_pk, B_pk,
    )
    if pk_in_hamiltonian is None:
        pk_in_hamiltonian = include_pk
    return build_composed_hamiltonian(
        spec, pk_in_hamiltonian=pk_in_hamiltonian, verbose=verbose,
    )


def build_composed_ch4(
    max_n_core: int = 2,
    max_n_val: int = 2,
    R_CH: float = 2.050,
    E_core: Optional[float] = None,
    include_pk: bool = True,
    pk_in_hamiltonian: Optional[bool] = None,
    A_pk: Optional[float] = None,
    B_pk: Optional[float] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Build the composed CH4 qubit Hamiltonian and count Pauli terms.

    Orbital basis (5 blocks):
        Core:    C hydrogenic (Z=6) up to max_n_core
        Bond 1-4: screened C (Z_eff=4) + H (Z=1), up to max_n_val

    Parameters
    ----------
    max_n_core : int
        Maximum n for core orbitals on C (Z=6).
    max_n_val : int
        Maximum n for valence orbitals.
    R_CH : float
        C-H bond distance in bohr (default: 2.050, expt 1.085 Å).
    E_core : float or None
        Core energy in Ha.
    include_pk : bool
        If True (default), compute PK on C-side valence blocks.
    pk_in_hamiltonian : bool or None
        If True, add PK to h1. If None, uses include_pk value.
    A_pk, B_pk : float or None
        PK parameters.
    verbose : bool
        Print progress and results.

    Returns
    -------
    dict
        Keys include M, Q, N_pauli, h1, h1_pk, eri, nuclear_repulsion,
        qubit_op, fermion_op, wall_time_s, blocks, spec_name.
    """
    spec = ch4_spec(max_n_core, max_n_val, R_CH, E_core, include_pk, A_pk, B_pk)
    if pk_in_hamiltonian is None:
        pk_in_hamiltonian = include_pk
    return build_composed_hamiltonian(
        spec, pk_in_hamiltonian=pk_in_hamiltonian, verbose=verbose,
    )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import sys
    if '--sweep' in sys.argv:
        print("Running composed LiH scaling sweep...\n")
        composed_lih_scaling_sweep()
    elif '--analysis' in sys.argv:
        print("Running composed LiH full analysis...\n")
        composed_lih_full_analysis()
    elif '--beh2' in sys.argv:
        print("Running composed BeH2 qubit analysis...\n")
        result = build_composed_beh2()
    elif '--compare' in sys.argv:
        print("Running LiH vs BeH2 scaling comparison...\n")
        compare_lih_beh2_scaling()
    elif '--h2o' in sys.argv:
        print("Running composed H2O qubit analysis...\n")
        result = build_composed_h2o()
    elif '--compare-all' in sys.argv:
        print("Running LiH / BeH2 / H2O scaling comparison...\n")
        compare_all_molecules_scaling()
    else:
        print("Running composed LiH qubit analysis at default parameters...\n")
        result = build_composed_lih()
