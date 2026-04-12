"""
Deuteron Hamiltonian Builder — Nuclear Qubit Hamiltonian via Jordan-Wigner.

Builds the nuclear Hamiltonian for deuterium (1 proton + 1 neutron) in the
harmonic oscillator shell model basis with the Minnesota NN potential, encodes
it via Jordan-Wigner transformation, and exports to Pauli representation.

Uses the uncoupled (n_r, l, m_l, m_s) basis for computing two-body matrix
elements. The spin-singlet and spin-triplet projectors of the Minnesota
potential are handled by explicit CG decomposition of the two-nucleon spin.

Tensor-product JW approach: separate JW encodings for proton and neutron
registers, combined via tensor product.

Units: MeV throughout.

Author: GeoVac Development Team
Date: April 2026
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from math import factorial, sqrt as msqrt

from geovac.nuclear.minnesota import (
    ho_length_parameter,
    minnesota_matrix_element_analytical,
    minnesota_params,
)
from geovac.nuclear.moshinsky import (
    moshinsky_bracket,
    lab_to_relative_matrix_element,
)
from geovac.nuclear.spin_orbit import ls_eigenvalue


# ---------------------------------------------------------------------------
# Single-particle state in (n_r, l, m_l, m_s) basis
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SPState:
    """Single-particle state in uncoupled (n_r, l, m_l, m_s) basis."""
    n_r: int       # radial quantum number
    l: int         # orbital angular momentum
    m_l: int       # orbital magnetic quantum number
    m_s: float     # spin projection (+0.5 or -0.5)
    species: str   # 'proton' or 'neutron'

    @property
    def N(self) -> int:
        """HO major shell number."""
        return 2 * self.n_r + self.l

    def __repr__(self) -> str:
        letters = 'spdfghijklmno'
        l_str = letters[self.l] if self.l < len(letters) else f'[{self.l}]'
        spin_str = 'u' if self.m_s > 0 else 'd'
        return f'{self.n_r}{l_str}(m_l={self.m_l},{spin_str},{self.species[0]})'


# ---------------------------------------------------------------------------
# DeuteronSpec
# ---------------------------------------------------------------------------

@dataclass
class DeuteronSpec:
    """Specification for the deuteron system."""
    N_shells: int = 2         # number of HO major shells (N=0..N_shells-1)
    hw: float = 10.0          # HO frequency in MeV

    @property
    def b(self) -> float:
        """HO length parameter in fm."""
        return ho_length_parameter(self.hw)


# ---------------------------------------------------------------------------
# Enumerate single-particle states
# ---------------------------------------------------------------------------

def enumerate_sp_states(spec: DeuteronSpec) -> Tuple[List[SPState], List[SPState]]:
    """
    Enumerate single-particle states in uncoupled basis for protons and neutrons.

    For N_shells=2: N=0 (0s, m_l=0, m_s=+/-1/2) + N=1 (0p, m_l=-1,0,1, m_s=+/-1/2)
    = 2 + 6 = 8 states per species.

    Returns (proton_states, neutron_states).
    """
    def _make_states(species: str) -> List[SPState]:
        states = []
        for N in range(spec.N_shells):
            for l in range(N, -1, -2):
                n_r = (N - l) // 2
                for m_l in range(-l, l + 1):
                    for m_s in [0.5, -0.5]:
                        states.append(SPState(n_r=n_r, l=l, m_l=m_l,
                                              m_s=m_s, species=species))
        return states

    return _make_states('proton'), _make_states('neutron')


# ---------------------------------------------------------------------------
# One-body Hamiltonian
# ---------------------------------------------------------------------------

def build_one_body(
    states: List[SPState],
    spec: DeuteronSpec,
) -> np.ndarray:
    """
    Build the one-body Hamiltonian matrix (diagonal: HO energies).

    H1[i,j] = delta_{ij} * hw * (N_i + 3/2)
    """
    n = len(states)
    h1 = np.zeros((n, n))
    for i, s in enumerate(states):
        h1[i, i] = spec.hw * (s.N + 1.5)
    return h1


# ---------------------------------------------------------------------------
# Two-body matrix elements for central, spin-dependent potential
# ---------------------------------------------------------------------------

def _spin_cg(ms1: float, ms2: float, S: int, MS: float) -> float:
    """
    Clebsch-Gordan coefficient <1/2 ms1, 1/2 ms2 | S MS> for two spin-1/2.

    Known analytically:
    S=0, MS=0: <+1/2,-1/2|0,0> = 1/sqrt(2), <-1/2,+1/2|0,0> = -1/sqrt(2)
    S=1, MS=+1: <+1/2,+1/2|1,1> = 1
    S=1, MS=0: <+1/2,-1/2|1,0> = 1/sqrt(2), <-1/2,+1/2|1,0> = 1/sqrt(2)
    S=1, MS=-1: <-1/2,-1/2|1,-1> = 1
    """
    if abs(ms1 + ms2 - MS) > 1e-10:
        return 0.0
    if S == 0 and abs(MS) < 1e-10:
        if abs(ms1 - 0.5) < 1e-10 and abs(ms2 + 0.5) < 1e-10:
            return 1.0 / msqrt(2.0)
        elif abs(ms1 + 0.5) < 1e-10 and abs(ms2 - 0.5) < 1e-10:
            return -1.0 / msqrt(2.0)
        return 0.0
    elif S == 1:
        if abs(MS - 1.0) < 1e-10:
            return 1.0 if (abs(ms1 - 0.5) < 1e-10 and abs(ms2 - 0.5) < 1e-10) else 0.0
        elif abs(MS + 1.0) < 1e-10:
            return 1.0 if (abs(ms1 + 0.5) < 1e-10 and abs(ms2 + 0.5) < 1e-10) else 0.0
        elif abs(MS) < 1e-10:
            if abs(ms1 - 0.5) < 1e-10 and abs(ms2 + 0.5) < 1e-10:
                return 1.0 / msqrt(2.0)
            elif abs(ms1 + 0.5) < 1e-10 and abs(ms2 - 0.5) < 1e-10:
                return 1.0 / msqrt(2.0)
            return 0.0
    return 0.0


def _orbital_cg(l1: int, m1: int, l2: int, m2: int, L: int, ML: int) -> float:
    """CG coefficient for integer angular momenta."""
    if m1 + m2 != ML:
        return 0.0
    from geovac.angular_integrals import wigner3j
    return (-1)**(l1 - l2 + ML) * msqrt(2*L + 1) * wigner3j(l1, l2, L, m1, m2, -ML)


def compute_tbme(
    states_p: List[SPState],
    states_n: List[SPState],
    spec: DeuteronSpec,
) -> Dict[Tuple[int, int, int, int], float]:
    """
    Compute proton-neutron two-body matrix elements.

    V_{ijkl} = <p_i n_j | V_NN | p_k n_l>

    For Minnesota (central, spin-dependent):
    <p_i n_j | V | p_k n_l>
    = sum_{L,S} CG(li,mi; lj,mj|L,ML) CG(lk,mk; ll,ml|L,ML)
              * CG(si,msi; sj,msj|S,MS) CG(sk,msk; sl,msl|S,MS)
              * <n_ri li, n_rj lj; L,S | V | n_rk lk, n_rl ll; L,S>_Moshinsky

    The last factor is computed via Moshinsky transformation to relative coordinates.

    Parameters
    ----------
    states_p : list of SPState (proton)
    states_n : list of SPState (neutron)
    spec : DeuteronSpec

    Returns
    -------
    dict (i_p, j_n, k_p, l_n) -> float (MeV)
    """
    b = spec.b
    tbme: Dict[Tuple[int, int, int, int], float] = {}

    n_p = len(states_p)
    n_n = len(states_n)

    # Define V_rel function for Moshinsky
    def V_rel(n_rel: int, l_rel: int, n_rel_p: int, l_rel_p: int,
              S: int, b_val: float) -> float:
        return minnesota_matrix_element_analytical(
            n_rel, l_rel, n_rel_p, l_rel_p, S, b_val)

    # Cache for lab-to-relative MEs: (n1,l1,n2,l2, n3,l3,n4,l4, L, S) -> float
    lab_rel_cache: Dict[Tuple, float] = {}

    for i in range(n_p):
        pi = states_p[i]
        for j in range(n_n):
            nj = states_n[j]
            for k in range(n_p):
                pk = states_p[k]
                for l in range(n_n):
                    nl = states_n[l]

                    # m_l conservation
                    ML = pi.m_l + nj.m_l
                    if pi.m_l + nj.m_l != pk.m_l + nl.m_l:
                        continue
                    # m_s conservation
                    MS = pi.m_s + nj.m_s
                    if abs(pi.m_s + nj.m_s - pk.m_s - nl.m_s) > 1e-10:
                        continue

                    me = 0.0

                    # Sum over L (total orbital ang. mom.)
                    L_min = max(abs(pi.l - nj.l), abs(pk.l - nl.l))
                    L_max = min(pi.l + nj.l, pk.l + nl.l)

                    for L in range(L_min, L_max + 1):
                        if abs(ML) > L:
                            continue

                        # Orbital CG coefficients
                        cg_orb_bra = _orbital_cg(pi.l, pi.m_l, nj.l, nj.m_l, L, ML)
                        if abs(cg_orb_bra) < 1e-15:
                            continue
                        cg_orb_ket = _orbital_cg(pk.l, pk.m_l, nl.l, nl.m_l, L, ML)
                        if abs(cg_orb_ket) < 1e-15:
                            continue

                        # Sum over S (total spin)
                        for S in [0, 1]:
                            if abs(MS) > S:
                                continue

                            # Spin CG coefficients
                            cg_spin_bra = _spin_cg(pi.m_s, nj.m_s, S, MS)
                            if abs(cg_spin_bra) < 1e-15:
                                continue
                            cg_spin_ket = _spin_cg(pk.m_s, nl.m_s, S, MS)
                            if abs(cg_spin_ket) < 1e-15:
                                continue

                            # Lab-to-relative ME via Moshinsky
                            cache_key = (pi.n_r, pi.l, nj.n_r, nj.l,
                                         pk.n_r, pk.l, nl.n_r, nl.l, L, S)
                            if cache_key not in lab_rel_cache:
                                lab_rel_cache[cache_key] = lab_to_relative_matrix_element(
                                    pi.n_r, pi.l, nj.n_r, nj.l,
                                    pk.n_r, pk.l, nl.n_r, nl.l,
                                    L, S, V_rel, b)
                            v_ls = lab_rel_cache[cache_key]

                            me += cg_orb_bra * cg_orb_ket * cg_spin_bra * cg_spin_ket * v_ls

                    if abs(me) > 1e-15:
                        tbme[(i, j, k, l)] = me

    return tbme


# ---------------------------------------------------------------------------
# Jordan-Wigner encoding (tensor product approach)
# ---------------------------------------------------------------------------

def _number_op_pauli(q: int, n_qubits: int) -> Dict[str, float]:
    """Number operator a^dag_q a_q = (I - Z_q)/2."""
    I_str = 'I' * n_qubits
    Z_str = list(I_str)
    Z_str[q] = 'Z'
    Z_str = ''.join(Z_str)
    return {I_str: 0.5, Z_str: -0.5}


def _one_body_pauli(p: int, q: int, n_qubits: int) -> Dict[str, float]:
    """
    JW encoding of a^dag_p a_q on a single-species register.

    For p == q: number operator (I - Z)/2.
    For p != q: (XX + YY)/4 * Z-string (real part only; imaginary parts
    cancel in Hermitian Hamiltonian when summing V_{ijkl} + V_{klij}).
    """
    if p == q:
        return _number_op_pauli(p, n_qubits)

    mn, mx = min(p, q), max(p, q)
    result: Dict[str, float] = {}

    for pauli_pair, coeff in [('XX', 0.25), ('YY', 0.25)]:
        s = ['I'] * n_qubits
        s[mn] = pauli_pair[0]
        s[mx] = pauli_pair[1]
        for z in range(mn + 1, mx):
            s[z] = 'Z'
        key = ''.join(s)
        result[key] = result.get(key, 0.0) + coeff

    # For p < q: add imaginary part XY - YX
    # For p > q: subtract (conjugate)
    sign = 1.0 if p < q else -1.0
    for pauli_pair, coeff in [('XY', sign * 0.25), ('YX', -sign * 0.25)]:
        s = ['I'] * n_qubits
        s[mn] = pauli_pair[0]
        s[mx] = pauli_pair[1]
        for z in range(mn + 1, mx):
            s[z] = 'Z'
        key = ''.join(s)
        result[key] = result.get(key, 0.0) + coeff

    return result


def _tensor_product_pauli(
    pauli_p: Dict[str, float],
    pauli_n: Dict[str, float],
) -> Dict[str, float]:
    """Tensor product of Pauli dicts from proton and neutron registers."""
    result: Dict[str, float] = {}
    for sp, cp in pauli_p.items():
        for sn, cn in pauli_n.items():
            key = sp + sn
            coeff = cp * cn
            if abs(coeff) > 1e-15:
                result[key] = result.get(key, 0.0) + coeff
    return result


def _add_pauli(
    target: Dict[str, float],
    source: Dict[str, float],
    scale: float = 1.0,
) -> None:
    """Add scaled Pauli terms from source to target dict."""
    for key, val in source.items():
        target[key] = target.get(key, 0.0) + scale * val


def _embed_single_species(
    pauli_dict: Dict[str, float],
    species: str,
    n_p: int,
    n_n: int,
) -> Dict[str, float]:
    """Embed a single-species Pauli dict into the full register."""
    result: Dict[str, float] = {}
    for pstr, coeff in pauli_dict.items():
        if species == 'proton':
            full_str = pstr + 'I' * n_n
        else:
            full_str = 'I' * n_p + pstr
        result[full_str] = result.get(full_str, 0.0) + coeff
    return result


# ---------------------------------------------------------------------------
# Build deuteron Hamiltonian
# ---------------------------------------------------------------------------

def build_deuteron_hamiltonian(
    N_shells: int = 2,
    hw: float = 10.0,
) -> Dict[str, Any]:
    """
    Build the deuteron qubit Hamiltonian.

    Parameters
    ----------
    N_shells : int
        Number of HO shells.
    hw : float
        HO frequency in MeV.

    Returns
    -------
    dict with keys: H_pauli, H_matrix, states_p, states_n, Q, Q_p, Q_n, metadata
    """
    spec = DeuteronSpec(N_shells=N_shells, hw=hw)
    states_p, states_n = enumerate_sp_states(spec)
    n_p = len(states_p)
    n_n = len(states_n)
    Q = n_p + n_n

    h1_p = build_one_body(states_p, spec)
    h1_n = build_one_body(states_n, spec)

    # Two-body matrix elements
    tbme = compute_tbme(states_p, states_n, spec)

    # Build Pauli Hamiltonian using tensor-product JW
    H_pauli: Dict[str, float] = {}

    # One-body: proton
    for i in range(n_p):
        if abs(h1_p[i, i]) > 1e-15:
            pauli_p = _number_op_pauli(i, n_p)
            embedded = _embed_single_species(pauli_p, 'proton', n_p, n_n)
            _add_pauli(H_pauli, embedded, h1_p[i, i])

    # One-body: neutron
    for j in range(n_n):
        if abs(h1_n[j, j]) > 1e-15:
            pauli_n = _number_op_pauli(j, n_n)
            embedded = _embed_single_species(pauli_n, 'neutron', n_p, n_n)
            _add_pauli(H_pauli, embedded, h1_n[j, j])

    # Two-body: V_{ijkl} a^dag_{p_i} a^dag_{n_j} a_{n_l} a_{p_k}
    # = V_{ijkl} (a^dag_{p_i} a_{p_k}) (a^dag_{n_j} a_{n_l})
    for (i, j, k, l), v in tbme.items():
        if abs(v) < 1e-15:
            continue
        pauli_p = _one_body_pauli(i, k, n_p)
        pauli_n = _one_body_pauli(j, l, n_n)
        pn_pauli = _tensor_product_pauli(pauli_p, pauli_n)
        _add_pauli(H_pauli, pn_pauli, v)

    # Clean up small terms
    H_pauli = {k: v for k, v in H_pauli.items() if abs(v) > 1e-12}

    # Build full matrix for diagonalization
    H_matrix = _build_fci_matrix(states_p, states_n, h1_p, h1_n, tbme)

    metadata = {
        'system': 'deuteron',
        'N_shells': N_shells,
        'hw_MeV': hw,
        'b_fm': spec.b,
        'energy_unit': 'MeV',
        'Q': Q, 'Q_p': n_p, 'Q_n': n_n,
    }

    return {
        'H_pauli': H_pauli,
        'H_matrix': H_matrix,
        'states_p': states_p,
        'states_n': states_n,
        'Q': Q, 'Q_p': n_p, 'Q_n': n_n,
        'metadata': metadata,
    }


# ---------------------------------------------------------------------------
# FCI matrix in the 1p+1n Hilbert space
# ---------------------------------------------------------------------------

def _build_fci_matrix(
    states_p: List[SPState],
    states_n: List[SPState],
    h1_p: np.ndarray,
    h1_n: np.ndarray,
    tbme: Dict[Tuple[int, int, int, int], float],
) -> np.ndarray:
    """
    Build the FCI matrix in (1 proton) x (1 neutron) space.
    Basis: |p_i, n_j>. Dimension: n_p * n_n.
    """
    n_p = len(states_p)
    n_n = len(states_n)
    dim = n_p * n_n
    H = np.zeros((dim, dim))

    for i in range(n_p):
        for j in range(n_n):
            idx_bra = i * n_n + j
            H[idx_bra, idx_bra] = h1_p[i, i] + h1_n[j, j]
            for k in range(n_p):
                for l in range(n_n):
                    idx_ket = k * n_n + l
                    v = tbme.get((i, j, k, l), 0.0)
                    if abs(v) > 1e-15:
                        H[idx_bra, idx_ket] += v

    return H


# ---------------------------------------------------------------------------
# Diagonalize and analyze
# ---------------------------------------------------------------------------

def diagonalize_deuteron(
    N_shells: int = 2,
    hw: float = 10.0,
) -> Dict[str, Any]:
    """Build and diagonalize the deuteron Hamiltonian."""
    data = build_deuteron_hamiltonian(N_shells=N_shells, hw=hw)
    H = data['H_matrix']
    evals, evecs = np.linalg.eigh(H)
    E_gs = evals[0]
    E_zpe = 3.0 * hw  # 3/2 * hw * 2 particles

    return {
        'E_gs': E_gs,
        'E_binding': E_gs,
        'eigenvalues': evals,
        'ground_state': evecs[:, 0],
        'H_data': data,
        'E_zpe': E_zpe,
    }


def analyze_deuteron_hamiltonian(
    H_pauli: Dict[str, float],
    Q: int,
) -> Dict[str, Any]:
    """Analyze the Pauli Hamiltonian."""
    identity_key = 'I' * Q
    non_identity = {k: v for k, v in H_pauli.items() if k != identity_key}

    n_terms = len(non_identity)
    one_norm = sum(abs(v) for v in non_identity.values())
    one_norm_total = sum(abs(v) for v in H_pauli.values())

    qwc_groups = _compute_qwc_groups(non_identity)

    n_z_only = sum(1 for k in non_identity if all(c in 'IZ' for c in k))
    n_with_xy = n_terms - n_z_only

    return {
        'n_pauli_terms': n_terms,
        'n_pauli_terms_total': len(H_pauli),
        'one_norm_ni': one_norm,
        'one_norm_total': one_norm_total,
        'n_qwc_groups': len(qwc_groups),
        'n_z_only_terms': n_z_only,
        'n_xy_terms': n_with_xy,
        'identity_coeff': H_pauli.get(identity_key, 0.0),
        'Q': Q,
    }


def _compute_qwc_groups(pauli_dict: Dict[str, float]) -> List[List[str]]:
    """Compute qubitwise commuting groups via greedy coloring."""
    terms = list(pauli_dict.keys())
    if not terms:
        return []

    def _qwc(a: str, b: str) -> bool:
        for ca, cb in zip(a, b):
            if ca == 'I' or cb == 'I':
                continue
            if ca != cb:
                return False
        return True

    groups: List[List[str]] = []
    for term in terms:
        placed = False
        for group in groups:
            if all(_qwc(term, other) for other in group):
                group.append(term)
                placed = True
                break
        if not placed:
            groups.append([term])
    return groups


def hw_scan(
    N_shells: int = 2,
    hw_values: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """Scan hw to find optimal oscillator frequency."""
    if hw_values is None:
        hw_values = list(np.arange(5.0, 27.5, 2.5))

    energies = []
    for hw in hw_values:
        try:
            result = diagonalize_deuteron(N_shells=N_shells, hw=hw)
            energies.append(result['E_gs'])
        except Exception:
            energies.append(float('nan'))

    energies = np.array(energies)
    valid = ~np.isnan(energies)
    if np.any(valid):
        idx_opt = int(np.argmin(energies[valid]))
        valid_hw = np.array(hw_values)[valid]
        return {
            'hw_values': hw_values,
            'energies': energies.tolist(),
            'optimal_hw': float(valid_hw[idx_opt]),
            'optimal_E': float(energies[valid][idx_opt]),
        }
    return {'hw_values': hw_values, 'energies': energies.tolist(),
            'optimal_hw': None, 'optimal_E': None}


# ===========================================================================
# He-4 (2 protons + 2 neutrons) — Track NF
# ===========================================================================

# Proton Coulomb constant: e^2 = 1.44 MeV·fm
_E2_MEV_FM = 1.44


def _coulomb_me_ho(
    n1: int, l1: int, n2: int, l2: int,
    n3: int, l3: int, n4: int, l4: int,
    L: int, b: float,
) -> float:
    """
    Coulomb matrix element <n1 l1, n2 l2; L | e^2/r_rel | n3 l3, n4 l4; L>
    via Moshinsky transformation to relative coordinates.

    The Coulomb interaction 1/r in the relative coordinate of two protons.
    Uses the HO basis in relative coordinates.

    Parameters
    ----------
    n1, l1, n2, l2 : int
        Bra lab-frame quantum numbers.
    n3, l3, n4, l4 : int
        Ket lab-frame quantum numbers.
    L : int
        Total orbital angular momentum.
    b : float
        HO length parameter in fm.

    Returns
    -------
    float
        Matrix element in MeV.
    """
    def V_coulomb_rel(n_rel: int, l_rel: int, n_rel_p: int, l_rel_p: int,
                      S: int, b_val: float) -> float:
        """<n' l' | e^2/r | n l> in relative HO basis (S-independent)."""
        if l_rel != l_rel_p:
            return 0.0
        l = l_rel
        # Use numerical integration
        r_max = 6.0 * b_val * np.sqrt(max(2*n_rel + l, 2*n_rel_p + l) + 3)
        n_grid = 4000
        r_grid = np.linspace(1e-10, r_max, n_grid)
        from geovac.nuclear.minnesota import _ho_radial_wf
        R_bra = _ho_radial_wf(n_rel_p, l, r_grid, b_val)
        R_ket = _ho_radial_wf(n_rel, l, r_grid, b_val)
        # Coulomb: e^2 / r_rel, but r_rel = r * sqrt(2) * b in relative coord
        # Actually, the relative coordinate uses b_rel = b (same oscillator length
        # for equal mass). The 1/r is just 1/r in fm.
        # Note: for equal-mass particles, the relative coordinate is
        # r_rel = (r1 - r2)/sqrt(2), so the actual distance is sqrt(2)*r_rel.
        # But the Moshinsky bracket already accounts for this coordinate change.
        # The Coulomb potential between the two protons is e^2/|r1-r2| = e^2/(sqrt(2)*r_rel).
        # However, the standard convention for Moshinsky brackets uses mass-weighted
        # Jacobi coordinates where the potential is expressed in terms of r_rel.
        # For equal masses: r_rel = (r1-r2)/sqrt(2), so |r1-r2| = sqrt(2)*r_rel.
        # The matrix element of e^2/|r1-r2| = (e^2/sqrt(2)) * <n'l'|1/r_rel|nl>
        V = _E2_MEV_FM / (np.sqrt(2.0) * r_grid)
        integrand = R_bra * V * R_ket * r_grid**2
        return float(np.trapezoid(integrand, r_grid))

    return lab_to_relative_matrix_element(
        n1, l1, n2, l2, n3, l3, n4, l4, L, 0, V_coulomb_rel, b)


def compute_same_species_tbme(
    states: List[SPState],
    spec: DeuteronSpec,
    include_coulomb: bool = False,
) -> Dict[Tuple[int, int, int, int], float]:
    """
    Compute antisymmetrized two-body matrix elements for same-species pairs.

    V_{ijkl}_AS = <i j | V | k l> - <i j | V | l k>  (direct - exchange)

    For identical fermions, antisymmetrization is required.

    Parameters
    ----------
    states : list of SPState
        Single-particle states for one species (all proton or all neutron).
    spec : DeuteronSpec
        System specification.
    include_coulomb : bool
        If True and species is 'proton', include Coulomb repulsion.

    Returns
    -------
    dict (i, j, k, l) -> float (MeV)
        Antisymmetrized matrix elements.
    """
    b = spec.b
    n = len(states)
    species = states[0].species if states else 'proton'
    use_coulomb = include_coulomb and species == 'proton'

    def V_rel_nuc(n_rel: int, l_rel: int, n_rel_p: int, l_rel_p: int,
                  S: int, b_val: float) -> float:
        return minnesota_matrix_element_analytical(
            n_rel, l_rel, n_rel_p, l_rel_p, S, b_val)

    # Coulomb in relative coords: e^2/(sqrt(2)*r), spin-independent
    _coulomb_rel_cache: Dict[Tuple, float] = {}
    def V_rel_coul(n_rel: int, l_rel: int, n_rel_p: int, l_rel_p: int,
                   S: int, b_val: float) -> float:
        if l_rel != l_rel_p:
            return 0.0
        key = (n_rel, l_rel, n_rel_p, l_rel_p)
        if key in _coulomb_rel_cache:
            return _coulomb_rel_cache[key]
        # <n'l|e^2/(sqrt(2)*r)|nl> in relative HO basis
        from geovac.nuclear.minnesota import _ho_radial_wf
        r_max = 6.0 * b_val * np.sqrt(
            max(2*n_rel + l_rel, 2*n_rel_p + l_rel_p) + 3)
        r_grid = np.linspace(1e-10, r_max, 4000)
        R_bra = _ho_radial_wf(n_rel_p, l_rel_p, r_grid, b_val)
        R_ket = _ho_radial_wf(n_rel, l_rel, r_grid, b_val)
        V = _E2_MEV_FM / (np.sqrt(2.0) * r_grid)
        val = float(np.trapezoid(R_bra * V * R_ket * r_grid**2, r_grid))
        _coulomb_rel_cache[key] = val
        return val

    def V_rel_total(n_rel, l_rel, n_rel_p, l_rel_p, S, b_val):
        me = V_rel_nuc(n_rel, l_rel, n_rel_p, l_rel_p, S, b_val)
        if use_coulomb:
            me += V_rel_coul(n_rel, l_rel, n_rel_p, l_rel_p, S, b_val)
        return me

    V_rel = V_rel_total

    # Cache for lab-to-relative MEs
    lab_rel_cache: Dict[Tuple, float] = {}

    tbme: Dict[Tuple[int, int, int, int], float] = {}

    for i in range(n):
        si = states[i]
        for j in range(n):
            sj = states[j]
            for k in range(n):
                sk = states[k]
                for l in range(n):
                    sl = states[l]

                    # Compute direct ME: <i j | V | k l>
                    direct = _compute_two_body_me(
                        si, sj, sk, sl, V_rel, b, lab_rel_cache)

                    # Compute exchange ME: <i j | V | l k>
                    exchange = _compute_two_body_me(
                        si, sj, sl, sk, V_rel, b, lab_rel_cache)

                    me_as = direct - exchange

                    if abs(me_as) > 1e-15:
                        tbme[(i, j, k, l)] = me_as

    return tbme


def _compute_two_body_me(
    s1: SPState, s2: SPState, s3: SPState, s4: SPState,
    V_rel_func, b: float,
    cache: Dict[Tuple, float],
) -> float:
    """
    Compute <s1 s2 | V | s3 s4> using CG decomposition and Moshinsky.

    This is the same logic as compute_tbme but generalized to any two states.
    """
    # Conservation laws
    ML = s1.m_l + s2.m_l
    if s1.m_l + s2.m_l != s3.m_l + s4.m_l:
        return 0.0
    MS = s1.m_s + s2.m_s
    if abs(s1.m_s + s2.m_s - s3.m_s - s4.m_s) > 1e-10:
        return 0.0

    me = 0.0
    L_min = max(abs(s1.l - s2.l), abs(s3.l - s4.l))
    L_max = min(s1.l + s2.l, s3.l + s4.l)

    for L in range(L_min, L_max + 1):
        if abs(ML) > L:
            continue
        cg_orb_bra = _orbital_cg(s1.l, s1.m_l, s2.l, s2.m_l, L, ML)
        if abs(cg_orb_bra) < 1e-15:
            continue
        cg_orb_ket = _orbital_cg(s3.l, s3.m_l, s4.l, s4.m_l, L, ML)
        if abs(cg_orb_ket) < 1e-15:
            continue

        for S in [0, 1]:
            if abs(MS) > S:
                continue
            cg_spin_bra = _spin_cg(s1.m_s, s2.m_s, S, MS)
            if abs(cg_spin_bra) < 1e-15:
                continue
            cg_spin_ket = _spin_cg(s3.m_s, s4.m_s, S, MS)
            if abs(cg_spin_ket) < 1e-15:
                continue

            cache_key = (s1.n_r, s1.l, s2.n_r, s2.l,
                         s3.n_r, s3.l, s4.n_r, s4.l, L, S)
            if cache_key not in cache:
                cache[cache_key] = lab_to_relative_matrix_element(
                    s1.n_r, s1.l, s2.n_r, s2.l,
                    s3.n_r, s3.l, s4.n_r, s4.l,
                    L, S, V_rel_func, b)
            v_ls = cache[cache_key]

            me += cg_orb_bra * cg_orb_ket * cg_spin_bra * cg_spin_ket * v_ls

    return me


def _compute_coulomb_me(
    s1: SPState, s2: SPState, s3: SPState, s4: SPState,
    b: float,
    cache: Dict[Tuple, float],
) -> float:
    """
    Compute <s1 s2 | e^2/r | s3 s4> Coulomb matrix element.

    Coulomb is spin-independent, so we only need orbital overlap for spin.
    """
    # Spin conservation: delta(ms1,ms3) * delta(ms2,ms4)
    if abs(s1.m_s - s3.m_s) > 1e-10 or abs(s2.m_s - s4.m_s) > 1e-10:
        return 0.0
    # m_l conservation
    if s1.m_l + s2.m_l != s3.m_l + s4.m_l:
        return 0.0

    ML = s1.m_l + s2.m_l
    me = 0.0
    L_min = max(abs(s1.l - s2.l), abs(s3.l - s4.l))
    L_max = min(s1.l + s2.l, s3.l + s4.l)

    for L in range(L_min, L_max + 1):
        if abs(ML) > L:
            continue
        cg_orb_bra = _orbital_cg(s1.l, s1.m_l, s2.l, s2.m_l, L, ML)
        if abs(cg_orb_bra) < 1e-15:
            continue
        cg_orb_ket = _orbital_cg(s3.l, s3.m_l, s4.l, s4.m_l, L, ML)
        if abs(cg_orb_ket) < 1e-15:
            continue

        cache_key = (s1.n_r, s1.l, s2.n_r, s2.l,
                     s3.n_r, s3.l, s4.n_r, s4.l, L)
        if cache_key not in cache:
            cache[cache_key] = _coulomb_me_ho(
                s1.n_r, s1.l, s2.n_r, s2.l,
                s3.n_r, s3.l, s4.n_r, s4.l,
                L, b)
        v_coul = cache[cache_key]

        me += cg_orb_bra * cg_orb_ket * v_coul

    return me


def compute_cross_species_tbme(
    states_a: List[SPState],
    states_b: List[SPState],
    spec: DeuteronSpec,
) -> Dict[Tuple[int, int, int, int], float]:
    """
    Compute two-body matrix elements between different species (no exchange).

    V_{ijkl} = <a_i b_j | V_NN | a_k b_l>  (direct only, no antisymmetrization)

    This is the same as compute_tbme from the deuteron builder.

    Parameters
    ----------
    states_a : list of SPState
        First species (proton).
    states_b : list of SPState
        Second species (neutron).
    spec : DeuteronSpec

    Returns
    -------
    dict (i, j, k, l) -> float (MeV)
    """
    return compute_tbme(states_a, states_b, spec)


# ---------------------------------------------------------------------------
# Slater determinant enumeration for multi-particle species
# ---------------------------------------------------------------------------

def _enumerate_slater_dets(n_orbitals: int, n_particles: int) -> List[Tuple[int, ...]]:
    """
    Enumerate all Slater determinants for n_particles in n_orbitals.

    Each SD is a sorted tuple of occupied orbital indices.

    Returns
    -------
    list of tuples
        Each tuple contains n_particles orbital indices in ascending order.
    """
    from itertools import combinations
    return [sd for sd in combinations(range(n_orbitals), n_particles)]


def _sd_phase(sd: Tuple[int, ...], p: int, q: int) -> Tuple[int, Tuple[int, ...]]:
    """
    Compute the phase and new SD from applying a^dag_p a_q to |sd>.

    Returns (phase, new_sd) where phase is +1 or -1, or (0, ()) if annihilation
    fails (q not in sd or p already in sd after removal).
    """
    if q not in sd:
        return 0, ()
    sd_list = list(sd)
    # Annihilate q: count swaps to bring q to position
    pos_q = sd_list.index(q)
    phase = (-1) ** pos_q
    sd_list.pop(pos_q)

    if p in sd_list:
        return 0, ()  # Pauli exclusion

    # Create p: find insertion position and count swaps
    sd_list.append(p)
    sd_list.sort()
    pos_p = sd_list.index(p)
    phase *= (-1) ** pos_p

    return phase, tuple(sd_list)


def _sd_phase_two_body(
    sd: Tuple[int, ...], p: int, q: int, r: int, s: int,
) -> Tuple[int, Tuple[int, ...]]:
    """
    Compute phase and new SD from a^dag_p a^dag_q a_s a_r |sd>.

    Note operator ordering: a_s a_r (destroy r first, then s).

    Returns (phase, new_sd) or (0, ()) if the operation fails.
    """
    # Annihilate r
    if r not in sd:
        return 0, ()
    sd_list = list(sd)
    pos_r = sd_list.index(r)
    phase = (-1) ** pos_r
    sd_list.pop(pos_r)

    # Annihilate s
    if s not in sd_list:
        return 0, ()
    pos_s = sd_list.index(s)
    phase *= (-1) ** pos_s
    sd_list.pop(pos_s)

    # Create q
    if q in sd_list:
        return 0, ()
    sd_list.append(q)
    sd_list.sort()
    pos_q = sd_list.index(q)
    phase *= (-1) ** pos_q

    # Create p
    if p in sd_list:
        return 0, ()
    sd_list.append(p)
    sd_list.sort()
    pos_p = sd_list.index(p)
    phase *= (-1) ** pos_p

    return phase, tuple(sd_list)


# ---------------------------------------------------------------------------
# He-4 FCI matrix builder
# ---------------------------------------------------------------------------

def _build_he4_fci_matrix(
    states_p: List[SPState],
    states_n: List[SPState],
    h1_p: np.ndarray,
    h1_n: np.ndarray,
    tbme_pp: Dict[Tuple[int, int, int, int], float],
    tbme_nn: Dict[Tuple[int, int, int, int], float],
    tbme_pn: Dict[Tuple[int, int, int, int], float],
) -> np.ndarray:
    """
    Build the FCI matrix for 2 protons + 2 neutrons.

    Basis: |SD_p> x |SD_n> where SD_p and SD_n are Slater determinants
    for 2 protons in n_p orbitals and 2 neutrons in n_n orbitals.
    """
    n_p = len(states_p)
    n_n = len(states_n)

    sds_p = _enumerate_slater_dets(n_p, 2)
    sds_n = _enumerate_slater_dets(n_n, 2)

    dim_p = len(sds_p)
    dim_n = len(sds_n)
    dim = dim_p * dim_n

    # Index lookup
    sd_p_idx = {sd: i for i, sd in enumerate(sds_p)}
    sd_n_idx = {sd: j for j, sd in enumerate(sds_n)}

    H = np.zeros((dim, dim))

    # One-body: protons
    for ip, sd_p in enumerate(sds_p):
        for a in range(n_p):
            for c in range(n_p):
                if abs(h1_p[a, c]) < 1e-15:
                    continue
                phase, new_sd = _sd_phase(sd_p, a, c)
                if phase == 0:
                    continue
                jp = sd_p_idx.get(new_sd)
                if jp is None:
                    continue
                for jn in range(dim_n):
                    bra = ip * dim_n + jn
                    ket = jp * dim_n + jn
                    H[bra, ket] += phase * h1_p[a, c]

    # One-body: neutrons
    for jn, sd_n in enumerate(sds_n):
        for b in range(n_n):
            for d in range(n_n):
                if abs(h1_n[b, d]) < 1e-15:
                    continue
                phase, new_sd = _sd_phase(sd_n, b, d)
                if phase == 0:
                    continue
                kn = sd_n_idx.get(new_sd)
                if kn is None:
                    continue
                for ip in range(dim_p):
                    bra = ip * dim_n + jn
                    ket = ip * dim_n + kn
                    H[bra, ket] += phase * h1_n[b, d]

    # Two-body: proton-proton (already antisymmetrized in tbme_pp)
    # H_pp = (1/4) sum_{ijkl} <ij||kl> a^dag_i a^dag_j a_l a_k
    # where <ij||kl> = <ij|kl> - <ij|lk> is the antisymmetrized ME.
    # The 1/4 factor accounts for the double-counting in the sum:
    # swapping (ij) or (kl) gives the same matrix element with a sign flip
    # that's absorbed by the antisymmetry of the operator product.
    for ip, sd_p in enumerate(sds_p):
        for (i, j, k, l), v in tbme_pp.items():
            if abs(v) < 1e-15:
                continue
            # a^dag_i a^dag_j a_l a_k |sd_p>
            phase, new_sd = _sd_phase_two_body(sd_p, i, j, k, l)
            if phase == 0:
                continue
            jp = sd_p_idx.get(new_sd)
            if jp is None:
                continue
            for jn in range(dim_n):
                bra = ip * dim_n + jn
                ket = jp * dim_n + jn
                H[bra, ket] += 0.25 * phase * v

    # Two-body: neutron-neutron (already antisymmetrized in tbme_nn)
    for jn, sd_n in enumerate(sds_n):
        for (i, j, k, l), v in tbme_nn.items():
            if abs(v) < 1e-15:
                continue
            phase, new_sd = _sd_phase_two_body(sd_n, i, j, k, l)
            if phase == 0:
                continue
            kn = sd_n_idx.get(new_sd)
            if kn is None:
                continue
            for ip in range(dim_p):
                bra = ip * dim_n + jn
                ket = ip * dim_n + kn
                H[bra, ket] += 0.25 * phase * v

    # Two-body: proton-neutron (no exchange between species)
    # V_{ijkl} (a^dag_{p_i} a_{p_k}) (a^dag_{n_j} a_{n_l})
    for ip, sd_p in enumerate(sds_p):
        for jn, sd_n in enumerate(sds_n):
            for (i, j, k, l), v in tbme_pn.items():
                if abs(v) < 1e-15:
                    continue
                # Proton: a^dag_i a_k
                phase_p, new_sd_p = _sd_phase(sd_p, i, k)
                if phase_p == 0:
                    continue
                kp = sd_p_idx.get(new_sd_p)
                if kp is None:
                    continue
                # Neutron: a^dag_j a_l
                phase_n, new_sd_n = _sd_phase(sd_n, j, l)
                if phase_n == 0:
                    continue
                kn = sd_n_idx.get(new_sd_n)
                if kn is None:
                    continue

                bra = ip * dim_n + jn
                ket = kp * dim_n + kn
                H[bra, ket] += phase_p * phase_n * v

    return H


# ---------------------------------------------------------------------------
# He-4 JW Pauli Hamiltonian builder
# ---------------------------------------------------------------------------

def _two_body_pauli_same_species(
    i: int, j: int, k: int, l: int,
    n_qubits: int,
) -> Dict[str, float]:
    """
    JW encoding of (1/2) * V_{ijkl}_AS * a^dag_i a^dag_j a_l a_k
    on a single-species register.

    Uses the identity:
    a^dag_i a^dag_j a_l a_k = (a^dag_i a_k)(a^dag_j a_l) - delta_{jk}(a^dag_i a_l)

    But it is simpler to encode directly using the Pauli decomposition of
    the two-body operator.

    For the antisymmetrized matrix element, we encode:
    (1/2) V_{ijkl}_AS * [n_i n_j delta_{ik} delta_{jl}
                          - various one-body and two-body Pauli terms]

    Actually, the standard approach is:
    a^dag_p a^dag_q a_s a_r = (a^dag_p a_r)(a^dag_q a_s) - delta_{qr}(a^dag_p a_s)
    """
    result: Dict[str, float] = {}

    # a^dag_i a^dag_j a_l a_k = (a^dag_i a_k)(a^dag_j a_l) - delta_{jk}(a^dag_i a_l)
    # First term: product of one-body operators
    pauli_ik = _one_body_pauli(i, k, n_qubits)
    pauli_jl = _one_body_pauli(j, l, n_qubits)

    # Multiply Pauli dictionaries (within same register)
    for s1, c1 in pauli_ik.items():
        for s2, c2 in pauli_jl.items():
            # Multiply Pauli strings on same register
            coeff = c1 * c2
            product_str, phase = _multiply_pauli_strings(s1, s2)
            coeff *= phase
            if abs(coeff) > 1e-15:
                result[product_str] = result.get(product_str, 0.0) + coeff

    # Second term: -delta_{jk} * a^dag_i a_l
    if j == k:
        pauli_il = _one_body_pauli(i, l, n_qubits)
        for s, c in pauli_il.items():
            result[s] = result.get(s, 0.0) - c

    return result


def _multiply_pauli_strings(s1: str, s2: str) -> Tuple[str, complex]:
    """
    Multiply two Pauli strings on the same register.

    Returns (result_string, phase) where phase includes factors of i from
    Pauli multiplication rules.
    """
    result = []
    phase = 1.0 + 0.0j
    for a, b in zip(s1, s2):
        p, ph = _multiply_single_pauli(a, b)
        result.append(p)
        phase *= ph
    return ''.join(result), phase


def _multiply_single_pauli(a: str, b: str) -> Tuple[str, complex]:
    """Multiply two single-qubit Paulis. Returns (result, phase)."""
    if a == 'I':
        return b, 1.0
    if b == 'I':
        return a, 1.0
    if a == b:
        return 'I', 1.0

    # Non-trivial products with phases
    table = {
        ('X', 'Y'): ('Z', 1j),
        ('Y', 'X'): ('Z', -1j),
        ('Y', 'Z'): ('X', 1j),
        ('Z', 'Y'): ('X', -1j),
        ('Z', 'X'): ('Y', 1j),
        ('X', 'Z'): ('Y', -1j),
    }
    return table[(a, b)]


def build_he4_hamiltonian(
    N_shells: int = 2,
    hw: float = 15.0,
    include_coulomb: bool = True,
) -> Dict[str, Any]:
    """
    Build the He-4 qubit Hamiltonian for 2 protons + 2 neutrons.

    Parameters
    ----------
    N_shells : int
        Number of HO shells.
    hw : float
        HO frequency in MeV.
    include_coulomb : bool
        Whether to include pp Coulomb repulsion.

    Returns
    -------
    dict with keys: H_pauli, H_matrix, states_p, states_n,
        sds_p, sds_n, Q, Q_p, Q_n, metadata
    """
    spec = DeuteronSpec(N_shells=N_shells, hw=hw)
    states_p, states_n = enumerate_sp_states(spec)
    n_p = len(states_p)
    n_n = len(states_n)
    Q = n_p + n_n

    h1_p = build_one_body(states_p, spec)
    h1_n = build_one_body(states_n, spec)

    # Two-body matrix elements
    # pp: antisymmetrized, optionally with Coulomb
    tbme_pp = compute_same_species_tbme(states_p, spec,
                                         include_coulomb=include_coulomb)
    # nn: antisymmetrized, no Coulomb
    tbme_nn = compute_same_species_tbme(states_n, spec,
                                         include_coulomb=False)
    # pn: direct only (no antisymmetrization between species)
    tbme_pn = compute_cross_species_tbme(states_p, states_n, spec)

    # Build FCI matrix
    H_matrix = _build_he4_fci_matrix(
        states_p, states_n, h1_p, h1_n, tbme_pp, tbme_nn, tbme_pn)

    # Build Pauli Hamiltonian
    H_pauli: Dict[str, float] = {}

    # One-body: proton
    for i in range(n_p):
        if abs(h1_p[i, i]) > 1e-15:
            pauli_p = _number_op_pauli(i, n_p)
            embedded = _embed_single_species(pauli_p, 'proton', n_p, n_n)
            _add_pauli(H_pauli, embedded, h1_p[i, i])

    # One-body: neutron
    for j in range(n_n):
        if abs(h1_n[j, j]) > 1e-15:
            pauli_n = _number_op_pauli(j, n_n)
            embedded = _embed_single_species(pauli_n, 'neutron', n_p, n_n)
            _add_pauli(H_pauli, embedded, h1_n[j, j])

    # Two-body: pp (antisymmetrized), factor of 1/4 for double counting
    for (i, j, k, l), v in tbme_pp.items():
        if abs(v) < 1e-15:
            continue
        pauli_2b = _two_body_pauli_same_species(i, j, k, l, n_p)
        # Embed in full register (proton part only)
        for pstr, coeff in pauli_2b.items():
            full_coeff = 0.25 * v * coeff
            if abs(full_coeff) < 1e-15:
                continue
            # Only keep real part (imaginary parts cancel for Hermitian H)
            if abs(full_coeff.imag if isinstance(full_coeff, complex) else 0) > abs(full_coeff.real if isinstance(full_coeff, complex) else full_coeff) * 0.01:
                continue
            real_coeff = full_coeff.real if isinstance(full_coeff, complex) else full_coeff
            full_str = pstr + 'I' * n_n
            H_pauli[full_str] = H_pauli.get(full_str, 0.0) + real_coeff

    # Two-body: nn (antisymmetrized), factor of 1/4
    for (i, j, k, l), v in tbme_nn.items():
        if abs(v) < 1e-15:
            continue
        pauli_2b = _two_body_pauli_same_species(i, j, k, l, n_n)
        for pstr, coeff in pauli_2b.items():
            full_coeff = 0.25 * v * coeff
            if abs(full_coeff) < 1e-15:
                continue
            if abs(full_coeff.imag if isinstance(full_coeff, complex) else 0) > abs(full_coeff.real if isinstance(full_coeff, complex) else full_coeff) * 0.01:
                continue
            real_coeff = full_coeff.real if isinstance(full_coeff, complex) else full_coeff
            full_str = 'I' * n_p + pstr
            H_pauli[full_str] = H_pauli.get(full_str, 0.0) + real_coeff

    # Two-body: pn (direct, cross-species tensor product)
    for (i, j, k, l), v in tbme_pn.items():
        if abs(v) < 1e-15:
            continue
        pauli_p = _one_body_pauli(i, k, n_p)
        pauli_n = _one_body_pauli(j, l, n_n)
        pn_pauli = _tensor_product_pauli(pauli_p, pauli_n)
        _add_pauli(H_pauli, pn_pauli, v)

    # Clean up small terms and imaginary residuals
    H_pauli = {k: v for k, v in H_pauli.items() if abs(v) > 1e-12}

    sds_p = _enumerate_slater_dets(n_p, 2)
    sds_n = _enumerate_slater_dets(n_n, 2)

    metadata = {
        'system': 'he4',
        'N_shells': N_shells,
        'hw_MeV': hw,
        'b_fm': spec.b,
        'energy_unit': 'MeV',
        'Q': Q, 'Q_p': n_p, 'Q_n': n_n,
        'include_coulomb': include_coulomb,
        'n_protons': 2, 'n_neutrons': 2,
        'dim_p': len(sds_p), 'dim_n': len(sds_n),
        'dim_total': len(sds_p) * len(sds_n),
    }

    return {
        'H_pauli': H_pauli,
        'H_matrix': H_matrix,
        'states_p': states_p,
        'states_n': states_n,
        'sds_p': sds_p,
        'sds_n': sds_n,
        'Q': Q, 'Q_p': n_p, 'Q_n': n_n,
        'metadata': metadata,
    }


def diagonalize_he4(
    N_shells: int = 2,
    hw: float = 15.0,
    include_coulomb: bool = True,
) -> Dict[str, Any]:
    """
    Build and diagonalize the He-4 Hamiltonian.

    Parameters
    ----------
    N_shells : int
        Number of HO shells.
    hw : float
        HO frequency in MeV.
    include_coulomb : bool
        Whether to include pp Coulomb repulsion.

    Returns
    -------
    dict with E_gs, eigenvalues, ground_state, H_data, metadata
    """
    data = build_he4_hamiltonian(N_shells=N_shells, hw=hw,
                                  include_coulomb=include_coulomb)
    H = data['H_matrix']
    evals, evecs = np.linalg.eigh(H)
    E_gs = evals[0]

    return {
        'E_gs': E_gs,
        'eigenvalues': evals,
        'ground_state': evecs[:, 0],
        'H_data': data,
        'metadata': data['metadata'],
    }


def analyze_he4_hamiltonian(
    H_pauli: Dict[str, float],
    Q: int,
) -> Dict[str, Any]:
    """Analyze the He-4 Pauli Hamiltonian (same interface as deuteron)."""
    return analyze_deuteron_hamiltonian(H_pauli, Q)


def he4_hw_scan(
    N_shells: int = 2,
    hw_values: Optional[List[float]] = None,
    include_coulomb: bool = True,
) -> Dict[str, Any]:
    """Scan hw to find optimal oscillator frequency for He-4."""
    if hw_values is None:
        hw_values = list(np.arange(10.0, 35.0, 2.5))

    energies = []
    for hw in hw_values:
        try:
            result = diagonalize_he4(N_shells=N_shells, hw=hw,
                                      include_coulomb=include_coulomb)
            energies.append(result['E_gs'])
        except Exception:
            energies.append(float('nan'))

    energies = np.array(energies)
    valid = ~np.isnan(energies)
    if np.any(valid):
        idx_opt = int(np.argmin(energies[valid]))
        valid_hw = np.array(hw_values)[valid]
        return {
            'hw_values': hw_values,
            'energies': energies.tolist(),
            'optimal_hw': float(valid_hw[idx_opt]),
            'optimal_E': float(energies[valid][idx_opt]),
        }
    return {'hw_values': hw_values, 'energies': energies.tolist(),
            'optimal_hw': None, 'optimal_E': None}
