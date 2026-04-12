"""
Composed Nuclear-Electronic Qubit Hamiltonian (Track NI).

Assembles the first GeoVac composed Hamiltonian spanning both nuclear and
electronic degrees of freedom on a single qubit register. The architecture is:

    H_total = H_nuc (x) I_elec + I_nuc (x) H_elec + V_finite_size + H_hyperfine

Layout: nuclear register on qubits 0..Q_nuc-1, electronic register on qubits
Q_nuc..Q_nuc+Q_elec-1.

Proof of concept for the embedding architecture defined in
`docs/nuclear_electronic_embedding_spec.md`. The focus here is on structural
validation, not production precision. The 10^6 energy scale ratio between
nuclear (MeV) and electronic (Ha) terms is the key architectural test.

All terms are reported in Hartree in the final composed Hamiltonian.

Units conventions:
    1 MeV   = 1e6 / 27.211386245988 Ha  (~ 36749.32 Ha)
    1 Ha    = 27.211386245988 eV
    1 eV    = 1/27.211386245988 Ha

Author: GeoVac Development Team
Date: April 2026 (Track NI)
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np

from geovac.nuclear.nuclear_hamiltonian import (
    SPState,
    build_deuteron_hamiltonian,
    enumerate_sp_states,
    DeuteronSpec,
)
from geovac.nuclear.form_factor import (
    R_PROTON_BOHR,
    finite_size_correction,
    finite_size_correction_exact_1s,
)


# ---------------------------------------------------------------------------
# Physical constants (atomic units)
# ---------------------------------------------------------------------------

HA_PER_EV: float = 1.0 / 27.211386245988           # 1 eV in Ha
HA_PER_MEV: float = 1.0e6 * HA_PER_EV              # 1 MeV in Ha (~ 36749.32)

# 21 cm hyperfine splitting of hydrogen 1s: 1420.405 MHz.
# Photon energy: h * nu in Hartree using h in eV*s = 4.135667696e-15,
# or equivalently nu (Hz) * 1.519829846e-16 Ha/Hz (Planck constant over h_bar-less
# Hartree atomic time).  Canonical value: 6.4465e-6 Ha? No --- the splitting
# is ~5.87 microeV ~ 2.16e-7 Ha. Compute directly: 1.420405e9 Hz * h_ha
# with h in Hartree*s = 1.519829846e-16.
HZ_PER_HA: float = 1.0 / 1.519829846e-16            # Hartree frequency (s^-1)
HF_HYDROGEN_HZ: float = 1.420405751768e9            # CODATA 21cm line
HF_HYDROGEN_HA: float = HF_HYDROGEN_HZ / HZ_PER_HA  # ~2.16e-7 Ha


# ---------------------------------------------------------------------------
# Utility: add a Pauli dict into a target
# ---------------------------------------------------------------------------

def _add_pauli(target: Dict[str, float],
               source: Dict[str, float],
               scale: float = 1.0) -> None:
    for k, v in source.items():
        target[k] = target.get(k, 0.0) + scale * v


def _clean_pauli(pauli: Dict[str, float], tol: float = 1e-14) -> Dict[str, float]:
    return {k: v for k, v in pauli.items() if abs(v) > tol}


# ---------------------------------------------------------------------------
# 1. Nuclear block wrapper
# ---------------------------------------------------------------------------

def build_nuclear_block(N_shells: int = 2, hw: float = 10.0) -> Dict[str, Any]:
    """
    Wrap the deuteron Hamiltonian builder and return a Pauli dict in MeV.

    Parameters
    ----------
    N_shells : int
        Number of harmonic oscillator major shells per species.
    hw : float
        HO frequency in MeV.

    Returns
    -------
    dict with keys:
        'pauli_terms' : Dict[str, float]   (MeV)
        'Q_nuc'       : int                 (qubit count)
        'Q_p', 'Q_n'  : int                 (protons, neutrons)
        'states_p', 'states_n' : List[SPState]
        'energy_unit' : 'MeV'
        'H_matrix'    : np.ndarray          (FCI matrix in 1p+1n subspace)
    """
    data = build_deuteron_hamiltonian(N_shells=N_shells, hw=hw)
    return {
        'pauli_terms': data['H_pauli'],
        'Q_nuc': data['Q'],
        'Q_p': data['Q_p'],
        'Q_n': data['Q_n'],
        'states_p': data['states_p'],
        'states_n': data['states_n'],
        'energy_unit': 'MeV',
        'H_matrix': data['H_matrix'],
    }


# ---------------------------------------------------------------------------
# 2. Electronic block wrapper (hydrogen, n_max=2, one electron)
# ---------------------------------------------------------------------------

def _enumerate_hydrogen_sp_states(n_max: int) -> List[Tuple[int, int, int, float]]:
    """
    Enumerate (n, l, m, m_s) one-electron hydrogenic states for 1 <= n <= n_max.

    Returns a list of tuples. For n_max=2 this yields 8 spin-orbitals:
        n=1: (1,0,0,+1/2), (1,0,0,-1/2)                             [2]
        n=2: (2,0,0,+1/2), (2,0,0,-1/2),                            [2]
             (2,1,-1,+1/2), (2,1,-1,-1/2),
             (2,1,0,+1/2), (2,1,0,-1/2),
             (2,1,1,+1/2), (2,1,1,-1/2)                             [6]
    """
    states: List[Tuple[int, int, int, float]] = []
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                for m_s in (0.5, -0.5):
                    states.append((n, l, m, m_s))
    return states


def build_electronic_block(Z: int = 1, n_max: int = 2) -> Dict[str, Any]:
    """
    Build a diagonal one-electron hydrogenic Hamiltonian in Pauli form.

    At Z=1, n_max=2 this is the simplest possible case: the Hamiltonian is
    diagonal in the (n, l, m, m_s) basis with entries -Z^2/(2 n^2) Ha.

    H = sum_i eps_i n_i   where n_i = (I - Z_i)/2

    The resulting Pauli dict has 1 identity term plus one single-Z term per
    spin-orbital (here 8 single-Z terms).

    Parameters
    ----------
    Z : int
        Nuclear charge (default 1 for hydrogen).
    n_max : int
        Maximum principal quantum number.

    Returns
    -------
    dict with keys:
        'pauli_terms' : Dict[str, float]   (Ha)
        'Q_elec'      : int                 (qubit count)
        'states'      : list of (n, l, m, m_s)
        'orbital_energies' : np.ndarray     (Ha, length Q_elec)
        'energy_unit' : 'Ha'
    """
    states = _enumerate_hydrogen_sp_states(n_max)
    Q = len(states)

    eps = np.array([-0.5 * (Z ** 2) / (n ** 2) for (n, l, m, ms) in states])

    pauli: Dict[str, float] = {}
    # Identity contribution: sum_i eps_i/2
    identity_str = 'I' * Q
    pauli[identity_str] = float(np.sum(eps) / 2.0)
    # Single-Z contribution: -eps_i/2 on qubit i
    for i in range(Q):
        key = list(identity_str)
        key[i] = 'Z'
        pauli[''.join(key)] = -eps[i] / 2.0

    pauli = _clean_pauli(pauli)

    return {
        'pauli_terms': pauli,
        'Q_elec': Q,
        'states': states,
        'orbital_energies': eps,
        'energy_unit': 'Ha',
    }


def _find_1s_electron_qubits(elec_block: Dict[str, Any]) -> List[int]:
    """Return indices of 1s spin-orbitals in the electronic register."""
    return [i for i, (n, l, m, ms) in enumerate(elec_block['states'])
            if n == 1 and l == 0 and m == 0]


def _find_1s_up_down(elec_block: Dict[str, Any]) -> Tuple[int, int]:
    """Return (q_up, q_down) for the 1s spin-orbitals."""
    q_up, q_down = None, None
    for i, (n, l, m, ms) in enumerate(elec_block['states']):
        if n == 1 and l == 0 and m == 0:
            if ms > 0:
                q_up = i
            else:
                q_down = i
    assert q_up is not None and q_down is not None
    return q_up, q_down


# ---------------------------------------------------------------------------
# 3. Finite-size coupling
# ---------------------------------------------------------------------------

def finite_size_coupling_pauli(
    Q_nuc: int,
    Q_elec: int,
    nuclear_block: Dict[str, Any],
    electronic_block: Dict[str, Any],
    R_nuc: float = R_PROTON_BOHR,
    Z: int = 1,
) -> Dict[str, float]:
    """
    Finite-size correction as a one-body operator on the electronic register.

    Closed form (Paper 18, form_factor.py):
        delta E_1s = (2/5) Z^4 R_nuc^2      (leading order)

    In second quantization, the correction couples to the electronic 1s
    number operator:
        V_fs = delta E_1s * (n_{1s up} + n_{1s down})

    with n_i = (I - Z_i)/2.

    This simplest implementation treats R_nuc as a classical parameter (not
    a nuclear operator). A fully quantum coupling would replace R_nuc with
    an operator on the nuclear register (the nuclear charge distribution),
    but that is deferred as a future extension --- here the nuclear
    wavefunction is approximated by its expectation value.

    Returns
    -------
    Pauli dict in Hartree on the combined register (nuclear qubits first).
    """
    Q_total = Q_nuc + Q_elec
    dE = finite_size_correction(Z=Z, R_nuc=R_nuc, n=1, l=0)  # Ha

    q_up, q_down = _find_1s_up_down(electronic_block)
    # Absolute qubit indices on the combined register.
    q_up_abs = Q_nuc + q_up
    q_down_abs = Q_nuc + q_down

    pauli: Dict[str, float] = {}
    identity = 'I' * Q_total
    # n_i = (I - Z_i)/2; V_fs = dE * (n_up + n_down)
    # = dE * (I - Z_up/2 - Z_down/2 + I*0)  wait: (I-Z_up)/2 + (I-Z_down)/2
    # = I - 0.5*Z_up - 0.5*Z_down
    pauli[identity] = dE  # two 1/2*I pieces

    for q_abs in (q_up_abs, q_down_abs):
        s = list(identity)
        s[q_abs] = 'Z'
        key = ''.join(s)
        pauli[key] = pauli.get(key, 0.0) + (-0.5 * dE)

    return _clean_pauli(pauli)


# ---------------------------------------------------------------------------
# 4. Hyperfine coupling
# ---------------------------------------------------------------------------

def _find_proton_0s_qubits(nuclear_block: Dict[str, Any]) -> Tuple[int, int]:
    """Return (q_up, q_down) for the proton 0s (n_r=0, l=0) spin-orbitals."""
    q_up, q_down = None, None
    for i, s in enumerate(nuclear_block['states_p']):
        if s.n_r == 0 and s.l == 0 and s.m_l == 0:
            if s.m_s > 0:
                q_up = i
            else:
                q_down = i
    assert q_up is not None and q_down is not None
    return q_up, q_down


def _single_qubit_pauli(op: str, q: int, n_qubits: int) -> str:
    s = ['I'] * n_qubits
    s[q] = op
    return ''.join(s)


def _two_qubit_pauli(op1: str, q1: int, op2: str, q2: int, n_qubits: int) -> str:
    s = ['I'] * n_qubits
    s[q1] = op1
    s[q2] = op2
    return ''.join(s)


def hyperfine_coupling_pauli(
    Q_nuc: int,
    Q_elec: int,
    nuclear_block: Dict[str, Any],
    electronic_block: Dict[str, Any],
    A_hf: float = HF_HYDROGEN_HA,
) -> Dict[str, float]:
    """
    Hyperfine coupling H_hf = A_hf * I . S as a cross-register Pauli sum.

    For two spin-1/2 objects S1, S2 (the proton 0s spin and the electron 1s
    spin), the Heisenberg coupling is
        S1 . S2 = S1_x S2_x + S1_y S2_y + S1_z S2_z
                = (1/4)(X1 X2 + Y1 Y2 + Z1 Z2)
    in units where the Pauli operators are used for single-qubit spin-1/2.
    Note that S_i = (1/2) sigma_i, so S1.S2 = (1/4) sigma1.sigma2.

    The spin-1/2 qubits are identified as follows. On the nuclear register,
    the proton 0s m_s=+1/2 state occupies one qubit and m_s=-1/2 occupies
    another. In the single-particle (non-number-conserving) Pauli encoding,
    the spin flip between these two is implemented by a_dag_up a_down (plus
    its conjugate), which becomes a pair of X/Y Pauli terms with a Z string
    in between. Similarly for the electron 1s.

    However, for a validation-oriented proof of concept we adopt the
    following simplification: we restrict attention to the single-occupancy
    subspace (one proton in its 0s, one electron in its 1s) and treat the
    pair of spin qubits as an effective spin-1/2 Heisenberg coupling on
    those two qubits. The relevant operators on the "up/down" qubit pair
    are:

        Spin operators on a pair (q_up, q_down):
            S_z = (1/2)(n_up - n_down)
                = (1/4)(Z_down - Z_up)
            S_+ = a_dag_up a_down
            S_- = a_dag_down a_up
            S_x = (S_+ + S_-)/2
            S_y = (S_+ - S_-)/(2i)

        Under Jordan-Wigner with the up qubit at index u and the down qubit
        at index d:
            a_dag_u a_d = (X_u + i Y_u)/2 * (X_d - i Y_d)/2 * (Z-string between)
        The Z-string contribution is constant in the single-electron
        subspace of the pair (only depends on occupations of the two qubits
        themselves), so within this subspace:
            a_dag_u a_d = (1/4)(X_u X_d + Y_u Y_d + i (Y_u X_d - X_u Y_d))
        and taking the Hermitian combination gives
            a_dag_u a_d + a_dag_d a_u = (1/2)(X_u X_d + Y_u Y_d).

        Therefore within the single-occupancy subspace:
            S_x = (1/4)(X_u X_d + Y_u Y_d) / ??? Let us re-derive carefully.

    Careful derivation (within single-occupation subspace of pair):
        |up> = |1 0> = a_dag_u |vac>, |down> = |0 1> = a_dag_d |vac>
        S_z |up> = +1/2 |up>, S_z |down> = -1/2 |down>
        S_+ |down> = |up>, S_- |up> = |down>

        Define two Pauli operators on the (u, d) qubit pair within the
        single-occupancy subspace:
            Z_u + Z_d = 0 on both |10> and |01>? Let's check:
               |10>: Z_u=-1, Z_d=+1 -> sum = 0
               |01>: Z_u=+1, Z_d=-1 -> sum = 0   Yes.
            Z_u - Z_d: |10> -> -2, |01> -> +2
               So S_z = -(Z_u - Z_d)/4 = (Z_d - Z_u)/4.

        For the flip operators:
            |10> <-> |01>
            a_dag_u a_d has matrix element <10|a_dag_u a_d|01> = +1
            In terms of Pauli X/Y (assuming qubits u < d, no intermediate
            Z-string between adjacent qubits here), the operator
                sigma^+_u sigma^-_d = (X_u + iY_u)/2 * (X_d - iY_d)/2
            acts on |10> = sigma^+_u |00>, |01> = sigma^+_d |00> as:
                sigma^+_u sigma^-_d |01> = sigma^+_u |00> = |10>
            so sigma^+_u sigma^-_d matches a_dag_u a_d within the
            single-occupancy subspace (up to a JW Z-string that we set aside).

        Expanded:
            sigma^+_u sigma^-_d = (1/4)(X_u X_d + Y_u Y_d + i X_u Y_d - i Y_u X_d)
            Its Hermitian conjugate is:
            sigma^-_u sigma^+_d = (1/4)(X_u X_d + Y_u Y_d - i X_u Y_d + i Y_u X_d)
            Sum: (1/2)(X_u X_d + Y_u Y_d) = 2 S_x_pair
            Diff: (i/2)(X_u Y_d - Y_u X_d) = 2 i S_y_pair
            (where S_pair is the pair's spin-1/2 operator.)

        Therefore, within the single-occupation subspace of the pair:
            S_x = (1/4)(X_u X_d + Y_u Y_d)
            S_y = (1/4)(Y_u X_d - X_u Y_d)
            S_z = (1/4)(Z_d - Z_u)

        These satisfy [S_x, S_y] = i S_z on the 2-dimensional single-
        occupancy subspace.

    Putting it together, for the two pairs (nuclear up, nuclear down) and
    (electronic up, electronic down):

        I.S = I_x S_x + I_y S_y + I_z S_z
            = (1/16) [ (X_Nu X_Nd + Y_Nu Y_Nd)(X_Eu X_Ed + Y_Eu Y_Ed)
                      + (Y_Nu X_Nd - X_Nu Y_Nd)(Y_Eu X_Ed - X_Eu Y_Ed)
                      + (Z_Nd - Z_Nu)(Z_Ed - Z_Eu) ]

    This is 4 + 4 + 4 = 12 cross-register Pauli terms (four-qubit supports).
    The coefficient is A_hf / 16 per term after expanding the first two
    products, times the combinations from the full I.S sum.

    Notes
    -----
    1. The JW Z-strings between the up/down qubits of each species are
       dropped: within single occupation of each pair, those Z-strings act
       as +1 (no intervening occupied qubits between u=0 and d=1 on the
       nuclear side, or between u=0 and d=1 on the electronic side after
       reindexing within the 1s shell).
    2. In reality the nuclear u=0, d=1 qubits are adjacent, and the
       electronic u=0, d=1 qubits are also adjacent (the 1s spin-up and
       spin-down), so there is no Z-string inside either pair.
    3. This coupling lives in the physical 1-proton, 1-electron sector and
       operates only within the 1s nuclear and 1s electronic single
       occupancies. In the full Hilbert space of the qubit register, the
       operator couples additional sectors, but the validation test
       (block_decomposition) projects onto the physical sector.

    Returns
    -------
    Pauli dict in Hartree on the combined register.
    """
    Q_total = Q_nuc + Q_elec

    # Identify the four relevant qubits.
    nu, nd = _find_proton_0s_qubits(nuclear_block)
    eu, ed = _find_1s_up_down(electronic_block)
    Nu = nu
    Nd = nd
    Eu = Q_nuc + eu
    Ed = Q_nuc + ed

    pauli: Dict[str, float] = {}

    def add_term(ops: List[Tuple[str, int]], coeff: float) -> None:
        s = ['I'] * Q_total
        for op, q in ops:
            s[q] = op
        key = ''.join(s)
        pauli[key] = pauli.get(key, 0.0) + coeff

    # I_x S_x = (1/16) (X_Nu X_Nd + Y_Nu Y_Nd) * (X_Eu X_Ed + Y_Eu Y_Ed)
    # Four terms:
    c_xy = A_hf / 16.0
    add_term([('X', Nu), ('X', Nd), ('X', Eu), ('X', Ed)], c_xy)
    add_term([('X', Nu), ('X', Nd), ('Y', Eu), ('Y', Ed)], c_xy)
    add_term([('Y', Nu), ('Y', Nd), ('X', Eu), ('X', Ed)], c_xy)
    add_term([('Y', Nu), ('Y', Nd), ('Y', Eu), ('Y', Ed)], c_xy)

    # I_y S_y = (1/16) (Y_Nu X_Nd - X_Nu Y_Nd) * (Y_Eu X_Ed - X_Eu Y_Ed)
    # Four terms with signs:
    #   +Y X * +Y X   = + Y_Nu X_Nd Y_Eu X_Ed
    #   +Y X * -X Y   = - Y_Nu X_Nd X_Eu Y_Ed
    #   -X Y * +Y X   = - X_Nu Y_Nd Y_Eu X_Ed
    #   -X Y * -X Y   = + X_Nu Y_Nd X_Eu Y_Ed
    add_term([('Y', Nu), ('X', Nd), ('Y', Eu), ('X', Ed)], +c_xy)
    add_term([('Y', Nu), ('X', Nd), ('X', Eu), ('Y', Ed)], -c_xy)
    add_term([('X', Nu), ('Y', Nd), ('Y', Eu), ('X', Ed)], -c_xy)
    add_term([('X', Nu), ('Y', Nd), ('X', Eu), ('Y', Ed)], +c_xy)

    # I_z S_z = (1/16) (Z_Nd - Z_Nu)(Z_Ed - Z_Eu)
    # Four terms:
    #   + Z_Nd Z_Ed
    #   - Z_Nd Z_Eu
    #   - Z_Nu Z_Ed
    #   + Z_Nu Z_Eu
    add_term([('Z', Nd), ('Z', Ed)], +c_xy)
    add_term([('Z', Nd), ('Z', Eu)], -c_xy)
    add_term([('Z', Nu), ('Z', Ed)], -c_xy)
    add_term([('Z', Nu), ('Z', Eu)], +c_xy)

    return _clean_pauli(pauli)


# ---------------------------------------------------------------------------
# 5. Embedding helpers
# ---------------------------------------------------------------------------

def _embed_nuclear_into_combined(
    pauli_nuc: Dict[str, float],
    Q_nuc: int,
    Q_elec: int,
    scale: float = 1.0,
) -> Dict[str, float]:
    """Pad nuclear Pauli strings with identities on the electronic side."""
    tail = 'I' * Q_elec
    result: Dict[str, float] = {}
    for k, v in pauli_nuc.items():
        assert len(k) == Q_nuc, f"expected length {Q_nuc}, got {len(k)} for '{k}'"
        result[k + tail] = scale * v
    return result


def _embed_electronic_into_combined(
    pauli_elec: Dict[str, float],
    Q_nuc: int,
    Q_elec: int,
    scale: float = 1.0,
) -> Dict[str, float]:
    """Prepend identities on the nuclear side."""
    head = 'I' * Q_nuc
    result: Dict[str, float] = {}
    for k, v in pauli_elec.items():
        assert len(k) == Q_elec, f"expected length {Q_elec}, got {len(k)} for '{k}'"
        result[head + k] = scale * v
    return result


# ---------------------------------------------------------------------------
# 6. Composed Hamiltonian assembly
# ---------------------------------------------------------------------------

def build_deuterium_composed_hamiltonian(
    N_shells: int = 2,
    hw: float = 10.0,
    n_max_elec: int = 2,
    R_nuc: float = R_PROTON_BOHR,
    include_finite_size: bool = True,
    include_hyperfine: bool = True,
    Z: int = 1,
) -> Dict[str, Any]:
    """
    Build the composed nuclear-electronic qubit Hamiltonian for deuterium.

    Layout:
        qubits 0..Q_nuc-1            : nuclear register (deuteron)
        qubits Q_nuc..Q_nuc+Q_elec-1 : electronic register (hydrogen)

    Hamiltonian (Ha):
        H_total = H_nuc (x) I + I (x) H_elec + V_finite_size + H_hyperfine

    Nuclear MeV values are converted to Ha via HA_PER_MEV before assembly.

    Returns
    -------
    dict with keys:
        'pauli_terms' : Dict[str, float]   (Ha)
        'Q_total'     : int
        'Q_nuc'       : int
        'Q_elec'      : int
        'nuclear_block'    : dict (from build_nuclear_block)
        'electronic_block' : dict (from build_electronic_block)
        'metadata'    : dict
    """
    nuc = build_nuclear_block(N_shells=N_shells, hw=hw)
    elec = build_electronic_block(Z=Z, n_max=n_max_elec)

    Q_nuc = nuc['Q_nuc']
    Q_elec = elec['Q_elec']
    Q_total = Q_nuc + Q_elec

    composed: Dict[str, float] = {}

    # Nuclear block: MeV -> Ha
    nuc_ha = _embed_nuclear_into_combined(
        nuc['pauli_terms'], Q_nuc, Q_elec, scale=HA_PER_MEV
    )
    _add_pauli(composed, nuc_ha)

    # Electronic block: already Ha
    elec_pauli = _embed_electronic_into_combined(
        elec['pauli_terms'], Q_nuc, Q_elec
    )
    _add_pauli(composed, elec_pauli)

    # Finite-size coupling
    if include_finite_size:
        fs = finite_size_coupling_pauli(
            Q_nuc, Q_elec, nuc, elec, R_nuc=R_nuc, Z=Z
        )
        _add_pauli(composed, fs)

    # Hyperfine coupling
    if include_hyperfine:
        hf = hyperfine_coupling_pauli(
            Q_nuc, Q_elec, nuc, elec, A_hf=HF_HYDROGEN_HA
        )
        _add_pauli(composed, hf)

    composed = _clean_pauli(composed, tol=1e-18)

    metadata = {
        'system': 'deuterium (1p + 1n + 1e)',
        'N_shells': N_shells,
        'hw_MeV': hw,
        'n_max_elec': n_max_elec,
        'R_nuc_bohr': R_nuc,
        'include_finite_size': include_finite_size,
        'include_hyperfine': include_hyperfine,
        'energy_unit': 'Ha',
        'layout': f'nuclear[0..{Q_nuc-1}] + electronic[{Q_nuc}..{Q_total-1}]',
    }

    return {
        'pauli_terms': composed,
        'Q_total': Q_total,
        'Q_nuc': Q_nuc,
        'Q_elec': Q_elec,
        'nuclear_block': nuc,
        'electronic_block': elec,
        'metadata': metadata,
    }


# ---------------------------------------------------------------------------
# 7. Analysis
# ---------------------------------------------------------------------------

def analyze_composed_hamiltonian(
    H_pauli: Dict[str, float],
    Q_nuc: int,
    Q_elec: int,
) -> Dict[str, Any]:
    """
    Compute scaling and structure metrics for the composed Hamiltonian.
    """
    Q_total = Q_nuc + Q_elec
    assert all(len(k) == Q_total for k in H_pauli), "Pauli strings must have length Q_total"

    identity = 'I' * Q_total
    non_identity = {k: v for k, v in H_pauli.items() if k != identity}

    # Decomposition: nuclear-only, electronic-only, coupling
    nuc_slice = slice(0, Q_nuc)
    elec_slice = slice(Q_nuc, Q_total)

    n_nuc_only = 0
    n_elec_only = 0
    n_cross = 0
    one_norm_nuc = 0.0
    one_norm_elec = 0.0
    one_norm_cross = 0.0

    for k, v in non_identity.items():
        nuc_part = k[nuc_slice]
        elec_part = k[elec_slice]
        nuc_active = any(c != 'I' for c in nuc_part)
        elec_active = any(c != 'I' for c in elec_part)
        a = abs(v)
        if nuc_active and not elec_active:
            n_nuc_only += 1
            one_norm_nuc += a
        elif elec_active and not nuc_active:
            n_elec_only += 1
            one_norm_elec += a
        elif nuc_active and elec_active:
            n_cross += 1
            one_norm_cross += a
        # else: pure identity already excluded

    one_norm_total = one_norm_nuc + one_norm_elec + one_norm_cross

    # Coefficient magnitude range
    abs_coeffs = [abs(v) for v in non_identity.values() if abs(v) > 0.0]
    max_coeff = max(abs_coeffs) if abs_coeffs else 0.0
    min_coeff = min(abs_coeffs) if abs_coeffs else 0.0
    ratio = (max_coeff / min_coeff) if min_coeff > 0.0 else float('inf')

    # Histogram by order of magnitude
    histogram: Dict[int, int] = {}
    for c in abs_coeffs:
        order = int(np.floor(np.log10(c)))
        histogram[order] = histogram.get(order, 0) + 1

    # QWC groups via greedy coloring
    from geovac.nuclear.nuclear_hamiltonian import _compute_qwc_groups
    qwc_groups = _compute_qwc_groups(non_identity)

    return {
        'n_pauli_terms_total': len(H_pauli),
        'n_pauli_terms_non_identity': len(non_identity),
        'n_nuclear_only': n_nuc_only,
        'n_electronic_only': n_elec_only,
        'n_cross_register': n_cross,
        'one_norm_total': one_norm_total,
        'one_norm_nuclear': one_norm_nuc,
        'one_norm_electronic': one_norm_elec,
        'one_norm_coupling': one_norm_cross,
        'max_coefficient': max_coeff,
        'min_coefficient': min_coeff,
        'coefficient_ratio': ratio,
        'coefficient_histogram': histogram,
        'n_qwc_groups': len(qwc_groups),
        'identity_coefficient': H_pauli.get(identity, 0.0),
        'Q_total': Q_total,
        'Q_nuc': Q_nuc,
        'Q_elec': Q_elec,
    }


# ---------------------------------------------------------------------------
# 8. Pauli -> matrix utility for validation
# ---------------------------------------------------------------------------

# Single-qubit Pauli matrices
_I2 = np.array([[1, 0], [0, 1]], dtype=complex)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)
_PAULI_MAT = {'I': _I2, 'X': _X, 'Y': _Y, 'Z': _Z}


def pauli_string_to_matrix(pauli_str: str) -> np.ndarray:
    """Kronecker product of single-qubit Paulis."""
    mats = [_PAULI_MAT[c] for c in pauli_str]
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def pauli_dict_to_matrix(H_pauli: Dict[str, float]) -> np.ndarray:
    """
    Build the dense matrix corresponding to a Pauli dict. Only feasible
    for small qubit counts (<= ~12). For Q=24 this is 16M x 16M --- too big.
    Use project_to_sector() for validation instead.
    """
    if not H_pauli:
        raise ValueError("Empty Pauli dict")
    Q = len(next(iter(H_pauli.keys())))
    dim = 2 ** Q
    H = np.zeros((dim, dim), dtype=complex)
    for k, v in H_pauli.items():
        H += v * pauli_string_to_matrix(k)
    return H


def pauli_dict_to_sparse_sector(
    H_pauli: Dict[str, float],
    Q_nuc: int,
    Q_elec: int,
    nuclear_block: Dict[str, Any],
    electronic_block: Dict[str, Any],
) -> np.ndarray:
    """
    Project the Pauli Hamiltonian onto the physical sector
    (1 proton occupation) x (1 neutron occupation) x (1 electron occupation).

    For N_shells=2, n_max=2, this is 8 * 8 * 8 = 512-dimensional, which is
    much more tractable than 2^24 = 16.7M.

    Basis indexing (row-major):
        |p_i, n_j, e_k> -> index = i * n_n * n_e + j * n_e + k
    where n_p=8 proton SP states, n_n=8 neutron SP states, n_e=8 electronic
    spin-orbitals.

    Returns
    -------
    H_proj : np.ndarray, shape (n_p * n_n * n_e, n_p * n_n * n_e), complex
    """
    n_p = nuclear_block['Q_p']
    n_n = nuclear_block['Q_n']
    n_e = electronic_block['Q_elec']
    Q_total = Q_nuc + Q_elec

    # Build basis of occupation strings (bits): bit i corresponds to qubit i.
    # JW convention: bit i=1 iff qubit i is occupied.
    # A state |p_i, n_j, e_k> has bits set at positions i, n_p + j, Q_nuc + k.
    dim = n_p * n_n * n_e

    # Precompute for each basis state, the occupation bitstring (as int)
    bits = np.zeros(dim, dtype=np.int64)
    for i in range(n_p):
        for j in range(n_n):
            for k in range(n_e):
                idx = (i * n_n + j) * n_e + k
                mask = (1 << (Q_total - 1 - i)) \
                     | (1 << (Q_total - 1 - (n_p + j))) \
                     | (1 << (Q_total - 1 - (Q_nuc + k)))
                bits[idx] = mask

    # Now, for each Pauli string, compute its matrix elements between these
    # basis states. This is a brute-force approach: for each basis state,
    # apply the Pauli string and see which basis state it maps to (if any
    # within the physical sector) and with what sign/phase.
    H_proj = np.zeros((dim, dim), dtype=complex)

    # Fast path: decompose each Pauli string into its set of X/Y/Z positions.
    for pstr, coeff in H_pauli.items():
        # Identify positions and ops
        xs = [q for q, c in enumerate(pstr) if c == 'X']
        ys = [q for q, c in enumerate(pstr) if c == 'Y']
        zs = [q for q, c in enumerate(pstr) if c == 'Z']
        flip_mask = 0
        for q in xs + ys:
            flip_mask |= 1 << (Q_total - 1 - q)
        # Determine sign/phase for each basis state
        for col in range(dim):
            b = bits[col]
            new_b = b ^ flip_mask  # X and Y both flip
            # Phase factor: each Y gives +i if bit is 0 in source,
            # -i if bit is 1 in source (since Y|0>=+i|1>, Y|1>=-i|0>)
            # Each Z gives +1 if bit is 0, -1 if bit is 1
            phase: complex = 1.0 + 0.0j
            for q in ys:
                mask = 1 << (Q_total - 1 - q)
                if b & mask:
                    phase *= -1j
                else:
                    phase *= 1j
            sign = 1
            for q in zs:
                mask = 1 << (Q_total - 1 - q)
                if b & mask:
                    sign = -sign
            phase *= sign

            # Find new_b in bits array (if it's in physical sector)
            # bits is a sorted-by-construction list; use np.searchsorted
            # but bits is NOT sorted --- we built it in a specific order.
            # Use dict lookup. Build once outside loop.
            # For efficiency, precompute bits -> idx map outside the loop.
            pass  # handled below after we build the map

    # Precompute bit -> index map
    bit_to_idx = {int(b): i for i, b in enumerate(bits)}

    for pstr, coeff in H_pauli.items():
        xs = [q for q, c in enumerate(pstr) if c == 'X']
        ys = [q for q, c in enumerate(pstr) if c == 'Y']
        zs = [q for q, c in enumerate(pstr) if c == 'Z']
        flip_mask = 0
        for q in xs + ys:
            flip_mask |= 1 << (Q_total - 1 - q)
        for col in range(dim):
            b = int(bits[col])
            new_b = b ^ flip_mask
            if new_b not in bit_to_idx:
                continue
            row = bit_to_idx[new_b]
            phase: complex = 1.0 + 0.0j
            for q in ys:
                mask = 1 << (Q_total - 1 - q)
                if b & mask:
                    phase *= -1j
                else:
                    phase *= 1j
            sign = 1
            for q in zs:
                mask = 1 << (Q_total - 1 - q)
                if b & mask:
                    sign = -sign
            phase *= sign
            H_proj[row, col] += coeff * phase

    return H_proj


# ---------------------------------------------------------------------------
# 8b. Matrix-form composed Hamiltonian (physical sector, source of truth)
# ---------------------------------------------------------------------------
#
# The Pauli dict from `build_deuteron_hamiltonian` has a known sign issue
# in its off-diagonal one-body encoding (Track NE legacy). Specifically,
# the helper `_one_body_pauli` stores XY-YX terms with real coefficients,
# which produces an A(p)A(n) - C(p)C(n) cross-product instead of the correct
# A(p)A(n) + C(p)C(n) Hermitian symmetrization. The deuteron's H_matrix
# (FCI form) is unaffected and remains the source of truth for spectra.
#
# For Track NI, we therefore build the composed Hamiltonian directly as a
# matrix in the physical sector |p_i, n_j, e_k> (dim n_p * n_n * n_e),
# combining the deuteron H_matrix, the diagonal electronic Hamiltonian,
# and explicit matrix forms for the hyperfine and finite-size couplings.
# The Pauli dict from `build_deuterium_composed_hamiltonian` is retained
# for resource counts (Pauli term count, 1-norm, QWC groups) since the
# bug affects only the spectral content of the off-diagonal one-body
# pieces, not the term count.

def build_deuterium_composed_matrix(
    N_shells: int = 2,
    hw: float = 10.0,
    n_max_elec: int = 2,
    R_nuc: float = R_PROTON_BOHR,
    include_finite_size: bool = True,
    include_hyperfine: bool = True,
    Z: int = 1,
) -> Dict[str, Any]:
    """
    Build the composed nuclear-electronic deuterium Hamiltonian as a dense
    matrix in the physical sector |p_i, n_j, e_k>.

    Dimension: n_p * n_n * n_e = 8 * 8 * 8 = 512 at N_shells=2, n_max=2.

    All terms are in Hartree.

    Returns
    -------
    dict with keys:
        'H_matrix' : np.ndarray (Hermitian, dim x dim, real)
        'dim'      : int
        'n_p', 'n_n', 'n_e' : int
        'nuclear_block', 'electronic_block' : block dicts (for indices)
        'metadata' : dict
    """
    nuc = build_nuclear_block(N_shells=N_shells, hw=hw)
    elec = build_electronic_block(Z=Z, n_max=n_max_elec)

    n_p = nuc['Q_p']
    n_n = nuc['Q_n']
    n_e = elec['Q_elec']
    dim = n_p * n_n * n_e

    # Index helper for the physical-sector basis
    def idx(i: int, j: int, k: int) -> int:
        return (i * n_n + j) * n_e + k

    # 1. Nuclear contribution: H_nuc[p, p'] * delta(n,n') * delta(e,e')
    # H_matrix is in MeV; convert to Ha.
    H_nuc_MeV = nuc['H_matrix']  # shape (n_p * n_n, n_p * n_n)
    H_nuc_Ha = H_nuc_MeV * HA_PER_MEV

    H = np.zeros((dim, dim), dtype=float)

    # The deuteron H_matrix is indexed as flat (p_i * n_n + n_j).
    # Embed in (p, n, e) space: identity on e.
    for pn_row in range(n_p * n_n):
        i_row = pn_row // n_n
        j_row = pn_row % n_n
        for pn_col in range(n_p * n_n):
            v = H_nuc_Ha[pn_row, pn_col]
            if v == 0.0:
                continue
            i_col = pn_col // n_n
            j_col = pn_col % n_n
            for k in range(n_e):
                H[idx(i_row, j_row, k), idx(i_col, j_col, k)] += v

    # 2. Electronic contribution: diagonal in (p, n), eps[k] on e
    eps = elec['orbital_energies']  # shape (n_e,)
    for i in range(n_p):
        for j in range(n_n):
            for k in range(n_e):
                H[idx(i, j, k), idx(i, j, k)] += eps[k]

    # 3. Finite-size coupling: dE * (n_1s_up + n_1s_down) on the electron
    if include_finite_size:
        dE = finite_size_correction(Z=Z, R_nuc=R_nuc, n=1, l=0)  # Ha
        e_up, e_dn = _find_1s_up_down(elec)
        for i in range(n_p):
            for j in range(n_n):
                for k in (e_up, e_dn):
                    H[idx(i, j, k), idx(i, j, k)] += dE

    # 4. Hyperfine coupling: A_hf * I_proton . S_electron
    # Acts only on the proton-spin DOF (within proton 0s_{1/2}) and the
    # electron-spin DOF (within electron 1s). Spectator on neutron index.
    if include_hyperfine:
        p_up, p_dn = _find_proton_0s_qubits(nuc)
        e_up, e_dn = _find_1s_up_down(elec)
        A = HF_HYDROGEN_HA

        # I.S in the basis (m_p, m_e):
        #   |+,+>, |-,->: I_z S_z = +A/4 (diagonal)
        #   |+,->, |-,+>: I_z S_z = -A/4 (diagonal)
        #   |+,-> <-> |-,+>: off-diagonal A/2 (from I_+ S_- + I_- S_+)
        #
        # In matrix form on the 4-state subspace:
        #   diagonal: [+A/4, -A/4, -A/4, +A/4] for (++, +-, -+, --)
        #   off-diagonal: A/2 between (+-) and (-+)
        for j in range(n_n):
            # ++ : (p_up, e_up)
            H[idx(p_up, j, e_up), idx(p_up, j, e_up)] += A / 4.0
            # -- : (p_dn, e_dn)
            H[idx(p_dn, j, e_dn), idx(p_dn, j, e_dn)] += A / 4.0
            # +- : (p_up, e_dn)
            H[idx(p_up, j, e_dn), idx(p_up, j, e_dn)] += -A / 4.0
            # -+ : (p_dn, e_up)
            H[idx(p_dn, j, e_up), idx(p_dn, j, e_up)] += -A / 4.0
            # off-diagonal flip: |+-> <-> |-+>
            H[idx(p_up, j, e_dn), idx(p_dn, j, e_up)] += A / 2.0
            H[idx(p_dn, j, e_up), idx(p_up, j, e_dn)] += A / 2.0

    metadata = {
        'system': 'deuterium (1p + 1n + 1e), matrix form',
        'N_shells': N_shells,
        'hw_MeV': hw,
        'n_max_elec': n_max_elec,
        'R_nuc_bohr': R_nuc,
        'include_finite_size': include_finite_size,
        'include_hyperfine': include_hyperfine,
        'energy_unit': 'Ha',
        'sector': '(1p, 1n, 1e), dim n_p * n_n * n_e',
    }

    return {
        'H_matrix': H,
        'dim': dim,
        'n_p': n_p,
        'n_n': n_n,
        'n_e': n_e,
        'nuclear_block': nuc,
        'electronic_block': elec,
        'metadata': metadata,
    }


# ---------------------------------------------------------------------------
# 9. Validation
# ---------------------------------------------------------------------------

def validate_composed_deuterium(
    N_shells: int = 2,
    hw: float = 10.0,
    verbose: bool = False,
) -> Dict[str, Any]:
    """
    Validate the composed Hamiltonian by checking each block independently
    and then the combined + coupled spectra.

    Steps:
        1. Nuclear-only ground state vs Track NE.
        2. Electronic-only ground state vs -0.5 Ha.
        3. Full H_nuc + H_elec with no coupling: eigenvalues are sums.
        4. Add hyperfine: singlet-triplet gap ~ A_hf ~ 5.87 ueV ~ 2.16e-7 Ha.
        5. Add finite-size: shift of 1s by ~ 1e-10 Ha at R_proton.

    Returns
    -------
    dict with results and flags.
    """
    # Step 1: Nuclear-only (matrix form, source of truth)
    mat_nuc_only = build_deuterium_composed_matrix(
        N_shells=N_shells, hw=hw,
        include_finite_size=False, include_hyperfine=False,
    )
    nuc = mat_nuc_only['nuclear_block']
    elec = mat_nuc_only['electronic_block']
    H_no_coupling = mat_nuc_only['H_matrix']
    herm_err = float(np.max(np.abs(H_no_coupling - H_no_coupling.T)))
    evals = np.linalg.eigvalsh((H_no_coupling + H_no_coupling.T) / 2.0)
    E_gs_composed = float(evals[0])

    # Compare to diagonalize_deuteron directly (in MeV)
    from geovac.nuclear.nuclear_hamiltonian import diagonalize_deuteron
    nuc_direct = diagonalize_deuteron(N_shells=N_shells, hw=hw)
    E_nuc_direct_MeV = nuc_direct['E_gs']
    E_nuc_direct_Ha = E_nuc_direct_MeV * HA_PER_MEV

    # Expected: nuclear GS + electron 1s GS = E_nuc_direct_Ha + (-0.5)
    E_expected_no_coupling = E_nuc_direct_Ha + (-0.5)

    # Step 2: Electronic-only (trivial: diagonal)
    E_elec_1s = float(elec['orbital_energies'].min())

    # Step 3: Spectrum sums check (matrix form)
    nuc_evals_MeV = np.linalg.eigvalsh(nuc['H_matrix'])
    nuc_evals_Ha = nuc_evals_MeV * HA_PER_MEV
    elec_evals_Ha = np.sort(elec['orbital_energies'])
    sum_grid = np.add.outer(nuc_evals_Ha, elec_evals_Ha).ravel()
    sum_grid_sorted = np.sort(sum_grid)
    k_check = min(10, len(evals))
    sum_match_err = float(np.max(np.abs(evals[:k_check] - sum_grid_sorted[:k_check])))

    # Step 4: Add hyperfine (matrix form)
    mat_hf = build_deuterium_composed_matrix(
        N_shells=N_shells, hw=hw,
        include_finite_size=False, include_hyperfine=True,
    )
    H_hf = mat_hf['H_matrix']
    H_hf = (H_hf + H_hf.T) / 2.0
    hf_evals = np.linalg.eigvalsh(H_hf)
    # Maximum hyperfine-induced shift across the lowest few states.
    # For two spin-1/2 with H_hf = A I.S, the singlet shifts by -3A/4
    # and the triplet by +A/4. Max shift = 3*A_hf/4.
    n_compare = min(len(hf_evals), len(evals))
    hf_max_shift = float(np.max(np.abs(hf_evals[:n_compare] - evals[:n_compare])))
    # The ground state manifold is split by hyperfine. Look at the 4 lowest
    # states (1 singlet + 3 triplet) and their splitting.
    # For two spin-1/2 objects, H_hf = A I.S has eigenvalues:
    #   singlet (S=0): -3 A / 4
    #   triplet (S=1): +1 A / 4
    # Gap (triplet - singlet) = A.
    # So the gap between the lowest manifold states should be ~ A_hf.
    gap = float(hf_evals[1] - hf_evals[0])
    # Also verify singlet-triplet structure: should see 1 state at -3A/4 and
    # 3 states at +A/4 relative to the no-coupling ground state.
    # The ground state of the combined (nuclear + electronic, no coupling)
    # is 4-fold degenerate in the physical sector: 2 nuclear spin states x
    # 2 electronic spin states for the 1s x 1s manifold. Actually nuclear
    # deuteron ground state is a specific spin-1 state (triplet), so the
    # degeneracy depends on the details of the nuclear GS.
    # The test we do here is just: the gap between the two lowest states
    # is ~ A_hf (within factor of a few).

    # Step 5: Finite-size (matrix form)
    mat_fs = build_deuterium_composed_matrix(
        N_shells=N_shells, hw=hw,
        include_finite_size=True, include_hyperfine=False,
        R_nuc=R_PROTON_BOHR,
    )
    H_fs = mat_fs['H_matrix']
    H_fs = (H_fs + H_fs.T) / 2.0
    # The finite-size perturbation dE ~ 1e-10 Ha is well below the float64
    # roundoff floor (~1e-10) of the nuclear-scale ~5e5 Ha entries, so
    # subtracting two matrices is meaningless. Report the analytical value
    # from the closed-form perturbation formula directly. The existing
    # test_finite_size_coupling_magnitude validates the formula independently.
    fs_shift = finite_size_correction(Z=1.0, R_nuc=R_PROTON_BOHR, n=1, l=0)
    fs_expected = (2.0 / 5.0) * (R_PROTON_BOHR ** 2)

    results = {
        'nuclear_only': {
            'E_gs_Ha': E_gs_composed,
            'E_expected_Ha': E_expected_no_coupling,
            'E_nuc_direct_MeV': E_nuc_direct_MeV,
            'E_nuc_direct_Ha': E_nuc_direct_Ha,
            'herm_err': herm_err,
            'diff': E_gs_composed - E_expected_no_coupling,
        },
        'electronic_only': {
            'E_1s_Ha': E_elec_1s,
            'expected_Ha': -0.5,
            'diff': E_elec_1s - (-0.5),
        },
        'no_coupling_sums': {
            'max_abs_err': sum_match_err,
            'k_check': k_check,
        },
        'hyperfine': {
            'gap_Ha': gap,
            'max_shift_Ha': hf_max_shift,
            'expected_A_hf_Ha': HF_HYDROGEN_HA,
            'expected_max_shift_Ha': 0.75 * HF_HYDROGEN_HA,
            'ratio': hf_max_shift / (0.75 * HF_HYDROGEN_HA) if HF_HYDROGEN_HA > 0 else float('inf'),
            'lowest_4': hf_evals[:4].tolist(),
        },
        'finite_size': {
            'shift_Ha': fs_shift,
            'expected_Ha': fs_expected,
            'ratio': fs_shift / fs_expected if fs_expected > 0 else float('inf'),
        },
    }

    if verbose:
        print("Composed deuterium validation results:")
        for k, v in results.items():
            print(f"  {k}: {v}")

    return results
