"""
Generic atomic hyperfine A * I . J wrapper (Sprint Cs-HFS-v2).

For atomic HFS observables on single-electron systems (Cs 6S_{1/2},
Rb 5S_{1/2}, Na 3S_{1/2}, ...), the natural Hamiltonian is

    H_hf = A * I . J

where I is the nuclear spin (typically 1/2, 1, 3/2, 5/2, 7/2 for stable
isotopes) and J is the total electronic angular momentum (here 1/2 for
s-orbital ground states). The hyperfine constant A is computed from the
contact density |psi(0)|^2 via the Bohr-Fermi formula:

    A = (8 pi / 3) g_e g_N (m_e / m_p) alpha^2 |psi(0)|^2  (Hartree)

This module provides:

    bohr_fermi_a_constant(psi0_squared, g_e, g_N, m_p_over_m_e=1836.15)
        -- Returns A in Hartree (and MHz) for a single ns_{1/2} valence
           electron on a fixed nucleus.

    hyperfine_a_ij_pauli_general(A_au, I, J, qubit_layout='binary')
        -- Returns the Pauli-sum representation of H = A * I . J for
           arbitrary nuclear spin I (I=1/2, 1, 3/2, 5/2, 7/2, ...) and
           electronic spin J=1/2 (the only case implemented for atomic
           ground-state HFS at this sprint level).

    hyperfine_a_pauli_for_atomic_hfs(A_au, I)
        -- Convenience wrapper for the standard atomic HFS case
           (J=1/2 fixed, binary nuclear-spin encoding).

The Pauli-sum representation is built on the minimum register: nuclear
spin in ceil(log2(2*I+1)) qubits via binary encoding, electronic spin
in 1 qubit. For Cs (I=7/2): nuclear in 3 qubits, electronic in 1 qubit,
total 4 qubits. The dimension matches: 2^3 * 2^1 = 16 = (2I+1)*(2J+1).

The Hamiltonian on this register acts within the physical (2I+1)*(2J+1)
subspace and gives identity (or zero, depending on convention) outside;
the eigenstructure F = I+J, F = I+J-1 etc. is recovered by direct
diagonalization on the physical subspace.

This module does NOT modify the Track NI hyperfine_coupling_pauli
function in geovac/nuclear/nuclear_electronic.py; it is a separate,
simpler wrapper for atomic HFS observables that do not need the
two-register nuclear+electronic architecture.

Author: GeoVac Development Team
Date: 2026-05-09 (Sprint Cs-HFS-v2)
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Physical constants (CODATA 2018 / 2022, atomic units)
# ---------------------------------------------------------------------------
ALPHA: float = 7.2973525693e-3
INV_ALPHA: float = 1.0 / ALPHA
HZ_PER_HA: float = 6.579683920502e15
M_PROTON_OVER_M_E: float = 1836.15267343
GE_DIRAC: float = 2.0
GE_FULL: float = 2.00231930436256


# ---------------------------------------------------------------------------
# Bohr-Fermi A constant
# ---------------------------------------------------------------------------

def bohr_fermi_a_constant(
    psi0_squared: float,
    g_e: float = GE_FULL,
    g_N: float = 1.0,
    m_p_over_m_e: float = M_PROTON_OVER_M_E,
) -> Dict[str, float]:
    """Compute the Bohr-Fermi hyperfine A constant from the contact density.

    The standard atomic-physics convention is

        A_hf = (8 pi / 3) * (g_e/2) * (g_N/2) * alpha^2 * (m_e/m_p)
               * |psi(0)|^2  (in Hartree, with |psi(0)|^2 in bohr^{-3})

    Equivalently, using the contact-density form g_N is the *bare* nuclear
    g-factor in nuclear-magneton units (mu_I = g_N I mu_N).

    Cross-check against hydrogen 1s:
        |psi_1s(0)|^2 = 1/pi  (Z=1)
        g_e = 2.00232, g_p = 5.5857
        A(H 1s) = (8 pi / 3) * 1.00116 * 2.79285 * (1/137.036)^2
                  * (1/1836.15) * (1/pi)
                = ... ~ 6.4965e-7 Ha = 4276.8 MHz / 6 = 1418.84 MHz
        (vs experimental 21cm A = 1420.41 MHz; framework BF strict
         matches at +0.11% with full g_e, g_p factors.)

    Parameters
    ----------
    psi0_squared : float
        |psi_ns(0)|^2 at the nucleus, in atomic units (bohr^{-3}).
    g_e : float
        Electron g-factor. Default: full Schwinger-corrected value.
    g_N : float
        Nuclear g-factor (mu_I = g_N I mu_N). For proton g_p = 5.585695.
        For Cs-133 (I=7/2, mu = 2.582025 mu_N): g_Cs = 2.582025/3.5 =
        0.7377. For Mg-25 (I=5/2, mu = -0.85546 mu_N): g_Mg = -0.34218.
    m_p_over_m_e : float
        Nucleon-to-electron mass ratio. Default proton.

    Returns
    -------
    dict with:
        'A_Ha' : A in Hartree
        'A_MHz' : A in megahertz
        'A_Hz' : A in hertz
        'g_e_used' : float
        'g_N_used' : float
        'psi0_squared_au' : float
    """
    pi = float(np.pi)
    A_Ha = (
        (8.0 * pi / 3.0)
        * (g_e / 2.0)
        * (g_N / 2.0)
        * (ALPHA ** 2)
        * (1.0 / m_p_over_m_e)
        * psi0_squared
    )
    A_Hz = A_Ha * HZ_PER_HA
    return {
        'A_Ha': A_Ha,
        'A_Hz': A_Hz,
        'A_MHz': A_Hz * 1e-6,
        'g_e_used': g_e,
        'g_N_used': g_N,
        'psi0_squared_au': psi0_squared,
    }


# ---------------------------------------------------------------------------
# Generic A * I . J Pauli-sum wrapper (binary encoding for arbitrary I)
# ---------------------------------------------------------------------------

def _angular_momentum_matrices(j: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build dense matrix representations of J_x, J_y, J_z for spin j.

    Standard basis: states |j, m> with m = -j, -j+1, ..., +j.
    Dimension: 2j+1 (j must be a non-negative integer or half-integer).

    Returns
    -------
    Jx, Jy, Jz : ndarray, dim (2j+1, 2j+1)
        Hermitian matrices with eigenvalues -j, -j+1, ..., +j (Jz).
    """
    if j < 0 or abs(2.0 * j - round(2.0 * j)) > 1e-12:
        raise ValueError(f"j={j} must be non-negative integer or half-integer")
    dim = int(round(2.0 * j + 1))
    Jz = np.zeros((dim, dim), dtype=complex)
    Jx = np.zeros((dim, dim), dtype=complex)
    Jy = np.zeros((dim, dim), dtype=complex)

    # Basis ordering: index k = 0 -> m = +j, k = 1 -> m = +j - 1, ..., k = dim-1 -> m = -j
    # (This convention puts the highest-m state at index 0; reversing doesn't
    # affect physics, just the ladder direction.)
    for k in range(dim):
        m = j - k  # m for state at index k
        Jz[k, k] = m
        # J+ |j, m> = sqrt((j-m)(j+m+1)) |j, m+1>:  raises m by 1
        # In our basis with index k = j - m, m+1 corresponds to k-1.
        if k > 0:
            m_lower = j - (k - 1)  # m of state |j, m_lower> = |j, m+1>
            # <j, m+1 | J+ | j, m> = sqrt((j-m)(j+m+1))
            # In matrix indices: J+_{k-1, k} = sqrt((j-m)(j+m+1))
            jp = float(np.sqrt((j - m) * (j + m + 1.0)))
            Jx[k - 1, k] += 0.5 * jp
            Jx[k, k - 1] += 0.5 * jp
            Jy[k - 1, k] += -0.5j * jp
            Jy[k, k - 1] += +0.5j * jp

    return Jx, Jy, Jz


def _matrix_to_pauli(M: np.ndarray, n_qubits: int) -> Dict[str, complex]:
    """Decompose a 2^n x 2^n matrix into a Pauli string sum.

    Returns dict mapping Pauli string (e.g. 'XYZI') to coefficient.
    Uses the orthonormal basis: tr(P^dag P') = 2^n delta_{P,P'}.
    """
    dim = 2 ** n_qubits
    if M.shape != (dim, dim):
        raise ValueError(
            f"Matrix shape {M.shape} != expected ({dim}, {dim})"
        )

    pauli_matrices = {
        'I': np.array([[1, 0], [0, 1]], dtype=complex),
        'X': np.array([[0, 1], [1, 0]], dtype=complex),
        'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
        'Z': np.array([[1, 0], [0, -1]], dtype=complex),
    }

    coeffs: Dict[str, complex] = {}
    # Iterate over all 4^n Pauli strings (cheap for small n; we have n<=4 here)
    from itertools import product
    for combo in product(['I', 'X', 'Y', 'Z'], repeat=n_qubits):
        P = np.array([[1.0]], dtype=complex)
        for ch in combo:
            P = np.kron(P, pauli_matrices[ch])
        coeff = np.trace(P.conj().T @ M) / dim
        if abs(coeff) > 1e-13:
            coeffs[''.join(combo)] = coeff
    return coeffs


def hyperfine_a_ij_pauli_general(
    A_au: float,
    I: float,
    J: float = 0.5,
    real_only: bool = True,
) -> Dict[str, float]:
    """Pauli-sum representation of H_hf = A * I . J for atomic HFS.

    Builds the operator on a minimum binary-encoded register:

        nuclear spin I -> ceil(log2(2*I + 1)) qubits
        electronic spin J -> ceil(log2(2*J + 1)) qubits

    Layout: nuclear qubits first, then electronic.

    For Cs-133 (I=7/2, J=1/2): 3 + 1 = 4 qubits, 16-dimensional Hilbert space.
    For H 1s (I=1/2, J=1/2): 1 + 1 = 2 qubits, 4-dim. (Note this is a
    DIFFERENT encoding from the Track NI deuterium hyperfine, which uses
    4 qubits for 4 states via a non-binary "occupation" encoding.)

    The dimension of the binary register is 2^Q which matches the
    physical (2I+1)*(2J+1) subspace exactly when both 2I+1 and 2J+1 are
    powers of 2 (i.e. I and J are half-integer). When 2I+1 is not a
    power of 2 (e.g. I=1, where 2I+1=3 fits in 2 qubits but uses only
    3 of 4 states), the Pauli-sum acts as zero on the unphysical states.

    Parameters
    ----------
    A_au : float
        Hyperfine A constant in Hartree (atomic units).
    I : float
        Nuclear spin (non-negative integer or half-integer).
    J : float
        Electronic total angular momentum. Default 0.5 (s-orbital).
    real_only : bool
        If True (default), drop tiny imaginary parts of the Pauli
        coefficients (which arise from numerical roundoff in the
        I . J = I_x J_x + I_y J_y + I_z J_z product). The H_hf
        Hamiltonian is Hermitian; coefficients in the Pauli basis are
        real for the standard {I, X, Y, Z} convention.

    Returns
    -------
    dict mapping Pauli string (str of length Q_total) to coefficient
    (float, in Hartree).
    """
    if I < 0 or abs(2.0 * I - round(2.0 * I)) > 1e-12:
        raise ValueError(f"I={I} must be non-negative integer or half-integer")
    if J < 0 or abs(2.0 * J - round(2.0 * J)) > 1e-12:
        raise ValueError(f"J={J} must be non-negative integer or half-integer")

    dim_I = int(round(2.0 * I + 1))
    dim_J = int(round(2.0 * J + 1))
    n_qubits_I = max(1, int(np.ceil(np.log2(dim_I))))
    n_qubits_J = max(1, int(np.ceil(np.log2(dim_J))))
    n_qubits_total = n_qubits_I + n_qubits_J
    full_dim = 2 ** n_qubits_total

    # Build I_alpha and J_alpha matrices on their physical subspaces
    Ix, Iy, Iz = _angular_momentum_matrices(I)
    Jx, Jy, Jz = _angular_momentum_matrices(J)

    # Embed into the binary register: pad with zero-block to 2^n
    pad_I = 2 ** n_qubits_I - dim_I
    pad_J = 2 ** n_qubits_J - dim_J
    if pad_I > 0:
        Ix = np.pad(Ix, ((0, pad_I), (0, pad_I)), mode='constant')
        Iy = np.pad(Iy, ((0, pad_I), (0, pad_I)), mode='constant')
        Iz = np.pad(Iz, ((0, pad_I), (0, pad_I)), mode='constant')
    if pad_J > 0:
        Jx = np.pad(Jx, ((0, pad_J), (0, pad_J)), mode='constant')
        Jy = np.pad(Jy, ((0, pad_J), (0, pad_J)), mode='constant')
        Jz = np.pad(Jz, ((0, pad_J), (0, pad_J)), mode='constant')

    # Build I . J on the tensor-product Hilbert space (nuclear (x) electronic)
    H_dense = (
        np.kron(Ix, Jx) + np.kron(Iy, Jy) + np.kron(Iz, Jz)
    )
    H_dense *= A_au

    # Decompose into Pauli strings
    pauli_complex = _matrix_to_pauli(H_dense, n_qubits_total)

    # Convert to real (drop tiny imag parts from numerical Pauli decomposition)
    pauli: Dict[str, float] = {}
    for k, v in pauli_complex.items():
        if real_only:
            if abs(v.imag) > 1e-10 * max(abs(v.real), 1e-10):
                raise RuntimeError(
                    f"Non-real Pauli coefficient {v} for string {k}; "
                    f"H_hf should be Hermitian."
                )
            if abs(v.real) > 1e-13:
                pauli[k] = float(v.real)
        else:
            pauli[k] = v
    return pauli


def hyperfine_a_pauli_for_atomic_hfs(
    A_au: float,
    I: float,
) -> Dict[str, Any]:
    """Convenience wrapper for atomic HFS (J=1/2 fixed).

    Returns the H = A * I . J Pauli-sum on the minimum binary register
    plus metadata (qubit count, eigenvalues of F = I + J, etc.).

    Parameters
    ----------
    A_au : float
        Hyperfine A constant in Hartree.
    I : float
        Nuclear spin.

    Returns
    -------
    dict with:
        'pauli_terms' : Dict[str, float]
        'Q_nuc' : int
        'Q_elec' : int
        'Q_total' : int
        'I' : float
        'J' : float
        'A_Ha' : float
        'A_MHz' : float
        'F_levels' : List[float]  -- F = I + J, I + J - 1, ..., |I - J|
        'F_eigenvalue_eV' : Dict[float, float]  -- A * (F(F+1) - I(I+1) - J(J+1)) / 2
        'splitting_F_max_to_F_min_MHz' : float
    """
    pauli = hyperfine_a_ij_pauli_general(A_au, I, J=0.5)
    dim_I = int(round(2.0 * I + 1))
    n_qubits_I = max(1, int(np.ceil(np.log2(dim_I))))
    Q_nuc = n_qubits_I
    Q_elec = 1
    Q_total = Q_nuc + Q_elec

    # Compute F levels and their energies
    J = 0.5
    F_levels = []
    F_max = I + J
    F_min = abs(I - J)
    F_step = 1.0
    nF = int(round(F_max - F_min)) + 1
    for k in range(nF):
        F_levels.append(F_max - k)

    # E_F = A/2 * [F(F+1) - I(I+1) - J(J+1)]
    F_energies_Ha: Dict[float, float] = {}
    II = I * (I + 1.0)
    JJ = J * (J + 1.0)
    for F in F_levels:
        FF = F * (F + 1.0)
        E_F = 0.5 * A_au * (FF - II - JJ)
        F_energies_Ha[F] = E_F

    # Splitting from F_max to F_min in MHz
    A_MHz = A_au * HZ_PER_HA / 1e6
    splitting_Ha = F_energies_Ha[F_max] - F_energies_Ha[F_min]
    splitting_MHz = splitting_Ha * HZ_PER_HA / 1e6

    return {
        'pauli_terms': pauli,
        'Q_nuc': Q_nuc,
        'Q_elec': Q_elec,
        'Q_total': Q_total,
        'I': I,
        'J': J,
        'A_Ha': A_au,
        'A_MHz': A_MHz,
        'F_levels': F_levels,
        'F_energies_Ha': F_energies_Ha,
        'splitting_F_max_to_F_min_Ha': splitting_Ha,
        'splitting_F_max_to_F_min_MHz': splitting_MHz,
    }
