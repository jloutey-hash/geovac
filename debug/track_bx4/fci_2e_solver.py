"""
Correct 2-electron FCI solver using direct second-quantization algebra.
Works for non-Hermitian (TC) Hamiltonians.
"""
import numpy as np
from itertools import combinations
from typing import Dict, List, Tuple


def apply_op_string(
    ops: Tuple[Tuple[int, int], ...],
    det: Tuple[int, ...],
) -> Tuple[float, Tuple[int, ...]]:
    """
    Apply a string of creation/annihilation operators to a Fock state.

    ops: tuple of (index, type) where type=1 for creation, 0 for annihilation.
         Applied RIGHT-TO-LEFT (last element first).
    det: sorted tuple of occupied spin-orbital indices.

    Returns (phase, new_det) or (0.0, ()) if result is zero.
    """
    state = list(det)
    phase = 1.0

    # Apply operators right-to-left
    for op_idx, op_type in reversed(ops):
        if op_type == 0:  # annihilation
            if op_idx not in state:
                return 0.0, ()
            pos = state.index(op_idx)
            phase *= (-1) ** pos
            state.pop(pos)
        else:  # creation
            if op_idx in state:
                return 0.0, ()
            # Find insertion position (maintain sorted order)
            pos = 0
            while pos < len(state) and state[pos] < op_idx:
                pos += 1
            phase *= (-1) ** pos
            state.insert(pos, op_idx)

    return phase, tuple(state)


def build_fci_matrix_2e(
    fermion_op,
    n_qubits: int,
) -> np.ndarray:
    """
    Build FCI matrix in the 2-electron sector from a FermionOperator.

    This is correct for non-Hermitian operators.
    """
    n_so = n_qubits
    dets = list(combinations(range(n_so), 2))
    N_SD = len(dets)
    det_map = {d: i for i, d in enumerate(dets)}

    H = np.zeros((N_SD, N_SD), dtype=complex)

    for term, coeff in fermion_op.terms.items():
        if abs(coeff) < 1e-15:
            continue

        if len(term) == 0:
            # Identity: contributes to diagonal
            for I in range(N_SD):
                H[I, I] += coeff
            continue

        # For each ket determinant, apply the operator string
        for J in range(N_SD):
            phase, new_det = apply_op_string(term, dets[J])
            if phase == 0.0:
                continue
            if new_det in det_map:
                I = det_map[new_det]
                H[I, J] += coeff * phase

    return H


def solve_fci_2e(
    h1: np.ndarray,
    eri: np.ndarray,
    nuclear_repulsion: float,
    n_qubits: int,
) -> Tuple[np.ndarray, int]:
    """
    Build and solve FCI in the 2-electron sector.

    Parameters
    ----------
    h1 : (M, M) one-body integrals (spatial)
    eri : (M, M, M, M) two-body integrals in chemist notation (spatial)
    nuclear_repulsion : float
    n_qubits : int (= 2*M)

    Returns
    -------
    eigenvalues : np.ndarray (sorted by real part)
    N_SD : int
    """
    from geovac.qubit_encoding import build_fermion_op_from_integrals

    fermion_op = build_fermion_op_from_integrals(h1, eri, nuclear_repulsion)
    H = build_fci_matrix_2e(fermion_op, n_qubits)

    vals = np.linalg.eigvals(H)
    return np.sort(vals.real), H.shape[0]


if __name__ == "__main__":
    # Quick test
    from geovac.tc_integrals import compute_tc_integrals_block, tc_eri_to_chemist

    EXACT_HE = -2.9037
    Z = 2.0

    for max_n in [1, 2]:
        states = []
        for n in range(1, max_n + 1):
            for l in range(n):
                for m in range(-l, l + 1):
                    states.append((n, l, m))
        M = len(states)
        Q = 2 * M

        h1 = np.zeros((M, M))
        for i, (n, l, m) in enumerate(states):
            h1[i, i] = -Z**2 / (2.0 * n**2)

        for ang in [False, True]:
            tc_eri = compute_tc_integrals_block(Z, states, 2000, include_angular=ang)
            eri = tc_eri_to_chemist(tc_eri, M)
            vals, N_SD = solve_fci_2e(h1, eri, -0.25, Q)
            E = vals[0]
            err = abs(E - EXACT_HE) / abs(EXACT_HE) * 100
            label = "full" if ang else "rad "
            print(f"max_n={max_n} {label}: E={E:.6f}, err={err:.3f}%, N_SD={N_SD}")
