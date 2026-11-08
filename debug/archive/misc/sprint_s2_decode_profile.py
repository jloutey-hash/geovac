"""S2 — Decode the universal {4,16,16,9,9,9,6,3,3} interior bond-rank profile.

For each cut k inside a 5-orbital sub-block ({1s, 2s, 2p_{-1,0,+1}} -> 10
qubits under JW with alpha/beta interleaving), enumerate which Pauli
strings contribute on the LEFT side and count their distinct types.

This gives a structural reading of where each rank number comes from.

JW qubit ordering for one 5-orbital sub-block (alpha/beta interleaved):
  qubit 0: 1s_alpha
  qubit 1: 1s_beta
  qubit 2: 2s_alpha
  qubit 3: 2s_beta
  qubit 4: 2p_{-1}_alpha
  qubit 5: 2p_{-1}_beta
  qubit 6: 2p_0_alpha
  qubit 7: 2p_0_beta
  qubit 8: 2p_{+1}_alpha
  qubit 9: 2p_{+1}_beta

Expected interior cuts (within one sub-block):
  cut 1: chi = 4   (left has 1s_alpha only -- 1 qubit, 4 Pauli types)
  cut 2: chi = 16  (left has 1s, both spins -- 2 qubits, saturate Pauli)
  cut 3: chi = 16  (left has 1s + 2s_alpha -- 3 qubits, NOT saturated)
  cut 4: chi = 9   (left has 1s + 2s -- 2s pair closed)
  cuts 5-6: chi = 9 (active 2p_{-1} pair)
  cut 7: chi = 6   (2p_{-1} pair closed, 2p_0 active)
  cuts 8-9: chi = 3 (2p_{-1, 0} closed, only 2p_{+1} active)
  cut 10: chi = 2  (boundary, all 1s/2s/2p closed)

The pattern suggests chi tracks the "number of active orbital pairs"
that have been opened but not yet closed. Let's verify by enumerating.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Tuple

import numpy as np
from openfermion import QubitOperator

from geovac.composed_qubit import build_composed_hamiltonian
from geovac.molecular_spec import lih_spec


def enumerate_left_paulis_at_cut(
    qop: QubitOperator, cut: int,
) -> Dict[Tuple, List[float]]:
    """For each distinct LEFT Pauli string at cut k, list the coefficients
    of all Pauli terms that have it on the left."""
    left_paulis: Dict[Tuple, List[float]] = defaultdict(list)
    for pauli_tuple, c in qop.terms.items():
        left = tuple((q, op) for (q, op) in pauli_tuple if q < cut)
        left_paulis[left].append(float(c.real if hasattr(c, 'real') else c))
    return left_paulis


def classify_left_pauli(left_pauli: Tuple, cut: int) -> Dict:
    """Classify a Pauli string by its 'shape' — qubit support and Pauli type."""
    if not left_pauli:
        return {'support': [], 'pattern': 'I', 'paulis': [], 'n_X': 0, 'n_Y': 0,
                'n_Z': 0, 'n_fermion_ops': 0}
    qubits = [q for q, _ in left_pauli]
    paulis = [op for _, op in left_pauli]
    # Count Pauli types
    n_X = paulis.count('X')
    n_Y = paulis.count('Y')
    n_Z = paulis.count('Z')
    # Identify spin parity (X/Y count = fermion operator support; Z = number)
    n_fermion_ops = n_X + n_Y
    return {
        'support': qubits,
        'paulis': paulis,
        'n_X': n_X,
        'n_Y': n_Y,
        'n_Z': n_Z,
        'n_fermion_ops': n_fermion_ops,
        'pattern': ''.join(paulis),
    }


def main():
    # Build LiH for the Pauli operator
    spec = lih_spec()
    res = build_composed_hamiltonian(spec, verbose=False)
    qop = res['qubit_op']
    print(f"LiH composed: {len(qop.terms)} Pauli terms")

    # Analyze cuts 1 through 9 (interior of first sub-block)
    print(f"\n{'='*70}")
    print(f"Left-Pauli enumeration at each interior cut of sub-block 1")
    print(f"{'='*70}")

    for cut in range(1, 11):
        left_paulis = enumerate_left_paulis_at_cut(qop, cut)
        n_distinct = len(left_paulis)

        # For each distinct left Pauli, classify
        # Group by (n_fermion_ops, n_Z) to see structure
        by_shape: Dict[Tuple, int] = defaultdict(int)
        for left_p, _coefs in left_paulis.items():
            cls = classify_left_pauli(left_p, cut)
            shape_key = (cls['n_fermion_ops'], cls['n_Z'])
            by_shape[shape_key] += 1

        print(f"\nCut {cut:2d}: {n_distinct} distinct left Pauli strings")
        for shape, count in sorted(by_shape.items()):
            n_fermion, n_Z = shape
            print(f"   {count:3d} strings with {n_fermion} X/Y + {n_Z} Z operators")

    # Now check: is the rank ALSO 16 at cut 2 because we saturate or because
    # the coefficient matrix has full rank?
    print(f"\n{'='*70}")
    print(f"Coefficient matrix rank at each cut")
    print(f"{'='*70}")

    for cut in range(1, 11):
        coef = {}
        for pauli_tuple, c in qop.terms.items():
            left = tuple((q, op) for (q, op) in pauli_tuple if q < cut)
            right = tuple((q, op) for (q, op) in pauli_tuple if q >= cut)
            coef[(left, right)] = float(c.real if hasattr(c, 'real') else c)

        left_keys = sorted(set(k[0] for k in coef))
        right_keys = sorted(set(k[1] for k in coef))
        if not left_keys or not right_keys:
            continue

        left_idx = {k: i for i, k in enumerate(left_keys)}
        right_idx = {k: i for i, k in enumerate(right_keys)}
        M = np.zeros((len(left_keys), len(right_keys)))
        for (l, r), v in coef.items():
            M[left_idx[l], right_idx[r]] = v
        rank = int(np.linalg.matrix_rank(M))
        print(f"  cut {cut:2d}: {len(left_keys):4d} left x {len(right_keys):4d} right, "
              f"matrix rank = {rank}")


if __name__ == '__main__':
    main()
