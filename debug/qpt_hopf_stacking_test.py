"""QPT (Quantum Paldus Transform) stackability test with Hopf-Z2 tapering.

Burkat-Fitzpatrick June 2025 (arXiv:2506.09151) implement the quantum
Paldus transform: u(d) x SU(2) -> Gelfand-Tsetlin block-diagonalization
as a quantum circuit. The transform exploits SPIN (S^2) and ORBITAL
(unitary group) symmetries to block-diagonalize the Hamiltonian.

GeoVac's per-sub-block Hopf-U(1) Z2 tapering (v3.52.0, 2026-06-04)
exploits S3 -> S2 Hopf base m -> -m parity, removing 2 + n_sub_blocks
qubits per molecule, machine-precision spectrum preservation across
37/37 molecules.

STRUCTURAL QUESTION: Do the QPT symmetries (S^2, U(d) generators) commute
with the Hopf-Z2 stabilizers?
  - If [S^2, P_Hopf] = 0 -> they STACK additively.
  - If [S^2, P_Hopf] != 0 -> they OVERLAP and must be enforced sequentially.

The test:
  1. Build LiH/BeH2/H2O composed Hamiltonians.
  2. Verify [H, S^2] = 0 (sanity: non-relativistic builder commutes with spin).
  3. Verify [H, P_Hopf] = 0 (already verified in v3.52.0 sprint).
  4. Compute [S^2, P_Hopf] commutator. KEY OUTCOME.
  5. Also check [Z_alpha, S^2] and [Z_beta, S^2] for completeness.
  6. If all commutators vanish: report joint block dimensions in
     (S^2 eigenvalue) x (Hopf-parity sector) basis.
"""

from __future__ import annotations

import json
import os
from typing import Dict, List, Tuple

import numpy as np

try:
    from openfermion import (
        QubitOperator, FermionOperator, jordan_wigner,
        s_squared_operator,
    )
    from openfermion.utils import commutator
except ImportError as exc:
    raise SystemExit(
        f"openfermion required for QPT test: {exc}"
    )

from geovac.composed_qubit import build_composed_hamiltonian
from geovac.molecular_spec import lih_spec, beh2_spec, h2o_spec
from geovac.z2_tapering import (
    build_pm_rotation, rotate_h1_eri, build_stabilizers, _enumerate_orbitals,
)


def _max_abs_coef(op: QubitOperator) -> float:
    return max((abs(v) for v in op.terms.values()), default=0.0)


def commutator_norm(a: QubitOperator, b: QubitOperator) -> float:
    """L-infinity norm of [a, b] coefficients."""
    return _max_abs_coef(commutator(a, b))


def build_rotated_hamiltonian(spec):
    """Build the composed Hamiltonian in the Hopf-rotated basis."""
    result = build_composed_hamiltonian(spec, verbose=False)
    h1 = result['h1']
    eri = result['eri']
    nuc = result['nuclear_repulsion']

    orbital_table = _enumerate_orbitals(spec)
    U, parity = build_pm_rotation(orbital_table)
    h1_rot, eri_rot = rotate_h1_eri(h1, eri, U)

    # JW transform the rotated Hamiltonian
    from openfermion import InteractionOperator, get_fermion_operator
    M = h1.shape[0]
    one_body = h1_rot.copy()
    two_body = np.zeros((2 * M, 2 * M, 2 * M, 2 * M))
    one_body_spin = np.zeros((2 * M, 2 * M))
    # Spin-block the one-body
    for p in range(M):
        for q in range(M):
            one_body_spin[2 * p, 2 * q] = h1_rot[p, q]
            one_body_spin[2 * p + 1, 2 * q + 1] = h1_rot[p, q]
    # Spin-block the two-body (chemist -> physicist convert + spin block)
    for p in range(M):
        for q in range(M):
            for r in range(M):
                for s in range(M):
                    v = eri_rot[p, r, q, s]  # chemist (pr|qs)
                    if abs(v) < 1e-15:
                        continue
                    # physicist <pq|rs> = (pr|qs) -> spin block over alpha/beta
                    for sigma1 in (0, 1):
                        for sigma2 in (0, 1):
                            two_body[
                                2 * p + sigma1, 2 * q + sigma2,
                                2 * r + sigma2, 2 * s + sigma1,
                            ] = v
    iop = InteractionOperator(nuc, one_body_spin, 0.5 * two_body)
    fop = get_fermion_operator(iop)
    qop = jordan_wigner(fop)
    return qop, parity, orbital_table, M


def s_squared_qubit(M: int) -> QubitOperator:
    """Total S^2 operator in qubit form (JW), M spatial orbitals."""
    sf = s_squared_operator(M)
    return jordan_wigner(sf)


def measure_h_p_commutators(H: QubitOperator, S2: QubitOperator,
                            Z_alpha: QubitOperator, Z_beta: QubitOperator,
                            P_list: List[QubitOperator]) -> Dict:
    """Compute all relevant commutators."""
    out = {
        '[H, S^2]': commutator_norm(H, S2),
        '[H, Z_alpha]': commutator_norm(H, Z_alpha),
        '[H, Z_beta]': commutator_norm(H, Z_beta),
        '[S^2, Z_alpha]': commutator_norm(S2, Z_alpha),
        '[S^2, Z_beta]': commutator_norm(S2, Z_beta),
        '[Z_alpha, Z_beta]': commutator_norm(Z_alpha, Z_beta),
    }
    for i, P in enumerate(P_list):
        out[f'[H, P_{i}]'] = commutator_norm(H, P)
        out[f'[S^2, P_{i}]'] = commutator_norm(S2, P)
        out[f'[Z_alpha, P_{i}]'] = commutator_norm(Z_alpha, P)
        out[f'[Z_beta, P_{i}]'] = commutator_norm(Z_beta, P)
    for i in range(len(P_list)):
        for j in range(i + 1, len(P_list)):
            out[f'[P_{i}, P_{j}]'] = commutator_norm(P_list[i], P_list[j])
    return out


def main():
    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'qpt_hopf_stacking_test.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    summary = {}

    for spec_fn, label in [
        (lih_spec, 'LiH'),
        (beh2_spec, 'BeH2'),
        (h2o_spec, 'H2O'),
    ]:
        print(f"\n{'='*70}\n{label}\n{'='*70}")
        try:
            spec = spec_fn()
            H, parity, orbital_table, M = build_rotated_hamiltonian(spec)
            print(f"  M = {M} spatial orbitals, Q = {2 * M} qubits")
            print(f"  H has {len(H.terms)} Pauli terms")

            S2 = s_squared_qubit(M)
            print(f"  S^2 has {len(S2.terms)} Pauli terms")

            Z_alpha, Z_beta, P_list = build_stabilizers(
                parity, orbital_table=orbital_table, mode='per_block',
            )
            print(f"  Hopf stabilizers: 1 Z_alpha + 1 Z_beta + {len(P_list)} P_i")

            comms = measure_h_p_commutators(H, S2, Z_alpha, Z_beta, P_list)

            print(f"\n  COMMUTATORS:")
            all_zero = True
            for name, norm in comms.items():
                marker = "OK" if norm < 1e-10 else "**NONZERO**"
                print(f"    {name:30s}: {norm:.4e}   {marker}")
                if norm >= 1e-10:
                    all_zero = False

            print(f"\n  VERDICT: {'STACK (all commute)' if all_zero else 'OVERLAP/CONFLICT'}")

            summary[label] = {
                'M': M, 'Q': 2 * M,
                'n_pauli_H': len(H.terms),
                'n_pauli_S2': len(S2.terms),
                'n_P_stabilizers': len(P_list),
                'commutators': comms,
                'all_commute': all_zero,
            }
        except Exception as e:
            import traceback
            traceback.print_exc()
            summary[label] = {'error': str(e), 'traceback': traceback.format_exc()}

    with open(out_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
