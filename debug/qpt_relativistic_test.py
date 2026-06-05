"""QPT compatibility test on the relativistic Dirac-on-S^3 Tier 2 builder.

The non-relativistic QPT test confirmed [H, S^2] = [H, P_Hopf] = [S^2, P_Hopf] = 0
bit-exactly. With spin-orbit coupling (Tier 2 relativistic), [H_rel, S^2] is
nonzero in general because spin-orbit mixes (l, s) couplings. The conserved
quantum number is total J (with m_j eigenvalue). QPT-style adaptation must
use J^2 instead of S^2 in the relativistic regime.

This test characterizes the relativistic regime:
  1. Build H_rel(alpha=CODATA) and H_rel(alpha=0) for LiH.
  2. Compute || H_rel(alpha) - H_rel(0) ||_op_norm via term-count and 1-norm.
  3. Show alpha^2 scaling of the spin-orbit perturbation by varying alpha.
  4. Construct J_z in qubit form and verify [H_rel, J_z] = 0.
  5. Construct the m_j -> -m_j parity operator P_spinor and check
     [H_rel, P_spinor] = 0.

The expected pattern:
  - H_rel = H_rel(alpha=0) + alpha^2 * V_SO + O(alpha^4)
  - [H_rel(alpha=0), S^2] = 0 (non-rel limit)
  - [H_rel(alpha), S^2] = O(alpha^2) - finite spin-orbit term
  - [H_rel, J^2] = [H_rel, J_z] = 0 to all orders
  - [H_rel, P_spinor (m_j -> -m_j)] = 0 - structural

If the last point holds, then a J^2-adapted QPT *would* still stack with a
spinor-extended Hopf-Z2 tapering. This is the scoping result.
"""

from __future__ import annotations

import json
import os
from typing import Dict, List, Tuple

import numpy as np

from openfermion import QubitOperator
from openfermion.utils import commutator

from geovac.composed_qubit_relativistic import build_composed_hamiltonian_relativistic
from geovac.molecular_spec import lih_spec_relativistic, beh_spec_relativistic


CODATA_ALPHA = 7.2973525693e-3


def _max_abs(op: QubitOperator) -> float:
    return max((abs(v) for v in op.terms.values()), default=0.0)


def _l1_norm(op: QubitOperator) -> float:
    """Sum of absolute coefficients (Pauli 1-norm)."""
    return sum(abs(v) for v in op.terms.values())


def commutator_max(a: QubitOperator, b: QubitOperator) -> float:
    return _max_abs(commutator(a, b))


def operator_diff_metrics(op_a: QubitOperator, op_b: QubitOperator) -> Dict:
    """Compute ||A - B|| in several senses."""
    diff = op_a - op_b
    return {
        'n_terms_diff': len(diff.terms),
        'max_abs_coef': _max_abs(diff),
        'l1_norm': _l1_norm(diff),
    }


def build_spinor_orbital_table(spec):
    """Build orbital table for the relativistic builder.

    Each entry: (sub_block_key, n, kappa, m_j_doubled)
    where m_j_doubled = 2 * m_j (so values are integer odd e.g. 1, -1, 3, -3).
    """
    from geovac.composed_qubit_relativistic import enumerate_dirac_labels
    orbital_table = []
    for blk_idx, blk in enumerate(spec.blocks):
        l_min = getattr(blk, 'l_min', 0)
        center_labels = enumerate_dirac_labels(blk.max_n, l_min=l_min)
        for lab in center_labels:
            sb_key = (blk_idx, blk.label, 'center')
            orbital_table.append((sb_key, lab.n_fock, lab.kappa, lab.two_m_j))
        if blk.has_h_partner:
            partner_max_n = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            partner_labels = enumerate_dirac_labels(partner_max_n, l_min=0)
            for lab in partner_labels:
                sb_key = (blk_idx, blk.label, 'partner')
                orbital_table.append((sb_key, lab.n_fock, lab.kappa, lab.two_m_j))
    return orbital_table


def build_Jz_operator(orbital_table) -> QubitOperator:
    """Construct J_z = sum_i m_j(i) * n_i as a JW qubit operator.

    Each spinor orbital is its own JW qubit. n_i = (1 - Z_i) / 2.
    """
    Jz = QubitOperator()
    for i, entry in enumerate(orbital_table):
        if len(entry) >= 4:
            _sb, _n, _kappa, mj2 = entry[0], entry[1], entry[2], entry[3]
        else:
            continue
        # m_j is stored as 2*m_j in the relativistic builder (integer)
        mj = mj2 / 2.0
        # n_i = 1/2 - Z_i/2
        Jz += QubitOperator((), mj * 0.5)
        Jz += QubitOperator(((i, 'Z'),), -mj * 0.5)
    return Jz


def build_total_N_operator(n_qubits: int) -> QubitOperator:
    """Total particle number sum_i n_i as a JW qubit operator."""
    N = QubitOperator()
    for i in range(n_qubits):
        N += QubitOperator((), 0.5)
        N += QubitOperator(((i, 'Z'),), -0.5)
    return N


def build_mj_parity_projection(orbital_table) -> QubitOperator:
    """Build the m_j -> -m_j permutation projector as a Pauli Z-string.

    The permutation P_spinor flips m_j -> -m_j at fixed (n, kappa). On the
    symmetric/antisymmetric eigenbasis of P_spinor, this becomes a Z-string
    counting electrons on antisymmetric spinor combinations.

    For this DIAGNOSTIC test we use a simpler proxy: the global Z-string
    over all spinor qubits, similar to the global Hopf-Z2.
    """
    Z_global = QubitOperator()
    pauli = tuple((i, 'Z') for i in range(len(orbital_table)))
    Z_global += QubitOperator(pauli, 1.0)
    return Z_global


def test_alpha_scaling(spec, alphas: List[float]) -> Dict:
    """Build H_rel at multiple alpha values, characterize alpha-scaling."""
    Hs = {}
    for a in alphas:
        result = build_composed_hamiltonian_relativistic(
            spec, alpha_num=a, verbose=False,
        )
        Hs[a] = result['qubit_op']

    # Reference: alpha = 0 (non-relativistic)
    H0 = Hs[0.0]
    out = {
        'alphas': alphas,
        'H0_n_terms': len(H0.terms),
        'H0_l1_norm': _l1_norm(H0),
        'differences': {},
    }
    for a in alphas:
        if a == 0:
            continue
        diff = Hs[a] - H0
        out['differences'][a] = {
            'n_terms_diff': len(diff.terms),
            'max_abs_coef': _max_abs(diff),
            'l1_norm': _l1_norm(diff),
            'scaling_vs_alpha_squared': _l1_norm(diff) / (a ** 2),
        }
    return out, Hs


def main():
    results = {}

    for spec_fn, label in [(lih_spec_relativistic, 'LiH')]:
        print(f"\n{'='*70}\n{label} (relativistic)\n{'='*70}")
        spec = spec_fn()
        orbital_table = build_spinor_orbital_table(spec)
        Q = len(orbital_table)
        print(f"  Q = {Q} qubits (spinor JW)")

        # --- Alpha scaling test ---
        alphas = [0.0, CODATA_ALPHA / 2, CODATA_ALPHA, CODATA_ALPHA * 2]
        alpha_results, Hs = test_alpha_scaling(spec, alphas)
        print(f"\n  H_rel Pauli term counts and 1-norms:")
        print(f"    alpha=0:        n_terms={alpha_results['H0_n_terms']}, "
              f"1-norm={alpha_results['H0_l1_norm']:.4e}")
        for a, info in alpha_results['differences'].items():
            print(f"    alpha={a:.6f}: ||H - H_0||_1 = {info['l1_norm']:.4e}, "
                  f"||/alpha^2 = {info['scaling_vs_alpha_squared']:.4e}, "
                  f"max_coef = {info['max_abs_coef']:.4e}")

        # Determine if alpha-scaling is consistent with alpha^2
        ratios = [
            info['scaling_vs_alpha_squared']
            for info in alpha_results['differences'].values()
        ]
        scaling_consistent = (
            max(ratios) / min(ratios) < 1.05 if ratios else False
        )
        print(f"\n  Alpha^2 scaling consistency: {scaling_consistent} "
              f"(ratio max/min = {max(ratios)/min(ratios):.4f})")

        # --- Conserved-quantity tests at CODATA alpha ---
        H = Hs[CODATA_ALPHA]
        Jz = build_Jz_operator(orbital_table)
        N = build_total_N_operator(Q)
        P_global = build_mj_parity_projection(orbital_table)

        print(f"\n  Conserved-quantity commutators at alpha=CODATA:")
        c_HJz = commutator_max(H, Jz)
        c_HN = commutator_max(H, N)
        c_HP = commutator_max(H, P_global)
        c_JzN = commutator_max(Jz, N)
        print(f"    [H_rel, J_z]                = {c_HJz:.4e}   {'OK' if c_HJz < 1e-10 else 'NONZERO'}")
        print(f"    [H_rel, N_total]            = {c_HN:.4e}   {'OK' if c_HN < 1e-10 else 'NONZERO'}")
        print(f"    [H_rel, P_global_spinor]    = {c_HP:.4e}   {'OK' if c_HP < 1e-10 else 'NONZERO'}")
        print(f"    [J_z, N_total]              = {c_JzN:.4e}   {'OK' if c_JzN < 1e-10 else 'NONZERO'}")

        results[label] = {
            'Q': Q,
            'alpha_scaling': alpha_results,
            'alpha_squared_scaling_consistent': scaling_consistent,
            'commutators_at_codata': {
                '[H_rel, J_z]': c_HJz,
                '[H_rel, N]': c_HN,
                '[H_rel, P_global_spinor]': c_HP,
                '[J_z, N]': c_JzN,
            },
        }

    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'qpt_relativistic_test.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
