"""
Track BD: H2O composed 1-norm computation and analysis.

Computes:
1. H2O composed qubit Hamiltonian at n_max=1 (Q=14) and n_max=2 (Q=70)
2. 1-norm, Pauli terms, qubits for each
3. PK decomposition: how much 1-norm from PK vs electronic terms
4. BeH2 verification at n_max=1 (Q=10)
5. Scaling analysis across systems
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

from openfermion import QubitOperator
from geovac.composed_qubit import build_composed_h2o, build_composed_beh2, build_composed_lih
from geovac.trotter_bounds import pauli_1norm


def compute_1norm_decomposition(qubit_op):
    """Compute total 1-norm and decompose by Pauli weight."""
    total = 0.0
    by_weight = {}  # weight -> sum |c_i|
    max_coeff = 0.0
    n_terms = 0

    for term, coeff in qubit_op.terms.items():
        c = abs(coeff)
        total += c
        n_terms += 1
        if c > max_coeff:
            max_coeff = c
        weight = len(term)  # number of non-identity Paulis
        by_weight[weight] = by_weight.get(weight, 0.0) + c

    return {
        'total_1norm': total,
        'max_coeff': max_coeff,
        'n_terms': n_terms,
        'by_weight': {str(k): v for k, v in sorted(by_weight.items())},
    }


def compute_pk_decomposition(system_builder, builder_kwargs_with_pk, builder_kwargs_no_pk):
    """Compare 1-norm with and without PK to isolate PK contribution."""
    print("\n--- With PK ---")
    result_pk = system_builder(**builder_kwargs_with_pk)
    norm_pk = pauli_1norm(result_pk['qubit_op'])

    print("\n--- Without PK ---")
    result_no_pk = system_builder(**builder_kwargs_no_pk)
    norm_no_pk = pauli_1norm(result_no_pk['qubit_op'])

    return {
        'with_pk': norm_pk,
        'without_pk': norm_no_pk,
        'pk_contribution': norm_pk - norm_no_pk,
        'pk_fraction': (norm_pk - norm_no_pk) / norm_pk if norm_pk > 0 else 0,
        'n_pauli_with_pk': result_pk['N_pauli'],
        'n_pauli_without_pk': result_no_pk['N_pauli'],
    }


def main():
    results = {}

    # ================================================================
    # 1. H2O at n_max=1 (Q=14)
    # ================================================================
    print("=" * 70)
    print("H2O n_max=1")
    print("=" * 70)

    t0 = time.perf_counter()
    h2o_n1 = build_composed_h2o(max_n_core=1, max_n_val=1, verbose=True)
    t1 = time.perf_counter()

    decomp_n1 = compute_1norm_decomposition(h2o_n1['qubit_op'])
    norm_n1 = decomp_n1['total_1norm']

    print(f"\nH2O n_max=1: Q={h2o_n1['Q']}, N_pauli={h2o_n1['N_pauli']}, "
          f"1-norm={norm_n1:.4f} Ha, time={t1-t0:.1f}s")

    results['h2o_nmax1'] = {
        'Q': h2o_n1['Q'],
        'N_pauli': h2o_n1['N_pauli'],
        'one_norm': norm_n1,
        'max_coeff': decomp_n1['max_coeff'],
        'by_weight': decomp_n1['by_weight'],
        'M': h2o_n1['M'],
        'wall_time_s': t1 - t0,
        'nuclear_repulsion': h2o_n1['nuclear_repulsion'],
    }

    # PK decomposition for n_max=1
    print("\n--- PK decomposition (n_max=1) ---")
    pk_decomp_n1 = compute_pk_decomposition(
        build_composed_h2o,
        {'max_n_core': 1, 'max_n_val': 1, 'verbose': False, 'include_pk': True},
        {'max_n_core': 1, 'max_n_val': 1, 'verbose': False, 'include_pk': False},
    )
    print(f"  With PK:    1-norm = {pk_decomp_n1['with_pk']:.4f} Ha")
    print(f"  Without PK: 1-norm = {pk_decomp_n1['without_pk']:.4f} Ha")
    print(f"  PK contribution:     {pk_decomp_n1['pk_contribution']:.4f} Ha "
          f"({pk_decomp_n1['pk_fraction']:.1%})")
    results['h2o_nmax1']['pk_decomposition'] = pk_decomp_n1

    # ================================================================
    # 2. H2O at n_max=2 (Q=70)
    # ================================================================
    print("\n" + "=" * 70)
    print("H2O n_max=2")
    print("=" * 70)

    t0 = time.perf_counter()
    h2o_n2 = build_composed_h2o(max_n_core=2, max_n_val=2, verbose=True)
    t1 = time.perf_counter()

    decomp_n2 = compute_1norm_decomposition(h2o_n2['qubit_op'])
    norm_n2 = decomp_n2['total_1norm']

    print(f"\nH2O n_max=2: Q={h2o_n2['Q']}, N_pauli={h2o_n2['N_pauli']}, "
          f"1-norm={norm_n2:.4f} Ha, time={t1-t0:.1f}s")

    results['h2o_nmax2'] = {
        'Q': h2o_n2['Q'],
        'N_pauli': h2o_n2['N_pauli'],
        'one_norm': norm_n2,
        'max_coeff': decomp_n2['max_coeff'],
        'by_weight': decomp_n2['by_weight'],
        'M': h2o_n2['M'],
        'wall_time_s': t1 - t0,
        'nuclear_repulsion': h2o_n2['nuclear_repulsion'],
    }

    # PK decomposition for n_max=2
    print("\n--- PK decomposition (n_max=2) ---")
    pk_decomp_n2 = compute_pk_decomposition(
        build_composed_h2o,
        {'max_n_core': 2, 'max_n_val': 2, 'verbose': False, 'include_pk': True},
        {'max_n_core': 2, 'max_n_val': 2, 'verbose': False, 'include_pk': False},
    )
    print(f"  With PK:    1-norm = {pk_decomp_n2['with_pk']:.4f} Ha")
    print(f"  Without PK: 1-norm = {pk_decomp_n2['without_pk']:.4f} Ha")
    print(f"  PK contribution:     {pk_decomp_n2['pk_contribution']:.4f} Ha "
          f"({pk_decomp_n2['pk_fraction']:.1%})")
    results['h2o_nmax2']['pk_decomposition'] = pk_decomp_n2

    # ================================================================
    # 3. H2O scaling exponent
    # ================================================================
    if norm_n1 > 0 and norm_n2 > 0:
        Q1, Q2 = h2o_n1['Q'], h2o_n2['Q']
        exponent = np.log(norm_n2 / norm_n1) / np.log(Q2 / Q1)
        print(f"\nH2O 1-norm scaling: Q^{exponent:.2f}")
        print(f"  (from Q={Q1}: {norm_n1:.2f} to Q={Q2}: {norm_n2:.2f})")
        results['h2o_scaling'] = {
            'exponent': exponent,
            'Q_values': [Q1, Q2],
            'norm_values': [norm_n1, norm_n2],
        }

    # ================================================================
    # 4. BeH2 verification at n_max=1 (Q=10), expected 198.20 Ha
    # ================================================================
    print("\n" + "=" * 70)
    print("BeH2 verification (n_max=1, Q=10)")
    print("=" * 70)

    beh2_n1 = build_composed_beh2(max_n_core=1, max_n_val=1, verbose=True)
    norm_beh2_n1 = pauli_1norm(beh2_n1['qubit_op'])
    print(f"\nBeH2 n_max=1: Q={beh2_n1['Q']}, N_pauli={beh2_n1['N_pauli']}, "
          f"1-norm={norm_beh2_n1:.4f} Ha")
    print(f"  Expected: 198.20 Ha -> {'MATCH' if abs(norm_beh2_n1 - 198.20) < 0.1 else 'MISMATCH'}")

    results['beh2_nmax1_verify'] = {
        'Q': beh2_n1['Q'],
        'N_pauli': beh2_n1['N_pauli'],
        'one_norm': norm_beh2_n1,
        'expected': 198.20,
        'match': abs(norm_beh2_n1 - 198.20) < 0.1,
    }

    # Also compute BeH2 at n_max=2 for comparison
    print("\n" + "=" * 70)
    print("BeH2 n_max=2 (Q=50)")
    print("=" * 70)

    beh2_n2 = build_composed_beh2(max_n_core=2, max_n_val=2, verbose=True)
    norm_beh2_n2 = pauli_1norm(beh2_n2['qubit_op'])
    print(f"\nBeH2 n_max=2: Q={beh2_n2['Q']}, N_pauli={beh2_n2['N_pauli']}, "
          f"1-norm={norm_beh2_n2:.4f} Ha")
    print(f"  Expected: 354.89 Ha -> {'MATCH' if abs(norm_beh2_n2 - 354.89) < 0.1 else 'MISMATCH'}")

    results['beh2_nmax2_verify'] = {
        'Q': beh2_n2['Q'],
        'N_pauli': beh2_n2['N_pauli'],
        'one_norm': norm_beh2_n2,
        'expected': 354.89,
        'match': abs(norm_beh2_n2 - 354.89) < 0.1,
    }

    # ================================================================
    # 5. LiH at n_max=1 and n_max=2 for comparison
    # ================================================================
    print("\n" + "=" * 70)
    print("LiH n_max=1 and n_max=2 (for comparison)")
    print("=" * 70)

    lih_n1 = build_composed_lih(max_n_core=1, max_n_val=1, verbose=False)
    norm_lih_n1 = pauli_1norm(lih_n1['qubit_op'])
    print(f"LiH n_max=1: Q={lih_n1['Q']}, N_pauli={lih_n1['N_pauli']}, "
          f"1-norm={norm_lih_n1:.4f} Ha")

    lih_n2 = build_composed_lih(max_n_core=2, max_n_val=2, verbose=False)
    norm_lih_n2 = pauli_1norm(lih_n2['qubit_op'])
    print(f"LiH n_max=2: Q={lih_n2['Q']}, N_pauli={lih_n2['N_pauli']}, "
          f"1-norm={norm_lih_n2:.4f} Ha")
    print(f"  Expected n_max=2: 37.33 Ha -> {'MATCH' if abs(norm_lih_n2 - 37.33) < 0.1 else 'MISMATCH'}")

    # LiH PK decomposition at n_max=2
    pk_decomp_lih = compute_pk_decomposition(
        build_composed_lih,
        {'max_n_core': 2, 'max_n_val': 2, 'verbose': False, 'include_pk': True},
        {'max_n_core': 2, 'max_n_val': 2, 'verbose': False, 'include_pk': False},
    )
    print(f"LiH PK decomposition (n_max=2):")
    print(f"  With PK:    {pk_decomp_lih['with_pk']:.4f} Ha")
    print(f"  Without PK: {pk_decomp_lih['without_pk']:.4f} Ha")
    print(f"  PK fraction: {pk_decomp_lih['pk_fraction']:.1%}")

    results['lih_nmax1'] = {
        'Q': lih_n1['Q'],
        'N_pauli': lih_n1['N_pauli'],
        'one_norm': norm_lih_n1,
    }
    results['lih_nmax2'] = {
        'Q': lih_n2['Q'],
        'N_pauli': lih_n2['N_pauli'],
        'one_norm': norm_lih_n2,
        'pk_decomposition': pk_decomp_lih,
    }

    # LiH scaling
    lih_exp = np.log(norm_lih_n2 / norm_lih_n1) / np.log(lih_n2['Q'] / lih_n1['Q'])
    print(f"LiH 1-norm scaling: Q^{lih_exp:.2f}")
    results['lih_scaling'] = {'exponent': lih_exp}

    # BeH2 scaling
    beh2_exp = np.log(norm_beh2_n2 / norm_beh2_n1) / np.log(beh2_n2['Q'] / beh2_n1['Q'])
    print(f"BeH2 1-norm scaling: Q^{beh2_exp:.2f}")
    results['beh2_scaling'] = {'exponent': beh2_exp}

    # ================================================================
    # 6. Cross-system comparison table
    # ================================================================
    print("\n" + "=" * 70)
    print("CROSS-SYSTEM 1-NORM COMPARISON")
    print("=" * 70)
    print(f"{'System':<10} {'n_max':<6} {'Q':<6} {'N_pauli':<10} {'1-norm (Ha)':<14} {'Scaling'}")
    print("-" * 60)
    print(f"{'LiH':<10} {'1':<6} {lih_n1['Q']:<6} {lih_n1['N_pauli']:<10} {norm_lih_n1:<14.2f}")
    print(f"{'LiH':<10} {'2':<6} {lih_n2['Q']:<6} {lih_n2['N_pauli']:<10} {norm_lih_n2:<14.2f} Q^{lih_exp:.2f}")
    print(f"{'BeH2':<10} {'1':<6} {beh2_n1['Q']:<6} {beh2_n1['N_pauli']:<10} {norm_beh2_n1:<14.2f}")
    print(f"{'BeH2':<10} {'2':<6} {beh2_n2['Q']:<6} {beh2_n2['N_pauli']:<10} {norm_beh2_n2:<14.2f} Q^{beh2_exp:.2f}")
    print(f"{'H2O':<10} {'1':<6} {h2o_n1['Q']:<6} {h2o_n1['N_pauli']:<10} {norm_n1:<14.2f}")
    print(f"{'H2O':<10} {'2':<6} {h2o_n2['Q']:<6} {h2o_n2['N_pauli']:<10} {norm_n2:<14.2f} Q^{results.get('h2o_scaling', {}).get('exponent', 0):.2f}")

    print("\n" + "=" * 70)
    print("PK DECOMPOSITION COMPARISON (n_max=2)")
    print("=" * 70)
    print(f"{'System':<10} {'With PK':<14} {'Without PK':<14} {'PK contrib':<14} {'PK %'}")
    print("-" * 60)
    print(f"{'LiH':<10} {pk_decomp_lih['with_pk']:<14.2f} {pk_decomp_lih['without_pk']:<14.2f} "
          f"{pk_decomp_lih['pk_contribution']:<14.2f} {pk_decomp_lih['pk_fraction']:.1%}")
    print(f"{'H2O':<10} {pk_decomp_n2['with_pk']:<14.2f} {pk_decomp_n2['without_pk']:<14.2f} "
          f"{pk_decomp_n2['pk_contribution']:<14.2f} {pk_decomp_n2['pk_fraction']:.1%}")

    # ================================================================
    # 7. Save results
    # ================================================================
    out_dir = Path(__file__).resolve().parent

    # JSON data
    json_path = out_dir / 'h2o_1norm_data.json'

    # Make all values JSON-serializable
    def make_serializable(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, dict):
            return {k: make_serializable(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [make_serializable(v) for v in obj]
        return obj

    with open(json_path, 'w') as f:
        json.dump(make_serializable(results), f, indent=2)
    print(f"\nData saved to {json_path}")

    return results


if __name__ == '__main__':
    main()
