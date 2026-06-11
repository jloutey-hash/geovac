"""
Track DF Sprint 3: Be atom qubit Hamiltonian comparison.

Builds two Be Hamiltonians:
1. "Nested" (full-atom): single block, all orbitals at Z=4, all ERIs, no PK
2. "Composed": core (1s² at Z=4) + valence (2s² at Z_eff=2 with PK)

Both are JW-encoded and compared on: Q, Pauli terms, ERI density, 1-norm.
Also compares at max_n=2 and max_n=3 to measure scaling.
"""

import sys
import json
import time
import numpy as np

sys.path.insert(0, '.')

from geovac.molecular_spec import MolecularSpec, OrbitalBlock
from geovac.composed_qubit import build_composed_hamiltonian


# ============================================================================
# Spec factories
# ============================================================================

def be_nested_spec(max_n: int = 2) -> MolecularSpec:
    """Full-atom Be: single block, all orbitals at Z=4, no PK.

    This is the "nested" approach: all 4 electrons in a single Hilbert space,
    all ERIs computed with Gaunt selection rules, no block decomposition.
    """
    return MolecularSpec(
        name=f'Be_nested_n{max_n}',
        blocks=[
            OrbitalBlock(
                label='Be_full',
                block_type='core',  # arbitrary label
                Z_center=4.0,
                n_electrons=4,
                max_n=max_n,
                pk_A=0.0,  # no PK
                pk_B=0.0,
            ),
        ],
        nuclear_repulsion_constant=0.0,  # atom, no V_NN
    )


def be_composed_spec(max_n: int = 2) -> MolecularSpec:
    """Composed Be: core (1s²) + valence (2s²) with PK.

    Mimics the composed architecture: core block at Z=4 with n_max=1,
    valence block at Z_eff=2 (Z - n_core) with PK barrier from core.
    """
    Z = 4
    n_core = 2
    Z_eff_val = Z - n_core  # 2.0

    # PK parameters from Paper 17 Table 1 (Be³⁺ inv2)
    A_pk = 13.01
    B_pk = 12.53

    return MolecularSpec(
        name=f'Be_composed_n{max_n}',
        blocks=[
            OrbitalBlock(
                label='Be_core',
                block_type='core',
                Z_center=float(Z),
                n_electrons=2,
                max_n=1,  # only 1s
                pk_A=0.0,
                pk_B=0.0,
            ),
            OrbitalBlock(
                label='Be_valence',
                block_type='bond_pair',
                Z_center=float(Z_eff_val),
                n_electrons=2,
                max_n=max_n,
                pk_A=A_pk,
                pk_B=B_pk,
            ),
        ],
        nuclear_repulsion_constant=0.0,
    )


# ============================================================================
# Metrics computation
# ============================================================================

def compute_metrics(result: dict) -> dict:
    """Extract key metrics from build_composed_hamiltonian result."""
    Q = result['Q']
    N_pauli = result['N_pauli']
    eri = result['eri']
    h1 = result['h1']
    qubit_op = result['qubit_op']

    # ERI density
    M = result['M']
    n_eri_nz = int(np.count_nonzero(np.abs(eri) > 1e-15))
    eri_density = n_eri_nz / max(1, M**4)

    # h1 density
    n_h1_nz = int(np.count_nonzero(np.abs(h1) > 1e-15))
    h1_density = n_h1_nz / max(1, M**2)

    # 1-norm: sum of absolute values of Pauli coefficients
    one_norm = sum(abs(coeff) for coeff in qubit_op.terms.values())

    # Non-identity 1-norm (exclude identity term)
    from openfermion import QubitOperator
    identity = QubitOperator(())
    id_coeff = qubit_op.terms.get((), 0.0)
    one_norm_ni = one_norm - abs(id_coeff)

    # Pauli/Q ratio
    pauli_per_q = N_pauli / Q if Q > 0 else 0

    return {
        'M': M,
        'Q': Q,
        'N_pauli': N_pauli,
        'Pauli_per_Q': round(pauli_per_q, 2),
        'n_eri_nonzero': n_eri_nz,
        'ERI_density_pct': round(eri_density * 100, 2),
        'n_h1_nonzero': n_h1_nz,
        'h1_density_pct': round(h1_density * 100, 2),
        'one_norm_Ha': round(one_norm, 2),
        'one_norm_ni_Ha': round(one_norm_ni, 2),
        'wall_time_s': round(result['wall_time_s'], 2),
    }


# ============================================================================
# FCI energy (optional, for validation)
# ============================================================================

def fci_energy_from_qubit_op(qubit_op, n_qubits: int, n_electrons: int) -> float:
    """Compute FCI ground state energy via exact diagonalization of qubit Hamiltonian.

    Only feasible for small Q (< ~24).
    """
    from openfermion import get_sparse_operator
    import scipy.sparse.linalg as spla

    sparse_H = get_sparse_operator(qubit_op, n_qubits=n_qubits)
    # Find ground state
    if sparse_H.shape[0] <= 64:
        # Small enough for dense diag
        H_dense = sparse_H.toarray()
        evals = np.linalg.eigvalsh(H_dense)
        return float(evals[0])
    else:
        evals, _ = spla.eigsh(sparse_H, k=1, which='SA')
        return float(evals[0])


# ============================================================================
# Main
# ============================================================================

def main():
    all_results = {}

    for max_n in [2, 3]:
        print(f"\n{'='*70}")
        print(f"  Be COMPARISON: max_n = {max_n}")
        print(f"{'='*70}")

        # --- Nested (full-atom) ---
        print(f"\n--- Nested (full-atom, no PK) ---")
        spec_nested = be_nested_spec(max_n=max_n)
        t0 = time.perf_counter()
        res_nested = build_composed_hamiltonian(spec_nested, pk_in_hamiltonian=False, verbose=True)
        metrics_nested = compute_metrics(res_nested)
        print(f"  Q={metrics_nested['Q']}, Pauli={metrics_nested['N_pauli']}, "
              f"ERI density={metrics_nested['ERI_density_pct']}%, "
              f"1-norm={metrics_nested['one_norm_Ha']} Ha")

        # FCI energy for small systems
        if metrics_nested['Q'] <= 24:
            try:
                E_fci = fci_energy_from_qubit_op(
                    res_nested['qubit_op'], metrics_nested['Q'], 4)
                metrics_nested['FCI_energy_Ha'] = round(E_fci, 6)
                # Be exact: -14.6674 Ha
                err = abs(E_fci - (-14.6674)) / 14.6674 * 100
                metrics_nested['FCI_error_pct'] = round(err, 3)
                print(f"  FCI energy: {E_fci:.6f} Ha (error: {err:.3f}%)")
            except Exception as e:
                print(f"  FCI failed: {e}")

        # --- Composed (core + valence + PK) ---
        print(f"\n--- Composed (core + valence + PK) ---")
        spec_composed = be_composed_spec(max_n=max_n)
        res_composed = build_composed_hamiltonian(spec_composed, pk_in_hamiltonian=True, verbose=True)
        metrics_composed = compute_metrics(res_composed)
        print(f"  Q={metrics_composed['Q']}, Pauli={metrics_composed['N_pauli']}, "
              f"ERI density={metrics_composed['ERI_density_pct']}%, "
              f"1-norm={metrics_composed['one_norm_Ha']} Ha")

        # Composed without PK for comparison
        print(f"\n--- Composed (core + valence, NO PK) ---")
        res_composed_nopk = build_composed_hamiltonian(spec_composed, pk_in_hamiltonian=False, verbose=True)
        metrics_composed_nopk = compute_metrics(res_composed_nopk)
        print(f"  Q={metrics_composed_nopk['Q']}, Pauli={metrics_composed_nopk['N_pauli']}, "
              f"ERI density={metrics_composed_nopk['ERI_density_pct']}%, "
              f"1-norm={metrics_composed_nopk['one_norm_Ha']} Ha")

        # FCI for composed (small enough?)
        if metrics_composed['Q'] <= 24:
            try:
                E_fci_comp = fci_energy_from_qubit_op(
                    res_composed['qubit_op'], metrics_composed['Q'], 4)
                metrics_composed['FCI_energy_Ha'] = round(E_fci_comp, 6)
                err = abs(E_fci_comp - (-14.6674)) / 14.6674 * 100
                metrics_composed['FCI_error_pct'] = round(err, 3)
                print(f"  Composed FCI energy: {E_fci_comp:.6f} Ha (error: {err:.3f}%)")
            except Exception as e:
                print(f"  FCI failed: {e}")

        # --- Summary table ---
        print(f"\n--- Summary (max_n={max_n}) ---")
        print(f"{'Metric':<30} {'Nested':>12} {'Composed':>12} {'Composed(noPK)':>15}")
        print('-' * 72)
        for key in ['Q', 'N_pauli', 'Pauli_per_Q', 'ERI_density_pct', 'h1_density_pct',
                     'one_norm_Ha', 'one_norm_ni_Ha']:
            v_n = metrics_nested.get(key, 'N/A')
            v_c = metrics_composed.get(key, 'N/A')
            v_cn = metrics_composed_nopk.get(key, 'N/A')
            print(f"{key:<30} {str(v_n):>12} {str(v_c):>12} {str(v_cn):>15}")

        all_results[f'max_n_{max_n}'] = {
            'nested': metrics_nested,
            'composed': metrics_composed,
            'composed_nopk': metrics_composed_nopk,
        }

    # --- Scaling analysis ---
    print(f"\n{'='*70}")
    print(f"  SCALING ANALYSIS")
    print(f"{'='*70}")
    if 'max_n_2' in all_results and 'max_n_3' in all_results:
        for arch in ['nested', 'composed']:
            Q2 = all_results['max_n_2'][arch]['Q']
            Q3 = all_results['max_n_3'][arch]['Q']
            P2 = all_results['max_n_2'][arch]['N_pauli']
            P3 = all_results['max_n_3'][arch]['N_pauli']
            if Q2 > 0 and Q3 > 0 and P2 > 0 and P3 > 0:
                import math
                alpha = math.log(P3 / P2) / math.log(Q3 / Q2)
                print(f"{arch}: Q={Q2}->{Q3}, Pauli={P2}->{P3}, scaling exponent alpha={alpha:.2f}")

    # Save
    output_path = 'debug/data/track_df_be_comparison.json'
    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to {output_path}")

    # Go/no-go
    print(f"\n{'='*70}")
    print(f"  GO/NO-GO (Sprint 3)")
    print(f"{'='*70}")
    n2_nested = all_results.get('max_n_2', {}).get('nested', {})
    n2_composed = all_results.get('max_n_2', {}).get('composed', {})
    pq_nested = n2_nested.get('Pauli_per_Q', 999)
    pq_composed = n2_composed.get('Pauli_per_Q', 0)
    ratio = pq_nested / pq_composed if pq_composed > 0 else 999
    print(f"Nested Pauli/Q = {pq_nested}, Composed Pauli/Q = {pq_composed}")
    print(f"Ratio = {ratio:.2f} (gate: < 2.0)")
    if ratio < 2.0:
        print("PASS: Nested Pauli/Q within 2x of composed")
    else:
        print("FAIL: Nested Pauli/Q exceeds 2x of composed")


if __name__ == '__main__':
    main()
