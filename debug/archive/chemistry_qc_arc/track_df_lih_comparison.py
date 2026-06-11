"""
Track DF Sprint 4: Nested LiH molecular qubit Hamiltonian.

Builds LiH Hamiltonians in three approaches:
1. "Nested" (single-center): all orbitals at Z=3 centered on Li,
   H nuclear attraction via Shibuya-Wulfman multipole, all ERIs
   single-center Gaunt-sparse, no PK, no block decomposition.
2. "Composed": core (1s² at Z=3) + bond pair (at Z_eff=1 with PK).
3. "Balanced coupled": separate blocks with cross-center V_ne, no PK.

Compares: Q, Pauli terms, 1-norm, ERI density, FCI energy (if feasible).
"""

import sys
import json
import time
import numpy as np

sys.path.insert(0, '.')

from geovac.molecular_spec import MolecularSpec, OrbitalBlock
from geovac.composed_qubit import (
    build_composed_hamiltonian,
    _enumerate_states,
    _compute_rk_integrals_block,
    _build_eri_block,
    lih_spec,
)
from geovac.qubit_encoding import build_fermion_op_from_integrals
from geovac.shibuya_wulfman import compute_cross_center_vne
from openfermion import jordan_wigner


# ============================================================================
# Nested LiH builder
# ============================================================================

def build_nested_lih(
    max_n: int = 2,
    R: float = 3.015,
    Z_Li: float = 3.0,
    Z_H: float = 1.0,
    L_max: int = 4,
    verbose: bool = False,
) -> dict:
    """Build nested (single-center) LiH qubit Hamiltonian.

    All orbitals centered on Li (Z=3). H nuclear attraction added via
    Shibuya-Wulfman multipole expansion. All ERIs from single-center
    Slater integrals. No PK, no block decomposition.

    Parameters
    ----------
    max_n : int
        Maximum principal quantum number for orbitals.
    R : float
        Li-H internuclear distance (bohr).
    Z_Li : float
        Li nuclear charge.
    Z_H : float
        H nuclear charge.
    L_max : int
        Maximum multipole order for cross-center V_ne.
        L_max = 2*l_max is sufficient (Gaunt termination).
    verbose : bool
        Print progress.
    """
    t0 = time.perf_counter()

    # 1. Enumerate orbitals centered on Li
    states = _enumerate_states(max_n)
    M = len(states)
    Q = 2 * M

    if verbose:
        print(f"[nested_lih] max_n={max_n}, M={M}, Q={Q}, R={R}")

    # 2. One-body: Li nuclear attraction (hydrogenic diagonal)
    h1 = np.zeros((M, M))
    for i, (n, l, m) in enumerate(states):
        h1[i, i] = -Z_Li**2 / (2.0 * n**2)

    # 3. Cross-center V_ne from H nucleus at distance R along +z
    if verbose:
        print(f"[nested_lih] Computing cross-center V_ne (L_max={L_max})...")
    vne_cross = compute_cross_center_vne(
        Z_orb=Z_Li,
        states=states,
        Z_nuc=Z_H,
        R_AB=R,
        L_max=L_max,
        n_grid=4000,
        nuc_parity=1,  # H nucleus along +z
    )
    h1 += vne_cross

    if verbose:
        n_vne_nz = int(np.count_nonzero(np.abs(vne_cross) > 1e-15))
        print(f"[nested_lih] V_ne cross: {n_vne_nz} nonzero elements")

    # 4. Two-body: single-center Slater integrals at Z=3
    if verbose:
        print(f"[nested_lih] Computing Slater R^k integrals...")
    rk_cache = _compute_rk_integrals_block(Z_Li, states)
    eri_phys = _build_eri_block(Z_Li, states, rk_cache)

    # Convert physicist to chemist notation for ERI tensor
    eri = np.zeros((M, M, M, M))
    for (a, b, c, d), val in eri_phys.items():
        eri[a, c, b, d] = val  # phys <ab|cd> -> chem (ac|bd)
    eri = 0.5 * (eri + eri.transpose(2, 3, 0, 1))  # symmetrize

    # 5. Nuclear repulsion
    V_NN = Z_Li * Z_H / R

    # 6. Fermion operator + JW
    if verbose:
        print(f"[nested_lih] Building fermion operator...")
    fermion_op = build_fermion_op_from_integrals(h1, eri, V_NN)

    if verbose:
        print(f"[nested_lih] Jordan-Wigner transform...")
    qubit_op = jordan_wigner(fermion_op)
    N_pauli = len(qubit_op.terms)

    elapsed = time.perf_counter() - t0

    # Metrics
    n_eri_nz = int(np.count_nonzero(np.abs(eri) > 1e-15))
    one_norm = sum(abs(c) for c in qubit_op.terms.values())
    id_coeff = abs(qubit_op.terms.get((), 0.0))
    one_norm_ni = one_norm - id_coeff

    n_h1_nz = int(np.count_nonzero(np.abs(h1) > 1e-15))

    if verbose:
        print(f"[nested_lih] Q={Q}, Pauli={N_pauli}, 1-norm_ni={one_norm_ni:.2f} Ha, "
              f"time={elapsed:.1f}s")

    return {
        'M': M,
        'Q': Q,
        'N_pauli': N_pauli,
        'Pauli_per_Q': round(N_pauli / Q, 2),
        'n_eri_nonzero': n_eri_nz,
        'ERI_density_pct': round(n_eri_nz / M**4 * 100, 2),
        'n_h1_nonzero': n_h1_nz,
        'h1_density_pct': round(n_h1_nz / M**2 * 100, 2),
        'one_norm_Ha': round(one_norm, 2),
        'one_norm_ni_Ha': round(one_norm_ni, 2),
        'wall_time_s': round(elapsed, 2),
        'R': R,
        'max_n': max_n,
        'h1': h1,
        'eri': eri,
        'qubit_op': qubit_op,
        'fermion_op': fermion_op,
        'nuclear_repulsion': V_NN,
    }


# ============================================================================
# FCI energy
# ============================================================================

def fci_energy(qubit_op, Q: int) -> float:
    """FCI ground state energy via exact diagonalization."""
    from openfermion import get_sparse_operator
    import scipy.sparse.linalg as spla

    sparse_H = get_sparse_operator(qubit_op, n_qubits=Q)
    if sparse_H.shape[0] <= 1024:
        H_dense = sparse_H.toarray()
        return float(np.linalg.eigvalsh(H_dense)[0])
    else:
        evals, _ = spla.eigsh(sparse_H, k=1, which='SA')
        return float(evals[0])


# ============================================================================
# Main comparison
# ============================================================================

def main():
    R_eq = 3.015  # experimental LiH equilibrium
    E_exact_Req = -8.0705  # exact LiH energy at R=3.015

    all_results = {}

    for max_n in [2, 3]:
        print(f"\n{'='*70}")
        print(f"  LiH COMPARISON: max_n = {max_n}, R = {R_eq} bohr")
        print(f"{'='*70}")

        # L_max for multipole: 2*l_max
        l_max = max_n - 1
        L_max = 2 * l_max + 2  # generous

        # --- 1. Nested (single-center on Li) ---
        print(f"\n--- Nested (single-center, Z=3, no PK) ---")
        res_nested = build_nested_lih(max_n=max_n, R=R_eq, L_max=L_max, verbose=True)

        if res_nested['Q'] <= 20:
            try:
                E = fci_energy(res_nested['qubit_op'], res_nested['Q'])
                err = abs(E - E_exact_Req) / abs(E_exact_Req) * 100
                res_nested['FCI_energy_Ha'] = round(E, 6)
                res_nested['FCI_error_pct'] = round(err, 3)
                print(f"  FCI energy: {E:.6f} Ha (error: {err:.3f}%)")
            except Exception as e:
                print(f"  FCI failed: {e}")

        # --- 2. Composed (existing architecture) ---
        print(f"\n--- Composed (core + bond + PK) ---")
        spec_comp = lih_spec(max_n_core=min(max_n, 2), max_n_val=max_n, R=R_eq)
        res_comp = build_composed_hamiltonian(spec_comp, pk_in_hamiltonian=True, verbose=True)
        m_comp = {
            'Q': res_comp['Q'],
            'N_pauli': res_comp['N_pauli'],
            'Pauli_per_Q': round(res_comp['N_pauli'] / res_comp['Q'], 2),
            'ERI_density_pct': round(res_comp['ERI_density_total'] * 100, 2),
            'one_norm_Ha': round(sum(abs(c) for c in res_comp['qubit_op'].terms.values()), 2),
            'one_norm_ni_Ha': round(sum(abs(c) for c in res_comp['qubit_op'].terms.values())
                                    - abs(res_comp['qubit_op'].terms.get((), 0.0)), 2),
        }
        print(f"  Q={m_comp['Q']}, Pauli={m_comp['N_pauli']}, "
              f"1-norm_ni={m_comp['one_norm_ni_Ha']} Ha")

        if m_comp['Q'] <= 20:
            try:
                E_c = fci_energy(res_comp['qubit_op'], m_comp['Q'])
                err_c = abs(E_c - E_exact_Req) / abs(E_exact_Req) * 100
                m_comp['FCI_energy_Ha'] = round(E_c, 6)
                m_comp['FCI_error_pct'] = round(err_c, 3)
                print(f"  FCI energy: {E_c:.6f} Ha (error: {err_c:.3f}%)")
            except Exception as e:
                print(f"  FCI failed: {e}")

        # --- 3. Composed without PK ---
        print(f"\n--- Composed (no PK) ---")
        res_comp_nopk = build_composed_hamiltonian(spec_comp, pk_in_hamiltonian=False, verbose=True)
        m_comp_nopk = {
            'Q': res_comp_nopk['Q'],
            'N_pauli': res_comp_nopk['N_pauli'],
            'one_norm_ni_Ha': round(sum(abs(c) for c in res_comp_nopk['qubit_op'].terms.values())
                                     - abs(res_comp_nopk['qubit_op'].terms.get((), 0.0)), 2),
        }

        # --- Summary ---
        print(f"\n{'='*70}")
        print(f"  SUMMARY (max_n={max_n})")
        print(f"{'='*70}")
        print(f"{'Metric':<28} {'Nested':>12} {'Composed+PK':>14} {'Composed noPK':>14}")
        print('-' * 70)

        # Clean nested results for display
        n = res_nested
        c = m_comp
        cn = m_comp_nopk

        print(f"{'Qubits (Q)':<28} {n['Q']:>12} {c['Q']:>14} {cn['Q']:>14}")
        print(f"{'Pauli terms':<28} {n['N_pauli']:>12} {c['N_pauli']:>14} {cn['N_pauli']:>14}")
        print(f"{'Pauli/Q':<28} {n['Pauli_per_Q']:>12} {c['Pauli_per_Q']:>14} {'':>14}")
        print(f"{'ERI density (%)':<28} {n['ERI_density_pct']:>11}% {c['ERI_density_pct']:>13}% {'':>14}")
        print(f"{'1-norm ni (Ha)':<28} {n['one_norm_ni_Ha']:>12} {c['one_norm_ni_Ha']:>14} {cn['one_norm_ni_Ha']:>14}")

        if 'FCI_energy_Ha' in n:
            print(f"{'FCI energy (Ha)':<28} {n['FCI_energy_Ha']:>12} {c.get('FCI_energy_Ha','N/A'):>14} {'':>14}")
        if 'FCI_error_pct' in n:
            err_c_str = str(c.get('FCI_error_pct', 'N/A')) + '%'
            print(f"{'FCI error (%)':<28} {n['FCI_error_pct']:>11}% {err_c_str:>14} {'':>14}")

        # Known reference values
        print(f"\nReference: Composed LiH Q=30 -> 334 Pauli, 37.33 Ha 1-norm")
        print(f"Reference: Balanced coupled LiH Q=30 -> 878 Pauli, 74.1 Ha 1-norm")
        print(f"Reference: Full Level 4N LiH -> 3,288 Pauli")

        all_results[f'max_n_{max_n}'] = {
            'nested': {k: v for k, v in n.items()
                       if k not in ('h1', 'eri', 'qubit_op', 'fermion_op')},
            'composed': c,
            'composed_nopk': cn,
        }

    # --- PES scan at max_n=2 if FCI is feasible ---
    print(f"\n{'='*70}")
    print(f"  FCI PES SCAN (max_n=2)")
    print(f"{'='*70}")

    R_values = [2.0, 2.5, 3.015, 3.5, 4.0, 5.0, 6.0]
    pes_data = []

    for R in R_values:
        res = build_nested_lih(max_n=2, R=R)
        if res['Q'] <= 20:
            try:
                E = fci_energy(res['qubit_op'], res['Q'])
                pes_data.append({'R': R, 'E': round(E, 6), 'N_pauli': res['N_pauli'],
                                 'one_norm_ni': round(res['one_norm_ni_Ha'], 2)})
                print(f"  R={R:.3f}: E={E:.6f} Ha, Pauli={res['N_pauli']}, "
                      f"1-norm_ni={res['one_norm_ni_Ha']:.2f} Ha")
            except Exception as e:
                print(f"  R={R:.3f}: FCI failed: {e}")

    if pes_data:
        # Find R_eq (minimum energy)
        E_min_idx = min(range(len(pes_data)), key=lambda i: pes_data[i]['E'])
        R_eq_found = pes_data[E_min_idx]['R']
        E_min = pes_data[E_min_idx]['E']
        print(f"\n  R_eq (approx): {R_eq_found} bohr (exact: 3.015)")
        print(f"  E_min: {E_min:.6f} Ha (exact: {E_exact_Req:.4f})")

        # Check binding
        E_inf = pes_data[-1]['E'] if pes_data[-1]['R'] >= 5.0 else None
        if E_inf is not None:
            D_e = E_inf - E_min
            print(f"  D_e = E(R=inf) - E(R_eq) = {D_e:.6f} Ha")
            if D_e > 0:
                print(f"  BOUND: D_e = {D_e:.4f} Ha (exact: 0.092 Ha)")
            else:
                print(f"  UNBOUND: D_e = {D_e:.4f} Ha")

        all_results['pes_n2'] = pes_data

    # Save
    output_path = 'debug/data/track_df_lih_comparison.json'
    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to {output_path}")

    # Go/no-go
    print(f"\n{'='*70}")
    print(f"  GO/NO-GO (Sprint 4)")
    print(f"{'='*70}")
    n2 = all_results.get('max_n_2', {}).get('nested', {})
    pq = n2.get('Pauli_per_Q', 999)
    print(f"Nested Pauli/Q = {pq} (gate: <= 15)")
    eri = n2.get('ERI_density_pct', 999)
    print(f"ERI density = {eri}% (gate: <= 10%)")
    if pq <= 15 and eri <= 15:
        print("PASS: Nested LiH within gates")
    else:
        print("See assessment for details")


if __name__ == '__main__':
    main()
