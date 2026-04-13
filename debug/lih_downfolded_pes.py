"""
LiH PES comparison: PK vs downfolded core potential.

Scans R from 2.0 to 5.0 bohr, computing FCI energy at each point
with both the PK Gaussian barrier and the exact mean-field (2J-K)
downfolded core potential.

Author: GeoVac Development Team
Date: April 2026
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import json
import time
import numpy as np
from geovac.downfolding import compute_cross_rk
from geovac.composed_qubit import (
    _enumerate_states, _ck_coefficient, _radial_wf_grid,
    _compute_rk_integrals_block, _build_eri_block,
    _compute_pk_matrix_elements,
    build_fermion_op_from_integrals,
)
from openfermion import jordan_wigner, get_sparse_operator
from scipy.sparse.linalg import eigsh


def compute_downfolded_potential(
    states_val: list, Z_val: float, Z_core: float,
    n_grid: int = 8000,
) -> np.ndarray:
    """Compute exact mean-field (2J - K) core potential on valence orbitals.

    Assumes a single occupied 1s core orbital at Z_core.
    """
    Mv = len(states_val)
    V_df = np.zeros((Mv, Mv))

    for p, (np_, lp, mp) in enumerate(states_val):
        for q, (nq, lq, mq) in enumerate(states_val):
            if q < p:
                V_df[p, q] = V_df[q, p]
                continue

            # Coulomb: J(p,q) = (pq|cc) — pair (p,q) at Z_val, pair (c,c) at Z_core
            J_val = 0.0
            k_max_J = min(lp + lq, 0)  # core is 1s → l_c=0
            for k in range(k_max_J + 1):
                c_pq = _ck_coefficient(lp, mp, lq, mq, k)
                c_cc = _ck_coefficient(0, 0, 0, 0, k)
                if abs(c_pq) < 1e-15 or abs(c_cc) < 1e-15:
                    continue
                rk = compute_cross_rk(
                    np_, lp, Z_val, nq, lq, Z_val,
                    1, 0, Z_core, 1, 0, Z_core, k, n_grid=n_grid,
                )
                J_val += c_pq * c_cc * rk

            # Exchange: K(p,q) = (pc|cq) — mixed pairs
            K_val = 0.0
            k_max_K = min(lp, lq)
            for k in range(k_max_K + 1):
                c_pc = _ck_coefficient(lp, mp, 0, 0, k)
                c_cq = _ck_coefficient(0, 0, lq, mq, k)
                if abs(c_pc) < 1e-15 or abs(c_cq) < 1e-15:
                    continue
                rk = compute_cross_rk(
                    np_, lp, Z_val, 1, 0, Z_core,
                    1, 0, Z_core, nq, lq, Z_val, k, n_grid=n_grid,
                )
                K_val += c_pc * c_cq * rk

            V_df[p, q] = 2 * J_val - K_val
            V_df[q, p] = 2 * J_val - K_val

    return V_df


def build_lih_valence_hamiltonian(
    R: float,
    max_n: int = 2,
    core_method: str = 'pk',  # 'pk', 'downfolded', or 'none'
) -> dict:
    """Build LiH VALENCE-ONLY Hamiltonian (no core orbitals, 2 valence electrons).

    Valence space: screened Li (Z_eff=1) + H (Z=1) orbitals.
    Core effect enters via PK or downfolded potential on Li valence orbitals.
    FCI-tractable: Q=20 for max_n=2 (2 electrons in 10 spatial orbitals).
    """
    Z_A = 3.0   # Li
    Z_B = 1.0   # H
    Z_eff = 1.0  # screened Li valence
    E_core = -7.2799  # Li²⁺ core energy

    states_val_li = _enumerate_states(max_n)
    states_val_h = _enumerate_states(max_n)
    Mv_li = len(states_val_li)
    Mv_h = len(states_val_h)
    M = Mv_li + Mv_h
    Q = 2 * M

    off_h = Mv_li

    # --- h1 ---
    h1 = np.zeros((M, M))
    for i, (n, l, m) in enumerate(states_val_li):
        h1[i, i] = -Z_eff ** 2 / (2 * n ** 2)
    for i, (n, l, m) in enumerate(states_val_h):
        h1[off_h + i, off_h + i] = -Z_B ** 2 / (2 * n ** 2)

    # --- Core potential on valence-Li ---
    if core_method == 'pk':
        pk_mat = _compute_pk_matrix_elements(Z_eff, states_val_li, 6.93, 7.00)
        for i in range(Mv_li):
            for j in range(Mv_li):
                h1[i, j] += pk_mat[i, j]
        h1_core_pot = pk_mat
    elif core_method == 'downfolded':
        V_df = compute_downfolded_potential(states_val_li, Z_eff, Z_A)
        for i in range(Mv_li):
            for j in range(Mv_li):
                h1[i, j] += V_df[i, j]
        h1_core_pot = V_df
    else:
        h1_core_pot = np.zeros((Mv_li, Mv_li))

    # --- ERI (block diagonal, valence only) ---
    eri = np.zeros((M, M, M, M))

    rk_vli = _compute_rk_integrals_block(Z_eff, states_val_li)
    eri_vli = _build_eri_block(Z_eff, states_val_li, rk_vli)
    for (a, b, c, d), val in eri_vli.items():
        eri[a, c, b, d] += val

    rk_vh = _compute_rk_integrals_block(Z_B, states_val_h)
    eri_vh = _build_eri_block(Z_B, states_val_h, rk_vh)
    for (a, b, c, d), val in eri_vh.items():
        eri[off_h + a, off_h + c, off_h + b, off_h + d] += val

    eri = 0.5 * (eri + eri.transpose(2, 3, 0, 1))

    # --- Nuclear repulsion ---
    V_NN = Z_A * Z_B / R
    from geovac.composed_qubit import _v_cross_nuc_1s
    V_cross = _v_cross_nuc_1s(Z_A, 2, Z_B, R)
    nuc_rep = V_NN + V_cross + E_core

    # --- JW transform ---
    fermion_op = build_fermion_op_from_integrals(h1, eri, nuc_rep)
    qubit_op = jordan_wigner(fermion_op)

    return {
        'qubit_op': qubit_op,
        'Q': Q,
        'M': M,
        'R': R,
        'nuc_rep': nuc_rep,
        'core_method': core_method,
        'h1_core_pot_diag': np.diag(h1_core_pot).tolist(),
        'n_electrons': 2,
    }


def fci_energy(qubit_op: 'QubitOperator', n_qubits: int,
               n_electrons: int = 4) -> float:
    """Compute FCI ground state energy from qubit operator."""
    sparse_H = get_sparse_operator(qubit_op, n_qubits=n_qubits)
    evals, _ = eigsh(sparse_H, k=1, which='SA')
    return float(evals[0])


def main():
    R_values = np.arange(2.0, 5.5, 0.25)
    R_exact = 3.015  # experimental R_eq

    print("LiH PES Comparison: PK vs Downfolded vs No-Core-Potential")
    print("=" * 70)
    print(f"{'R (bohr)':>10s}  {'E_PK (Ha)':>12s}  {'E_DF (Ha)':>12s}  {'E_none (Ha)':>12s}")
    print("-" * 70)

    results = {'R': [], 'E_pk': [], 'E_df': [], 'E_none': []}

    for R in R_values:
        energies = {}
        for method in ['pk', 'downfolded', 'none']:
            t0 = time.perf_counter()
            ham = build_lih_valence_hamiltonian(R, max_n=2, core_method=method)
            E = fci_energy(ham['qubit_op'], ham['Q'])
            wall = time.perf_counter() - t0
            energies[method] = E

        results['R'].append(float(R))
        results['E_pk'].append(energies['pk'])
        results['E_df'].append(energies['downfolded'])
        results['E_none'].append(energies['none'])

        print(f"{R:10.3f}  {energies['pk']:12.6f}  {energies['downfolded']:12.6f}  {energies['none']:12.6f}")

    # Find R_eq for each method
    for method, key in [('PK', 'E_pk'), ('Downfolded', 'E_df'), ('None', 'E_none')]:
        Es = results[key]
        Rs = results['R']
        i_min = np.argmin(Es)
        if 0 < i_min < len(Rs) - 1:
            # Parabolic interpolation for R_eq
            r1, r2, r3 = Rs[i_min - 1], Rs[i_min], Rs[i_min + 1]
            e1, e2, e3 = Es[i_min - 1], Es[i_min], Es[i_min + 1]
            denom = 2 * ((r2 - r1) * (e2 - e3) - (r2 - r3) * (e2 - e1))
            if abs(denom) > 1e-12:
                R_eq = r2 - ((r2 - r1) ** 2 * (e2 - e3) - (r2 - r3) ** 2 * (e2 - e1)) / denom
            else:
                R_eq = r2
            E_min = min(Es)
            err = abs(R_eq - R_exact) / R_exact * 100
            print(f"\n{method:12s}: R_eq = {R_eq:.3f} bohr ({err:.1f}% error), E_min = {E_min:.6f} Ha")
        else:
            print(f"\n{method:12s}: minimum at edge of scan (R={Rs[i_min]:.3f})")

    # Exact reference
    print(f"\n{'Experiment':12s}: R_eq = {R_exact:.3f} bohr, E_exact = -8.0706 Ha")

    # Save data
    out_path = Path('debug/data/lih_downfolded_pes.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nData saved to {out_path}")


if __name__ == '__main__':
    main()
