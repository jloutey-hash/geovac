"""
Track EP-2b: Two-fermion entanglement on the Bargmann-Segal HO lattice.

Paper 27 §VII.A Prediction 1 — empirical test.

For the HO (Minnesota NN singlet channel) at N_max in {2, 3}:
  (a) Build H_kin (HO diagonal only) and H_full = H_kin + V_NN,
  (b) Diagonalize the spatial-singlet FCI (M_L = 0 block),
  (c) Report ground-state 1-RDM von Neumann entropy,
  (d) Report Frobenius-mass diagonality of V_NN in the H_kin eigenbasis
      and the hottest pair-state diagonal / off-diagonal edge.
  (e) Compare to He Coulomb EP-1 values.

Expected (from Paper 27 §VII.A):
  - S_HO << S_He at matched basis size (Minnesota's Gaussian kernel is
    short-range, no 1/r_12 singularity to drive hot-node concentration).
  - No single hot pair-state dominates (distributed Frobenius mass).

Refutation would be: S_HO ~ S_He or hot-node concentration on the HO
ground state.

Data saved to: debug/data/ep2b_ho_two_fermion.json
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from debug.entanglement_geometry import (
    build_1rdm_from_singlet_ci,
    compute_entanglement_measures,
)
from geovac.nuclear.ho_two_fermion import build_decomposed_ho_hamiltonians


def _frobenius_diagonality(H_vee: np.ndarray, H_kin: np.ndarray):
    """Return (unsquared Frobenius ratio, relative commutator norm)."""
    eH, U = np.linalg.eigh(H_kin)
    V_in = U.T @ H_vee @ U
    V_frob = float(np.linalg.norm(H_vee))
    diag_norm = float(np.sqrt(np.sum(np.diag(V_in) ** 2)))
    C = H_kin @ H_vee - H_vee @ H_kin
    kin_frob = float(np.linalg.norm(H_kin))
    return {
        'V_diag_fraction_unsquared': diag_norm / V_frob if V_frob > 0 else 0.0,
        'V_frobenius': V_frob,
        'commutator_norm': float(np.linalg.norm(C)),
        'relative_commutator_norm_over_kin':
            float(np.linalg.norm(C)) / kin_frob if kin_frob > 0 else 0.0,
    }


def _hottest(H_vee: np.ndarray, configs, orbitals):
    N = H_vee.shape[0]
    diag = np.diag(H_vee).copy()
    diag_idx = np.argsort(-np.abs(diag))
    diag_top = []
    for idx in diag_idx[:5]:
        (i, j) = configs[idx]
        diag_top.append({
            'pair': f'({orbitals[i]},{orbitals[j]})',
            'Vii': float(diag[idx]),
        })
    off = H_vee.copy()
    np.fill_diagonal(off, 0.0)
    i_hot, j_hot = np.unravel_index(np.argmax(np.abs(off)), off.shape)
    (a1, a2) = configs[i_hot]
    (b1, b2) = configs[j_hot]
    return {
        'diagonal_top5': diag_top,
        'hottest_offdiag': {
            'a': f'({orbitals[a1]},{orbitals[a2]})',
            'b': f'({orbitals[b1]},{orbitals[b2]})',
            'value': float(off[i_hot, j_hot]),
        },
        'concentration_ratio':
            float(np.abs(diag[diag_idx[0]])) / float(np.sum(np.abs(diag))),
    }


def run_ho(N_max: int, hw: float = 10.0):
    t0 = time.time()
    data = build_decomposed_ho_hamiltonians(N_max, hw)

    configs = data['configs']
    n_spatial = data['n_spatial']
    H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
    H_vee = data['H_vee_full']
    H_full = data['H_full']

    # Diagonalize full + kin-only.
    eigs_full, vecs_full = np.linalg.eigh(H_full)
    eigs_kin, vecs_kin = np.linalg.eigh(H_kin)
    ci_full = vecs_full[:, 0]
    ci_kin = vecs_kin[:, 0]

    rho_full = build_1rdm_from_singlet_ci(ci_full, configs, n_spatial)
    rho_kin = build_1rdm_from_singlet_ci(ci_kin, configs, n_spatial)
    ent_full = compute_entanglement_measures(rho_full)
    ent_kin = compute_entanglement_measures(rho_kin)

    frob = _frobenius_diagonality(H_vee, H_kin)
    hot = _hottest(H_vee, configs, data['orbitals'])

    wall = time.time() - t0
    return {
        'N_max': N_max,
        'hw_MeV': hw,
        'b_length_fm': data['b_length'],
        'n_spatial': n_spatial,
        'n_configs': data['n_configs'],
        'E_full_MeV': float(eigs_full[0]),
        'E_kin_MeV': float(eigs_kin[0]),
        'S_full_nats': float(ent_full['von_neumann_entropy']),
        'S_kin_nats': float(ent_kin['von_neumann_entropy']),
        'occ_full_top4': [float(x) for x in ent_full['occupation_numbers'][:4]],
        'occ_kin_top4': [float(x) for x in ent_kin['occupation_numbers'][:4]],
        'frobenius': frob,
        'hot_node': hot,
        'wall_seconds': wall,
    }


def main():
    results = {}
    for N_max in (2, 3):
        print(f"\n=== Bargmann-Segal HO, N_max={N_max} ===")
        r = run_ho(N_max, hw=10.0)
        results[f'HO_Nmax{N_max}_hw10'] = r
        print(f"  n_spatial = {r['n_spatial']}, n_configs = {r['n_configs']}")
        print(f"  E_full = {r['E_full_MeV']:.4f} MeV   "
              f"E_kin = {r['E_kin_MeV']:.4f} MeV")
        print(f"  S_full = {r['S_full_nats']:.4e} nats  "
              f"S_kin  = {r['S_kin_nats']:.4e} nats")
        print(f"  V_diag_fraction (H_kin eigenbasis) = "
              f"{r['frobenius']['V_diag_fraction_unsquared']:.4f}")
        print(f"  rel [H_kin,V]/H_kin = "
              f"{r['frobenius']['relative_commutator_norm_over_kin']:.4e}")
        print(f"  hottest diag : {r['hot_node']['diagonal_top5'][0]}")
        print(f"  concentration: {r['hot_node']['concentration_ratio']:.4f}")
        print(f"  wall = {r['wall_seconds']:.1f} s")

    # He-Coulomb reference from EP-1 (for the summary panel).
    results['_He_coulomb_reference_from_EP1'] = {
        'n_max2_S_full_nats': 0.03954502764024503,
        'n_max3_S_full_nats': 0.04081105136647098,
        'n_max3_V_diag_fraction_unsquared': 0.9202580563910042,
        'n_max3_hot_pair': '((1,0,0),(1,0,0))',
        'n_max3_hot_Vii': 0.625,
    }

    # EP-1 He concentration: |V_1s1s| / sum |V_ii| on the pair graph.
    # From debug/data/energy_graph_nmax3.json we can compute this offline.
    # Provide qualitative verdict here:
    for key, r in results.items():
        if not key.startswith('HO_'):
            continue
        s_he = (0.03954502764024503 if r['N_max'] == 2
                else 0.04081105136647098)
        ratio = r['S_full_nats'] / s_he if s_he > 0 else float('nan')
        print(f"\n  [verdict N_max={r['N_max']}] S_HO / S_He = {ratio:.3f}")
        if ratio < 0.5:
            print("    CONFIRMS Paper 27 §VII.A (S_HO << S_He)")
        elif ratio < 1.2:
            print("    COMPARABLE — weak evidence either way")
        else:
            print("    REFUTES Paper 27 §VII.A (S_HO >= S_He)")

    out = os.path.join(os.path.dirname(__file__), 'data',
                       'ep2b_ho_two_fermion.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved -> {out}")


if __name__ == '__main__':
    main()
