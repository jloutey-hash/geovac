"""
Track EP-1: Entropy partition for whole-subsystem entanglement.

Thesis: Entanglement entropy in GeoVac is generated entirely by V_ee,
not by the graph Laplacian kinetic piece. Paper 26 showed this at
single-orbital granularity. EP-1 extends to whole-subsystem entanglement
for He (and possibly LiH).

Method:
  (a) psi_full  = ground state of H_full = h1 + V_ee
  (b) psi_kin   = ground state of H_kin  = h1 alone (V_ee = 0)
  (c) S(rho_A) via 1-RDM von Neumann entropy (Paper 26 machinery)
  (d) Report S_full, S_kin, S_kin / S_full

Prediction: S_kin / S_full < 0.05.
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from debug.energy_entanglement_decoupling import (
    build_decomposed_hamiltonians,
    solve_and_entangle,
)


def run_he(Z_float: float, n_max: int):
    """Compute S_full, S_kin for He at given Z, n_max."""
    data = build_decomposed_hamiltonians(Z_float, n_max)

    # Full Hamiltonian: diagonal h1 + off-diagonal h1 + full V_ee
    H_full = data['H_h1_diag'] + data['H_h1_offdiag'] + data['H_vee_full']

    # Kinetic-only: diagonal h1 + off-diagonal h1 (V_ee = 0)
    H_kin = data['H_h1_diag'] + data['H_h1_offdiag']

    configs = data['configs']
    n_spatial = data['n_spatial']

    E_full, ci_full, rho_full, ent_full = solve_and_entangle(H_full, configs, n_spatial)
    E_kin, ci_kin, rho_kin, ent_kin = solve_and_entangle(H_kin, configs, n_spatial)

    S_full = ent_full['von_neumann_entropy']
    S_kin = ent_kin['von_neumann_entropy']

    # Also report occupation spectrum (for sanity)
    occ_full = ent_full['occupation_numbers']
    occ_kin = ent_kin['occupation_numbers']

    return {
        'Z': Z_float,
        'n_max': n_max,
        'n_configs': data['n_configs'],
        'n_spatial': n_spatial,
        'E_full': float(E_full),
        'E_kin': float(E_kin),
        'S_full': float(S_full),
        'S_kin': float(S_kin),
        'ratio_S_kin_over_S_full': float(S_kin / S_full) if S_full > 1e-14 else None,
        'occ_full_top4': [float(x) for x in occ_full[:4]],
        'occ_kin_top4': [float(x) for x in occ_kin[:4]],
    }


def main():
    results = {}

    print("=" * 70)
    print("Track EP-1: Entropy partition (kinetic-only vs full Hamiltonian)")
    print("=" * 70)

    for n_max in (2, 3):
        print(f"\n--- He (Z=2), n_max={n_max} ---")
        r = run_he(2.0, n_max)
        print(f"  n_configs = {r['n_configs']}, n_spatial = {r['n_spatial']}")
        print(f"  E_full = {r['E_full']:.6f} Ha")
        print(f"  E_kin  = {r['E_kin']:.6f} Ha  (V_ee = 0)")
        print(f"  S_full = {r['S_full']:.6e} nats")
        print(f"  S_kin  = {r['S_kin']:.6e} nats")
        if r['ratio_S_kin_over_S_full'] is not None:
            print(f"  ratio  = {r['ratio_S_kin_over_S_full']:.6e}")
        else:
            print(f"  ratio  = undefined (S_full near zero)")
        print(f"  top-4 occ (full): {r['occ_full_top4']}")
        print(f"  top-4 occ (kin):  {r['occ_kin_top4']}")
        results[f'He_nmax{n_max}'] = r

    # Save
    out_path = os.path.join(os.path.dirname(__file__), 'data',
                             'ep1_entropy_partition.json')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")

    # Verdict
    print("\n" + "=" * 70)
    print("Verdict summary:")
    print("=" * 70)
    for key, r in results.items():
        ratio = r['ratio_S_kin_over_S_full']
        if ratio is None:
            status = "CONFIRMED (S_kin = 0 exactly)"
        elif ratio < 0.05:
            status = "CONFIRMED"
        elif ratio < 0.30:
            status = "WEAKENED"
        else:
            status = "REFUTED"
        ratio_str = f"{ratio:.3e}" if ratio is not None else "0 (exact)"
        print(f"  {key}: S_kin/S_full = {ratio_str} -> {status}")


if __name__ == '__main__':
    main()
