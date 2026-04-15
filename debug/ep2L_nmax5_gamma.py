"""
Track EP-2L: extend the gamma(Z, n_max) table to n_max=5.

EP-2i tested n_max=2,3,4 and showed local slope (Z=15->30) converges to
2 from above with a mild monotone n_max downshift. n_max=4 gave 1.966.
Push to n_max=5 to test whether the trend continues toward gamma=2 or
plateaus below it.

For He-like at n_max=5: 55 spatial orbitals, ~few hundred M_L=0 singlet
configs. Should run in under 60s per Z-point.

Output: debug/data/ep2L_nmax5_gamma.json
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from debug.energy_entanglement_decoupling import (
    build_decomposed_hamiltonians, solve_and_entangle,
)


def _row(Z, n_max):
    t0 = time.time()
    data = build_decomposed_hamiltonians(Z, n_max)
    H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
    H_vee = data['H_vee_full']
    _, _, _, ent = solve_and_entangle(
        H_kin + H_vee, data['configs'], data['n_spatial'])
    e, U = np.linalg.eigh(H_kin)
    V_in = U.T @ H_vee @ U
    V_frob = np.linalg.norm(V_in)
    V_diag = np.sqrt(np.sum(np.diag(V_in) ** 2))
    V_off = float(np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2)))
    kin = float(np.linalg.norm(H_kin))
    e_sort = np.sort(e)
    return {
        'Z': Z, 'n_max': n_max,
        'S_B': float(ent['von_neumann_entropy']),
        'w_B_dimless': V_off / kin,
        'delta_B': float(abs(e_sort[1] - e_sort[0]) / kin),
        'n_configs': data['n_configs'],
        'wall_seconds': time.time() - t0,
    }


def main():
    print("EP-2L: gamma local slope at n_max=5")
    print("=" * 60)

    rows = []
    Zs = [2.0, 3.0, 4.0, 6.0, 10.0, 15.0, 30.0]
    for n_max in (2, 3, 4, 5):
        print(f"\n--- n_max = {n_max} ---")
        for Z in Zs:
            r = _row(Z, n_max)
            rows.append(r)
            ratio = r['w_B_dimless'] / r['delta_B']
            print(f"  Z={Z:5.1f}  S={r['S_B']:.3e}  w/d={ratio:.3e}  "
                  f"({r['wall_seconds']:.2f}s, {r['n_configs']} cfg)")

    # Local slope at large Z
    summary = {}
    print("\n" + "-" * 60)
    print("Local slope (Z=15 -> 30):")
    for n_max in (2, 3, 4, 5):
        sub = [r for r in rows if r['n_max'] == n_max]
        sub_15 = [r for r in sub if r['Z'] == 15.0][0]
        sub_30 = [r for r in sub if r['Z'] == 30.0][0]
        x1 = sub_15['w_B_dimless'] / sub_15['delta_B']
        x2 = sub_30['w_B_dimless'] / sub_30['delta_B']
        slope = (np.log(sub_30['S_B']) - np.log(sub_15['S_B'])) / \
                (np.log(x2) - np.log(x1))
        summary[n_max] = float(slope)
        print(f"  n_max={n_max}: gamma_loc = {slope:.4f}")

    diffs = [summary[n+1] - summary[n] for n in (2, 3, 4)]
    print(f"\nMonotonicity (n_max+1 minus n_max): {diffs}")
    print(f"All decreasing? {all(d < 0 for d in diffs)}")

    out = {'rows': rows, 'local_slopes_Z15_to_30': summary}
    out_path = os.path.join(os.path.dirname(__file__), 'data',
                            'ep2L_nmax5_gamma.json')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved -> {out_path}")


if __name__ == '__main__':
    main()
