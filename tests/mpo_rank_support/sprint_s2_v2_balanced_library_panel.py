"""Sprint S2-v2 — Extended negative control: balanced-coupled library panel.

Sprint S1 showed that LiH composed has chi_k = 2 at every sub-block
boundary, and that LiH balanced-coupled lifts this to chi_k = 9 at both
boundaries (bit-identical). The structural reading: chi=2 is tied to
F4 (cross-block ERI vanishing) in composed; balanced-coupled has
non-zero cross-block ERIs (cross-center V_ne via multipole) which lift
the boundary bond rank.

This panel extends the negative control to BeH2 and H2O balanced-coupled
to confirm:
  (i)  chi_k > 2 at sub-block boundaries (structural reading robust).
  (ii) The lift (chi_boundary_balanced - 2) is consistent across
       different molecules with the same orbital basis.

For each molecule we record the chi_k profile across all cuts and
flag the sub-block boundaries.
"""

from __future__ import annotations

import json
import os
from typing import List, Dict
import numpy as np

from openfermion import QubitOperator

from geovac.composed_qubit import (
    build_composed_hamiltonian, _enumerate_states,
)
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import lih_spec, beh2_spec, h2o_spec


def operator_schmidt_rank(qop: QubitOperator, cut: int, n_qubits: int,
                          rel_thr: float = 1e-12) -> int:
    coef = {}
    for pauli_tuple, c in qop.terms.items():
        left = tuple((q, op) for (q, op) in pauli_tuple if q < cut)
        right = tuple((q, op) for (q, op) in pauli_tuple if q >= cut)
        coef[(left, right)] = c
    left_keys = sorted(set(k[0] for k in coef))
    right_keys = sorted(set(k[1] for k in coef))
    if not left_keys or not right_keys:
        return 0
    left_idx = {k: i for i, k in enumerate(left_keys)}
    right_idx = {k: i for i, k in enumerate(right_keys)}
    M = np.zeros((len(left_keys), len(right_keys)))
    for (l, r), v in coef.items():
        M[left_idx[l], right_idx[r]] = float(v.real if hasattr(v, 'real') else v)
    sv = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(sv > rel_thr * max(sv[0], 1e-30)))


def get_subblock_boundaries(spec) -> List[int]:
    boundaries = []
    cumulative = 0
    for blk in spec.blocks:
        l_min = getattr(blk, 'l_min', 0)
        cs = _enumerate_states(blk.max_n, l_min=l_min)
        cumulative += len(cs)
        boundaries.append(2 * cumulative)
        if blk.has_h_partner:
            pn = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            ps = _enumerate_states(pn)
            cumulative += len(ps)
            boundaries.append(2 * cumulative)
    return boundaries[:-1]


def analyze_pair(name: str, spec_fn, R: float):
    spec = spec_fn()
    print(f"\n{'='*70}\n{name}: R = {R} bohr\n{'='*70}")

    composed = build_composed_hamiltonian(spec, verbose=False)
    H_c = composed['qubit_op']
    M = composed['h1'].shape[0]
    Q = 2 * M
    print(f"  composed: M={M}, Q={Q}, Pauli={len(H_c.terms)}")

    balanced = build_balanced_hamiltonian(spec, R=R, verbose=False)
    H_b = balanced['qubit_op']
    print(f"  balanced: Pauli={len(H_b.terms)}")

    boundaries = get_subblock_boundaries(spec)
    print(f"  sub-block boundaries: {boundaries}")

    boundary_data = []
    for b in boundaries:
        chi_c = operator_schmidt_rank(H_c, b, Q)
        chi_b = operator_schmidt_rank(H_b, b, Q)
        lift = chi_b - chi_c
        boundary_data.append({
            'cut': b, 'chi_composed': chi_c, 'chi_balanced': chi_b,
            'lift': lift,
        })
        print(f"    cut {b}: chi_composed={chi_c}, chi_balanced={chi_b}, "
              f"lift={lift}")

    return {
        'molecule': name, 'M': M, 'Q': Q, 'R': R,
        'n_pauli_composed': len(H_c.terms),
        'n_pauli_balanced': len(H_b.terms),
        'sub_block_boundaries': boundaries,
        'boundary_data': boundary_data,
    }


def main():
    results = {}
    cases = [
        ('LiH', lih_spec, 3.015),
        ('BeH2', beh2_spec, 2.539),
        ('H2O', h2o_spec, 1.808),
    ]
    for name, fn, R in cases:
        try:
            results[name] = analyze_pair(name, fn, R)
        except Exception as e:
            import traceback; traceback.print_exc()
            results[name] = {'error': str(e)}

    # Print summary table
    print(f"\n{'='*70}\nNEGATIVE CONTROL SUMMARY\n{'='*70}")
    print(f"{'mol':>6s} {'cut':>4s} {'chi_comp':>9s} {'chi_bal':>9s} "
          f"{'lift':>6s} {'lift/sb':>9s}")
    print('-' * 60)
    for name, r in results.items():
        if 'error' in r:
            print(f"  {name}: ERROR")
            continue
        n_sb_bd = len(r['boundary_data'])
        for b in r['boundary_data']:
            sb_idx = r['boundary_data'].index(b)
            print(f"{name:>6s} {b['cut']:>4d} {b['chi_composed']:>9d} "
                  f"{b['chi_balanced']:>9d} {b['lift']:>6d} "
                  f"{'.' if sb_idx == 0 else '':>9s}")

    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'sprint_s2_v2_balanced_library_panel.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved to {out_path}")


if __name__ == '__main__':
    main()
