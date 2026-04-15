"""
Track EP-2f: Full-library sweep of the dimensionless entropy power law.

Collects every UNIQUE (Z_eff, n_max) single-center 2e block across the
full composed-geometry molecular library (first, second, third rows
plus TM hydrides).  Deduplicates so each He-like CI is solved once.
Fits the combined dimensionless power law and tabulates the
log-log local slope d(log S)/d(log w̃) across the sampled Z.

Output: debug/data/ep2f_full_library.json
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from debug.energy_entanglement_decoupling import (
    build_decomposed_hamiltonians,
    solve_and_entangle,
)
from geovac import molecular_spec as ms


# All single-center spec factories in molecular_spec.py (skip
# multi-center factories in composed_qubit.py which have different API).
SPEC_FACTORIES = [
    ('LiH', ms.lih_spec),
    ('BeH2', ms.beh2_spec),
    ('CH4', ms.ch4_spec),
    ('NH3', ms.nh3_spec),
    ('H2O', ms.h2o_spec),
    ('HF',  ms.hf_spec),
    ('NaH', ms.nah_spec),
    ('MgH2', ms.mgh2_spec),
    ('SiH4', ms.sih4_spec),
    ('PH3', ms.ph3_spec),
    ('H2S', ms.h2s_spec),
    ('HCl', ms.hcl_spec),
    ('KH',  ms.kh_spec),
    ('CaH2', ms.cah2_spec),
    ('GeH4', ms.geh4_spec),
    ('AsH3', ms.ash3_spec),
    ('H2Se', ms.h2se_spec),
    ('HBr',  ms.hbr_spec),
]


def _row(Z, n_max=3):
    data = build_decomposed_hamiltonians(Z, n_max)
    H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
    H_vee = data['H_vee_full']
    _, _, _, ent = solve_and_entangle(
        H_kin + H_vee, data['configs'], data['n_spatial'])
    _e, U = np.linalg.eigh(H_kin)
    V_in = U.T @ H_vee @ U
    V_frob = np.linalg.norm(V_in)
    V_diag = np.sqrt(np.sum(np.diag(V_in) ** 2))
    V_off = np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2))
    return {
        'S_B': float(ent['von_neumann_entropy']),
        'w_B_dimless': float(V_off / np.linalg.norm(H_kin)),
        'V_frobenius': float(V_frob),
        'V_diag_fraction': float(V_diag / V_frob) if V_frob > 0 else 0.0,
    }


def _fit(Zs, Ss, Ws):
    Ws = np.asarray(Ws, float); Ss = np.asarray(Ss, float)
    if len(Ws) < 2:
        return {'alpha': None, 'A': None, 'r2': None, 'n': int(len(Ws))}
    lx = np.log(Ws); ly = np.log(Ss)
    alpha, intercept = np.polyfit(lx, ly, 1)
    pred = alpha * lx + intercept
    ss = np.sum((ly - pred) ** 2)
    st = np.sum((ly - ly.mean()) ** 2)
    return {
        'alpha': float(alpha),
        'A': float(np.exp(intercept)),
        'r2': float(1.0 - ss / st) if st > 0 else 0.0,
        'n': int(len(Ws)),
    }


def main():
    n_max_analysis = 3
    print("=" * 72)
    print(f"EP-2f full-library sweep (n_max={n_max_analysis})")
    print("=" * 72)

    # Collect (Z, mol, block) for every 2e single-center block.
    block_index = []
    unique_Zs = set()
    for mol_name, factory in SPEC_FACTORIES:
        try:
            spec = factory()
        except Exception as e:
            print(f"  skip {mol_name}: {e}")
            continue
        for blk in spec.blocks:
            if blk.n_electrons != 2 or blk.has_h_partner:
                continue
            Z = float(blk.Z_center)
            block_index.append({
                'molecule': mol_name,
                'block_label': blk.label,
                'block_type': blk.block_type,
                'Z_eff': Z,
                'n_max_native': blk.max_n,
                'l_min': blk.l_min,
            })
            unique_Zs.add(Z)

    # Include He-like atomic reference + high-Z extension to probe the
    # asymptotic alpha -> 2 limit (second-order RS).
    for Z in (2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
              12.0, 15.0, 20.0, 30.0, 50.0, 100.0):
        unique_Zs.add(Z)

    # Evaluate each unique Z once.
    Z_to_point = {}
    print(f"\nUnique Z_eff values: {sorted(unique_Zs)}")
    for Z in sorted(unique_Zs):
        t0 = time.time()
        Z_to_point[Z] = _row(Z, n_max=n_max_analysis)
        Z_to_point[Z]['wall_seconds'] = time.time() - t0
        print(f"  Z={Z:5.2f}  "
              f"S={Z_to_point[Z]['S_B']:.4e}  "
              f"w={Z_to_point[Z]['w_B_dimless']:.4e}  "
              f"({Z_to_point[Z]['wall_seconds']:.2f}s)")

    # Emit full block table by mapping each block to its Z point.
    block_rows = []
    for b in block_index:
        p = Z_to_point[b['Z_eff']]
        row = {**b, **p}
        block_rows.append(row)

    # Combined fit over unique Z points (not weighted by block count).
    Zs_sorted = sorted(unique_Zs)
    Ss_u = [Z_to_point[z]['S_B'] for z in Zs_sorted]
    Ws_u = [Z_to_point[z]['w_B_dimless'] for z in Zs_sorted]
    combined_fit = _fit(Zs_sorted, Ss_u, Ws_u)

    # Local log-log slope at each Z.
    local_slopes = []
    lx = np.log(Ws_u); ly = np.log(Ss_u)
    for i in range(len(Zs_sorted)):
        if 0 < i < len(Zs_sorted) - 1:
            slope = (ly[i + 1] - ly[i - 1]) / (lx[i + 1] - lx[i - 1])
        elif i == 0:
            slope = (ly[1] - ly[0]) / (lx[1] - lx[0])
        else:
            slope = (ly[-1] - ly[-2]) / (lx[-1] - lx[-2])
        local_slopes.append((float(Zs_sorted[i]), float(slope)))

    print("\n" + "-" * 72)
    print(f"Combined fit (unique Z, {combined_fit['n']} points): "
          f"alpha={combined_fit['alpha']:.4f}, "
          f"A={combined_fit['A']:.4f}, "
          f"R2={combined_fit['r2']:.6f}")

    print("\nLocal log-log slope alpha_loc(Z):")
    print("  Z         alpha_loc")
    for Z, sl in local_slopes:
        print(f"  {Z:5.2f}     {sl:.4f}")

    # Split into low/high Z for the convergence narrative
    low_fit = _fit(
        [z for z in Zs_sorted if z <= 4.5],
        [Z_to_point[z]['S_B'] for z in Zs_sorted if z <= 4.5],
        [Z_to_point[z]['w_B_dimless'] for z in Zs_sorted if z <= 4.5])
    mid_fit = _fit(
        [z for z in Zs_sorted if 5 <= z <= 10],
        [Z_to_point[z]['S_B'] for z in Zs_sorted if 5 <= z <= 10],
        [Z_to_point[z]['w_B_dimless'] for z in Zs_sorted if 5 <= z <= 10])
    high_fit = _fit(
        [z for z in Zs_sorted if z > 10],
        [Z_to_point[z]['S_B'] for z in Zs_sorted if z > 10],
        [Z_to_point[z]['w_B_dimless'] for z in Zs_sorted if z > 10])

    print("\nBand fits:")
    for tag, f in [('Z<=4.5', low_fit), ('5<=Z<=10', mid_fit), ('Z>10', high_fit)]:
        if f['n'] < 2 or f['alpha'] is None:
            print(f"  {tag:10s}: n={f['n']:2d}  (too few)")
            continue
        print(f"  {tag:10s}: n={f['n']:2d}  "
              f"alpha={f['alpha']:.4f}  A={f['A']:.4f}  R2={f['r2']:.4f}")

    # Library block coverage stats
    by_mol = {}
    for b in block_index:
        by_mol.setdefault(b['molecule'], 0)
        by_mol[b['molecule']] += 1
    print(f"\nLibrary coverage: {len(block_index)} single-center 2e blocks "
          f"across {len(by_mol)} molecules")

    out = os.path.join(os.path.dirname(__file__), 'data',
                       'ep2f_full_library.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump({
            'unique_Z_points': [{'Z': z, **Z_to_point[z]} for z in Zs_sorted],
            'block_rows': block_rows,
            'combined_fit_unique_Z': combined_fit,
            'band_fits': {
                'Z<=4.5': low_fit, '5<=Z<=10': mid_fit, 'Z>10': high_fit,
            },
            'local_log_log_slope': local_slopes,
            'library_coverage': by_mol,
        }, f, indent=2)
    print(f"\nSaved -> {out}")


if __name__ == '__main__':
    main()
