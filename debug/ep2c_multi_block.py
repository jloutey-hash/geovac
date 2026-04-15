"""
Track EP-2c multi-block: extend the dimensionless entropy power law
across composed-molecule orbital blocks.

For each 2-electron block of each molecule:
  - Treat as an effective He-like system at Z = Z_center, max_n.
  - Build H_1 = h1(graph) + (Z_eff offdiag), V_ee = Slater integrals.
  - Diagonalize singlet FCI; compute S_B, w_B, w̃_B = w_B / ||H_1||_F.

Test whether A ≈ 8.16, α ≈ 2.383 from the He-like isoelectronic pilot
are block-universal constants or vary by block type.

Molecules sampled:
  - He-like Z sweep (reference, 6 points from EP-2c pilot)
  - LiH cores/bonds  (2 blocks)
  - HF  cores + lone pairs  (up to 4 blocks, skipping bond)
  - H2O cores + lone pairs  (up to 3 blocks)
  - NH3 cores + lone pairs  (up to 2 blocks)

Output: debug/data/ep2c_multi_block.json
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
from geovac.molecular_spec import (
    lih_spec, h2o_spec, nh3_spec, hf_spec,
)


def _w_B_and_S(Z_float: float, n_max: int) -> dict:
    """Build He-like CI at (Z, n_max) and return S_B, w̃_B, diagnostics."""
    data = build_decomposed_hamiltonians(Z_float, n_max)
    H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
    H_vee = data['H_vee_full']
    H_full = H_kin + H_vee
    _, _, _, ent = solve_and_entangle(
        H_full, data['configs'], data['n_spatial'])
    S_B = float(ent['von_neumann_entropy'])

    _eH, U = np.linalg.eigh(H_kin)
    V_in = U.T @ H_vee @ U
    V_frob = float(np.linalg.norm(V_in))
    V_diag_frob = float(np.sqrt(np.sum(np.diag(V_in) ** 2)))
    V_off = float(np.sqrt(max(0.0, V_frob ** 2 - V_diag_frob ** 2)))
    kin_frob = float(np.linalg.norm(H_kin))
    return {
        'S_B': S_B,
        'w_B': V_off,
        'w_B_dimless': V_off / kin_frob if kin_frob > 0 else 0.0,
        'V_frobenius': V_frob,
        'V_diag_fraction_unsquared': V_diag_frob / V_frob if V_frob > 0 else 0.0,
        'H_kin_frobenius': kin_frob,
        'n_configs': data['n_configs'],
    }


def _fit(xs, ys):
    lx = np.log(np.asarray(xs, float))
    ly = np.log(np.asarray(ys, float))
    alpha, intercept = np.polyfit(lx, ly, 1)
    pred = alpha * lx + intercept
    ss = np.sum((ly - pred) ** 2)
    st = np.sum((ly - ly.mean()) ** 2)
    return {
        'alpha': float(alpha),
        'A': float(np.exp(intercept)),
        'r2': float(1.0 - ss / st) if st > 0 else 0.0,
    }


def collect_molecule_blocks(spec, molecule_name, n_max_override=None):
    """Return list of (block_label, Z_eff, n_max, block_type) for 2e blocks.

    Skip 'bond' type blocks that have a hydrogen partner (not a pure
    single-center He-like subsystem).
    """
    rows = []
    for blk in spec.blocks:
        if blk.n_electrons != 2:
            continue
        # Only accept purely single-center blocks for this pilot.
        if blk.has_h_partner:
            continue
        n_max_use = n_max_override if n_max_override else blk.max_n
        # Need n_max >= 2 for cross-shell V_ee off-diagonal.
        if n_max_use < 2:
            n_max_use = 2
        rows.append({
            'molecule': molecule_name,
            'block_label': blk.label,
            'block_type': blk.block_type,
            'Z_eff': float(blk.Z_center),
            'n_max': n_max_use,
            'l_min': blk.l_min,
        })
    return rows


def main():
    print("=" * 72)
    print("Track EP-2c multi-block: dimensionless power law across blocks")
    print("=" * 72)

    # Layer 1: the EP-2c He-like reference.
    rows = []
    print("\nHe-like isoelectronic reference (n_max=3):")
    for Z in (2.0, 3.0, 4.0, 6.0, 8.0, 10.0):
        t0 = time.time()
        r = _w_B_and_S(Z, n_max=3)
        r.update({
            'molecule': 'He-like',
            'block_label': f'He-like_Z{int(Z)}',
            'block_type': 'core',
            'Z_eff': Z,
            'n_max': 3,
            'l_min': 0,
            'wall_seconds': time.time() - t0,
        })
        rows.append(r)
        print(f"  Z={Z:4.1f}  S_B={r['S_B']:.3e}  "
              f"w̃_B={r['w_B_dimless']:.3e}")

    # Layer 2: composed-molecule blocks (use block's max_n; n_max=3 for
    # core blocks where computationally tractable).
    mols = [
        ('LiH', lih_spec()),
        ('HF', hf_spec()),
        ('H2O', h2o_spec()),
        ('NH3', nh3_spec()),
    ]

    print("\nComposed-molecule 2e single-center blocks (n_max=3):")
    for name, spec in mols:
        meta = collect_molecule_blocks(spec, name, n_max_override=3)
        for m in meta:
            t0 = time.time()
            try:
                r = _w_B_and_S(m['Z_eff'], m['n_max'])
            except Exception as e:
                print(f"  SKIP {name}/{m['block_label']}: {e}")
                continue
            r.update(m)
            r['wall_seconds'] = time.time() - t0
            rows.append(r)
            print(f"  {name:4s}/{m['block_label']:12s}  "
                  f"type={m['block_type']:10s}  "
                  f"Z={m['Z_eff']:4.1f}  "
                  f"S_B={r['S_B']:.3e}  "
                  f"w̃_B={r['w_B_dimless']:.3e}")

    # Deduplicate: (Z, n_max) pairs are effectively the same CI.
    # For the fit, use all rows (many composed blocks duplicate the
    # He-like pilot at the same Z).
    S = [r['S_B'] for r in rows]
    w_dim = [r['w_B_dimless'] for r in rows]

    fit_all = _fit(w_dim, S)
    print("\n" + "-" * 72)
    print(f"Combined fit (all {len(rows)} blocks): "
          f"alpha = {fit_all['alpha']:.4f}, "
          f"A = {fit_all['A']:.4f}, "
          f"R² = {fit_all['r2']:.6f}")

    # Subset fits by block type
    by_type = {}
    for r in rows:
        by_type.setdefault(r['block_type'], []).append(r)
    subset_fits = {}
    print("\nBy block type:")
    for btype, group in by_type.items():
        if len(group) < 2:
            print(f"  {btype:10s}: n={len(group)}  (too few for fit)")
            continue
        f = _fit([r['w_B_dimless'] for r in group],
                 [r['S_B'] for r in group])
        subset_fits[btype] = f
        print(f"  {btype:10s}: n={len(group):2d}  "
              f"alpha={f['alpha']:.4f}  "
              f"A={f['A']:.4f}  "
              f"R²={f['r2']:.4f}")

    # Reference pilot alpha=2.383, A=8.16
    print("\nReference (He-like pilot, 6 Z values): alpha=2.383, A=8.16, R²=0.998")

    # Residual scatter test: do composed-block points fall on the
    # He-like line at their dimensionless w̃_B?
    A_ref = 8.16
    alpha_ref = 2.383
    max_rel_resid = 0.0
    for r in rows:
        if r['molecule'] == 'He-like':
            continue
        S_pred = A_ref * (r['w_B_dimless']) ** alpha_ref
        rel = abs(S_pred - r['S_B']) / r['S_B']
        max_rel_resid = max(max_rel_resid, rel)
    print(f"\nMax relative deviation of composed blocks from "
          f"He-like line: {max_rel_resid:.4f}")

    summary = {
        'rows': rows,
        'combined_fit': fit_all,
        'subset_fits': subset_fits,
        'he_like_reference': {'alpha': alpha_ref, 'A': A_ref, 'r2': 0.9982},
        'max_relative_residual_from_ref_line': max_rel_resid,
    }
    out = os.path.join(os.path.dirname(__file__), 'data',
                       'ep2c_multi_block.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved -> {out}")


if __name__ == '__main__':
    main()
