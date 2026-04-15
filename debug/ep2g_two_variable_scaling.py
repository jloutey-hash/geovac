"""
Track EP-2g: Two-variable scaling S_B = f(w̃_B, δ_B).

EP-2e showed that the single-variable predictor w̃_B = ||V_ee_off||_F
/ ||H_1||_F fails for bond blocks by 20x because near-degeneracy of
the two-center H_1 eigenstates pushes V_ee mixing into the
quasi-degenerate regime.

Hypothesis: adding the smallest H_1 gap
    δ_B = ΔE_1 / ||H_1||_F
as a second predictor collapses single-center and bond data onto a
universal surface.  Two candidate forms:
  (A) S_B = A * (w̃_B / δ_B)^γ        (Brillouin-Wigner-like ratio)
  (B) S_B = A * w̃_B^α * f(w̃_B / δ_B)  (regime interpolation)

Sample set:
  - He-like Z ∈ {2..10, 15, 30, 100} (single-center, varying δ)
  - Composed single-center 2e blocks (cores + lone pairs, matched to Z)
  - LiH bond block at R ∈ {2, 3, 4, 6} bohr — VARIES δ systematically
    at fixed orbital structure.  This is the key diagnostic: if
    (w̃, δ) works, all R-points should fall on the universal surface.

Output: debug/data/ep2g_two_variable_scaling.json
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
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import lih_spec


def _single_center_row(Z, n_max=3):
    data = build_decomposed_hamiltonians(Z, n_max)
    H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
    H_vee = data['H_vee_full']
    H_full = H_kin + H_vee
    _, _, _, ent = solve_and_entangle(
        H_full, data['configs'], data['n_spatial'])
    S_B = float(ent['von_neumann_entropy'])
    e, U = np.linalg.eigh(H_kin)
    V_in = U.T @ H_vee @ U
    V_frob = float(np.linalg.norm(V_in))
    V_diag = float(np.sqrt(np.sum(np.diag(V_in) ** 2)))
    V_off = float(np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2)))
    kin_frob = float(np.linalg.norm(H_kin))
    # Config-space H_kin gap
    e_sorted = np.sort(e)
    dE1 = float(e_sorted[1] - e_sorted[0])
    return {
        'S_B': S_B,
        'w_B_dimless': V_off / kin_frob,
        'delta_B': abs(dE1) / kin_frob,
        'gap_absolute': abs(dE1),
        'kin_frobenius': kin_frob,
    }


def _enumerate_states(max_n, l_min=0):
    states = []
    for n in range(1, max_n + 1):
        for l in range(l_min, n):
            for m in range(-l, l + 1):
                states.append((n, l, m))
    return states


def _bond_block_row(R):
    """Extract LiH bond-block 2e FCI at bond length R."""
    spec = lih_spec(R=R)
    result = build_balanced_hamiltonian(spec, R=R, verbose=False)
    h1 = result['h1']
    eri = result['eri']
    blocks = spec.blocks

    # Compute orbital offsets and m quantum numbers.
    orb_idx, orb_m, is_bond = [], [], []
    offset = 0
    for blk in blocks:
        center_states = _enumerate_states(blk.max_n, l_min=blk.l_min)
        for (n, l, m) in center_states:
            orb_idx.append(offset); orb_m.append(m)
            is_bond.append('bond' in blk.label)
            offset += 1
        if blk.has_h_partner:
            partner_n = blk.max_n_partner or blk.max_n
            for (n, l, m) in _enumerate_states(partner_n):
                orb_idx.append(offset); orb_m.append(m)
                is_bond.append('bond' in blk.label)
                offset += 1

    bond_orb = [i for i, b in zip(orb_idx, is_bond) if b]
    bond_m = [orb_m[i] for i, b in enumerate(is_bond) if b]
    h1_sub = h1[np.ix_(bond_orb, bond_orb)]
    eri_sub = eri[np.ix_(bond_orb, bond_orb, bond_orb, bond_orb)]

    # 2e singlet FCI on M_L=0 configs.
    configs = [(i, j) for i in range(len(bond_orb))
               for j in range(i, len(bond_orb))
               if bond_m[i] + bond_m[j] == 0]
    n_cfg = len(configs)
    H_full = np.zeros((n_cfg, n_cfg))
    H_kin = np.zeros((n_cfg, n_cfg))
    for I, (i, j) in enumerate(configs):
        for J in range(I, n_cfg):
            p, q = configs[J]
            I_p = [(i, j)] + ([(j, i)] if i != j else [])
            J_p = [(p, q)] + ([(q, p)] if p != q else [])
            N_I = np.sqrt(float(len(I_p)))
            N_J = np.sqrt(float(len(J_p)))
            me_full = 0.0; me_kin = 0.0
            for a, b in I_p:
                for c, d in J_p:
                    if b == d:
                        me_full += h1_sub[a, c]; me_kin += h1_sub[a, c]
                    if a == c:
                        me_full += h1_sub[b, d]; me_kin += h1_sub[b, d]
                    me_full += eri_sub[a, c, b, d]
            me_full /= (N_I * N_J); me_kin /= (N_I * N_J)
            H_full[I, J] = H_full[J, I] = me_full
            H_kin[I, J] = H_kin[J, I] = me_kin
    V_ee = H_full - H_kin

    eigs_full, vecs_full = np.linalg.eigh(H_full)
    eigs_kin, _ = np.linalg.eigh(H_kin)

    # 1-RDM and entropy
    from debug.entanglement_geometry import (
        build_1rdm_from_singlet_ci, compute_entanglement_measures,
    )
    rho = build_1rdm_from_singlet_ci(vecs_full[:, 0], configs, len(bond_orb))
    ent = compute_entanglement_measures(rho)
    S_B = float(ent['von_neumann_entropy'])

    _e, U = np.linalg.eigh(H_kin)
    V_in = U.T @ V_ee @ U
    V_frob = float(np.linalg.norm(V_in))
    V_diag = float(np.sqrt(np.sum(np.diag(V_in) ** 2)))
    V_off = float(np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2)))
    kin_frob = float(np.linalg.norm(H_kin))
    e_sorted = np.sort(_e)
    dE1 = float(e_sorted[1] - e_sorted[0])
    return {
        'S_B': S_B,
        'w_B_dimless': V_off / kin_frob,
        'delta_B': abs(dE1) / kin_frob,
        'gap_absolute': abs(dE1),
        'kin_frobenius': kin_frob,
        'R': float(R),
        'E_full': float(eigs_full[0]),
    }


def _fit_power(xs, ys):
    lx = np.log(np.asarray(xs, float))
    ly = np.log(np.asarray(ys, float))
    a, b = np.polyfit(lx, ly, 1)
    pred = a * lx + b
    ss = float(np.sum((ly - pred) ** 2))
    st = float(np.sum((ly - ly.mean()) ** 2))
    return float(a), float(np.exp(b)), float(1 - ss / st) if st > 0 else 0.0


def main():
    rows = []
    print("=" * 72)
    print("EP-2g: Two-variable scaling -- S_B vs (w_B, delta_B)")
    print("=" * 72)

    print("\nSingle-center He-like Z-sweep:")
    for Z in (2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
              12.0, 15.0, 20.0, 30.0, 50.0, 100.0):
        r = _single_center_row(Z, n_max=3)
        r['system'] = f'He-like_Z{Z:g}'
        r['kind'] = 'single_center'
        rows.append(r)
        print(f"  Z={Z:6.1f}  S={r['S_B']:.3e}  w={r['w_B_dimless']:.3e}  "
              f"d={r['delta_B']:.3e}  w/d={r['w_B_dimless']/r['delta_B']:.3f}")

    print("\nLiH bond block at varying R (two-center H_1):")
    for R in (1.5, 2.0, 3.015, 4.0, 5.0, 6.0, 8.0):
        try:
            r = _bond_block_row(R)
        except Exception as e:
            print(f"  R={R}: FAILED {e}")
            continue
        r['system'] = f'LiH_bond_R{R:g}'
        r['kind'] = 'bond'
        rows.append(r)
        print(f"  R={R:5.2f}  S={r['S_B']:.3e}  w={r['w_B_dimless']:.3e}  "
              f"d={r['delta_B']:.3e}  w/d={r['w_B_dimless']/r['delta_B']:.3f}")

    # Test three scaling laws
    sc_rows = [r for r in rows if r['kind'] == 'single_center']
    bond_rows = [r for r in rows if r['kind'] == 'bond']

    print("\n" + "-" * 72)
    print("Fit 1: S = A * w^alpha  (single-variable, EP-2c/f)")
    for name, grp in [('single-center', sc_rows), ('bond', bond_rows),
                      ('all', rows)]:
        if len(grp) < 2:
            continue
        a, A, r2 = _fit_power([r['w_B_dimless'] for r in grp],
                              [r['S_B'] for r in grp])
        print(f"  {name:14s} n={len(grp):2d}  "
              f"alpha={a:.4f}  A={A:.4f}  R^2={r2:.4f}")

    print("\nFit 2: S = A * (w/delta)^gamma  (Brillouin-Wigner ratio)")
    for name, grp in [('single-center', sc_rows), ('bond', bond_rows),
                      ('all', rows)]:
        if len(grp) < 2:
            continue
        ratio = [r['w_B_dimless'] / r['delta_B'] for r in grp]
        a, A, r2 = _fit_power(ratio, [r['S_B'] for r in grp])
        print(f"  {name:14s} n={len(grp):2d}  "
              f"gamma={a:.4f}  A={A:.4f}  R^2={r2:.4f}")

    # Fit 3: Multilinear in log-log with both features.
    # log S = alpha * log w̃ - beta * log δ + const
    print("\nFit 3: log S = alpha*log(w) + beta*log(delta) + const (2D)")
    for name, grp in [('single-center', sc_rows), ('bond', bond_rows),
                      ('all', rows)]:
        if len(grp) < 3:
            continue
        lw = np.log([r['w_B_dimless'] for r in grp])
        ld = np.log([r['delta_B'] for r in grp])
        ls = np.log([r['S_B'] for r in grp])
        X = np.column_stack([lw, ld, np.ones_like(lw)])
        coef, *_ = np.linalg.lstsq(X, ls, rcond=None)
        pred = X @ coef
        ss = float(np.sum((ls - pred) ** 2))
        st = float(np.sum((ls - ls.mean()) ** 2))
        r2 = 1 - ss / st if st > 0 else 0.0
        print(f"  {name:14s} n={len(grp):2d}  "
              f"alpha={coef[0]:.4f}  beta={coef[1]:.4f}  "
              f"logA={coef[2]:.4f}  R^2={r2:.4f}")

    # Does the bond block land on the (w̃/δ) universal line?
    ratio_all = [r['w_B_dimless'] / r['delta_B'] for r in rows]
    S_all = [r['S_B'] for r in rows]
    g_all, A_all, r2_all = _fit_power(ratio_all, S_all)
    print("\n" + "=" * 72)
    print(f"UNIVERSAL fit (single-center + bond, {len(rows)} points):")
    print(f"  S = {A_all:.4f} * (w/delta)^{g_all:.4f}, R^2 = {r2_all:.6f}")
    max_rel = 0.0
    for r in rows:
        S_pred = A_all * (r['w_B_dimless'] / r['delta_B']) ** g_all
        rel = abs(S_pred - r['S_B']) / r['S_B']
        if rel > max_rel:
            max_rel = rel
    print(f"  Max relative deviation: {max_rel:.4f}")
    print("=" * 72)

    out = os.path.join(os.path.dirname(__file__), 'data',
                       'ep2g_two_variable_scaling.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump({
            'rows': rows,
            'universal_w_over_delta': {
                'gamma': g_all, 'A': A_all, 'r2': r2_all,
                'max_rel_deviation': max_rel,
            },
        }, f, indent=2)
    print(f"\nSaved -> {out}")


if __name__ == '__main__':
    main()
