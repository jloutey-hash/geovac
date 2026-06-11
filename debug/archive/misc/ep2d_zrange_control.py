"""
Track EP-2d: Disambiguate the lone-pair A offset.

The EP-2c multi-block fit reported cores (α=2.37, A=7.78) vs lone-pairs
(α=2.24, A=4.65). BUT: both block-types drive the same He-like CI at
their own Z_eff. So the only difference is *which Z range* they sample:
  - cores:      Z ∈ {2, 3, 4, 6, 7, 8, 9, 10}
  - lone pairs: Z ∈ {5, 6, 7}

Test: fit an equally dense Z-grid on the He-like CI, then take the
same fit restricted to (a) Z ≤ 4, (b) Z ∈ [5,7], (c) Z ≥ 8. If A
varies across these sub-ranges by the same magnitude as the EP-2c
core-vs-lone-pair split, the offset is a log-log fit-range artifact
of the finite-basis higher-order residue, not a block-type signature.

A flat (α, A) across sub-ranges would mean the power law is truly
scale-free and the lone-pair offset reflects a real block-type
invariant.

Output: debug/data/ep2d_zrange_control.json
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from debug.energy_entanglement_decoupling import (
    build_decomposed_hamiltonians,
    solve_and_entangle,
)


def _row(Z, n_max=3):
    data = build_decomposed_hamiltonians(Z, n_max)
    H_kin = data['H_h1_diag'] + data['H_h1_offdiag']
    H_vee = data['H_vee_full']
    _, _, _, ent = solve_and_entangle(
        H_kin + H_vee, data['configs'], data['n_spatial'])
    S = float(ent['von_neumann_entropy'])
    _e, U = np.linalg.eigh(H_kin)
    V_in = U.T @ H_vee @ U
    V_frob = np.linalg.norm(V_in)
    V_diag = np.sqrt(np.sum(np.diag(V_in) ** 2))
    V_off = np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2))
    return S, float(V_off / np.linalg.norm(H_kin))


def _fit(w, S):
    lx = np.log(np.asarray(w, float))
    ly = np.log(np.asarray(S, float))
    alpha, intercept = np.polyfit(lx, ly, 1)
    pred = alpha * lx + intercept
    ss = np.sum((ly - pred) ** 2)
    st = np.sum((ly - ly.mean()) ** 2)
    return float(alpha), float(np.exp(intercept)), float(1.0 - ss / st)


def main():
    # Dense Z grid, same Hamiltonian family.
    Z_grid = list(range(2, 13))  # 2..12
    data = []
    for Z in Z_grid:
        S, w = _row(float(Z), n_max=3)
        data.append((Z, S, w))
        print(f"  Z={Z:2d}  S={S:.4e}  w_dim={w:.4e}")

    Zs = np.array([r[0] for r in data])
    Ss = np.array([r[1] for r in data])
    Ws = np.array([r[2] for r in data])

    print("\nFit on full range Z=2..12:")
    a, A, r2 = _fit(Ws, Ss)
    print(f"  alpha={a:.4f}  A={A:.4f}  R^2={r2:.6f}")

    bins = [
        ('low  Z<=4',   Zs <= 4),
        ('mid  5<=Z<=7', (Zs >= 5) & (Zs <= 7)),
        ('high Z>=8',   Zs >= 8),
        ('EP-2c cores mask (Z in {2,3,4,6,7,8,9,10})',
         np.isin(Zs, [2, 3, 4, 6, 7, 8, 9, 10])),
        ('EP-2c lone-pair mask (Z in {5,6,7})',
         np.isin(Zs, [5, 6, 7])),
    ]

    print("\nSub-range fits (identical Hamiltonian family):")
    bin_fits = {}
    for label, mask in bins:
        if mask.sum() < 2:
            continue
        a, A, r2 = _fit(Ws[mask], Ss[mask])
        bin_fits[label] = {'alpha': a, 'A': A, 'r2': r2, 'n': int(mask.sum())}
        print(f"  {label:50s}: n={int(mask.sum()):2d}  "
              f"alpha={a:.4f}  A={A:.4f}  R^2={r2:.4f}")

    # Headline test: does the EP-2c lone-pair mask reproduce A=4.65
    # on the pure He-like family (which is uncontaminated by block-type)?
    ep2c_lp = bin_fits['EP-2c lone-pair mask (Z in {5,6,7})']
    ep2c_co = bin_fits['EP-2c cores mask (Z in {2,3,4,6,7,8,9,10})']
    print("\n" + "=" * 72)
    print("Verdict on the EP-2c core/lone-pair offset:")
    print(f"  He-like at EP-2c core    Z set: A = {ep2c_co['A']:.4f}")
    print(f"  He-like at EP-2c lone-pair Z set: A = {ep2c_lp['A']:.4f}")
    print(f"  Difference: {100 * (ep2c_co['A'] - ep2c_lp['A']) / ep2c_co['A']:.1f}%")
    print("\nEP-2c multi-block reported core A=7.78 vs lone-pair A=4.65")
    print(f"    (40% offset). If this sub-range fit reproduces that same")
    print(f"    offset WITHOUT any block-type variation, the offset is a")
    print(f"    Z-range fit artifact, not a block-type invariant.")
    print("=" * 72)

    out = os.path.join(os.path.dirname(__file__), 'data',
                       'ep2d_zrange_control.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump({
            'Z_grid': Z_grid,
            'points': [{'Z': int(z), 'S_B': float(s), 'w_B_dimless': float(w)}
                       for z, s, w in data],
            'full_fit': {'alpha': float(_fit(Ws, Ss)[0]),
                         'A': float(_fit(Ws, Ss)[1]),
                         'r2': float(_fit(Ws, Ss)[2])},
            'bin_fits': bin_fits,
        }, f, indent=2)
    print(f"\nSaved -> {out}")


if __name__ == '__main__':
    main()
