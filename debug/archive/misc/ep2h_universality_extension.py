"""
Track EP-2h: three simultaneous tests of the EP-2g universal curve.

Q1 (other bonds): does the same (w/delta) universal curve describe
    BeH2 and H2O bonds?  Same gamma ~ 2.23?

Q2 (n_max convergence): does the single-center exponent drift toward
    2 as n_max increases (more orbitals → better non-degenerate RS)?

Q3 (bounded form): is S = log(2) * tanh^2(c * (w/delta)^gamma) a
    better fit than the pure power law, given the physical log(2)
    ceiling on two-orbital mixing?

Output: debug/data/ep2h_universality_extension.json
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from debug.energy_entanglement_decoupling import (
    build_decomposed_hamiltonians, solve_and_entangle,
)
from debug.entanglement_geometry import (
    build_1rdm_from_singlet_ci, compute_entanglement_measures,
)
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import lih_spec, beh2_spec, h2o_spec, hf_spec


def _single_center_row(Z, n_max):
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
    kin_frob = float(np.linalg.norm(H_kin))
    e_sort = np.sort(e)
    return {
        'S_B': float(ent['von_neumann_entropy']),
        'w_B_dimless': V_off / kin_frob,
        'delta_B': float(abs(e_sort[1] - e_sort[0]) / kin_frob),
        'kin_frobenius': kin_frob,
        'n_max': n_max,
    }


def _enumerate_states(max_n, l_min=0):
    return [(n, l, m) for n in range(1, max_n + 1)
            for l in range(l_min, n) for m in range(-l, l + 1)]


def _bond_block_row(spec_fn, R, bond_block_index=0, max_n=None, verbose=False):
    """Extract the 2e FCI for one bond block of a composed molecule."""
    kwargs = {'R': R} if R is not None else {}
    if max_n is not None:
        kwargs['max_n'] = max_n
    spec = spec_fn(**kwargs)
    result = build_balanced_hamiltonian(spec, R=R if R is not None else 3.015,
                                        verbose=verbose)
    h1 = result['h1']
    eri = result['eri']
    blocks = spec.blocks

    # Identify target bond block.
    bond_indices = [i for i, blk in enumerate(blocks) if blk.block_type == 'bond']
    target_block_idx = bond_indices[bond_block_index]
    target_label = blocks[target_block_idx].label

    orb_m = []
    orb_in_bond = []
    offset = 0
    target_orb_idx = []
    for i_blk, blk in enumerate(blocks):
        center_states = _enumerate_states(blk.max_n, l_min=blk.l_min)
        for (n, l, m) in center_states:
            if i_blk == target_block_idx:
                target_orb_idx.append(offset)
                orb_m.append(m)
            offset += 1
        if blk.has_h_partner:
            partner_n = blk.max_n_partner or blk.max_n
            for (n, l, m) in _enumerate_states(partner_n):
                if i_blk == target_block_idx:
                    target_orb_idx.append(offset)
                    orb_m.append(m)
                offset += 1

    h1_sub = h1[np.ix_(target_orb_idx, target_orb_idx)]
    eri_sub = eri[np.ix_(target_orb_idx, target_orb_idx,
                         target_orb_idx, target_orb_idx)]

    # 2e singlet FCI, M_L=0.
    configs = [(i, j) for i in range(len(target_orb_idx))
               for j in range(i, len(target_orb_idx))
               if orb_m[i] + orb_m[j] == 0]
    n_cfg = len(configs)
    H_full = np.zeros((n_cfg, n_cfg)); H_kin = np.zeros((n_cfg, n_cfg))
    for I, (i, j) in enumerate(configs):
        for J in range(I, n_cfg):
            p, q = configs[J]
            I_p = [(i, j)] + ([(j, i)] if i != j else [])
            J_p = [(p, q)] + ([(q, p)] if p != q else [])
            N_I = np.sqrt(float(len(I_p))); N_J = np.sqrt(float(len(J_p)))
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
            H_kin[I, J]  = H_kin[J, I]  = me_kin

    V_ee = H_full - H_kin
    eigs, vecs = np.linalg.eigh(H_full)
    rho = build_1rdm_from_singlet_ci(vecs[:, 0], configs, len(target_orb_idx))
    S_B = float(compute_entanglement_measures(rho)['von_neumann_entropy'])

    e, U = np.linalg.eigh(H_kin)
    V_in = U.T @ V_ee @ U
    V_frob = np.linalg.norm(V_in)
    V_diag = np.sqrt(np.sum(np.diag(V_in) ** 2))
    V_off = float(np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2)))
    kin_frob = float(np.linalg.norm(H_kin))
    e_sort = np.sort(e)
    return {
        'S_B': S_B,
        'w_B_dimless': V_off / kin_frob,
        'delta_B': float(abs(e_sort[1] - e_sort[0]) / kin_frob),
        'kin_frobenius': kin_frob,
        'n_orb': len(target_orb_idx),
        'n_cfg': n_cfg,
        'block_label': target_label,
        'R': float(R) if R is not None else None,
        'max_n': max_n,
    }


def _fit_power(xs, ys):
    lx = np.log(np.asarray(xs, float)); ly = np.log(np.asarray(ys, float))
    a, b = np.polyfit(lx, ly, 1)
    pred = a * lx + b
    ss = float(np.sum((ly - pred) ** 2))
    st = float(np.sum((ly - ly.mean()) ** 2))
    return {'gamma': float(a), 'A': float(np.exp(b)),
            'r2': float(1 - ss / st) if st > 0 else 0.0, 'n': len(xs)}


def _fit_bounded(xs, ys, log2=np.log(2.0)):
    """Fit S = log(2) * tanh^2(c * x^gamma)."""
    def model(x, c, gamma):
        return log2 * np.tanh(c * np.asarray(x) ** gamma) ** 2
    try:
        popt, _ = curve_fit(model, xs, ys, p0=[0.4, 2.2],
                            maxfev=10000, bounds=([1e-6, 0.5], [100, 10]))
        pred = model(xs, *popt)
        ss = float(np.sum((np.array(ys) - pred) ** 2))
        st = float(np.sum((np.array(ys) - np.mean(ys)) ** 2))
        return {'c': float(popt[0]), 'gamma': float(popt[1]),
                'r2': float(1 - ss / st) if st > 0 else 0.0}
    except Exception as e:
        return {'error': str(e)}


def main():
    print("=" * 72)
    print("EP-2h: bond universality + n_max convergence + bounded form")
    print("=" * 72)

    rows = []

    # Q2: single-center at n_max=2,3,4 (alpha drift test).
    print("\nQ2: single-center n_max scan, Z in {2,3,4,6,10}")
    for n_max in (2, 3, 4):
        for Z in (2.0, 3.0, 4.0, 6.0, 10.0):
            r = _single_center_row(Z, n_max)
            r['system'] = f'HeLike_Z{Z:g}_nmax{n_max}'
            r['kind'] = 'single_center'
            r['Z'] = Z
            rows.append(r)
            print(f"  Z={Z:4.1f} n_max={n_max}  "
                  f"S={r['S_B']:.3e}  "
                  f"w={r['w_B_dimless']:.3e}  "
                  f"delta={r['delta_B']:.3e}  "
                  f"w/d={r['w_B_dimless']/r['delta_B']:.3f}")

    # Q1: LiH bonds from EP-2g, plus BeH2 and H2O bonds.
    print("\nQ1: LiH bond R-sweep (max_n=2) + BeH2 + H2O + HF")
    for R in (2.0, 3.015, 4.0, 5.0):
        r = _bond_block_row(lih_spec, R=R, max_n=2)
        r['system'] = f'LiH_bond_R{R:g}_nmax2'; r['kind'] = 'bond_LiH'
        rows.append(r)
        print(f"  LiH  R={R:5.2f} nmax=2  S={r['S_B']:.3e}  "
              f"w={r['w_B_dimless']:.3e}  d={r['delta_B']:.3e}  "
              f"w/d={r['w_B_dimless']/r['delta_B']:.3f}")

    # BeH2: use default R, vary max_n=2
    try:
        r = _bond_block_row(beh2_spec, R=None, bond_block_index=0, max_n=2)
        r['system'] = 'BeH2_bond1_nmax2'; r['kind'] = 'bond_BeH2'
        rows.append(r)
        print(f"  BeH2 default         S={r['S_B']:.3e}  "
              f"w={r['w_B_dimless']:.3e}  d={r['delta_B']:.3e}  "
              f"w/d={r['w_B_dimless']/r['delta_B']:.3f}")
    except Exception as e:
        print(f"  BeH2: FAILED {e}")

    # H2O: default geometry
    try:
        r = _bond_block_row(h2o_spec, R=None, bond_block_index=0, max_n=2)
        r['system'] = 'H2O_bond1_nmax2'; r['kind'] = 'bond_H2O'
        rows.append(r)
        print(f"  H2O  default         S={r['S_B']:.3e}  "
              f"w={r['w_B_dimless']:.3e}  d={r['delta_B']:.3e}  "
              f"w/d={r['w_B_dimless']/r['delta_B']:.3f}")
    except Exception as e:
        print(f"  H2O: FAILED {e}")

    # HF: single bond
    try:
        r = _bond_block_row(hf_spec, R=None, bond_block_index=0, max_n=2)
        r['system'] = 'HF_bond_nmax2'; r['kind'] = 'bond_HF'
        rows.append(r)
        print(f"  HF   default         S={r['S_B']:.3e}  "
              f"w={r['w_B_dimless']:.3e}  d={r['delta_B']:.3e}  "
              f"w/d={r['w_B_dimless']/r['delta_B']:.3f}")
    except Exception as e:
        print(f"  HF: FAILED {e}")

    # --- Analysis ---
    print("\n" + "-" * 72)
    print("Q2 answer: single-center alpha(Z) by n_max (pure power-law fit)")
    for n_max in (2, 3, 4):
        sub = [r for r in rows if r['kind'] == 'single_center'
               and r['n_max'] == n_max]
        if len(sub) < 2:
            continue
        f = _fit_power([x['w_B_dimless'] / x['delta_B'] for x in sub],
                       [x['S_B'] for x in sub])
        print(f"  n_max={n_max}: gamma={f['gamma']:.4f}  "
              f"A={f['A']:.4f}  R^2={f['r2']:.4f}  n={f['n']}")

    print("\nQ1 answer: per-bond gamma (pure power-law fit)")
    for kind in ('bond_LiH', 'bond_BeH2', 'bond_H2O', 'bond_HF'):
        sub = [r for r in rows if r['kind'] == kind]
        if len(sub) < 2:
            # Single point: just report values
            for x in sub:
                print(f"  {kind:11s}: 1 pt  "
                      f"S={x['S_B']:.3e} w/d={x['w_B_dimless']/x['delta_B']:.3f}")
            continue
        f = _fit_power([x['w_B_dimless'] / x['delta_B'] for x in sub],
                       [x['S_B'] for x in sub])
        print(f"  {kind:11s}: gamma={f['gamma']:.4f}  "
              f"A={f['A']:.4f}  R^2={f['r2']:.4f}  n={f['n']}")

    # Power-law fit on all data
    all_r = rows
    f_all = _fit_power([x['w_B_dimless'] / x['delta_B'] for x in all_r],
                       [x['S_B'] for x in all_r])
    print(f"\nAll {len(all_r)} points, pure power law:")
    print(f"  gamma={f_all['gamma']:.4f}  A={f_all['A']:.4f}  "
          f"R^2={f_all['r2']:.4f}")

    # Q3: bounded form
    xs_all = np.array([x['w_B_dimless'] / x['delta_B'] for x in all_r])
    ys_all = np.array([x['S_B'] for x in all_r])
    f_bounded = _fit_bounded(xs_all, ys_all)
    print(f"\nAll {len(all_r)} points, bounded form "
          f"S = log(2) * tanh^2(c * (w/d)^gamma):")
    if 'error' in f_bounded:
        print(f"  FAILED: {f_bounded['error']}")
    else:
        print(f"  c = {f_bounded['c']:.4f}  "
              f"gamma = {f_bounded['gamma']:.4f}  "
              f"R^2 = {f_bounded['r2']:.6f}")
        # Residuals
        def model(x, c, g):
            return np.log(2.0) * np.tanh(c * x ** g) ** 2
        resid = ys_all - model(xs_all, f_bounded['c'], f_bounded['gamma'])
        max_resid = float(np.max(np.abs(resid)))
        rms_resid = float(np.sqrt(np.mean(resid ** 2)))
        print(f"  max abs residual = {max_resid:.4e}  "
              f"RMS = {rms_resid:.4e}")

    out = os.path.join(os.path.dirname(__file__), 'data',
                       'ep2h_universality_extension.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump({
            'rows': rows,
            'n_max_fits': {
                f'n_max_{n}': _fit_power(
                    [x['w_B_dimless'] / x['delta_B'] for x in rows
                     if x['kind'] == 'single_center' and x['n_max'] == n],
                    [x['S_B'] for x in rows
                     if x['kind'] == 'single_center' and x['n_max'] == n])
                for n in (2, 3, 4)
            },
            'combined_power_law': f_all,
            'combined_bounded': f_bounded,
        }, f, indent=2)
    print(f"\nSaved -> {out}")


if __name__ == '__main__':
    main()
