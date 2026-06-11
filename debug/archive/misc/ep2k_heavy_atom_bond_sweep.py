"""
Track EP-2k: HF and H2O bond R-sweeps.

EP-2g established the universal curve on LiH bond R-sweep + single-center
Z-sweep. EP-2h added single-point checks for BeH2/H2O. This sprint adds
R-sweeps for HF and H2O — heavier atoms, smaller bonds, larger Z mismatch
between center and partner.

Question: do HF and H2O bonds follow the same gamma_loc as LiH at
matched (w/delta) ratio, or does heavy-atom Z_eff change the local slope?

Output: debug/data/ep2k_heavy_atom_bond_sweep.json
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from debug.entanglement_geometry import (
    build_1rdm_from_singlet_ci, compute_entanglement_measures,
)
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import hf_spec, h2o_spec, lih_spec


def _enumerate_states(max_n, l_min=0):
    return [(n, l, m) for n in range(1, max_n + 1)
            for l in range(l_min, n) for m in range(-l, l + 1)]


def _hf_nuclei(R):
    return [
        {'Z': 9.0, 'position': (0.0, 0.0, 0.0), 'label': 'F'},
        {'Z': 1.0, 'position': (0.0, 0.0, R), 'label': 'H'},
    ]


def _h2o_nuclei(R, angle_HOH_deg=104.5):
    half = np.radians(angle_HOH_deg / 2.0)
    return [
        {'Z': 8.0, 'position': (0.0, 0.0, 0.0), 'label': 'O'},
        {'Z': 1.0, 'position': (R * np.sin(half), 0.0, R * np.cos(half)),
         'label': 'H1'},
        {'Z': 1.0, 'position': (-R * np.sin(half), 0.0, R * np.cos(half)),
         'label': 'H2'},
    ]


def _bond_block_row(spec_fn, R, nuclei_fn, bond_block_index=0):
    spec = spec_fn(R=R)
    result = build_balanced_hamiltonian(spec, R=R, nuclei=nuclei_fn(R),
                                        verbose=False)
    h1 = result['h1']; eri = result['eri']
    blocks = spec.blocks

    bond_indices = [i for i, b in enumerate(blocks) if b.block_type == 'bond']
    target_idx = bond_indices[bond_block_index]
    target_label = blocks[target_idx].label

    orb_m, target_orb = [], []
    offset = 0
    for i_blk, blk in enumerate(blocks):
        for (n, l, m) in _enumerate_states(blk.max_n, l_min=blk.l_min):
            if i_blk == target_idx:
                target_orb.append(offset); orb_m.append(m)
            offset += 1
        if blk.has_h_partner:
            partner_n = blk.max_n_partner or blk.max_n
            for (n, l, m) in _enumerate_states(partner_n):
                if i_blk == target_idx:
                    target_orb.append(offset); orb_m.append(m)
                offset += 1

    h1_sub = h1[np.ix_(target_orb, target_orb)]
    eri_sub = eri[np.ix_(target_orb, target_orb, target_orb, target_orb)]
    M = len(target_orb)
    configs = [(i, j) for i in range(M) for j in range(i, M)
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
            H_kin[I, J] = H_kin[J, I] = me_kin
    V_ee = H_full - H_kin
    eigs, vecs = np.linalg.eigh(H_full)
    rho = build_1rdm_from_singlet_ci(vecs[:, 0], configs, M)
    S_B = float(compute_entanglement_measures(rho)['von_neumann_entropy'])
    e, U = np.linalg.eigh(H_kin)
    V_in = U.T @ V_ee @ U
    V_frob = np.linalg.norm(V_in)
    V_diag = np.sqrt(np.sum(np.diag(V_in) ** 2))
    V_off = float(np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2)))
    kin = float(np.linalg.norm(H_kin))
    e_sort = np.sort(e)
    return {
        'R': float(R),
        'block': target_label,
        'n_orb': M,
        'n_cfg': n_cfg,
        'E_full': float(eigs[0]),
        'S_B': S_B,
        'w_B_dimless': V_off / kin,
        'delta_B': float(abs(e_sort[1] - e_sort[0]) / kin),
        'kin_frobenius': kin,
    }


def _local_slope(rs):
    """Compute local log-log slope between consecutive R points."""
    rs = sorted(rs, key=lambda r: r['R'])
    slopes = []
    for i in range(len(rs) - 1):
        if rs[i]['S_B'] <= 0 or rs[i+1]['S_B'] <= 0:
            continue
        x1 = rs[i]['w_B_dimless'] / rs[i]['delta_B']
        x2 = rs[i+1]['w_B_dimless'] / rs[i+1]['delta_B']
        if x1 == x2: continue
        slope = (np.log(rs[i+1]['S_B']) - np.log(rs[i]['S_B'])) / \
                (np.log(x2) - np.log(x1))
        slopes.append((float(0.5 * (rs[i]['R'] + rs[i+1]['R'])), float(slope)))
    return slopes


def main():
    t0 = time.time()
    print("=" * 72)
    print("EP-2k: Heavy-atom bond R-sweeps (HF, H2O, LiH reference)")
    print("=" * 72)

    out = {}

    # LiH reference (already in EP-2g, recompute for self-consistency)
    print("\nLiH reference (R = 2..6 bohr):")
    lih_rows = []
    for R in (2.0, 3.015, 4.0, 5.0, 6.0):
        try:
            r = _bond_block_row(lih_spec, R,
                                lambda Rv: [
                                    {'Z': 3.0, 'position': (0,0,0), 'label': 'Li'},
                                    {'Z': 1.0, 'position': (0,0,Rv), 'label': 'H'},
                                ])
            lih_rows.append(r)
            ratio = r['w_B_dimless'] / r['delta_B'] if r['delta_B'] > 1e-10 else float('inf')
            print(f"  R={R:5.2f}  S={r['S_B']:.3e}  "
                  f"w={r['w_B_dimless']:.3e}  d={r['delta_B']:.3e}  w/d={ratio:.3f}")
        except Exception as e:
            print(f"  R={R}: FAILED {e}")
    out['LiH'] = {'rows': lih_rows, 'local_slopes': _local_slope(lih_rows)}

    # HF (R = 1.5..3 bohr; equilibrium 1.733)
    print("\nHF (R = 1.3..3.0 bohr; eq 1.733):")
    hf_rows = []
    for R in (1.3, 1.733, 2.0, 2.5, 3.0):
        try:
            r = _bond_block_row(hf_spec, R, _hf_nuclei)
            hf_rows.append(r)
            ratio = r['w_B_dimless'] / r['delta_B'] if r['delta_B'] > 1e-10 else float('inf')
            print(f"  R={R:5.2f}  S={r['S_B']:.3e}  "
                  f"w={r['w_B_dimless']:.3e}  d={r['delta_B']:.3e}  w/d={ratio:.3f}")
        except Exception as e:
            print(f"  R={R}: FAILED {e}")
    out['HF'] = {'rows': hf_rows, 'local_slopes': _local_slope(hf_rows)}

    # H2O (R = 1.5..3 bohr; equilibrium 1.809)
    print("\nH2O (R_OH = 1.5..3.0 bohr; eq 1.809):")
    h2o_rows = []
    for R in (1.5, 1.809, 2.2, 2.6, 3.0):
        try:
            r = _bond_block_row(h2o_spec, R, _h2o_nuclei)
            h2o_rows.append(r)
            ratio = r['w_B_dimless'] / r['delta_B'] if r['delta_B'] > 1e-10 else float('inf')
            print(f"  R={R:5.2f}  S={r['S_B']:.3e}  "
                  f"w={r['w_B_dimless']:.3e}  d={r['delta_B']:.3e}  w/d={ratio:.3f}")
        except Exception as e:
            print(f"  R={R}: FAILED {e}")
    out['H2O'] = {'rows': h2o_rows, 'local_slopes': _local_slope(h2o_rows)}

    # Summary: compare local slopes across molecules at matched w/d
    print("\n" + "-" * 72)
    print("Local log-log slope (S vs w/d) at consecutive R-points:")
    for mol in ('LiH', 'HF', 'H2O'):
        slopes = out[mol]['local_slopes']
        print(f"  {mol}:")
        for R_mid, sl in slopes:
            print(f"    R~{R_mid:5.2f}: gamma_loc = {sl:.4f}")

    # Universal-curve check: ratio S_actual / S_predicted using
    # A=0.157, gamma=2.228 from EP-2g.
    A_uni, g_uni = 0.157, 2.228
    print("\nUniversal-curve check (A=0.157, gamma=2.228):")
    for mol in ('LiH', 'HF', 'H2O'):
        for r in out[mol]['rows']:
            if r['delta_B'] < 1e-10:
                continue
            ratio = r['w_B_dimless'] / r['delta_B']
            S_pred = A_uni * ratio ** g_uni
            dev = r['S_B'] / S_pred if S_pred > 0 else float('inf')
            print(f"  {mol} R={r['R']:5.2f}: S/S_pred = {dev:.3f}")

    out['_meta'] = {
        'wall_seconds': time.time() - t0,
        'reference_universal': {'A': A_uni, 'gamma': g_uni},
    }
    out_path = os.path.join(os.path.dirname(__file__), 'data',
                            'ep2k_heavy_atom_bond_sweep.json')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved -> {out_path}")
    print(f"Wall: {time.time() - t0:.1f}s")


if __name__ == '__main__':
    main()
