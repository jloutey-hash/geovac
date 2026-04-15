"""
Track EP-2i: follow-on diagnostics for EP-2h.

Q1 (BeH2 saturation): is BeH2 bond saturation at equilibrium driven
    by Z_eff asymmetry (Be_eff=2 vs H=1, bigger gap than LiH's 1:1)?
    Diagnose the bond-block H_kin eigenvalue spectrum.

Q2 (H2O factor-2): is the H2O bond-block offset from the single-center
    curve due to cross-block ERIs (lone-pair influence via the
    balanced builder) or intrinsic to the bond block itself?
    Compare balanced H2O bond (with cross-block) vs composed bond
    (within-block only).

Q3 (gamma(Z, n_max) surface): map gamma(Z, n_max) on a 5x3 grid
    and check whether the two asymptotic limits (Z->inf in EP-2f
    gives alpha=1.92; n_max->inf in EP-2h gives alpha>2.3) reconcile
    on a 2D crossover surface.

Output: debug/data/ep2i_deep_diagnostics.json
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from debug.energy_entanglement_decoupling import (
    build_decomposed_hamiltonians, solve_and_entangle,
)
from debug.entanglement_geometry import (
    build_1rdm_from_singlet_ci, compute_entanglement_measures,
)
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.composed_qubit import build_composed_hamiltonian
from geovac.molecular_spec import lih_spec, beh2_spec, h2o_spec


def _enumerate_states(max_n, l_min=0):
    return [(n, l, m) for n in range(1, max_n + 1)
            for l in range(l_min, n) for m in range(-l, l + 1)]


def _extract_bond_subblock(result, spec, bond_block_index=0):
    """Extract h1/eri sub-matrix, m-labels, and target label for a bond block."""
    h1 = result['h1']; eri = result['eri']
    blocks = spec.blocks
    bond_indices = [i for i, blk in enumerate(blocks)
                    if blk.block_type == 'bond']
    target_idx = bond_indices[bond_block_index]
    target_label = blocks[target_idx].label

    orb_m = []
    offset = 0
    target_orb_idx = []
    for i_blk, blk in enumerate(blocks):
        for (n, l, m) in _enumerate_states(blk.max_n, l_min=blk.l_min):
            if i_blk == target_idx:
                target_orb_idx.append(offset); orb_m.append(m)
            offset += 1
        if blk.has_h_partner:
            partner_n = blk.max_n_partner or blk.max_n
            for (n, l, m) in _enumerate_states(partner_n):
                if i_blk == target_idx:
                    target_orb_idx.append(offset); orb_m.append(m)
                offset += 1

    h1_sub = h1[np.ix_(target_orb_idx, target_orb_idx)]
    eri_sub = eri[np.ix_(target_orb_idx, target_orb_idx,
                         target_orb_idx, target_orb_idx)]
    return h1_sub, eri_sub, orb_m, target_label


def _solve_bond_block(h1_sub, eri_sub, orb_m):
    """2e singlet M_L=0 FCI on a sub-block."""
    M = h1_sub.shape[0]
    configs = [(i, j) for i in range(M) for j in range(i, M)
               if orb_m[i] + orb_m[j] == 0]
    n = len(configs)
    H_full = np.zeros((n, n)); H_kin = np.zeros((n, n))
    for I, (i, j) in enumerate(configs):
        for J in range(I, n):
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
    e_kin, U = np.linalg.eigh(H_kin)
    V_in = U.T @ V_ee @ U
    V_frob = np.linalg.norm(V_in)
    V_diag = np.sqrt(np.sum(np.diag(V_in) ** 2))
    V_off = float(np.sqrt(max(0.0, V_frob ** 2 - V_diag ** 2)))
    kin = float(np.linalg.norm(H_kin))
    e_sort = np.sort(e_kin)
    return {
        'S_B': S_B,
        'w_B_dimless': V_off / kin,
        'delta_B': float(abs(e_sort[1] - e_sort[0]) / kin),
        'kin_spectrum_first6': [float(x) for x in e_sort[:6]],
        'kin_gaps_first5': [float(e_sort[i+1] - e_sort[i]) for i in range(5)],
        'n_orb': M, 'n_cfg': n,
    }


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
    kin = float(np.linalg.norm(H_kin))
    e_sort = np.sort(e)
    return {
        'S_B': float(ent['von_neumann_entropy']),
        'w_B_dimless': V_off / kin,
        'delta_B': float(abs(e_sort[1] - e_sort[0]) / kin),
        'Z': Z, 'n_max': n_max,
    }


def _fit(ratios, Ss):
    lx = np.log(np.asarray(ratios, float))
    ly = np.log(np.asarray(Ss, float))
    a, b = np.polyfit(lx, ly, 1)
    pred = a * lx + b
    ss = float(np.sum((ly - pred) ** 2))
    st = float(np.sum((ly - ly.mean()) ** 2))
    return {'gamma': float(a), 'A': float(np.exp(b)),
            'r2': float(1 - ss / st) if st > 0 else 0.0}


def main():
    print("=" * 72)
    print("EP-2i: BeH2 diagnosis + H2O cross-block + gamma(Z, n_max) surface")
    print("=" * 72)

    out = {'Q1_BeH2': {}, 'Q2_H2O': {}, 'Q3_gamma_surface': {}}

    # --- Q1: BeH2 saturation diagnosis ---
    print("\nQ1: LiH vs BeH2 bond-block diagnostics (default geometry)")
    for mol_name, spec_fn in (('LiH', lih_spec), ('BeH2', beh2_spec)):
        spec = spec_fn()
        # Z_eff on bond block
        bond_blks = [b for b in spec.blocks if b.block_type == 'bond']
        Z_center = bond_blks[0].Z_center; Z_partner = bond_blks[0].Z_partner
        # Balanced (with cross-block ERIs)
        r_bal = build_balanced_hamiltonian(spec, R=3.015 if mol_name == 'LiH' else 2.54,
                                           verbose=False)
        h1_sub, eri_sub, m_list, lbl = _extract_bond_subblock(r_bal, spec)
        diag_bal = _solve_bond_block(h1_sub, eri_sub, m_list)

        # Composed (no cross-block ERIs)
        r_comp = build_composed_hamiltonian(spec, pk_in_hamiltonian=False,
                                            verbose=False)
        h1_sub_c, eri_sub_c, _, _ = _extract_bond_subblock(r_comp, spec)
        diag_comp = _solve_bond_block(h1_sub_c, eri_sub_c, m_list)

        print(f"  {mol_name} bond block:  Z_center={Z_center}  Z_partner={Z_partner}")
        print(f"    balanced:  S={diag_bal['S_B']:.3e}  "
              f"w={diag_bal['w_B_dimless']:.3e}  d={diag_bal['delta_B']:.3e}  "
              f"w/d={diag_bal['w_B_dimless']/diag_bal['delta_B']:.3f}")
        wd_str = (f"{diag_comp['w_B_dimless']/diag_comp['delta_B']:.3f}"
                  if diag_comp['delta_B'] > 1e-10 else 'DEGENERATE')
        print(f"    composed:  S={diag_comp['S_B']:.3e}  "
              f"w={diag_comp['w_B_dimless']:.3e}  d={diag_comp['delta_B']:.3e}  "
              f"w/d={wd_str}")
        print(f"    H_kin first-5 gaps: "
              f"{[f'{g:.3f}' for g in diag_bal['kin_gaps_first5']]}")

        out['Q1_BeH2'][mol_name] = {
            'Z_center': Z_center, 'Z_partner': Z_partner,
            'balanced': diag_bal, 'composed': diag_comp,
        }

    # --- Q2: H2O bond cross-block vs within-block ---
    print("\nQ2: H2O bond block with and without cross-block ERIs")
    spec = h2o_spec()
    r_bal = build_balanced_hamiltonian(spec, R=3.015, verbose=False)
    r_comp = build_composed_hamiltonian(spec, pk_in_hamiltonian=False,
                                        verbose=False)
    h1_sub_b, eri_sub_b, m_list, _ = _extract_bond_subblock(r_bal, spec)
    h1_sub_c, eri_sub_c, _, _ = _extract_bond_subblock(r_comp, spec)
    diag_bal = _solve_bond_block(h1_sub_b, eri_sub_b, m_list)
    diag_comp = _solve_bond_block(h1_sub_c, eri_sub_c, m_list)
    print(f"  balanced (w/ cross-block):  S={diag_bal['S_B']:.3e}  "
          f"w/d={diag_bal['w_B_dimless']/diag_bal['delta_B']:.3f}")
    wd_str2 = (f"{diag_comp['w_B_dimless']/diag_comp['delta_B']:.3f}"
               if diag_comp['delta_B'] > 1e-10 else 'DEGENERATE')
    print(f"  composed (within-block):    S={diag_comp['S_B']:.3e}  "
          f"w/d={wd_str2}")
    S_ratio = diag_bal['S_B'] / diag_comp['S_B'] if diag_comp['S_B'] > 0 else float('nan')
    print(f"  S_balanced / S_composed = {S_ratio:.3f}")
    out['Q2_H2O'] = {
        'balanced': diag_bal, 'composed': diag_comp,
        'S_ratio': float(S_ratio),
    }

    # --- Q3: gamma(Z, n_max) surface ---
    print("\nQ3: gamma(Z, n_max) surface on Z in {2,3,4,6,10,15,30}, n_max in {2,3,4}")
    Zs = [2.0, 3.0, 4.0, 6.0, 10.0, 15.0, 30.0]
    rows = []
    for n_max in (2, 3, 4):
        for Z in Zs:
            r = _single_center_row(Z, n_max)
            rows.append(r)
    # Fit gamma per n_max across full Z range
    gamma_by_nmax = {}
    for n_max in (2, 3, 4):
        sub = [r for r in rows if r['n_max'] == n_max]
        f = _fit([r['w_B_dimless'] / r['delta_B'] for r in sub],
                 [r['S_B'] for r in sub])
        gamma_by_nmax[n_max] = f
        print(f"  n_max={n_max}: gamma={f['gamma']:.4f}  "
              f"A={f['A']:.4f}  R^2={f['r2']:.4f}  n={len(sub)}")

    # Fit gamma per Z across n_max (less well-conditioned but informative)
    gamma_by_Z = {}
    for Z in Zs:
        sub = [r for r in rows if abs(r['Z'] - Z) < 1e-6]
        if len(sub) < 2:
            continue
        f = _fit([r['w_B_dimless'] / r['delta_B'] for r in sub],
                 [r['S_B'] for r in sub])
        gamma_by_Z[Z] = f
        print(f"  Z={Z:5.1f}:  gamma={f['gamma']:.4f}  "
              f"A={f['A']:.4f}  R^2={f['r2']:.4f}  n={len(sub)}")

    # Local slopes at large Z, large n_max
    print("\n  Large-(Z, n_max) local slope:")
    for n_max in (2, 3, 4):
        sub = sorted([r for r in rows if r['n_max'] == n_max],
                     key=lambda x: x['Z'])
        if len(sub) >= 2:
            # Last two points
            ly = [np.log(r['S_B']) for r in sub[-2:]]
            lx = [np.log(r['w_B_dimless'] / r['delta_B']) for r in sub[-2:]]
            slope = (ly[1] - ly[0]) / (lx[1] - lx[0])
            print(f"    n_max={n_max}, Z={sub[-2]['Z']}->{sub[-1]['Z']}: "
                  f"local slope = {slope:.4f}")

    out['Q3_gamma_surface'] = {
        'grid_rows': rows,
        'gamma_by_nmax': gamma_by_nmax,
        'gamma_by_Z': gamma_by_Z,
    }

    # --- Verdict summaries ---
    print("\n" + "=" * 72)
    print("Verdict summaries:")
    print("=" * 72)

    # Q1
    lih_kin = out['Q1_BeH2']['LiH']['balanced']['kin_gaps_first5'][0]
    beh2_kin = out['Q1_BeH2']['BeH2']['balanced']['kin_gaps_first5'][0]
    print(f"  Q1 (BeH2 saturation):")
    print(f"     LiH  first H_kin gap (abs): {lih_kin:.4f}")
    print(f"     BeH2 first H_kin gap (abs): {beh2_kin:.4f}")
    if beh2_kin < lih_kin:
        print(f"     -> BeH2 has SMALLER H_kin gap ({beh2_kin/lih_kin:.2f}x), "
              "consistent with earlier saturation at equilibrium.")
    else:
        print(f"     -> BeH2 has LARGER H_kin gap, saturation driven by something else.")

    # Q2
    print(f"  Q2 (H2O offset): "
          f"S_balanced/S_composed = {out['Q2_H2O']['S_ratio']:.3f}")
    if abs(out['Q2_H2O']['S_ratio'] - 1) < 0.2:
        print(f"     -> Cross-block ERIs don't move S much; offset is intrinsic.")
    else:
        print(f"     -> Cross-block ERIs drive the offset significantly.")

    # Q3
    print(f"  Q3 (gamma(Z, n_max) surface): see per-Z and per-n_max fits above.")

    out_path = os.path.join(os.path.dirname(__file__), 'data',
                            'ep2i_deep_diagnostics.json')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved -> {out_path}")


if __name__ == '__main__':
    main()
