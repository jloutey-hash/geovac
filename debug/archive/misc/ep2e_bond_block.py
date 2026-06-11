"""
Track EP-2e: Does the LiH bond block (two-center H_1) fall on the
same dimensionless entropy power law as single-center blocks?

The bond block has orbitals on two centers (Li-side + H-side). Its
H_1 contains kinetic energy on each sub-block plus the full two-center
Coulomb attraction (Li-orbital feels H-nucleus, H-orbital feels
Li-nucleus — cross-center V_ne from the balanced builder). Its V_ee
is the two-body Coulomb projected onto bond orbitals only (within
the bond block).

We extract the bond sub-matrix from the balanced LiH Hamiltonian and
perform a 2-electron singlet FCI on the M_L=0 bond-orbital subspace,
then plot its (w_B_dim, S_B) point against the single-center curve from
EP-2f.

Output: debug/data/ep2e_bond_block.json
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import lih_spec


def _enumerate_states(max_n, l_min=0):
    """Enumerate (n, l, m) with n_min<=n<=max_n, l_min<=l<n, |m|<=l."""
    states = []
    for n in range(1, max_n + 1):
        for l in range(l_min, n):
            for m in range(-l, l + 1):
                states.append((n, l, m))
    return states


def _bond_orbital_metadata(spec, n_max=2):
    """Return orbital metadata for the LiH bond block (center + partner).

    The composed builder orders orbitals as:
        [core_center, bond_center, bond_partner]
    with the core block having no partner. Each sub-block enumerates
    (n, l, m) for n=1..max_n, l=0..n-1, m=-l..l.

    Returns
    -------
    list of dicts:
        {'global_idx': int, 'sub_block': str, 'nlm': (n,l,m), 'center': 'Li'/'H'}
    """
    blocks = spec.blocks
    orb_info = []
    global_idx = 0

    for blk in blocks:
        center_states = _enumerate_states(blk.max_n, l_min=blk.l_min)
        for nlm in center_states:
            orb_info.append({
                'global_idx': global_idx,
                'block_label': blk.label,
                'sub_block': blk.label + '_center',
                'nlm': nlm,
                'center': 'Li' if 'core' in blk.label or blk.Z_center > 1.5 else 'X',
                'parent_block': blk.label,
            })
            global_idx += 1

        if blk.has_h_partner:
            partner_n = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            partner_states = _enumerate_states(partner_n)
            for nlm in partner_states:
                orb_info.append({
                    'global_idx': global_idx,
                    'block_label': blk.label,
                    'sub_block': blk.label + '_partner',
                    'nlm': nlm,
                    'center': 'H',
                    'parent_block': blk.label,
                })
                global_idx += 1

    return orb_info


def _solve_2e_singlet_fci(h1_sub: np.ndarray, eri_sub: np.ndarray,
                          orb_m: list):
    """Build and diagonalize 2-electron singlet CI in an orbital subspace.

    Uses the chemistry convention ERI: eri[p,r,q,s] = (pr|qs) = <pq|rs>.
    Singlet spatial configs: (i <= j), spatial symmetric.
    Two-electron singlet matrix element:
        <(ij)|H|(kl)> = h1 terms + (1 + P_spatial)/N * <ij|V|kl>_direct
    For singlet spin-adapted configs this is standard (see
    debug/energy_entanglement_decoupling).
    """
    M = h1_sub.shape[0]
    configs = []
    for i in range(M):
        for j in range(i, M):
            if orb_m[i] + orb_m[j] == 0:
                configs.append((i, j))
    n_configs = len(configs)
    if n_configs == 0:
        raise RuntimeError('no M_L=0 configs')

    H = np.zeros((n_configs, n_configs))
    for I, (i, j) in enumerate(configs):
        for J in range(I, n_configs):
            p, q = configs[J]
            I_perms = [(i, j)] + ([(j, i)] if i != j else [])
            J_perms = [(p, q)] + ([(q, p)] if p != q else [])
            N_I = np.sqrt(float(len(I_perms)))
            N_J = np.sqrt(float(len(J_perms)))
            me = 0.0
            for a_i, b_i in I_perms:
                for c_i, d_i in J_perms:
                    # h1 contributions
                    if b_i == d_i:
                        me += h1_sub[a_i, c_i]
                    if a_i == c_i:
                        me += h1_sub[b_i, d_i]
                    # V_ee contribution: <a b | V | c d> = eri[a,c,b,d]
                    me += eri_sub[a_i, c_i, b_i, d_i]
            me /= (N_I * N_J)
            H[I, J] = me
            H[J, I] = me
    return H, configs


def _build_1rdm_singlet(ci_coeffs, configs, n_spatial):
    rho = np.zeros((n_spatial, n_spatial))
    for I, (i, j) in enumerate(configs):
        for J, (p, q) in enumerate(configs):
            coeff = ci_coeffs[I] * ci_coeffs[J]
            if abs(coeff) < 1e-16:
                continue
            N_I = np.sqrt(2.0) if i != j else 1.0
            N_J = np.sqrt(2.0) if p != q else 1.0
            I_perms = [(i, j)] + ([(j, i)] if i != j else [])
            J_perms = [(p, q)] + ([(q, p)] if p != q else [])
            for a, b in I_perms:
                for c, d in J_perms:
                    if b == d:
                        rho[a, c] += coeff / (N_I * N_J)
    rho *= 2.0
    return rho


def _von_neumann_entropy(rho):
    occ = np.linalg.eigvalsh(rho)
    occ = np.clip(occ, 0.0, None)
    total = occ.sum()
    if total <= 1e-14:
        return 0.0
    p = occ / total
    p = p[p > 1e-14]
    return float(-np.sum(p * np.log(p)))


def main():
    t0 = time.time()
    spec = lih_spec()
    print("Building LiH balanced Hamiltonian (no PK)...")
    result = build_balanced_hamiltonian(spec, R=3.015, verbose=True)
    h1 = result['h1']
    eri = result['eri']
    M = h1.shape[0]
    print(f"Full Hamiltonian: M = {M}")

    orb_info = _bond_orbital_metadata(spec)
    bond_idx = [o['global_idx'] for o in orb_info if 'bond' in o['parent_block']]
    bond_m = [o['nlm'][2] for o in orb_info if 'bond' in o['parent_block']]
    print(f"Bond-block orbitals: {len(bond_idx)}  indices {bond_idx}")
    for o in orb_info:
        if 'bond' in o['parent_block']:
            print(f"  idx={o['global_idx']:2d}  "
                  f"{o['sub_block']:20s}  nlm={o['nlm']}  center={o['center']}")

    # Extract sub-matrices.
    h1_sub = h1[np.ix_(bond_idx, bond_idx)]
    eri_sub = eri[np.ix_(bond_idx, bond_idx, bond_idx, bond_idx)]

    # Solve 2-electron singlet FCI in the M_L=0 bond subspace.
    H_full, configs = _solve_2e_singlet_fci(h1_sub, eri_sub, bond_m)

    # Kinetic-only restriction: zero V_ee (within-bond-block ERIs).
    H_kin_config, _ = _solve_2e_singlet_fci(h1_sub, np.zeros_like(eri_sub), bond_m)
    # V_ee in config-space = full - kin.
    V_ee_config = H_full - H_kin_config

    # Diagonalize both.
    eigs_full, vecs_full = np.linalg.eigh(H_full)
    eigs_kin, vecs_kin = np.linalg.eigh(H_kin_config)

    rho_full = _build_1rdm_singlet(vecs_full[:, 0], configs, len(bond_idx))
    rho_kin = _build_1rdm_singlet(vecs_kin[:, 0], configs, len(bond_idx))
    S_full = _von_neumann_entropy(rho_full)
    S_kin = _von_neumann_entropy(rho_kin)

    # Dimensionless off-diagonal mass w_B_dim.
    _e, U = np.linalg.eigh(H_kin_config)
    V_in = U.T @ V_ee_config @ U
    V_frob = np.linalg.norm(V_in)
    V_diag_frob = np.sqrt(np.sum(np.diag(V_in) ** 2))
    V_off = np.sqrt(max(0.0, V_frob ** 2 - V_diag_frob ** 2))
    kin_frob = np.linalg.norm(H_kin_config)
    w_B_dim = float(V_off / kin_frob) if kin_frob > 0 else 0.0

    # Reference-line prediction using EP-2f local fit.
    # Bond Z_eff mix: Li side Z_eff=1 + H side Z=1. Expected w_B_dim range
    # puts it in the low-Z regime.  Use EP-2c reference alpha=2.383,
    # A=8.16 AND the local slope at small Z for comparison.
    A_ref, alpha_ref = 8.16, 2.383
    S_pred_ref = A_ref * w_B_dim ** alpha_ref

    # Also use EP-2f Z<=4.5 band fit: alpha=2.588, A=14.02
    A_lowZ, alpha_lowZ = 14.0241, 2.588
    S_pred_lowZ = A_lowZ * w_B_dim ** alpha_lowZ

    print("\n" + "=" * 72)
    print(f"LiH bond block (2-electron FCI, two-center H_1):")
    print(f"  n_orbitals        = {len(bond_idx)}")
    print(f"  n_configs (M_L=0) = {len(configs)}")
    print(f"  E_full (GS)       = {eigs_full[0]:.6f} Ha")
    print(f"  E_kin  (GS)       = {eigs_kin[0]:.6f} Ha")
    print(f"  S_full            = {S_full:.4e} nats")
    print(f"  S_kin             = {S_kin:.4e} nats")
    print(f"  w_B_dim (dimless)     = {w_B_dim:.4e}")
    print(f"  ||H_kin||_F       = {kin_frob:.4f}")
    print(f"  V_diag_fraction   = {V_diag_frob / V_frob:.4f}")
    print(f"\nSingle-center reference predictions at this w_B_dim:")
    print(f"  S_pred (EP-2c He-like, alpha=2.383, A=8.16)  = {S_pred_ref:.4e}")
    print(f"  S_pred (EP-2f low-Z, alpha=2.588, A=14.02)   = {S_pred_lowZ:.4e}")
    ratio_ref = S_full / S_pred_ref if S_pred_ref > 0 else float('nan')
    ratio_lowZ = S_full / S_pred_lowZ if S_pred_lowZ > 0 else float('nan')
    print(f"\nRatio S_bond / S_pred:")
    print(f"  vs reference   = {ratio_ref:.3f}")
    print(f"  vs low-Z fit   = {ratio_lowZ:.3f}")

    if 0.5 < ratio_lowZ < 2.0:
        verdict = ('The bond block lies on the same single-center curve '
                   'within factor-of-2 — two-center H_1 does not produce '
                   'a qualitatively new regime.')
    else:
        verdict = ('The bond block is off the single-center curve by '
                   f'{ratio_lowZ:.2f}x — two-center H_1 produces distinct '
                   'behavior.')
    print(f"\nVerdict: {verdict}")

    summary = {
        'n_bond_orbitals': len(bond_idx),
        'n_bond_configs': len(configs),
        'bond_orbital_indices': bond_idx,
        'bond_orbital_info': [
            {'global_idx': o['global_idx'], 'sub_block': o['sub_block'],
             'nlm': list(o['nlm']), 'center': o['center']}
            for o in orb_info if 'bond' in o['parent_block']
        ],
        'E_full': float(eigs_full[0]),
        'E_kin': float(eigs_kin[0]),
        'S_full': S_full,
        'S_kin': S_kin,
        'w_B_dimless': w_B_dim,
        'V_frobenius': float(V_frob),
        'V_diag_fraction_unsquared': float(V_diag_frob / V_frob),
        'kin_frobenius': float(kin_frob),
        'reference_predictions': {
            'S_pred_he_like_ref': float(S_pred_ref),
            'S_pred_low_Z_band': float(S_pred_lowZ),
            'ratio_S_bond_over_S_pred_ref': float(ratio_ref),
            'ratio_S_bond_over_S_pred_low_Z': float(ratio_lowZ),
        },
        'wall_seconds': time.time() - t0,
    }
    out = os.path.join(os.path.dirname(__file__), 'data',
                       'ep2e_bond_block.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved -> {out}")


if __name__ == '__main__':
    main()
