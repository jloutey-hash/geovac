"""Species-II discriminator: does which-site entanglement tell a BINDER from a
NON-BINDER, and is the 'fission aperture stuck open' reading real?

Control logic (the point of this script):
  LiH  binds in the balanced framework (known success, Paper 17/19).
  NaH/KH do NOT bind (the W1c-W1e wall, monotonic descent).
If the species-II reading is right ("the wall = fission aperture jammed at
MAXIMAL site-entanglement"), then:
  - NaH/KH: S_mode / S_mode_ceiling ~ 1 (pinned at ceiling) at ALL R.
  - LiH:    S_mode / S_mode_ceiling < 1 somewhere (aperture partially closed
            = a real bond), distinguishing it from the non-binders.
If LiH ALSO sits pinned at the ceiling, the entanglement-pinning is NOT the
bind/no-bind discriminator and the 'stuck aperture' reading is WEAKENED.
Either outcome is informative; we LOOK and report.

Two observables per (molecule, R), both degeneracy-ROBUST (ensemble reduced
density matrix averaged over any degenerate ground manifold, so no arbitrary-
eigenvector artifact — the trap the deleted fabrication fell into):
  S_mode : von Neumann entropy of the heavy|H mode-reduced density matrix.
  S_occ  : Shannon entropy of the electrons-per-site distribution (classical
           which-site; ceiling ln(N_sites)).
Also reported: ground-state degeneracy d (states within 1e-6 Ha), Schmidt rank
(effective dim of the heavy|H cut), and ln(rank) as the mode-entropy ceiling.

Architecture: screened_cross_center + cross_block_h1, multi_zeta OFF (Na-only),
so ALL molecules use the same architecture. FCI is particle-number-projected;
ground energy validated bit-exact against library coupled_fci_energy.
"""

from __future__ import annotations

import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import (coupled_fci_energy, _excitation_phase,
                                        _double_excitation_phase)
from geovac.molecular_spec import nah_spec, kh_spec, lih_spec


def build_fci(result, n_electrons, k=8):
    """Particle-number-projected FCI. Returns (eigvals, eigvecs, a_str, b_str, nb, M).

    Mirrors coupled_composition.coupled_fci_energy term-for-term.
    """
    M = result['M']; h1 = result['h1']; eri = result['eri']
    nuc = result['nuclear_repulsion']
    n_up = n_electrons // 2; n_down = n_electrons - n_up
    a_str = list(combinations(range(M), n_up))
    b_str = list(combinations(range(M), n_down))
    na, nb = len(a_str), len(b_str)
    a_idx = {s: i for i, s in enumerate(a_str)}
    b_idx = {s: i for i, s in enumerate(b_str)}
    n_det = na * nb
    H = lil_matrix((n_det, n_det))

    def di(ai, bi):
        return ai * nb + bi

    for ai, al in enumerate(a_str):
        for bi, be in enumerate(b_str):
            E = nuc
            for p in al:
                E += h1[p, p]
            for p in be:
                E += h1[p, p]
            for i in range(n_up):
                for j in range(i + 1, n_up):
                    p, q = al[i], al[j]; E += eri[p, p, q, q] - eri[p, q, q, p]
            for i in range(n_down):
                for j in range(i + 1, n_down):
                    p, q = be[i], be[j]; E += eri[p, p, q, q] - eri[p, q, q, p]
            for p in al:
                for q in be:
                    E += eri[p, p, q, q]
            H[di(ai, bi), di(ai, bi)] = E

    for ai, al in enumerate(a_str):
        aset = set(al)
        for p in al:
            for r in range(M):
                if r in aset:
                    continue
                new = tuple(sorted((aset - {p}) | {r}))
                if new not in a_idx:
                    continue
                ain = a_idx[new]; ph = _excitation_phase(al, p, r)
                for bi, be in enumerate(b_str):
                    val = ph * h1[r, p]
                    for q in al:
                        if q == p:
                            continue
                        val += ph * (eri[r, p, q, q] - eri[r, q, q, p])
                    for q in be:
                        val += ph * eri[r, p, q, q]
                    if abs(val) > 1e-14:
                        H[di(ai, bi), di(ain, bi)] += val

    for bi, be in enumerate(b_str):
        bset = set(be)
        for p in be:
            for r in range(M):
                if r in bset:
                    continue
                new = tuple(sorted((bset - {p}) | {r}))
                if new not in b_idx:
                    continue
                bin_ = b_idx[new]; ph = _excitation_phase(be, p, r)
                for ai, al in enumerate(a_str):
                    val = ph * h1[r, p]
                    for q in be:
                        if q == p:
                            continue
                        val += ph * (eri[r, p, q, q] - eri[r, q, q, p])
                    for q in al:
                        val += ph * eri[r, p, q, q]
                    if abs(val) > 1e-14:
                        H[di(ai, bi), di(ai, bin_)] += val

    for ai, al in enumerate(a_str):
        aset = set(al); occ = list(al)
        for i1 in range(n_up):
            for i2 in range(i1 + 1, n_up):
                p, q = occ[i1], occ[i2]
                for r in range(M):
                    if r in aset:
                        continue
                    for s in range(r + 1, M):
                        if s in aset:
                            continue
                        new = tuple(sorted((aset - {p, q}) | {r, s}))
                        if new not in a_idx:
                            continue
                        ain = a_idx[new]
                        ph = _double_excitation_phase(al, p, q, r, s)
                        val = ph * (eri[r, p, s, q] - eri[r, q, s, p])
                        if abs(val) > 1e-14:
                            for bi in range(nb):
                                H[di(ai, bi), di(ain, bi)] += val

    for bi, be in enumerate(b_str):
        bset = set(be); occ = list(be)
        for i1 in range(n_down):
            for i2 in range(i1 + 1, n_down):
                p, q = occ[i1], occ[i2]
                for r in range(M):
                    if r in bset:
                        continue
                    for s in range(r + 1, M):
                        if s in bset:
                            continue
                        new = tuple(sorted((bset - {p, q}) | {r, s}))
                        if new not in b_idx:
                            continue
                        bin_ = b_idx[new]
                        ph = _double_excitation_phase(be, p, q, r, s)
                        val = ph * (eri[r, p, s, q] - eri[r, q, s, p])
                        if abs(val) > 1e-14:
                            for ai in range(na):
                                H[di(ai, bi), di(ai, bin_)] += val

    for ai, al in enumerate(a_str):
        aset = set(al)
        for pa in al:
            for ra in range(M):
                if ra in aset:
                    continue
                newa = tuple(sorted((aset - {pa}) | {ra}))
                if newa not in a_idx:
                    continue
                ain = a_idx[newa]; pha = _excitation_phase(al, pa, ra)
                for bi, be in enumerate(b_str):
                    bset = set(be)
                    for pb in be:
                        for rb in range(M):
                            if rb in bset:
                                continue
                            newb = tuple(sorted((bset - {pb}) | {rb}))
                            if newb not in b_idx:
                                continue
                            bin_ = b_idx[newb]; phb = _excitation_phase(be, pb, rb)
                            val = pha * phb * eri[ra, pa, rb, pb]
                            if abs(val) > 1e-14:
                                H[di(ai, bi), di(ain, bin_)] += val

    H = H.tocsr(); H = 0.5 * (H + H.T)
    kk = min(k, n_det - 1)
    vals, vecs = eigsh(H, k=kk, which='SA')
    order = np.argsort(vals)
    return vals[order], vecs[:, order], a_str, b_str, nb, M


def site_of_spatial(result):
    site = []
    for b in result['blocks']:
        is_H = b['label'].endswith('_partner')
        site.extend([1 if is_H else 0] * b['n_orbitals'])
    return np.array(site, dtype=int)


def reshape_psi(vec, a_str, b_str, nb, site_spatial, M):
    """Return dict (keyA,keyB)->signed amp for one eigenvector (heavy|H cut)."""
    def sg(g):
        return site_spatial[g if g < M else g - M]
    amps = {}
    for ai, al in enumerate(a_str):
        for bi, be in enumerate(b_str):
            amp = vec[ai * nb + bi]
            if abs(amp) < 1e-13:
                continue
            occ = [p for p in al] + [M + q for q in be]
            seenA = 0; inv_before = 0; A = []; B = []
            for g in occ:
                if sg(g) == 0:
                    seenA += 1; A.append(g)
                else:
                    inv_before += seenA; B.append(g)
            inv = len(A) * len(B) - inv_before
            sign = -1.0 if (inv % 2) else 1.0
            amps[(tuple(sorted(A)), tuple(sorted(B)))] = \
                amps.get((tuple(sorted(A)), tuple(sorted(B))), 0.0) + sign * amp
    return amps


def ensemble_mode_entropy(vecs_manifold, a_str, b_str, nb, site_spatial, M):
    """von Neumann entropy of rho_A = (1/d) sum_i Tr_B|psi_i><psi_i|.

    Degeneracy-robust: basis-independent over the degenerate manifold.
    Returns (entropy_nats, schmidt_rank).
    """
    d = len(vecs_manifold)
    psis = [reshape_psi(v, a_str, b_str, nb, site_spatial, M) for v in vecs_manifold]
    keysA = sorted({k[0] for ps in psis for k in ps})
    keysB = sorted({k[1] for ps in psis for k in ps})
    iA = {k: i for i, k in enumerate(keysA)}
    iB = {k: i for i, k in enumerate(keysB)}
    rhoA = np.zeros((len(keysA), len(keysA)))
    for ps in psis:
        mat = np.zeros((len(keysA), len(keysB)))
        for (kA, kB), v in ps.items():
            mat[iA[kA], iB[kB]] = v
        rhoA += mat @ mat.T
    rhoA /= d
    ev = np.linalg.eigvalsh(rhoA)
    ev = ev[ev > 1e-15]
    ev = ev / ev.sum()
    S = float(-np.sum(ev * np.log(ev)))
    return S, int(len(ev))


def ensemble_occ_entropy(vecs_manifold, a_str, b_str, nb, site_spatial, M, n_sites=2):
    """Shannon entropy of electrons-per-site distribution, ensemble-averaged."""
    def sg(g):
        return site_spatial[g if g < M else g - M]
    d = len(vecs_manifold)
    dist = {}
    for v in vecs_manifold:
        for ai, al in enumerate(a_str):
            for bi, be in enumerate(b_str):
                w = abs(v[ai * nb + bi]) ** 2
                if w < 1e-14:
                    continue
                c = [0] * n_sites
                for p in al:
                    c[sg(p)] += 1
                for q in be:
                    c[sg(M + q)] += 1
                key = tuple(c)
                dist[key] = dist.get(key, 0.0) + w / d
    p = np.array(list(dist.values()))
    p = p[p > 1e-15]; p = p / p.sum()
    return float(-np.sum(p * np.log(p))), {str(k): float(w) for k, w in dist.items()}


def measure(spec_fn, name, Z, n_electrons, R_grid):
    spec = spec_fn()
    rows = []
    print(f"\n{name} (Z={Z}, {n_electrons}e)", flush=True)
    for j, R in enumerate(R_grid):
        result = build_balanced_hamiltonian(
            spec, R=R, screened_cross_center=True,
            multi_zeta_basis=False, cross_block_h1=True)
        M = result['M']; site_spatial = site_of_spatial(result)
        vals, vecs, a_str, b_str, nb, MM = build_fci(result, n_electrons, k=8)
        # degeneracy: states within 1e-6 Ha of ground
        gap = float(vals[1] - vals[0]) if len(vals) > 1 else float('inf')
        d = int(np.sum(vals - vals[0] < 1e-6))
        manifold = [vecs[:, i] for i in range(d)]
        S_mode, rank = ensemble_mode_entropy(manifold, a_str, b_str, nb, site_spatial, M)
        S_occ, dist = ensemble_occ_entropy(manifold, a_str, b_str, nb, site_spatial, M)
        S_mode_ceiling = math.log(rank) if rank > 1 else 0.0
        ratio = S_mode / S_mode_ceiling if S_mode_ceiling > 0 else 0.0
        rows.append({'R_bohr': float(R), 'E_Ha': float(vals[0]), 'gap_Ha': gap,
                     'degeneracy': d, 'S_mode': S_mode, 'schmidt_rank': rank,
                     'S_mode_ceiling_lnrank': S_mode_ceiling,
                     'S_mode_over_ceiling': ratio, 'S_occ': S_occ,
                     'occ_dist': dist})
        print(f"  R={R:5.2f}  E={vals[0]:11.4f}  gap={gap:8.5f}  d={d}  "
              f"S_mode={S_mode:.4f}/lnrank({rank})={S_mode_ceiling:.4f}"
              f"={ratio:.3f}  S_occ={S_occ:.4f}", flush=True)
    return {'molecule': name, 'Z': Z, 'n_electrons': n_electrons, 'rows': rows}


def validate(spec_fn, n_el):
    spec = spec_fn()
    res = build_balanced_hamiltonian(spec, R=3.566, screened_cross_center=True,
                                     multi_zeta_basis=False, cross_block_h1=True)
    ref = coupled_fci_energy(res, n_el)['E_coupled']
    vals, *_ = build_fci(res, n_el, k=4)
    return abs(vals[0] - ref)


def main():
    print("=" * 72)
    print("Species-II discriminator: binder (LiH) vs non-binders (NaH, KH)")
    print("Q: does S_mode/ceiling distinguish LiH (binds) from NaH/KH (wall)?")
    print("=" * 72)
    dv = validate(nah_spec, 2)
    dv2 = validate(lih_spec, 4)
    print(f"FCI validation vs library: NaH diff={dv:.2e}  LiH diff={dv2:.2e}  "
          f"{'OK' if max(dv, dv2) < 1e-6 else 'MISMATCH'}")

    R_grid = [2.0, 2.5, 3.0, 3.566, 4.0, 5.0, 7.0, 10.0]
    out = {}
    out['NaH'] = measure(nah_spec, 'NaH', 11, 2, R_grid)
    out['KH'] = measure(kh_spec, 'KH', 19, 2, R_grid)
    out['LiH'] = measure(lih_spec, 'LiH', 3, 4, R_grid)

    print("\n" + "=" * 72)
    print("DISCRIMINATOR: S_mode / ceiling (=1 means pinned at maximal)")
    print(f"{'R':>6} {'NaH':>8} {'KH':>8} {'LiH':>8}   "
          f"{'NaH_Socc':>9} {'KH_Socc':>9} {'LiH_Socc':>9}")
    for i, R in enumerate(R_grid):
        rn = out['NaH']['rows'][i]; rk = out['KH']['rows'][i]; rl = out['LiH']['rows'][i]
        print(f"{R:>6.2f} {rn['S_mode_over_ceiling']:>8.3f} "
              f"{rk['S_mode_over_ceiling']:>8.3f} {rl['S_mode_over_ceiling']:>8.3f}   "
              f"{rn['S_occ']:>9.4f} {rk['S_occ']:>9.4f} {rl['S_occ']:>9.4f}")
    print(f"\nln2={math.log(2):.4f} (S_occ ceiling for 2 sites)")
    print("If LiH ratio dips below NaH/KH near R_eq -> entanglement discriminates")
    print("bind/no-bind, species-II 'stuck aperture' reading SUPPORTED.")
    print("If all three pinned at ratio~1 -> ceiling trivial, reading WEAKENED.")

    out['_meta'] = {'R_grid': R_grid, 'ln2': math.log(2),
                    'fci_validation_NaH': float(dv), 'fci_validation_LiH': float(dv2),
                    'architecture': 'screened+cross_block_h1, multi_zeta=False'}
    Path('debug/data').mkdir(parents=True, exist_ok=True)
    Path('debug/data/species_ii_discriminator.json').write_text(json.dumps(out, indent=2))
    print("\nWrote debug/data/species_ii_discriminator.json")


if __name__ == '__main__':
    main()
