"""Species-II aperture signature — which-site entanglement curve S_site(R).

Charter (debug/confinement_reframing_charter.md) §7c/§7d. The fission-aperture
reading says molecular dissociation is a SPECIES-II aperture: which-site
entanglement should SATURATE at ln(N_sites) as R opens, NOT diverge. The W1e
chemistry wall (bonding pair collapsing without limit) is read as spectral
machinery (diverges) wrongly applied to a fission aperture (should saturate).

This script does NOT pre-commit to a scalar strain metric. It MEASURES THE
CURVE S_site(R) and asks two qualitative, pre-committed questions:
  (1) Does S_site(R) saturate at ln(N_sites) as R opens?  (fission signature)
  (2) Is the curve Z-invariant?  NaH (Z=11) vs KH (Z=19) — both 2 sites, both
      2 electrons, isostructural (identical Pauli=239, Q=20). Topology-locked
      => Z-invariant (NaH ~ KH). Charge-locked => NaH != KH => reframe WRONG.

The curve is the data; any scalar is chosen afterward as description.

Method (feedback_tc_correction.md): particle-number-projected FCI ONLY. We
build the FCI in the (n_up, n_down) determinant sector mirroring
geovac.coupled_composition.coupled_fci_energy EXACTLY (same h1/eri/phase
conventions), but return the ground EIGENVECTOR + strings. Site entanglement
is the fermionic mode-entanglement across the heavy|H spatial cut, computed in
the determinant basis with proper creation-operator reordering signs
(alpha-block ascending, then beta-block ascending — the builder's convention).
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
from geovac.coupled_composition import _excitation_phase, _double_excitation_phase
from geovac.molecular_spec import nah_spec, kh_spec


def fci_ground(result, n_electrons):
    """Particle-number-projected FCI ground state (faithful to coupled_fci_energy).

    Returns (E_gs, eigvec, alpha_strings, beta_strings, n_beta).
    """
    M = result['M']
    h1 = result['h1']
    eri = result['eri']
    nuc = result['nuclear_repulsion']
    n_up = n_electrons // 2
    n_down = n_electrons - n_up

    alpha_strings = list(combinations(range(M), n_up))
    beta_strings = list(combinations(range(M), n_down))
    na, nb = len(alpha_strings), len(beta_strings)
    n_det = na * nb
    a_idx = {s: i for i, s in enumerate(alpha_strings)}
    b_idx = {s: i for i, s in enumerate(beta_strings)}
    H = lil_matrix((n_det, n_det))

    def di(ai, bi):
        return ai * nb + bi

    # diagonal
    for ai, alpha in enumerate(alpha_strings):
        for bi, beta in enumerate(beta_strings):
            E = nuc
            for p in alpha:
                E += h1[p, p]
            for p in beta:
                E += h1[p, p]
            for i in range(n_up):
                for j in range(i + 1, n_up):
                    p, q = alpha[i], alpha[j]
                    E += eri[p, p, q, q] - eri[p, q, q, p]
            for i in range(n_down):
                for j in range(i + 1, n_down):
                    p, q = beta[i], beta[j]
                    E += eri[p, p, q, q] - eri[p, q, q, p]
            for p in alpha:
                for q in beta:
                    E += eri[p, p, q, q]
            H[di(ai, bi), di(ai, bi)] = E

    # single excitations within alpha
    for ai, alpha in enumerate(alpha_strings):
        aset = set(alpha)
        for p in alpha:
            for r in range(M):
                if r in aset:
                    continue
                new = tuple(sorted((aset - {p}) | {r}))
                if new not in a_idx:
                    continue
                ain = a_idx[new]
                ph = _excitation_phase(alpha, p, r)
                for bi, beta in enumerate(beta_strings):
                    val = ph * h1[r, p]
                    for q in alpha:
                        if q == p:
                            continue
                        val += ph * (eri[r, p, q, q] - eri[r, q, q, p])
                    for q in beta:
                        val += ph * eri[r, p, q, q]
                    if abs(val) > 1e-14:
                        H[di(ai, bi), di(ain, bi)] += val

    # single excitations within beta
    for bi, beta in enumerate(beta_strings):
        bset = set(beta)
        for p in beta:
            for r in range(M):
                if r in bset:
                    continue
                new = tuple(sorted((bset - {p}) | {r}))
                if new not in b_idx:
                    continue
                bin_ = b_idx[new]
                ph = _excitation_phase(beta, p, r)
                for ai, alpha in enumerate(alpha_strings):
                    val = ph * h1[r, p]
                    for q in beta:
                        if q == p:
                            continue
                        val += ph * (eri[r, p, q, q] - eri[r, q, q, p])
                    for q in alpha:
                        val += ph * eri[r, p, q, q]
                    if abs(val) > 1e-14:
                        H[di(ai, bi), di(ai, bin_)] += val

    # double excitations within alpha
    for ai, alpha in enumerate(alpha_strings):
        aset = set(alpha)
        occ = list(alpha)
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
                        ph = _double_excitation_phase(alpha, p, q, r, s)
                        val = ph * (eri[r, p, s, q] - eri[r, q, s, p])
                        if abs(val) > 1e-14:
                            for bi in range(nb):
                                H[di(ai, bi), di(ain, bi)] += val

    # double excitations within beta
    for bi, beta in enumerate(beta_strings):
        bset = set(beta)
        occ = list(beta)
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
                        ph = _double_excitation_phase(beta, p, q, r, s)
                        val = ph * (eri[r, p, s, q] - eri[r, q, s, p])
                        if abs(val) > 1e-14:
                            for ai in range(na):
                                H[di(ai, bi), di(ai, bin_)] += val

    # alpha-beta double excitations
    for ai, alpha in enumerate(alpha_strings):
        aset = set(alpha)
        for pa in alpha:
            for ra in range(M):
                if ra in aset:
                    continue
                newa = tuple(sorted((aset - {pa}) | {ra}))
                if newa not in a_idx:
                    continue
                ain = a_idx[newa]
                pha = _excitation_phase(alpha, pa, ra)
                for bi, beta in enumerate(beta_strings):
                    bset = set(beta)
                    for pb in beta:
                        for rb in range(M):
                            if rb in bset:
                                continue
                            newb = tuple(sorted((bset - {pb}) | {rb}))
                            if newb not in b_idx:
                                continue
                            bin_ = b_idx[newb]
                            phb = _excitation_phase(beta, pb, rb)
                            val = pha * phb * eri[ra, pa, rb, pb]
                            if abs(val) > 1e-14:
                                H[di(ai, bi), di(ain, bin_)] += val

    H = H.tocsr()
    H = 0.5 * (H + H.T)  # symmetrize (builder fills upper/lower asymmetrically)
    k = min(4, n_det - 1)
    if k < 1:
        vals = np.array([H[0, 0]])
        vecs = np.array([[1.0]])
    else:
        vals, vecs = eigsh(H, k=k, which='SA')
    order = np.argsort(vals)
    E_gs = float(vals[order[0]])
    vec = np.asarray(vecs[:, order[0]]).ravel()
    return E_gs, vec, alpha_strings, beta_strings, nb


def site_of_spatial(result):
    site = []
    labels = []
    for b in result['blocks']:
        is_H = b['label'].endswith('_partner')
        site.extend([1 if is_H else 0] * b['n_orbitals'])
        labels.append(b['label'])
    return np.array(site, dtype=int), labels


def site_entanglement(vec, alpha_strings, beta_strings, nb, site_spatial, M):
    """Fermionic mode-entanglement (nats) across heavy(site 0)|H(site 1) cut.

    Determinant = (a†_{alpha,up} ascending)(a†_{beta,down} ascending)|0>.
    Global spin-orbital index: up spatial p -> p; down spatial q -> M+q.
    Site of global g = site_spatial[g if g<M else g-M].
    Reorder operator product to (A ascending)(B ascending); sign = (-1)^(inv),
    inv = #pairs (a in A, b in B) with b before a in the builder order.
    """
    def site_g(g):
        return site_spatial[g if g < M else g - M]

    amps = {}
    for ai, alpha in enumerate(alpha_strings):
        for bi, beta in enumerate(beta_strings):
            amp = vec[ai * nb + bi]
            if abs(amp) < 1e-13:
                continue
            occ = [p for p in alpha] + [M + q for q in beta]
            seenA = 0
            inv_before = 0  # sum over B of (#A before it)
            A = []
            B = []
            for g in occ:
                if site_g(g) == 0:
                    seenA += 1
                    A.append(g)
                else:
                    inv_before += seenA
                    B.append(g)
            inv = len(A) * len(B) - inv_before  # #A after each B
            sign = -1.0 if (inv % 2) else 1.0
            kA = tuple(sorted(A))
            kB = tuple(sorted(B))
            amps[(kA, kB)] = amps.get((kA, kB), 0.0) + sign * amp

    if not amps:
        return 0.0
    keysA = sorted({k[0] for k in amps})
    keysB = sorted({k[1] for k in amps})
    iA = {k: i for i, k in enumerate(keysA)}
    iB = {k: i for i, k in enumerate(keysB)}
    psi = np.zeros((len(keysA), len(keysB)))
    for (kA, kB), v in amps.items():
        psi[iA[kA], iB[kB]] = v
    sv = np.linalg.svd(psi, compute_uv=False)
    p = sv ** 2
    p = p[p > 1e-15]
    p = p / p.sum()
    return float(-np.sum(p * np.log(p)))


def measure(spec_fn, name, Z, n_electrons, n_sites, R_grid):
    spec = spec_fn()
    rows = []
    print(f"\n{name} (Z={Z}): {n_electrons}e, {n_sites} sites, "
          f"ceiling ln{n_sites}={math.log(n_sites):.4f}", flush=True)
    for j, R in enumerate(R_grid):
        result = build_balanced_hamiltonian(
            spec, R=R, screened_cross_center=True,
            multi_zeta_basis=False, cross_block_h1=True,
        )
        M = result['M']
        site_spatial, labels = site_of_spatial(result)
        if j == 0:
            nH = int((site_spatial == 1).sum())
            print(f"  M={M}  heavy={M-nH} H={nH}  blocks={labels}", flush=True)
        E, vec, a_s, b_s, nb = fci_ground(result, n_electrons)
        S = site_entanglement(vec, a_s, b_s, nb, site_spatial, M)
        rows.append({'R_bohr': float(R), 'E_Ha': float(E), 'S_site_nats': float(S)})
        print(f"  R={R:5.2f}  E={E:13.5f}  S_site={S:.5f}", flush=True)
    return {'molecule': name, 'Z': Z, 'n_electrons': n_electrons,
            'n_sites': n_sites, 'ln_n_sites': math.log(n_sites), 'rows': rows}


def main():
    R_grid = [2.0, 2.5, 3.0, 3.566, 4.0, 5.0, 7.0, 10.0]
    out = {}
    print("=" * 70)
    print("Species-II signature: which-site entanglement curve S_site(R)")
    print("Decisive pair: NaH (Z=11) vs KH (Z=19), isostructural 2-site 2e")
    print("Architecture: screened + cross_block_h1, multi_zeta OFF (Na-only)")
    print("Pre-committed Q1: does S_site(R) plateau? Q2: NaH ~ KH (Z-invariant)?")
    print("=" * 70)
    out['NaH'] = measure(nah_spec, 'NaH', 11, 2, 2, R_grid)
    out['KH'] = measure(kh_spec, 'KH', 19, 2, 2, R_grid)

    na = {r['R_bohr']: r['S_site_nats'] for r in out['NaH']['rows']}
    kh = {r['R_bohr']: r['S_site_nats'] for r in out['KH']['rows']}
    diffs = [abs(na[R] - kh[R]) for R in R_grid]
    ln2 = math.log(2)
    print("\n" + "=" * 70)
    print("CURVE COMPARISON (the data — interpret ceiling AFTER seeing it)")
    print(f"{'R':>6} {'S(NaH)':>10} {'S(KH)':>10} {'|diff|':>9}")
    for R, d in zip(R_grid, diffs):
        print(f"{R:>6.2f} {na[R]:>10.5f} {kh[R]:>10.5f} {d:>9.5f}")
    print(f"\nreference: ln2={ln2:.5f}  ln4=2ln2={2*ln2:.5f}")
    print(f"large-R: S(NaH)={na[R_grid[-1]]:.5f}  S(KH)={kh[R_grid[-1]]:.5f}")
    print(f"max|S(NaH)-S(KH)| = {max(diffs):.5f}   mean = {np.mean(diffs):.5f}")
    print("\nQ1 saturation: do both curves plateau at large R?")
    print("Q2 Z-invariance: max|diff| small (curves coincide)?")

    out['_meta'] = {
        'R_grid': R_grid, 'ln2': ln2, 'ln4': 2 * ln2,
        'max_abs_diff_NaH_KH': float(max(diffs)),
        'mean_abs_diff_NaH_KH': float(np.mean(diffs)),
        'S_NaH_largeR': na[R_grid[-1]], 'S_KH_largeR': kh[R_grid[-1]],
        'architecture': 'screened+cross_block_h1, multi_zeta=False',
    }
    Path('debug/data').mkdir(parents=True, exist_ok=True)
    Path('debug/data/species_ii_ordering.json').write_text(json.dumps(out, indent=2))
    print("\nWrote debug/data/species_ii_ordering.json")


if __name__ == '__main__':
    main()
