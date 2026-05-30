"""ln3 replication: MgH2 (Z=12) vs CaH2 (Z=20), 3-site 4e isostructural pair.

The clincher for the species-II saturation prediction. NaH/KH pinned at
ln2=0.6931 (2 sites). If MgH2/CaH2 pin at ln3=1.0986 (3 sites) — a DIFFERENT
ceiling value the curve cannot fake — the "saturates at ln(N_sites)" reading
survives a second, harder test, and Z-invariance is replicated at a new
topology. 4 electrons, M=20, FCI dim C(20,2)^2 = 36100; coarse R grid.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

import debug.species_ii_ordering as s
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import mgh2_spec, cah2_spec


def site_of_spatial_3way(result):
    """3-way site label per spatial orbital: 0=heavy, 1=H_a, 2=H_b.

    Partner sub-blocks are labelled '..._bond_1_partner', '..._bond_2_partner';
    each distinct partner index is its own H site. Everything else = heavy.
    """
    site = []
    h_index = {}
    for b in result['blocks']:
        lbl = b['label']
        if lbl.endswith('_partner'):
            # extract bond index for distinct H sites
            key = lbl
            if key not in h_index:
                h_index[key] = len(h_index) + 1  # 1, 2, ...
            site.extend([h_index[key]] * b['n_orbitals'])
        else:
            site.extend([0] * b['n_orbitals'])
    return np.array(site, dtype=int), sorted(h_index.keys())


def site_occupation_entropy(vec, alpha_strings, beta_strings, nb, site3, M, n_sites):
    """Shannon entropy (nats) of the electron-count-per-site distribution.

    Unambiguous, gauge-free, no fermion signs needed: for each determinant,
    count how many electrons sit on each site (heavy=0, H_a=1, H_b=2), weight
    by |amp|^2, accumulate the joint distribution over occupation patterns, and
    return its Shannon entropy. Max for n_sites equally-likely single-site
    patterns at dissociation is ln(n_sites). This is the coarse-grained
    'which-site' classical entropy — the species-II saturation observable.
    """
    def site_g(g):
        return site3[g if g < M else g - M]

    dist = {}
    for ai, alpha in enumerate(alpha_strings):
        for bi, beta in enumerate(beta_strings):
            w = abs(vec[ai * nb + bi]) ** 2
            if w < 1e-14:
                continue
            counts = [0] * n_sites
            for p in alpha:
                counts[site_g(p)] += 1
            for q in beta:
                counts[M + q if False else site_g(M + q)] += 1
            key = tuple(counts)
            dist[key] = dist.get(key, 0.0) + w
    p = np.array(list(dist.values()))
    p = p[p > 1e-15]
    p = p / p.sum()
    return float(-np.sum(p * np.log(p))), dist


def measure(spec_fn, name, Z, n_electrons, n_sites, R_grid):
    spec = spec_fn()
    rows = []
    print(f"\n{name} (Z={Z}): {n_electrons}e, {n_sites} sites, "
          f"ceiling ln{n_sites}={math.log(n_sites):.4f}", flush=True)
    for j, R in enumerate(R_grid):
        result = build_balanced_hamiltonian(
            spec, R=R, screened_cross_center=True,
            multi_zeta_basis=True, cross_block_h1=True,
        )
        M = result['M']
        site_spatial, labels = s.site_of_spatial(result)
        if j == 0:
            nH = int((site_spatial == 1).sum())
            print(f"  M={M}  heavy={M-nH} H={nH}  blocks={labels}", flush=True)
        E, vec, a_s, b_s, nb = s.fci_ground(result, n_electrons)
        S_bipart = s.site_entanglement(vec, a_s, b_s, nb, site_spatial, M)
        site3, hkeys = site_of_spatial_3way(result)
        n_site3 = int(site3.max()) + 1
        S_occ, dist = site_occupation_entropy(vec, a_s, b_s, nb, site3, M, n_site3)
        rows.append({'R_bohr': float(R), 'E_Ha': float(E),
                     'S_bipartite_nats': float(S_bipart),
                     'S_occ_3way_nats': float(S_occ),
                     'n_site3': n_site3})
        print(f"  R={R:5.2f}  E={E:13.5f}  S_bipart(heavy|H)={S_bipart:.5f}  "
              f"S_occ_3way={S_occ:.5f}  (ln3={math.log(3):.4f})", flush=True)
    return {'molecule': name, 'Z': Z, 'n_electrons': n_electrons,
            'n_sites': n_sites, 'ln_n_sites': math.log(n_sites), 'rows': rows}


def main():
    # Note: site partition counts ALL H-partner orbitals as "site 1" (lumped).
    # For 3 physical sites (Mg + 2 H) the heavy|H bipartition gives a 2-way cut,
    # whose saturation ceiling is still ln2 unless we cut 3 ways. So we ALSO
    # report a 3-way site entropy below. First the heavy|H bipartite S_site.
    R_grid = [2.5, 3.5, 5.0, 8.0]
    out = {}
    print("=" * 70)
    print("ln(N_sites) replication: MgH2 (Z=12) vs CaH2 (Z=20), 4e")
    print("=" * 70)
    out['MgH2'] = measure(mgh2_spec, 'MgH2', 12, 4, 3, R_grid)
    out['CaH2'] = measure(cah2_spec, 'CaH2', 20, 4, 3, R_grid)

    # PRIMARY: bipartite heavy|H mode entanglement (clean, ceiling ln2 for this
    # 2-way cut). Tests Z-invariance replication at 4 electrons.
    mg = {r['R_bohr']: r['S_bipartite_nats'] for r in out['MgH2']['rows']}
    ca = {r['R_bohr']: r['S_bipartite_nats'] for r in out['CaH2']['rows']}
    diffs = [abs(mg[R] - ca[R]) for R in R_grid]
    print("\n" + "=" * 70)
    print("PRIMARY: bipartite heavy|H mode entanglement (ceiling ln2; 4e pair)")
    print(f"{'R':>6} {'S(MgH2)':>10} {'S(CaH2)':>10} {'|diff|':>9}")
    for R, d in zip(R_grid, diffs):
        print(f"{R:>6.2f} {mg[R]:>10.5f} {ca[R]:>10.5f} {d:>9.5f}")
    print(f"\nln2={math.log(2):.5f}  ln3={math.log(3):.5f}")
    print(f"max|MgH2-CaH2| bipartite = {max(diffs):.5f}  mean = {np.mean(diffs):.5f}")
    print("Q2-replication: MgH2 ~ CaH2 (Z=12 vs 20) on a 4-electron pair?")
    # SECONDARY (diagnostic only, NOT a ceiling test): 3-way occupation entropy
    mg3 = {r['R_bohr']: r['S_occ_3way_nats'] for r in out['MgH2']['rows']}
    ca3 = {r['R_bohr']: r['S_occ_3way_nats'] for r in out['CaH2']['rows']}
    print("\nSECONDARY (correlation diagnostic, not a saturation ceiling):")
    print(f"{'R':>6} {'occ(MgH2)':>11} {'occ(CaH2)':>11}")
    for R in R_grid:
        print(f"{R:>6.2f} {mg3[R]:>11.5f} {ca3[R]:>11.5f}")

    out['_meta'] = {'R_grid': R_grid, 'ln2': math.log(2), 'ln3': math.log(3),
                    'max_abs_diff_bipartite': float(max(diffs)),
                    'mean_abs_diff_bipartite': float(np.mean(diffs)),
                    'S_MgH2_largeR_bipartite': mg[R_grid[-1]],
                    'S_CaH2_largeR_bipartite': ca[R_grid[-1]]}
    Path('debug/data/species_ii_ln3.json').write_text(json.dumps(out, indent=2))
    print("\nWrote debug/data/species_ii_ln3.json")


if __name__ == '__main__':
    main()
