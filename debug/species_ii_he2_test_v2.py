"""He2 responsiveness test v2 — CORRECTED build (two coupled single-center He shells).

v1 (debug/species_ii_he2_test.py) used a single bond-block and the cross_h1 GUARD
caught it: cross_h1 == 0 at ALL R (even 1.5 bohr) => the two centers never coupled =>
product state, artifact, NOT a real frozen signal. Same latent issue our memo flagged
for the H2 single-bond-block control.

Fix (instrument, NOT threshold): build He2 as TWO single-center He blocks (one closed
1s^2 shell per nucleus), like N2 puts a block on each center. cross_block_h1 (s-s) +
cross-center V_ne now actually couple them. Partition by NUCLEUS, not by '_partner' label.

Pre-registered falsifier UNCHANGED (debug/species_ii_he2_prediction_registered.md):
  FROZEN     if range(S_mode, R>=4) < 0.02  -> interaction-mechanism (B); my bet
  RESPONSIVE if range(S_mode, R>=4) > 0.05  -> pair-count (A); my bet falsified

GUARD retained: if cross_h1 is still ~0 at small R, the build is STILL decoupled and
the result is invalid (report, do not spin).
"""

from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import debug.species_ii_discriminator as disc
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import MolecularSpec, OrbitalBlock


def he2_spec_v2(R=5.6, max_n=2):
    """He2 as two single-center closed 1s^2 shells on two nuclei (coupled like N2)."""
    blocks = [
        OrbitalBlock(label='HeA', block_type='lone_pair', Z_center=2.0,
                     n_electrons=2, max_n=max_n, center_nucleus_idx=0,
                     Z_nuc_center=2.0),
        OrbitalBlock(label='HeB', block_type='lone_pair', Z_center=2.0,
                     n_electrons=2, max_n=max_n, center_nucleus_idx=1,
                     Z_nuc_center=2.0),
    ]
    nuclei = [
        {'Z': 2.0, 'position': (0.0, 0.0, 0.0), 'label': 'HeA'},
        {'Z': 2.0, 'position': (0.0, 0.0, float(R)), 'label': 'HeB'},
    ]
    return MolecularSpec(
        name='He2', blocks=blocks, nuclear_repulsion_constant=4.0 / R,
        description='He2: two coupled single-center 1s^2 shells; count-vs-interaction',
        nuclei=nuclei,
    )


def he2_site(result):
    """Partition spatial orbitals by NUCLEUS: HeA-derived -> 0, HeB-derived -> 1."""
    site = []
    for b in result['blocks']:
        lbl = b['label']
        s = 1 if ('HeB' in lbl or lbl.endswith('_partner')) else 0
        site.extend([s] * b['n_orbitals'])
    return np.array(site, dtype=int)


def cross_site_coupling(result, site):
    h1 = result['h1']; M = len(site); mx = 0.0
    for p in range(M):
        for q in range(M):
            if site[p] != site[q]:
                mx = max(mx, abs(float(h1[p, q])))
    return mx


def main():
    R_grid = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.6, 6.0, 7.0, 8.0]
    print("=" * 74)
    print("He2 v2 (CORRECTED): two coupled single-center He shells")
    print("GUARD: cross_h1 must be NONZERO at small R or build is still decoupled")
    print("=" * 74)

    # one-shot block-structure check at the smallest R
    spec0 = he2_spec_v2(R=1.5)
    r0 = build_balanced_hamiltonian(spec0, R=1.5, nuclei=spec0.nuclei,
                                    screened_cross_center=True,
                                    multi_zeta_basis=False, cross_block_h1=True)
    s0 = he2_site(r0)
    print(f"blocks: {[(b['label'], b['n_orbitals']) for b in r0['blocks']]}")
    print(f"site vector: {s0.tolist()}   (n0={int(np.sum(s0==0))}, n1={int(np.sum(s0==1))})")
    print(f"cross_h1 at R=1.5: {cross_site_coupling(r0, s0):.4f}  "
          f"(MUST be > 0 for a valid coupled build)\n")

    rows = []
    print("He2 v2 sweep:", flush=True)
    for R in R_grid:
        spec = he2_spec_v2(R=R)
        result = build_balanced_hamiltonian(spec, R=R, nuclei=spec.nuclei,
                                            screened_cross_center=True,
                                            multi_zeta_basis=False, cross_block_h1=True)
        M = result['M']; site = he2_site(result)
        vals, vecs, a_str, b_str, nb, MM = disc.build_fci(result, 4, k=8)
        gap = float(vals[1] - vals[0]) if len(vals) > 1 else float('inf')
        d = int(np.sum(vals - vals[0] < 1e-6))
        manifold = [vecs[:, i] for i in range(d)]
        S_mode, rank = disc.ensemble_mode_entropy(manifold, a_str, b_str, nb, site, M)
        S_occ, _ = disc.ensemble_occ_entropy(manifold, a_str, b_str, nb, site, M)
        coup = cross_site_coupling(result, site)
        rows.append({'R_bohr': float(R), 'E_Ha': float(vals[0]), 'gap_Ha': gap,
                     'degeneracy': d, 'S_mode': S_mode, 'schmidt_rank': rank,
                     'S_occ': S_occ, 'cross_site_coupling': coup})
        print(f"  R={R:5.2f}  E={vals[0]:10.4f}  d={d}  S_mode={S_mode:.5f}  "
              f"S_occ={S_occ:.5f}  rank={rank}  cross_h1={coup:.5f}", flush=True)

    smode = np.array([r['S_mode'] for r in rows])
    Rs = np.array([r['R_bohr'] for r in rows])
    coup = np.array([r['cross_site_coupling'] for r in rows])
    full_range = float(smode.max() - smode.min())
    large = smode[Rs >= 4.0]; large_range = float(large.max() - large.min())
    small = smode[Rs < 4.0]; small_range = float(small.max() - small.min())

    # validity guard
    coup_small = coup[Rs < 3.0]
    build_valid = bool(coup_small.max() > 1e-4)

    if not build_valid:
        verdict = "INVALID: cross_h1 ~0 at small R, centers STILL decoupled"
        head = "build-invalid"
    elif large_range < 0.02:
        verdict = "FROZEN at large R -> interaction-mechanism (B) -> MY BET HOLDS"
        head = "frozen-interaction-mechanism"
    elif large_range > 0.05:
        verdict = "RESPONSIVE at large R -> pair-count (A) -> MY BET FALSIFIED"
        head = "responsive-pair-count"
    else:
        verdict = "AMBIGUOUS (large-R range 0.02-0.05) -> no spin"
        head = "ambiguous"

    if smode.std() > 1e-9 and coup.std() > 1e-9:
        corr_sc = float(np.corrcoef(smode, coup)[0, 1])
    else:
        corr_sc = None

    print("\n" + "=" * 74)
    print(f"build_valid (cross_h1>0 at small R) = {build_valid}  "
          f"[cross_h1(R<3) max = {coup_small.max():.5f}]")
    print(f"S_mode full range           = {full_range:.5f}")
    print(f"S_mode range, LARGE R (>=4) = {large_range:.5f}   <- PRE-REGISTERED TEST")
    print(f"S_mode range, small R (<4)  = {small_range:.5f}")
    print(f"corr(S_mode, cross_h1)      = {corr_sc}   (my bet: response tracks coupling)")
    print("Reference: 1-pair frozen ~0.000-0.005; 2-pair INTERACTING (LiH) 0.20+")
    print("-" * 74)
    print("VERDICT:", verdict)
    print("=" * 74)

    out = {'system': 'He2_v2_two_center', 'n_electrons': 4, 'pairs': 2,
           'pair_interaction': 'NON-interacting closed 1s^2 shells (now coupled in build)',
           'R_grid': [float(x) for x in R_grid], 'rows': rows,
           'build_valid': build_valid, 'cross_h1_small_R_max': float(coup_small.max()),
           'S_mode_full_range': full_range,
           'S_mode_large_R_range_R>=4': large_range,
           'S_mode_small_R_range_R<4': small_range,
           'corr_Smode_cross_h1': corr_sc,
           'verdict': verdict, 'headline': head,
           'prereg_file': 'debug/species_ii_he2_prediction_registered.md'}
    Path('debug/data').mkdir(parents=True, exist_ok=True)
    Path('debug/data/species_ii_he2_test_v2.json').write_text(json.dumps(out, indent=2))
    print("Wrote debug/data/species_ii_he2_test_v2.json")


if __name__ == '__main__':
    main()
