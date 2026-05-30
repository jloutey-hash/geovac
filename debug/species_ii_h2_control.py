"""De-confounder: H2, a 2-ELECTRON BINDER, vs NaH (2e non-binder) and LiH (4e binder).

The discriminator sprint found: non-binders (NaH/KH) have R-INDEPENDENT (frozen)
which-site entanglement; the binder (LiH) has strongly R-responsive entanglement.
But LiH differs from NaH in 3 ways: 4e vs 2e, explicit core vs frozen core, AND
it binds. So "responsive vs frozen" is confounded.

H2 is the clean control: 2 electrons (like NaH), NO core, and it BINDS (like LiH).
  - If H2 entanglement is RESPONSIVE (varies with R like LiH) -> responsiveness
    tracks BINDING, not electron count. Discriminator survives. SUPPORTED.
  - If H2 entanglement is FROZEN (flat like NaH) -> responsiveness tracks
    electron-count/core-treatment, NOT binding. Discriminator confounded. WEAKENED.

H2 spec: one bond block, Z_center=Z_partner=1, 2 electrons, no frozen/explicit
core, nuclear_repulsion = Z_A*Z_B/R = 1/R. Same balanced pipeline, same FCI,
same entanglement measure as the discriminator sprint.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

import debug.species_ii_discriminator as disc
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import MolecularSpec, OrbitalBlock


def h2_spec(R=1.4, max_n=2):
    """Minimal H2: single bond block, both centers Z=1, 2 electrons, no core."""
    blocks = [OrbitalBlock(
        label='H2_bond', block_type='bond', Z_center=1.0, n_electrons=2,
        max_n=max_n, has_h_partner=True, Z_partner=1.0, max_n_partner=max_n,
        Z_nuc_center=1.0,
    )]
    return MolecularSpec(
        name='H2', blocks=blocks,
        nuclear_repulsion_constant=1.0 / R,  # placeholder; rebuilt per R below
        description='H2: 2e, no core, 2-electron binder control',
    )


def measure_h2(R_grid):
    rows = []
    print("\nH2 (2e binder, no core)", flush=True)
    for j, R in enumerate(R_grid):
        spec = h2_spec(R=R)
        # set nuclear repulsion for this R (1/R, Z_A=Z_B=1)
        spec.nuclear_repulsion_constant = 1.0 / R
        nuclei = [
            {"Z": 1.0, "position": (0.0, 0.0, 0.0), "label": "Ha"},
            {"Z": 1.0, "position": (0.0, 0.0, float(R)), "label": "Hb"},
        ]
        result = build_balanced_hamiltonian(
            spec, R=R, nuclei=nuclei, screened_cross_center=True,
            multi_zeta_basis=False, cross_block_h1=True)
        # balanced may overwrite nuclear_repulsion; ensure 1/R for a no-core diatomic
        result['nuclear_repulsion'] = 1.0 / R
        M = result['M']; site = disc.site_of_spatial(result)
        vals, vecs, a_str, b_str, nb, MM = disc.build_fci(result, 2, k=8)
        gap = float(vals[1] - vals[0]) if len(vals) > 1 else float('inf')
        d = int(np.sum(vals - vals[0] < 1e-6))
        manifold = [vecs[:, i] for i in range(d)]
        S_mode, rank = disc.ensemble_mode_entropy(manifold, a_str, b_str, nb, site, M)
        S_occ, distn = disc.ensemble_occ_entropy(manifold, a_str, b_str, nb, site, M)
        ceil = math.log(rank) if rank > 1 else 0.0
        ratio = S_mode / ceil if ceil > 0 else 0.0
        rows.append({'R_bohr': float(R), 'E_Ha': float(vals[0]), 'gap_Ha': gap,
                     'degeneracy': d, 'S_mode': S_mode, 'schmidt_rank': rank,
                     'S_mode_ceiling_lnrank': ceil, 'S_mode_over_ceiling': ratio,
                     'S_occ': S_occ, 'occ_dist': distn})
        print(f"  R={R:5.2f}  E={vals[0]:9.4f}  gap={gap:7.4f}  d={d}  "
              f"S_mode={S_mode:.4f}/ln{rank}={ceil:.4f}={ratio:.3f}  S_occ={S_occ:.4f}",
              flush=True)
    return rows


def main():
    R_grid = [1.0, 1.4, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0]
    print("=" * 72)
    print("H2 de-confounder: 2e binder. Frozen(NaH-like) or responsive(LiH-like)?")
    print("=" * 72)
    rows = measure_h2(R_grid)
    ratios = [r['S_mode_over_ceiling'] for r in rows]
    smodes = [r['S_mode'] for r in rows]
    rng_ratio = max(ratios) - min(ratios)
    rng_smode = max(smodes) - min(smodes)
    print("\n" + "=" * 72)
    print(f"H2 ratio range over PES = {rng_ratio:.4f}  (S_mode range = {rng_smode:.4f})")
    print("Reference from discriminator sprint:")
    print("  NaH (2e, no bind): ratio range 0.004   (FROZEN)")
    print("  KH  (2e, no bind): ratio range 0.002   (FROZEN)")
    print("  LiH (4e, BINDS):   ratio range 0.203   (RESPONSIVE)")
    print(f"  H2  (2e, BINDS):   ratio range {rng_ratio:.4f}   "
          f"-> {'RESPONSIVE (tracks BINDING)' if rng_ratio > 0.05 else 'FROZEN (tracks e-count, confounded)'}")
    # also find energy minimum (does H2 bind in this framework?)
    Es = [r['E_Ha'] for r in rows]
    imin = int(np.argmin(Es))
    print(f"\nH2 PES min at R={R_grid[imin]:.2f} bohr (E={Es[imin]:.4f}); "
          f"E(R=10)={Es[-1]:.4f}; binds_in_framework={imin not in (0, len(Es)-1)}")
    out = {'H2': {'rows': rows, 'ratio_range': rng_ratio, 'smode_range': rng_smode},
           '_ref': {'NaH_ratio_range': 0.004, 'KH_ratio_range': 0.002,
                    'LiH_ratio_range': 0.203}}
    Path('debug/data/species_ii_h2_control.json').write_text(json.dumps(out, indent=2))
    print("\nWrote debug/data/species_ii_h2_control.json")


if __name__ == '__main__':
    main()
