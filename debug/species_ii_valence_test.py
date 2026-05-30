"""Valence-freedom test (Josh's hypothesis): MgH2 splits the hypotheses.

PRE-REGISTERED predictions (locked before running, per audit-claim discipline).
Measured so far (debug/data/species_ii_discriminator.json + h2_control.json):
  H2  : 2e, 1 valence pair, binds(spurious)  -> entanglement FROZEN (range 0.000)
  NaH : 2e, 1 valence pair, no bind          -> FROZEN (range 0.004)
  KH  : 2e, 1 valence pair, no bind          -> FROZEN (range 0.002)
  LiH : 4e, 1 valence pair + 1 explicit core -> RESPONSIVE (range 0.203)

MgH2: 4 active electrons, frozen [Ne] core (NO explicit core pair), 2 valence
pairs (Mg 3s^2 + 2 H), does NOT bind in-framework (W1c wall). It splits the
four hypotheses for "what makes entanglement responsive":
  H1 total electron count        -> MgH2 has 4    -> predict RESPONSIVE
  H3a valence freedom (Josh)     -> 2 valence pairs-> predict RESPONSIVE
  H2  binding                    -> doesn't bind  -> predict FROZEN
  H3b explicit core-valence corr -> no expl. core -> predict FROZEN

Threshold (pre-set): RESPONSIVE if ratio-range > 0.05; FROZEN if < 0.02.
(LiH ref = 0.203 responsive; 2e set = 0.002-0.004 frozen; 0.02-0.05 = ambiguous.)

  -> MgH2 RESPONSIVE  : supports {electron-count, valence-freedom}; kills {binding, core-valence}
  -> MgH2 FROZEN      : supports {binding, core-valence}; kills {electron-count, valence-freedom}

Either way: cuts 4 hypotheses to 2. Honest caveat baked in: MgH2 is 3-site
(lumped heavy|H bipartite cut), a geometry difference from the 2-site set; note
it, don't hide it.

Architecture identical to discriminator (screened + cross_block_h1, multi_zeta
OFF). FCI particle-number-projected, validated bit-exact vs library.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

import debug.species_ii_discriminator as disc
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import mgh2_spec


def main():
    R_grid = [2.0, 2.5, 3.0, 3.566, 4.0, 5.0, 7.0, 10.0]
    print("=" * 72)
    print("VALENCE-FREEDOM TEST: MgH2 (4e, all-valence, 2 valence pairs, no bind)")
    print("Predictions: count->RESP  valence->RESP  binding->FROZEN  core-val->FROZEN")
    print("=" * 72)

    # validate FCI vs library at one point
    spec = mgh2_spec()
    res0 = build_balanced_hamiltonian(spec, R=3.566, screened_cross_center=True,
                                      multi_zeta_basis=False, cross_block_h1=True)
    ref = coupled_fci_energy(res0, 4)['E_coupled']
    vv, *_ = disc.build_fci(res0, 4, k=4)
    print(f"FCI validation vs library: diff={abs(vv[0]-ref):.2e} "
          f"{'OK' if abs(vv[0]-ref) < 1e-6 else 'MISMATCH'}")

    rows = []
    print("\nMgH2 (Z=12, 4e)")
    for R in R_grid:
        result = build_balanced_hamiltonian(spec, R=R, screened_cross_center=True,
                                            multi_zeta_basis=False, cross_block_h1=True)
        M = result['M']; site = disc.site_of_spatial(result)
        vals, vecs, a_str, b_str, nb, MM = disc.build_fci(result, 4, k=8)
        gap = float(vals[1] - vals[0]) if len(vals) > 1 else float('inf')
        d = int(np.sum(vals - vals[0] < 1e-6))
        manifold = [vecs[:, i] for i in range(d)]
        S_mode, rank = disc.ensemble_mode_entropy(manifold, a_str, b_str, nb, site, M)
        ceil = math.log(rank) if rank > 1 else 0.0
        ratio = S_mode / ceil if ceil > 0 else 0.0
        nH = int((site == 1).sum())
        rows.append({'R_bohr': float(R), 'E_Ha': float(vals[0]), 'gap_Ha': gap,
                     'degeneracy': d, 'S_mode': S_mode, 'schmidt_rank': rank,
                     'S_mode_over_ceiling': ratio, 'n_heavy': M - nH, 'n_H': nH})
        print(f"  R={R:5.2f}  E={vals[0]:10.4f}  gap={gap:7.4f}  d={d}  "
              f"S_mode={S_mode:.4f}/ln{rank}={ceil:.4f}={ratio:.3f}", flush=True)

    ratios = [r['S_mode_over_ceiling'] for r in rows]
    rng = max(ratios) - min(ratios)
    Es = [r['E_Ha'] for r in rows]
    imin = int(np.argmin(Es))
    binds = imin not in (0, len(Es) - 1)

    if rng > 0.05:
        verdict = "RESPONSIVE -> supports {electron-count, VALENCE-FREEDOM}; kills {binding, core-valence}"
    elif rng < 0.02:
        verdict = "FROZEN -> supports {binding, core-valence}; kills {electron-count, valence-freedom}"
    else:
        verdict = "AMBIGUOUS (0.02-0.05) -> inconclusive"

    print("\n" + "=" * 72)
    print(f"MgH2 ratio range over PES = {rng:.4f}")
    print(f"  ref: LiH(4e,binds)=0.203 RESPONSIVE | 2e set=0.002-0.004 FROZEN")
    print(f"  MgH2 binds in-framework: {binds} (PES min at R={R_grid[imin]:.2f})")
    print(f"VERDICT: {verdict}")
    print("CAVEAT: MgH2 is 3-site (heavy|H bipartite lumps 2 H); geometry differs")
    print("        from the 2-site set. Real but note it.")

    out = {'MgH2': {'Z': 12, 'n_electrons': 4, 'rows': rows,
                    'ratio_range': rng, 'binds': binds},
           '_predictions': {
               'electron_count': 'RESPONSIVE', 'valence_freedom_josh': 'RESPONSIVE',
               'binding': 'FROZEN', 'core_valence': 'FROZEN'},
           '_ref': {'LiH': 0.203, 'NaH': 0.004, 'KH': 0.002, 'H2': 0.000},
           '_verdict': verdict}
    Path('debug/data/species_ii_valence_test.json').write_text(json.dumps(out, indent=2))
    print("\nWrote debug/data/species_ii_valence_test.json")


if __name__ == '__main__':
    main()
