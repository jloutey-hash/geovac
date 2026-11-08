"""Compute joint (Hopf-Z2 sector) x (S^2 eigenvalue) block dimensions.

After the QPT stackability test (qpt_hopf_stacking_test.py) confirmed
all commutators are bit-exact zero, this script makes the resource
savings concrete by counting state-space dimensions per joint sector.

For LiH (the smallest molecule), tabulate:
  dim(N_alpha = a, N_beta = b, P_0 = s0, P_1 = s1, ..., S^2 = j(j+1))

Compare:
  - Full Hilbert (2^Q)
  - Particle-number sector (combinatorial)
  - + Hopf-Z2 per-sub-block sector (factor 2^n_sb)
  - + S^2 sector (combinatorial spin coupling)

The cleanest savings metric: dim(target sector) / dim(particle-number-only)
gives the additional reduction from QPT + Hopf above the standard
Bravyi 2017 particle-number/spin-parity tapering.
"""

from __future__ import annotations

import json
import os
from typing import Dict
from itertools import combinations

import numpy as np

from geovac.composed_qubit import build_composed_hamiltonian
from geovac.molecular_spec import lih_spec
from geovac.z2_tapering import build_pm_rotation, _enumerate_orbitals


def count_sector_dimensions(M: int, parity: np.ndarray,
                            orbital_table, n_electrons: int,
                            target_S2: float = 0.0) -> Dict:
    """Enumerate Slater determinants and bucket by symmetry sectors.

    For each SD |I, J> with alpha-set I and beta-set J:
      N_alpha = |I|, N_beta = |J|
      P_i = sign over antisym orbitals in sub-block i

    Note: a single SD is not an S^2 eigenstate in general. We report SD
    counts in each (N, P) sector and the dimension of the spin-adapted
    CSF subspace.
    """
    # Group antisym orbitals by sub-block
    sb_to_antisym = {}
    for p, (sb, n, l, m) in enumerate(orbital_table):
        sb_to_antisym.setdefault(sb, [])
        if parity[p] == -1:
            sb_to_antisym[sb].append(p)
    # Sub-blocks with at least one antisym orbital
    active_sb = [sb for sb, anti in sb_to_antisym.items() if anti]
    n_P = len(active_sb)

    # Enumerate (N_alpha, N_beta) decompositions of n_electrons.
    n_alpha_range = range(0, min(n_electrons, M) + 1)

    sector_counts = {}
    for n_alpha in n_alpha_range:
        n_beta = n_electrons - n_alpha
        if n_beta < 0 or n_beta > M:
            continue
        sz = (n_alpha - n_beta) / 2.0

        # Enumerate alpha occupations
        for alpha_set in combinations(range(M), n_alpha):
            alpha_mask = np.zeros(M, dtype=int)
            for p in alpha_set:
                alpha_mask[p] = 1
            for beta_set in combinations(range(M), n_beta):
                beta_mask = np.zeros(M, dtype=int)
                for p in beta_set:
                    beta_mask[p] = 1

                # P_i sign per sub-block
                P_signs = []
                for sb in active_sb:
                    anti = sb_to_antisym[sb]
                    n_anti_occ = sum(int(alpha_mask[p]) + int(beta_mask[p])
                                     for p in anti)
                    P_signs.append(+1 if n_anti_occ % 2 == 0 else -1)
                P_key = tuple(P_signs)

                key = (n_alpha, n_beta, P_key)
                sector_counts[key] = sector_counts.get(key, 0) + 1

    return {
        'sector_counts': {str(k): v for k, v in sector_counts.items()},
        'total_SD': sum(sector_counts.values()),
        'n_active_sub_blocks': n_P,
        'active_sub_blocks': [str(sb) for sb in active_sb],
    }


def main():
    spec = lih_spec()
    result = build_composed_hamiltonian(spec, verbose=False)
    h1 = result['h1']
    M = h1.shape[0]
    Q = 2 * M

    orbital_table = _enumerate_orbitals(spec)
    _, parity = build_pm_rotation(orbital_table)

    n_electrons = 4  # LiH

    print(f"LiH composed: M = {M}, Q = {Q}, n_electrons = {n_electrons}")
    print(f"Full Hilbert: 2^{Q} = {2 ** Q}")

    # Particle-number sector
    sector_data = count_sector_dimensions(M, parity, orbital_table, n_electrons)
    n_SD_total = sector_data['total_SD']
    n_active_sb = sector_data['n_active_sub_blocks']
    print(f"Particle-number sector (N={n_electrons}): {n_SD_total} Slater determinants")
    print(f"Active sub-blocks (Hopf-Z2): {n_active_sb}")

    # Singlet (N_alpha = N_beta = 2 for LiH)
    n_alpha_singlet = n_electrons // 2
    n_beta_singlet = n_electrons - n_alpha_singlet
    singlet_count = sum(
        v for k, v in sector_data['sector_counts'].items()
        if eval(k)[0] == n_alpha_singlet and eval(k)[1] == n_beta_singlet
    )
    print(f"S_z = 0 sector (N_alpha=N_beta={n_alpha_singlet}): {singlet_count} SDs")

    # Joint with all P_i = +1 sector (typical ground state)
    target_P = tuple([1] * n_active_sb)
    target_count = sector_data['sector_counts'].get(
        str((n_alpha_singlet, n_beta_singlet, target_P)), 0
    )
    print(f"S_z = 0, all P = +1 sector: {target_count} SDs")
    print(f"Reduction factor (full / target): {n_SD_total / max(target_count, 1):.1f}x")

    # Print per-(N_alpha, N_beta) breakdown
    print(f"\n  Joint (N_alpha, N_beta, P_signs) sector counts:")
    items = sorted(sector_data['sector_counts'].items(), key=lambda x: -x[1])
    for k, v in items[:20]:
        print(f"    {k}: {v}")

    # Theoretical savings comparison
    print(f"\n=== RESOURCE COMPARISON ===")
    print(f"  Full Slater determinant space (N={n_electrons}): {n_SD_total}")
    print(f"  + Spin parity (S_z = 0 sub-sector):              {singlet_count}")
    print(f"  + Hopf-Z2 per-sub-block ({n_active_sb} Z2 sectors): {2 ** n_active_sb}x finer")
    print(f"  Target sector (S_z=0, ground-state Hopf):        {target_count}")
    print(f"  Singlet S^2 = 0 sub-subspace (CSF count):")
    print(f"    typical = singlet SDs / spin-degeneracy ~ {singlet_count // 3} (rough)")

    out_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'data', 'qpt_joint_block_dim.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump({
            'M': M, 'Q': Q,
            'n_electrons': n_electrons,
            'full_hilbert': 2 ** Q,
            'particle_number_sector_SDs': n_SD_total,
            'sz_zero_sector_SDs': singlet_count,
            'all_plus_hopf_sz_zero': target_count,
            'n_active_sub_blocks': n_active_sb,
            'sector_counts_topk': {k: v for k, v in items[:30]},
        }, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
