"""Sprint F6 Step 1 — feasibility test for NaH max_n=4 with full F3 stack.

Examines:
  - NaH max_n=4 basis structure (M, Q, orbital labels, sub-block layout)
  - Multi-zeta registry coverage at max_n=4 (which orbitals get substituted)
  - FCI dimension at 2-electron singlet sector
  - Single-point build timing estimate (cross-V_ne + multi-zeta + cross-block h1)

Gate decision (per sprint plan):
  - FCI dim < 10k AND build < 1h per single-point -> proceed to Step 2
  - FCI dim > 50k OR build > 3h -> STOP, infeasibility
  - Intermediate -> proceed with caveats
"""

from __future__ import annotations

import json
import sys
import time
from itertools import combinations
from pathlib import Path
from math import comb

import numpy as np

sys.stdout.reconfigure(line_buffering=True)

from geovac.molecular_spec import nah_spec
from geovac.multi_zeta_orbitals import get_physical_valence_orbitals


def main():
    out_path = Path("debug/data/sprint_f6_step1_feasibility.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Sprint F6 Step 1 — NaH max_n=4 feasibility")
    print("=" * 70)

    # --- (a) Spec structure ---
    spec = nah_spec(max_n=4)
    print(f"\nNaH spec at max_n=4:")
    print(f"  name = {spec.name}")
    print(f"  num blocks = {len(spec.blocks)}")
    for i, b in enumerate(spec.blocks):
        print(f"  Block[{i}]: label={b.label}, type={b.block_type}, "
              f"max_n={b.max_n}, Z_center={b.Z_center}, n_electrons={b.n_electrons}, "
              f"has_h_partner={b.has_h_partner}, max_n_partner={b.max_n_partner}, "
              f"n_val_offset={b.n_val_offset}, Z_nuc_center={b.Z_nuc_center}")

    # --- (b) Count orbitals (M) without doing full build ---
    # For NaH (a single bond block with has_h_partner), the orbital structure:
    #   Na side: sum over block_n in [1..max_n] of (orbital count at each block_n)
    #     block_n=1: l=0 (1 orbital)
    #     block_n=2: l in {0, 1} = 1 + 3 = 4 orbitals
    #     block_n=3: l in {0, 1, 2} = 1 + 3 + 5 = 9 orbitals
    #     block_n=4: l in {0, 1, 2, 3} = 1 + 3 + 5 + 7 = 16 orbitals
    #   Total Na side at max_n=4: 1 + 4 + 9 + 16 = 30 orbitals
    #   Same for H side: 30 orbitals
    #   Total M = 60 spatial orbitals, Q = 120 spinors

    def count_orbitals_per_side(max_n):
        # Hydrogenic enumeration: for each block_n in 1..max_n, l ranges 0..(block_n-1)
        # with (2l+1) m-values per l
        total = 0
        per_n = []
        for n in range(1, max_n + 1):
            n_orbs = sum(2 * l + 1 for l in range(n))
            per_n.append((n, n_orbs))
            total += n_orbs
        return total, per_n

    M_Na, per_n_Na = count_orbitals_per_side(4)
    M_H, per_n_H = count_orbitals_per_side(4)
    M_total = M_Na + M_H
    Q_total = 2 * M_total

    print(f"\nOrbital counts:")
    print(f"  Na side per block_n: {per_n_Na}")
    print(f"  H side per block_n:  {per_n_H}")
    print(f"  M (spatial) = {M_total}, Q (spinors) = {Q_total}")

    # --- (c) Multi-zeta registry coverage ---
    physical_orbitals = get_physical_valence_orbitals(11)  # Na
    print(f"\nMulti-zeta registry for Na (Z=11):")
    print(f"  Number of entries: {len(physical_orbitals)}")
    for orb in physical_orbitals:
        print(f"    n_orbital={orb.n_orbital}, l_orbital={orb.l_orbital}")

    # Dispatch: physical_n = block_n + n_val_offset = block_n + 2
    # So block_n=1 -> physical n=3, block_n=2 -> physical n=4, etc.
    print(f"\nMulti-zeta dispatch at max_n=4 (n_val_offset=2):")
    coverage_table = []
    for block_n in range(1, 5):
        for l in range(block_n):
            physical_n = block_n + 2
            in_registry = any(
                (orb.n_orbital == physical_n and orb.l_orbital == l)
                for orb in physical_orbitals
            )
            label = f"Na {physical_n}{'spdf'[l]}"
            coverage_table.append({
                "block_n": block_n,
                "l": l,
                "physical_n": physical_n,
                "label": label,
                "in_registry": in_registry,
            })
            mark = "[YES] SUBSTITUTED" if in_registry else "      hydrogenic placeholder"
            print(f"  block_n={block_n}, l={l}: {label}  ({mark})")

    # --- (d) FCI dim for n_electrons=2 singlet ---
    # 2-electron singlet sector: n_alpha = n_beta = 1
    # FCI dim = C(M, 1) * C(M, 1) = M^2
    n_e = 2
    n_alpha = n_e // 2
    n_beta = n_e // 2
    fci_dim = comb(M_total, n_alpha) * comb(M_total, n_beta)
    print(f"\n2-electron singlet FCI dim at M={M_total}:")
    print(f"  C({M_total},{n_alpha}) * C({M_total},{n_beta}) = {fci_dim}")

    # Compare to previous sprints
    fci_dim_max_n_2 = comb(10, 1) ** 2  # M=10 at max_n=2
    fci_dim_max_n_3 = comb(28, 1) ** 2  # M=28 at max_n=3
    print(f"  (for context: max_n=2 had fci_dim={fci_dim_max_n_2}; "
          f"max_n=3 had fci_dim={fci_dim_max_n_3})")

    # --- (e) Multipole content at max_n=4 ---
    # Cross-V_ne L_max = 2 * l_max = 2 * 3 = 6 at max_n=4 (vs L_max=4 at max_n=3, L_max=2 at max_n=2)
    l_max_n2 = 1
    l_max_n3 = 2
    l_max_n4 = 3
    L_max_xvne_n2 = 2 * l_max_n2
    L_max_xvne_n3 = 2 * l_max_n3
    L_max_xvne_n4 = 2 * l_max_n4
    print(f"\nCross-V_ne multipole content:")
    print(f"  max_n=2: l_max={l_max_n2}, L_max={L_max_xvne_n2}")
    print(f"  max_n=3: l_max={l_max_n3}, L_max={L_max_xvne_n3}")
    print(f"  max_n=4: l_max={l_max_n4}, L_max={L_max_xvne_n4}")

    # --- (f) Cross-block h1 architectural scope ---
    # cross_block_h1.py module is s-s only currently!
    # At max_n=4, the s-s pairs are (1s,1s), (1s,2s), (1s,3s), (1s,4s) on each side combined
    # cross-block s-s pairs: 4 (Na) x 4 (H) = 16 pairs
    print(f"\nCross-block h1 (s-s only) scope at max_n=4:")
    s_orbs_Na = sum(1 for (n, l) in [(n, l) for n in range(1, 5) for l in range(n)] if l == 0)
    s_orbs_H = s_orbs_Na
    n_xb_h1_ss_pairs = s_orbs_Na * s_orbs_H
    print(f"  s-orbitals per side: {s_orbs_Na}")
    print(f"  cross-block s-s pairs: {n_xb_h1_ss_pairs}")
    print(f"  WARNING: l>0 cross-block h1 pairs raise NotImplementedError")
    print(f"  -> Architecture handles {n_xb_h1_ss_pairs} of total possible {(M_total//2)**2} cross-block pairs")
    print(f"  -> s-only coverage = {n_xb_h1_ss_pairs / (M_total//2)**2 * 100:.1f}%")

    # --- (g) Build-time estimate ---
    # max_n=3 took ~10 s per build (per F1 max_n=3 memo §1.3).
    # max_n=4 has 60/28 = 2.14x more orbitals.
    # h1 build scales O(M^2) per matrix element -> 4.6x more elements
    # ERI build scales O(M^4) -> 21x more elements
    # cross-V_ne scales O(M^2) x O(L_max+1) -> 4.6x x 1.4x = 6.5x
    # multipole content: L_max=6 vs L_max=4 -> ~1.4x more multipole terms
    # So roughly 21x ERI cost dominates -> ~210 s ~= 3.5 min per build
    M_ratio = M_total / 28
    eri_ratio_estimate = M_ratio ** 4
    xvne_ratio_estimate = M_ratio ** 2 * (L_max_xvne_n4 + 1) / (L_max_xvne_n3 + 1)
    print(f"\nBuild-time estimate (from max_n=3 baseline ~10 s):")
    print(f"  M_ratio = {M_ratio:.2f}")
    print(f"  ERI scaling estimate (~M^4): {eri_ratio_estimate:.1f}x -> {eri_ratio_estimate * 10:.0f}s = {eri_ratio_estimate * 10 / 60:.1f}min")
    print(f"  cross-V_ne scaling estimate: {xvne_ratio_estimate:.1f}x -> {xvne_ratio_estimate * 10:.0f}s = {xvne_ratio_estimate * 10 / 60:.1f}min")
    print(f"  Combined (dominant ERI): ~{eri_ratio_estimate * 10:.0f}s = {eri_ratio_estimate * 10 / 60:.1f}min per single-point build")

    # FCI assembly: O(fci_dim^2 * M^2) for matrix construction, but with sparsity ~M^2 per row
    fci_assembly_estimate = (fci_dim / 100) ** 2 * (M_total / 10) ** 2 * 0.1
    print(f"  FCI assembly estimate: ~{fci_assembly_estimate:.0f}s = {fci_assembly_estimate / 60:.1f}min")

    total_estimate_min = (eri_ratio_estimate * 10 + fci_assembly_estimate) / 60
    print(f"  TOTAL estimate per single-point: ~{total_estimate_min:.1f}min")
    print(f"  For 2 single-points (R=3.5, R=10): ~{2 * total_estimate_min:.1f}min")

    # --- (h) Gate decision ---
    print(f"\n{'=' * 70}")
    print(f"GATE DECISION:")
    print(f"  FCI dim = {fci_dim} (< 10k: {fci_dim < 10000}, < 50k: {fci_dim < 50000})")
    print(f"  Build estimate per single-point ~= {total_estimate_min:.1f}min")

    if fci_dim < 10_000 and total_estimate_min < 60:
        gate = "PROCEED_TO_STEP_2"
        gate_reason = "Both FCI dim < 10k and build < 1h satisfied"
    elif fci_dim > 50_000 or total_estimate_min > 180:
        gate = "STOP_INFEASIBILITY"
        gate_reason = "FCI dim > 50k or build > 3h"
    else:
        gate = "PROCEED_WITH_CAVEATS"
        gate_reason = "Intermediate regime"

    print(f"  -> Gate verdict: {gate}")
    print(f"  -> Reason: {gate_reason}")
    print(f"{'=' * 70}")

    out = {
        "sprint": "F6 Step 1 — NaH max_n=4 feasibility",
        "date": "2026-05-23",
        "spec": {
            "name": spec.name,
            "num_blocks": len(spec.blocks),
            "blocks": [
                {
                    "label": b.label,
                    "type": b.block_type,
                    "max_n": b.max_n,
                    "Z_center": b.Z_center,
                    "n_electrons": b.n_electrons,
                    "has_h_partner": b.has_h_partner,
                    "max_n_partner": b.max_n_partner,
                    "n_val_offset": b.n_val_offset,
                    "Z_nuc_center": b.Z_nuc_center,
                }
                for b in spec.blocks
            ],
        },
        "orbital_counts": {
            "M_Na": M_Na,
            "M_H": M_H,
            "M_total": M_total,
            "Q_total": Q_total,
            "per_n_Na": [{"block_n": n, "n_orbs": k} for n, k in per_n_Na],
        },
        "multi_zeta_registry": {
            "Z": 11,
            "num_entries": len(physical_orbitals),
            "entries": [
                {"n_orbital": orb.n_orbital, "l_orbital": orb.l_orbital}
                for orb in physical_orbitals
            ],
            "coverage_at_max_n_4": coverage_table,
        },
        "fci": {
            "n_electrons": n_e,
            "M_total": M_total,
            "fci_dim_singlet": fci_dim,
            "context_max_n_2": fci_dim_max_n_2,
            "context_max_n_3": fci_dim_max_n_3,
        },
        "multipole": {
            "l_max_max_n_2": l_max_n2,
            "l_max_max_n_3": l_max_n3,
            "l_max_max_n_4": l_max_n4,
            "L_max_xvne_max_n_4": L_max_xvne_n4,
        },
        "cross_block_h1_scope": {
            "s_orbs_per_side": s_orbs_Na,
            "n_ss_pairs": n_xb_h1_ss_pairs,
            "total_xb_pairs": (M_total // 2) ** 2,
            "s_only_coverage_pct": n_xb_h1_ss_pairs / (M_total // 2) ** 2 * 100,
            "l_gt_0_pairs_raise": True,
        },
        "build_time_estimate_min": total_estimate_min,
        "gate": gate,
        "gate_reason": gate_reason,
    }
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote: {out_path}")


if __name__ == "__main__":
    main()
