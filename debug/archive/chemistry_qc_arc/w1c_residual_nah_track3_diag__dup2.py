"""
Diagnostic experiment for Track 3: estimate the cross-V_ne reduction
that would result from replacing the hydrogenic Z=1 1s on Na with the
physical-n hydrogenic 3s shape (still hydrogenic, but at the actual
principal quantum number).

This is a *one-shot diagnostic*, NOT a production change to the
balanced-coupled architecture. It tests whether the wavefunction-shape
correction (replacing block_n=1 → physical_n=3 for the Na valence)
moves cross-V_ne in the direction needed to close the wall.
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np

from geovac.shibuya_wulfman import compute_cross_center_vne


def main():
    out_path = Path('debug/data/w1c_residual_nah_track3_diag.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Na valence states under FRAMEWORK convention (block_n labels):
    framework_states = [(1, 0, 0), (2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, 1)]
    # Same orbital labels but with PHYSICAL n labels for [Ne] core valence:
    # block_n=1 (l=0) -> Na 3s (physical n=3)
    # block_n=2 (l=0) -> Na 4s (physical n=4)
    # block_n=2 (l=1) -> Na 3p (physical n=3, lowest p outside core)
    # We test BOTH conventions for block_n=2,l=1: physical 4p (matches framework's
    # n_val_offset = period - 1 rule) and physical 3p (chemically natural).
    physical_states_4p = [(3, 0, 0), (4, 0, 0), (4, 1, -1), (4, 1, 0), (4, 1, 1)]
    physical_states_3p = [(3, 0, 0), (4, 0, 0), (3, 1, -1), (3, 1, 0), (3, 1, 1)]

    R_grid = [2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0]
    Z_nuc_H = 1.0
    Z_orb = 1.0  # asymptotic Na valence Z

    print("=== Diagnostic: cross-V_ne shape correction ===")
    print("Trace of cross-V_ne matrix (Z_orb=1, Z_nuc_H=1) for 3 conventions")
    print()
    print(f"{'R':>6s}  {'framework (block_n)':>20s}  {'physical (4p)':>14s}  "
          f"{'physical (3p)':>14s}  {'frame-vs-3p diff':>16s}")
    print('-' * 80)
    rows = []
    for R in R_grid:
        v_frame = compute_cross_center_vne(
            Z_orb, framework_states, Z_nuc_H, R,
            L_max=4, n_grid=8000, direction=(0, 0, 1),
        )
        v_4p = compute_cross_center_vne(
            Z_orb, physical_states_4p, Z_nuc_H, R,
            L_max=4, n_grid=8000, direction=(0, 0, 1),
        )
        v_3p = compute_cross_center_vne(
            Z_orb, physical_states_3p, Z_nuc_H, R,
            L_max=4, n_grid=8000, direction=(0, 0, 1),
        )
        tf = float(np.trace(v_frame))
        t4 = float(np.trace(v_4p))
        t3 = float(np.trace(v_3p))
        diff = tf - t3
        print(f"{R:6.3f}  {tf:>20.4f}  {t4:>14.4f}  {t3:>14.4f}  {diff:>16.4f}")
        rows.append({
            'R': float(R),
            'framework_trace': tf,
            'physical_4p_trace': t4,
            'physical_3p_trace': t3,
            'framework_minus_physical_3p': diff,
        })

    print()
    # The framework over-attracts most at small R. The "physical" convention
    # would *under-attract* the Na-side electrons at small R relative to
    # framework. Let's check whether the differential shape change closes
    # the W1c-residual gap.
    diff_smallR = rows[0]['framework_minus_physical_3p']
    diff_largeR = rows[-1]['framework_minus_physical_3p']
    descent_correction = diff_smallR - diff_largeR
    print(f"Differential (small R - large R) = {diff_smallR:.4f} - {diff_largeR:.4f} "
          f"= {descent_correction:.4f} Ha")
    print(f"Compare to W1c-residual descent = 0.357 Ha (Track 2 baseline) "
          f"and PK addition = -0.052 Ha")

    out = {
        'description': 'Cross-V_ne diagnostic: compare framework block_n labels '
                      'vs physical n labels for the Na valence orbitals',
        'rows': rows,
        'descent_correction_estimate': float(descent_correction),
    }
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\n[saved] {out_path}")


if __name__ == '__main__':
    main()
