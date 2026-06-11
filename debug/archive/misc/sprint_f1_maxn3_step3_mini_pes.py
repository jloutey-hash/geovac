"""Sprint F1 max_n=3 Step 3: mini-PES scan.

Per the Step 2 decision gate: PROCEED_TO_MINI_PES because D_e^combined > 0
even though natural occupations are still (1, 1) — substantive new content
is the dominant natural orbital shows a true bonding combination
(Na 3s ± H 1s, 50/50 amplitude split).

This sprint runs a 4-point PES at R ∈ {3.0, 3.5, 4.0, 5.0} bohr to
locate R_eq if there's an internal minimum, plus dissociation reference
at R = 10.0 bohr. Optional extension to R ∈ {2.0, 2.5, 7.0} to test
the over-attraction regime if needed.

Verifies:
  - Internal minimum location vs P1 range [3.0, 4.5] bohr
  - Concavity near R_min (genuine internal minimum, not flat)
  - D_e vs P2 range [0.0375, 0.150] Ha
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(line_buffering=True)

# Import the FCI machinery from Step 2
sys.path.insert(0, str(Path('debug').resolve()))
from sprint_f1_maxn3_step2_fci import run_arch

from geovac.molecular_spec import nah_spec


def main():
    out_path = Path('debug/data/sprint_f1_maxn3_step3_pes.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Sprint F1 max_n=3 Step 3: Mini-PES at NaH max_n=3")
    print("=" * 70)
    print()

    spec = nah_spec(max_n=3)

    # Two phases: first the canonical 4-point + diss reference; then if
    # an internal minimum appears we ALSO probe shorter R to verify the
    # PES doesn't keep descending below R_min.
    canonical_R = [3.0, 3.5, 4.0, 5.0, 10.0]
    extension_R = [2.0, 2.5, 7.0]  # extra points for shape verification

    print(f"Canonical R-points: {canonical_R}")
    print(f"Extension R-points: {extension_R}")
    print()

    # Run combined architecture only (the headline)
    print("--- COMBINED W1c + multi-zeta ---")
    pes_combined = {}
    for R in canonical_R + extension_R:
        r = run_arch(spec, R, screened=True, multi_zeta=True,
                     label='combined_w1c_mz', verbose=False)
        pes_combined[R] = r
        print(f"  R={R:6.3f}  E={r['E_gs']:+.6f}  top_no={r['natural_occupations_top_8'][0]:.4f}  "
              f"dom_NO Na/H = {r['dominant_no_na_amp2_total']:.3f}/{r['dominant_no_h_amp2_total']:.3f}  "
              f"[{r['total_time_s']:.0f}s]", flush=True)
    print()

    # Also run W1c-alone for the differential
    print("--- W1c alone (no multi-zeta) ---")
    pes_w1c_alone = {}
    for R in canonical_R + extension_R:
        r = run_arch(spec, R, screened=True, multi_zeta=False,
                     label='w1c_alone', verbose=False)
        pes_w1c_alone[R] = r
        print(f"  R={R:6.3f}  E={r['E_gs']:+.6f}  top_no={r['natural_occupations_top_8'][0]:.4f}  "
              f"[{r['total_time_s']:.0f}s]", flush=True)
    print()

    # Find R_min for combined
    R_sorted = sorted(pes_combined.keys())
    E_combined = [pes_combined[R]['E_gs'] for R in R_sorted]
    R_min_idx = int(np.argmin(E_combined))
    R_min = R_sorted[R_min_idx]
    E_min = E_combined[R_min_idx]
    E_diss = pes_combined[10.0]['E_gs']
    D_e_pes = E_diss - E_min

    print("--- PES summary (COMBINED) ---")
    print(f"  {'R':>6}  {'E_combined':>13}  {'E_w1c_alone':>13}  {'mz_diff':>10}")
    for R in R_sorted:
        E_c = pes_combined[R]['E_gs']
        E_w = pes_w1c_alone[R]['E_gs']
        diff = E_w - E_c
        marker = '  <-- R_min' if R == R_min else ''
        print(f"  {R:6.3f}  {E_c:+13.6f}  {E_w:+13.6f}  {diff:+10.6f}{marker}")
    print()

    # PES verification
    print(f"  R_min (smallest tested in canonical): {R_min} bohr")
    print(f"  E_min: {E_min:+.6f} Ha")
    print(f"  D_e (PES): {D_e_pes:+.6f} Ha")
    print(f"  Smallest tested R: {R_sorted[0]}")
    print()

    # Internal minimum check
    in_canonical = R_min in canonical_R
    is_internal_min = (R_min not in [min(R_sorted), max(R_sorted)] and R_min < 10.0)
    is_extension_min = R_min in extension_R
    smallest_R_min = R_min == min(R_sorted)

    print(f"  Internal minimum (not at boundary)? {'YES' if is_internal_min else 'NO'}")
    print(f"  Min at smallest tested R (over-attraction)? "
          f"{'YES' if smallest_R_min else 'NO'}")
    print()

    # Concavity check around R_min
    if 0 < R_min_idx < len(R_sorted) - 1:
        R_left = R_sorted[R_min_idx - 1]
        R_right = R_sorted[R_min_idx + 1]
        E_left = E_combined[R_min_idx - 1]
        E_right = E_combined[R_min_idx + 1]
        is_concave = (E_left > E_min) and (E_right > E_min)
        print(f"  Concavity check: E({R_left})={E_left:+.4f}, E({R_min})={E_min:+.4f}, "
              f"E({R_right})={E_right:+.4f}")
        print(f"  Concave (E_min < both neighbors)? {'YES' if is_concave else 'NO'}")
    else:
        is_concave = False
        print(f"  Min at boundary, cannot check concavity")
    print()

    # Predictions verification
    print("=" * 70)
    print("PREDICTIONS VERIFICATION")
    print("=" * 70)
    print()

    # P1: internal PES minimum at R_eq in [3.0, 4.5] bohr
    p1_pass = is_internal_min and (3.0 <= R_min <= 4.5)
    print(f"P1: Internal PES minimum in R in [3.0, 4.5] bohr")
    print(f"    Actual R_min = {R_min} bohr")
    print(f"    Falsifier: R_min at smallest tested R (e.g. 2.0)")
    print(f"    Result: {'PASS' if p1_pass else 'FAIL'}")
    print()

    # P2: D_e in [0.0375, 0.150] Ha
    p2_pass = 0.0375 <= D_e_pes <= 0.150 if p1_pass else False
    print(f"P2: D_e in [0.0375, 0.150] Ha (within 2x of experimental)")
    print(f"    Actual D_e = {D_e_pes:+.6f} Ha")
    print(f"    Result: {'PASS' if p2_pass else 'FAIL (conditional on P1)'}")
    print()

    # P3: |mz_diff| ∈ [0.02, 0.30] Ha at R_eq
    mz_diff_at_R_min = pes_w1c_alone[R_min]['E_gs'] - pes_combined[R_min]['E_gs']
    p3_pass = 0.02 <= abs(mz_diff_at_R_min) <= 0.30
    print(f"P3: |E_W1c - E_W1c+mz| in [0.02, 0.30] Ha at R_eq")
    print(f"    Actual at R={R_min}: |{abs(mz_diff_at_R_min):.4f}| Ha")
    print(f"    Result: {'PASS' if p3_pass else 'FAIL'}")
    print()

    # Verdict
    if p1_pass and p2_pass and p3_pass:
        verdict = "THREE-BUCKET-PARTITION-VERIFIED"
    elif p1_pass and p3_pass:
        verdict = "PARTIAL (P1 + P3 PASS, P2 outside range)"
    elif p1_pass:
        verdict = "PARTIAL (P1 PASS only)"
    elif not p1_pass and smallest_R_min:
        verdict = "CLEAN NEGATIVE (P1 FAIL — over-attraction wall)"
    else:
        verdict = "MIXED"

    print("=" * 70)
    print(f"STRUCTURAL VERDICT: {verdict}")
    print("=" * 70)

    # Save
    results = {
        'sprint': 'F1 max_n=3 Step 3',
        'date': '2026-05-23',
        'canonical_R': canonical_R,
        'extension_R': extension_R,
        'all_R_sorted': R_sorted,
        'E_combined_by_R': {f'{R:.3f}': pes_combined[R]['E_gs'] for R in R_sorted},
        'E_w1c_alone_by_R': {f'{R:.3f}': pes_w1c_alone[R]['E_gs'] for R in R_sorted},
        'top_no_combined_by_R': {f'{R:.3f}': pes_combined[R]['natural_occupations_top_8'][0] for R in R_sorted},
        'dom_no_NaH_split_combined': {
            f'{R:.3f}': [pes_combined[R]['dominant_no_na_amp2_total'],
                         pes_combined[R]['dominant_no_h_amp2_total']] for R in R_sorted
        },
        'R_min': R_min,
        'E_min': E_min,
        'E_diss_R_10': pes_combined[10.0]['E_gs'],
        'D_e_pes_Ha': D_e_pes,
        'is_internal_min': is_internal_min,
        'is_concave_at_min': is_concave,
        'smallest_tested_R': min(R_sorted),
        'predictions': {
            'P1_internal_min_in_3_to_4p5': {
                'predicted_range_bohr': [3.0, 4.5],
                'actual_R_min': R_min,
                'PASS': p1_pass,
            },
            'P2_D_e_in_0p0375_to_0p150': {
                'predicted_range_Ha': [0.0375, 0.150],
                'actual_D_e_Ha': D_e_pes,
                'PASS': p2_pass,
            },
            'P3_mz_diff_in_0p02_to_0p30': {
                'predicted_range_Ha': [0.02, 0.30],
                'actual_mz_diff_Ha': abs(mz_diff_at_R_min),
                'PASS': p3_pass,
            },
        },
        'verdict': verdict,
    }

    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == '__main__':
    main()
