"""Sprint F1 max_n=3 Step 1: algebraic kernel differential.

Computes the cross-V_ne kernel matrix element on the Na 3s diagonal
(orbital index 0 in the NaH spec) at R_eq = 3.566 bohr with both
hydrogenic Z=1 placeholder and physical multi-zeta basis.

Sanity check confirming the multi-zeta machinery works at max_n=3.
Should give similar magnitude to the max_n=2 result (~-0.135 Ha range).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(line_buffering=True)

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import nah_spec


def main():
    out_path = Path('debug/data/sprint_f1_maxn3_step1_kernel.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Sprint F1 max_n=3 Step 1: Cross-V_ne kernel differential at Na 3s")
    print("=" * 70)
    print()

    R = 3.566  # NaH experimental R_eq
    spec = nah_spec(max_n=3)

    results = {
        'sprint': 'F1 max_n=3 Step 1: kernel differential',
        'date': '2026-05-23',
        'R_bohr': R,
        'spec': 'NaH max_n=3',
        'M': None,
        'tests': [],
    }

    # Run 4 architectures: bare/no-mz, bare/mz, W1c/no-mz, W1c/mz
    archs = [
        ('bare_no_mz', False, False),
        ('bare_mz',    False, True),
        ('w1c_no_mz',  True,  False),
        ('w1c_mz',     True,  True),
    ]

    print("Building 4 architectures at R = {:.3f} bohr...".format(R))
    print()

    h1_results = {}
    for label, screened, mz in archs:
        t0 = time.perf_counter()
        ham = build_balanced_hamiltonian(
            spec, R=R, n_grid_vne=4000, L_max=4,
            screened_cross_center=screened,
            multi_zeta_basis=mz,
            verbose=False,
        )
        dt = time.perf_counter() - t0
        M = ham['M']
        h1 = ham['h1']
        # Na 3s is at orbital index 0
        h1_na3s_diag = float(h1[0, 0])
        # h1_cross_vne is the cross-V_ne contribution
        h1_cv = ham['h1_cross_vne']
        h1_cv_na3s = float(h1_cv[0, 0])
        # The cross-V_ne diagonal at Na 3s = ⟨Na 3s | sum_off-nuclei -Z_n/|r-R_n| | Na 3s⟩
        # In NaH the only off-center nucleus seen by Na is H, so this is
        # ⟨Na 3s | -1/|r-R_H| | Na 3s⟩
        # Z_H = 1
        results['tests'].append({
            'label': label,
            'screened': screened,
            'multi_zeta': mz,
            'h1_na3s_diag_Ha': h1_na3s_diag,
            'h1_cross_vne_na3s_diag_Ha': h1_cv_na3s,
            'h1_diag_first_8': h1.diagonal()[:8].tolist(),
            'h1_cross_vne_diag_first_8': h1_cv.diagonal()[:8].tolist(),
            'multi_zeta_diagnostics': ham.get('multi_zeta_diagnostics', []),
            'build_time_s': dt,
        })
        print(f"  {label:12}: h1[0,0]={h1_na3s_diag:+.6f}, "
              f"h1_cv[0,0]={h1_cv_na3s:+.6f} ({dt:.1f}s)")
        h1_results[label] = (h1_na3s_diag, h1_cv_na3s, h1)
        results['M'] = M

    print()
    print("Differentials (mz vs no-mz on the same screening level):")
    print()

    diff_bare = h1_results['bare_mz'][1] - h1_results['bare_no_mz'][1]
    diff_w1c = h1_results['w1c_mz'][1] - h1_results['w1c_no_mz'][1]
    print(f"  bare:  h1_cv_na3s(mz) - h1_cv_na3s(no_mz) = {diff_bare:+.6f} Ha")
    print(f"  W1c:   h1_cv_na3s(mz) - h1_cv_na3s(no_mz) = {diff_w1c:+.6f} Ha")
    print()
    print("Comparison with alpha-PES Step 1 (max_n=2 baseline):")
    print(f"  alpha-PES Step 1 reported -0.135 Ha for bare cross-V_ne mz-vs-no_mz "
          f"diagonal differential at R=3.566.")
    print(f"  At max_n=3 (this sprint), bare: {diff_bare:+.6f} Ha; "
          f"W1c: {diff_w1c:+.6f} Ha.")
    print()

    results['differential'] = {
        'bare_mz_minus_no_mz': diff_bare,
        'w1c_mz_minus_no_mz': diff_w1c,
        'alpha_pes_step1_max_n_2_reference': -0.135,
        'consistent_with_alpha_pes': abs(diff_bare - (-0.135)) < 0.1,
    }

    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Saved to {out_path}")


if __name__ == '__main__':
    main()
