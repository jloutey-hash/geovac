"""Sprint F2 Step 1b — verify 3D quadrature convergence for off-diagonal
Na 3s -> Na 4s matrix element.

Step 1's 1% relative error on the offdiag matrix element looked suspicious;
verify by sweeping n_radial and r_max.
"""

from __future__ import annotations

import json

import numpy as np

from geovac.composed_qubit import _radial_wf_grid
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from sprint_f2_step1_kernel_diagnostic import (
    cross_vne_3d_quadrature_same_center,
    cross_vne_converged_multipole,
)


def main():
    R_AB = 3.566

    def R_Na3s(r): return _radial_wf_grid(1.0, 3, 0, r)
    def R_Na4s(r): return _radial_wf_grid(1.0, 4, 0, r)
    def R_H1s(r): return _radial_wf_grid(1.0, 1, 0, r)

    # Reference multipole at L_max=30 (super-high) to be sure
    mp_ref_B = cross_vne_converged_multipole(
        Z_orb=1.0,
        n1=3, l1=0, m1=0,
        n2=4, l2=0, m2=0,
        Z_C=1.0, R_AB=R_AB, L_max_high=30,
    )

    print(f"Reference multipole (L_max=30): {mp_ref_B:+.12f} Ha")
    print()
    print("3D quadrature convergence: <Na 3s | V_ne(H) | Na 4s>")
    print(f"{'n_radial':<10} {'r_max':<8} {'integral_Ha':<20} {'diff_from_mp':<15}")

    rows = []
    for n_radial in [100, 200, 400, 800, 1600, 3200]:
        for r_max in [40.0, 80.0, 120.0]:
            v_3d = cross_vne_3d_quadrature_same_center(
                R_bra=R_Na3s, l_bra=0, m_bra=0,
                R_ket=R_Na4s, l_ket=0, m_ket=0,
                R_C=R_AB, Z_C=1.0,
                n_radial=n_radial, r_max=r_max,
            )
            diff = v_3d - mp_ref_B
            print(f"{n_radial:<10} {r_max:<8} {v_3d:+.12f}    {diff:+.4e}")
            rows.append({
                'n_radial': n_radial,
                'r_max': r_max,
                'val_Ha': v_3d,
                'diff_from_mp_Ha': diff,
            })

    out = {
        'sprint': 'F2 Step 1b — 3D quadrature convergence check',
        'R_AB_bohr': R_AB,
        'multipole_reference_Lmax30_Ha': mp_ref_B,
        'convergence_scan': rows,
    }
    with open('debug/data/sprint_f2_step1b_convergence.json', 'w') as f:
        json.dump(out, f, indent=2)
    print()
    print('Wrote: debug/data/sprint_f2_step1b_convergence.json')


if __name__ == '__main__':
    main()
