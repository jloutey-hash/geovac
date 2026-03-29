"""
Fast Sections 4-5: minimal R points for R_eq comparison.
Section 3 result: Z_eff(1.81) = 6.0000. Z_eff injection is not the issue.
"""

import json
import time
import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.level4_multichannel import solve_level4_h2_multichannel
from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK

L_MAX = 2
N_ALPHA = 100
N_RE = 300


def main():
    t0 = time.time()
    results = {}

    # ===== SECTION 5: PK effect (most informative) =====
    print("=" * 60)
    print("SECTION 5: PK effect on Z_A=6, Z_B=1 bond pair")
    print("=" * 60)

    core = CoreScreening(Z=8)
    core.solve()
    pk = AbInitioPK(core, n_core=2)
    print(f"PK: A={pk.A:.4f}, B={pk.B:.4f}, r_core={pk.r_core:.4f}")
    pk_potentials = [pk.pk_dict(atom='A')]

    # Coarse scan for R_eq: 0.5 to 4.0 in steps of 0.2
    R_scan = np.arange(0.5, 4.1, 0.2)

    E_nopk = []
    E_pk = []

    print(f"\n{'R':>6s}  {'E_nopk':>12s}  {'E_pk':>12s}  {'dE_pk':>10s}")
    print(f"{'-'*6}  {'-'*12}  {'-'*12}  {'-'*10}")

    for R in R_scan:
        sys.stdout.flush()
        res_n = solve_level4_h2_multichannel(
            R, l_max=L_MAX, Z_A=6.0, Z_B=1.0,
            n_alpha=N_ALPHA, n_Re=N_RE, verbose=False,
            origin='charge_center'
        )
        res_p = solve_level4_h2_multichannel(
            R, l_max=L_MAX, Z_A=6.0, Z_B=1.0,
            n_alpha=N_ALPHA, n_Re=N_RE, verbose=False,
            origin='charge_center', pk_potentials=pk_potentials
        )
        E_n = res_n['E_total']
        E_p = res_p['E_total']
        dE = E_p - E_n
        E_nopk.append(E_n)
        E_pk.append(E_p)
        print(f"{R:6.2f}  {E_n:12.6f}  {E_p:12.6f}  {dE:+10.6f}")
        sys.stdout.flush()

    idx_n = np.argmin(E_nopk)
    idx_p = np.argmin(E_pk)
    R_eq_nopk = float(R_scan[idx_n])
    R_eq_pk = float(R_scan[idx_p])

    print(f"\nR_eq (no PK):  {R_eq_nopk:.1f} bohr")
    print(f"R_eq (PK):     {R_eq_pk:.1f} bohr")
    print(f"R_eq shift:    {R_eq_pk - R_eq_nopk:+.1f} bohr")

    results['section5_pk'] = {
        'pk_A': pk.A, 'pk_B': pk.B, 'pk_r_core': pk.r_core,
        'R_eq_nopk': R_eq_nopk, 'R_eq_pk': R_eq_pk,
        'R_eq_shift': R_eq_pk - R_eq_nopk,
    }

    t_elapsed = time.time() - t0
    print(f"\nTotal time: {t_elapsed:.0f}s")

    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'data', 'h2o_bond_pk_effect.json')
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Saved to {outpath}")


if __name__ == '__main__':
    main()
