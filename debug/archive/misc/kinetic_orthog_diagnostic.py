"""
Diagnostic: Löwdin kinetic orthogonalization correction for BeH2.

Computes DeltaT(R) at 5 R-points and compares to full exchange energy.
Also runs full PES scans with full_exchange and full_exchange_kinetic
to compare R_eq.
"""

import json
import os
import numpy as np

from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK
from geovac.level4_multichannel import solve_level4_h2_multichannel
from geovac.inter_fiber_coupling import (
    full_exchange_inter_fiber_energy,
    kinetic_orthogonalization_energy,
)
from geovac.composed_triatomic import ComposedTriatomicSolver


def main():
    print("=" * 64)
    print("Kinetic Orthogonalization Diagnostic — BeH2")
    print("=" * 64)

    # Setup
    print("\nSetting up core + PK...")
    core = CoreScreening(Z=4, l_max=2, n_alpha=200)
    core.solve(verbose=False)
    pk = AbInitioPK(core, n_core=2)
    pk_potentials = [pk.pk_dict(atom='A')]

    Z_eff, Z_B, l_max, n_alpha = 2.0, 1.0, 2, 100

    # Part 1: Point-by-point diagnostic
    print("\n--- Point-by-point DeltaT(R) diagnostic ---")
    print(f"  {'R':>5s}  {'DeltaT':>10s}  {'E_exch':>10s}  "
          f"{'|DT/Eex|':>8s}  {'||S_AB||':>8s}  {'Re2_avg':>8s}")
    print(f"  {'-'*5}  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*8}  {'-'*8}")

    diag_data = {}
    R_points = [2.0, 2.5, 3.0, 4.0, 5.0]

    for R in R_points:
        l4 = solve_level4_h2_multichannel(
            R=R, Z_A=Z_eff, Z_B=Z_B,
            l_max=l_max, n_alpha=n_alpha, n_Re=300,
            verbose=False, pk_potentials=pk_potentials,
        )

        # Full exchange
        exch = full_exchange_inter_fiber_energy(
            l4, R, Z_A=Z_eff, Z_B=Z_B,
            l_max=l_max, n_alpha=n_alpha,
            pk_potentials=pk_potentials,
            n_sample_Re=10,
        )
        E_exch = exch['E_exchange']

        # Kinetic correction
        kin = kinetic_orthogonalization_energy(
            l4, R, Z_A=Z_eff, Z_B=Z_B,
            l_max=l_max, n_alpha=n_alpha,
            pk_potentials=pk_potentials,
            n_sample_Re=10,
        )
        delta_T = kin['delta_T']
        ratio = abs(delta_T / E_exch) if abs(E_exch) > 1e-15 else 0.0

        print(f"  {R:5.1f}  {delta_T:10.6f}  {E_exch:10.6f}  "
              f"{ratio:8.4f}  {kin['S_AB_frobenius']:8.4f}  "
              f"{kin['Re2_avg']:8.4f}")

        diag_data[str(R)] = {
            'delta_T': delta_T,
            'E_exchange': E_exch,
            'ratio_abs': ratio,
            'S_AB_frobenius': kin['S_AB_frobenius'],
            'Re2_avg': kin['Re2_avg'],
            'T_l': kin['T_l'].tolist(),
            'l_values': kin['l_values'].tolist(),
            'E_val': kin['E_val'],
            'S_AB_matrix': kin['S_AB_matrix'].tolist(),
        }

    # Part 2: PES comparison
    print("\n--- PES comparison: full_exchange vs full_exchange_kinetic ---")

    R_ref = 2.507  # bohr

    for mode in ['full_exchange', 'full_exchange_kinetic']:
        print(f"\n  Mode: {mode}")
        solver = ComposedTriatomicSolver.BeH2(
            l_max=2, interbond_mode=mode, verbose=False)
        solver.solve_core()
        R_grid = np.arange(1.8, 4.5, 0.2)
        pes = solver.scan_pes(R_grid=R_grid, n_Re=300)

        R_eq = pes['R_eq']
        D_e = pes['D_e']
        err_pct = abs(R_eq - R_ref) / R_ref * 100.0

        print(f"  R_eq = {R_eq:.3f} bohr (ref {R_ref}), "
              f"error = {err_pct:.1f}%")
        print(f"  D_e  = {D_e:.6f} Ha")
        print(f"  E_min = {pes['E_min']:.6f} Ha")

        diag_data[f'pes_{mode}'] = {
            'R_eq': R_eq,
            'D_e': D_e,
            'E_min': pes['E_min'],
            'R_eq_error_pct': err_pct,
        }

    # Save diagnostic
    os.makedirs('debug/data', exist_ok=True)
    with open('debug/data/kinetic_orthog_diagnostic.json', 'w') as f:
        json.dump(diag_data, f, indent=2)
    print("\nDiagnostic saved to debug/data/kinetic_orthog_diagnostic.json")


if __name__ == '__main__':
    main()
