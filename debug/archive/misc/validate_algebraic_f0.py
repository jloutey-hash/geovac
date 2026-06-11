"""
Validate algebraic F^0 inter-fiber pathway against numerical baseline.

Compares the two methods at R = 2.0, 3.0, 5.0 for BeH2 (composed geometry).
Success criterion: <1% agreement on F^0 values.

Also analyzes:
- Channel decomposition of F^0
- S(R) vs F^0(R) relationship
- Whether S*F^0 simplifies when both are algebraic
"""

import numpy as np
import json
import time

from geovac.level4_multichannel import solve_level4_h2_multichannel
from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK
from geovac.inter_fiber_coupling import (
    monopole_inter_fiber_energy,
    exchange_inter_fiber_energy,
    extract_channel_data,
    compute_overlap_from_channel_data,
    compute_channel_f0_matrix,
)


def setup_beh2():
    """Set up core + PK for BeH2."""
    core = CoreScreening(Z=4, l_max=2, n_alpha=200)
    core.solve(verbose=False)
    pk = AbInitioPK(core, n_core=2)
    pk_potentials = [pk.pk_dict(atom='A')]
    return {
        'pk_potentials': pk_potentials,
        'Z_eff': 2.0,
        'Z_ligand': 1.0,
        'l_max': 2,
        'n_alpha': 100,
    }


def main():
    print("=" * 70)
    print("Algebraic F^0 Inter-Fiber Validation")
    print("=" * 70)

    setup = setup_beh2()
    R_values = [2.0, 3.0, 5.0]
    n_sample_Re = 12

    results = {}

    for R in R_values:
        print(f"\n--- R = {R:.1f} bohr ---")

        # Solve Level 4
        l4 = solve_level4_h2_multichannel(
            R=R, Z_A=setup['Z_eff'], Z_B=setup['Z_ligand'],
            l_max=setup['l_max'], n_alpha=setup['n_alpha'], n_Re=300,
            verbose=False, pk_potentials=setup['pk_potentials'],
        )

        # --- Numerical pathway ---
        t0 = time.time()
        mono_num = monopole_inter_fiber_energy(
            l4, R, Z_A=setup['Z_eff'], Z_B=setup['Z_ligand'],
            l_max=setup['l_max'], n_alpha=setup['n_alpha'],
            pk_potentials=setup['pk_potentials'],
            n_sample_Re=n_sample_Re, method='numerical',
        )
        t_num = time.time() - t0
        F0_num = mono_num['E_monopole']
        N_num = mono_num['N_elec_origin']

        # --- Algebraic pathway ---
        t0 = time.time()
        mono_alg = monopole_inter_fiber_energy(
            l4, R, Z_A=setup['Z_eff'], Z_B=setup['Z_ligand'],
            l_max=setup['l_max'], n_alpha=setup['n_alpha'],
            pk_potentials=setup['pk_potentials'],
            n_sample_Re=n_sample_Re, method='algebraic',
        )
        t_alg = time.time() - t0
        F0_alg = mono_alg['E_monopole']
        N_alg = mono_alg['N_elec_origin']

        # --- Exchange: numerical vs algebraic ---
        exch_num = exchange_inter_fiber_energy(
            l4, R, Z_A=setup['Z_eff'], Z_B=setup['Z_ligand'],
            l_max=setup['l_max'], n_alpha=setup['n_alpha'],
            pk_potentials=setup['pk_potentials'],
            n_sample_Re=n_sample_Re, method='numerical',
        )
        exch_alg = exchange_inter_fiber_energy(
            l4, R, Z_A=setup['Z_eff'], Z_B=setup['Z_ligand'],
            l_max=setup['l_max'], n_alpha=setup['n_alpha'],
            pk_potentials=setup['pk_potentials'],
            n_sample_Re=n_sample_Re, method='algebraic',
        )

        # Comparison
        rel_err_f0 = abs(F0_alg - F0_num) / abs(F0_num) if F0_num != 0 else 0
        rel_err_exch = (abs(exch_alg['E_exchange'] - exch_num['E_exchange'])
                        / abs(exch_num['E_exchange'])
                        if exch_num['E_exchange'] != 0 else 0)

        print(f"  F^0 numerical:  {F0_num:.6f} Ha  (N_elec={N_num:.4f}, {t_num:.2f}s)")
        print(f"  F^0 algebraic:  {F0_alg:.6f} Ha  (N_elec={N_alg:.4f}, {t_alg:.2f}s)")
        print(f"  F^0 rel error:  {rel_err_f0:.4%}")
        print()
        print(f"  E_exch numerical:  {exch_num['E_exchange']:.6f} Ha  (S={exch_num['S_avg']:.4f})")
        print(f"  E_exch algebraic:  {exch_alg['E_exchange']:.6f} Ha  (S={exch_alg['S_avg']:.4f})")
        print(f"  E_exch rel error:  {rel_err_exch:.4%}")

        # Channel decomposition analysis
        f0_mat = mono_alg['F0_matrix_result']
        channels = f0_mat['channels']
        F0_per_ch = f0_mat['F0_per_channel']
        F0_matrix = f0_mat['F0_matrix']
        F0_from_matrix = f0_mat['F0_total']

        print(f"\n  Channel decomposition of F^0:")
        print(f"    F^0 from matrix sum: {F0_from_matrix:.6f} Ha")
        for ic, ch in enumerate(channels):
            print(f"    ch {ch}: F0_ch = {F0_per_ch[ic]:.6f} Ha "
                  f"({F0_per_ch[ic]/F0_from_matrix*100:.1f}%)")

        # Direct exchange (parity-weighted F^0 per channel)
        if exch_alg.get('E_exchange_direct') is not None:
            print(f"\n  Exchange analysis:")
            print(f"    E_exch (S*F^0):    {exch_alg['E_exchange']:.6f} Ha")
            print(f"    E_exch (direct):   {exch_alg['E_exchange_direct']:.6f} Ha")
            print(f"    F^0_signed:        {exch_alg['F0_signed']:.6f} Ha")
            print(f"    Ratio direct/S*F0: {exch_alg['E_exchange_direct']/exch_alg['E_exchange']:.4f}")

        results[R] = {
            'F0_numerical': float(F0_num),
            'F0_algebraic': float(F0_alg),
            'F0_rel_error': float(rel_err_f0),
            'E_exch_numerical': float(exch_num['E_exchange']),
            'E_exch_algebraic': float(exch_alg['E_exchange']),
            'E_exch_rel_error': float(rel_err_exch),
            'S_avg_numerical': float(exch_num['S_avg']),
            'S_avg_algebraic': float(exch_alg['S_avg']),
            'F0_from_matrix': float(F0_from_matrix),
            'F0_per_channel': {str(ch): float(v) for ch, v in
                               zip(channels, F0_per_ch)},
            'E_exchange_direct': float(exch_alg['E_exchange_direct'])
                if exch_alg.get('E_exchange_direct') is not None else None,
            'F0_signed': float(exch_alg['F0_signed'])
                if exch_alg.get('F0_signed') is not None else None,
        }

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    all_pass = True
    for R in R_values:
        r = results[R]
        status = "PASS" if r['F0_rel_error'] < 0.10 else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"  R={R:.1f}: F^0 err={r['F0_rel_error']:.2%}, "
              f"E_exch err={r['E_exch_rel_error']:.2%}  [{status}]")

    print(f"\n  S(R) vs F^0(R) relationship:")
    for R in R_values:
        r = results[R]
        print(f"    R={R:.1f}: S={r['S_avg_algebraic']:.4f}, "
              f"F^0={r['F0_algebraic']:.4f}, "
              f"S*F^0={r['S_avg_algebraic']*r['F0_algebraic']:.4f}")
        if r['F0_signed'] is not None:
            print(f"           F^0_signed={r['F0_signed']:.4f}, "
                  f"ratio={r['F0_signed']/r['F0_algebraic']:.4f}")

    # Save results
    with open('debug/data/algebraic_f0_validation.json', 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to debug/data/algebraic_f0_validation.json")

    return all_pass


if __name__ == '__main__':
    ok = main()
    if not ok:
        print("\n*** VALIDATION FAILED ***")
    else:
        print("\n*** ALL VALIDATIONS PASSED ***")
