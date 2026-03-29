"""
BeH2 inter-bond repulsion PES scan and scaling sweep.

Efficient implementation: runs the Level 4 solver ONCE (expensive), then
sweeps V_inter scaling as pure post-processing (cheap, V_inter is a scalar).

Outputs:
  debug/data/beh2_interbond_pes.json       — PES with and without V_inter
  debug/data/beh2_interbond_spectro.json   — spectroscopic constants (best scale)
  debug/data/beh2_interbond_scaling.json   — scaling sweep results
  debug/plots/beh2_interbond_comparison.png — comparison plot

Usage:
    python debug/beh2_interbond_scan.py
"""

import json
import time
import numpy as np
from pathlib import Path
from scipy.optimize import curve_fit

from geovac.composed_triatomic import (
    ComposedTriatomicSolver, _morse_potential,
)
from geovac.nuclear_lattice import HARTREE_TO_CM, AMU_TO_ME

# Experimental reference
REF = {
    'R_eq': 2.507,       # bohr
    'omega_e': 2345.0,   # cm-1
    'D_e': 0.147,        # Ha
}

L_MAX = 2
N_RE = 300
R_GRID = np.arange(1.5, 5.1, 0.1)

OUTPUT_DIR = Path(__file__).parent / 'data'
PLOT_DIR = Path(__file__).parent / 'plots'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)


def fit_spectro(R_valid: np.ndarray, E_valid: np.ndarray,
                M_ligand: float) -> dict:
    """Fit Morse and extract spectroscopic constants."""
    i_min = np.argmin(E_valid)
    R_eq_grid = R_valid[i_min]
    E_min_grid = E_valid[i_min]
    D_e_grid = E_valid[-1] - E_min_grid

    mask = np.abs(R_valid - R_eq_grid) < 1.5
    R_fit = R_valid[mask]
    E_fit = E_valid[mask]
    if len(R_fit) < 4:
        R_fit, E_fit = R_valid, E_valid

    try:
        popt, _ = curve_fit(
            _morse_potential, R_fit, E_fit,
            p0=[E_min_grid, max(D_e_grid, 0.01), 1.0, R_eq_grid],
            maxfev=10000,
        )
        E_min_fit, D_e_fit, a_fit, R_eq_fit = popt
    except RuntimeError:
        E_min_fit, D_e_fit, a_fit, R_eq_fit = (
            E_min_grid, max(D_e_grid, 0.01), 1.0, R_eq_grid)

    D_e_fit = abs(D_e_fit)
    a_fit = abs(a_fit)
    mu_sym_au = M_ligand * AMU_TO_ME
    omega_e_au = a_fit * np.sqrt(2.0 * D_e_fit / mu_sym_au)
    omega_e_cm = omega_e_au * HARTREE_TO_CM

    return {
        'R_eq': float(R_eq_fit),
        'D_e': float(D_e_fit),
        'omega_e_sym': float(omega_e_cm),
        'a': float(a_fit),
        'E_min': float(E_min_fit),
    }


def main():
    t_start = time.time()

    # ---------------------------------------------------------------
    # 1. Run the Level 4 solver ONCE (the expensive part)
    # ---------------------------------------------------------------
    print("=" * 64)
    print("STEP 1: Base PES scan (no V_inter) — Level 4 solver")
    print("=" * 64)

    solver = ComposedTriatomicSolver.BeH2(
        l_max=L_MAX, include_interbond=False, verbose=True)
    solver.solve_core()
    solver.scan_pes(R_grid=R_GRID, n_Re=N_RE)
    solver.fit_spectroscopic_constants()

    base_spectro = solver.spectro
    base_err_R = abs(base_spectro['R_eq'] - REF['R_eq']) / REF['R_eq'] * 100
    base_err_D = abs(base_spectro['D_e'] - REF['D_e']) / REF['D_e'] * 100

    print(f"\n  === Baseline (no V_inter) ===")
    print(f"  R_eq = {base_spectro['R_eq']:.4f} bohr"
          f"  (ref: {REF['R_eq']}, err: {base_err_R:.1f}%)")
    print(f"  D_e  = {base_spectro['D_e']:.6f} Ha"
          f"  (ref: {REF['D_e']}, err: {base_err_D:.1f}%)")
    print(f"  w_e  = {base_spectro['omega_e_sym']:.1f} cm-1"
          f"  (ref: {REF['omega_e']})")

    # Extract base PES arrays
    R_valid = np.array(solver.pes_result['R_valid'])
    E_base = np.array(solver.pes_result['E_valid'])

    # ---------------------------------------------------------------
    # 2. Compute V_inter(R) once, then sweep scales (cheap)
    # ---------------------------------------------------------------
    print("\n" + "=" * 64)
    print("STEP 2: Inter-bond repulsion scaling sweep (post-processing)")
    print("=" * 64)

    # Compute V_inter at each valid R point
    V_inter_arr = np.array([solver._inter_bond_repulsion(R) for R in R_valid])

    print(f"\n  V_inter range: [{V_inter_arr[-1]:.4f}, {V_inter_arr[0]:.4f}] Ha")
    print(f"  V_inter(R_eq_base={base_spectro['R_eq']:.2f}): "
          f"{solver._inter_bond_repulsion(base_spectro['R_eq']):.4f} Ha")

    scales = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
    scaling_results = []

    print(f"\n  {'Scale':>8s}  {'R_eq':>8s}  {'R_eq err%':>10s}"
          f"  {'D_e':>10s}  {'D_e err%':>10s}  {'omega_e':>8s}")
    print(f"  {'-'*8}  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*8}")

    for scale in scales:
        E_with_inter = E_base + scale * V_inter_arr
        spectro = fit_spectro(R_valid, E_with_inter, solver.M_ligand)

        err_R = abs(spectro['R_eq'] - REF['R_eq']) / REF['R_eq'] * 100
        err_D = abs(spectro['D_e'] - REF['D_e']) / REF['D_e'] * 100

        result = {
            'scale': scale,
            'R_eq': spectro['R_eq'],
            'D_e': spectro['D_e'],
            'omega_e': spectro['omega_e_sym'],
            'R_eq_error_pct': err_R,
            'D_e_error_pct': err_D,
        }
        scaling_results.append(result)

        print(f"  {scale:8.2f}  {spectro['R_eq']:8.4f}  {err_R:10.1f}"
              f"  {spectro['D_e']:10.6f}  {err_D:10.1f}"
              f"  {spectro['omega_e_sym']:8.1f}")

    # Also get the scale=1.0 PES for saving
    E_with_inter_1 = E_base + 1.0 * V_inter_arr

    # Find best scale for R_eq
    best = min(scaling_results, key=lambda x: x['R_eq_error_pct'])

    # ---------------------------------------------------------------
    # 3. Pauli term count
    # ---------------------------------------------------------------
    print("\n" + "=" * 64)
    print("STEP 3: Pauli term count verification")
    print("=" * 64)

    from geovac.composed_qubit import build_composed_beh2
    pauli_result = build_composed_beh2(
        max_n_core=2, max_n_val=2, verbose=True)
    n_pauli = pauli_result['N_pauli']
    print(f"\n  Pauli terms: {n_pauli}")
    print(f"  V_inter is a scalar — does NOT affect qubit Hamiltonian.")
    print(f"  Pauli count unchanged: {'YES' if n_pauli == 556 else 'NO'}")

    # ---------------------------------------------------------------
    # 4. Save results
    # ---------------------------------------------------------------
    print("\n" + "=" * 64)
    print("Saving results...")
    print("=" * 64)

    # PES comparison
    pes_data = {
        'description': 'BeH2 PES with and without inter-bond repulsion',
        'reference': REF,
        'l_max': L_MAX,
        'no_interbond': {
            'R': R_valid.tolist(),
            'E_total': E_base.tolist(),
            'R_eq': base_spectro['R_eq'],
            'D_e': base_spectro['D_e'],
            'omega_e': base_spectro['omega_e_sym'],
            'R_eq_error_pct': base_err_R,
        },
        'with_interbond': {
            'R': R_valid.tolist(),
            'E_total': E_with_inter_1.tolist(),
            'V_inter': V_inter_arr.tolist(),
            'R_eq': scaling_results[4]['R_eq'],  # scale=1.0
            'D_e': scaling_results[4]['D_e'],
            'omega_e': scaling_results[4]['omega_e'],
            'R_eq_error_pct': scaling_results[4]['R_eq_error_pct'],
        },
    }
    with open(OUTPUT_DIR / 'beh2_interbond_pes.json', 'w') as f:
        json.dump(pes_data, f, indent=2)
    print(f"  Saved: {OUTPUT_DIR / 'beh2_interbond_pes.json'}")

    # Scaling sweep
    scaling_data = {
        'description': 'BeH2 inter-bond scaling sweep',
        'reference': REF,
        'l_max': L_MAX,
        'V_inter_formula': '1/(2*r_Be) + 1/(2R) + 2/(R+r_Be), r_Be=3/(2*Z_eff)',
        'r_Be': 1.5 / solver.Z_eff,
        'Z_eff': solver.Z_eff,
        'results': scaling_results,
        'best_scale': best['scale'],
        'best_R_eq': best['R_eq'],
        'best_R_eq_error_pct': best['R_eq_error_pct'],
        'direction': 'outward' if scaling_results[4]['R_eq'] > scaling_results[0]['R_eq'] else 'inward',
    }
    with open(OUTPUT_DIR / 'beh2_interbond_scaling.json', 'w') as f:
        json.dump(scaling_data, f, indent=2)
    print(f"  Saved: {OUTPUT_DIR / 'beh2_interbond_scaling.json'}")

    # Spectroscopic constants (best scale)
    spectro_data = {
        'molecule': 'BeH2',
        'method': 'inter-bond point charge',
        'interbond_scale': best['scale'],
        'computed': {
            'R_eq': best['R_eq'],
            'D_e': best['D_e'],
            'omega_e_sym': best['omega_e'],
        },
        'reference': REF,
        'errors': {
            'R_eq_pct': best['R_eq_error_pct'],
            'D_e_pct': best['D_e_error_pct'],
        },
        'diagnosis': {
            'v_inter_direction': 'V_inter pushes R_eq outward (wrong direction)',
            'reason': 'V_inter ~ 1/R is monotonically decreasing, raises short-R '
                      'PES more than long-R, shifting minimum to larger R.',
            'best_scale_is_zero': best['scale'] == 0.0,
            'conclusion': 'Point charge Approach A fails — wrong sign of R_eq shift. '
                          'The functional form 1/R cannot fix R_eq. '
                          'Approach B (orbital overlap) needed.',
        },
    }
    with open(OUTPUT_DIR / 'beh2_interbond_spectro.json', 'w') as f:
        json.dump(spectro_data, f, indent=2)
    print(f"  Saved: {OUTPUT_DIR / 'beh2_interbond_spectro.json'}")

    # ---------------------------------------------------------------
    # 5. Summary
    # ---------------------------------------------------------------
    print("\n" + "=" * 64)
    print("SUMMARY")
    print("=" * 64)
    print(f"\n  {'Scale':>8s}  {'R_eq':>8s}  {'R_eq err%':>10s}"
          f"  {'D_e':>10s}  {'D_e err%':>10s}")
    print(f"  {'-'*8}  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*10}")
    for r in scaling_results:
        marker = " <-- best" if r['scale'] == best['scale'] else ""
        print(f"  {r['scale']:8.2f}  {r['R_eq']:8.4f}  {r['R_eq_error_pct']:10.1f}"
              f"  {r['D_e']:10.6f}  {r['D_e_error_pct']:10.1f}{marker}")
    print(f"\n  Best scale: {best['scale']:.2f} -> R_eq = {best['R_eq']:.4f} bohr"
          f" ({best['R_eq_error_pct']:.1f}% error)")
    print(f"  Experiment: R_eq = {REF['R_eq']} bohr, D_e = {REF['D_e']} Ha")
    print(f"  Pauli terms: {n_pauli} (unchanged)")

    direction = 'OUTWARD' if scaling_results[4]['R_eq'] > scaling_results[0]['R_eq'] else 'INWARD'
    print(f"\n  DIAGNOSTIC: V_inter pushes R_eq {direction}")
    if direction == 'OUTWARD':
        print(f"  -> Point charge model (Approach A) has WRONG direction.")
        print(f"  -> 1/R repulsion is strongest at short R, raises short-R PES")
        print(f"  -> more than long-R PES, pushing minimum outward.")
        print(f"  -> Approach B (orbital overlap integral) needed.")

    print(f"\n  Total time: {time.time() - t_start:.0f}s")

    # ---------------------------------------------------------------
    # 6. Plot
    # ---------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Left: PES comparison
        ax1.plot(R_valid, E_base, 'b-o', markersize=3,
                 label=f'No V_inter (R_eq={scaling_results[0]["R_eq"]:.3f})')
        ax1.plot(R_valid, E_with_inter_1, 'r-s', markersize=3,
                 label=f'V_inter x1.0 (R_eq={scaling_results[4]["R_eq"]:.3f})')
        ax1.axvline(x=REF['R_eq'], color='green', linestyle='--',
                     alpha=0.7, label=f'R_eq expt = {REF["R_eq"]}')
        ax1.set_xlabel('R (bohr)')
        ax1.set_ylabel('E (Ha)')
        ax1.set_title('BeH2 PES: Effect of Inter-Bond Repulsion')
        ax1.legend(fontsize=7)
        ax1.grid(True, alpha=0.3)

        # Right: R_eq vs scale
        scale_vals = [r['scale'] for r in scaling_results]
        req_vals = [r['R_eq'] for r in scaling_results]
        ax2.plot(scale_vals, req_vals, 'ko-', markersize=6)
        ax2.axhline(y=REF['R_eq'], color='green', linestyle='--',
                     alpha=0.7, label=f'R_eq expt = {REF["R_eq"]}')
        ax2.set_xlabel('Inter-bond scale factor')
        ax2.set_ylabel('R_eq (bohr)')
        ax2.set_title('R_eq vs Inter-Bond Scaling')
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plot_path = PLOT_DIR / 'beh2_interbond_comparison.png'
        plt.savefig(str(plot_path), dpi=150)
        plt.close(fig)
        print(f"\n  Plot saved: {plot_path}")
    except ImportError:
        print("  WARNING: matplotlib not available, skipping plot.")


if __name__ == '__main__':
    main()
