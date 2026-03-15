"""
Run the full ab initio nuclear lattice pipeline using pre-computed PES data.

This script uses PES data from the prolate CI computation (Eckart + sigma_g/sigma_u)
and completes the Morse fit, nuclear lattice construction, and spectrum comparison.
"""

import numpy as np
import sys
import os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from geovac.nuclear_lattice import NuclearLattice
from morse_fit import fit_morse, compare_to_experiment, morse_potential, EXPERIMENTAL_H2


def main() -> None:
    t0 = time.time()

    print("=" * 70)
    print("AB INITIO NUCLEAR LATTICE FOR H2")
    print("Electron graph -> PES -> Morse fit -> Nuclear lattice -> Spectrum")
    print("Zero experimental input (charges and masses only)")
    print("=" * 70)

    # ================================================================
    # Step 1: PES data (from prolate CI computation)
    # ================================================================
    # Try to load from file first, fall back to hardcoded CI results
    pes_file = os.path.join(script_dir, 'h2_pes_data.txt')
    loaded = False

    if os.path.exists(pes_file):
        try:
            data = np.loadtxt(pes_file)
            if data.shape[0] >= 6:
                R = data[:, 0]
                E_total = data[:, 1]
                loaded = True
                print(f"\nLoaded PES data from {pes_file} ({len(R)} points)")
        except Exception:
            pass

    if not loaded:
        # Hardcoded from prolate CI runs (Eckart + sigma_g/sigma_u)
        print("\nUsing pre-computed prolate CI data")
        R = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
        E_total = np.array([
            -0.961887, -1.066655, -1.107189,
            -1.116793, -1.110971, -1.097598, -1.080932,
        ])

    print(f"\n{'R (bohr)':>10s}  {'E_total (Ha)':>14s}")
    print("-" * 28)
    for i in range(len(R)):
        print(f"{R[i]:10.2f}  {E_total[i]:14.6f}")

    # ================================================================
    # Step 2: Morse fit
    # ================================================================
    print("\n--- STEP 2: Morse Fit ---")
    morse_params = fit_morse(R, E_total, verbose=True)

    if not morse_params:
        print("ERROR: Morse fit failed.")
        return

    # ================================================================
    # Step 3: Compare to experiment
    # ================================================================
    print("\n--- STEP 3: Compare to Experiment ---")
    errors = compare_to_experiment(morse_params, verbose=True)

    # ================================================================
    # Step 4: Build nuclear lattices
    # ================================================================
    print("\n--- STEP 4: Build Nuclear Lattices ---")

    nuc_ai = NuclearLattice(
        D_e=morse_params['D_e'],
        omega_e=morse_params['omega_e_cm'],
        B_e=morse_params['B_e_cm'],
        alpha_e=morse_params['alpha_e_cm'],
        J_max=10,
    )

    exp = EXPERIMENTAL_H2
    nuc_exp = NuclearLattice(
        D_e=exp['D_e'],
        omega_e=exp['omega_e_cm'],
        B_e=exp['B_e_cm'],
        alpha_e=exp['alpha_e_cm'],
        J_max=10,
    )

    print(f"\nAb initio lattice: {nuc_ai.vib.n_states} vibrational states, J_max=10")
    print(f"Experimental lattice: {nuc_exp.vib.n_states} vibrational states")

    # ================================================================
    # Step 5: Compare transitions
    # ================================================================
    print("\n--- STEP 5: Rovibrational Spectrum ---")
    HARTREE_TO_CM = 219474.63
    transitions = {}

    if nuc_ai.vib.n_states >= 2 and nuc_exp.vib.n_states >= 2:
        E01_ai = (nuc_ai.vib.morse_energy(1) - nuc_ai.vib.morse_energy(0)) * HARTREE_TO_CM
        E01_exp = (nuc_exp.vib.morse_energy(1) - nuc_exp.vib.morse_energy(0)) * HARTREE_TO_CM
        transitions['v01'] = {
            'ai': E01_ai, 'exp': E01_exp,
            'err_pct': (E01_ai - E01_exp) / E01_exp * 100,
            'label': 'v=0->1 (fundamental)',
            'exp_nist': 4161.0,
        }

    if nuc_ai.vib.n_states >= 3 and nuc_exp.vib.n_states >= 3:
        E02_ai = (nuc_ai.vib.morse_energy(2) - nuc_ai.vib.morse_energy(0)) * HARTREE_TO_CM
        E02_exp = (nuc_exp.vib.morse_energy(2) - nuc_exp.vib.morse_energy(0)) * HARTREE_TO_CM
        transitions['v02'] = {
            'ai': E02_ai, 'exp': E02_exp,
            'err_pct': (E02_ai - E02_exp) / E02_exp * 100,
            'label': 'v=0->2 (first overtone)',
            'exp_nist': 8087.0,
        }

    E_rot_ai = (nuc_ai.rovibrational_energy(0, 1)
                - nuc_ai.rovibrational_energy(0, 0)) * HARTREE_TO_CM
    E_rot_exp = (nuc_exp.rovibrational_energy(0, 1)
                 - nuc_exp.rovibrational_energy(0, 0)) * HARTREE_TO_CM
    transitions['J01'] = {
        'ai': E_rot_ai, 'exp': E_rot_exp,
        'err_pct': (E_rot_ai - E_rot_exp) / E_rot_exp * 100,
        'label': 'J=0->1 (rotational)',
        'exp_nist': 118.0,
    }

    if nuc_ai.vib.n_states >= 2:
        E_R0_ai = (nuc_ai.rovibrational_energy(1, 1)
                   - nuc_ai.rovibrational_energy(0, 0)) * HARTREE_TO_CM
        E_R0_exp = (nuc_exp.rovibrational_energy(1, 1)
                    - nuc_exp.rovibrational_energy(0, 0)) * HARTREE_TO_CM
        transitions['R0'] = {
            'ai': E_R0_ai, 'exp': E_R0_exp,
            'err_pct': (E_R0_ai - E_R0_exp) / E_R0_exp * 100,
            'label': 'R(0): (0,0)->(1,1)',
            'exp_nist': 4497.8,
        }

    print(f"\n{'Transition':<25s}  {'Ab initio':>10s}  {'Expt':>10s}  "
          f"{'NIST':>10s}  {'Err %':>8s}")
    print("-" * 68)
    for key, t in transitions.items():
        nist_str = f"{t['exp_nist']:.1f}" if 'exp_nist' in t else "---"
        print(f"  {t['label']:<25s}  {t['ai']:10.1f}  {t['exp']:10.1f}  "
              f"{nist_str:>10s}  {t['err_pct']:+8.2f}%")

    # ================================================================
    # Step 6: Plot
    # ================================================================
    print("\n--- STEP 6: Plot ---")
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    valid = ~np.isnan(E_total)
    ax.plot(R[valid], E_total[valid], 'ko', markersize=8,
            label='Ab initio (Prolate CI)', zorder=5)

    R_fine = np.linspace(0.6, 7.0, 500)
    D_e = morse_params['D_e']
    a = morse_params['a']
    R_e = morse_params['R_e']
    E_atoms = -1.0
    V_fit = morse_potential(R_fine, D_e, a, R_e) + E_atoms
    ax.plot(R_fine, V_fit, 'b-', linewidth=2,
            label=f'Morse fit (D_e={D_e:.4f} Ha)')

    E_min = E_atoms - D_e
    ax.axhline(y=E_atoms, color='gray', linestyle='--', alpha=0.5,
               label='E(H+H) = -1.0 Ha')
    ax.plot(R_e, E_min, 'r*', markersize=15, zorder=10,
            label=f'R_e = {R_e:.3f} bohr')
    ax.annotate(f'D_e = {D_e:.4f} Ha\n= {D_e*27.2114:.3f} eV',
                xy=(R_e, E_min), xytext=(R_e + 1.5, E_min + 0.02),
                fontsize=10, arrowprops=dict(arrowstyle='->', color='red'),
                color='red')

    ax.set_xlabel('R (bohr)', fontsize=13)
    ax.set_ylabel('E_total (Hartree)', fontsize=13)
    ax.set_title('H$_2$ Potential Energy Surface: Ab Initio from Graph Laplacian',
                 fontsize=14)
    ax.legend(fontsize=10, loc='lower right')
    ax.set_xlim(0.5, 7.0)
    y_lo = min(E_min - 0.05, np.nanmin(E_total[valid]) - 0.05)
    ax.set_ylim(y_lo, E_atoms + 0.15)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    plot_file = os.path.join(script_dir, 'h2_pes_morse.png')
    plt.savefig(plot_file, dpi=150)
    plt.close()
    print(f"Plot saved to {plot_file}")

    # ================================================================
    # Step 7: Write results.md
    # ================================================================
    results_file = os.path.join(script_dir, 'results.md')
    with open(results_file, 'w') as f:
        f.write("# Ab Initio Nuclear Lattice: H2 Rovibrational Spectrum\n\n")
        f.write("**Date:** 2026-03-14\n\n")
        f.write("**Method:** Prolate spheroidal CI (Eckart Z_eff + sigma_g/sigma_u "
                "2x2 CI) -> Morse fit -> Nuclear Lattice\n\n")
        f.write("**Input:** Nuclear charges (Z=1) and hydrogen mass only. "
                "Zero experimental spectroscopic input.\n\n")

        f.write("## PES Data\n\n")
        f.write("| R (bohr) | E_total (Ha) |\n")
        f.write("|----------:|-------------:|\n")
        for i in range(len(R)):
            if not np.isnan(E_total[i]):
                f.write(f"| {R[i]:.2f} | {E_total[i]:.6f} |\n")

        f.write("\n## Morse Parameters\n\n")
        f.write("| Parameter | Ab Initio | Experiment | Error |\n")
        f.write("|:----------|----------:|-----------:|------:|\n")

        params_table = [
            ('D_e (Ha)', 'D_e', '.6f', '.4f'),
            ('R_e (bohr)', 'R_e', '.4f', '.3f'),
            ('omega_e (cm-1)', 'omega_e_cm', '.2f', '.2f'),
            ('omega_e*x_e (cm-1)', 'omega_e_xe_cm', '.2f', '.2f'),
            ('B_e (cm-1)', 'B_e_cm', '.3f', '.3f'),
            ('alpha_e (cm-1)', 'alpha_e_cm', '.3f', '.3f'),
        ]

        for label, key, fmt_ai, fmt_exp in params_table:
            ai_val = morse_params.get(key, float('nan'))
            exp_val = EXPERIMENTAL_H2.get(key, float('nan'))
            err = errors.get(key, float('nan'))
            f.write(f"| {label} | {ai_val:{fmt_ai}} | {exp_val:{fmt_exp}} | "
                    f"{err:+.2f}% |\n")

        f.write("\n## Rovibrational Transitions\n\n")
        f.write("| Transition | Ab Initio (cm-1) | Expt (cm-1) | "
                "NIST (cm-1) | Error |\n")
        f.write("|:-----------|------------------:|------------:|"
                "------------:|------:|\n")

        for key, t in transitions.items():
            nist = t.get('exp_nist', float('nan'))
            nist_str = f"{nist:.1f}" if not np.isnan(nist) else "---"
            f.write(f"| {t['label']} | {t['ai']:.1f} | {t['exp']:.1f} | "
                    f"{nist_str} | {t['err_pct']:+.2f}% |\n")

        f.write("\n## Notes\n\n")
        f.write("- PES computed via prolate spheroidal CI with Eckart Z_eff "
                "optimization\n")
        f.write("- Two orbitals: sigma_g (bonding) + sigma_u (antibonding) from "
                "exact H2+ solutions\n")
        f.write("- 2x2 CI matrix captures ionic-covalent mixing\n")
        f.write("- V_ee computed via azimuthal averaging (elliptic K integral)\n")
        f.write("- Morse V(R) = D_e[1-exp(-a(R-R_e))]^2 - D_e fitted via "
                "least squares\n")
        f.write("- Spectroscopic constants derived analytically from Morse "
                "parameters\n")
        f.write("- Nuclear lattice: Morse SU(2) vibrational chain (finite "
                "representation) + SO(3) rotational paraboloid\n")
        f.write("- Vibration-rotation coupling alpha_e included\n")
        f.write("- D_e underestimated by ~33% (missing dynamical correlation "
                "beyond 2x2 CI)\n")
        f.write("- R_e, omega_e, and B_e accurate to ~1-3% (shape of PES well "
                "preserved)\n")

    print(f"Results written to {results_file}")

    dt = time.time() - t0
    print(f"\n{'='*70}")
    print(f"Total pipeline time: {dt:.1f}s")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
