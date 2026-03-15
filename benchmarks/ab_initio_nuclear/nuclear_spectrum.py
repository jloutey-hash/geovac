"""
Ab initio nuclear lattice construction and rovibrational spectrum for H2.

Full pipeline:
1. Compute H2 PES using prolate spheroidal Eckart HF
2. Fit Morse potential to extract D_e, R_e, a
3. Derive spectroscopic constants omega_e, omega_e_xe, B_e, alpha_e
4. Build nuclear lattice (Morse SU(2) + SO(3))
5. Compare rovibrational transitions to experiment
6. Plot PES with Morse fit

Zero experimental input: only nuclear charges (Z=1) and masses (m_H).
"""

import numpy as np
import sys
import os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

from typing import Dict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from geovac.nuclear_lattice import NuclearLattice
from h2_pes import compute_pes_hf, compute_pes_ci, save_pes, R_VALUES
from morse_fit import fit_morse, compare_to_experiment, morse_potential, EXPERIMENTAL_H2


def build_nuclear_lattices(
    ai_params: Dict[str, float],
    J_max: int = 10,
    verbose: bool = True,
) -> Dict:
    """Build nuclear lattices from ab initio and experimental parameters."""
    nuc_ai = NuclearLattice(
        D_e=ai_params['D_e'],
        omega_e=ai_params['omega_e_cm'],
        B_e=ai_params['B_e_cm'],
        alpha_e=ai_params['alpha_e_cm'],
        J_max=J_max,
    )

    exp = EXPERIMENTAL_H2
    nuc_exp = NuclearLattice(
        D_e=exp['D_e'],
        omega_e=exp['omega_e_cm'],
        B_e=exp['B_e_cm'],
        alpha_e=exp['alpha_e_cm'],
        J_max=J_max,
    )

    if verbose:
        print(f"\nAb initio lattice: {nuc_ai.vib.n_states} vibrational states, "
              f"J_max={J_max}")
        print(f"Experimental lattice: {nuc_exp.vib.n_states} vibrational states")

    return {'ai': nuc_ai, 'exp': nuc_exp}


def compare_transitions(
    nuc_ai: NuclearLattice,
    nuc_exp: NuclearLattice,
    verbose: bool = True,
) -> Dict[str, Dict[str, float]]:
    """Compare rovibrational transitions between ab initio and experiment."""
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

    if verbose:
        print("\n" + "=" * 78)
        print("ROVIBRATIONAL TRANSITION COMPARISON")
        print("=" * 78)
        print(f"  {'Transition':<25s}  {'Ab initio':>10s}  {'Expt':>10s}  "
              f"{'NIST':>10s}  {'Err %':>8s}")
        print("-" * 68)
        for key, t in transitions.items():
            nist_str = f"{t['exp_nist']:.1f}" if 'exp_nist' in t else "---"
            print(f"  {t['label']:<25s}  {t['ai']:10.1f}  {t['exp']:10.1f}  "
                  f"{nist_str:>10s}  {t['err_pct']:+8.2f}%")

    return transitions


def plot_pes_and_fit(
    R: np.ndarray,
    E_total: np.ndarray,
    morse_params: Dict[str, float],
    outfile: str,
    verbose: bool = True,
) -> None:
    """Plot PES with Morse fit."""
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    valid = ~np.isnan(E_total)
    ax.plot(R[valid], E_total[valid], 'ko', markersize=8,
            label='Ab initio (Eckart HF)', zorder=5)

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
    plt.savefig(outfile, dpi=150)
    plt.close()

    if verbose:
        print(f"\nPlot saved to {outfile}")


def write_results_md(
    morse_params: Dict[str, float],
    errors: Dict[str, float],
    transitions: Dict[str, Dict[str, float]],
    pes_data: Dict[str, np.ndarray],
    outfile: str,
) -> None:
    """Write comparison tables to results.md."""
    exp = EXPERIMENTAL_H2

    with open(outfile, 'w') as f:
        f.write("# Ab Initio Nuclear Lattice: H2 Rovibrational Spectrum\n\n")
        f.write("**Date:** 2026-03-14\n\n")
        f.write("**Method:** Prolate spheroidal Eckart HF (Z_eff optimized) "
                "-> Morse fit -> Nuclear Lattice\n\n")
        f.write("**Input:** Nuclear charges (Z=1) and hydrogen mass only. "
                "Zero experimental spectroscopic input.\n\n")

        f.write("## PES Data\n\n")
        f.write("| R (bohr) | E_total (Ha) |\n")
        f.write("|----------:|-------------:|\n")
        for i in range(len(pes_data['R'])):
            if not np.isnan(pes_data['E_total'][i]):
                f.write(f"| {pes_data['R'][i]:.2f} | "
                        f"{pes_data['E_total'][i]:.6f} |\n")

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
            exp_val = exp.get(key, float('nan'))
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
        f.write("- PES computed with Eckart variational HF in prolate spheroidal "
                "coordinates\n")
        f.write("- Z_eff optimized at each R: solves H2+ at R*Z_eff, then "
                "evaluates physical energy via scaling\n")
        f.write("- V_ee (Coulomb integral J_gg) computed via azimuthal averaging "
                "(elliptic K kernel)\n")
        f.write("- Morse V(R) = D_e[1-exp(-a(R-R_e))]^2 - D_e fitted via "
                "least squares\n")
        f.write("- Spectroscopic constants derived analytically from Morse "
                "parameters\n")
        f.write("- Nuclear lattice: Morse SU(2) vibrational chain + SO(3) "
                "rotational paraboloid\n")
        f.write("- Vibration-rotation coupling alpha_e included\n")
        f.write("- HF misses ~20% of D_e (electron correlation); spectroscopic "
                "constants inherit this systematic shift\n")


def main() -> None:
    """Run the full ab initio nuclear lattice pipeline."""
    t0_total = time.time()

    print("=" * 70)
    print("AB INITIO NUCLEAR LATTICE FOR H2")
    print("Electron graph -> PES -> Morse fit -> Nuclear lattice -> Spectrum")
    print("Zero experimental input (charges and masses only)")
    print("=" * 70)

    # Step 1: Compute PES
    print("\n--- STEP 1: Compute PES (Eckart HF) ---")
    pes = compute_pes_hf(verbose=True)
    save_pes(pes, os.path.join(script_dir, 'h2_pes_data.txt'))

    # Step 2: Morse fit
    print("\n--- STEP 2: Morse Fit ---")
    morse_params = fit_morse(pes['R'], pes['E_total'], verbose=True)
    if not morse_params:
        print("ERROR: Morse fit failed.")
        return

    # Step 3: Compare to experiment
    print("\n--- STEP 3: Compare to Experiment ---")
    errors = compare_to_experiment(morse_params, verbose=True)

    # Step 4: Build nuclear lattices
    print("\n--- STEP 4: Build Nuclear Lattices ---")
    lattices = build_nuclear_lattices(morse_params, J_max=10, verbose=True)

    # Step 5: Compare transitions
    print("\n--- STEP 5: Rovibrational Spectrum ---")
    transitions = compare_transitions(
        lattices['ai'], lattices['exp'], verbose=True)

    # Step 6: Plot
    print("\n--- STEP 6: Plot ---")
    plot_file = os.path.join(script_dir, 'h2_pes_morse.png')
    plot_pes_and_fit(pes['R'], pes['E_total'], morse_params,
                     plot_file, verbose=True)

    # Step 7: Write results
    results_file = os.path.join(script_dir, 'results.md')
    write_results_md(morse_params, errors, transitions, pes, results_file)
    print(f"Results written to {results_file}")

    dt_total = time.time() - t0_total
    print(f"\n{'='*70}")
    print(f"Total pipeline time: {dt_total:.1f}s")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
