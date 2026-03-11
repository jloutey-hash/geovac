"""
Shell-Dependent Cross-Nuclear Attenuation — LiH PES Diagnostic (v0.9.23)

Compares PES with use_shell_radius=True (new) vs use_shell_radius=False (old/v0.9.18).
Runs at nmax=3 with use_dmatrix='hybrid' to match v0.9.18 architecture.

Output: debug/data/shell_radius_lih_results.txt
"""

import warnings
import sys
import numpy as np

sys.path.insert(0, '.')

from geovac.lattice_index import (
    MolecularLatticeIndex, LatticeIndex, compute_bsse_correction,
)


def compute_atomic_energies(nmax: int = 3):
    """Compute separated atom energies at given nmax."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        li = LatticeIndex(
            n_electrons=3, max_n=nmax, nuclear_charge=3,
            vee_method='slater_full', h1_method='exact', fci_method='auto',
        )
        E_li = li.compute_ground_state(n_states=1)[0][0]

        h = LatticeIndex(
            n_electrons=1, max_n=nmax, nuclear_charge=1,
            vee_method='slater_full', h1_method='exact', fci_method='auto',
        )
        E_h = h.compute_ground_state(n_states=1)[0][0]
    return E_li, E_h


def run_pes(R_vals, nmax: int = 3, use_shell_radius: bool = True,
            use_dmatrix='hybrid'):
    """Run PES scan and return energies."""
    energies = []
    for R in R_vals:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mol = MolecularLatticeIndex(
                Z_A=3, Z_B=1, nmax_A=nmax, nmax_B=nmax,
                R=R, n_electrons=4,
                n_bridges=20, vee_method='slater_full',
                fci_method='auto', use_dmatrix=use_dmatrix,
                use_shell_radius=use_shell_radius,
            )
            E, _ = mol.compute_ground_state(n_states=1)
            energies.append(E[0])
    return energies


def main():
    nmax = 3
    R_vals = [1.5, 2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]

    print("=" * 70)
    print("Shell-Dependent Cross-Nuclear Attenuation — LiH PES (v0.9.23)")
    print("=" * 70)
    print(f"nmax = {nmax}")
    print(f"Architecture: hybrid (Fourier diag + SW off-diag)")
    print()

    # Atomic energies
    E_li, E_h = compute_atomic_energies(nmax)
    E_atoms = E_li + E_h
    print(f"Atomic energies:")
    print(f"  E_Li = {E_li:.6f} Ha")
    print(f"  E_H  = {E_h:.6f} Ha")
    print(f"  E_atoms = {E_atoms:.6f} Ha")
    print()

    # BSSE — old (no shell radius)
    print("Computing BSSE (use_shell_radius=False)...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bsse_old = compute_bsse_correction(
            Z_A=3, Z_B=1, nmax_A=nmax, nmax_B=nmax, R=3.015,
            n_electrons_A=3, n_electrons_B=1,
            vee_method='slater_full', fci_method='auto',
            use_shell_radius=False,
        )
    BSSE_old = bsse_old['BSSE']

    # BSSE — new (with shell radius)
    print("Computing BSSE (use_shell_radius=True)...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bsse_new = compute_bsse_correction(
            Z_A=3, Z_B=1, nmax_A=nmax, nmax_B=nmax, R=3.015,
            n_electrons_A=3, n_electrons_B=1,
            vee_method='slater_full', fci_method='auto',
            use_shell_radius=True,
        )
    BSSE_new = bsse_new['BSSE']

    print(f"\nBSSE (old, no shell radius) = {BSSE_old:.6f} Ha")
    print(f"BSSE (new, shell radius)    = {BSSE_new:.6f} Ha")
    print()

    # PES — old (no shell radius, same as v0.9.18)
    print("Running PES (use_shell_radius=False, v0.9.18 baseline)...")
    E_old = run_pes(R_vals, nmax=nmax, use_shell_radius=False)

    # PES — new (with shell radius)
    print("Running PES (use_shell_radius=True, v0.9.23)...")
    E_new = run_pes(R_vals, nmax=nmax, use_shell_radius=True)

    # Output table
    lines = []
    lines.append("Shell-Dependent Cross-Nuclear Attenuation — LiH PES (v0.9.23)")
    lines.append("=" * 70)
    lines.append(f"nmax = {nmax}")
    lines.append(f"Architecture: hybrid (Fourier diag + SW off-diag)")
    lines.append("")
    lines.append(f"Atomic energies:")
    lines.append(f"  E_Li = {E_li:.6f} Ha")
    lines.append(f"  E_H  = {E_h:.6f} Ha")
    lines.append(f"  E_atoms = {E_atoms:.6f} Ha")
    lines.append("")
    lines.append(f"BSSE (old, no shell radius) = {BSSE_old:.6f} Ha")
    lines.append(f"BSSE (new, shell radius)    = {BSSE_new:.6f} Ha")
    lines.append("")

    # Side-by-side comparison
    header = (f"{'R':>8s}  {'E_old':>12s}  {'D_raw_old':>10s}  "
              f"{'D_CP_old':>10s}  {'E_new':>12s}  {'D_raw_new':>10s}  "
              f"{'D_CP_new':>10s}")
    lines.append("PES comparison (old = v0.9.18 baseline, new = v0.9.23):")
    lines.append(header)

    best_De_old = -999.0
    best_R_old = R_vals[0]
    best_De_new = -999.0
    best_R_new = R_vals[0]

    for i, R in enumerate(R_vals):
        D_raw_old = E_atoms - E_old[i]
        D_CP_old = D_raw_old + BSSE_old
        D_raw_new = E_atoms - E_new[i]
        D_CP_new = D_raw_new + BSSE_new

        line = (f"{R:8.3f}  {E_old[i]:12.6f}  {D_raw_old:10.6f}  "
                f"{D_CP_old:10.6f}  {E_new[i]:12.6f}  {D_raw_new:10.6f}  "
                f"{D_CP_new:10.6f}")
        lines.append(line)

        if D_CP_old > best_De_old:
            best_De_old = D_CP_old
            best_R_old = R
        if D_CP_new > best_De_new:
            best_De_new = D_CP_new
            best_R_new = R

    lines.append("")
    lines.append("Summary:")
    lines.append(f"  Old (v0.9.18): R_eq ~ {best_R_old:.3f}, "
                 f"D_e_CP = {best_De_old:.6f} Ha")
    lines.append(f"  New (v0.9.23): R_eq ~ {best_R_new:.3f}, "
                 f"D_e_CP = {best_De_new:.6f} Ha")
    lines.append(f"  Experiment:    R_eq = 3.015,   D_e = 0.092000 Ha")
    lines.append("")

    # Diagnostics
    D_CP_R6_new = (E_atoms - E_new[-1]) + BSSE_new
    D_CP_R6_old = (E_atoms - E_old[-1]) + BSSE_old
    lines.append("Diagnostics:")
    lines.append(f"  1. Dissociation (old): |D_e_CP(R=6)| = "
                 f"{abs(D_CP_R6_old):.6f} {'PASS' if abs(D_CP_R6_old) <= 0.001 else 'FAIL'}")
    lines.append(f"     Dissociation (new): |D_e_CP(R=6)| = "
                 f"{abs(D_CP_R6_new):.6f} {'PASS' if abs(D_CP_R6_new) <= 0.001 else 'FAIL'}")
    lines.append(f"  2. Bound minimum (old): {'PASS' if best_De_old > 0 else 'FAIL'} "
                 f"(D_e_CP = {best_De_old:.6f} at R = {best_R_old})")
    lines.append(f"     Bound minimum (new): {'PASS' if best_De_new > 0 else 'FAIL'} "
                 f"(D_e_CP = {best_De_new:.6f} at R = {best_R_new})")
    lines.append(f"  3. R_eq in [2.5, 3.5] (old): "
                 f"{'PASS' if 2.5 <= best_R_old <= 3.5 else 'FAIL'} "
                 f"(R_eq = {best_R_old})")
    lines.append(f"     R_eq in [2.5, 3.5] (new): "
                 f"{'PASS' if 2.5 <= best_R_new <= 3.5 else 'FAIL'} "
                 f"(R_eq = {best_R_new})")
    lines.append(f"  4. BSSE changed: "
                 f"{'YES' if abs(BSSE_new - BSSE_old) > 0.001 else 'NO'} "
                 f"(delta = {BSSE_new - BSSE_old:.6f})")
    lines.append("")

    # Per-shell correction magnitudes
    lines.append("Per-shell cross-nuclear corrections at R=3.015:")
    R_eq = 3.015
    for n_val in [1, 2, 3]:
        R_eff = MolecularLatticeIndex._shell_effective_R(n_val, 3, R_eq)
        v_old = MolecularLatticeIndex._fourier_cross_attraction(
            n_val, 0, 3, 1, R_eq)
        v_new = MolecularLatticeIndex._fourier_cross_attraction(
            n_val, 0, 3, 1, R_eff)
        pct = 100.0 * (v_new - v_old) / abs(v_old) if abs(v_old) > 1e-15 else 0
        lines.append(f"  Li {n_val}s: R_eff={R_eff:.3f}, V_old={v_old:.6f}, "
                     f"V_new={v_new:.6f} ({pct:+.1f}%)")

    for n_val in [1, 2, 3]:
        R_eff = MolecularLatticeIndex._shell_effective_R(n_val, 1, R_eq)
        v_old = MolecularLatticeIndex._fourier_cross_attraction(
            n_val, 0, 1, 3, R_eq)
        v_new = MolecularLatticeIndex._fourier_cross_attraction(
            n_val, 0, 1, 3, R_eff)
        pct = 100.0 * (v_new - v_old) / abs(v_old) if abs(v_old) > 1e-15 else 0
        lines.append(f"  H  {n_val}s: R_eff={R_eff:.3f}, V_old={v_old:.6f}, "
                     f"V_new={v_new:.6f} ({pct:+.1f}%)")

    output = "\n".join(lines)
    print()
    print(output)

    # Save to file
    with open("debug/data/shell_radius_lih_results.txt", "w") as f:
        f.write(output + "\n")
    print(f"\nSaved to debug/data/shell_radius_lih_results.txt")


if __name__ == "__main__":
    main()
