#!/usr/bin/env python3
"""
Bridge Weight Correction Validation
=====================================

Validates the STO-overlap bridge weight correction by:
1. Scanning the H2 Full CI PES to find R_eq (should be 1.40 ± 0.02 Bohr)
2. Verifying correlation energy error remains < 0.5%
3. Comparing old (pure exp) vs new (STO overlap) weight profiles

Usage:
    python debug/debug_bridge_correction.py

Output:
    debug/plots/bridge_correction_validation.png
"""

import sys
import os
import io
import time

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from geovac import GeometricLattice, MoleculeHamiltonian


# -- Parameters --
MAX_N = 4
N_BRIDGES = 20
KINETIC_SCALE = -1 / 16
BRIDGE_DECAY_RATE = 1.0

# Experimental references
R_EQ_EXP = 1.401      # Bohr
E_H2_EXP = -1.1745    # Hartree (electronic energy, no V_NN)
R_EQ_TOL = 0.02       # Bohr tolerance on R_eq
ENERGY_TOL = 0.005     # 0.5% relative error on correlation energy


def compute_energy_silent(R: float) -> float:
    """Compute H2 Full CI energy at bond length R, suppressing output."""
    atom_A = GeometricLattice(
        max_n=MAX_N, nucleus_position=(0.0, 0.0, 0.0), nuclear_charge=1,
    )
    atom_B = GeometricLattice(
        max_n=MAX_N, nucleus_position=(R, 0.0, 0.0), nuclear_charge=1,
    )
    h2 = MoleculeHamiltonian(
        lattices=[atom_A, atom_B],
        connectivity=[(0, 1, N_BRIDGES)],
        kinetic_scale=KINETIC_SCALE,
        bridge_decay_rate=BRIDGE_DECAY_RATE,
    )

    _stdout = sys.stdout
    sys.stdout = io.StringIO()
    try:
        E_ci, _ = h2.compute_ground_state(n_states=1, method="full_ci")
    finally:
        sys.stdout = _stdout

    return E_ci[0]


def run_validation() -> None:
    """Full validation of bridge weight correction."""
    print("=" * 70)
    print("  Bridge Weight Correction — STO Overlap Validation")
    print("=" * 70)

    # ----------------------------------------------------------------
    # 1. PES scan: find R_eq
    # ----------------------------------------------------------------
    print("\n[1/3] Scanning Full CI PES around equilibrium...")
    R_values = np.arange(0.8, 3.01, 0.05)
    E_values = []

    t0 = time.time()
    for i, R in enumerate(R_values):
        E = compute_energy_silent(R)
        E_values.append(E)
        if (i + 1) % 10 == 0 or i == 0:
            print(f"       R = {R:.2f} Bohr  ->  E = {E:.6f} Ha")

    E_values = np.array(E_values)
    elapsed_pes = time.time() - t0
    print(f"       PES scan: {len(R_values)} points in {elapsed_pes:.1f}s")

    # Find minimum
    idx_min = np.argmin(E_values)
    R_eq = R_values[idx_min]
    E_min = E_values[idx_min]

    # Refine with parabolic interpolation around minimum
    if 1 <= idx_min <= len(R_values) - 2:
        R_a, R_b, R_c = R_values[idx_min - 1], R_values[idx_min], R_values[idx_min + 1]
        E_a, E_b, E_c = E_values[idx_min - 1], E_values[idx_min], E_values[idx_min + 1]
        denom = 2.0 * ((R_b - R_a) * (E_b - E_c) - (R_b - R_c) * (E_b - E_a))
        if abs(denom) > 1e-15:
            R_eq_refined = R_b - (
                (R_b - R_a)**2 * (E_b - E_c) - (R_b - R_c)**2 * (E_b - E_a)
            ) / denom
            E_min_refined = compute_energy_silent(R_eq_refined)
            if E_min_refined < E_min:
                R_eq = R_eq_refined
                E_min = E_min_refined

    print(f"\n       R_eq (computed)  = {R_eq:.4f} Bohr")
    print(f"       R_eq (expt)     = {R_EQ_EXP:.3f} Bohr")
    print(f"       Delta R         = {abs(R_eq - R_EQ_EXP):.4f} Bohr")
    print(f"       E_min           = {E_min:.6f} Ha")

    r_eq_pass = abs(R_eq - R_EQ_EXP) <= R_EQ_TOL
    print(f"       R_eq test:      {'PASS' if r_eq_pass else 'FAIL'} "
          f"(tolerance ±{R_EQ_TOL} Bohr)")

    # ----------------------------------------------------------------
    # 2. Energy precision check at R_eq
    # ----------------------------------------------------------------
    print("\n[2/3] Checking correlation energy precision...")

    # Separated atom energy (2 × H atom ground state)
    E_atom = compute_energy_silent(20.0)  # Large R → separated atoms
    E_binding = E_min - E_atom
    print(f"       E_atoms (separated) = {E_atom:.6f} Ha")
    print(f"       E_binding           = {E_binding:.6f} Ha")
    print(f"       E_binding (expt)    ~ -0.174 Ha")

    # Relative error (use binding energy for correlation check)
    E_binding_exp = -0.1745  # Approximate experimental binding energy
    if abs(E_binding_exp) > 1e-10:
        rel_error = abs((E_binding - E_binding_exp) / E_binding_exp)
    else:
        rel_error = 0.0
    energy_pass = rel_error < ENERGY_TOL
    print(f"       Relative error      = {rel_error:.4f} ({rel_error*100:.2f}%)")
    print(f"       Energy test:        {'PASS' if energy_pass else 'FAIL'} "
          f"(tolerance < {ENERGY_TOL*100:.1f}%)")

    # ----------------------------------------------------------------
    # 3. Bridge weight profile comparison
    # ----------------------------------------------------------------
    print("\n[3/3] Generating weight profile comparison...")

    R_profile = np.linspace(0.5, 5.0, 200)
    zeta = BRIDGE_DECAY_RATE

    W_exp = np.exp(-zeta * R_profile)
    zR = zeta * R_profile
    W_sto = (1.0 + zR + zR**2 / 3.0) * np.exp(-zR)

    # ----------------------------------------------------------------
    # Plot
    # ----------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        r"H$_2$ Bridge Weight Correction — STO Overlap Validation",
        fontsize=14, fontweight="bold",
    )

    # (a) PES
    ax = axes[0, 0]
    ax.plot(R_values, E_values, "o-", color="darkblue", markersize=3, linewidth=1.5,
            label="Full CI PES")
    ax.axvline(R_eq, color="red", linestyle="--", linewidth=1.2,
               label=f"$R_{{eq}}$ = {R_eq:.3f} Bohr")
    ax.axvline(R_EQ_EXP, color="gray", linestyle=":", linewidth=1.2,
               label=f"Expt. = {R_EQ_EXP:.3f} Bohr")
    ax.set_xlabel("R (Bohr)")
    ax.set_ylabel("Energy (Hartree)")
    ax.set_title("(a) Full CI Potential Energy Surface")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (b) Weight profiles
    ax = axes[0, 1]
    ax.plot(R_profile, W_exp, "-", color="red", linewidth=2,
            label=r"Old: $e^{-\zeta R}$")
    ax.plot(R_profile, W_sto, "-", color="blue", linewidth=2,
            label=r"New: $(1+\zeta R + (\zeta R)^2/3)\,e^{-\zeta R}$")
    ax.axvline(R_EQ_EXP, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("R (Bohr)")
    ax.set_ylabel("Bridge Weight W(R)")
    ax.set_title(r"(b) Bridge Weight: Exp vs STO Overlap ($\zeta$=1.0)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (c) Fractional gradient |dW/dR| / W
    ax = axes[1, 0]
    frac_exp = np.ones_like(R_profile) * zeta  # d/dR [exp(-zR)] / exp(-zR) = zeta
    frac_sto = zR * (1.0 + zR) / (3.0 * (1.0 + zR + zR**2 / 3.0))
    ax.plot(R_profile, frac_exp, "-", color="red", linewidth=2,
            label=r"Old: $\zeta$ (constant)")
    ax.plot(R_profile, frac_sto, "-", color="blue", linewidth=2,
            label=r"New: STO fractional decay")
    ax.axvline(R_EQ_EXP, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("R (Bohr)")
    ax.set_ylabel(r"$|dW/dR| \,/\, W$")
    ax.set_title("(c) Fractional Gradient (Compressive Bias)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (d) Summary table
    ax = axes[1, 1]
    ax.axis("off")
    table_data = [
        ["Metric", "Value", "Target", "Status"],
        ["R_eq (Bohr)", f"{R_eq:.4f}", f"{R_EQ_EXP:.3f} ± {R_EQ_TOL}",
         "PASS" if r_eq_pass else "FAIL"],
        ["E_min (Ha)", f"{E_min:.6f}", "—", "—"],
        ["E_binding (Ha)", f"{E_binding:.6f}", "~ -0.174", "—"],
        ["Rel. error", f"{rel_error*100:.2f}%", f"< {ENERGY_TOL*100:.1f}%",
         "PASS" if energy_pass else "FAIL"],
        ["Bridge form", "STO overlap", "—", "—"],
    ]
    colors = []
    for i, row in enumerate(table_data):
        if i == 0:
            colors.append(["lightsteelblue"] * 4)
        elif "FAIL" in row:
            colors.append(["mistyrose"] * 4)
        else:
            colors.append(["honeydew"] * 4)

    table = ax.table(cellText=table_data, cellColours=colors,
                     loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.5)
    ax.set_title("(d) Validation Summary", fontsize=12)

    plt.tight_layout()

    plot_dir = os.path.join(os.path.dirname(__file__), "plots")
    os.makedirs(plot_dir, exist_ok=True)
    plot_path = os.path.join(plot_dir, "bridge_correction_validation.png")
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"       Plot saved: {plot_path}")
    plt.close(fig)

    # ----------------------------------------------------------------
    # Final verdict
    # ----------------------------------------------------------------
    print("\n" + "=" * 70)
    all_pass = r_eq_pass and energy_pass
    if all_pass:
        print("  ALL TESTS PASSED")
    else:
        print("  SOME TESTS FAILED")
        if not r_eq_pass:
            print(f"    - R_eq = {R_eq:.4f} outside [{R_EQ_EXP - R_EQ_TOL:.3f}, "
                  f"{R_EQ_EXP + R_EQ_TOL:.3f}]")
        if not energy_pass:
            print(f"    - Relative error {rel_error*100:.2f}% exceeds {ENERGY_TOL*100:.1f}%")
    print("=" * 70)


if __name__ == "__main__":
    run_validation()
