#!/usr/bin/env python3
"""
H2 Potential Energy Surface (PES) -- Bond Dissociation Curve
=============================================================

Maps the Full CI ground state energy of H2 as a function of
internuclear distance R, producing the classic Morse-like potential
energy curve.

For each distance R:
  1. Build two hydrogen lattices with nucleus_position separated by R
  2. Stitch with distance-dependent bridge weights: W = exp(-R)
  3. Solve for the exact 2-electron ground state via Full CI
  4. Compute binding energy: D_e(R) = 2*E(H_atom) - E_CI(R)

In the graph framework, the CI electronic energy E_CI(R) directly
captures the binding landscape: the node weights (-Z/n^2) encode
the nuclear-electron attraction, while the Full CI tensor product
adds cross-nuclear attraction and electron-electron repulsion.
The balance of these terms produces the equilibrium geometry.

Usage:
    python benchmarks/scripts/h2_dissociation.py

Output:
    benchmarks/figures/h2_dissociation_curve.png
"""

import sys
import os
import time
import io

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from geovac import GeometricLattice, MoleculeHamiltonian


def compute_h2_energy_ci(
    R: float,
    max_n: int,
    n_bridges: int,
    kinetic_scale: float,
    bridge_decay_rate: float,
) -> dict:
    """Compute H2 Full CI energy at internuclear distance R.

    Parameters
    ----------
    R : float
        Internuclear distance in Bohr.
    max_n : int
        Principal quantum number cutoff per atom.
    n_bridges : int
        Number of topological bridge edges.
    kinetic_scale : float
        Universal kinetic scale constant.
    bridge_decay_rate : float
        Exponential decay rate for bridge weights (1/Bohr).

    Returns
    -------
    dict with keys: R, E_ci, V_nn, E_total, bridge_weight
    """
    atom_A = GeometricLattice(
        max_n=max_n,
        nucleus_position=(0.0, 0.0, 0.0),
        nuclear_charge=1,
    )
    atom_B = GeometricLattice(
        max_n=max_n,
        nucleus_position=(R, 0.0, 0.0),
        nuclear_charge=1,
    )

    h2 = MoleculeHamiltonian(
        lattices=[atom_A, atom_B],
        connectivity=[(0, 1, n_bridges)],
        kinetic_scale=kinetic_scale,
        bridge_decay_rate=bridge_decay_rate,
    )

    # Suppress Full CI solver's verbose output
    _stdout = sys.stdout
    sys.stdout = io.StringIO()
    try:
        E_ci, _ = h2.compute_ground_state(n_states=1, method="full_ci")
    finally:
        sys.stdout = _stdout

    V_nn = 1.0 / R
    bw = h2.bridge_info[0]["bridge_weight"] if h2.bridge_info else 0.0

    return {
        "R": R,
        "E_ci": E_ci[0],
        "V_nn": V_nn,
        "E_total": E_ci[0] + V_nn,
        "bridge_weight": bw,
    }


def run_dissociation_curve() -> None:
    """Compute and plot the H2 dissociation curve using Full CI."""

    # Parameters -- max_n=4 gives 30 states/atom, 3600-dim CI space (~0.03s/point)
    max_n = 4
    n_bridges = 20
    kinetic_scale = -1 / 16
    bridge_decay_rate = 1.0
    HA_TO_EV = 27.2114

    # Dense grid near equilibrium, coarser outward
    R_values = np.array([
        0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
        1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8,
        2.0, 2.5, 3.0, 4.0, 5.0, 6.0,
    ])

    print("=" * 72)
    print("  H2 Potential Energy Surface -- Full CI Dissociation Curve")
    print("=" * 72)
    print(f"  Method:         Full CI (exact 2-electron)")
    print(f"  Basis:          max_n = {max_n} per atom (30 states, 3600-dim CI)")
    print(f"  Bridges:        {n_bridges}")
    print(f"  Kinetic scale:  {kinetic_scale}")
    print(f"  Bridge decay:   W = exp(-{bridge_decay_rate} * R)")
    print(f"  R range:        [{R_values[0]:.1f}, {R_values[-1]:.1f}] Bohr")
    print(f"  Points:         {len(R_values)}")
    print()

    # Isolated atom reference
    atom_iso = GeometricLattice(max_n=max_n, nuclear_charge=1)
    h_iso = MoleculeHamiltonian(
        lattices=[atom_iso], connectivity=[], kinetic_scale=kinetic_scale,
    )
    E_atom, _ = eigsh(h_iso.hamiltonian, k=1, which="SA")
    E_2H = 2.0 * E_atom[0]
    print(f"  Isolated H atom:    E = {E_atom[0]:.6f} Ha")
    print(f"  Dissociation limit: 2E(H) = {E_2H:.6f} Ha")
    print()

    # Sweep
    print(f"  {'R':>6} {'E_CI':>14} {'V_nn':>10} {'E_total':>14} "
          f"{'D_e (eV)':>10} {'W_bridge':>10}")
    print(f"  {'-'*6} {'-'*14} {'-'*10} {'-'*14} {'-'*10} {'-'*10}")

    results = []
    t0 = time.time()
    for R in R_values:
        res = compute_h2_energy_ci(
            R, max_n, n_bridges, kinetic_scale, bridge_decay_rate,
        )
        results.append(res)
        # Binding energy: positive = bound (CI energy below 2H limit)
        D_e_ci = E_2H - res["E_ci"]
        print(
            f"  {res['R']:>6.2f} {res['E_ci']:>14.6f} {res['V_nn']:>10.6f} "
            f"{res['E_total']:>14.6f} {D_e_ci * HA_TO_EV:>10.3f} "
            f"{res['bridge_weight']:>10.6f}"
        )

    elapsed = time.time() - t0
    print()
    print(f"  Sweep time: {elapsed:.2f}s  ({elapsed / len(R_values):.2f}s per point)")
    print()

    # Extract arrays
    R_arr = np.array([r["R"] for r in results])
    E_ci_arr = np.array([r["E_ci"] for r in results])
    V_nn_arr = np.array([r["V_nn"] for r in results])
    E_total_arr = np.array([r["E_total"] for r in results])

    # Binding energy curve: D_e(R) = 2*E(H) - E_CI(R)
    # Positive D_e means the CI energy is below the 2H limit (bound)
    D_e_arr = E_2H - E_ci_arr

    # Find equilibrium from E_CI minimum (most negative = most bound)
    idx_min_ci = np.argmin(E_ci_arr)
    R_eq_ci = R_arr[idx_min_ci]
    E_min_ci = E_ci_arr[idx_min_ci]
    D_e_ci = E_2H - E_min_ci

    print("  CI Equilibrium (E_CI minimum):")
    print(f"    R_eq    = {R_eq_ci:.2f} Bohr")
    print(f"    E_CI    = {E_min_ci:.6f} Ha")
    print(f"    D_e     = {D_e_ci:.6f} Ha ({D_e_ci * HA_TO_EV:.3f} eV)")
    print()

    # Asymptotic behavior
    E_asym = E_ci_arr[-1]
    asym_diff = abs(E_asym - E_2H)
    print("  Asymptotic behavior:")
    print(f"    E_CI(R={R_arr[-1]:.0f}) = {E_asym:.6f} Ha")
    print(f"    2E(H)         = {E_2H:.6f} Ha")
    print(f"    Difference:     {asym_diff:.6f} Ha ({asym_diff / abs(E_2H) * 100:.2f}%)")
    print()

    # Experimental comparison
    R_exp_eq = 1.401   # Bohr
    D_exp = 0.1745     # Ha (4.75 eV)
    print("  Experimental comparison:")
    print(f"    {'':>16} {'GeoVac':>12} {'Experiment':>12}")
    print(f"    {'R_eq (Bohr)':>16} {R_eq_ci:>12.2f} {R_exp_eq:>12.3f}")
    print(f"    {'D_e (eV)':>16} {D_e_ci * HA_TO_EV:>12.3f} {D_exp * HA_TO_EV:>12.3f}")
    print()

    # ------------------------------------------------------------------
    # Plot: two panels
    # ------------------------------------------------------------------
    print("Generating plot...")

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 10), height_ratios=[1.2, 1],
    )
    fig.suptitle(
        r"H$_2$ Potential Energy Surface — Full CI Dissociation Curve"
        f"\n(max_n={max_n}, {n_bridges} bridges, "
        r"$W_{\mathrm{bridge}} = e^{-R}$"
        ")",
        fontsize=13, fontweight="bold",
    )

    # --- Top panel: Energy vs R ---
    ax1.plot(
        R_arr, E_ci_arr,
        "o-", color="darkblue", linewidth=2, markersize=5,
        label=r"$E_{\mathrm{CI}}(R)$ (2-electron, Full CI)",
    )
    ax1.plot(
        R_arr, E_total_arr,
        "s--", color="steelblue", linewidth=1.2, markersize=4, alpha=0.7,
        label=r"$E_{\mathrm{CI}}(R) + 1/R$ (with nuclear repulsion)",
    )
    ax1.axhline(
        E_2H, color="red", linestyle=":", linewidth=1.5, alpha=0.7,
        label=f"2E(H) = {E_2H:.4f} Ha (dissociation limit)",
    )
    ax1.plot(
        R_eq_ci, E_min_ci, "*", color="gold", markersize=15, zorder=5,
        markeredgecolor="black", markeredgewidth=0.5,
        label=f"CI minimum: R={R_eq_ci:.2f} Bohr",
    )
    ax1.axvline(
        R_exp_eq, color="gray", linestyle="--", linewidth=1, alpha=0.5,
        label=f"Expt. R_eq = {R_exp_eq:.3f} Bohr",
    )

    ax1.set_xlabel("Internuclear Distance R (Bohr)", fontsize=11)
    ax1.set_ylabel("Energy (Hartree)", fontsize=11)
    ax1.set_title("Energy Curves", fontsize=11)
    ax1.legend(loc="lower right", fontsize=8)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, R_arr[-1] + 0.3)
    y_pad = 0.15
    ax1.set_ylim(
        min(E_ci_arr.min(), E_total_arr.min()) - y_pad,
        max(E_2H, E_total_arr.max()) + y_pad,
    )

    # --- Bottom panel: Binding energy D_e(R) ---
    ax2.plot(
        R_arr, D_e_arr * HA_TO_EV,
        "o-", color="darkgreen", linewidth=2, markersize=5,
        label=r"$D_e(R) = 2E(\mathrm{H}) - E_{\mathrm{CI}}(R)$",
    )
    ax2.axhline(
        0, color="black", linestyle="-", linewidth=0.8, alpha=0.5,
    )
    ax2.axhline(
        D_exp * HA_TO_EV, color="red", linestyle=":", linewidth=1.5, alpha=0.7,
        label=f"Expt. D_e = {D_exp * HA_TO_EV:.2f} eV",
    )
    ax2.plot(
        R_eq_ci, D_e_ci * HA_TO_EV, "*", color="gold", markersize=15, zorder=5,
        markeredgecolor="black", markeredgewidth=0.5,
        label=f"Max binding: {D_e_ci * HA_TO_EV:.2f} eV at R={R_eq_ci:.2f}",
    )
    ax2.axvline(
        R_exp_eq, color="gray", linestyle="--", linewidth=1, alpha=0.5,
    )

    ax2.set_xlabel("Internuclear Distance R (Bohr)", fontsize=11)
    ax2.set_ylabel("Binding Energy (eV)", fontsize=11)
    ax2.set_title("Binding Energy Curve", fontsize=11)
    ax2.legend(loc="upper right", fontsize=8)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, R_arr[-1] + 0.3)

    plt.tight_layout()

    # Save
    fig_dir = os.path.join(os.path.dirname(__file__), "..", "figures")
    os.makedirs(fig_dir, exist_ok=True)
    plot_path = os.path.join(fig_dir, "h2_dissociation_curve.png")
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved to: {plot_path}")
    plt.close(fig)

    print()
    print("=" * 72)
    print(f"  PES complete. {len(R_values)} points in {elapsed:.2f}s.")
    print(f"  Equilibrium: R = {R_eq_ci:.2f} Bohr")
    print(f"  Binding:     D_e = {D_e_ci * HA_TO_EV:.3f} eV")
    print("=" * 72)


if __name__ == "__main__":
    run_dissociation_curve()
