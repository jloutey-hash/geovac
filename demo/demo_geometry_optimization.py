#!/usr/bin/env python3
"""
H2 Automatic Geometry Optimization
====================================

Finds the equilibrium bond length of H2 by gradient descent on the
Full CI potential energy surface.

At each step:
  1. Compute E(R) via Full CI
  2. Compute E(R + dR) for numerical gradient
  3. Force = -dE/dR
  4. Update: R_new = R + alpha * Force
  5. Converge when |Force| < tolerance

Starting far from equilibrium (R_init = 2.5 Bohr), the optimizer
should converge to the PES minimum at R ~ 1.30 Bohr identified
in the dissociation curve benchmark.

Usage:
    python demo/demo_geometry_optimization.py

Output:
    debug/plots/h2_geometry_optimization.png
"""

import sys
import os
import io
import time

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from geovac import GeometricLattice, MoleculeHamiltonian


# -- Molecular parameters (fixed throughout optimization) --
MAX_N = 4
N_BRIDGES = 20
KINETIC_SCALE = -1 / 16
BRIDGE_DECAY_RATE = 1.0


def compute_energy(R: float) -> float:
    """Compute H2 Full CI electronic energy at internuclear distance R.

    The Full CI tensor product space includes cross-nuclear attraction
    and electron-electron repulsion, so E_CI(R) directly captures the
    binding landscape with a well-defined minimum.
    """
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

    # Suppress Full CI verbose output
    _stdout = sys.stdout
    sys.stdout = io.StringIO()
    try:
        E_ci, _ = h2.compute_ground_state(n_states=1, method="full_ci")
    finally:
        sys.stdout = _stdout

    return E_ci[0]


def run_geometry_optimization(
    R_init: float = 2.5,
    alpha: float = 0.5,
    tol: float = 1e-4,
    dR: float = 0.001,
    max_steps: int = 200,
    R_min_clamp: float = 0.4,
    R_max_clamp: float = 10.0,
) -> None:
    """Optimize H2 bond length by gradient descent on the Full CI PES.

    Parameters
    ----------
    R_init : float
        Starting internuclear distance (Bohr).
    alpha : float
        Learning rate (Bohr^2 / Hartree).
    tol : float
        Convergence tolerance on |Force| (Hartree/Bohr).
    dR : float
        Finite difference step for numerical gradient.
    max_steps : int
        Maximum optimization steps.
    R_min_clamp : float
        Minimum allowed R to prevent singularities.
    R_max_clamp : float
        Maximum allowed R.
    """
    HA_TO_EV = 27.2114

    print("=" * 65)
    print("  H2 Geometry Optimization -- Gradient Descent on Full CI PES")
    print("=" * 65)
    print(f"  Basis:        max_n = {MAX_N} per atom")
    print(f"  Bridges:      {N_BRIDGES}")
    print(f"  R_init:       {R_init:.2f} Bohr")
    print(f"  Learning rate: {alpha}")
    print(f"  Tolerance:    |F| < {tol} Ha/Bohr")
    print(f"  Finite diff:  dR = {dR}")
    print()

    # Optimization loop
    R = R_init
    history = []

    print(f"  {'Step':>4} {'R (Bohr)':>10} {'E_total (Ha)':>14} "
          f"{'Force (Ha/B)':>14} {'|dR|':>10}")
    print(f"  {'-'*4} {'-'*10} {'-'*14} {'-'*14} {'-'*10}")

    t0 = time.time()
    converged = False

    for step in range(max_steps):
        # Energy at current R and perturbed R
        E = compute_energy(R)
        E_plus = compute_energy(R + dR)

        # Numerical force: F = -dE/dR
        gradient = (E_plus - E) / dR
        force = -gradient

        # Store history
        history.append({"step": step, "R": R, "E": E, "force": force})

        # Compute step size
        delta_R = alpha * force
        step_size = abs(delta_R)

        # Print every step (fast enough at ~0.06s/step)
        print(
            f"  {step:>4} {R:>10.5f} {E:>14.6f} "
            f"{force:>14.6f} {step_size:>10.6f}"
        )

        # Check convergence
        if abs(force) < tol:
            converged = True
            break

        # Update R with clamping
        R_new = R + delta_R
        R = np.clip(R_new, R_min_clamp, R_max_clamp)

    elapsed = time.time() - t0
    n_steps = len(history)

    print()
    if converged:
        print(f"  CONVERGED in {n_steps} steps ({elapsed:.2f}s)")
    else:
        print(f"  Reached max steps ({max_steps}), elapsed {elapsed:.2f}s")
    print()

    # Final state
    final = history[-1]
    print("  Optimization Result:")
    print(f"    R_init      = {R_init:.4f} Bohr")
    print(f"    R_final     = {final['R']:.5f} Bohr")
    print(f"    E_final     = {final['E']:.6f} Ha")
    print(f"    Final Force  = {final['force']:.2e} Ha/Bohr")
    print(f"    Steps        = {n_steps}")
    print(f"    Time         = {elapsed:.2f}s ({elapsed / n_steps:.2f}s per step)")
    print()

    # Comparison
    R_exp = 1.401
    R_pes = 1.40  # Target from STO-overlap bridge correction
    print("  Comparison:")
    print(f"    {'':>20} {'Optimizer':>12} {'PES Sweep':>12} {'Experiment':>12}")
    print(f"    {'R_eq (Bohr)':>20} {final['R']:>12.4f} {R_pes:>12.2f} {R_exp:>12.3f}")
    print()

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    print("Generating plot...")

    steps_arr = np.array([h["step"] for h in history])
    R_arr = np.array([h["R"] for h in history])
    E_arr = np.array([h["E"] for h in history])
    F_arr = np.array([h["force"] for h in history])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), height_ratios=[1, 1])
    fig.suptitle(
        r"H$_2$ Geometry Optimization — Gradient Descent on Full CI PES",
        fontsize=13, fontweight="bold",
    )

    # --- Top: Energy vs Step ---
    ax1.plot(steps_arr, E_arr, "o-", color="darkblue", linewidth=2, markersize=4)
    ax1.axhline(
        E_arr[-1], color="red", linestyle=":", linewidth=1, alpha=0.6,
        label=f"E_final = {E_arr[-1]:.6f} Ha",
    )
    ax1.set_xlabel("Optimization Step")
    ax1.set_ylabel("Total Energy (Hartree)")
    ax1.set_title("Energy Convergence")
    ax1.legend(loc="upper right", fontsize=9)
    ax1.grid(True, alpha=0.3)

    # --- Bottom: R vs Step ---
    ax2.plot(steps_arr, R_arr, "o-", color="darkgreen", linewidth=2, markersize=4)
    ax2.axhline(
        R_arr[-1], color="red", linestyle=":", linewidth=1, alpha=0.6,
        label=f"R_final = {R_arr[-1]:.4f} Bohr",
    )
    ax2.axhline(
        R_exp, color="gray", linestyle="--", linewidth=1, alpha=0.5,
        label=f"Expt. R_eq = {R_exp:.3f} Bohr",
    )
    ax2.set_xlabel("Optimization Step")
    ax2.set_ylabel("Internuclear Distance R (Bohr)")
    ax2.set_title("Bond Length Convergence")
    ax2.legend(loc="upper right", fontsize=9)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save
    plot_dir = os.path.join(os.path.dirname(__file__), "..", "debug", "plots")
    os.makedirs(plot_dir, exist_ok=True)
    plot_path = os.path.join(plot_dir, "h2_geometry_optimization.png")
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved to: {plot_path}")
    plt.close(fig)

    print()
    print("=" * 65)
    print(f"  Optimization complete: R = {final['R']:.5f} Bohr in {n_steps} steps")
    print("=" * 65)


if __name__ == "__main__":
    run_geometry_optimization()
