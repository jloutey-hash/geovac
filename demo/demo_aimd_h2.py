#!/usr/bin/env python3
"""
H2 Ab Initio Molecular Dynamics (AIMD)
========================================

Simulates real-time nuclear vibrations of H2 by coupling the exact
quantum electronic energy (Full CI) to classical nuclear kinematics
via the Velocity Verlet integrator.

At each step:
  1. Compute the Full CI electronic energy E_CI(R) as the potential
  2. Evaluate the force F = -dE_CI/dR by finite difference
  3. Advance nuclei with Velocity Verlet (symplectic, time-reversible)

Starting from a compressed geometry (R_0 = 1.0 Bohr, v_0 = 0),
the molecule oscillates around its equilibrium R_eq ~ 1.29 Bohr.
Total energy (electronic + nuclear kinetic) is conserved by the
symplectic integrator.

Usage:
    python demo/demo_aimd_h2.py

Output:
    debug/plots/h2_aimd.png
"""

import sys
import os
import io
import time

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from geovac import GeometricLattice, MoleculeHamiltonian


# -- Molecular parameters (same as geometry optimizer) --
MAX_N = 4
N_BRIDGES = 20
KINETIC_SCALE = -1 / 16
BRIDGE_DECAY_RATE = 1.0


def compute_energy(R: float) -> float:
    """Compute H2 Full CI electronic energy at internuclear distance R.

    Returns E_CI(R) which serves as the potential energy surface
    for nuclear dynamics.
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

    _stdout = sys.stdout
    sys.stdout = io.StringIO()
    try:
        E_ci, _ = h2.compute_ground_state(n_states=1, method="full_ci")
    finally:
        sys.stdout = _stdout

    return E_ci[0]


def compute_force(R: float, dR: float = 0.001) -> tuple:
    """Compute energy and force at internuclear distance R.

    Returns (E, F) where F = -dE/dR via central difference.
    """
    E_minus = compute_energy(R - dR)
    E_plus = compute_energy(R + dR)
    E_center = compute_energy(R)
    force = -(E_plus - E_minus) / (2 * dR)
    return E_center, force


def run_aimd(
    R_init: float = 1.0,
    v_init: float = 0.0,
    dt: float = 1.0,
    n_steps: int = 600,
    dR: float = 0.001,
    R_min_clamp: float = 0.3,
) -> None:
    """Run Ab Initio Molecular Dynamics for H2.

    Parameters
    ----------
    R_init : float
        Initial internuclear distance (Bohr). Compressed from R_eq ~ 1.29.
    v_init : float
        Initial relative velocity of nuclei (Bohr / a.u. time).
    dt : float
        Nuclear time step (atomic units, ~0.024 fs).
    n_steps : int
        Number of Velocity Verlet steps.
    dR : float
        Finite difference step for force evaluation.
    R_min_clamp : float
        Minimum allowed R to prevent singularities.
    """
    # Reduced mass of H2 in atomic units (m_proton / 2)
    mu = 1836.15 / 2.0
    HA_TO_EV = 27.2114
    AU_TO_FS = 0.02419  # 1 a.u. of time ~ 0.024 fs

    T_total = dt * n_steps
    T_total_fs = T_total * AU_TO_FS

    print("=" * 65)
    print("  H2 Ab Initio Molecular Dynamics -- Velocity Verlet")
    print("=" * 65)
    print(f"  Basis:         max_n = {MAX_N} per atom")
    print(f"  Bridges:       {N_BRIDGES}")
    print(f"  Reduced mass:  mu = {mu:.2f} a.u. (m_p / 2)")
    print(f"  R_init:        {R_init:.2f} Bohr (compressed from R_eq ~ 1.29)")
    print(f"  v_init:        {v_init} Bohr/a.u.")
    print(f"  dt:            {dt:.1f} a.u. ({dt * AU_TO_FS:.4f} fs)")
    print(f"  Steps:         {n_steps}")
    print(f"  Total time:    {T_total:.0f} a.u. ({T_total_fs:.2f} fs)")
    print()

    # ------------------------------------------------------------------
    # Velocity Verlet integration
    # ------------------------------------------------------------------
    print("Running AIMD simulation...")
    print()
    print(f"  {'Step':>4} {'t (a.u.)':>10} {'R (Bohr)':>10} "
          f"{'E_pot (Ha)':>12} {'E_kin (Ha)':>12} {'E_tot (Ha)':>12} "
          f"{'Force':>12}")
    print(f"  {'-'*4} {'-'*10} {'-'*10} "
          f"{'-'*12} {'-'*12} {'-'*12} {'-'*12}")

    t0 = time.time()

    R = R_init
    v = v_init

    # Initial force
    E_pot, F = compute_force(R, dR)
    E_kin = 0.5 * mu * v**2
    E_tot = E_pot + E_kin

    # Storage
    history = [{
        "step": 0, "t": 0.0, "R": R, "v": v,
        "E_pot": E_pot, "E_kin": E_kin, "E_tot": E_tot, "force": F,
    }]

    print(
        f"  {0:>4} {0.0:>10.1f} {R:>10.5f} "
        f"{E_pot:>12.6f} {E_kin:>12.6f} {E_tot:>12.6f} "
        f"{F:>12.6f}"
    )

    for step in range(1, n_steps + 1):
        # Velocity Verlet: half-step velocity
        a = F / mu
        v_half = v + 0.5 * a * dt

        # Full-step position
        R_new = R + v_half * dt

        # Clamp to prevent singularity
        if R_new < R_min_clamp:
            R_new = R_min_clamp
            v_half = 0.0  # Reflect

        # Force at new position
        E_pot_new, F_new = compute_force(R_new, dR)

        # Complete velocity step
        a_new = F_new / mu
        v_new = v_half + 0.5 * a_new * dt

        # Energies
        E_kin_new = 0.5 * mu * v_new**2
        E_tot_new = E_pot_new + E_kin_new

        # Store
        t_now = step * dt
        history.append({
            "step": step, "t": t_now, "R": R_new, "v": v_new,
            "E_pot": E_pot_new, "E_kin": E_kin_new, "E_tot": E_tot_new,
            "force": F_new,
        })

        # Print every 10 steps + first 5 + last step
        if step <= 5 or step % 10 == 0 or step == n_steps:
            print(
                f"  {step:>4} {t_now:>10.1f} {R_new:>10.5f} "
                f"{E_pot_new:>12.6f} {E_kin_new:>12.6f} {E_tot_new:>12.6f} "
                f"{F_new:>12.6f}"
            )

        # Advance state
        R = R_new
        v = v_new
        F = F_new

    elapsed = time.time() - t0
    print()
    print(f"  Simulation complete: {elapsed:.2f}s "
          f"({elapsed / n_steps:.2f}s per step)")
    print()

    # ------------------------------------------------------------------
    # Analysis
    # ------------------------------------------------------------------
    t_arr = np.array([h["t"] for h in history])
    R_arr = np.array([h["R"] for h in history])
    E_pot_arr = np.array([h["E_pot"] for h in history])
    E_kin_arr = np.array([h["E_kin"] for h in history])
    E_tot_arr = np.array([h["E_tot"] for h in history])

    R_min = np.min(R_arr)
    R_max = np.max(R_arr)
    R_mean = np.mean(R_arr)

    E_tot_mean = np.mean(E_tot_arr)
    E_tot_drift = np.max(np.abs(E_tot_arr - E_tot_arr[0]))
    E_tot_rel_drift = E_tot_drift / abs(E_tot_arr[0]) * 100

    # Estimate vibration period from R(t) oscillation
    # Find zero-crossings of R(t) - R_mean (ascending)
    R_centered = R_arr - R_mean
    crossings = []
    for i in range(1, len(R_centered)):
        if R_centered[i - 1] < 0 and R_centered[i] >= 0:
            # Linear interpolation for crossing time
            frac = -R_centered[i - 1] / (R_centered[i] - R_centered[i - 1])
            t_cross = t_arr[i - 1] + frac * dt
            crossings.append(t_cross)

    if len(crossings) >= 2:
        periods = np.diff(crossings)
        T_vib = np.mean(periods)
        T_vib_fs = T_vib * AU_TO_FS
        freq_cm1 = 1.0 / (T_vib_fs * 1e-15 * 2.998e10) if T_vib_fs > 0 else 0
    else:
        T_vib = 0.0
        T_vib_fs = 0.0
        freq_cm1 = 0.0

    print("  AIMD Results:")
    print(f"    R_min       = {R_min:.5f} Bohr")
    print(f"    R_max       = {R_max:.5f} Bohr")
    print(f"    R_mean      = {R_mean:.5f} Bohr")
    print(f"    Amplitude   = {(R_max - R_min) / 2:.5f} Bohr")
    print()
    print("  Energy Conservation:")
    print(f"    E_tot(0)    = {E_tot_arr[0]:.8f} Ha")
    print(f"    E_tot(end)  = {E_tot_arr[-1]:.8f} Ha")
    print(f"    Max drift   = {E_tot_drift:.2e} Ha ({E_tot_rel_drift:.4f}%)")
    print()
    print("  Vibration:")
    if T_vib > 0:
        print(f"    Period      = {T_vib:.1f} a.u. ({T_vib_fs:.2f} fs)")
        print(f"    Frequency   = {freq_cm1:.0f} cm^-1")
        print(f"    Crossings   = {len(crossings)}")
        print(f"    Expt. H2:   T ~ 8.1 fs, v ~ 4161 cm^-1")
    else:
        print("    Less than 2 zero-crossings detected; "
              "increase n_steps for period estimate.")
    print()

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    print("Generating plot...")

    t_fs = t_arr * AU_TO_FS

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), height_ratios=[1, 1])
    fig.suptitle(
        r"H$_2$ Ab Initio Molecular Dynamics — Velocity Verlet on Full CI PES",
        fontsize=13, fontweight="bold",
    )

    # --- Top: Energy vs Time ---
    ax1.plot(t_fs, E_tot_arr, color="black", linewidth=1.5,
             label=f"Total Energy (drift = {E_tot_drift:.2e} Ha)")
    ax1.plot(t_fs, E_pot_arr, color="darkblue", linewidth=1.0, alpha=0.8,
             label="Potential Energy (E_CI)")
    ax1.plot(t_fs, E_kin_arr, color="red", linewidth=1.0, alpha=0.8,
             label="Kinetic Energy (nuclear)")
    ax1.set_xlabel("Time (fs)")
    ax1.set_ylabel("Energy (Hartree)")
    ax1.set_title("Energy Conservation")
    ax1.legend(loc="right", fontsize=8)
    ax1.grid(True, alpha=0.3)

    # --- Bottom: R vs Time ---
    ax2.plot(t_fs, R_arr, color="darkgreen", linewidth=1.5)
    ax2.axhline(
        R_mean, color="gray", linestyle="--", linewidth=1, alpha=0.6,
        label=f"R_mean = {R_mean:.4f} Bohr",
    )
    ax2.axhline(
        1.293, color="red", linestyle=":", linewidth=1, alpha=0.5,
        label="R_eq = 1.293 Bohr (optimizer)",
    )
    ax2.set_xlabel("Time (fs)")
    ax2.set_ylabel("Internuclear Distance R (Bohr)")
    ax2.set_title("Nuclear Vibration")
    ax2.legend(loc="upper right", fontsize=8)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save
    plot_dir = os.path.join(os.path.dirname(__file__), "..", "debug", "plots")
    os.makedirs(plot_dir, exist_ok=True)
    plot_path = os.path.join(plot_dir, "h2_aimd.png")
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved to: {plot_path}")
    plt.close(fig)

    print()
    print("=" * 65)
    print(f"  AIMD complete: {n_steps} steps in {elapsed:.2f}s")
    print(f"  Vibration: R = [{R_min:.4f}, {R_max:.4f}] Bohr "
          f"around {R_mean:.4f}")
    if T_vib > 0:
        print(f"  Period: {T_vib_fs:.2f} fs ({freq_cm1:.0f} cm^-1)")
    print(f"  Energy conservation: {E_tot_drift:.2e} Ha drift")
    print("=" * 65)


if __name__ == "__main__":
    run_aimd()
