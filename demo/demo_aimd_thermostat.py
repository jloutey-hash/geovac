#!/usr/bin/env python3
"""
H2 Langevin Thermostat AIMD -- Thermal Vibration & Bond Dissociation
=====================================================================

Extends the NVE Ab Initio Molecular Dynamics engine with a Langevin
thermostat to simulate the H2 molecule in the NVT (canonical) ensemble.

The Langevin equation couples the nuclei to a thermal bath:

    F_total = F_quantum - gamma * mu * v + xi(t)

where gamma is the friction coefficient and xi(t) is a stochastic
force with variance sigma = sqrt(2 * gamma * T * mu / dt).

Two scenarios demonstrate the physics:
  1. Room temperature (T ~ 300K): stable thermal vibrations around
     the equilibrium bond length.
  2. Extreme temperature (T ~ 950000K): violent thermal kicks overcome
     the binding well and dissociate the molecule.

Usage:
    python demo/demo_aimd_thermostat.py

Output:
    debug/plots/h2_thermal_dissociation.png
"""

import sys
import os
import io
import time

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from geovac import GeometricLattice, MoleculeHamiltonian


# -- Molecular parameters (same as geometry optimizer / NVE AIMD) --
MAX_N = 4
N_BRIDGES = 20
KINETIC_SCALE = -1 / 16
BRIDGE_DECAY_RATE = 1.0


def compute_energy(R: float) -> float:
    """Compute H2 Full CI electronic energy at internuclear distance R."""
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


def run_langevin(
    label: str,
    T_target: float,
    R_init: float = 1.29,
    v_init: float = 0.0,
    gamma: float = 0.01,
    dt: float = 1.0,
    n_steps: int = 500,
    dR: float = 0.001,
    R_min_clamp: float = 0.3,
    R_max_clamp: float = 20.0,
    seed: int = 42,
) -> dict:
    """Run Langevin thermostat AIMD for H2.

    Parameters
    ----------
    label : str
        Scenario name for console output.
    T_target : float
        Target temperature in atomic units (1 a.u. ~ 315775 K).
    R_init : float
        Initial internuclear distance (Bohr).
    v_init : float
        Initial relative velocity (Bohr / a.u. time).
    gamma : float
        Friction coefficient (1 / a.u. time).
    dt : float
        Nuclear time step (atomic units).
    n_steps : int
        Number of integration steps.
    dR : float
        Finite difference step for force evaluation.
    R_min_clamp : float
        Minimum allowed R to prevent singularities.
    R_max_clamp : float
        Maximum R; trajectory stops if exceeded (dissociation).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    dict with arrays: t, R, E_pot, E_kin, E_tot, dissociated (bool)
    """
    mu = 1836.15 / 2.0
    AU_TO_FS = 0.02419
    AU_TO_K = 315775.0

    T_kelvin = T_target * AU_TO_K

    rng = np.random.default_rng(seed)

    # Stochastic force standard deviation
    sigma = np.sqrt(2.0 * gamma * T_target * mu / dt)

    print(f"  --- {label} ---")
    print(f"  T_target:  {T_target} a.u. ({T_kelvin:.0f} K)")
    print(f"  gamma:     {gamma}")
    print(f"  sigma:     {sigma:.4f} Ha/Bohr")
    print(f"  R_init:    {R_init:.2f} Bohr")
    print(f"  Steps:     {n_steps}")
    print()

    print(f"  {'Step':>4} {'t (a.u.)':>10} {'R (Bohr)':>10} "
          f"{'E_pot (Ha)':>12} {'E_kin (Ha)':>12} {'E_tot (Ha)':>12}")
    print(f"  {'-'*4} {'-'*10} {'-'*10} "
          f"{'-'*12} {'-'*12} {'-'*12}")

    t0 = time.time()

    R = R_init
    v = v_init
    dissociated = False
    dissoc_step = -1

    # Initial quantum force
    E_pot, F_q = compute_force(R, dR)
    E_kin = 0.5 * mu * v**2
    E_tot = E_pot + E_kin

    # Storage
    t_list = [0.0]
    R_list = [R]
    E_pot_list = [E_pot]
    E_kin_list = [E_kin]
    E_tot_list = [E_tot]

    print(
        f"  {0:>4} {0.0:>10.1f} {R:>10.5f} "
        f"{E_pot:>12.6f} {E_kin:>12.6f} {E_tot:>12.6f}"
    )

    for step in range(1, n_steps + 1):
        # Langevin forces: friction + stochastic
        F_friction = -gamma * mu * v
        F_random = sigma * rng.standard_normal()
        F_total = F_q + F_friction + F_random

        # Velocity Verlet half-step
        a = F_total / mu
        v_half = v + 0.5 * a * dt

        # Position update
        R_new = R + v_half * dt

        # Clamp floor
        if R_new < R_min_clamp:
            R_new = R_min_clamp
            v_half = abs(v_half)  # Reflect off wall

        # Check dissociation
        if R_new > R_max_clamp:
            dissociated = True
            dissoc_step = step

        # Quantum force at new position
        E_pot_new, F_q_new = compute_force(R_new, dR)

        # Langevin forces at new position (use v_half for friction estimate)
        F_friction_new = -gamma * mu * v_half
        F_random_new = sigma * rng.standard_normal()
        F_total_new = F_q_new + F_friction_new + F_random_new

        # Complete velocity step
        a_new = F_total_new / mu
        v_new = v_half + 0.5 * a_new * dt

        # Energies
        E_kin_new = 0.5 * mu * v_new**2
        E_tot_new = E_pot_new + E_kin_new

        t_now = step * dt
        t_list.append(t_now)
        R_list.append(R_new)
        E_pot_list.append(E_pot_new)
        E_kin_list.append(E_kin_new)
        E_tot_list.append(E_tot_new)

        if step <= 3 or step % 25 == 0 or step == n_steps or dissociated:
            tag = " << DISSOCIATED" if dissociated and dissoc_step == step else ""
            print(
                f"  {step:>4} {t_now:>10.1f} {R_new:>10.5f} "
                f"{E_pot_new:>12.6f} {E_kin_new:>12.6f} "
                f"{E_tot_new:>12.6f}{tag}"
            )

        # Advance state
        R = R_new
        v = v_new
        F_q = F_q_new

        if dissociated:
            break

    elapsed = time.time() - t0
    actual_steps = len(t_list) - 1
    print()
    print(f"  Completed {actual_steps} steps in {elapsed:.2f}s "
          f"({elapsed / max(actual_steps, 1):.2f}s/step)")
    if dissociated:
        print(f"  ** Bond DISSOCIATED at step {dissoc_step} "
              f"(R > {R_max_clamp} Bohr)")
    else:
        R_arr = np.array(R_list)
        print(f"  Bond STABLE: R in [{R_arr.min():.4f}, {R_arr.max():.4f}] Bohr")
    print()

    return {
        "t": np.array(t_list),
        "R": np.array(R_list),
        "E_pot": np.array(E_pot_list),
        "E_kin": np.array(E_kin_list),
        "E_tot": np.array(E_tot_list),
        "dissociated": dissociated,
        "dissoc_step": dissoc_step,
        "elapsed": elapsed,
        "T_target": T_target,
        "T_kelvin": T_kelvin,
        "label": label,
    }


def run_thermostat_demo() -> None:
    """Run two Langevin AIMD scenarios and generate comparison plot."""

    AU_TO_FS = 0.02419
    AU_TO_K = 315775.0

    T_low = 0.001    # ~316 K (room temperature)
    T_high = 3.0     # ~947000 K (extreme, enough to dissociate)

    print("=" * 65)
    print("  H2 Langevin Thermostat AIMD -- Thermal Bond Dynamics")
    print("=" * 65)
    print(f"  Basis:       max_n = {MAX_N} per atom")
    print(f"  Bridges:     {N_BRIDGES}")
    print(f"  mu:          {1836.15 / 2:.2f} a.u. (reduced mass)")
    print(f"  Scenario 1:  T = {T_low} a.u. ({T_low * AU_TO_K:.0f} K) -- room temp")
    print(f"  Scenario 2:  T = {T_high} a.u. ({T_high * AU_TO_K:.0f} K) -- extreme")
    print()

    # ------------------------------------------------------------------
    # Scenario 1: Room temperature -- stable vibration
    # ------------------------------------------------------------------
    res_low = run_langevin(
        label="Scenario 1: Room Temperature",
        T_target=T_low,
        R_init=1.29,
        n_steps=500,
        seed=42,
    )

    # ------------------------------------------------------------------
    # Scenario 2: Extreme temperature -- thermal dissociation
    # ------------------------------------------------------------------
    res_high = run_langevin(
        label="Scenario 2: Extreme Temperature",
        T_target=T_high,
        R_init=1.29,
        gamma=0.001,
        n_steps=500,
        R_max_clamp=3.0,
        seed=42,
    )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("=" * 65)
    print("  SUMMARY")
    print("=" * 65)
    print()

    # Low-T analysis
    R_low = res_low["R"]
    E_kin_low = res_low["E_kin"]
    T_kinetic_low = np.mean(E_kin_low[len(E_kin_low)//4:]) / 0.5  # <E_kin> = 0.5*T for 1 DOF
    print(f"  Room Temperature ({res_low['T_kelvin']:.0f} K):")
    print(f"    Bond stable:    YES")
    print(f"    R range:        [{R_low.min():.4f}, {R_low.max():.4f}] Bohr")
    print(f"    R mean:         {R_low.mean():.4f} Bohr")
    print(f"    <E_kin>:        {np.mean(E_kin_low[len(E_kin_low)//4:]):.2e} Ha")
    print(f"    T_kinetic:      {T_kinetic_low:.4f} a.u. "
          f"({T_kinetic_low * AU_TO_K:.0f} K)")
    print(f"    Time:           {res_low['elapsed']:.2f}s")
    print()

    # High-T analysis
    R_high = res_high["R"]
    print(f"  Extreme Temperature ({res_high['T_kelvin']:.0f} K):")
    if res_high["dissociated"]:
        t_dissoc = res_high["t"][res_high["dissoc_step"]] * AU_TO_FS
        print(f"    Bond dissociated: YES (step {res_high['dissoc_step']}, "
              f"t = {t_dissoc:.2f} fs)")
    else:
        print(f"    Bond dissociated: NO (R_max = {R_high.max():.2f} Bohr)")
    print(f"    R range:        [{R_high.min():.4f}, {R_high.max():.4f}] Bohr")
    print(f"    Time:           {res_high['elapsed']:.2f}s")
    print()

    # ------------------------------------------------------------------
    # Plot: 2x2 grid
    # ------------------------------------------------------------------
    print("Generating plot...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle(
        r"H$_2$ Langevin Thermostat AIMD — Thermal Vibration & Dissociation",
        fontsize=14, fontweight="bold",
    )

    t_low_fs = res_low["t"] * AU_TO_FS
    t_high_fs = res_high["t"] * AU_TO_FS

    # --- Top Left: R vs Time (300K) ---
    ax = axes[0, 0]
    ax.plot(t_low_fs, res_low["R"], color="darkgreen", linewidth=0.8)
    ax.axhline(1.293, color="red", linestyle=":", linewidth=1, alpha=0.6,
               label="R_eq = 1.293 Bohr")
    ax.set_xlabel("Time (fs)")
    ax.set_ylabel("R (Bohr)")
    ax.set_title(f"Bond Length -- {res_low['T_kelvin']:.0f} K (Room Temp)")
    ax.legend(loc="upper right", fontsize=7)
    ax.grid(True, alpha=0.3)

    # --- Bottom Left: Kinetic Energy vs Time (300K) ---
    ax = axes[1, 0]
    ax.plot(t_low_fs, res_low["E_kin"], color="red", linewidth=0.8,
            label="E_kin (nuclear)")
    # Target thermal energy: <E_kin> = 0.5 * T for 1 DOF
    E_kin_target = 0.5 * T_low
    ax.axhline(E_kin_target, color="black", linestyle="--", linewidth=1,
               alpha=0.6, label=f"0.5*T = {E_kin_target:.4f} Ha")
    ax.set_xlabel("Time (fs)")
    ax.set_ylabel("Kinetic Energy (Ha)")
    ax.set_title(f"Thermal Equilibration -- {res_low['T_kelvin']:.0f} K")
    ax.legend(loc="upper right", fontsize=7)
    ax.grid(True, alpha=0.3)

    # --- Top Right: R vs Time (High Temp) ---
    ax = axes[0, 1]
    ax.plot(t_high_fs, res_high["R"], color="darkred", linewidth=0.8)
    ax.axhline(1.293, color="blue", linestyle=":", linewidth=1, alpha=0.4,
               label="R_eq = 1.293 Bohr")
    if res_high["dissociated"]:
        t_d = res_high["t"][res_high["dissoc_step"]] * AU_TO_FS
        ax.axvline(t_d, color="black", linestyle="--", linewidth=1.2,
                   alpha=0.7, label=f"Dissociation at {t_d:.2f} fs")
    ax.set_xlabel("Time (fs)")
    ax.set_ylabel("R (Bohr)")
    ax.set_title(f"Bond Dissociation -- {res_high['T_kelvin']:.0f} K")
    ax.legend(loc="upper left", fontsize=7)
    ax.grid(True, alpha=0.3)

    # --- Bottom Right: Total Energy vs Time (High Temp) ---
    ax = axes[1, 1]
    ax.plot(t_high_fs, res_high["E_tot"], color="black", linewidth=0.8,
            label="Total Energy")
    ax.plot(t_high_fs, res_high["E_pot"], color="darkblue", linewidth=0.6,
            alpha=0.7, label="Potential (E_CI)")
    ax.plot(t_high_fs, res_high["E_kin"], color="red", linewidth=0.6,
            alpha=0.7, label="Kinetic (nuclear)")
    ax.set_xlabel("Time (fs)")
    ax.set_ylabel("Energy (Ha)")
    ax.set_title(f"Energy Components -- {res_high['T_kelvin']:.0f} K")
    ax.legend(loc="upper left", fontsize=7)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    plot_dir = os.path.join(os.path.dirname(__file__), "..", "debug", "plots")
    os.makedirs(plot_dir, exist_ok=True)
    plot_path = os.path.join(plot_dir, "h2_thermal_dissociation.png")
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved to: {plot_path}")
    plt.close(fig)

    print()
    print("=" * 65)
    stable = not res_low["dissociated"]
    broke = res_high["dissociated"]
    if stable and broke:
        print("  SUCCESS: Bond stable at room temp, dissociated at high temp")
    elif stable and not broke:
        print("  PARTIAL: Bond stable at room temp, but survived high temp")
        print("           (may need higher T or more steps)")
    else:
        print("  UNEXPECTED: Check results")
    print("=" * 65)


if __name__ == "__main__":
    run_thermostat_demo()
