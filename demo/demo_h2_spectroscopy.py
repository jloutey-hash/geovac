#!/usr/bin/env python3
"""
H2 Molecular Delta-Kick Spectroscopy
=====================================

Applies the delta-kick method to the H2 molecule to extract
the molecular absorption spectrum from a single time propagation.

The molecular dipole operator is constructed as a block-diagonal
matrix from the atomic dipole operators of each hydrogen atom,
embedded in the full molecular basis.

This demonstrates that the TimePropagator handles multi-center
molecular topologies with the same exactness as single atoms:
perfect unitarity over thousands of steps.

Usage:
    python demo/demo_h2_spectroscopy.py

Output:
    debug/plots/h2_spectrum.png
"""

import sys
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.sparse import lil_matrix, block_diag
from scipy.signal import find_peaks

# Ensure package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from geovac import GeometricLattice, MoleculeHamiltonian
from geovac.dynamics import TimePropagator


def build_molecular_dipole_z(
    lattices: list,
    n_total: int,
) -> "csr_matrix":
    """Build molecular Z-dipole operator as block-diagonal from atomic dipoles.

    Each atomic dipole operator is placed at the correct block offset
    in the full molecular basis. Cross-atom dipole couplings are neglected
    (they vanish for well-separated atoms and are small for bonded ones).

    Parameters
    ----------
    lattices : list of GeometricLattice
        The atomic lattices composing the molecule.
    n_total : int
        Total number of states in the molecular system.

    Returns
    -------
    V_z : csr_matrix
        Molecular dipole operator (sparse, real, symmetric).
    """
    atomic_dipoles = [TimePropagator.build_dipole_z(lat) for lat in lattices]
    return block_diag(atomic_dipoles, format="csr")


def run_h2_spectroscopy(
    max_n: int = 10,
    n_bridges: int = 40,
    kick_strength: float = 0.001,
    dt: float = 0.05,
    n_steps: int = 20000,
) -> None:
    """Run delta-kick spectroscopy on the H2 molecule.

    Parameters
    ----------
    max_n : int
        Principal quantum number cutoff per atom.
    n_bridges : int
        Number of topological bridge edges connecting the two atoms.
    kick_strength : float
        Momentum kick amplitude k.
    dt : float
        Time step in atomic units.
    n_steps : int
        Number of propagation steps.
    """
    T_total = dt * n_steps
    freq_res = 2 * np.pi / T_total
    HA_TO_EV = 27.2114

    print("=" * 65)
    print("  H2 Molecular Delta-Kick Spectroscopy")
    print("=" * 65)
    print(f"  Basis:       max_n = {max_n} per atom")
    print(f"  Bridges:     {n_bridges} topological edges")
    print(f"  Kick:        k = {kick_strength}")
    print(f"  Time step:   dt = {dt} a.u.")
    print(f"  Total time:  T = {T_total:.1f} a.u.  ({n_steps} steps)")
    print(f"  Freq. resolution: dw = {freq_res:.5f} Ha")
    print()

    # ------------------------------------------------------------------
    # Step 1: Build H2 molecular Hamiltonian
    # ------------------------------------------------------------------
    print("[1/5] Building H2 molecular Hamiltonian...")
    t0 = time.time()

    kinetic_scale = -1 / 16
    atom_A = GeometricLattice(max_n=max_n)
    atom_B = GeometricLattice(max_n=max_n)

    h2 = MoleculeHamiltonian(
        lattices=[atom_A, atom_B],
        connectivity=[(0, 1, n_bridges)],
        kinetic_scale=kinetic_scale,
    )

    n_mol = h2.n_total_states
    H_mol = h2.hamiltonian

    print(f"       States per atom: {atom_A.num_states}")
    print(f"       Total molecular states: {n_mol}")
    print(f"       Bridges created: {h2.bridge_info[0]['n_bridges_actual']}")
    print(f"       Done in {time.time() - t0:.2f}s")
    print()

    # ------------------------------------------------------------------
    # Step 2: Full diagonalization for reference spectrum
    # ------------------------------------------------------------------
    print("[2/5] Solving eigenstates and building dipole operator...")
    t0 = time.time()

    H_dense = H_mol.toarray()
    E_all, psi_all = eigh(H_dense)

    psi_gs = psi_all[:, 0]
    E_gs = E_all[0]

    # Build molecular dipole operator
    V_z = build_molecular_dipole_z(h2.lattices, n_mol)
    V_z_dense = V_z.toarray()

    # Find dipole-active transitions
    dipole_transitions = []
    for i in range(1, len(E_all)):
        d_i = abs(float(psi_gs @ V_z_dense @ psi_all[:, i]))
        if d_i > 0.01:
            dE = E_all[i] - E_gs
            dipole_transitions.append({
                "state": i,
                "dE": dE,
                "mu": d_i,
                "strength": d_i**2 * dE,
            })

    dipole_transitions.sort(key=lambda x: x["dE"])

    print(f"       E_gs = {E_gs:.6f} Ha")
    print(f"       Eigenvalue range: [{E_all[0]:.6f}, {E_all[-1]:.6f}] Ha")
    print(f"       {len(dipole_transitions)} dipole-active transitions")
    print(f"       Done in {time.time() - t0:.2f}s")
    print()

    # Print top transitions by oscillator strength
    top_by_strength = sorted(dipole_transitions, key=lambda x: -x["strength"])[:10]
    print("       Top 10 transitions by oscillator strength:")
    for tr in sorted(top_by_strength, key=lambda x: x["dE"]):
        print(
            f"         dE = {tr['dE']:.6f} Ha "
            f"({tr['dE'] * HA_TO_EV:.3f} eV), "
            f"|mu| = {tr['mu']:.4f}, "
            f"f = {tr['strength']:.4f}"
        )
    print()

    # ------------------------------------------------------------------
    # Step 3: Apply delta-kick and propagate
    # ------------------------------------------------------------------
    print("[3/5] Applying delta-kick and propagating...")
    t0 = time.time()

    # Kicked state: psi(0) ~ psi_0 + i*k*(Z*psi_0)
    z_psi = V_z.dot(psi_gs)
    psi0 = psi_gs.astype(complex) + 1j * kick_strength * z_psi
    psi0 /= np.linalg.norm(psi0)

    prop = TimePropagator(H_mol, dt=dt)

    times = np.zeros(n_steps + 1)
    dipole_signal = np.zeros(n_steps + 1)

    # Record mu(t=0)
    dipole_signal[0] = np.real(np.vdot(psi0, V_z.dot(psi0)))

    # Evolve
    psi_t = psi0.copy()
    for step in range(1, n_steps + 1):
        psi_t = prop.step(psi_t)
        times[step] = step * dt
        dipole_signal[step] = np.real(np.vdot(psi_t, V_z.dot(psi_t)))

    # Subtract DC offset
    mu_dc = dipole_signal[0]
    dipole_osc = dipole_signal - mu_dc

    elapsed = time.time() - t0
    norm_final = np.abs(np.vdot(psi_t, psi_t))
    print(f"       Done in {elapsed:.2f}s  ({n_steps / elapsed:.0f} steps/s)")
    print(f"       Norm at end: {norm_final:.12f}")
    print(f"       Max |mu_osc|: {np.max(np.abs(dipole_osc)):.2e}")
    print()

    # ------------------------------------------------------------------
    # Step 4: FFT -> absorption spectrum
    # ------------------------------------------------------------------
    print("[4/5] Computing absorption spectrum...")

    # Hann window + zero-padding
    window = np.hanning(len(dipole_osc))
    signal_windowed = dipole_osc * window

    n_padded = 4 * len(signal_windowed)
    spectrum = np.fft.rfft(signal_windowed, n=n_padded)
    freqs = np.fft.rfftfreq(n_padded, d=dt)
    energies_fft = 2 * np.pi * freqs

    # Amplitude spectrum (compressed dynamic range)
    amplitude = np.abs(spectrum)
    amp_max = np.max(amplitude)
    amp_norm = amplitude / amp_max if amp_max > 0 else amplitude

    # Log-scale peak detection
    log_amp = np.log10(amp_norm + 1e-12)
    log_amp_shifted = log_amp - np.min(log_amp)
    log_max = np.max(log_amp_shifted)
    log_norm = log_amp_shifted / log_max if log_max > 0 else log_amp_shifted

    dE_bin = energies_fft[1] - energies_fft[0] if len(energies_fft) > 1 else 1.0
    min_distance = max(3, int(0.003 / dE_bin))

    peak_indices, _ = find_peaks(
        log_norm,
        height=0.15,
        distance=min_distance,
        prominence=0.05,
    )

    # Sort by amplitude
    peak_amps = amp_norm[peak_indices]
    sorted_order = np.argsort(peak_amps)[::-1]
    peak_indices = peak_indices[sorted_order]
    peak_energies = energies_fft[peak_indices]

    print(f"       Found {len(peak_indices)} spectral peaks")
    print()

    # ------------------------------------------------------------------
    # Step 5: Report top molecular absorption peaks
    # ------------------------------------------------------------------
    print("[5/5] Molecular absorption spectrum -- top peaks:")
    print()

    # Match FFT peaks to exact eigenvalue gaps
    exact_resolvable = [
        tr for tr in dipole_transitions if tr["dE"] > 1.5 * freq_res
    ]

    matched = []
    used_peaks = set()
    for tr in exact_resolvable:
        if len(peak_energies) == 0:
            break
        dists = np.abs(peak_energies - tr["dE"])
        for attempt in np.argsort(dists):
            if attempt not in used_peaks and dists[attempt] < 2 * freq_res:
                dE_fft = peak_energies[attempt]
                err_pct = abs(dE_fft - tr["dE"]) / tr["dE"] * 100
                used_peaks.add(attempt)
                matched.append({
                    "exact": tr["dE"],
                    "fft": dE_fft,
                    "error_pct": err_pct,
                    "mu": tr["mu"],
                    "strength": tr["strength"],
                })
                break

    # Sort matched by strength and show top 5
    matched_sorted = sorted(matched, key=lambda m: -m["strength"])
    n_show = min(5, len(matched_sorted))

    print(f"  {'#':<4} {'Exact (Ha)':>12} {'FFT (Ha)':>12} {'(eV)':>10} "
          f"{'Error':>8} {'|mu|':>8} {'Strength':>10}")
    print(f"  {'-'*4} {'-'*12} {'-'*12} {'-'*10} "
          f"{'-'*8} {'-'*8} {'-'*10}")

    for rank, m in enumerate(matched_sorted[:n_show], 1):
        print(
            f"  {rank:<4} {m['exact']:>12.6f} {m['fft']:>12.6f} "
            f"{m['fft'] * HA_TO_EV:>10.3f} "
            f"{m['error_pct']:>7.2f}% {m['mu']:>8.4f} "
            f"{m['strength']:>10.4f}"
        )

    print()
    n_resolvable = len(exact_resolvable)
    print(f"  Total resolvable transitions: {n_resolvable}")
    print(f"  Peaks matched: {len(matched)}/{n_resolvable}")
    if matched:
        errors = [m["error_pct"] for m in matched]
        print(f"  Mean error: {np.mean(errors):.2f}%")
        print(f"  Max error:  {np.max(errors):.2f}%")
    print()

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    print("Generating plot...")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 9), height_ratios=[1, 1.5])
    fig.suptitle(
        r"H$_2$ Molecular Delta-Kick Spectroscopy",
        fontsize=14,
        fontweight="bold",
    )

    # --- Top panel: Time-domain dipole ringing ---
    ax1.plot(times, dipole_signal, color="steelblue", linewidth=0.3, alpha=0.8)
    ax1.set_xlabel("Time (a.u.)")
    ax1.set_ylabel(r"Dipole moment $\mu(t)$ (a.u.)")
    ax1.set_title(
        r"Induced Dipole Ringing After Delta-Kick (H$_2$ molecule)"
    )
    ax1.set_xlim(0, T_total)
    ax1.grid(True, alpha=0.3)

    # --- Bottom panel: Absorption spectrum (log scale) ---
    # Set plot range to cover the actual transition energies
    if dipole_transitions:
        max_dE = max(tr["dE"] for tr in dipole_transitions)
        E_max_plot = min(max_dE * 1.15, 1.5)
        E_max_plot = max(E_max_plot, 0.15)
    else:
        E_max_plot = 0.35
    mask = (energies_fft > 0) & (energies_fft <= E_max_plot)

    ax2.semilogy(
        energies_fft[mask],
        amp_norm[mask],
        color="darkblue",
        linewidth=0.8,
        label="FFT Amplitude",
    )

    # Mark exact eigenvalue gaps
    for i, tr in enumerate(exact_resolvable[:15]):
        lbl = f"Eigengap {tr['dE']:.4f} Ha" if i < 5 else None
        if tr["dE"] <= E_max_plot:
            ax2.axvline(
                tr["dE"],
                color="red",
                linestyle="--",
                linewidth=0.8,
                alpha=0.4,
                label=lbl,
            )

    # Mark matched FFT peaks
    for i, m in enumerate(matched_sorted[:n_show]):
        idx = np.argmin(np.abs(energies_fft - m["fft"]))
        ax2.plot(
            m["fft"],
            amp_norm[idx],
            "v",
            color="limegreen",
            markersize=8,
            zorder=5,
            label=f"Peak {i+1}: {m['fft']:.4f} Ha ({m['fft']*HA_TO_EV:.2f} eV)"
            if i < 5 else None,
        )

    ax2.set_xlabel("Energy (Hartree)")
    ax2.set_ylabel("Spectral Amplitude (log scale)")
    ax2.set_title(r"H$_2$ Absorption Spectrum vs Exact Eigenvalue Gaps")
    ax2.set_xlim(0, E_max_plot)
    ax2.set_ylim(bottom=1e-5, top=2.0)
    ax2.legend(loc="upper right", fontsize=7, ncol=2)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save
    plot_dir = os.path.join(os.path.dirname(__file__), "..", "debug", "plots")
    os.makedirs(plot_dir, exist_ok=True)
    plot_path = os.path.join(plot_dir, "h2_spectrum.png")
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved to: {plot_path}")
    plt.close(fig)

    print()
    print("=" * 65)
    if matched:
        print(f"  SUCCESS: {len(matched)} molecular transitions recovered")
        print(f"  Mean accuracy: {np.mean([m['error_pct'] for m in matched]):.2f}%")
        print(f"  Norm conservation: {norm_final:.12f}")
    print("  H2 spectroscopy complete.")
    print("=" * 65)


if __name__ == "__main__":
    run_h2_spectroscopy()
