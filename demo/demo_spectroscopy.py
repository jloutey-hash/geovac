#!/usr/bin/env python3
"""
Delta-Kick Spectroscopy: Real-Time Excited State Spectrometer
=============================================================

Extracts the full absorption spectrum of the graph Hamiltonian from a single
time propagation using the delta-kick method:

    1. Prepare ground state psi_0
    2. Apply weak momentum kick: psi(0) = psi_0 + i*k*(Z*psi_0)
    3. Propagate under static H -- the kicked state "rings" at all
       dipole-active transition frequencies
    4. Record induced dipole mu(t) = <psi(t)|Z|psi(t)>
    5. Fourier transform mu(t) -> absorption spectrum

The peak positions reveal transition energies; peak heights reveal
oscillator strengths. No driving field needed -- one simulation gives
the entire spectrum.

Validation: FFT peaks are compared against exact eigenvalue gaps from
full diagonalization of the Hamiltonian. Lyman series energies are
shown as reference markers for the analytical (continuous) hydrogen.

Usage:
    python demo/demo_spectroscopy.py

Output:
    debug/plots/hydrogen_spectrum.png
"""

import sys
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.signal import find_peaks

# Ensure package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from geovac import AtomicSolver
from geovac.dynamics import TimePropagator


def run_spectroscopy(
    max_n: int = 10,
    kick_strength: float = 0.001,
    dt: float = 0.05,
    n_steps: int = 20000,
) -> None:
    """Run delta-kick spectroscopy on atomic hydrogen.

    Parameters
    ----------
    max_n : int
        Principal quantum number cutoff (basis size).
    kick_strength : float
        Momentum kick amplitude k (weak-field regime).
    dt : float
        Time step in atomic units.
    n_steps : int
        Number of propagation steps (T = dt * n_steps).
    """
    T_total = dt * n_steps
    freq_res = 2 * np.pi / T_total
    print("=" * 65)
    print("  Delta-Kick Spectroscopy -- Hydrogen Absorption Spectrum")
    print("=" * 65)
    print(f"  Basis:       max_n = {max_n}")
    print(f"  Kick:        k = {kick_strength}")
    print(f"  Time step:   dt = {dt} a.u.")
    print(f"  Total time:  T = {T_total:.1f} a.u.  ({n_steps} steps)")
    print(f"  Freq. resolution: dw = {freq_res:.5f} Ha")
    print()

    # ------------------------------------------------------------------
    # Step 1: Build system and solve for eigenstates
    # ------------------------------------------------------------------
    print("[1/5] Building Hamiltonian and solving eigenstates...")
    t0 = time.time()

    solver = AtomicSolver(max_n=max_n, Z=1)
    n_basis = solver.n_states

    # Full diagonalization for reference spectrum
    H_dense = solver.H.toarray()
    E_all, psi_all = eigh(H_dense)

    psi_gs = psi_all[:, 0]
    E_gs = E_all[0]

    print(f"       {n_basis} basis states, E_gs = {E_gs:.6f} Ha")
    print(f"       Eigenvalue range: [{E_all[0]:.6f}, {E_all[-1]:.6f}] Ha")
    print(f"       Done in {time.time() - t0:.2f}s")
    print()

    # ------------------------------------------------------------------
    # Step 2: Build dipole operator and catalog transitions
    # ------------------------------------------------------------------
    print("[2/5] Building dipole operator and applying delta-kick...")

    V_z = TimePropagator.build_dipole_z(solver.lattice)
    V_z_dense = V_z.toarray()

    # Catalog all dipole-active transitions from ground state
    dipole_transitions = []
    for i in range(1, len(E_all)):
        d_i = abs(float(psi_gs @ V_z_dense @ psi_all[:, i]))
        if d_i > 0.01:
            dE = E_all[i] - E_gs
            dipole_transitions.append({
                "state": i,
                "dE": dE,
                "mu": d_i,
                "strength": d_i**2 * dE,  # oscillator strength proxy
            })

    # Sort by transition energy
    dipole_transitions.sort(key=lambda x: x["dE"])

    # Merge degenerate transitions (same energy within tolerance)
    merged = []
    tol = 1e-4
    for tr in dipole_transitions:
        if merged and abs(tr["dE"] - merged[-1]["dE"]) < tol:
            # Sum oscillator strengths for degenerate levels
            merged[-1]["mu_total"] = np.sqrt(
                merged[-1]["mu_total"] ** 2 + tr["mu"] ** 2
            )
            merged[-1]["strength"] += tr["strength"]
            merged[-1]["count"] += 1
        else:
            merged.append({
                "dE": tr["dE"],
                "mu_total": tr["mu"],
                "strength": tr["strength"],
                "count": 1,
            })

    print(f"       {len(dipole_transitions)} dipole-active eigenstates")
    print(f"       {len(merged)} distinct transition energies")
    print()
    print("       Strongest transitions (exact eigenvalue gaps):")
    top_transitions = sorted(merged, key=lambda x: -x["strength"])[:10]
    for tr in sorted(top_transitions, key=lambda x: x["dE"]):
        resolvable = "Y" if tr["dE"] > 2 * freq_res else "~"
        print(
            f"         dE = {tr['dE']:.6f} Ha, "
            f"|mu| = {tr['mu_total']:.4f}, "
            f"f = {tr['strength']:.4f}, "
            f"degen = {tr['count']}, "
            f"resolvable: {resolvable}"
        )
    print()

    # Apply delta-kick: psi(0) ~ psi_0 + i*k*(Z*psi_0)
    z_psi = V_z.dot(psi_gs)
    psi0 = psi_gs.astype(complex) + 1j * kick_strength * z_psi
    psi0 /= np.linalg.norm(psi0)

    # ------------------------------------------------------------------
    # Step 3: Time propagation -- record induced dipole mu(t)
    # ------------------------------------------------------------------
    print("[3/5] Propagating (static Hamiltonian, no driving field)...")
    t0 = time.time()

    prop = TimePropagator(solver.H, dt=dt)

    times = np.zeros(n_steps + 1)
    dipole_signal = np.zeros(n_steps + 1)

    # Record mu(t=0)
    dipole_signal[0] = np.real(np.vdot(psi0, V_z.dot(psi0)))

    # Evolve step-by-step with the static Hamiltonian
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
    # Step 4: Fourier transform -> absorption spectrum
    # ------------------------------------------------------------------
    print("[4/5] Computing absorption spectrum...")

    # Apply Hann window to reduce spectral leakage
    window = np.hanning(len(dipole_osc))
    signal_windowed = dipole_osc * window

    # Zero-pad to 4x length for sharper peak resolution
    n_padded = 4 * len(signal_windowed)
    spectrum = np.fft.rfft(signal_windowed, n=n_padded)
    freqs = np.fft.rfftfreq(n_padded, d=dt)
    energies_fft = 2 * np.pi * freqs  # omega = 2*pi*f

    # Use amplitude spectrum (not power) to compress dynamic range
    amplitude = np.abs(spectrum)

    # Normalize
    amp_max = np.max(amplitude)
    amp_norm = amplitude / amp_max if amp_max > 0 else amplitude

    # Find peaks on log-scale amplitude to detect weak transitions
    # Add small floor to avoid log(0)
    log_amp = np.log10(amp_norm + 1e-12)
    log_amp_shifted = log_amp - np.min(log_amp)  # shift to positive
    log_max = np.max(log_amp_shifted)
    log_norm = log_amp_shifted / log_max if log_max > 0 else log_amp_shifted

    # Adaptive minimum distance between peaks
    dE_bin = energies_fft[1] - energies_fft[0] if len(energies_fft) > 1 else 1.0
    min_distance = max(3, int(0.003 / dE_bin))

    peak_indices, peak_props = find_peaks(
        log_norm,
        height=0.15,
        distance=min_distance,
        prominence=0.05,
    )

    # Sort by amplitude (strongest first)
    peak_amps = amp_norm[peak_indices]
    sorted_order = np.argsort(peak_amps)[::-1]
    peak_indices = peak_indices[sorted_order]

    peak_energies = energies_fft[peak_indices]

    print(f"       Found {len(peak_indices)} spectral peaks")
    print()

    # ------------------------------------------------------------------
    # Step 5: Validate FFT peaks against exact eigenvalue gaps
    # ------------------------------------------------------------------
    print("[5/5] Validating FFT peaks against exact Hamiltonian eigenvalues:")
    print()

    # Build list of exact transition energies (only resolvable ones)
    exact_transitions = [
        tr for tr in merged if tr["dE"] > 1.5 * freq_res
    ]

    print(f"  {'Exact dE (Ha)':>14} {'FFT Peak (Ha)':>14} {'Error':>10} "
          f"{'|mu|':>8} {'Strength':>10}")
    print(f"  {'-'*14} {'-'*14} {'-'*10} {'-'*8} {'-'*10}")

    matched = []
    used_peaks = set()  # Avoid double-matching
    for tr in exact_transitions[:15]:
        dE_exact = tr["dE"]
        best_match = None
        if len(peak_energies) > 0:
            dists = np.abs(peak_energies - dE_exact)
            # Find closest *unused* peak
            for attempt in np.argsort(dists):
                if attempt not in used_peaks and dists[attempt] < 2 * freq_res:
                    best_match = attempt
                    break
        if best_match is not None:
            dE_fft = peak_energies[best_match]
            err_pct = abs(dE_fft - dE_exact) / dE_exact * 100
            used_peaks.add(best_match)
            print(
                f"  {dE_exact:>14.6f} {dE_fft:>14.6f} "
                f"{err_pct:>9.2f}% {tr['mu_total']:>8.4f} "
                f"{tr['strength']:>10.4f}"
            )
            matched.append({
                "exact": dE_exact,
                "fft": dE_fft,
                "error_pct": err_pct,
            })
        else:
            print(
                f"  {dE_exact:>14.6f} {'(not found)':>14} "
                f"{'':>10} {tr['mu_total']:>8.4f} "
                f"{tr['strength']:>10.4f}"
            )

    print()
    if matched:
        errors = [m["error_pct"] for m in matched]
        print(f"  Peaks matched:  {len(matched)}/{len(exact_transitions[:15])}")
        print(f"  Mean error:     {np.mean(errors):.2f}%")
        print(f"  Max error:      {np.max(errors):.2f}%")
    else:
        print("  No peaks matched. Try increasing T_total for better resolution.")
    print()

    # Also show Lyman series comparison
    lyman_n = np.arange(2, 7)
    lyman_dE = 0.5 * (1.0 - 1.0 / lyman_n**2)
    print("  Reference: Analytical Lyman series (continuous hydrogen):")
    for n, dE in zip(lyman_n, lyman_dE):
        print(f"    1s -> {n}p:  dE = {dE:.6f} Ha")
    print("  (Graph Hamiltonian eigenvalues differ from the continuous limit)")
    print()

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    print("Generating plot...")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 9), height_ratios=[1, 1.5])
    fig.suptitle(
        "Delta-Kick Spectroscopy \u2014 Graph Hydrogen Absorption Spectrum",
        fontsize=14,
        fontweight="bold",
    )

    # --- Top panel: Time-domain dipole ringing ---
    ax1.plot(times, dipole_signal, color="steelblue", linewidth=0.3, alpha=0.8)
    ax1.set_xlabel("Time (a.u.)")
    ax1.set_ylabel(r"Dipole moment $\mu(t)$ (a.u.)")
    ax1.set_title("Induced Dipole Ringing After Delta-Kick")
    ax1.set_xlim(0, T_total)
    ax1.grid(True, alpha=0.3)

    # --- Bottom panel: Absorption spectrum (log scale) ---
    # Plot range: up to the max transition energy + margin
    E_max_plot = min(0.35, max(tr["dE"] for tr in merged) * 1.2)
    mask = (energies_fft > 0) & (energies_fft <= E_max_plot)

    ax2.semilogy(
        energies_fft[mask],
        amp_norm[mask],
        color="darkblue",
        linewidth=0.8,
        label="FFT Amplitude",
    )

    # Mark exact eigenvalue gaps as vertical lines
    for i, tr in enumerate(exact_transitions[:12]):
        lbl = f"Eigengap {tr['dE']:.4f} Ha" if i < 6 else None
        ax2.axvline(
            tr["dE"],
            color="red",
            linestyle="--",
            linewidth=1.0,
            alpha=0.5,
            label=lbl,
        )

    # Mark Lyman series as faint reference
    for n, dE in zip(lyman_n, lyman_dE):
        if dE <= E_max_plot:
            ax2.axvline(
                dE,
                color="gray",
                linestyle=":",
                linewidth=1.0,
                alpha=0.4,
                label=f"Lyman 1s->{n}p" if n == 2 else None,
            )

    # Mark matched FFT peaks
    for m in matched:
        idx = np.argmin(np.abs(energies_fft - m["fft"]))
        ax2.plot(
            m["fft"],
            amp_norm[idx],
            "v",
            color="limegreen",
            markersize=8,
            zorder=5,
            label="FFT peak" if m is matched[0] else None,
        )

    ax2.set_xlabel("Energy (Hartree)")
    ax2.set_ylabel("Spectral Amplitude (log scale)")
    ax2.set_title("Absorption Spectrum vs Exact Eigenvalue Gaps")
    ax2.set_xlim(0, E_max_plot)
    ax2.set_ylim(bottom=1e-5, top=2.0)
    ax2.legend(loc="upper right", fontsize=7, ncol=2)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save
    plot_dir = os.path.join(os.path.dirname(__file__), "..", "debug", "plots")
    os.makedirs(plot_dir, exist_ok=True)
    plot_path = os.path.join(plot_dir, "hydrogen_spectrum.png")
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved to: {plot_path}")
    plt.close(fig)

    print()
    print("=" * 65)
    if matched:
        print(f"  SUCCESS: {len(matched)} transition energies recovered from FFT")
        print(f"  Mean accuracy: {np.mean([m['error_pct'] for m in matched]):.2f}%")
    print("  Spectroscopy complete.")
    print("=" * 65)


if __name__ == "__main__":
    run_spectroscopy()
