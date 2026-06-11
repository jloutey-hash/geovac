"""G4-6a analytical apex correction via Bessel-zero eigenvalues.

The substrate's radial eigenvalues (from tridiagonal FD) are wrong for high-m
azimuthal modes because the centrifugal barrier m^2/rho^2 at the first grid
point (rho=a) is unresolved. The CONTINUUM eigenvalues are (j_{j,m}/R)^2
where j_{j,m} are zeros of J_m.

Strategy:
1. For each azimuthal mode m = (k+1/2)/alpha:
   - Compute continuum eigenvalues from Bessel zeros
   - Compute FD eigenvalues from the tridiagonal matrix (same as substrate)
   - Per-mode correction DK_rad(t,m) = K_rad^cont - K_rad^FD
2. Sum over modes with alpha-derivative to get tip correction
3. tip_corrected = tip_substrate + correction
"""

import numpy as np
from scipy.special import jv
from scipy.optimize import brentq
import sys
import os
import json

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from geovac.gravity.warped_dirac import (
    DiscreteDiskDiracSpectral, DiscreteWedgeDiracSpectral
)


def bessel_zeros(nu: float, n_zeros: int) -> np.ndarray:
    """Compute first n_zeros positive zeros of J_nu(x) for real nu >= 0.

    Uses McMahon asymptotic for initial guesses, then brentq polish.
    McMahon: j_{nu,s} ~ beta - (mu-1)/(8*beta) - ...
    where beta = (s + nu/2 - 1/4)*pi, mu = 4*nu^2.
    """
    mu = 4 * nu**2
    zeros = []
    for s in range(1, n_zeros + 5):  # overshoot slightly
        beta = (s + nu / 2.0 - 0.25) * np.pi
        # McMahon first correction
        x_guess = beta - (mu - 1) / (8 * beta)
        if x_guess <= 0:
            x_guess = nu + s * np.pi  # fallback

        # Bracket: search around guess
        dx = np.pi * 0.45
        x_lo = max(1e-10, x_guess - dx)
        x_hi = x_guess + dx

        # Ensure bracket contains a sign change
        f_lo = jv(nu, x_lo)
        f_hi = jv(nu, x_hi)
        if f_lo * f_hi > 0:
            # Widen bracket
            for expand in range(1, 10):
                x_lo = max(1e-10, x_guess - dx * (1 + expand * 0.5))
                x_hi = x_guess + dx * (1 + expand * 0.5)
                f_lo = jv(nu, x_lo)
                f_hi = jv(nu, x_hi)
                if f_lo * f_hi < 0:
                    break
            if f_lo * f_hi > 0:
                continue

        try:
            z = brentq(lambda x: jv(nu, x), x_lo, x_hi, xtol=1e-12)
            if z > 0 and (len(zeros) == 0 or z > zeros[-1] + 0.5):
                zeros.append(z)
        except (ValueError, RuntimeError):
            continue

        if len(zeros) >= n_zeros:
            break

    return np.array(zeros[:n_zeros])


def continuum_radial_eigenvalues(m: float, R: float, n_zeros: int = 100) -> np.ndarray:
    """Compute continuum eigenvalues (j_{j,m}/R)^2 for the disk-Dirac radial problem."""
    if m < 0:
        m = abs(m)
    zeros = bessel_zeros(m, n_zeros)
    return (zeros / R) ** 2


def fd_radial_eigenvalues(m: float, N_rho: int, a: float) -> np.ndarray:
    """Compute FD radial eigenvalues from the tridiagonal Hermitian polar Laplacian."""
    H = np.zeros((N_rho, N_rho))
    for i in range(N_rho):
        k = i + 1
        rho = k * a
        H[i, i] = 2.0 / a**2 + (m**2 - 0.25) / rho**2
        if i > 0:
            H[i, i - 1] = -1.0 / a**2
        if i < N_rho - 1:
            H[i, i + 1] = -1.0 / a**2
    return np.linalg.eigvalsh(H)


def per_mode_heat_trace_correction(m: float, t: float, R: float,
                                    N_rho: int, a: float,
                                    n_bessel_zeros: int = 200) -> float:
    """Compute DK_rad(t, m) = K_rad^continuum - K_rad^FD for one azimuthal mode."""
    # Continuum: sum over Bessel zeros (truncated when contribution < eps)
    cont_eigs = continuum_radial_eigenvalues(m, R, n_bessel_zeros)
    K_cont = np.sum(np.exp(-t * cont_eigs))

    # FD: sum over tridiagonal eigenvalues
    fd_eigs = fd_radial_eigenvalues(m, N_rho, a)
    K_fd = np.sum(np.exp(-t * fd_eigs))

    return K_cont - K_fd


def compute_tip_with_correction(a, N_rho, N_phi, t_values, alpha_eps=0.1,
                                 n_bessel=100):
    """Compute tip_corrected = tip_substrate + analytical apex correction."""
    R = N_rho * a
    alpha_plus = 1.0 + alpha_eps
    alpha_minus = 1.0 - alpha_eps

    # Precompute substrate heat traces for all t
    print("    Building substrate objects...", flush=True)
    disk = DiscreteDiskDiracSpectral(N_rho, a, N_phi)
    wp = DiscreteWedgeDiracSpectral(N_rho, a, N_phi, alpha_plus)
    wm = DiscreteWedgeDiracSpectral(N_rho, a, N_phi, alpha_minus)

    tips_raw = []
    for t in t_values:
        K_disk = disk.heat_trace(t)
        K_wp = wp.heat_trace(t)
        K_wm = wm.heat_trace(t)
        tip_raw = (K_wp - K_wm) / (2 * alpha_eps) - K_disk
        tips_raw.append(tip_raw)

    # Precompute per-mode corrections (expensive: Bessel zeros + FD eigenvalues)
    print(f"    Computing per-mode Bessel corrections for {N_phi} modes...", flush=True)
    corrections = np.zeros(len(t_values))

    modes_done = 0
    for k_idx in range(N_phi):
        if k_idx <= N_phi // 2:
            k = k_idx
        else:
            k = k_idx - N_phi

        m_disk = abs(k + 0.5)
        m_wp = abs(k + 0.5) / alpha_plus
        m_wm = abs(k + 0.5) / alpha_minus

        # Continuum eigenvalues (Bessel zeros)
        cont_eigs_disk = continuum_radial_eigenvalues(m_disk, R, n_bessel)
        cont_eigs_wp = continuum_radial_eigenvalues(m_wp, R, n_bessel)
        cont_eigs_wm = continuum_radial_eigenvalues(m_wm, R, n_bessel)

        # FD eigenvalues
        fd_eigs_disk = fd_radial_eigenvalues(m_disk, N_rho, a)
        fd_eigs_wp = fd_radial_eigenvalues(m_wp, N_rho, a)
        fd_eigs_wm = fd_radial_eigenvalues(m_wm, N_rho, a)

        for i, t in enumerate(t_values):
            # Per-mode correction = (continuum - FD), with spinor factor 2
            dK_disk = np.sum(np.exp(-t * cont_eigs_disk)) - np.sum(np.exp(-t * fd_eigs_disk))
            dK_wp = np.sum(np.exp(-t * cont_eigs_wp)) - np.sum(np.exp(-t * fd_eigs_wp))
            dK_wm = np.sum(np.exp(-t * cont_eigs_wm)) - np.sum(np.exp(-t * fd_eigs_wm))

            # tip correction for this mode (factor 2 for spinor doubling)
            mode_tip_corr = 2 * ((dK_wp - dK_wm) / (2 * alpha_eps) - dK_disk)
            corrections[i] += mode_tip_corr

        modes_done += 1
        if modes_done % 20 == 0:
            print(f"      {modes_done}/{N_phi} modes done", flush=True)

    tips_corrected = np.array(tips_raw) + corrections
    return np.array(tips_raw), corrections, tips_corrected


def main():
    print("=" * 72)
    print("G4-6a BESSEL CORRECTION: analytical apex UV from continuum eigenvalues")
    print("=" * 72)

    A_cont = 1.0 / (24 * np.pi)
    B_cont = 1.0 / 6.0

    a = 0.05
    N_rho = 200
    N_phi = 120  # Production substrate
    R = N_rho * a
    # Also add tail modes beyond N_phi that the substrate doesn't have
    N_tail = 80  # Extra modes beyond substrate (m from N_phi/2+0.5 to N_phi/2+N_tail+0.5)

    t_test = np.array([0.5, 1.0, 2.0, 5.0, 10.0])

    print(f"\n  Parameters: a={a}, N_rho={N_rho}, N_phi={N_phi}, R={R}, N_tail={N_tail}")
    print(f"  A_cont = {A_cont:.8f}, B_cont = {B_cont:.8f}")
    n_bessel = N_rho  # Match Bessel zeros to FD count for 1:1 replacement
    print(f"  Bessel zeros per mode: {n_bessel} (matched to N_rho)")

    print(f"\n  Computing substrate + per-mode correction...", flush=True)
    tips_raw, corrections, tips_corr = compute_tip_with_correction(
        a, N_rho, N_phi, t_test, alpha_eps=0.1, n_bessel=n_bessel
    )

    # TAIL CORRECTION: add continuum-only modes beyond the substrate
    # These modes have m = N_phi/2 + 0.5 + k for k = 0, 1, ..., N_tail-1
    # The substrate doesn't have them at all, so full continuum contribution is added
    print(f"    Computing tail correction ({N_tail} extra modes)...", flush=True)
    alpha_eps_val = 0.1
    alpha_plus = 1.0 + alpha_eps_val
    alpha_minus = 1.0 - alpha_eps_val
    tail_corrections = np.zeros(len(t_test))

    m_base = N_phi // 2 + 0.5  # First mode beyond substrate
    for k_tail in range(N_tail):
        m_tail = m_base + k_tail
        # Full continuum contribution (no FD to subtract — mode absent from substrate)
        cont_eigs_d = continuum_radial_eigenvalues(m_tail, R, N_rho)
        cont_eigs_p = continuum_radial_eigenvalues(m_tail / alpha_plus, R, N_rho)
        cont_eigs_m_val = continuum_radial_eigenvalues(m_tail / alpha_minus, R, N_rho)

        for i, t in enumerate(t_test):
            K_d = np.sum(np.exp(-t * cont_eigs_d))
            K_p = np.sum(np.exp(-t * cont_eigs_p))
            K_m_v = np.sum(np.exp(-t * cont_eigs_m_val))
            # Factor 2 for spinor, factor 2 for +/- k symmetry
            mode_tip = 4 * ((K_p - K_m_v) / (2 * alpha_eps_val) - K_d)
            tail_corrections[i] += mode_tip

        if (k_tail + 1) % 20 == 0:
            print(f"      {k_tail+1}/{N_tail} tail modes", flush=True)

    tips_corr_full = tips_corr + tail_corrections

    print(f"\n  {'t':>6} | {'tip_raw':>10} {'mode_corr':>10} {'tail_corr':>10} "
          f"{'tip_full':>10} | {'continuum':>10} {'raw%':>6} {'full%':>6}")
    print("  " + "-" * 95)

    for i, t in enumerate(t_test):
        cont = A_cont / t + B_cont
        raw_rec = tips_raw[i] / cont * 100
        full_rec = tips_corr_full[i] / cont * 100
        print(f"  {t:6.3f} | {tips_raw[i]:+10.5f} {corrections[i]:+10.5f} "
              f"{tail_corrections[i]:+10.5f} {tips_corr_full[i]:+10.5f} | "
              f"{cont:10.5f} {raw_rec:6.1f}% {full_rec:6.1f}%")

    # A extraction
    print(f"\n  A EXTRACTION: A_est = t * (tip - B_cont)")
    print(f"  {'t':>6} | {'A_raw':>10} {'A_full_corr':>12} | {'target':>10} {'recovery':>10}")
    print("  " + "-" * 65)
    for i, t in enumerate(t_test):
        A_raw = t * (tips_raw[i] - B_cont)
        A_full = t * (tips_corr_full[i] - B_cont)
        rec = A_full / A_cont * 100
        print(f"  {t:6.3f} | {A_raw:+10.6f} {A_full:+12.8f} | {A_cont:10.8f} {rec:+10.1f}%")

    # Save
    output = {
        'parameters': {'a': a, 'N_rho': N_rho, 'N_phi': N_phi, 'R': R},
        'A_continuum': A_cont, 'B_continuum': B_cont,
        't_values': t_test.tolist(),
        'tips_raw': tips_raw.tolist(),
        'corrections': corrections.tolist(),
        'tips_corrected': tips_corr.tolist(),
    }
    os.makedirs('debug/data', exist_ok=True)
    with open('debug/data/g4_6a_bessel_correction.json', 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\n  Saved to debug/data/g4_6a_bessel_correction.json")


if __name__ == '__main__':
    main()
