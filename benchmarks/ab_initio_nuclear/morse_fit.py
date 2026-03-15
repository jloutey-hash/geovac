"""
Morse potential fit to ab initio PES data.

Extracts spectroscopic constants: D_e, R_e, a, omega_e, omega_e_xe, B_e, alpha_e.
"""

import numpy as np
from scipy.optimize import curve_fit
from typing import Dict, Tuple

# Physical constants
HARTREE_TO_CM = 219474.63       # cm^-1 per Hartree
BOHR_TO_M = 5.29177210903e-11   # meters per bohr
M_H_AMU = 1.00782503207        # hydrogen atomic mass in amu
AMU_TO_KG = 1.66053906660e-27  # kg per amu
HBAR_SI = 1.054571817e-34      # J*s
C_CM = 2.99792458e10            # cm/s
AMU_TO_ME = 1822.888486209     # electron masses per amu


def morse_potential(R: np.ndarray, D_e: float, a: float, R_e: float) -> np.ndarray:
    """Morse potential V(R) = D_e * [1 - exp(-a*(R - R_e))]^2 - D_e.

    The zero of energy is at the dissociation limit (R -> infinity).
    So at R = R_e, V = -D_e (the well bottom).

    Parameters
    ----------
    R : ndarray
        Bond lengths (bohr).
    D_e : float
        Well depth (Hartree, positive).
    a : float
        Range parameter (1/bohr).
    R_e : float
        Equilibrium bond length (bohr).

    Returns
    -------
    V : ndarray
        Potential energy relative to dissociation limit.
    """
    return D_e * (1.0 - np.exp(-a * (R - R_e)))**2 - D_e


def fit_morse(
    R: np.ndarray,
    E_total: np.ndarray,
    verbose: bool = True,
) -> Dict[str, float]:
    """Fit Morse potential to PES data.

    The PES data E_total(R) includes nuclear repulsion. We fit to the
    Morse form shifted so that V(infinity) = E_atoms = -1.0 Ha for H2.

    Parameters
    ----------
    R : ndarray
        Bond lengths (bohr).
    E_total : ndarray
        Total energies (Hartree).
    verbose : bool
        Print results.

    Returns
    -------
    dict with D_e, a, R_e, and spectroscopic constants.
    """
    # Remove any NaN values
    valid = ~np.isnan(E_total)
    R = R[valid]
    E_total = E_total[valid]

    # E_atoms for H2 = 2 * (-0.5) = -1.0 Ha
    E_atoms = -1.0

    # Shift so dissociation limit is at 0
    # V(R) = E_total(R) - E_atoms
    V = E_total - E_atoms

    # Initial guess: find minimum
    idx_min = np.argmin(E_total)
    R_e_guess = R[idx_min]
    D_e_guess = -V[idx_min]  # V at minimum is -D_e
    a_guess = 1.0

    try:
        popt, pcov = curve_fit(
            morse_potential, R, V,
            p0=[D_e_guess, a_guess, R_e_guess],
            bounds=([0.01, 0.1, 0.5], [1.0, 5.0, 5.0]),
            maxfev=10000,
        )
        D_e, a, R_e = popt
        perr = np.sqrt(np.diag(pcov))
    except Exception as e:
        if verbose:
            print(f"Morse fit failed: {e}")
        return {}

    # Compute spectroscopic constants
    # Reduced mass of H2: mu = m_H / 2
    mu_amu = M_H_AMU / 2.0
    mu_au = mu_amu * AMU_TO_ME  # in electron masses (atomic units)

    # omega_e = a * sqrt(2 * D_e / mu)  (in atomic units, then convert)
    omega_e_au = a * np.sqrt(2.0 * D_e / mu_au)
    omega_e_cm = omega_e_au * HARTREE_TO_CM

    # omega_e * x_e = omega_e^2 / (4 * D_e)  (in Hartree, then convert)
    omega_e_xe_au = omega_e_au**2 / (4.0 * D_e)
    omega_e_xe_cm = omega_e_xe_au * HARTREE_TO_CM

    # B_e = hbar^2 / (2 * mu * R_e^2) in atomic units, then convert
    # In atomic units: hbar = 1, so B_e = 1 / (2 * mu * R_e^2)
    B_e_au = 1.0 / (2.0 * mu_au * R_e**2)
    B_e_cm = B_e_au * HARTREE_TO_CM

    # alpha_e (vibration-rotation coupling)
    # alpha_e = -6 * B_e^2 / omega_e * (1 + a * R_e * B_e / omega_e)  (rough)
    # More accurate: alpha_e = 6 * sqrt(omega_e_xe * B_e^3) / omega_e - 6 * B_e^2 / omega_e
    # Use the standard Morse expression:
    # alpha_e = -6 * B_e^2 / omega_e + 6 * B_e * sqrt(B_e * omega_e_xe) / omega_e
    alpha_e_au = -6.0 * B_e_au**2 / omega_e_au + 6.0 * B_e_au * np.sqrt(B_e_au * omega_e_xe_au) / omega_e_au
    alpha_e_cm = abs(alpha_e_au) * HARTREE_TO_CM

    result = {
        'D_e': D_e,
        'a': a,
        'R_e': R_e,
        'omega_e_cm': omega_e_cm,
        'omega_e_xe_cm': omega_e_xe_cm,
        'B_e_cm': B_e_cm,
        'alpha_e_cm': alpha_e_cm,
        'mu_amu': mu_amu,
        'fit_error': perr,
    }

    if verbose:
        print("\n" + "=" * 60)
        print("MORSE FIT RESULTS")
        print("=" * 60)
        print(f"  D_e   = {D_e:.6f} Ha  ({D_e * 27.2114:.4f} eV)")
        print(f"  R_e   = {R_e:.4f} bohr")
        print(f"  a     = {a:.4f} bohr^-1")
        print(f"  omega_e    = {omega_e_cm:.2f} cm^-1")
        print(f"  omega_e_xe = {omega_e_xe_cm:.2f} cm^-1")
        print(f"  B_e        = {B_e_cm:.3f} cm^-1")
        print(f"  alpha_e    = {alpha_e_cm:.3f} cm^-1")

    return result


# Experimental H2 values
EXPERIMENTAL_H2 = {
    'D_e': 0.1745,
    'R_e': 1.401,
    'omega_e_cm': 4401.21,
    'omega_e_xe_cm': 121.34,
    'B_e_cm': 60.853,
    'alpha_e_cm': 3.062,
}


def compare_to_experiment(
    ai_params: Dict[str, float],
    verbose: bool = True,
) -> Dict[str, float]:
    """Compare ab initio Morse parameters to experiment.

    Returns
    -------
    dict with percentage errors.
    """
    exp = EXPERIMENTAL_H2
    errors = {}

    keys = ['D_e', 'R_e', 'omega_e_cm', 'omega_e_xe_cm', 'B_e_cm', 'alpha_e_cm']
    labels = ['D_e (Ha)', 'R_e (bohr)', 'omega_e (cm-1)', 'omega_e*x_e (cm-1)',
              'B_e (cm-1)', 'alpha_e (cm-1)']

    if verbose:
        print("\n" + "=" * 70)
        print("COMPARISON TO EXPERIMENT")
        print("=" * 70)
        print(f"  {'Parameter':<20s}  {'Ab initio':>12s}  {'Expt':>12s}  {'Error %':>10s}")
        print("-" * 60)

    for key, label in zip(keys, labels):
        if key in ai_params and key in exp:
            ai = ai_params[key]
            ex = exp[key]
            err_pct = (ai - ex) / ex * 100.0
            errors[key] = err_pct
            if verbose:
                print(f"  {label:<20s}  {ai:12.4f}  {ex:12.4f}  {err_pct:+10.2f}%")

    return errors
