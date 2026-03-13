"""
Coupled Electron-Nuclear Lattice
=================================

Couples the electronic LCAO-FCI graph (lattice_index.py) to the nuclear
Morse vibrational graph (nuclear_lattice.py) to compute:

    E_total(v) = E_elec(R(v)) + E_nuc(v)

where R(v) = <R>_v is the Morse expectation value at vibrational state v.

Key result: the nuclear graph's vibrational frequency omega_e encodes the
PES curvature at equilibrium. Requiring that the Fock-weighted kinetic
correction produces this curvature DERIVES lambda from first principles:

    lambda = (k_Morse - k_bare) / k_correction_per_unit_lambda

where k_Morse = 2 * a^2 * D_e is the Morse force constant.

Author: GeoVac Development Team
Date: March 2026
"""

import numpy as np
from scipy.special import digamma
from typing import Tuple, Optional, Dict, List

from geovac.nuclear_lattice import (
    MorseVibrationLattice, NuclearLattice, build_diatomic,
    DIATOMIC_CONSTANTS, HARTREE_TO_CM, AMU_TO_ME,
)


# ======================================================================
# Morse Oscillator Analytical Tools
# ======================================================================

def morse_range_parameter(omega_e_hartree: float, mu_me: float,
                          D_e_hartree: float) -> float:
    """
    Morse range parameter a = omega_e * sqrt(mu / (2 * D_e)).

    Parameters
    ----------
    omega_e_hartree : float
        Harmonic frequency in Hartree.
    mu_me : float
        Reduced mass in electron masses.
    D_e_hartree : float
        Dissociation energy in Hartree.

    Returns
    -------
    float
        Range parameter a in bohr^{-1}.
    """
    return omega_e_hartree * np.sqrt(mu_me / (2.0 * D_e_hartree))


def morse_lambda_parameter(mu_me: float, D_e_hartree: float,
                           a: float) -> float:
    """
    Morse dimensionless parameter lambda = sqrt(2 * mu * D_e) / a.

    This equals j + 1/2 where j is the SU(2) representation label.

    Parameters
    ----------
    mu_me : float
        Reduced mass in electron masses.
    D_e_hartree : float
        Dissociation energy in Hartree.
    a : float
        Range parameter in bohr^{-1}.

    Returns
    -------
    float
        Dimensionless Morse parameter lambda.
    """
    return np.sqrt(2.0 * mu_me * D_e_hartree) / a


def morse_force_constant(a: float, D_e_hartree: float) -> float:
    """
    Morse force constant k = 2 * a^2 * D_e (second derivative at r_e).

    This is the curvature of the Morse PES at equilibrium, in Hartree/bohr^2.

    Parameters
    ----------
    a : float
        Range parameter in bohr^{-1}.
    D_e_hartree : float
        Dissociation energy in Hartree.

    Returns
    -------
    float
        Force constant in Hartree / bohr^2.
    """
    return 2.0 * a**2 * D_e_hartree


def morse_expectation_R(v: int, r_e: float, a: float,
                        lam: float) -> float:
    """
    Exact expectation value <R>_v for Morse oscillator state v.

    Uses the analytical result (Dahl & Springborg 1988):

        <r>_v = r_e - (1/a) * [psi(2*lambda - 2*v) - ln(2*lambda)]

    where psi is the digamma function.

    The anharmonicity of the Morse potential shifts the average position
    to larger R for higher v (bond stretches with excitation).

    Parameters
    ----------
    v : int
        Vibrational quantum number.
    r_e : float
        Equilibrium bond length in bohr.
    a : float
        Morse range parameter in bohr^{-1}.
    lam : float
        Morse dimensionless parameter lambda.

    Returns
    -------
    float
        Expectation value <R>_v in bohr.
    """
    s = 2.0 * lam - 2.0 * v
    if s <= 0:
        return np.inf  # Beyond dissociation
    return r_e - (1.0 / a) * (digamma(s) - np.log(2.0 * lam))


def morse_classical_turning_points(v: int, r_e: float, a: float,
                                   D_e: float,
                                   E_v: float) -> Tuple[float, float]:
    """
    Classical inner and outer turning points for Morse state v.

    At energy E_v, the Morse potential V(r) = D_e[1-exp(-a(r-r_e))]^2
    equals E_v at two points:

        r_inner = r_e - (1/a) * ln(1 + sqrt(E_v/D_e))
        r_outer = r_e - (1/a) * ln(1 - sqrt(E_v/D_e))

    Parameters
    ----------
    v : int
        Vibrational quantum number (for documentation).
    r_e : float
        Equilibrium bond length in bohr.
    a : float
        Morse range parameter in bohr^{-1}.
    D_e : float
        Dissociation energy in Hartree.
    E_v : float
        Morse energy of state v in Hartree.

    Returns
    -------
    Tuple[float, float]
        (r_inner, r_outer) in bohr.
    """
    if E_v >= D_e or E_v <= 0:
        return (0.0, np.inf)
    ratio = np.sqrt(E_v / D_e)
    r_inner = r_e - (1.0 / a) * np.log(1.0 + ratio)
    r_outer = r_e - (1.0 / a) * np.log(1.0 - ratio)
    return (r_inner, r_outer)


def compute_morse_parameters(molecule: str) -> Dict[str, float]:
    """
    Compute all Morse analytical parameters for a diatomic molecule.

    Parameters
    ----------
    molecule : str
        Molecule name (e.g., 'LiH', 'H2').

    Returns
    -------
    dict
        Dictionary with keys: r_e, D_e, omega_e_hartree, omega_e_xe_hartree,
        mu_me, a, lam, j, v_max, k_morse.
    """
    if molecule not in DIATOMIC_CONSTANTS:
        raise ValueError(f"Unknown molecule '{molecule}'. "
                         f"Available: {list(DIATOMIC_CONSTANTS.keys())}")

    c = DIATOMIC_CONSTANTS[molecule]
    D_e = c['D_e']
    omega_e_h = c['omega_e'] / HARTREE_TO_CM
    omega_e_xe_h = c['omega_e_xe'] / HARTREE_TO_CM
    mu_me = c['mu_amu'] * AMU_TO_ME
    r_e = c['r_e']

    a = morse_range_parameter(omega_e_h, mu_me, D_e)
    lam = morse_lambda_parameter(mu_me, D_e, a)
    j = omega_e_h / (2.0 * omega_e_xe_h) - 0.5
    v_max = int(np.floor(j))
    k = morse_force_constant(a, D_e)

    return {
        'r_e': r_e,
        'D_e': D_e,
        'omega_e_hartree': omega_e_h,
        'omega_e_xe_hartree': omega_e_xe_h,
        'mu_me': mu_me,
        'a': a,
        'lam': lam,
        'j': j,
        'v_max': v_max,
        'k_morse': k,
    }


def compute_R_expectation_table(molecule: str,
                                v_max_override: Optional[int] = None
                                ) -> np.ndarray:
    """
    Compute <R>_v for all vibrational states of a molecule.

    Parameters
    ----------
    molecule : str
        Molecule name.
    v_max_override : int, optional
        Override maximum v (default: use Morse v_max).

    Returns
    -------
    np.ndarray
        Array of shape (n_states, 3): columns are [v, <R>_v, E_nuc(v)].
    """
    params = compute_morse_parameters(molecule)
    r_e = params['r_e']
    a = params['a']
    lam = params['lam']
    omega_h = params['omega_e_hartree']
    omegax_h = params['omega_e_xe_hartree']

    v_max = v_max_override if v_max_override is not None else params['v_max']

    table = []
    for v in range(v_max + 1):
        R_v = morse_expectation_R(v, r_e, a, lam)
        E_v = omega_h * (v + 0.5) - omegax_h * (v + 0.5)**2
        table.append([v, R_v, E_v])

    return np.array(table)


# ======================================================================
# Force Constant Matching: Deriving lambda
# ======================================================================

def numerical_force_constant(R_values: np.ndarray,
                             E_values: np.ndarray,
                             R_eq: float) -> float:
    """
    Estimate force constant k = d^2E/dR^2 at R_eq from discrete PES data.

    Uses polynomial fitting near the equilibrium point.

    Parameters
    ----------
    R_values : np.ndarray
        Bond lengths in bohr.
    E_values : np.ndarray
        Energies in Hartree at each R.
    R_eq : float
        Equilibrium bond length for evaluating the second derivative.

    Returns
    -------
    float
        Force constant in Hartree/bohr^2. Positive means a minimum exists.
    """
    # Fit polynomial of degree 4 near equilibrium
    # Use points within ±2 bohr of R_eq
    mask = np.abs(R_values - R_eq) < 2.0
    if mask.sum() < 4:
        mask = np.ones(len(R_values), dtype=bool)

    R_fit = R_values[mask]
    E_fit = E_values[mask]

    # Fit quartic polynomial
    coeffs = np.polyfit(R_fit - R_eq, E_fit, min(4, len(R_fit) - 1))
    poly = np.poly1d(coeffs)

    # Second derivative of polynomial at R_eq (= 0 in shifted coords)
    poly_dd = poly.deriv(2)
    return float(poly_dd(0.0))


def derive_lambda_from_force_constant(
    k_morse: float,
    R_values: np.ndarray,
    E_bare: np.ndarray,
    E_correction_per_lambda: np.ndarray,
    R_eq: float,
) -> Dict[str, float]:
    """
    Derive lambda by matching PES curvature to Morse force constant.

    The corrected PES is:
        E_total(R; lambda) = E_bare(R) + lambda * Delta_E(R)

    The force constant is:
        k(lambda) = k_bare + lambda * k_correction

    Setting k(lambda) = k_morse and solving:
        lambda = (k_morse - k_bare) / k_correction

    Parameters
    ----------
    k_morse : float
        Target force constant from nuclear spectroscopy (Hartree/bohr^2).
    R_values : np.ndarray
        Bond lengths used in PES computation.
    E_bare : np.ndarray
        Bare electronic energies (no kinetic correction) at each R.
    E_correction_per_lambda : np.ndarray
        Kinetic correction energy per unit lambda at each R.
    R_eq : float
        Equilibrium bond length (bohr).

    Returns
    -------
    dict
        Dictionary with keys: lambda_derived, k_morse, k_bare,
        k_correction, k_total.
    """
    k_bare = numerical_force_constant(R_values, E_bare, R_eq)
    k_corr = numerical_force_constant(R_values, E_correction_per_lambda, R_eq)

    if abs(k_corr) < 1e-12:
        lambda_derived = np.inf
    else:
        lambda_derived = (k_morse - k_bare) / k_corr

    k_total = k_bare + lambda_derived * k_corr if np.isfinite(lambda_derived) else np.nan

    return {
        'lambda_derived': lambda_derived,
        'k_morse': k_morse,
        'k_bare': k_bare,
        'k_correction_per_lambda': k_corr,
        'k_total': k_total,
    }


# ======================================================================
# Coupled E_total(v) Computation
# ======================================================================

def coupled_energy_table(
    v_R_table: np.ndarray,
    E_elec_at_R: Dict[float, float],
    E_separated: float,
) -> np.ndarray:
    """
    Compute E_total(v) = E_elec(R(v)) + E_nuc(v) for each vibrational state.

    Parameters
    ----------
    v_R_table : np.ndarray
        Output from compute_R_expectation_table: columns [v, <R>_v, E_nuc(v)].
    E_elec_at_R : dict
        Mapping R -> E_elec (electronic energy including V_NN).
        R values should be close to the <R>_v values. Linear interpolation
        is used for R values not exactly in the dict.
    E_separated : float
        Energy of separated atoms (for binding energy reference).

    Returns
    -------
    np.ndarray
        Array of shape (n_states, 6): columns are
        [v, R_v, E_nuc, E_elec, E_total, D_e_effective].
    """
    # Build interpolation from E_elec_at_R
    R_pts = np.array(sorted(E_elec_at_R.keys()))
    E_pts = np.array([E_elec_at_R[r] for r in R_pts])

    result = []
    for row in v_R_table:
        v, R_v, E_nuc = row[0], row[1], row[2]

        if np.isinf(R_v) or R_v > R_pts[-1]:
            E_elec = E_separated
        elif R_v < R_pts[0]:
            E_elec = E_pts[0]  # Extrapolation (clamp)
        else:
            E_elec = np.interp(R_v, R_pts, E_pts)

        E_total = E_elec + E_nuc
        D_eff = E_separated - E_elec  # Binding at this geometry

        result.append([v, R_v, E_nuc, E_elec, E_total, D_eff])

    return np.array(result)


def find_equilibrium_v(E_total: np.ndarray) -> Tuple[int, float]:
    """
    Find the vibrational state v* that minimizes E_total.

    Parameters
    ----------
    E_total : np.ndarray
        Total energies for each v (column 4 of coupled_energy_table output).

    Returns
    -------
    Tuple[int, float]
        (v_eq, E_min) — equilibrium v and its energy.
    """
    idx = np.argmin(E_total)
    return int(idx), float(E_total[idx])
