"""
Perturbative cusp energy correction for Level 3 (He) and Level 4 (H2).

The Kato cusp condition requires the exact wavefunction to satisfy:
    d(Psi_bar)/d(r12)|_{r12=0} = (1/2) Psi_bar(r12=0)

where Psi_bar is the spherically averaged wavefunction around the coalescence
point.  A finite partial-wave expansion (l_max < infinity) cannot exactly
satisfy this derivative condition because the Coulomb cusp requires all l.

The energy error from the missing cusp is estimated from the known
partial-wave convergence law (Schwartz 1962, Hill 1985, Kutzelnigg & Morgan
1992):

    Delta E_l = E_l - E_{l-1} ~ -A / (l + 1/2)^4

with leading coefficient (Hill 1985)

    A = 6 pi^2 * int_0^inf |Psi(r, r, r12 = 0)|^2 r^5 dr .

For the exact He ground state A = 0.074226 Ha (Bromley & Mitroy,
arXiv:physics/0601049, Eqs. 8-9, written there as A_E = -0.074226 with the
opposite sign convention).  The repo's own 2D variational increments
(level3_variational) give effective coefficients 0.0831 ... 0.0782 at
l = 3..7, approaching this value from above as the higher-order
(l + 1/2)^-5 term dies away.  A is NOT (10/pi) <delta^3(r12)>; that
relation, used here previously, overstates the He coefficient by ~4.6x.

The cusp correction at finite l_max is the leading-order tail
    Delta E_cusp = -A * sum_{l=l_max+1}^{infinity} 1/(l + 1/2)^4 .

Only the He (Z = 2) coefficient is tabulated; no verified value is used for
other Z.  For Level 4 (H2) the module uses a MODEL: A(R) = A_He * delta(R) /
delta_He, i.e. the coefficient is taken proportional to the coalescence
density, with delta(R) from the hand-set interpolation in
coalescence_density_h2 (see its docstring).  Neither the proportionality nor
that interpolation is derived.

References:
  - Kato, Commun. Pure Appl. Math. 10, 151 (1957)
  - Schwartz, Phys. Rev. 126, 1015 (1962)
  - Hill, J. Chem. Phys. 83, 1173 (1985)
  - Kutzelnigg & Morgan, J. Chem. Phys. 96, 4484 (1992)
  - Thakkar & Smith, Phys. Rev. A 15, 1 (1977)
  - Bromley & Mitroy, arXiv:physics/0601049 (He partial-wave coefficients)
  - Paper 12, Section VII (cusp diagnosis)
  - Paper 13, Section III (angular eigenproblem)
"""

import numpy as np
from math import pi, sqrt
from typing import Tuple, Optional, Dict, List
from scipy.interpolate import CubicSpline


# ---------------------------------------------------------------------------
# Known exact coalescence densities for He-like systems
# <delta^3(r12)>_exact from Thakkar & Smith (1977) and related work
# ---------------------------------------------------------------------------

# For He (Z=2): <delta^3(r12)> = 0.1066 a.u.
# General Z: <delta^3(r12)> ~ Z^3 * 0.01333 for He-like (leading term)
# More precisely for neutral He: 0.1066 (Thakkar & Smith 1977)
_COALESCENCE_DENSITY_HE = 0.1066

# Schwartz exponent: E_{l_max} - E_exact ~ A / (l_max + 1/2)^4
# The 1/(l+1/2)^4 form is more accurate than 1/(l+1)^4 (Hill 1985)
_SCHWARTZ_EXPONENT = 4

# Leading partial-wave coefficient for the exact He ground state,
# A = 6 pi^2 int |Psi(r,r,0)|^2 r^5 dr (Hill 1985), value from
# Bromley & Mitroy, arXiv:physics/0601049, Eqs. 8-9 (A_E = -0.074226 there).
_SCHWARTZ_COEFF_HE = 0.074226


def _coefficient_from_delta(delta_r12: float) -> float:
    """Model coefficient A proportional to the coalescence density.

    A = A_He * delta / delta_He.  Exact at the He anchor by construction;
    elsewhere a model assumption (not derived), used for the H2 path and for
    the delta_r12 override.
    """
    return _SCHWARTZ_COEFF_HE * delta_r12 / _COALESCENCE_DENSITY_HE


def coalescence_density_he_like(Z: float) -> float:
    """Coalescence density <delta^3(r12)> for He-like atoms.

    Uses the exact value for He (Z=2) and Z^3 scaling for other charges.
    The Z^3 scaling comes from the hydrogenic cusp: Psi(r12=0) ~ Z^3/pi.

    Parameters
    ----------
    Z : float
        Nuclear charge.

    Returns
    -------
    float
        <delta^3(r12)> in atomic units.
    """
    return _COALESCENCE_DENSITY_HE * (Z / 2.0) ** 3


def schwartz_coefficient(Z: float) -> float:
    """Schwartz partial-wave convergence coefficient A for the He atom.

    The energy increment from partial wave l is approximately:
        Delta E_l ~ -A / (l + 1/2)^4

    with A = 6 pi^2 int |Psi(r,r,0)|^2 r^5 dr (Hill 1985); for the exact He
    ground state A = 0.074226 Ha (Bromley & Mitroy, arXiv:physics/0601049).
    Only Z = 2 is tabulated.

    Parameters
    ----------
    Z : float
        Nuclear charge (must be 2).

    Returns
    -------
    float
        Coefficient A in Hartree.
    """
    if abs(Z - 2.0) > 1e-12:
        raise ValueError(
            "schwartz_coefficient: only the He (Z=2) coefficient is "
            "tabulated; pass delta_r12 to use the proportional model"
        )
    return _SCHWARTZ_COEFF_HE


def cusp_correction_he(
    Z: float = 2.0,
    l_max: int = 0,
    delta_r12: Optional[float] = None,
) -> float:
    """Compute the cusp energy correction for He-like atoms at given l_max.

    Uses the Schwartz partial-wave convergence formula:
        Delta E_cusp = -sum_{l=l_max+1}^{infinity} A / (l + 1/2)^4

    The sum is computed exactly using the Hurwitz zeta function via direct
    summation (converges rapidly for the 1/n^4 tail).

    Parameters
    ----------
    Z : float
        Nuclear charge.
    l_max : int
        Maximum partial wave included in the calculation.
    delta_r12 : float or None
        Override coalescence density; A then follows the proportional
        model A = A_He * delta / delta_He.  If None, uses the He value.

    Returns
    -------
    float
        Energy correction Delta E_cusp (Hartree, negative).
    """
    if delta_r12 is not None:
        A = _coefficient_from_delta(delta_r12)
    else:
        A = schwartz_coefficient(Z)

    # Sum the tail: sum_{l=l_max+1}^{infinity} 1/(l + 1/2)^4
    # Converges rapidly; 200 terms is more than enough for 1/n^4
    tail = 0.0
    for l in range(l_max + 1, l_max + 201):
        tail += 1.0 / (l + 0.5) ** _SCHWARTZ_EXPONENT

    return -A * tail


def cusp_correction_he_extrapolated(
    energies_by_lmax: Dict[int, float],
    Z: float = 2.0,
) -> Tuple[float, float, float]:
    """Extrapolate the cusp correction from computed energies at multiple l_max.

    Fits the Schwartz formula E(l_max) = E_exact + A/(l_max + 1/2)^4
    to the computed energies and extracts E_exact.

    This is more reliable than the theoretical A coefficient because it
    uses the actual solver's convergence behavior.

    Parameters
    ----------
    energies_by_lmax : dict
        {l_max: energy} mapping.  Needs at least 2 entries.
    Z : float
        Nuclear charge (for reference only).

    Returns
    -------
    E_extrapolated : float
        Extrapolated CBS-limit energy.
    A_fitted : float
        Fitted Schwartz coefficient.
    residual : float
        RMS residual of the fit.
    """
    l_vals = sorted(energies_by_lmax.keys())
    if len(l_vals) < 2:
        raise ValueError("Need at least 2 l_max values for extrapolation")

    E_vals = np.array([energies_by_lmax[l] for l in l_vals])
    x_vals = np.array([1.0 / (l + 0.5) ** _SCHWARTZ_EXPONENT for l in l_vals])

    # Linear fit: E = E_exact + A * x
    # Using least squares: [1, x] @ [E_exact, A]^T = E
    design = np.column_stack([np.ones(len(l_vals)), x_vals])
    coeffs, residuals, _, _ = np.linalg.lstsq(design, E_vals, rcond=None)
    E_exact_fit = coeffs[0]
    A_fit = coeffs[1]

    E_predicted = design @ coeffs
    rms = np.sqrt(np.mean((E_vals - E_predicted) ** 2))

    return E_exact_fit, A_fit, rms


# ---------------------------------------------------------------------------
# Level 3: Cusp correction from coupled-channel solver output
# ---------------------------------------------------------------------------

def cusp_correction_from_solver(
    solver_result: dict,
    Z: float = 2.0,
    delta_r12: Optional[float] = None,
) -> dict:
    """Compute cusp correction from a Level 3 coupled-channel solver result.

    Parameters
    ----------
    solver_result : dict
        Output from solve_hyperspherical_algebraic_coupled().
        Must contain 'energy', 'l_max' keys.
    Z : float
        Nuclear charge.
    delta_r12 : float or None
        Override coalescence density.

    Returns
    -------
    dict
        Keys:
        - 'energy_uncorrected': original solver energy
        - 'delta_E_cusp': cusp correction (negative)
        - 'energy_corrected': corrected energy
        - 'error_uncorrected_pct': error vs exact He
        - 'error_corrected_pct': corrected error vs exact He
        - 'l_max': partial wave truncation
        - 'A_schwartz': Schwartz coefficient used
    """
    E_unc = solver_result['energy']
    l_max = solver_result['l_max']
    E_exact = -2.903724  # He exact (NIST)

    dE = cusp_correction_he(Z=Z, l_max=l_max, delta_r12=delta_r12)

    E_corr = E_unc + dE

    err_unc = abs(E_unc - E_exact) / abs(E_exact) * 100
    err_corr = abs(E_corr - E_exact) / abs(E_exact) * 100

    if delta_r12 is not None:
        A = _coefficient_from_delta(delta_r12)
    else:
        A = schwartz_coefficient(Z)

    return {
        'energy_uncorrected': E_unc,
        'delta_E_cusp': dE,
        'energy_corrected': E_corr,
        'error_uncorrected_pct': err_unc,
        'error_corrected_pct': err_corr,
        'l_max': l_max,
        'A_schwartz': A,
    }


# ---------------------------------------------------------------------------
# Level 4: R-dependent cusp correction for H2 PES
# ---------------------------------------------------------------------------

def coalescence_density_h2(R: float) -> float:
    """Estimate the coalescence density <delta^3(r12)>(R) for H2.

    Uses a simple interpolation between the united-atom (He) and
    separated-atom (2*H) limits:

    R -> 0:  <delta^3(r12)> -> <delta^3(r12)>_He (united atom, Z=2)
    R -> inf: <delta^3(r12)> -> 0 (two separated H atoms, no e-e coalescence
              in ground state singlet -- each atom has only 1 electron)

    The interpolation f(R) = exp(-beta R^2 / (1 + gamma R)) has two
    HAND-SET parameters, beta = 0.60 and gamma = 0.3, set by hand with the
    Gaussian-geminal H2 coalescence density (Komasa & Cencek, ~0.025-0.030
    a.u. near R = 1.4, i.e. f ~ 0.25) as the stated target.  They are not
    derived, and they do not reproduce that target: f(1.4) = 0.437,
    delta(1.4) = 0.0466 a.u.  They enter the H2 cusp correction directly.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).

    Returns
    -------
    float
        <delta^3(r12)>(R) in atomic units.
    """
    # United atom limit (He)
    delta_He = coalescence_density_he_like(Z=2.0)

    # At R=0 (united atom), both electrons on same nucleus -> He-like
    # At R->inf, electrons on separate atoms -> delta=0 for singlet H+H
    # The decay scale is set by the electron overlap, which for H2
    # goes roughly as exp(-R) for large R.
    # Model: delta(R) = delta_He * f(R) with f(0) = 1, f(inf) = 0,
    # f(R) = exp(-beta * R^2 / (1 + gamma*R)).  beta, gamma are HAND-SET
    # (target f(1.4) ~ 0.25 from Komasa & Cencek; actual f(1.4) = 0.437).
    beta = 0.60
    gamma = 0.3
    f_R = np.exp(-beta * R ** 2 / (1.0 + gamma * R))

    return delta_He * f_R


def cusp_correction_h2_point(
    R: float,
    l_max: int = 2,
    delta_r12: Optional[float] = None,
) -> float:
    """Compute the cusp energy correction for H2 at a single R.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    l_max : int
        Maximum partial wave in the Level 4 angular expansion.
    delta_r12 : float or None
        Override coalescence density at this R.

    Returns
    -------
    float
        Delta E_cusp(R) in Hartree (negative).
    """
    if delta_r12 is None:
        delta_r12 = coalescence_density_h2(R)

    # Schwartz formula with the model coefficient A(R) = A_He * delta(R)/delta_He
    A_R = _coefficient_from_delta(delta_r12)

    # Sum the tail
    tail = 0.0
    for l in range(l_max + 1, l_max + 201):
        tail += 1.0 / (l + 0.5) ** _SCHWARTZ_EXPONENT

    return -A_R * tail


def cusp_correction_h2_pes(
    R_grid: np.ndarray,
    E_pes: np.ndarray,
    l_max: int = 2,
    V_NN: Optional[np.ndarray] = None,
) -> dict:
    """Apply cusp correction to the H2 potential energy surface.

    Parameters
    ----------
    R_grid : ndarray
        Internuclear distances (bohr).
    E_pes : ndarray
        Total electronic + nuclear energy at each R.
    l_max : int
        Maximum partial wave used in the Level 4 solver.
    V_NN : ndarray or None
        Nuclear repulsion at each R.  If None, computes 1/R.

    Returns
    -------
    dict
        Keys:
        - 'R_grid': input R grid
        - 'E_uncorrected': original PES
        - 'delta_E_cusp': cusp correction at each R
        - 'E_corrected': corrected PES
        - 'D_e_uncorrected': uncorrected dissociation energy
        - 'D_e_corrected': corrected dissociation energy
        - 'D_e_exact': exact D_e for H2 (0.17447 Ha)
        - 'D_e_pct_uncorrected': % of exact D_e recovered (uncorrected)
        - 'D_e_pct_corrected': % of exact D_e recovered (corrected)
    """
    # Compute cusp correction at each R
    dE_cusp = np.array([cusp_correction_h2_point(R, l_max) for R in R_grid])

    E_corr = E_pes + dE_cusp

    # Dissociation energies
    # E(R->inf) for H2 is -1.0 Ha (two isolated H atoms)
    # Use the last few points to estimate dissociation limit
    E_inf = E_pes[-1] if R_grid[-1] > 8.0 else -1.0
    E_inf_corr = E_corr[-1] if R_grid[-1] > 8.0 else -1.0

    D_e_unc = E_inf - np.min(E_pes)
    D_e_corr = E_inf_corr - np.min(E_corr)

    D_e_exact = 0.17447  # Ha (Kolos & Wolniewicz)

    return {
        'R_grid': R_grid,
        'E_uncorrected': E_pes,
        'delta_E_cusp': dE_cusp,
        'E_corrected': E_corr,
        'D_e_uncorrected': D_e_unc,
        'D_e_corrected': D_e_corr,
        'D_e_exact': D_e_exact,
        'D_e_pct_uncorrected': D_e_unc / D_e_exact * 100,
        'D_e_pct_corrected': D_e_corr / D_e_exact * 100,
        'l_max': l_max,
    }


# ---------------------------------------------------------------------------
# Combined interface: correction from Level 4 solver output
# ---------------------------------------------------------------------------

def cusp_correction_from_level4(
    level4_result: dict,
    l_max: Optional[int] = None,
) -> dict:
    """Apply cusp correction to a Level 4 solver result.

    Parameters
    ----------
    level4_result : dict
        Output from solve_level4_h2_multichannel().
        Must contain 'E_total' and 'R' keys.
    l_max : int or None
        Override l_max. If None, uses level4_result['l_max'].

    Returns
    -------
    dict
        Keys:
        - 'E_uncorrected': original total energy
        - 'delta_E_cusp': cusp correction at this R
        - 'E_corrected': corrected total energy
        - 'R': internuclear distance
        - 'l_max': partial wave truncation
    """
    if l_max is None:
        l_max = level4_result['l_max']
    R = level4_result['R']
    E_unc = level4_result['E_total']

    dE = cusp_correction_h2_point(R, l_max)
    E_corr = E_unc + dE

    return {
        'E_uncorrected': E_unc,
        'delta_E_cusp': dE,
        'E_corrected': E_corr,
        'R': R,
        'l_max': l_max,
    }


# ---------------------------------------------------------------------------
# Utility: l_max convergence analysis
# ---------------------------------------------------------------------------

def analyze_cusp_convergence(
    Z: float = 2.0,
    l_max_range: Optional[List[int]] = None,
) -> dict:
    """Analyze how the cusp correction scales with l_max.

    Useful for verifying the 1/(l+1/2)^4 convergence and estimating
    how much error remains at each l_max.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    l_max_range : list of int or None
        l_max values to analyze. Default [0, 1, 2, 3, 4, 5].

    Returns
    -------
    dict
        Keys:
        - 'l_max_values': input l_max range
        - 'corrections': cusp correction at each l_max
        - 'A_schwartz': theoretical Schwartz coefficient
        - 'corrections_mHa': corrections in milli-Hartree
    """
    if l_max_range is None:
        l_max_range = [0, 1, 2, 3, 4, 5]

    A = schwartz_coefficient(Z)
    corrections = [cusp_correction_he(Z=Z, l_max=l) for l in l_max_range]

    return {
        'l_max_values': l_max_range,
        'corrections': corrections,
        'corrections_mHa': [c * 1000 for c in corrections],
        'A_schwartz': A,
    }
