"""
Log-Holonomy Module for GeoVac
================================

Computes discrete log-holonomies on the geometric lattice by circulating
around plaquettes in the (n, m) plane.  This measures the **curvature of
the edge weight function**, not a genuine Berry phase (which requires
complex-valued transition amplitudes).

A plaquette at (n, l, m) is the closed loop:

    |n,l,m> --T+--> |n+1,l,m> --L+--> |n+1,l,m+1> --T---> |n,l,m+1> --L---> |n,l,m>

For topological weights (w = 1/(n1*n2)), the T+ and T- edges share the same
weight 1/(n*(n+1)), so they cancel.  The net plaquette phase reduces to:

    theta(n) = ln(w_L(n+1)) - ln(w_L(n))
             = ln(1/(n+1)^2) - ln(1/n^2)
             = -2 * ln((n+1)/n)

The per-plaquette phase decays as ~2/n, giving exponent k = 1.0 exactly.

Note: This is distinct from Paper 1's claimed Berry phase (arg of complex
holonomy products from SU(2)/SU(1,1) transition operators).  The arg()
formula produces identically zero for real-valued CG coefficients.  See
debug/qa_sprint/berry_phase_reconciliation.md for full analysis.

For binary weights (all weights = 1), every plaquette phase is identically zero
since ln(1) - ln(1) = 0.

Author: GeoVac
Date: March 2026
"""

import numpy as np
from scipy.stats import linregress
from typing import Dict, List, Tuple

from geovac.lattice import GeometricLattice


def compute_plaquette_phase(lattice: GeometricLattice, n: int, l: int, m: int) -> float:
    """
    Compute the Berry phase around a single plaquette in the (n, m) plane.

    The plaquette loop is:
        |n,l,m> -> |n+1,l,m> -> |n+1,l,m+1> -> |n,l,m+1> -> |n,l,m>

    The phase is computed from the log of adjacency matrix weights:
        theta = ln(w_12) + ln(w_23) - ln(w_34) - ln(w_41)

    where the sign convention treats the first two edges as "forward" and the
    last two as "backward" (returning to the origin).

    Parameters
    ----------
    lattice : GeometricLattice
        The lattice to compute the plaquette phase on.
    n : int
        Principal quantum number of the plaquette origin.
    l : int
        Angular momentum quantum number.
    m : int
        Magnetic quantum number of the plaquette origin.
        The plaquette spans m and m+1.

    Returns
    -------
    float
        The plaquette phase (discrete holonomy) in units of radians-equivalent.
        Returns 0.0 if any of the four corner states do not exist in the lattice.

    Raises
    ------
    ValueError
        If n < 1.
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")

    # The four corners of the plaquette
    s1 = (n, l, m)
    s2 = (n + 1, l, m)
    s3 = (n + 1, l, m + 1)
    s4 = (n, l, m + 1)

    idx = lattice._state_index

    # Check all four states exist
    if any(s not in idx for s in (s1, s2, s3, s4)):
        return 0.0

    i1, i2, i3, i4 = idx[s1], idx[s2], idx[s3], idx[s4]

    # Extract adjacency weights for each edge
    w_12 = lattice.adjacency[i1, i2]  # T+: (n,l,m) -> (n+1,l,m)
    w_23 = lattice.adjacency[i2, i3]  # L+: (n+1,l,m) -> (n+1,l,m+1)
    w_34 = lattice.adjacency[i3, i4]  # T-: (n+1,l,m+1) -> (n,l,m+1)
    w_41 = lattice.adjacency[i4, i1]  # L-: (n,l,m+1) -> (n,l,m)

    # All weights must be positive for the log to be defined
    weights = [w_12, w_23, w_34, w_41]
    if any(w <= 0.0 for w in weights):
        return 0.0

    # Discrete holonomy: forward edges minus backward edges
    theta = np.log(w_12) + np.log(w_23) - np.log(w_34) - np.log(w_41)
    return float(theta)


def compute_berry_phase_at_n(lattice: GeometricLattice, n: int) -> float:
    """
    Compute the total Berry phase at principal quantum number n.

    Sums the plaquette phases over all valid (l, m) pairs where:
    - l ranges from 0 to n-1 (all angular momenta at shell n)
    - m ranges from -l to l-1 (both m and m+1 must be valid)
    - n+1 must be <= max_n (the plaquette needs the next shell)

    Parameters
    ----------
    lattice : GeometricLattice
        The lattice to compute the Berry phase on.
    n : int
        Principal quantum number.

    Returns
    -------
    float
        Total Berry phase (sum of all plaquette phases) at shell n.

    Raises
    ------
    ValueError
        If n < 1 or n >= lattice.max_n.
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    if n >= lattice.max_n:
        raise ValueError(
            f"n must be < max_n={lattice.max_n} (need shell n+1 for plaquettes), got {n}"
        )

    total = 0.0
    for l in range(n):  # l = 0, 1, ..., n-1
        for m in range(-l, l):  # m = -l, ..., l-1 (m+1 must be <= l)
            total += compute_plaquette_phase(lattice, n, l, m)
    return total


def fit_power_law(n_values: np.ndarray, phase_values: np.ndarray) -> Tuple[float, float]:
    """
    Fit a power law theta = A * n^(-k) using log-log linear regression.

    Parameters
    ----------
    n_values : np.ndarray
        Array of principal quantum numbers (must be > 0).
    phase_values : np.ndarray
        Array of Berry phase magnitudes (must be > 0 for log fitting).

    Returns
    -------
    tuple[float, float]
        (k, A) where theta ~ A * n^(-k).
        k is the power-law exponent (positive for decay).
        A is the amplitude prefactor.

    Raises
    ------
    ValueError
        If fewer than 2 data points, or if any values are non-positive.
    """
    n_arr = np.asarray(n_values, dtype=float)
    p_arr = np.asarray(phase_values, dtype=float)

    if len(n_arr) < 2:
        raise ValueError("Need at least 2 data points for power-law fit")

    # Use absolute values for fitting (phases may be negative)
    abs_phases = np.abs(p_arr)

    # Filter out zeros
    mask = abs_phases > 0
    if np.sum(mask) < 2:
        raise ValueError("Need at least 2 non-zero phase values for power-law fit")

    log_n = np.log(n_arr[mask])
    log_p = np.log(abs_phases[mask])

    slope, intercept, _, _, _ = linregress(log_n, log_p)

    k = -slope  # theta ~ n^(-k) => log(theta) = -k*log(n) + log(A)
    A = np.exp(intercept)

    return (k, A)


def compute_phase_convergence(
    max_n_values: List[int],
    Z: int = 1,
    topological_weights: bool = True,
) -> Dict[str, object]:
    """
    Study convergence of Berry phase power-law exponent with lattice size.

    For each max_n in max_n_values, builds a GeometricLattice, computes the
    Berry phase at each shell n = 1, ..., max_n-1, and fits a power law.

    Paper 1 predicts the exponent k converges to approximately 2.113.

    Parameters
    ----------
    max_n_values : List[int]
        List of max_n values to test (each must be >= 3).
    Z : int, optional
        Nuclear charge (default: 1).
    topological_weights : bool, optional
        Whether to use topological edge weights (default: True).
        With binary weights, all phases are zero and fitting will fail.

    Returns
    -------
    dict
        Keys:
        - 'max_n_values': list of max_n tested
        - 'phases': list of np.ndarray, each containing Berry phases at
          n = 1, ..., max_n-1
        - 'n_arrays': list of np.ndarray, the n values for each max_n
        - 'exponents': list of fitted k values (nan if fit fails)
        - 'amplitudes': list of fitted A values (nan if fit fails)
    """
    all_phases: List[np.ndarray] = []
    all_n_arrays: List[np.ndarray] = []
    exponents: List[float] = []
    amplitudes: List[float] = []

    for max_n in max_n_values:
        if max_n < 3:
            # Need at least n=1,2 for a 2-point fit
            all_phases.append(np.array([]))
            all_n_arrays.append(np.array([]))
            exponents.append(float('nan'))
            amplitudes.append(float('nan'))
            continue

        lattice = GeometricLattice(
            max_n=max_n,
            nuclear_charge=Z,
            topological_weights=topological_weights,
        )

        n_vals = np.arange(1, max_n)
        phases = np.array([compute_berry_phase_at_n(lattice, n) for n in n_vals])

        all_phases.append(phases)
        all_n_arrays.append(n_vals)

        try:
            k, A = fit_power_law(n_vals, phases)
            exponents.append(k)
            amplitudes.append(A)
        except ValueError:
            exponents.append(float('nan'))
            amplitudes.append(float('nan'))

    return {
        'max_n_values': list(max_n_values),
        'phases': all_phases,
        'n_arrays': all_n_arrays,
        'exponents': exponents,
        'amplitudes': amplitudes,
    }
