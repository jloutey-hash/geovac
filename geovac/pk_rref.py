"""
Ab initio R_ref for R-dependent Phillips-Kleinman pseudopotential.

The R-dependent PK weights the standard Gaussian PK barrier by a factor
that depends on the internuclear distance R:

    w_PK(R) = delta_{l,0} * min(cap, R / R_ref)

where R_ref is derived from atomic core properties (zero fitting).

This module provides three ab initio R_ref candidates:

(a) r_half_screen: distance where Z_eff(r) = Z - 1 (half the core
    electrons screened). This is where the core "edge" is — the natural
    length scale for core-valence separation.

(b) r_pk_width: 1/sqrt(B) from the ab initio PK Gaussian exponent B.
    This is the effective width of the PK barrier in real space.

(c) r_avg: mean core radius <r> = integral r * n(r) dr / integral n(r) dr.
    The first moment of the core density.

References:
    Paper 17, Sec VI.A (l_max divergence diagnosis)
    CLAUDE.md Section 2 (active frontier, R-dependent PK)
"""

import numpy as np
from typing import Dict, Optional, Tuple

from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK


def compute_rref_candidates(
    core: CoreScreening,
    pk: Optional[AbInitioPK] = None,
    n_core: int = 2,
) -> Dict[str, float]:
    """
    Compute ab initio R_ref candidates from solved core screening data.

    All candidates are derived from the atomic core wavefunction with
    zero fitted parameters.

    Parameters
    ----------
    core : CoreScreening
        Solved core screening object.
    pk : AbInitioPK or None
        Ab initio PK object. If None, created from core.
    n_core : int
        Number of core electrons.

    Returns
    -------
    dict
        Mapping candidate name -> R_ref value in bohr.
    """
    core._check_solved()

    if pk is None:
        pk = AbInitioPK(core, n_core=n_core)

    r_grid = core.r_grid
    n_r = core.n_r

    candidates: Dict[str, float] = {}

    # (a) r_half_screen: where Z_eff = Z - n_core/2
    # For a 2-electron core, this is where Z_eff = Z - 1.
    z_half_target = core.Z - n_core / 2.0
    r_test = np.linspace(0.05, 5.0, 500)
    z_vals = core.z_eff(r_test)

    # Find zero crossing of z_eff - z_half_target
    diff = z_vals - z_half_target
    sign_changes = np.where(np.diff(np.sign(diff)))[0]
    if len(sign_changes) > 0:
        from scipy.optimize import brentq
        from scipy.interpolate import CubicSpline
        z_spline = CubicSpline(r_test, z_vals)
        idx = sign_changes[0]
        r_half = float(brentq(
            lambda x: z_spline(x) - z_half_target,
            r_test[idx], r_test[idx + 1],
        ))
        candidates['r_half_screen'] = r_half

    # (b) r_pk_width: 1/sqrt(B) from the selected PK Gaussian exponent
    candidates['r_pk_width'] = 1.0 / np.sqrt(pk.B)

    # Also store the inv2-weighted core radius directly
    candidates['r_inv2'] = pk.r_core

    # (c) r_avg: mean core radius
    total = np.trapezoid(n_r, r_grid)
    if total > 0:
        r_avg = np.trapezoid(r_grid * n_r, r_grid) / total
        candidates['r_avg'] = r_avg

    # (d) r_rms: RMS core radius
    if total > 0:
        r2_avg = np.trapezoid(r_grid**2 * n_r, r_grid) / total
        candidates['r_rms'] = np.sqrt(r2_avg)

    # (e) r_median: where N_core = n_core/2
    N_core = core.N_core
    half = n_core / 2.0
    idx_cross = np.searchsorted(N_core, half)
    if 1 < idx_cross < len(r_grid) - 1:
        from scipy.optimize import brentq
        from scipy.interpolate import CubicSpline
        N_spline = CubicSpline(r_grid, N_core)
        r_med = float(brentq(
            lambda x: N_spline(x) - half,
            r_grid[idx_cross - 1], r_grid[idx_cross + 1],
        ))
        candidates['r_median'] = r_med

    return candidates


def select_rref(
    core: CoreScreening,
    pk: Optional[AbInitioPK] = None,
    n_core: int = 2,
    method: str = 'r_half_screen',
) -> float:
    """
    Select R_ref from core screening data.

    Parameters
    ----------
    core : CoreScreening
        Solved core screening object.
    pk : AbInitioPK or None
        Ab initio PK object.
    n_core : int
        Number of core electrons.
    method : str
        Which candidate to use: 'r_half_screen', 'r_pk_width', 'r_avg',
        'r_rms', 'r_median', 'r_inv2'.

    Returns
    -------
    float
        R_ref in bohr.
    """
    candidates = compute_rref_candidates(core, pk, n_core)
    if method not in candidates:
        raise ValueError(
            f"R_ref method '{method}' not available. "
            f"Available: {list(candidates.keys())}"
        )
    return candidates[method]


def pk_weight_r_dependent(
    R: float,
    R_ref: float,
    cap: float = 1.0,
) -> float:
    """
    R-dependent PK weight factor.

    w_PK(R) = min(cap, R / R_ref)

    At small R (R < R_ref * cap), the PK barrier is reduced because
    the core is more compressed and the valence-core overlap region
    shrinks. At large R (R > R_ref * cap), the weight saturates at cap.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    R_ref : float
        Reference distance from core screening (bohr).
    cap : float
        Maximum weight (default 1.0).

    Returns
    -------
    float
        PK weight factor.
    """
    return min(cap, R / R_ref)
