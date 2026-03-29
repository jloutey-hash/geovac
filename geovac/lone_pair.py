"""
Lone pair solver for composed natural geometries.

Wraps the Level 3 hyperspherical solver (Paper 13) for use as lone pair
blocks in polyatomic molecules like H₂O.  Each lone pair is a two-electron
system in the field of a screened nucleus (Z_eff from core screening).

Channel data adapter
--------------------
The inter-fiber coupling functions (inter_fiber_coupling.py) expect channel
data in Level 4 format: channels labelled by (l1, l2) tuples.  The Level 3
solver uses a single l index (l1 = l2 = l in the L=0 sector).  The adapter
``extract_channel_data_level3`` maps Level 3 channel l → (l, l) so that
overlap and exchange formulas work for both Level 3 and Level 4 blocks.

Inter-fiber coupling
--------------------
Generic overlap between any two channel_data dicts (Level 3 or Level 4):
    S(θ) = Σ_{(l1,l2)} phase(l1,l2,θ) · c^A_{l1,l2} · c^B_{l1,l2}
Only channels present in *both* fibers contribute.

References:
  - Paper 13: Hyperspherical lattice (He at 0.05%)
  - Paper 17: Composed natural geometries
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import eval_legendre
from typing import Any, Dict, Optional, Tuple

from geovac.hyperspherical_radial import solve_helium
from geovac.hyperspherical_angular import solve_angular


# ---------------------------------------------------------------------------
# 1. Lone pair solver wrapper
# ---------------------------------------------------------------------------

def solve_lone_pair(
    Z_eff: float,
    l_max: int = 3,
    n_alpha: int = 100,
    n_Re: int = 2000,
    N_R_angular: int = 200,
    R_max: float = 30.0,
    n_channels: int = 1,
    verbose: bool = False,
) -> dict:
    """
    Solve a lone pair as a two-electron system with effective nuclear charge.

    Thin wrapper around solve_helium (Level 3 hyperspherical solver).
    For He verification, call with Z_eff=2.

    Parameters
    ----------
    Z_eff : float
        Effective nuclear charge seen by the lone pair electrons.
        For O lone pairs: Z_eff ≈ 6 (Z_O - 2 core electrons).
    l_max : int
        Maximum partial wave in angular expansion.
    n_alpha : int
        Number of FD grid points for hyperangle alpha.
    n_Re : int
        Number of grid points for the radial (hyperradial) solve.
    N_R_angular : int
        Number of R points for computing adiabatic curves.
    R_max : float
        Maximum hyperradius (bohr).  Scale with 1/Z_eff for heavier atoms.
    n_channels : int
        Number of adiabatic channels.
    verbose : bool
        Print progress.

    Returns
    -------
    dict
        Same keys as solve_helium():
        - energy : float (ground state energy in Ha)
        - wavefunction : ndarray
        - R_grid_radial, R_grid_angular : ndarray
        - mu_curves, V_eff : ndarray
        - Z, l_max, n_channels, coupled : metadata
    """
    # For higher-Z atoms, the wavefunction is more compact.
    # Scale R_max down so the grid resolves the relevant region.
    scaled_R_max = R_max * (2.0 / Z_eff) if Z_eff > 2.0 else R_max

    return solve_helium(
        Z=Z_eff,
        l_max=l_max,
        n_alpha=n_alpha,
        N_R_angular=N_R_angular,
        R_max=scaled_R_max,
        N_R_radial=n_Re,
        n_channels=n_channels,
        coupled=False,
        verbose=verbose,
    )


# ---------------------------------------------------------------------------
# 2. Channel data extraction (Level 3 → common format)
# ---------------------------------------------------------------------------

def extract_channel_data_level3(
    helium_result: dict,
    Z: float,
    l_max: int,
    n_alpha: int,
    n_sample_R: int = 15,
) -> Dict[str, Any]:
    """
    Extract channel data from a Level 3 (hyperspherical) result.

    Produces the same dict format as ``extract_channel_data`` in
    inter_fiber_coupling.py, so that overlap, density, and exchange
    functions can consume Level 3 and Level 4 data interchangeably.

    Key convention mapping:
        Level 3 channel ``l`` → Level 4 channel ``(l, l)``
    Only diagonal (l, l) channels exist in Level 3 (L=0 sector: l1 = l2).

    Parameters
    ----------
    helium_result : dict
        Output of solve_helium() / solve_lone_pair().
    Z : float
        Nuclear charge used in the solve (Z_eff for lone pairs).
    l_max : int
        Maximum partial wave (must match the solve).
    n_alpha : int
        Number of alpha FD grid points (must match the solve).
    n_sample_R : int
        Number of hyperradius sample points within the significant
        wavefunction region.

    Returns
    -------
    data : dict
        Same keys as extract_channel_data():
        - Re_samples : ndarray (n_sample,)
        - F_values : ndarray (n_sample,)
        - alpha_grid : ndarray (n_alpha,)
        - h_alpha : float
        - channels : list of (l, l) tuples
        - n_ch : int (= l_max + 1)
        - vec_2d_list : list of ndarray (n_ch, n_alpha)
        - ch_weights : ndarray (n_sample, n_ch)
        - ang_density : ndarray (n_sample, n_alpha)
        - dRe : float
        - raw_norm : float
    """
    F = helium_result['wavefunction']       # shape (N_R_radial,)
    R_grid = helium_result['R_grid_radial']
    R_ang = helium_result['R_grid_angular']

    n_l = l_max + 1
    # Channel labels: (l, l) for l = 0, 1, ..., l_max
    channels = [(l, l) for l in range(n_l)]
    n_ch = n_l

    # Find significant region of |F|^2
    F2 = F ** 2
    threshold = 0.001 * F2.max()
    significant = np.where(F2 > threshold)[0]
    if len(significant) < 2:
        i_peak = np.argmax(F2)
        significant = np.arange(
            max(0, i_peak - 50), min(len(F2), i_peak + 50))

    R_lo = R_grid[significant[0]]
    R_hi = R_grid[significant[-1]]
    Re_samples = np.linspace(R_lo, R_hi, n_sample_R)

    # Alpha grid (must match solve_angular conventions)
    h_alpha = (np.pi / 2) / (n_alpha + 1)
    alpha_grid = (np.arange(n_alpha) + 1) * h_alpha

    # Interpolate F to sample points
    F_spline = CubicSpline(R_grid, F, extrapolate=True)
    F_values = F_spline(Re_samples)

    vec_2d_list = []
    ch_weights_all = []
    ang_density_all = []

    for k, R_s in enumerate(Re_samples):
        if R_s < 1e-8:
            vec_2d_list.append(np.zeros((n_ch, n_alpha)))
            ch_weights_all.append(np.zeros(n_ch))
            ang_density_all.append(np.zeros(n_alpha))
            continue

        # Solve angular problem at this R
        _, vecs = solve_angular(R_s, Z, l_max, n_alpha, n_channels=1)
        u = vecs[0]  # shape (n_l * n_alpha,)

        # Reshape: each Level 3 channel l maps to one row
        # u[l*n_alpha : (l+1)*n_alpha] is the alpha amplitude for channel l
        v2d = u.reshape(n_l, n_alpha)
        vec_2d_list.append(v2d)

        # Channel weights: |c_l|^2 = sum_j |v2d[l, j]|^2
        cw = np.sum(v2d ** 2, axis=1)
        ch_weights_all.append(cw)

        # Angular density at each alpha: sum_l |u_l(alpha)|^2
        ang_d = np.sum(v2d ** 2, axis=0)
        ang_density_all.append(ang_d)

    ch_weights_arr = np.array(ch_weights_all)   # (n_sample, n_ch)
    ang_density_arr = np.array(ang_density_all)  # (n_sample, n_alpha)

    # Raw normalization (same formula as inter_fiber_coupling)
    dRe = Re_samples[1] - Re_samples[0] if len(Re_samples) > 1 else 1.0
    raw_norm = 0.0
    for k in range(len(Re_samples)):
        F2_k = F_values[k] ** 2
        raw_norm += F2_k * np.sum(ang_density_all[k]) * h_alpha * dRe

    return {
        'Re_samples': Re_samples,
        'F_values': F_values,
        'alpha_grid': alpha_grid,
        'h_alpha': h_alpha,
        'channels': channels,
        'n_ch': n_ch,
        'vec_2d_list': vec_2d_list,
        'ch_weights': ch_weights_arr,
        'ang_density': ang_density_arr,
        'dRe': dRe,
        'raw_norm': raw_norm,
    }


# ---------------------------------------------------------------------------
# 3. Generic inter-fiber overlap (Level 3 ↔ Level 4, Level 3 ↔ Level 3)
# ---------------------------------------------------------------------------

def _common_channels(
    channels_A: list,
    channels_B: list,
) -> list:
    """Return list of (l1, l2) channels present in both A and B."""
    set_A = set(tuple(c) for c in channels_A)
    set_B = set(tuple(c) for c in channels_B)
    return sorted(set_A & set_B)


def compute_inter_fiber_overlap(
    channel_data_A: Dict[str, Any],
    channel_data_B: Dict[str, Any],
    bond_angle: float = np.pi / 2,
) -> Dict[str, Any]:
    """
    Compute inter-fiber channel overlap S(θ) between two fibers.

    Works for any combination of Level 3 and Level 4 channel data
    (produced by extract_channel_data or extract_channel_data_level3).

    The overlap at each sample point:
        S(θ) = Σ_{(l1,l2)∈common} P_{l1}(cos θ) · P_{l2}(cos θ)
               · c^A_{l1,l2} · c^B_{l1,l2}

    where c_{l1,l2} = |channel weight|^2, averaged over the wavefunction.

    For Level 3 ↔ Level 4: only the diagonal (l, l) channels are common.
    For Level 3 ↔ Level 3: all (l, l) channels are common.
    For Level 4 ↔ Level 4: all (l1, l2) channels are common.

    Parameters
    ----------
    channel_data_A, channel_data_B : dict
        Output of extract_channel_data() or extract_channel_data_level3().
    bond_angle : float
        Angle between the two fiber axes (radians).

    Returns
    -------
    result : dict
        S_avg : float — |F|²-weighted average overlap
        common_channels : list of (l1, l2) tuples that contributed
        per_channel_contribution : dict mapping (l1,l2) → weighted contrib
    """
    channels_A = channel_data_A['channels']
    channels_B = channel_data_B['channels']
    common = _common_channels(channels_A, channels_B)

    if not common:
        return {
            'S_avg': 0.0,
            'common_channels': [],
            'per_channel_contribution': {},
        }

    # Index lookups
    idx_A = {tuple(c): i for i, c in enumerate(channels_A)}
    idx_B = {tuple(c): i for i, c in enumerate(channels_B)}

    # Compute |F|²-weighted channel weights for each fiber
    def _weighted_channel_weights(cd: Dict[str, Any]) -> np.ndarray:
        """Return |F|²-weighted average channel weights, shape (n_ch,)."""
        F2 = cd['F_values'] ** 2
        total_F2 = np.sum(F2)
        if total_F2 < 1e-30:
            return np.zeros(cd['n_ch'])
        # ch_weights shape: (n_sample, n_ch)
        return np.sum(cd['ch_weights'] * F2[:, np.newaxis], axis=0) / total_F2

    w_A = _weighted_channel_weights(channel_data_A)
    w_B = _weighted_channel_weights(channel_data_B)

    cos_theta = np.cos(bond_angle)
    S_total = 0.0
    per_ch = {}

    for ch in common:
        l1, l2 = ch
        phase = float(
            eval_legendre(l1, cos_theta) * eval_legendre(l2, cos_theta)
        )
        contrib = phase * w_A[idx_A[ch]] * w_B[idx_B[ch]]
        per_ch[ch] = float(contrib)
        S_total += contrib

    return {
        'S_avg': float(S_total),
        'common_channels': common,
        'per_channel_contribution': per_ch,
    }


def compute_inter_fiber_overlap_detailed(
    channel_data_A: Dict[str, Any],
    channel_data_B: Dict[str, Any],
    bond_angle: float = np.pi / 2,
) -> Dict[str, Any]:
    """
    Compute per-sample-point overlap S(R_e, θ) between two fibers.

    More detailed than compute_inter_fiber_overlap: returns S at each
    sample R_e point of fiber A, using the nearest-R_e channel weights
    from fiber B.  This is useful for diagnostics.

    Parameters
    ----------
    channel_data_A, channel_data_B : dict
        Channel data from either Level 3 or Level 4.
    bond_angle : float
        Angle between fiber axes (radians).

    Returns
    -------
    result : dict
        S_avg : float
        S_samples_A : ndarray — S at each of fiber A's sample points
        Re_samples_A : ndarray — fiber A's R_e sample points
        common_channels : list
    """
    channels_A = channel_data_A['channels']
    channels_B = channel_data_B['channels']
    common = _common_channels(channels_A, channels_B)

    idx_A = {tuple(c): i for i, c in enumerate(channels_A)}
    idx_B = {tuple(c): i for i, c in enumerate(channels_B)}

    cos_theta = np.cos(bond_angle)
    phases = {
        ch: float(eval_legendre(ch[0], cos_theta) * eval_legendre(ch[1], cos_theta))
        for ch in common
    }

    # Use fiber B's |F|²-weighted average channel weights (single set)
    F2_B = channel_data_B['F_values'] ** 2
    total_F2_B = np.sum(F2_B)
    if total_F2_B > 1e-30:
        w_B = np.sum(
            channel_data_B['ch_weights'] * F2_B[:, np.newaxis], axis=0
        ) / total_F2_B
    else:
        w_B = np.zeros(channel_data_B['n_ch'])

    # Compute S at each of fiber A's sample points
    n_sample_A = len(channel_data_A['Re_samples'])
    S_samples = np.zeros(n_sample_A)

    for k in range(n_sample_A):
        cw_A = channel_data_A['ch_weights'][k]  # (n_ch_A,)
        for ch in common:
            S_samples[k] += phases[ch] * cw_A[idx_A[ch]] * w_B[idx_B[ch]]

    # |F_A|²-weighted average
    F2_A = channel_data_A['F_values'] ** 2
    total_F2_A = np.sum(F2_A)
    if total_F2_A > 1e-30:
        S_avg = np.sum(S_samples * F2_A) / total_F2_A
    else:
        S_avg = 0.0

    return {
        'S_avg': float(S_avg),
        'S_samples_A': S_samples,
        'Re_samples_A': channel_data_A['Re_samples'],
        'common_channels': common,
    }
