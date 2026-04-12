"""
Monopole inter-fiber coupling for polyatomic composed geometries.

Phase 1 implementation: k=0 (monopole) Coulomb coupling between bond-pair
fibers sharing a common central atom.

Physics: After solving each bond pair independently (Level 4), extract the
one-electron radial density from each fiber's wavefunction. Transform to
coordinates centered at the shared atom (Be). Compute the Slater F^0
(monopole) integral between the two fiber densities.

The monopole energy:
    E_monopole(R) = integral integral rho_A(d1) rho_B(d2) / max(d1, d2) dd1 dd2

where rho_A, rho_B are spherically-averaged electron densities at distance
d from the shared center. For identical fibers (D_inf_h), rho_A = rho_B.

Two computational pathways:

  method='numerical' (original):
    1. Re-run angular solver at sample R_e points to extract angular density
    2. Histogram-bin the one-electron density P(r) on a radial grid
    3. Shell theorem convolution to Be-centered coordinates
    4. Cumulative-charge F^0 integral

  method='algebraic' (v2.0.1):
    1. Single call to extract_channel_data() provides angular eigenvectors,
       channel weights, and angular densities at all sample R_e points
    2. Proper quadrature with spline-interpolated angular density replaces
       histogram binning for P(r)
    3. Channel decomposition: F^0 = sum_{ch,ch'} M(ch, ch') with per-channel
       densities P_ch(r) and F^0 matrix M(ch, ch')
    4. S(R) and F^0(R) computed from the same channel data (no duplicate
       angular solver calls)

Reference: design_inter_fiber_coupling.md, Sections 4.5 and 9 Phase 1.
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import eval_legendre
from typing import Dict, Any, Optional, List, Tuple

from geovac.level4_multichannel import solve_angular_multichannel


def _channel_rotation_phases(
    channels: list,
    bond_angle: float = np.pi,
) -> np.ndarray:
    """
    Compute the rotation phase for each (l1, l2) channel at bond angle theta.

    Derivation (general bond angle):
    ---------------------------------
    In Level 4 mol-frame hyperspherical coordinates, both electrons are
    described relative to the bond axis. Electron 1 has angular momentum l1
    and electron 2 has angular momentum l2, both measured from the bond axis.

    When rotating from fiber A's frame to fiber B's frame by the bond angle
    theta, the Wigner d-matrix for sigma channels (m=0) gives:

        d^l_{00}(theta) = P_l(cos theta)

    Both electrons share the bond-axis coordinate system, so BOTH transform
    under the rotation. The transformation matrix is diagonal in the (l1, l2)
    basis:

        <l'1, l'2 | R_theta | l1, l2> = P_{l1}(cos theta) * P_{l2}(cos theta)
                                          * delta_{l'1,l1} * delta_{l'2,l2}

    For linear molecules (theta = pi):
        P_l(cos pi) = P_l(-1) = (-1)^l
        => phase = (-1)^{l1} * (-1)^{l2} = (-1)^{l1+l2}
    which recovers the existing formula.

    For theta = 0 (aligned fibers): P_l(1) = 1, so phase = 1 for all channels.

    Parameters
    ----------
    channels : list of (l1, l2) tuples
        Channel labels from solve_angular_multichannel.
    bond_angle : float
        Angle between bond axes (radians). Default pi (linear molecule).

    Returns
    -------
    phases : ndarray of shape (n_ch,)
        P_{l1}(cos theta) * P_{l2}(cos theta) for each channel.
    """
    cos_theta = np.cos(bond_angle)
    return np.array([
        float(eval_legendre(ch[0], cos_theta) * eval_legendre(ch[1], cos_theta))
        for ch in channels
    ], dtype=float)


def extract_origin_density(
    level4_result: Dict[str, Any],
    R: float,
    Z_A: float,
    Z_B: float,
    l_max: int,
    n_alpha: int,
    pk_potentials: Optional[List[dict]] = None,
    n_r: int = 300,
    r_max: float = 10.0,
    n_sample_Re: int = 10,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract one-electron radial density P(r) at distance r from the
    Level 4 coordinate origin.

    The Level 4 solver does not directly expose the angular channel
    eigenvectors (they are computed during the adiabatic sweep but not
    stored). This function re-runs the angular solver at a set of
    representative R_e values to reconstruct the angular density profile,
    then integrates over the (R_e, alpha) grid to obtain P(r).

    Parameters
    ----------
    level4_result : dict
        Output of solve_level4_h2_multichannel().
    R : float
        Internuclear distance (bohr).
    Z_A, Z_B : float
        Nuclear charges used in the Level 4 solve.
    l_max : int
        Angular momentum limit.
    n_alpha : int
        Number of alpha grid points.
    pk_potentials : list or None
        Phillips-Kleinman pseudopotential specs.
    n_r : int
        Number of output radial grid points.
    r_max : float
        Maximum distance on output grid (bohr).
    n_sample_Re : int
        Number of R_e values at which to sample the angular eigenvectors.

    Returns
    -------
    r_grid : ndarray, shape (n_r,)
        Radial grid (bohr).
    P_r : ndarray, shape (n_r,)
        One-electron radial density, normalized so integral P(r) dr = 2.
    """
    F = level4_result['wavefunction']
    R_e_grid = level4_result['R_e_grid_radial']
    z0 = level4_result.get('z0', 0.0)

    # Find region where |F|^2 is significant (>0.1% of peak)
    F2 = F ** 2
    threshold = 0.001 * F2.max()
    significant = np.where(F2 > threshold)[0]
    if len(significant) < 2:
        i_peak = np.argmax(F2)
        significant = np.arange(
            max(0, i_peak - 50), min(len(F2), i_peak + 50))

    Re_lo = R_e_grid[significant[0]]
    Re_hi = R_e_grid[significant[-1]]

    # Sample angular eigenvectors at n_sample_Re points in the significant range
    Re_samples = np.linspace(Re_lo, Re_hi, n_sample_Re)
    h_alpha = (np.pi / 2) / (n_alpha + 1)
    alpha_grid = (np.arange(n_alpha) + 1) * h_alpha

    # Compute angular density |u(alpha)|^2 at each sample R_e
    u_sq_samples = np.zeros((n_sample_Re, n_alpha))
    for k, Re_s in enumerate(Re_samples):
        rho = R / (2.0 * Re_s)
        if rho < 1e-8:
            continue
        _, vecs, _, channels = solve_angular_multichannel(
            rho, Re_s, l_max, n_alpha=n_alpha,
            Z_A=Z_A, Z_B=Z_B, z0=z0,
            pk_potentials=pk_potentials,
        )
        vec = vecs[0]  # ground state eigenvector, shape (n_ch * n_alpha,)
        n_ch = len(channels)
        # Sum |amplitude|^2 over channels at each alpha point
        vec_2d = vec.reshape(n_ch, n_alpha)
        u_sq_samples[k, :] = np.sum(vec_2d ** 2, axis=0)

    # Interpolate F to sample R_e points
    F_spline = CubicSpline(R_e_grid, F, extrapolate=True)

    # Output grid
    dr = r_max / n_r
    r_grid = (np.arange(n_r) + 0.5) * dr
    P_r = np.zeros(n_r)

    # R_e sample spacing
    dRe = Re_samples[1] - Re_samples[0] if n_sample_Re > 1 else 1.0

    # Integrate over (R_e, alpha) grid to build P(r)
    # Electron 1 at r1 = R_e cos(alpha)
    # Electron 2 at r2 = R_e sin(alpha)
    for k in range(n_sample_Re):
        Re_s = Re_samples[k]
        F_val = F_spline(Re_s)
        F2_val = F_val ** 2

        for j in range(n_alpha):
            cos_a = np.cos(alpha_grid[j])
            sin_a = np.sin(alpha_grid[j])

            r1 = Re_s * cos_a
            r2 = Re_s * sin_a

            w = F2_val * u_sq_samples[k, j] * dRe * h_alpha

            # Bin electron 1 at distance r1 from origin
            i1 = int(r1 / dr)
            if 0 <= i1 < n_r:
                P_r[i1] += w

            # Bin electron 2 at distance r2 from origin
            i2 = int(r2 / dr)
            if 0 <= i2 < n_r:
                P_r[i2] += w

    # Normalize to 2 electrons
    total = np.sum(P_r) * dr
    if total > 1e-15:
        P_r *= 2.0 / total

    return r_grid, P_r


def extract_channel_data(
    level4_result: Dict[str, Any],
    R: float,
    Z_A: float,
    Z_B: float,
    l_max: int,
    n_alpha: int,
    pk_potentials: Optional[List[dict]] = None,
    n_sample_Re: int = 15,
) -> Dict[str, Any]:
    """
    Extract angular channel data from a Level 4 result (single solver pass).

    Runs the angular solver at sample R_e points within the significant
    region of |F(R_e)|^2 and returns structured channel coefficient data
    usable by both the algebraic density extraction and the overlap diagnostic,
    eliminating duplicate solver calls.

    Parameters
    ----------
    level4_result : dict
        Output of solve_level4_h2_multichannel().
    R : float
        Internuclear distance (bohr).
    Z_A, Z_B : float
        Nuclear charges.
    l_max : int
        Angular momentum limit.
    n_alpha : int
        Number of alpha grid points.
    pk_potentials : list or None
        Phillips-Kleinman pseudopotential specs.
    n_sample_Re : int
        Number of R_e values at which to sample angular eigenvectors.

    Returns
    -------
    data : dict
        Re_samples : ndarray (n_sample,)
        F_values : ndarray (n_sample,) — F(R_e) interpolated to samples
        alpha_grid : ndarray (n_alpha,)
        h_alpha : float — alpha grid spacing
        channels : list of (l1, l2) tuples
        n_ch : int — number of channels
        vec_2d_list : list of ndarray (n_ch, n_alpha) — eigenvectors
        ch_weights : ndarray (n_sample, n_ch) — |c_{l1l2}|^2 at each R_e
        ang_density : ndarray (n_sample, n_alpha) — sum_ch |u_ch(alpha)|^2
    """
    F = level4_result['wavefunction']
    R_e_grid = level4_result['R_e_grid_radial']
    z0 = level4_result.get('z0', 0.0)

    # Find region where |F|^2 is significant (>0.1% of peak)
    F2 = F ** 2
    threshold = 0.001 * F2.max()
    significant = np.where(F2 > threshold)[0]
    if len(significant) < 2:
        i_peak = np.argmax(F2)
        significant = np.arange(
            max(0, i_peak - 50), min(len(F2), i_peak + 50))

    Re_lo = R_e_grid[significant[0]]
    Re_hi = R_e_grid[significant[-1]]
    Re_samples = np.linspace(Re_lo, Re_hi, n_sample_Re)

    h_alpha = (np.pi / 2) / (n_alpha + 1)
    alpha_grid = (np.arange(n_alpha) + 1) * h_alpha

    F_spline = CubicSpline(R_e_grid, F, extrapolate=True)
    F_values = F_spline(Re_samples)

    channels_out = None
    n_ch = 0
    vec_2d_list = []
    ch_weights_all = []
    ang_density_all = []

    for k, Re_s in enumerate(Re_samples):
        rho = R / (2.0 * Re_s)
        if rho < 1e-8:
            if channels_out is not None:
                vec_2d_list.append(np.zeros((n_ch, n_alpha)))
                ch_weights_all.append(np.zeros(n_ch))
            else:
                vec_2d_list.append(None)
                ch_weights_all.append(None)
            ang_density_all.append(np.zeros(n_alpha))
            continue

        _, vecs, _, channels = solve_angular_multichannel(
            rho, Re_s, l_max, n_alpha=n_alpha,
            Z_A=Z_A, Z_B=Z_B, z0=z0,
            pk_potentials=pk_potentials,
        )

        if channels_out is None:
            channels_out = channels
            n_ch = len(channels)
            # Back-fill any skipped entries
            for idx in range(len(vec_2d_list)):
                if vec_2d_list[idx] is None:
                    vec_2d_list[idx] = np.zeros((n_ch, n_alpha))
                    ch_weights_all[idx] = np.zeros(n_ch)

        vec = vecs[0]  # ground state eigenvector
        v2d = vec.reshape(n_ch, n_alpha)
        vec_2d_list.append(v2d)

        # Channel weights: |c_{l1,l2}|^2 = sum_j |u_ch(alpha_j)|^2
        cw = np.sum(v2d ** 2, axis=1)
        ch_weights_all.append(cw)

        # Angular density at each alpha: sum_ch |u_ch(alpha)|^2
        ang_d = np.sum(v2d ** 2, axis=0)
        ang_density_all.append(ang_d)

    ch_weights_arr = np.array(ch_weights_all)  # (n_sample, n_ch)
    ang_density_arr = np.array(ang_density_all)  # (n_sample, n_alpha)

    # Compute normalization factor for density integration.
    # The quadrature P(r) = 2 * sum_k |F(Re_k)|^2 * rho_ang(alpha(r,Re)) * J * dRe
    # needs a global normalization so that integral P(r) dr = 2.
    # Pre-compute the raw integral so channel densities can be normalized consistently.
    dRe = Re_samples[1] - Re_samples[0] if len(Re_samples) > 1 else 1.0
    raw_norm = 0.0
    for k in range(len(Re_samples)):
        F2 = F_values[k] ** 2
        raw_norm += F2 * np.sum(ang_density_all[k]) * h_alpha * dRe
    # raw_norm approximates the total probability (should be ~1 for normalized F, Phi)

    return {
        'Re_samples': Re_samples,
        'F_values': F_values,
        'alpha_grid': alpha_grid,
        'h_alpha': h_alpha,
        'channels': [tuple(c) for c in channels_out] if channels_out else [],
        'n_ch': n_ch,
        'vec_2d_list': vec_2d_list,
        'ch_weights': ch_weights_arr,
        'ang_density': ang_density_arr,
        'dRe': dRe,
        'raw_norm': raw_norm,
    }


def extract_origin_density_algebraic(
    channel_data: Dict[str, Any],
    n_r: int = 300,
    r_max: float = 10.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Algebraic one-electron radial density from channel data.

    Replaces the histogram-binning approach with proper quadrature.
    The one-electron marginal density at distance r from the fiber origin:

        P(r) = 2 * integral_{r}^{inf} |F(R_e)|^2
               * rho_ang(arccos(r/R_e)) / sqrt(R_e^2 - r^2) * dR_e

    where rho_ang(alpha) = sum_ch |u_ch(alpha)|^2 is the Liouville-weighted
    angular density (interpolated via cubic spline), and the factor 2 counts
    both electrons symmetrically.

    Parameters
    ----------
    channel_data : dict
        Output of extract_channel_data().
    n_r : int
        Number of output radial grid points.
    r_max : float
        Maximum distance on output grid (bohr).

    Returns
    -------
    r_grid : ndarray (n_r,)
    P_r : ndarray (n_r,), normalized so integral P(r) dr = 2.
    """
    Re_samples = channel_data['Re_samples']
    F_values = channel_data['F_values']
    alpha_grid = channel_data['alpha_grid']
    ang_density = channel_data['ang_density']  # (n_sample, n_alpha)

    dr = r_max / n_r
    r_grid = (np.arange(n_r) + 0.5) * dr
    P_r = np.zeros(n_r)

    n_sample = len(Re_samples)
    dRe = Re_samples[1] - Re_samples[0] if n_sample > 1 else 1.0

    # Build angular density spline at each R_e sample
    ang_splines = []
    for k in range(n_sample):
        sp = CubicSpline(alpha_grid, ang_density[k], extrapolate=False)
        ang_splines.append(sp)

    # Quadrature: for each r, integrate over R_e > r
    for k in range(n_sample):
        Re_s = Re_samples[k]
        F2 = F_values[k] ** 2
        if F2 < 1e-30:
            continue

        for i_r in range(n_r):
            r = r_grid[i_r]
            if r < 1e-12 or Re_s <= r * 1.001:
                continue

            ratio = r / Re_s
            if ratio > 1.0:
                continue
            alpha_at_r = np.arccos(ratio)

            # Evaluate angular density at this alpha via spline
            rho_ang = ang_splines[k](alpha_at_r)
            if rho_ang is None or np.isnan(rho_ang):
                rho_ang = 0.0
            rho_ang = max(rho_ang, 0.0)

            jacobian = 1.0 / np.sqrt(Re_s**2 - r**2)
            P_r[i_r] += 2.0 * F2 * rho_ang * jacobian * dRe

    # Normalize to 2 electrons
    total = np.sum(P_r) * dr
    if total > 1e-15:
        P_r *= 2.0 / total

    return r_grid, P_r


def extract_channel_densities(
    channel_data: Dict[str, Any],
    n_r: int = 300,
    r_max: float = 10.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Per-channel one-electron radial densities from channel data.

    Decomposes the total density P(r) = sum_ch P_ch(r) where each
    P_ch uses only the angular amplitude |u_ch(alpha)|^2 from that channel.
    Normalized so that sum_ch integral P_ch(r) dr = 2 (two electrons).

    Parameters
    ----------
    channel_data : dict
        Output of extract_channel_data().
    n_r : int
        Output grid size.
    r_max : float
        Maximum distance (bohr).

    Returns
    -------
    r_grid : ndarray (n_r,)
    P_ch : ndarray (n_ch, n_r) — per-channel densities, normalized so
        that sum over all channels integrates to 2.
    """
    Re_samples = channel_data['Re_samples']
    F_values = channel_data['F_values']
    alpha_grid = channel_data['alpha_grid']
    vec_2d_list = channel_data['vec_2d_list']
    n_ch = channel_data['n_ch']

    dr = r_max / n_r
    r_grid = (np.arange(n_r) + 0.5) * dr
    P_ch = np.zeros((n_ch, n_r))

    n_sample = len(Re_samples)
    dRe = Re_samples[1] - Re_samples[0] if n_sample > 1 else 1.0

    # Build per-channel angular density splines
    ch_splines = []
    for k in range(n_sample):
        v2d = vec_2d_list[k]
        splines_k = []
        for ic in range(n_ch):
            ch_ang = v2d[ic, :] ** 2  # |u_ch(alpha)|^2
            sp = CubicSpline(alpha_grid, ch_ang, extrapolate=False)
            splines_k.append(sp)
        ch_splines.append(splines_k)

    for k in range(n_sample):
        Re_s = Re_samples[k]
        F2 = F_values[k] ** 2
        if F2 < 1e-30:
            continue

        for i_r in range(n_r):
            r = r_grid[i_r]
            if r < 1e-12 or Re_s <= r * 1.001:
                continue

            ratio = r / Re_s
            if ratio > 1.0:
                continue
            alpha_at_r = np.arccos(ratio)
            jacobian = 1.0 / np.sqrt(Re_s**2 - r**2)

            for ic in range(n_ch):
                rho_ch = ch_splines[k][ic](alpha_at_r)
                if rho_ch is None or np.isnan(rho_ch):
                    rho_ch = 0.0
                rho_ch = max(rho_ch, 0.0)
                P_ch[ic, i_r] += 2.0 * F2 * rho_ch * jacobian * dRe

    # Normalize: sum_ch integral P_ch(r) dr = 2 (two electrons)
    total = np.sum(P_ch) * dr
    if total > 1e-15:
        P_ch *= 2.0 / total

    return r_grid, P_ch


def compute_channel_f0_matrix(
    channel_data: Dict[str, Any],
    R: float,
    n_r: int = 300,
    r_max: float = 10.0,
) -> Dict[str, Any]:
    """
    Decompose F^0 into a channel-channel matrix M(ch, ch').

    The total F^0 = sum_{ch, ch'} M(ch, ch') where each entry is the
    Slater F^0 integral between the Be-centered per-channel densities:

        M(ch, ch') = integral integral P_ch_Be(d1) P_ch'_Be(d2)
                     / max(d1, d2) dd1 dd2

    Parameters
    ----------
    channel_data : dict
        Output of extract_channel_data().
    R : float
        Internuclear distance (bohr).
    n_r : int
        Radial grid size.
    r_max : float
        Maximum distance (bohr).

    Returns
    -------
    result : dict
        F0_total : float — total F^0 = sum M(ch,ch')
        F0_matrix : ndarray (n_ch, n_ch) — channel decomposition
        F0_per_channel : ndarray (n_ch,) — row sums (each channel's
            contribution to total F^0)
        channels : list of (l1, l2)
        P_ch_Be : list of ndarray — per-channel Be-centered densities
        d_grid : ndarray — Be-centered grid
    """
    # Step 1: Per-channel origin densities
    r_grid, P_ch_origin = extract_channel_densities(
        channel_data, n_r=n_r, r_max=r_max)

    n_ch = channel_data['n_ch']
    channels = channel_data['channels']

    # Step 2: Transform each channel density to Be center
    P_ch_Be_list = []
    d_grid_ref = None
    for ic in range(n_ch):
        d_grid, P_Be_ic = transform_to_center_density(
            r_grid, P_ch_origin[ic], R)
        if d_grid_ref is None:
            d_grid_ref = d_grid
        # Pad or truncate to common grid
        n_d = len(d_grid_ref)
        if len(P_Be_ic) < n_d:
            P_Be_ic = np.pad(P_Be_ic, (0, n_d - len(P_Be_ic)))
        elif len(P_Be_ic) > n_d:
            P_Be_ic = P_Be_ic[:n_d]
        P_ch_Be_list.append(P_Be_ic)

    # Step 3: F^0 matrix
    F0_matrix = np.zeros((n_ch, n_ch))
    for ic in range(n_ch):
        for jc in range(ic, n_ch):
            f0_ij = slater_f0_integral(d_grid_ref, P_ch_Be_list[ic],
                                       P_ch_Be_list[jc])
            F0_matrix[ic, jc] = f0_ij
            F0_matrix[jc, ic] = f0_ij

    F0_total = np.sum(F0_matrix)
    F0_per_channel = np.sum(F0_matrix, axis=1)

    return {
        'F0_total': float(F0_total),
        'F0_matrix': F0_matrix,
        'F0_per_channel': F0_per_channel,
        'channels': channels,
        'P_ch_Be': P_ch_Be_list,
        'd_grid': d_grid_ref,
    }


def compute_channel_fk_matrix(
    channel_data: Dict[str, Any],
    R: float,
    k: int = 0,
    n_r: int = 300,
    r_max: float = 10.0,
) -> Dict[str, Any]:
    """
    Decompose F^k into a channel-channel matrix M_k(ch, ch').

    Generalizes compute_channel_f0_matrix to arbitrary multipole order k.
    The total F^k = sum_{ch, ch'} M_k(ch, ch') where each entry is the
    Slater F^k integral between the Be-centered per-channel densities:

        M_k(ch, ch') = integral integral P_ch_Be(d1) P_ch'_Be(d2)
                       * d_<^k / d_>^{k+1} dd1 dd2

    Parameters
    ----------
    channel_data : dict
        Output of extract_channel_data().
    R : float
        Internuclear distance (bohr).
    k : int
        Multipole order (0=monopole, 1=dipole, 2=quadrupole, ...).
    n_r : int
        Radial grid size.
    r_max : float
        Maximum distance (bohr).

    Returns
    -------
    result : dict
        Fk_total : float — total F^k = sum M_k(ch,ch')
        Fk_matrix : ndarray (n_ch, n_ch) — channel decomposition
        Fk_per_channel : ndarray (n_ch,) — row sums
        channels : list of (l1, l2)
        k : int — multipole order
    """
    if k == 0:
        # Delegate to optimized F^0 path
        f0_result = compute_channel_f0_matrix(
            channel_data, R, n_r=n_r, r_max=r_max)
        return {
            'Fk_total': f0_result['F0_total'],
            'Fk_matrix': f0_result['F0_matrix'],
            'Fk_per_channel': f0_result['F0_per_channel'],
            'channels': f0_result['channels'],
            'k': 0,
            'P_ch_Be': f0_result['P_ch_Be'],
            'd_grid': f0_result['d_grid'],
        }

    # Step 1: Per-channel origin densities
    r_grid, P_ch_origin = extract_channel_densities(
        channel_data, n_r=n_r, r_max=r_max)

    n_ch = channel_data['n_ch']
    channels = channel_data['channels']

    # Step 2: Transform each channel density to Be center
    P_ch_Be_list = []
    d_grid_ref = None
    for ic in range(n_ch):
        d_grid, P_Be_ic = transform_to_center_density(
            r_grid, P_ch_origin[ic], R)
        if d_grid_ref is None:
            d_grid_ref = d_grid
        n_d = len(d_grid_ref)
        if len(P_Be_ic) < n_d:
            P_Be_ic = np.pad(P_Be_ic, (0, n_d - len(P_Be_ic)))
        elif len(P_Be_ic) > n_d:
            P_Be_ic = P_Be_ic[:n_d]
        P_ch_Be_list.append(P_Be_ic)

    # Step 3: F^k matrix
    Fk_matrix = np.zeros((n_ch, n_ch))
    for ic in range(n_ch):
        for jc in range(ic, n_ch):
            fk_ij = slater_fk_integral(d_grid_ref, P_ch_Be_list[ic],
                                        P_ch_Be_list[jc], k=k)
            Fk_matrix[ic, jc] = fk_ij
            Fk_matrix[jc, ic] = fk_ij

    Fk_total = np.sum(Fk_matrix)
    Fk_per_channel = np.sum(Fk_matrix, axis=1)

    return {
        'Fk_total': float(Fk_total),
        'Fk_matrix': Fk_matrix,
        'Fk_per_channel': Fk_per_channel,
        'channels': channels,
        'k': k,
        'P_ch_Be': P_ch_Be_list,
        'd_grid': d_grid_ref,
    }


def direct_exchange_inter_fiber_energy(
    level4_result: Dict[str, Any],
    R: float,
    Z_A: float,
    Z_B: float,
    l_max: int,
    n_alpha: int,
    pk_potentials: Optional[List[dict]] = None,
    n_r: int = 300,
    r_max: float = 10.0,
    n_sample_Re: int = 10,
    bond_angle: float = np.pi,
) -> Dict[str, Any]:
    """
    Compute direct exchange energy: E = -sum_ch phase_ch * F^0_ch.

    The rotation phase P_{l1}(cos θ) * P_{l2}(cos θ) generalizes the
    linear (-1)^{l1+l2} to arbitrary bond angles. Each channel's F^0
    contribution is weighted by its own rotation phase.

    Parameters
    ----------
    level4_result : dict
        Output of solve_level4_h2_multichannel().
    R, Z_A, Z_B, l_max, n_alpha, pk_potentials, n_r, r_max, n_sample_Re :
        Standard inter-fiber coupling parameters.
    bond_angle : float
        Angle between bond axes (radians). Default π (linear).

    Returns
    -------
    result : dict
        E_exchange : float — direct exchange energy (negative = attractive)
        F0_per_channel : ndarray — per-channel F^0 contributions
        channels : list of (l1, l2) tuples
        parity : ndarray — rotation phases for each channel
        F0_total : float — total (unsigned) F^0
    """
    channel_data = extract_channel_data(
        level4_result, R, Z_A, Z_B, l_max, n_alpha,
        pk_potentials=pk_potentials,
        n_sample_Re=n_sample_Re,
    )

    # Channel F^0 matrix
    f0_result = compute_channel_f0_matrix(
        channel_data, R, n_r=n_r, r_max=r_max)

    channels = f0_result['channels']
    F0_per_ch = f0_result['F0_per_channel']

    # Rotation phase: P_{l1}(cos θ) * P_{l2}(cos θ)
    phases = _channel_rotation_phases(channels, bond_angle)

    # Direct exchange: E = -sum_ch phase_ch * F0_ch
    F0_signed = float(np.dot(phases, F0_per_ch))
    E_exchange = -F0_signed

    return {
        'E_exchange': float(E_exchange),
        'F0_per_channel': F0_per_ch,
        'channels': channels,
        'parity': phases,
        'F0_total': f0_result['F0_total'],
        'F0_signed': F0_signed,
        'channel_data': channel_data,
    }


def multipole_exchange_inter_fiber_energy(
    level4_result: Dict[str, Any],
    R: float,
    Z_A: float,
    Z_B: float,
    l_max: int,
    n_alpha: int,
    pk_potentials: Optional[List[dict]] = None,
    n_r: int = 300,
    r_max: float = 10.0,
    n_sample_Re: int = 10,
    k_max: int = 2,
    bond_angle: float = np.pi,
) -> Dict[str, Any]:
    """
    Compute multipole exchange energy summing E_k from k=0 to k_max.

    At each multipole order k, the exchange energy is:

        E_k = -sum_{ch,ch'} phase_ch * G(l1, k, l1') * F^k(ch, ch')

    where phase_ch = P_{l1}(cos θ) * P_{l2}(cos θ) is the rotation phase,
    G(l, k, l') is the Gaunt integral, and F^k(ch, ch') is the Slater
    multipole integral. For θ=π: phase = (-1)^{l1+l2}.

    Parameters
    ----------
    level4_result : dict
        Output of solve_level4_h2_multichannel().
    R, Z_A, Z_B, l_max, n_alpha, pk_potentials, n_r, r_max, n_sample_Re :
        Standard inter-fiber coupling parameters.
    k_max : int
        Maximum multipole order (0=monopole only, 1=+dipole, 2=+quadrupole).
    bond_angle : float
        Angle between bond axes (radians). Default π (linear).

    Returns
    -------
    result : dict
        E_exchange : float — total multipole exchange energy
        E_per_k : dict — {k: E_k} for each multipole order
        Fk_per_k : dict — {k: Fk_total} for each multipole order
        channels : list of (l1, l2) tuples
        k_max : int
    """
    from geovac.hyperspherical_angular import gaunt_integral

    channel_data = extract_channel_data(
        level4_result, R, Z_A, Z_B, l_max, n_alpha,
        pk_potentials=pk_potentials,
        n_sample_Re=n_sample_Re,
    )

    channels = channel_data['channels']
    n_ch = channel_data['n_ch']

    # Rotation phase: P_{l1}(cos θ) * P_{l2}(cos θ)
    parity = _channel_rotation_phases(channels, bond_angle)

    E_total = 0.0
    E_per_k = {}
    Fk_per_k = {}
    Fk_matrices = {}

    for k in range(k_max + 1):
        # Compute F^k channel matrix
        fk_result = compute_channel_fk_matrix(
            channel_data, R, k=k, n_r=n_r, r_max=r_max)
        Fk_matrix = fk_result['Fk_matrix']
        Fk_per_k[k] = fk_result['Fk_total']
        Fk_matrices[k] = Fk_matrix

        if k == 0:
            # k=0 (monopole): no angular selection between channels.
            # All channel cross-terms contribute. This recovers the
            # direct exchange formula: E_0 = -sum_ch parity_ch * F0_ch.
            Fk_per_ch = np.sum(Fk_matrix, axis=1)
            E_k = -float(np.dot(parity, Fk_per_ch))
        else:
            # k >= 1: Gaunt angular coupling restricts which channels
            # interact. Selection rules: |l-l'| <= k <= l+l', l+k+l' even.
            #
            # The angular coupling factor G_factor(i,j,k) is the product
            # of Gaunt integrals for both electrons (Be-side and H-side),
            # normalized by the k=0 diagonal values so that the k>0
            # contributions are properly scaled relative to the monopole.
            E_k = 0.0
            for i in range(n_ch):
                l1_i, l2_i = channels[i]
                g0_l1i = gaunt_integral(l1_i, 0, l1_i)
                g0_l2i = gaunt_integral(l2_i, 0, l2_i)
                for j in range(n_ch):
                    l1_j, l2_j = channels[j]
                    G1 = gaunt_integral(l1_i, k, l1_j)
                    G2 = gaunt_integral(l2_i, k, l2_j)
                    if abs(G1) < 1e-15 or abs(G2) < 1e-15:
                        continue
                    g0_l1j = gaunt_integral(l1_j, 0, l1_j)
                    g0_l2j = gaunt_integral(l2_j, 0, l2_j)
                    norm_i = g0_l1i * g0_l2i
                    norm_j = g0_l1j * g0_l2j
                    norm = np.sqrt(abs(norm_i * norm_j))
                    if norm < 1e-15:
                        continue
                    G_factor = G1 * G2 / norm
                    E_k -= parity[i] * G_factor * Fk_matrix[i, j]

        E_per_k[k] = float(E_k)
        E_total += E_k

    return {
        'E_exchange': float(E_total),
        'E_per_k': E_per_k,
        'Fk_per_k': Fk_per_k,
        'Fk_matrices': Fk_matrices,
        'channels': channels,
        'parity': parity.tolist(),
        'k_max': k_max,
        'channel_data': channel_data,
    }


def compute_overlap_from_channel_data(
    channel_data: Dict[str, Any],
    bond_angle: float = np.pi,
) -> Dict[str, Any]:
    """
    Compute inter-fiber channel overlap S(R) from pre-extracted channel data.

    Identical physics to compute_overlap_diagnostic() but avoids
    re-running the angular solver (uses cached eigenvectors).

    S(rho, theta) = sum_{l1, l2} P_{l1}(cos theta) * P_{l2}(cos theta)
                    * |c^0_{l1 l2}(rho)|^2

    For theta=pi (linear): recovers (-1)^{l1+l2} formula.

    Parameters
    ----------
    channel_data : dict
        Output of extract_channel_data().
    bond_angle : float
        Angle between bond axes (radians). Default π (linear).

    Returns
    -------
    result : dict with keys S_avg, S_samples, Re_samples, channel_weights,
        channels, F0_signed_weights (parity-weighted channel contributions).
    """
    Re_samples = channel_data['Re_samples']
    F_values = channel_data['F_values']
    ch_weights = channel_data['ch_weights']  # (n_sample, n_ch)
    channels = channel_data['channels']
    n_ch = channel_data['n_ch']

    n_sample = len(Re_samples)
    F2_weights = F_values ** 2

    # Rotation phase: P_{l1}(cos θ) * P_{l2}(cos θ) for each channel
    phases = _channel_rotation_phases(channels, bond_angle)

    # S(rho, theta) at each sample R_e
    S_samples = np.dot(ch_weights, phases)  # (n_sample,)

    # |F|^2-weighted average
    total_weight = np.sum(F2_weights)
    if total_weight > 1e-15:
        S_avg = np.sum(S_samples * F2_weights) / total_weight
        avg_ch_weights = np.sum(
            ch_weights * F2_weights[:, np.newaxis], axis=0) / total_weight
    else:
        S_avg = 0.0
        avg_ch_weights = np.zeros(n_ch)

    ch_weight_dict = {}
    for i, ch in enumerate(channels):
        ch_weight_dict[tuple(ch)] = float(avg_ch_weights[i])

    return {
        'S_avg': float(S_avg),
        'S_samples': S_samples.tolist(),
        'Re_samples': Re_samples.tolist(),
        'F2_weights': F2_weights.tolist(),
        'channel_weights': ch_weight_dict,
        'channels': channels,
        'parity': phases,
        'avg_ch_weights': avg_ch_weights,
    }


def transform_to_center_density(
    r_grid: np.ndarray,
    P_origin: np.ndarray,
    R: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Transform density from fiber origin to shared center coordinates.

    For an s-wave electron at distance r from the fiber origin (located at
    distance s = R/2 from the shared center), the PDF of its distance d
    from the shared center is (shell theorem):

        p(d | r, s) = d / (2 r s)   for  |r - s| <= d <= r + s

    This preserves normalization: integral P_center(d) dd = integral P_origin(r) dr.

    Parameters
    ----------
    r_grid : ndarray
        Input radial grid (fiber-origin centered).
    P_origin : ndarray
        Radial density at distance from fiber origin.
    R : float
        Internuclear distance (bohr).

    Returns
    -------
    d_grid : ndarray
        Output grid (distance from shared center).
    P_center : ndarray
        Transformed density.
    """
    s = R / 2.0
    dr = r_grid[1] - r_grid[0]

    # Output grid extends to max(r) + s
    d_upper = r_grid[-1] + s + dr
    n_d = int(d_upper / dr) + 1
    d_grid = (np.arange(n_d) + 0.5) * dr
    P_center = np.zeros(n_d)

    for j in range(len(r_grid)):
        r = r_grid[j]
        if r < 1e-12 or P_origin[j] < 1e-30:
            continue

        d_min = abs(r - s)
        d_max_j = r + s

        # Output bins within the allowed range
        i_lo = max(0, int(d_min / dr))
        i_hi = min(n_d - 1, int(d_max_j / dr))

        for i in range(i_lo, i_hi + 1):
            d = d_grid[i]
            if d < d_min or d > d_max_j or d < 1e-12:
                continue
            kernel = d / (2.0 * r * s)
            P_center[i] += P_origin[j] * kernel * dr

    return d_grid, P_center


def slater_fk_integral(
    r_grid: np.ndarray,
    P_A: np.ndarray,
    P_B: np.ndarray,
    k: int = 0,
) -> float:
    """
    Compute Slater F^k integral between two radial densities.

        F^k = integral integral P_A(r1) P_B(r2) r_<^k / r_>^{k+1} dr1 dr2

    For k=0, this reduces to F^0 = integral integral P_A P_B / r_> dr1 dr2.

    Efficient O(N) computation using cumulative charge method:
        Q_k(r) = integral_0^r P_B(r') r'^k dr'
        T_k(r) = integral_r^inf P_B(r') / r'^{k+1} dr'
        phi_k(r) = Q_k(r) / r^{k+1} + T_k(r) * r^k
        F^k = integral P_A(r) phi_k(r) dr

    Parameters
    ----------
    r_grid : ndarray
        Radial grid (uniform spacing).
    P_A, P_B : ndarray
        Radial probability densities.
    k : int
        Multipole order (0=monopole, 1=dipole, 2=quadrupole, ...).

    Returns
    -------
    fk : float
        The F^k integral (Hartree).
    """
    dr = r_grid[1] - r_grid[0]

    if k == 0:
        # Optimized k=0 path (original slater_f0_integral)
        Q_B = np.cumsum(P_B) * dr
        P_over_r = np.where(r_grid > 1e-12, P_B / r_grid, 0.0)
        T_B = np.flip(np.cumsum(np.flip(P_over_r))) * dr
        phi_B = np.where(r_grid > 1e-12, Q_B / r_grid, 0.0) + T_B
        return float(np.sum(P_A * phi_B) * dr)

    # General k >= 1
    # Q_k(r) = integral_0^r P_B(r') r'^k dr'
    Q_k = np.cumsum(P_B * r_grid**k) * dr

    # T_k(r) = integral_r^inf P_B(r') / r'^{k+1} dr'
    P_over_rk1 = np.where(r_grid > 1e-12, P_B / r_grid**(k + 1), 0.0)
    T_k = np.flip(np.cumsum(np.flip(P_over_rk1))) * dr

    # phi_k(r) = Q_k(r) / r^{k+1} + T_k(r) * r^k
    phi_k = np.where(r_grid > 1e-12, Q_k / r_grid**(k + 1), 0.0) + T_k * r_grid**k

    return float(np.sum(P_A * phi_k) * dr)


def slater_f0_integral(
    r_grid: np.ndarray,
    P_A: np.ndarray,
    P_B: np.ndarray,
) -> float:
    """
    Compute Slater F^0 integral between two radial densities.

        F^0 = integral_0^inf integral_0^inf P_A(r1) P_B(r2) / max(r1, r2) dr1 dr2

    Convenience wrapper around slater_fk_integral(k=0).

    Parameters
    ----------
    r_grid : ndarray
        Radial grid (uniform spacing).
    P_A, P_B : ndarray
        Radial probability densities.

    Returns
    -------
    f0 : float
        The F^0 integral (Hartree).
    """
    return slater_fk_integral(r_grid, P_A, P_B, k=0)


def monopole_inter_fiber_energy(
    level4_result: Dict[str, Any],
    R: float,
    Z_A: float,
    Z_B: float,
    l_max: int,
    n_alpha: int,
    pk_potentials: Optional[List[dict]] = None,
    n_r: int = 300,
    r_max: float = 10.0,
    n_sample_Re: int = 10,
    method: str = 'numerical',
    channel_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Compute monopole (k=0) inter-fiber coupling energy from a Level 4 result.

    Pipeline:
    1. Extract one-electron density P(r) from the Level 4 wavefunction.
    2. Transform to shared-center (Be) coordinates via shell theorem.
    3. Compute Slater F^0 integral between the two fiber densities.

    For identical fibers (D_inf_h symmetry), both densities are the same.

    Parameters
    ----------
    level4_result : dict
        Output of solve_level4_h2_multichannel().
    R : float
        Internuclear distance (bohr).
    Z_A, Z_B : float
        Nuclear charges (Z_A = Z_eff of center atom, Z_B = ligand).
    l_max : int
        Angular momentum limit.
    n_alpha : int
        Alpha grid points.
    pk_potentials : list or None
        Phillips-Kleinman pseudopotential.
    method : str
        'numerical' (default): histogram-bin density extraction.
        'algebraic': quadrature-based density with channel decomposition.
    channel_data : dict or None
        Pre-computed channel data from extract_channel_data(). If None
        and method='algebraic', computed internally.

    Returns
    -------
    result : dict
        Keys: E_monopole, r_grid_origin, P_origin, d_grid, P_Be,
        N_elec_origin, N_elec_Be.
        If method='algebraic', also: channel_data, F0_matrix_result.
    """
    if method == 'algebraic':
        # Algebraic pathway: quadrature-based density + channel decomposition
        if channel_data is None:
            channel_data = extract_channel_data(
                level4_result, R, Z_A, Z_B, l_max, n_alpha,
                pk_potentials=pk_potentials,
                n_sample_Re=n_sample_Re,
            )

        # Step 1: Algebraic origin density (quadrature, no histogram)
        r_grid_origin, P_origin = extract_origin_density_algebraic(
            channel_data, n_r=n_r, r_max=r_max)

        # Step 2: Transform to Be-centered density
        d_grid, P_Be = transform_to_center_density(
            r_grid_origin, P_origin, R)

        # Step 3: F^0 integral
        E_mono = slater_f0_integral(d_grid, P_Be, P_Be)

        # Step 4: Channel decomposition of F^0
        f0_matrix_result = compute_channel_f0_matrix(
            channel_data, R, n_r=n_r, r_max=r_max)

        dr_origin = r_grid_origin[1] - r_grid_origin[0]
        dr_Be = d_grid[1] - d_grid[0]

        return {
            'E_monopole': E_mono,
            'r_grid_origin': r_grid_origin,
            'P_origin': P_origin,
            'd_grid': d_grid,
            'P_Be': P_Be,
            'N_elec_origin': float(np.sum(P_origin) * dr_origin),
            'N_elec_Be': float(np.sum(P_Be) * dr_Be),
            'channel_data': channel_data,
            'F0_matrix_result': f0_matrix_result,
        }

    # --- Numerical pathway (original) ---
    # Step 1: Extract origin-centered density
    r_grid_origin, P_origin = extract_origin_density(
        level4_result, R, Z_A, Z_B, l_max, n_alpha,
        pk_potentials=pk_potentials,
        n_r=n_r, r_max=r_max, n_sample_Re=n_sample_Re,
    )

    # Step 2: Transform to Be-centered density
    d_grid, P_Be = transform_to_center_density(r_grid_origin, P_origin, R)

    # Step 3: F^0 integral (identical fibers by D_inf_h symmetry)
    E_mono = slater_f0_integral(d_grid, P_Be, P_Be)

    dr_origin = r_grid_origin[1] - r_grid_origin[0]
    dr_Be = d_grid[1] - d_grid[0]

    return {
        'E_monopole': E_mono,
        'r_grid_origin': r_grid_origin,
        'P_origin': P_origin,
        'd_grid': d_grid,
        'P_Be': P_Be,
        'N_elec_origin': float(np.sum(P_origin) * dr_origin),
        'N_elec_Be': float(np.sum(P_Be) * dr_Be),
    }


def exchange_inter_fiber_energy(
    level4_result: Dict[str, Any],
    R: float,
    Z_A: float,
    Z_B: float,
    l_max: int,
    n_alpha: int,
    pk_potentials: Optional[List[dict]] = None,
    n_r: int = 300,
    r_max: float = 10.0,
    n_sample_Re: int = 10,
    S_avg: Optional[float] = None,
    method: str = 'numerical',
    bond_angle: float = np.pi,
) -> Dict[str, Any]:
    """
    Compute exchange inter-fiber coupling energy via E_exch = -S_avg * F^0.

    The exchange energy between two identical fibers is approximated as:

        E_exch(R) = -S_avg(R, θ) * F^0(R)

    where S_avg(R, θ) is the inter-fiber channel overlap at bond angle θ
    and F^0(R) is the monopole Coulomb integral.

    Parameters
    ----------
    level4_result : dict
        Output of solve_level4_h2_multichannel().
    R : float
        Internuclear distance (bohr).
    Z_A, Z_B : float
        Nuclear charges.
    l_max : int
        Angular momentum limit.
    n_alpha : int
        Alpha grid points.
    pk_potentials : list or None
        Phillips-Kleinman pseudopotential.
    n_r : int
        Radial grid points.
    r_max : float
        Maximum radius (bohr).
    n_sample_Re : int
        Number of R_e sample points.
    S_avg : float or None
        Pre-computed inter-fiber overlap. If None, computed internally.
    method : str
        'numerical' (default): original histogram-based pipeline.
        'algebraic': quadrature-based density + channel decomposition.
            Uses single extract_channel_data() call for both S and F^0.
    bond_angle : float
        Angle between bond axes (radians). Default π (linear).

    Returns
    -------
    result : dict
        Keys: E_exchange, S_avg, F0_monopole, E_monopole_raw.
        If method='algebraic', also: channel_data, F0_matrix,
        F0_signed (parity-weighted F^0), E_exchange_direct
        (channel-resolved exchange without S*F^0 factorization).
    """
    if method == 'algebraic':
        # Single channel data extraction for both S and F^0
        channel_data = extract_channel_data(
            level4_result, R, Z_A, Z_B, l_max, n_alpha,
            pk_potentials=pk_potentials,
            n_sample_Re=n_sample_Re,
        )

        # F^0 via algebraic monopole
        mono = monopole_inter_fiber_energy(
            level4_result, R, Z_A, Z_B, l_max, n_alpha,
            pk_potentials=pk_potentials,
            n_r=n_r, r_max=r_max, n_sample_Re=n_sample_Re,
            method='algebraic', channel_data=channel_data,
        )
        F0 = mono['E_monopole']
        f0_matrix_result = mono['F0_matrix_result']

        # S_avg from same channel data (no duplicate solver call)
        if S_avg is None:
            ovlp = compute_overlap_from_channel_data(
                channel_data, bond_angle=bond_angle)
            S_avg = ovlp['S_avg']
            phases = ovlp['parity']
            avg_ch_weights = ovlp['avg_ch_weights']
        else:
            channels = channel_data['channels']
            phases = _channel_rotation_phases(channels, bond_angle)
            avg_ch_weights = None

        # Standard exchange: E = -S * F^0
        E_exchange = -S_avg * F0

        # Channel-resolved exchange: E_direct = -sum_ch parity_ch * F0_ch
        # This avoids the S*F^0 factorization and weights each channel's
        # F^0 contribution by its own parity, potentially more accurate.
        F0_per_ch = f0_matrix_result['F0_per_channel']
        n_ch = channel_data['n_ch']
        if avg_ch_weights is not None and len(avg_ch_weights) == n_ch:
            # Weight F0_per_channel by rotation-phase-signed channel weights
            # F0_signed = sum_{ch} P_{l1}(cos θ) * P_{l2}(cos θ) * F0_ch
            F0_signed = float(np.dot(phases, F0_per_ch))
            E_exchange_direct = -F0_signed
        else:
            F0_signed = None
            E_exchange_direct = None

        return {
            'E_exchange': float(E_exchange),
            'S_avg': float(S_avg),
            'F0_monopole': float(F0),
            'E_monopole_raw': float(F0),
            'channel_data': channel_data,
            'F0_matrix': f0_matrix_result['F0_matrix'],
            'F0_per_channel': F0_per_ch,
            'F0_signed': F0_signed,
            'E_exchange_direct': E_exchange_direct,
        }

    # --- Numerical pathway (original) ---
    # Step 1: Compute monopole F^0 from full density at Be center
    mono = monopole_inter_fiber_energy(
        level4_result, R, Z_A, Z_B, l_max, n_alpha,
        pk_potentials=pk_potentials,
        n_r=n_r, r_max=r_max, n_sample_Re=n_sample_Re,
    )
    F0 = mono['E_monopole']

    # Step 2: Compute S_avg if not provided
    if S_avg is None:
        ovlp = compute_overlap_diagnostic(
            R=R, Z_A=Z_A, Z_B=Z_B,
            l_max=l_max, n_alpha=n_alpha,
            level4_result=level4_result,
            pk_potentials=pk_potentials,
            n_sample_Re=n_sample_Re,
            bond_angle=bond_angle,
        )
        S_avg = ovlp['S_avg']

    # Step 3: Exchange energy = -S * F^0
    E_exchange = -S_avg * F0

    return {
        'E_exchange': float(E_exchange),
        'S_avg': float(S_avg),
        'F0_monopole': float(F0),
        'E_monopole_raw': float(F0),
    }


def compute_overlap_diagnostic(
    R: float,
    Z_A: float,
    Z_B: float,
    l_max: int,
    n_alpha: int,
    level4_result: Dict[str, Any],
    pk_potentials: Optional[List[dict]] = None,
    n_sample_Re: int = 15,
    bond_angle: float = np.pi,
) -> Dict[str, Any]:
    """
    Compute inter-fiber channel overlap S(R) as a diagnostic for exchange.

    For identical fibers related by rotation of the bond axis by angle theta,
    the overlap between ground adiabatic channels of fibers A and B is:

        S(ρ, θ) = Σ_{l₁,l₂} P_{l₁}(cos θ) · P_{l₂}(cos θ) |c^0_{l₁l₂}(ρ)|²

    where P_l is the Legendre polynomial (Wigner d^l_{00}(θ) for σ channels).
    Both electrons transform because both are described relative to the bond
    axis in Level 4 coordinates.

    For θ = π (linear): P_l(-1) = (-1)^l, recovering S = Σ (-1)^{l₁+l₂} |c|².
    For θ = 0 (aligned): P_l(1) = 1, giving S = Σ |c|² = 1 (max overlap).

    The |F(R_e)|²-weighted average gives S_avg(R).

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    Z_A, Z_B : float
        Nuclear charges (Z_A = Z_eff center, Z_B = ligand).
    l_max : int
        Angular momentum limit.
    n_alpha : int
        Alpha grid points.
    level4_result : dict
        Output of solve_level4_h2_multichannel().
    pk_potentials : list or None
        Phillips-Kleinman pseudopotential.
    n_sample_Re : int
        Number of R_e sample points.
    bond_angle : float
        Angle between bond axes (radians). Default π (linear molecule).
        For H₂O: 1.824 rad (104.5°).

    Returns
    -------
    result : dict
        Keys:
        - S_avg: |F|²-weighted average overlap
        - S_samples: S(ρ) at each sample R_e
        - Re_samples: the sample R_e values
        - channel_weights: dict mapping (l1,l2) -> average |c|²
        - channels: list of (l1,l2) tuples
    """
    F = level4_result['wavefunction']
    R_e_grid = level4_result['R_e_grid_radial']
    z0 = level4_result.get('z0', 0.0)

    # Find significant region of |F|²
    F2 = F ** 2
    threshold = 0.001 * F2.max()
    significant = np.where(F2 > threshold)[0]
    if len(significant) < 2:
        i_peak = np.argmax(F2)
        significant = np.arange(
            max(0, i_peak - 50), min(len(F2), i_peak + 50))

    Re_lo = R_e_grid[significant[0]]
    Re_hi = R_e_grid[significant[-1]]
    Re_samples = np.linspace(Re_lo, Re_hi, n_sample_Re)

    # Interpolate F to sample points
    F_spline = CubicSpline(R_e_grid, F, extrapolate=True)

    S_samples = np.zeros(n_sample_Re)
    F2_weights = np.zeros(n_sample_Re)
    channels_out = None
    # Accumulate channel weights for reporting
    channel_weight_accum = None

    for k, Re_s in enumerate(Re_samples):
        rho = R / (2.0 * Re_s)
        if rho < 1e-8:
            continue

        _, vecs, _, channels = solve_angular_multichannel(
            rho, Re_s, l_max, n_alpha=n_alpha,
            Z_A=Z_A, Z_B=Z_B, z0=z0,
            pk_potentials=pk_potentials,
        )

        if channels_out is None:
            channels_out = channels
            n_ch = len(channels)
            channel_weight_accum = np.zeros(n_ch)

        vec = vecs[0]  # ground state, shape (n_ch * n_alpha,)
        n_ch = len(channels)
        vec_2d = vec.reshape(n_ch, n_alpha)

        # Channel weights: |c_{l1,l2}|² = sum_i |vec_2d[ch, i]|²
        # (eigenvector is unit-normalized, so these sum to 1)
        ch_weights = np.sum(vec_2d ** 2, axis=1)  # shape (n_ch,)

        # Rotation phase: P_{l1}(cos θ) * P_{l2}(cos θ) for each channel
        phases = _channel_rotation_phases(channels, bond_angle)

        # S(ρ, θ) = Σ P_{l1}(cos θ) · P_{l2}(cos θ) · |c_{l1,l2}|²
        S_samples[k] = np.sum(phases * ch_weights)

        F_val = F_spline(Re_s)
        F2_weights[k] = F_val ** 2

        channel_weight_accum += ch_weights * F_val ** 2

    # |F|²-weighted average
    total_weight = np.sum(F2_weights)
    if total_weight > 1e-15:
        S_avg = np.sum(S_samples * F2_weights) / total_weight
        channel_weight_accum /= total_weight
    else:
        S_avg = 0.0

    # Build channel weight dict
    ch_weight_dict = {}
    if channels_out is not None:
        for i, ch in enumerate(channels_out):
            ch_weight_dict[tuple(ch)] = float(channel_weight_accum[i])

    return {
        'S_avg': float(S_avg),
        'S_samples': S_samples.tolist(),
        'Re_samples': Re_samples.tolist(),
        'F2_weights': F2_weights.tolist(),
        'channel_weights': ch_weight_dict,
        'channels': [tuple(c) for c in channels_out] if channels_out else [],
    }


def extract_channel_1rdm(
    level4_result: Dict[str, Any],
    R: float,
    Z_A: float,
    Z_B: float,
    l_max: int,
    n_alpha: int,
    pk_potentials: Optional[List[dict]] = None,
    n_sample_Re: int = 15,
    channel_data: Optional[Dict[str, Any]] = None,
    bond_angle: float = np.pi,
) -> Dict[str, Any]:
    """
    Extract the one-particle reduced density matrix in the channel
    (partial-wave) representation from a Level 4 wavefunction.

    The 1-RDM is obtained by tracing out electron 2 from the two-electron
    wavefunction.  In the channel representation with channels (l1, l2):

        gamma_{l1, l1'}(R_e) = sum_{l2} sum_j v2d[(l1,l2), j] * v2d[(l1',l2), j]

    For general bond angle theta, the inter-fiber 1-RDM transform requires
    an l2-weighted version where each l2 contribution is weighted by
    P_{l2}(cos theta)^2.  When bond_angle != pi, both the standard gamma
    and the weighted gamma are returned. For bond_angle = pi, the weight
    is 1 for all l2, so gamma_weighted = gamma (backward compatible).

    Parameters
    ----------
    level4_result : dict
        Output of solve_level4_h2_multichannel().
    R : float
        Internuclear distance (bohr).
    Z_A, Z_B : float
        Nuclear charges.
    l_max : int
        Angular momentum limit.
    n_alpha : int
        Number of alpha grid points.
    pk_potentials : list or None
        Phillips-Kleinman pseudopotential specs.
    n_sample_Re : int
        Number of R_e sample points.
    channel_data : dict or None
        Pre-computed output of extract_channel_data().  If None, computed
        internally.
    bond_angle : float
        Angle between bond axes (radians). Default π (linear). Used to
        compute the l2-weighted gamma for the inter-fiber transform.

    Returns
    -------
    result : dict
        gamma_matrix : ndarray (n_l, n_l, n_sample) — 1-RDM at each R_e
        gamma_avg : ndarray (n_l, n_l) — |F|^2-weighted average 1-RDM
        gamma_weighted_avg : ndarray (n_l, n_l) — l2-weighted 1-RDM for
            inter-fiber transform at bond_angle (= gamma_avg when θ=π)
        l_values : ndarray — unique l1 values indexing the matrix
        Re_samples : ndarray — R_e sample points
        F_values : ndarray — F(R_e) at samples
        channels : list of (l1, l2) tuples
        ch_weights : ndarray (n_sample, n_ch) — per-channel weights
        offdiag_ratio : float — ||gamma_offdiag|| / ||gamma_diag|| for avg
    """
    if channel_data is None:
        channel_data = extract_channel_data(
            level4_result, R, Z_A, Z_B, l_max, n_alpha,
            pk_potentials=pk_potentials, n_sample_Re=n_sample_Re,
        )

    channels = channel_data['channels']
    n_ch = channel_data['n_ch']
    vec_2d_list = channel_data['vec_2d_list']
    Re_samples = channel_data['Re_samples']
    F_values = channel_data['F_values']
    ch_weights = channel_data['ch_weights']
    n_sample = len(Re_samples)

    # Build mapping: unique l1 values and channel grouping by l2
    l1_values_all = [ch[0] for ch in channels]
    l2_values_all = [ch[1] for ch in channels]
    l_values = np.array(sorted(set(l1_values_all)))
    n_l = len(l_values)
    l_to_idx = {l: i for i, l in enumerate(l_values)}

    # For each l2, collect list of (channel_index, l1) pairs
    l2_set = sorted(set(l2_values_all))
    l2_groups: Dict[int, List[Tuple[int, int]]] = {}
    for ic, (l1, l2) in enumerate(channels):
        if l2 not in l2_groups:
            l2_groups[l2] = []
        l2_groups[l2].append((ic, l1))

    # Precompute l2 rotation weights: P_{l2}(cos theta)^2
    # For theta=pi: P_{l2}(-1)^2 = 1 for all l2 (no weighting)
    cos_theta = np.cos(bond_angle)
    l2_weight = {l2: float(eval_legendre(l2, cos_theta)) ** 2
                 for l2 in l2_set}

    # Compute 1-RDM at each R_e sample (standard and l2-weighted)
    gamma_matrix = np.zeros((n_l, n_l, n_sample))
    gamma_weighted_matrix = np.zeros((n_l, n_l, n_sample))

    for k in range(n_sample):
        v2d = vec_2d_list[k]
        if v2d is None or np.all(v2d == 0):
            continue

        # Trace over l2: gamma_{l1, l1'} = sum_{l2} <v_{l1,l2} | v_{l1',l2}>
        # Weighted: gamma_w_{l1,l1'} = sum_{l2} P_{l2}^2 <v_{l1,l2}|v_{l1',l2}>
        for l2, group in l2_groups.items():
            w_l2 = l2_weight[l2]
            for ic_a, l1_a in group:
                i_a = l_to_idx[l1_a]
                va = v2d[ic_a, :]  # alpha profile
                for ic_b, l1_b in group:
                    i_b = l_to_idx[l1_b]
                    vb = v2d[ic_b, :]
                    dot_ab = np.dot(va, vb)
                    gamma_matrix[i_a, i_b, k] += dot_ab
                    gamma_weighted_matrix[i_a, i_b, k] += w_l2 * dot_ab

    # |F|^2-weighted average
    F2 = F_values ** 2
    total_F2 = np.sum(F2)
    if total_F2 > 1e-15:
        gamma_avg = np.sum(
            gamma_matrix * F2[np.newaxis, np.newaxis, :], axis=2
        ) / total_F2
        gamma_weighted_avg = np.sum(
            gamma_weighted_matrix * F2[np.newaxis, np.newaxis, :], axis=2
        ) / total_F2
    else:
        gamma_avg = np.zeros((n_l, n_l))
        gamma_weighted_avg = np.zeros((n_l, n_l))

    # Diagnostic: off-diagonal ratio
    diag_norm = np.linalg.norm(np.diag(gamma_avg))
    offdiag = gamma_avg - np.diag(np.diag(gamma_avg))
    offdiag_norm = np.linalg.norm(offdiag)
    offdiag_ratio = offdiag_norm / diag_norm if diag_norm > 1e-15 else 0.0

    return {
        'gamma_matrix': gamma_matrix,
        'gamma_avg': gamma_avg,
        'gamma_weighted_avg': gamma_weighted_avg,
        'l_values': l_values,
        'Re_samples': Re_samples,
        'F_values': F_values,
        'channels': channels,
        'ch_weights': ch_weights,
        'offdiag_ratio': float(offdiag_ratio),
    }


def _compute_l1_indexed_f0_matrix(
    channel_data: Dict[str, Any],
    R: float,
    n_r: int = 300,
    r_max: float = 10.0,
) -> Dict[str, Any]:
    """
    Compute F^0 matrix indexed by unique l1 values (not channels).

    Aggregates per-channel Be-centered densities P_{(l1,l2)}(d) over l2
    to get l1-marginal densities P_{l1}(d), then computes the Slater F^0
    integral between each pair (l1, l1').

    Parameters
    ----------
    channel_data : dict
        Output of extract_channel_data().
    R : float
        Internuclear distance (bohr).
    n_r : int
        Radial grid size.
    r_max : float
        Maximum distance (bohr).

    Returns
    -------
    result : dict
        F0_l1_matrix : ndarray (n_l, n_l) — F^0_{l1, l1'}
        l_values : ndarray — unique l1 values
        P_l1_Be : list of ndarray — l1-marginal Be-centered densities
        d_grid : ndarray — Be-centered radial grid
    """
    # Get per-channel F^0 data (includes Be-centered densities)
    f0_result = compute_channel_f0_matrix(
        channel_data, R, n_r=n_r, r_max=r_max)

    channels = f0_result['channels']
    P_ch_Be = f0_result['P_ch_Be']
    d_grid = f0_result['d_grid']

    # Build mapping: unique l1 values
    l1_all = [ch[0] for ch in channels]
    l_values = np.array(sorted(set(l1_all)))
    n_l = len(l_values)
    l_to_idx = {l: i for i, l in enumerate(l_values)}

    # Aggregate densities over l2 for each l1
    n_d = len(d_grid)
    P_l1_Be = [np.zeros(n_d) for _ in range(n_l)]
    for ic, ch in enumerate(channels):
        l1 = ch[0]
        idx = l_to_idx[l1]
        P_l1_Be[idx] += P_ch_Be[ic]

    # Compute F^0 matrix between l1-marginal densities
    F0_l1_matrix = np.zeros((n_l, n_l))
    for i in range(n_l):
        for j in range(i, n_l):
            f0_ij = slater_f0_integral(d_grid, P_l1_Be[i], P_l1_Be[j])
            F0_l1_matrix[i, j] = f0_ij
            F0_l1_matrix[j, i] = f0_ij

    return {
        'F0_l1_matrix': F0_l1_matrix,
        'l_values': l_values,
        'P_l1_Be': P_l1_Be,
        'd_grid': d_grid,
    }


def full_exchange_inter_fiber_energy(
    level4_result: Dict[str, Any],
    R: float,
    Z_A: float,
    Z_B: float,
    l_max: int,
    n_alpha: int,
    pk_potentials: Optional[List[dict]] = None,
    n_r: int = 300,
    r_max: float = 10.0,
    n_sample_Re: int = 10,
    use_approximate_f0: bool = False,
    bond_angle: float = np.pi,
) -> Dict[str, Any]:
    """
    Compute full exchange inter-fiber energy using the off-diagonal 1-RDM.

    The full exchange between fibers A and B in channel representation:

        K_AB(R, θ) = -sum_{l1,l1'} gamma^A_{l1,l1'}
                     * gamma^B_{l1',l1}(θ) * F^0_{l1,l1'}(R)

    General bond angle transform:
        gamma^B_{l1',l1}(θ) = P_{l1'}(cos θ) * P_{l1}(cos θ)
                              * gamma^A_weighted_{l1',l1}(θ)

    where gamma^A_weighted traces over l2 with weight P_{l2}(cos θ)^2:
        gamma^A_weighted_{l1,l1'}(θ) = sum_{l2} P_{l2}(cos θ)^2
                                        * sum_j v[(l1,l2),j] * v[(l1',l2),j]

    For θ = π: P_l(-1)^2 = 1 for all l, so gamma_weighted = gamma and
    gamma^B = (-1)^{l1+l1'} * gamma^A (the existing linear formula).

    Parameters
    ----------
    level4_result : dict
        Output of solve_level4_h2_multichannel().
    R : float
        Internuclear distance (bohr).
    Z_A, Z_B : float
        Nuclear charges.
    l_max : int
        Angular momentum limit.
    n_alpha : int
        Alpha grid points.
    pk_potentials : list or None
        Phillips-Kleinman pseudopotential.
    n_r : int
        Radial grid points.
    r_max : float
        Maximum radius (bohr).
    n_sample_Re : int
        Number of R_e sample points.
    use_approximate_f0 : bool
        If True, use F^0_avg for all (l1,l1') pairs. Default False.
    bond_angle : float
        Angle between bond axes (radians). Default π (linear).

    Returns
    -------
    result : dict
        E_exchange : float — full exchange energy (negative = attractive)
        E_exchange_diag : float — diagonal-only contribution
        E_exchange_offdiag : float — off-diagonal contribution
        gamma_avg : ndarray (n_l, n_l) — averaged 1-RDM
        parity_matrix : ndarray (n_l, n_l) — P_{l1}(cos θ) * P_{l1'}(cos θ)
        F0_l1_matrix : ndarray (n_l, n_l) — l1-indexed F^0
        l_values : ndarray — unique l1 values
        offdiag_ratio : float — ||gamma_offdiag|| / ||gamma_diag||
    """
    # Step 1: Extract channel data (shared for 1-RDM and F^0)
    channel_data = extract_channel_data(
        level4_result, R, Z_A, Z_B, l_max, n_alpha,
        pk_potentials=pk_potentials,
        n_sample_Re=n_sample_Re,
    )

    # Step 2: Extract 1-RDM (with l2-weighted version for general angles)
    rdm = extract_channel_1rdm(
        level4_result, R, Z_A, Z_B, l_max, n_alpha,
        pk_potentials=pk_potentials,
        n_sample_Re=n_sample_Re,
        channel_data=channel_data,
        bond_angle=bond_angle,
    )
    gamma_avg = rdm['gamma_avg']  # (n_l, n_l)
    gamma_weighted_avg = rdm['gamma_weighted_avg']  # l2-weighted
    l_values = rdm['l_values']
    n_l = len(l_values)

    # Step 3: l1-rotation matrix: P_{l1}(cos θ) * P_{l1'}(cos θ)
    # For θ=π: (-1)^{l1+l1'} (backward compatible)
    cos_theta = np.cos(bond_angle)
    parity_matrix = np.zeros((n_l, n_l))
    for i, l1 in enumerate(l_values):
        for j, l1p in enumerate(l_values):
            parity_matrix[i, j] = float(
                eval_legendre(l1, cos_theta) * eval_legendre(l1p, cos_theta))

    # Step 4: F^0 matrix
    if use_approximate_f0:
        # Approximate: use total F^0 scaled by density fractions
        mono = monopole_inter_fiber_energy(
            level4_result, R, Z_A, Z_B, l_max, n_alpha,
            pk_potentials=pk_potentials,
            n_r=n_r, r_max=r_max, n_sample_Re=n_sample_Re,
            method='algebraic', channel_data=channel_data,
        )
        F0_avg = mono['E_monopole']
        F0_l1_matrix = F0_avg * np.ones((n_l, n_l))
    else:
        # Exact: l1-indexed F^0 from aggregated channel densities
        f0_l1 = _compute_l1_indexed_f0_matrix(
            channel_data, R, n_r=n_r, r_max=r_max)
        F0_l1_matrix = f0_l1['F0_l1_matrix']

    # Step 5: Full exchange contraction (general bond angle)
    # K_AB(θ) = -sum_{l1,l1'} gamma^A_{l1,l1'}
    #           * P_{l1}(cos θ) * P_{l1'}(cos θ)
    #           * gamma^A_weighted_{l1',l1}(θ) * F^0_{l1,l1'}
    #
    # For θ=π: gamma_weighted = gamma and parity = (-1)^{l1+l1'},
    # recovering K_AB = -sum gamma^2 * (-1)^{l1+l1'} * F^0
    integrand = gamma_avg * gamma_weighted_avg * parity_matrix * F0_l1_matrix
    E_exchange = -np.sum(integrand)

    # Decompose into diagonal and off-diagonal contributions
    E_diag = -np.sum(np.diag(integrand))
    E_offdiag = E_exchange - E_diag

    return {
        'E_exchange': float(E_exchange),
        'E_exchange_diag': float(E_diag),
        'E_exchange_offdiag': float(E_offdiag),
        'gamma_avg': gamma_avg,
        'parity_matrix': parity_matrix,
        'F0_l1_matrix': F0_l1_matrix,
        'l_values': l_values,
        'offdiag_ratio': rdm['offdiag_ratio'],
        'channel_data': channel_data,
    }


def kinetic_orthogonalization_energy(
    level4_result: Dict[str, Any],
    R: float,
    Z_A: float,
    Z_B: float,
    l_max: int,
    n_alpha: int,
    pk_potentials: Optional[List[dict]] = None,
    n_sample_Re: int = 10,
    bond_angle: float = np.pi,
) -> Dict[str, Any]:
    """
    Compute the kinetic energy correction from Löwdin orthogonalization
    of overlapping fiber wavefunctions.

    The inter-fiber overlap in the channel basis for general bond angle θ:

        S^{AB}_{l1,l1'}(θ) = P_{l1'}(cos θ) * gamma^A_weighted_{l1,l1'}(θ)

    where gamma^A_weighted traces over l2 with weight P_{l2}(cos θ)^2.
    For θ=π: P_{l1'}(-1) = (-1)^{l1'} and gamma_weighted = gamma,
    recovering S^{AB} = (-1)^{l1'} * gamma^A.

    To leading order in S, the kinetic energy correction is:

        DeltaT(R) = +(1/2) * sum_{l1,l1'} |S^{AB}_{l1,l1'}|^2
                    * (T_{l1} + T_{l1'}) / 2

    Parameters
    ----------
    level4_result : dict
        Output of solve_level4_h2_multichannel().
    R : float
        Internuclear distance (bohr).
    Z_A, Z_B : float
        Nuclear charges.
    l_max : int
        Angular momentum limit.
    n_alpha : int
        Alpha grid points.
    pk_potentials : list or None
        Phillips-Kleinman pseudopotential.
    n_sample_Re : int
        Number of R_e sample points.
    bond_angle : float
        Angle between bond axes (radians). Default π (linear).

    Returns
    -------
    result : dict
        delta_T : float — kinetic orthogonalization correction (positive,
            repulsive).
        S_AB_matrix : ndarray (n_l, n_l) — inter-fiber overlap matrix.
        S_AB_frobenius : float — ||S^{AB}||_F (Frobenius norm of overlap).
        T_l : ndarray (n_l,) — centrifugal kinetic energy per l1.
        Re2_avg : float — <R_e^2> from |F|^2 weighting.
        E_val : float — valence eigenvalue (for scaling reference).
        l_values : ndarray — unique l1 values.
        gamma_avg : ndarray (n_l, n_l) — the 1-RDM used.
    """
    # Step 1: Extract channel data (shared for 1-RDM)
    channel_data = extract_channel_data(
        level4_result, R, Z_A, Z_B, l_max, n_alpha,
        pk_potentials=pk_potentials,
        n_sample_Re=n_sample_Re,
    )

    # Step 2: Extract 1-RDM (with l2-weighted version for general angles)
    rdm = extract_channel_1rdm(
        level4_result, R, Z_A, Z_B, l_max, n_alpha,
        pk_potentials=pk_potentials,
        n_sample_Re=n_sample_Re,
        channel_data=channel_data,
        bond_angle=bond_angle,
    )
    gamma_avg = rdm['gamma_avg']  # (n_l, n_l)
    gamma_weighted_avg = rdm['gamma_weighted_avg']
    l_values = rdm['l_values']
    n_l = len(l_values)

    # Step 3: Construct inter-fiber overlap matrix
    # S^{AB}_{l1,l1'}(θ) = P_{l1'}(cos θ) * gamma_weighted_{l1,l1'}(θ)
    # For θ=π: P_{l1'}(-1) = (-1)^{l1'} and gamma_weighted = gamma.
    cos_theta = np.cos(bond_angle)
    parity_col = np.array(
        [float(eval_legendre(l1p, cos_theta)) for l1p in l_values],
        dtype=float)
    S_AB = gamma_weighted_avg * parity_col[np.newaxis, :]

    # Step 4: Compute <R_e^2> from |F(R_e)|^2 weighting
    Re_samples = channel_data['Re_samples']
    F_values = channel_data['F_values']
    F2 = F_values ** 2
    total_F2 = np.sum(F2)
    if total_F2 > 1e-15:
        Re2_avg = float(np.sum(Re_samples ** 2 * F2) / total_F2)
    else:
        Re2_avg = 1.0  # fallback

    # Step 5: Centrifugal kinetic energy per channel
    # T_{l1} = l1*(l1+1) / (2 * <R_e^2>)
    T_l = np.array([l1 * (l1 + 1) / (2.0 * Re2_avg)
                     for l1 in l_values], dtype=float)

    # Step 6: Kinetic correction
    # DeltaT = +(1/2) * sum_{l1,l1'} |S^{AB}_{l1,l1'}|^2 * (T_{l1}+T_{l1'})/2
    S_AB_sq = S_AB ** 2
    T_avg_matrix = (T_l[:, np.newaxis] + T_l[np.newaxis, :]) / 2.0
    delta_T = 0.5 * float(np.sum(S_AB_sq * T_avg_matrix))

    # Valence eigenvalue for scaling reference
    E_val = float(level4_result['E_elec'])

    S_AB_frobenius = float(np.linalg.norm(S_AB, 'fro'))

    return {
        'delta_T': delta_T,
        'S_AB_matrix': S_AB,
        'S_AB_frobenius': S_AB_frobenius,
        'T_l': T_l,
        'Re2_avg': Re2_avg,
        'E_val': E_val,
        'l_values': l_values,
        'gamma_avg': gamma_avg,
    }
