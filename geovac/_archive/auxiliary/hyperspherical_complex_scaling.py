"""
Exterior Complex Scaling (ECS) for hyperspherical resonances.

Under the transformation R → R₀ + (R - R₀)e^{iθ} for R > R₀,
the Hamiltonian becomes non-Hermitian. Continuum states rotate
into the complex plane by angle -2θ, exposing resonances as
isolated complex eigenvalues:

    E_resonance = E_res - i Γ/2

The width Γ is the imaginary part of a complex eigenvalue — an
eigenvalue, not an integral. This is the natural endpoint of the
"eigenvalue not integral" philosophy applied to resonance physics.

References:
  - Simon, Phys. Lett. A 71, 211 (1979) — complex scaling theory
  - Moiseyev, Phys. Rep. 302, 212 (1998) — non-Hermitian QM review
  - Rescigno & McCurdy, Phys. Rev. A 62, 032706 (2000) — ECS method
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigs
from typing import Tuple, List, Dict, Optional

from geovac.hyperspherical_adiabatic import compute_adiabatic_curve, effective_potential
from geovac._archive.auxiliary.hyperspherical_resonances import (
    E_HE_EXACT, E_HE_PLUS, HARTREE_TO_EV, EXPERIMENTAL_RESONANCES,
)


def _smooth_theta(R: np.ndarray, R0: float, delta_R: float, theta: float) -> np.ndarray:
    """
    Smooth ECS scaling angle: ramps from 0 to theta over width delta_R.

    Uses a Fermi function: theta_j = theta / (1 + exp(-(R_j - R0) / (delta_R/6)))
    This gives ~0 for R << R0 and ~theta for R >> R0, with smooth transition.

    Parameters
    ----------
    R : ndarray
        Grid points.
    R0 : float
        ECS turning point.
    delta_R : float
        Transition width.
    theta : float
        Asymptotic rotation angle.

    Returns
    -------
    theta_arr : ndarray
        Local scaling angle at each grid point.
    """
    x = (R - R0) / (delta_R / 6.0)
    # Clip to avoid overflow in exp
    x = np.clip(x, -50.0, 50.0)
    return theta / (1.0 + np.exp(-x))


def solve_ecs_single_channel(
    V_eff_spline: CubicSpline,
    mu_spline: Optional[CubicSpline] = None,
    R_min: float = 0.05,
    R_max: float = 40.0,
    R0: float = 15.0,
    theta: float = 0.3,
    delta_R: float = 2.0,
    N_R: int = 3000,
    n_states: int = 20,
    sigma: complex = -2.5 + 0.0j,
) -> Dict:
    """
    Solve the complex-scaled single-channel hyperradial equation.

    The smooth ECS transformation applies a position-dependent rotation
    angle θ(R) that ramps from 0 (inner region) to θ (outer region)
    over a width δR centered at R₀.

    At each grid point, the kinetic energy picks up a factor e^{-2iθ_j}
    and the potential is analytically continued via the complex coordinate.

    Parameters
    ----------
    V_eff_spline : CubicSpline
        Real effective potential V_eff(R) from the adiabatic solver.
    mu_spline : CubicSpline, optional
        Angular eigenvalue mu(R). If provided, used for analytic continuation
        of V_eff in the complex region. Otherwise V_eff_spline is used directly.
    R_min : float
        Left boundary (bohr).
    R_max : float
        Right boundary (bohr).
    R0 : float
        ECS turning point (bohr). Should be beyond the wavefunction center
        but inside the grid.
    theta : float
        Rotation angle (radians). Typically 0.2-0.5.
    delta_R : float
        Smooth transition width (bohr).
    N_R : int
        Number of interior grid points.
    n_states : int
        Number of eigenvalues to find.
    sigma : complex
        Shift for shift-invert eigenvalue search.

    Returns
    -------
    result : dict
        Keys:
        - 'eigenvalues': ndarray of complex eigenvalues
        - 'R_grid': ndarray of real grid points
        - 'theta_grid': ndarray of local scaling angles
        - 'R0': float
        - 'theta': float
        - 'N_R': int
    """
    h = (R_max - R_min) / (N_R + 1)
    R_grid = R_min + (np.arange(N_R) + 1) * h

    # Compute local scaling angle at each grid point
    theta_grid = _smooth_theta(R_grid, R0, delta_R, theta)

    # Complex scaling factors
    exp_neg2i_theta = np.exp(-2j * theta_grid)

    # Build complex potential on the grid
    #
    # Key subtlety: V_eff(R) = μ(R)/R² + 15/(8R²) → E_threshold as R → ∞.
    # Under naive complex scaling V(z) = μ/z² + 15/(8z²), the threshold
    # also rotates, giving spurious positive imaginary parts. The correct
    # approach separates the asymptotic constant:
    #
    #   V_eff(R) = E_threshold + ΔV(R),  where ΔV → 0 as R → ∞
    #   V_ECS(j) = E_threshold + ΔV(z_j)  [threshold stays real]
    #
    # Since μ(R) + 15/8 = E_threshold·R² + Δμ(R) where Δμ → 0:
    #   ΔV(z) = Δμ(R_real)/z² = [μ(R) + 15/8 - E_threshold·R²]/z²
    #
    E_threshold = E_HE_PLUS  # -Z²/2 = -2.0 for He
    V_complex = np.zeros(N_R, dtype=complex)
    for j in range(N_R):
        if theta_grid[j] < 1e-10:
            # Inner region: real potential
            V_complex[j] = float(V_eff_spline(R_grid[j]))
        else:
            # Outer region: threshold-separated analytic continuation
            z_j = R0 + (R_grid[j] - R0) * np.exp(1j * theta_grid[j])
            R_j = R_grid[j]
            if mu_spline is not None:
                mu_val = float(mu_spline(R_j))
                # Δμ = μ + 15/8 - E_threshold·R²  (decaying part)
                delta_mu = mu_val + 15.0 / 8.0 - E_threshold * R_j**2
                V_complex[j] = E_threshold + delta_mu / z_j**2
            else:
                # Approximate: separate threshold from V_eff
                V_real = float(V_eff_spline(R_j))
                delta_V = V_real - E_threshold
                V_complex[j] = E_threshold + delta_V * (R_j / z_j)**2

    # Build complex tridiagonal Hamiltonian using sparse matrix
    H = lil_matrix((N_R, N_R), dtype=complex)

    for j in range(N_R):
        # Diagonal: scaled kinetic + complex potential
        H[j, j] = exp_neg2i_theta[j] / h**2 + V_complex[j]

        # Off-diagonal: scaled kinetic
        if j < N_R - 1:
            # Use geometric mean of scaling factors at j and j+1
            scale = np.sqrt(exp_neg2i_theta[j] * exp_neg2i_theta[j + 1])
            H[j, j + 1] = -0.5 * scale / h**2
            H[j + 1, j] = -0.5 * scale / h**2

    H_csr = csr_matrix(H)

    # Complex shift-invert eigenvalue solve
    evals, evecs = eigs(H_csr, k=n_states, sigma=sigma, which='LM')

    # Sort by real part
    idx = np.argsort(evals.real)
    evals = evals[idx]

    return {
        'eigenvalues': evals,
        'R_grid': R_grid,
        'theta_grid': theta_grid,
        'R0': R0,
        'theta': theta,
        'N_R': N_R,
    }


def solve_ecs_coupled(
    V_eff_splines: list,
    P_splines: list,
    mu_splines: Optional[list] = None,
    n_channels: int = 3,
    R_min: float = 0.05,
    R_max: float = 40.0,
    R0: float = 12.0,
    theta: float = 0.3,
    delta_R: float = 2.0,
    N_R: int = 2000,
    n_states: int = 50,
    sigma: complex = -2.5 + 0.0j,
) -> Dict:
    """
    Solve the complex-scaled coupled-channel hyperradial equation.

    Extends the single-channel ECS to coupled channels. The block-tridiagonal
    Hamiltonian from the coupled-channel solver gets ECS treatment:
    - Inner blocks: unchanged (real)
    - Outer blocks: kinetic scaled by e^{-2iθ}, potential analytically continued
    - P_μν coupling: scaled by e^{-iθ} (first derivative picks up one factor)

    Parameters
    ----------
    V_eff_splines : list of CubicSpline
        Real effective potentials for each channel.
    P_splines : list of lists of CubicSpline
        Non-adiabatic coupling P_μν(R) for each channel pair.
    mu_splines : list of CubicSpline, optional
        Angular eigenvalue splines for analytic continuation.
    n_channels : int
        Number of coupled channels.
    R_min, R_max, R0, theta, delta_R, N_R, n_states, sigma :
        Same as solve_ecs_single_channel.

    Returns
    -------
    result : dict
        Same keys as single-channel, plus 'n_channels'.
    """
    n_ch = n_channels
    h = (R_max - R_min) / (N_R + 1)
    R_grid = R_min + (np.arange(N_R) + 1) * h

    # Local scaling angles
    theta_grid = _smooth_theta(R_grid, R0, delta_R, theta)
    exp_neg2i_theta = np.exp(-2j * theta_grid)
    exp_neg1i_theta = np.exp(-1j * theta_grid)

    # Evaluate potentials on grid (threshold-separated analytic continuation)
    E_threshold = E_HE_PLUS  # -Z²/2 = -2.0 for He
    V_complex = np.zeros((n_ch, N_R), dtype=complex)
    for ch in range(n_ch):
        for j in range(N_R):
            if theta_grid[j] < 1e-10:
                V_complex[ch, j] = float(V_eff_splines[ch](R_grid[j]))
            else:
                z_j = R0 + (R_grid[j] - R0) * np.exp(1j * theta_grid[j])
                R_j = R_grid[j]
                if mu_splines is not None:
                    mu_val = float(mu_splines[ch](R_j))
                    delta_mu = mu_val + 15.0 / 8.0 - E_threshold * R_j**2
                    V_complex[ch, j] = E_threshold + delta_mu / z_j**2
                else:
                    V_real = float(V_eff_splines[ch](R_j))
                    delta_V = V_real - E_threshold
                    V_complex[ch, j] = E_threshold + delta_V * (R_j / z_j)**2

    # Evaluate P coupling on grid
    P_grid = np.zeros((n_ch, n_ch, N_R))
    for mu in range(n_ch):
        for nu in range(n_ch):
            P_grid[mu, nu] = P_splines[mu][nu](R_grid)

    # Build complex block-tridiagonal sparse matrix
    dim = N_R * n_ch
    H = lil_matrix((dim, dim), dtype=complex)

    for i in range(N_R):
        scale_2 = exp_neg2i_theta[i]
        scale_1 = exp_neg1i_theta[i]

        for mu in range(n_ch):
            row = i * n_ch + mu

            # Diagonal block: scaled kinetic + complex potential
            H[row, row] = scale_2 / h**2 + V_complex[mu, i]

        # Off-diagonal in R: coupling to i+1
        if i < N_R - 1:
            scale_2_avg = np.sqrt(exp_neg2i_theta[i] * exp_neg2i_theta[i + 1])
            scale_1_avg = np.sqrt(exp_neg1i_theta[i] * exp_neg1i_theta[i + 1])

            for mu in range(n_ch):
                row = i * n_ch + mu
                for nu in range(n_ch):
                    col = (i + 1) * n_ch + nu
                    val = 0.0 + 0.0j
                    if mu == nu:
                        val += -0.5 * scale_2_avg / h**2
                    # P coupling: scaled by e^{-iθ}
                    val -= scale_1_avg * P_grid[mu, nu, i] / (2.0 * h)
                    if abs(val) > 1e-16:
                        H[row, col] = val

        # Off-diagonal in R: coupling to i-1
        if i > 0:
            scale_2_avg = np.sqrt(exp_neg2i_theta[i] * exp_neg2i_theta[i - 1])
            scale_1_avg = np.sqrt(exp_neg1i_theta[i] * exp_neg1i_theta[i - 1])

            for mu in range(n_ch):
                row = i * n_ch + mu
                for nu in range(n_ch):
                    col = (i - 1) * n_ch + nu
                    val = 0.0 + 0.0j
                    if mu == nu:
                        val += -0.5 * scale_2_avg / h**2
                    # P coupling: opposite sign for backward derivative
                    val += scale_1_avg * P_grid[mu, nu, i] / (2.0 * h)
                    if abs(val) > 1e-16:
                        H[row, col] = val

    H_csr = csr_matrix(H)

    evals, evecs = eigs(H_csr, k=n_states, sigma=sigma, which='LM')

    idx = np.argsort(evals.real)
    evals = evals[idx]
    evecs = evecs[:, idx]

    # Compute channel norms: ||ψ_μ||² for each channel μ
    channel_norms = np.zeros((n_ch, n_states))
    for k in range(n_states):
        vec = evecs[:, k]
        for mu in range(n_ch):
            # Extract channel μ block: indices mu, mu+n_ch, mu+2*n_ch, ...
            ch_vec = vec[mu::n_ch]
            channel_norms[mu, k] = np.sum(np.abs(ch_vec)**2)

    return {
        'eigenvalues': evals,
        'eigenvectors': evecs,
        'channel_norms': channel_norms,
        'R_grid': R_grid,
        'theta_grid': theta_grid,
        'R0': R0,
        'theta': theta,
        'N_R': N_R,
        'n_channels': n_ch,
    }


def identify_resonances(
    eigenvalues: np.ndarray,
    theta: float,
    E_threshold: float = -2.0,
    angle_tolerance: float = 0.15,
    imag_threshold: float = -1e-5,
) -> List[Dict]:
    """
    Separate resonances from rotated continuum in the complex eigenvalue spectrum.

    Continuum states lie along the line arg(E - E_threshold) ≈ -2θ.
    Resonances deviate from this line and have small negative imaginary parts.

    Parameters
    ----------
    eigenvalues : ndarray of complex
        Complex eigenvalues from ECS solver.
    theta : float
        ECS rotation angle (radians).
    E_threshold : float
        Continuum threshold energy (Ha). Default: He+ at -2.0.
    angle_tolerance : float
        Max angular deviation from -2θ line to be considered continuum (radians).
    imag_threshold : float
        Maximum Im(E) for bound states (should be ~ 0).

    Returns
    -------
    resonances : list of dict
        Each dict has keys:
        - 'E_res': resonance position (Ha)
        - 'Gamma': width (Ha)
        - 'Gamma_eV': width (eV)
        - 'E_eV_above_gs': energy above He ground state (eV)
        - 'eigenvalue': the raw complex eigenvalue
        - 'angle': arg(E - E_threshold)
        - 'classification': 'bound', 'resonance', or 'continuum'
    """
    continuum_angle = -2.0 * theta
    resonances = []

    for E in eigenvalues:
        E_re = E.real
        E_im = E.imag

        # Classify
        if E_re < E_threshold and abs(E_im) < abs(imag_threshold):
            classification = 'bound'
        elif E_re > E_threshold and E_im < imag_threshold:
            # Above threshold with negative imaginary part
            dE = E - E_threshold
            angle = np.angle(dE)

            if abs(angle - continuum_angle) < angle_tolerance:
                classification = 'continuum'
            else:
                classification = 'resonance'
        else:
            classification = 'continuum'

        Gamma = -2.0 * E_im if E_im < 0 else 0.0

        resonances.append({
            'E_res': E_re,
            'Gamma': Gamma,
            'Gamma_eV': Gamma * HARTREE_TO_EV,
            'E_eV_above_gs': (E_re - E_HE_EXACT) * HARTREE_TO_EV,
            'eigenvalue': E,
            'angle': np.angle(E - E_threshold) if abs(E - E_threshold) > 1e-10 else 0.0,
            'classification': classification,
        })

    return resonances


def theta_stability_scan(
    V_eff_spline: CubicSpline,
    mu_spline: Optional[CubicSpline] = None,
    theta_values: Optional[List[float]] = None,
    R_min: float = 0.05,
    R_max: float = 40.0,
    R0: float = 15.0,
    delta_R: float = 2.0,
    N_R: int = 3000,
    n_states: int = 20,
    sigma: complex = -2.5 + 0.0j,
    stability_tol: float = 0.005,
) -> Dict:
    """
    Scan over multiple θ values to identify θ-stable resonances.

    True resonances have position and width independent of θ.
    Continuum states move linearly with θ.

    Parameters
    ----------
    V_eff_spline, mu_spline : solver inputs
    theta_values : list of float, optional
        Rotation angles to scan. Default: [0.15, 0.20, 0.25, 0.30, 0.35, 0.40].
    R_min, R_max, R0, delta_R, N_R, n_states, sigma : solver parameters
    stability_tol : float
        Max variation in Re(E) (Ha) for an eigenvalue to be θ-stable.

    Returns
    -------
    result : dict
        Keys:
        - 'theta_values': list of float
        - 'all_eigenvalues': list of ndarrays
        - 'stable_resonances': list of dicts with θ-stable resonance candidates
        - 'ground_state_energies': list of complex ground state at each θ
    """
    if theta_values is None:
        theta_values = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40]

    all_eigenvalues = []
    ground_states = []

    for th in theta_values:
        result = solve_ecs_single_channel(
            V_eff_spline, mu_spline,
            R_min=R_min, R_max=R_max, R0=R0,
            theta=th, delta_R=delta_R,
            N_R=N_R, n_states=n_states, sigma=sigma,
        )
        evals = result['eigenvalues']
        all_eigenvalues.append(evals)

        # Ground state: lowest real part with small imaginary part
        bound_mask = (evals.real < E_HE_PLUS) & (np.abs(evals.imag) < 0.01)
        if np.any(bound_mask):
            ground_states.append(evals[bound_mask][np.argmin(evals[bound_mask].real)])
        else:
            ground_states.append(evals[np.argmin(evals.real)])

    # Find θ-stable eigenvalues above threshold
    # For each eigenvalue at θ₀, check if a similar eigenvalue exists at all other θ
    ref_idx = len(theta_values) // 2  # Use middle θ as reference
    ref_evals = all_eigenvalues[ref_idx]
    ref_theta = theta_values[ref_idx]

    stable_resonances = []
    for E_ref in ref_evals:
        if E_ref.real < E_HE_PLUS:
            continue  # Skip bound states
        if E_ref.imag > -1e-6:
            continue  # Skip states with positive/zero imaginary part

        # Check stability across θ values
        matches = []
        for i, (th, evals) in enumerate(zip(theta_values, all_eigenvalues)):
            # Find closest eigenvalue to E_ref
            diffs = np.abs(evals - E_ref)
            best_idx = np.argmin(diffs)
            best_match = evals[best_idx]
            matches.append(best_match)

        # Check if Re(E) is stable
        re_parts = np.array([m.real for m in matches])
        im_parts = np.array([m.imag for m in matches])
        re_spread = np.max(re_parts) - np.min(re_parts)
        im_spread = np.max(im_parts) - np.min(im_parts)

        if re_spread < stability_tol:
            E_mean = np.mean(re_parts)
            Gamma_mean = -2.0 * np.mean(im_parts)

            # Check it's not on the -2θ line
            ref_angle = np.angle(E_ref - E_HE_PLUS)
            expected_continuum = -2.0 * ref_theta
            angle_dev = abs(ref_angle - expected_continuum)

            if angle_dev > 0.1:  # Not on continuum line
                stable_resonances.append({
                    'E_res': E_mean,
                    'Gamma': Gamma_mean,
                    'Gamma_eV': Gamma_mean * HARTREE_TO_EV,
                    'E_eV_above_gs': (E_mean - E_HE_EXACT) * HARTREE_TO_EV,
                    'Re_spread': re_spread,
                    'Im_spread': im_spread,
                    'matches': matches,
                    'n_theta': len(theta_values),
                })

    # Sort by energy
    stable_resonances.sort(key=lambda r: r['E_res'])

    return {
        'theta_values': theta_values,
        'all_eigenvalues': all_eigenvalues,
        'stable_resonances': stable_resonances,
        'ground_state_energies': ground_states,
    }


def run_ecs_analysis(
    Z: float = 2.0,
    l_max: int = 2,
    n_alpha: int = 100,
    n_channels: int = 5,
    N_R_angular: int = 200,
    R0: float = 15.0,
    theta: float = 0.3,
    delta_R: float = 2.0,
    N_R: int = 3000,
    n_states: int = 25,
    sigma: complex = -2.5 + 0.0j,
    verbose: bool = True,
) -> Dict:
    """
    Full ECS analysis pipeline: compute adiabatic curves, run ECS solver,
    identify resonances, compare with experiment.

    Parameters
    ----------
    Z, l_max, n_alpha, n_channels : angular solver parameters
    N_R_angular : int
        R grid for angular solve.
    R0, theta, delta_R, N_R, n_states, sigma : ECS solver parameters
    verbose : bool
        Print progress.

    Returns
    -------
    result : dict
        Full analysis results.
    """
    import time
    t0 = time.time()

    # Step 1: Compute adiabatic curves on dense grid
    R_ang_max = max(R0 + 5.0, 25.0)
    R_grid_ang = np.concatenate([
        np.linspace(0.1, 1.0, N_R_angular // 3),
        np.linspace(1.0, 5.0, N_R_angular // 3),
        np.linspace(5.0, R_ang_max, N_R_angular // 3 + 1),
    ])
    R_grid_ang = np.unique(R_grid_ang)

    if verbose:
        print(f"Computing adiabatic curves ({n_channels} channels)...")

    mu_curves = compute_adiabatic_curve(R_grid_ang, Z, l_max, n_alpha, n_channels)

    # Build splines for V_eff and mu
    V_eff_splines = []
    mu_splines = []
    for ch in range(n_channels):
        V_eff_ch = effective_potential(R_grid_ang, mu_curves[ch])
        V_eff_splines.append(CubicSpline(R_grid_ang, V_eff_ch, extrapolate=True))
        mu_splines.append(CubicSpline(R_grid_ang, mu_curves[ch], extrapolate=True))

    t1 = time.time()
    if verbose:
        print(f"  Angular solve: {t1 - t0:.1f}s")

    # Step 2: Single-channel ECS on lowest channel (validation)
    if verbose:
        print(f"Running single-channel ECS (R0={R0}, theta={theta:.2f})...")

    sc_result = solve_ecs_single_channel(
        V_eff_splines[0], mu_splines[0],
        R_min=0.05, R_max=40.0, R0=R0,
        theta=theta, delta_R=delta_R,
        N_R=N_R, n_states=n_states, sigma=sigma,
    )

    sc_resonances = identify_resonances(sc_result['eigenvalues'], theta)

    t2 = time.time()
    if verbose:
        n_bound = sum(1 for r in sc_resonances if r['classification'] == 'bound')
        n_res = sum(1 for r in sc_resonances if r['classification'] == 'resonance')
        n_cont = sum(1 for r in sc_resonances if r['classification'] == 'continuum')
        print(f"  Single-channel ECS: {t2 - t1:.1f}s")
        print(f"  Found: {n_bound} bound, {n_res} resonance, {n_cont} continuum")

        # Report ground state
        bound = [r for r in sc_resonances if r['classification'] == 'bound']
        if bound:
            gs = min(bound, key=lambda r: r['E_res'])
            print(f"  Ground state: {gs['E_res']:.6f} Ha "
                  f"(Im = {gs['eigenvalue'].imag:.2e})")

    # Step 3: θ-stability scan
    if verbose:
        print("Running theta-stability scan...")

    stability = theta_stability_scan(
        V_eff_splines[0], mu_splines[0],
        theta_values=[0.15, 0.20, 0.25, 0.30, 0.35, 0.40],
        R_min=0.05, R_max=40.0, R0=R0,
        delta_R=delta_R, N_R=N_R,
        n_states=n_states, sigma=sigma,
    )

    t3 = time.time()
    if verbose:
        print(f"  Stability scan: {t3 - t2:.1f}s")
        print(f"  θ-stable resonances: {len(stability['stable_resonances'])}")
        for sr in stability['stable_resonances']:
            print(f"    E = {sr['E_res']:.4f} Ha "
                  f"({sr['E_eV_above_gs']:.1f} eV above gs), "
                  f"Γ = {sr['Gamma_eV']:.4f} eV, "
                  f"Re spread = {sr['Re_spread']:.1e}")

    return {
        'single_channel': sc_result,
        'sc_resonances': sc_resonances,
        'stability': stability,
        'V_eff_splines': V_eff_splines,
        'mu_splines': mu_splines,
        'mu_curves': mu_curves,
        'R_grid_angular': R_grid_ang,
        'Z': Z,
        'l_max': l_max,
        'n_channels': n_channels,
    }
