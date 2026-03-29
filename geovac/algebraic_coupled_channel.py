"""
Coupled-channel radial solver using algebraic angular eigenvalues.

Connects the AlgebraicAngularSolver (which provides exact Hellmann-Feynman
P-matrix coupling) to the existing coupled-channel radial solver.

The key advantage over the FD-based coupling (hyperspherical_coupling.py):
dH/dR = V_coupling is precomputed and R-independent, so the Hellmann-Feynman
P-matrix is exact at every R without finite-difference noise.

Interface mapping:
    AlgebraicAngularSolver.solve(R, n_channels) → μ_ν(R), eigenvectors
    AlgebraicAngularSolver._coupling_full = dH/dR (R-independent)
    → V_eff_μ(R) = [μ_ν(R) + 15/8] / R²
    → P_μν(R) = ⟨Φ_μ|V_coupling|Φ_ν⟩ / (μ_μ - μ_ν)   (exact)

The coupled equations include the Q-matrix (second-derivative coupling)
via the closure approximation Q_μν ≈ Σ_κ P_μκ P_κν, summed over ALL
angular channels (not just the truncated set). This is essential because:
  - The P·d/dR coupling lowers the energy (channel mixing)
  - The Q term raises the energy (kinetic cost of non-adiabatic motion)
  - The Q diagonal (DBOC) and P coupling nearly cancel (~97% for He)
  - With only P and no Q, the energy always overshoots below exact

The sum over all channels in Q means the diagonal Q_μμ = 2*DBOC_μ(total),
which includes contributions from channels outside the truncated set.
These external contributions appear as a purely repulsive correction since
they have no corresponding P coupling to cancel them.

References:
  - Paper 13, Sections VII-VIII (coupled-channel structure)
  - debug/algebraic_dboc_results.md (97% DBOC cancellation)
  - Lin, Phys. Rep. 257, 1 (1995) — coupled-channel equations
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh
from typing import Tuple, Optional

from geovac.algebraic_angular import AlgebraicAngularSolver
from geovac.hyperspherical_adiabatic import effective_potential
from geovac.hyperspherical_radial import (
    solve_radial, solve_coupled_radial,
    solve_radial_spectral, solve_coupled_radial_spectral,
)


def _enforce_sign_consistency(
    vecs_prev: np.ndarray,
    vecs_curr: np.ndarray,
) -> np.ndarray:
    """Enforce eigenvector sign consistency between adjacent R points.

    Ensures ⟨Φ_μ(R_i)|Φ_μ(R_{i-1})⟩ > 0 for each channel μ.

    Parameters
    ----------
    vecs_prev : ndarray of shape (n_channels, dim)
        Sign-fixed eigenvectors at previous R point.
    vecs_curr : ndarray of shape (n_channels, dim)
        Eigenvectors at current R point.

    Returns
    -------
    vecs_fixed : ndarray of shape (n_channels, dim)
    """
    vecs_fixed = vecs_curr.copy()
    for mu in range(vecs_curr.shape[0]):
        if np.dot(vecs_prev[mu], vecs_curr[mu]) < 0:
            vecs_fixed[mu] *= -1
    return vecs_fixed


def _compute_dPdR(
    P_full: np.ndarray,
    evals_full: np.ndarray,
    V_coupling: np.ndarray,
    evecs_full: np.ndarray,
    n_ch: int,
    total_dim: int,
) -> np.ndarray:
    """Compute dP/dR for the truncated channel block, algebraically.

    Our storage convention is:
        P_stored[mu, nu] = V_ad[nu, mu] / (mu_mu - mu_nu)
    where V_ad = evecs^T @ V_coupling @ evecs is the coupling in the
    adiabatic basis. This equals -P_standard[mu, nu].

    Differentiating P_stored[mu, nu] = V_ad[nu, mu] / Delta_mu_nu:

        dP/dR = [dV_ad[nu,mu]/dR * Delta - V_ad[nu,mu] * dDelta/dR] / Delta^2

    where Delta = mu_mu - mu_nu and:
        dDelta/dR = V_ad[mu,mu] - V_ad[nu,nu]     (Hellmann-Feynman)
        dV_ad[nu,mu]/dR = sum_{k!=nu} P_s[nu,k] * V_ad[k,mu]
                        + sum_{k!=mu} P_s[mu,k] * V_ad[nu,k]

    where P_s[alpha,k] = V_ad[k,alpha] / (mu_alpha - mu_k) is the stored P
    for any pair (alpha, k) with alpha != k.

    This formula is verified to match finite-difference dP/dR to < 1e-9.

    Note: dP_stored[mu,mu]/dR = 0 identically (P_stored[mu,mu] = 0 for all R).
    This means the diagonal Q (DBOC) is exact under the closure approximation.
    The dP/dR correction only affects off-diagonal Q elements.

    Parameters
    ----------
    P_full : ndarray of shape (n_ch, total_dim)
        P_full[mu, kappa] = V_ad[kappa, mu] / (mu_mu - mu_kappa).
    evals_full : ndarray of shape (total_dim,)
        All eigenvalues at this R.
    V_coupling : ndarray of shape (total_dim, total_dim)
        R-independent coupling matrix in the original basis (dH/dR).
    evecs_full : ndarray of shape (total_dim, total_dim)
        All eigenvectors (columns) at this R.
    n_ch : int
        Number of truncated channels.
    total_dim : int
        Total number of angular channels.

    Returns
    -------
    dPdR : ndarray of shape (n_ch, n_ch)
        dP_stored[mu,nu]/dR for the truncated block.
    """
    # V in the adiabatic basis (symmetric)
    V_ad = evecs_full.T @ V_coupling @ evecs_full
    V_diag = np.diag(V_ad)

    dPdR = np.zeros((n_ch, n_ch))

    for mu_idx in range(n_ch):
        for nu_idx in range(n_ch):
            if mu_idx == nu_idx:
                continue

            Delta = evals_full[mu_idx] - evals_full[nu_idx]
            if abs(Delta) < 1e-10:
                continue

            V_nm = V_ad[nu_idx, mu_idx]

            # Sum 1: sum_{k!=nu} P_s[nu,k] * V_ad[k, mu]
            # P_s[nu, k] = V_ad[k, nu] / (mu_nu - mu_k)
            sum1 = 0.0
            for k in range(total_dim):
                if k == nu_idx:
                    continue
                gap_nk = evals_full[nu_idx] - evals_full[k]
                if abs(gap_nk) < 1e-10:
                    continue
                sum1 += (V_ad[k, nu_idx] / gap_nk) * V_ad[k, mu_idx]

            # Sum 2: sum_{k!=mu} P_s[mu,k] * V_ad[nu, k]
            # P_s[mu, k] = V_ad[k, mu] / (mu_mu - mu_k)
            sum2 = 0.0
            for k in range(total_dim):
                if k == mu_idx:
                    continue
                gap_mk = evals_full[mu_idx] - evals_full[k]
                if abs(gap_mk) < 1e-10:
                    continue
                sum2 += (V_ad[k, mu_idx] / gap_mk) * V_ad[nu_idx, k]

            dV_nm_dR = sum1 + sum2
            dDelta_dR = V_diag[mu_idx] - V_diag[nu_idx]

            dPdR[mu_idx, nu_idx] = (
                dV_nm_dR * Delta - V_nm * dDelta_dR
            ) / Delta**2

    return dPdR


def compute_algebraic_coupling(
    solver: AlgebraicAngularSolver,
    R_grid: np.ndarray,
    n_channels: int = 3,
    compute_exact_dPdR: bool = False,
) -> dict:
    """Compute adiabatic curves, P-matrix, and Q-matrix from algebraic solver.

    Uses the Hellmann-Feynman theorem with the precomputed R-independent
    V_coupling = dH/dR to get exact P-matrix elements at each R.

    The Q-matrix is computed via the closure approximation:
        Q_μν = Σ_{ALL κ} P_μκ P_κν
    summing over all angular channels (not just the truncated set).

    If compute_exact_dPdR=True, also computes the exact Q-matrix:
        Q_exact_μν = Σ_κ P_μκ P_κν + dP_μν/dR
    where dP/dR is computed algebraically from the Hellmann-Feynman quantities.
    The dP/dR correction only affects off-diagonal Q elements (dP_μμ/dR = 0
    by antisymmetry).

    Parameters
    ----------
    solver : AlgebraicAngularSolver
        Configured algebraic angular solver.
    R_grid : ndarray of shape (N_R,)
        Hyperradius grid.
    n_channels : int
        Number of adiabatic channels for the coupled equations.
    compute_exact_dPdR : bool
        If True, compute exact dP/dR and Q_exact in addition to closure Q.

    Returns
    -------
    result : dict
        Keys:
        - 'mu': ndarray (n_channels, N_R) -- adiabatic eigenvalues
        - 'P': ndarray (n_channels, n_channels, N_R) -- first-derivative coupling
        - 'Q': ndarray (n_channels, n_channels, N_R) -- closure Q (P*P over all)
        - 'DBOC_total': ndarray (n_channels, N_R) -- full DBOC (1/2 * Q diagonal)
        - 'DBOC_internal': ndarray (n_channels, N_R) -- DBOC from truncated set
        - 'DBOC_external': ndarray (n_channels, N_R) -- DBOC from excluded channels
        - 'R_grid': the input R grid
        If compute_exact_dPdR:
        - 'dPdR': ndarray (n_channels, n_channels, N_R) -- exact dP/dR
        - 'Q_exact': ndarray (n_channels, n_channels, N_R) -- Q + dP/dR
    """
    N_R = len(R_grid)
    V_coupling = solver._coupling_full
    casimir_all = np.concatenate(solver._channel_casimir)
    total_dim = len(casimir_all)
    n_ch = n_channels

    mu = np.zeros((n_ch, N_R))
    P = np.zeros((n_ch, n_ch, N_R))
    Q = np.zeros((n_ch, n_ch, N_R))
    DBOC_total = np.zeros((n_ch, N_R))
    DBOC_internal = np.zeros((n_ch, N_R))

    if compute_exact_dPdR:
        dPdR_arr = np.zeros((n_ch, n_ch, N_R))
        Q_exact = np.zeros((n_ch, n_ch, N_R))

    prev_vecs = None

    for i, R in enumerate(R_grid):
        # Full diagonalization to get ALL eigenvectors
        H = np.diag(casimir_all) + R * V_coupling
        evals_full, evecs_full = eigh(H)

        evals = evals_full[:n_ch]
        evecs = evecs_full[:, :n_ch].T  # (n_ch, dim)
        mu[:, i] = evals

        # Sign consistency
        if prev_vecs is not None:
            evecs = _enforce_sign_consistency(prev_vecs, evecs)
            # Also update evecs_full columns for consistency
            for c in range(n_ch):
                if np.dot(prev_vecs[c], evecs_full[:, c]) < 0:
                    evecs_full[:, c] *= -1
        prev_vecs = evecs

        # Compute P_μκ for ALL κ, for each truncated channel μ
        # P_full[mu_idx, kappa] = ⟨Φ_kappa|V|Φ_mu⟩ / (mu_mu - mu_kappa)
        P_full = np.zeros((n_ch, total_dim))
        for mu_idx in range(n_ch):
            V_phi = V_coupling @ evecs[mu_idx]
            for kappa in range(total_dim):
                if kappa == mu_idx:
                    continue
                gap = evals_full[mu_idx] - evals_full[kappa]
                if abs(gap) < 1e-10:
                    continue
                vec_k = evecs_full[:, kappa]
                P_full[mu_idx, kappa] = (vec_k @ V_phi) / gap

        # Extract truncated P matrix
        P[:, :, i] = P_full[:, :n_ch]

        # Q_μν via closure approximation over ALL channels
        for mu_idx in range(n_ch):
            for nu_idx in range(n_ch):
                q_val = 0.0
                for kappa in range(total_dim):
                    if kappa == mu_idx or kappa == nu_idx:
                        continue
                    gap_kn = evals_full[kappa] - evals_full[nu_idx]
                    if abs(gap_kn) < 1e-10:
                        continue
                    vec_nu = evecs[nu_idx] if nu_idx < n_ch else evecs_full[:, nu_idx]
                    vec_k = evecs_full[:, kappa]
                    p_kn = (vec_nu @ V_coupling @ vec_k) / gap_kn
                    q_val += P_full[mu_idx, kappa] * p_kn
                Q[mu_idx, nu_idx, i] = q_val

        # Compute exact dP/dR algebraically
        if compute_exact_dPdR:
            dPdR_i = _compute_dPdR(
                P_full, evals_full, V_coupling, evecs_full, n_ch, total_dim
            )
            dPdR_arr[:, :, i] = dPdR_i
            Q_exact[:, :, i] = Q[:, :, i] + dPdR_i

        # DBOC decomposition
        for mu_idx in range(n_ch):
            DBOC_total[mu_idx, i] = 0.5 * np.sum(P_full[mu_idx] ** 2)
            DBOC_internal[mu_idx, i] = 0.5 * np.sum(P[:, :, i][mu_idx] ** 2)

    DBOC_external = DBOC_total - DBOC_internal

    result = {
        'mu': mu,
        'P': P,
        'Q': Q,
        'DBOC_total': DBOC_total,
        'DBOC_internal': DBOC_internal,
        'DBOC_external': DBOC_external,
        'R_grid': R_grid,
    }

    if compute_exact_dPdR:
        result['dPdR'] = dPdR_arr
        result['Q_exact'] = Q_exact

    return result


def solve_hyperspherical_algebraic_coupled(
    Z: float = 2.0,
    n_basis: int = 15,
    l_max: int = 0,
    n_channels: int = 3,
    n_R: int = 200,
    R_min: float = 0.1,
    R_max: float = 30.0,
    N_R_radial: int = 3000,
    sigma: float = -3.0,
    q_mode: str = 'diagonal',
    verbose: bool = True,
    radial_method: str = 'fd',
    n_basis_radial: int = 25,
    alpha_radial: float = 2.0,
    matrix_method: str = 'quadrature',
) -> dict:
    """Full Level 3 coupled-channel solver using algebraic angular input.

    Combines the AlgebraicAngularSolver (exact Hellmann-Feynman P-matrix)
    with the existing coupled-channel radial solver.

    The coupled-channel equations are:
        [-½ d²/dR² δ_μν + V_eff_μ(R) δ_μν - P_μν(R) d/dR (- ½ Q_μν)] F_ν = E F_μ

    Parameters
    ----------
    Z : float
        Nuclear charge.
    n_basis : int
        Gegenbauer basis functions per l-channel.
    l_max : int
        Maximum partial wave.
    n_channels : int
        Number of coupled adiabatic channels.
    n_R : int
        R points for angular eigenvalue grid.
    R_min, R_max : float
        Hyperradial grid boundaries.
    N_R_radial : int
        Grid points for radial solver.
    sigma : float
        Shift-invert target energy.
    q_mode : str
        Treatment of second-derivative coupling:
        'none': P coupling only (no Q). Energy tends to overshoot below exact.
        'diagonal': DBOC_total on V_eff diagonal, no off-diagonal Q.
            Balances P coupling (lowers energy) with DBOC repulsion (raises).
        'full': Full closure Q = P·P over all channels. Tends to overcorrect.
        'internal': Q from truncated P only, plus external DBOC on diagonal.
        'exact': Exact Q = P·P + dP/dR over all channels. The dP/dR correction
            is computed algebraically from Hellmann-Feynman quantities, reducing
            the overcorrection from the closure approximation. The diagonal Q
            (DBOC) is unchanged because dP_μμ/dR = 0 by antisymmetry; only
            off-diagonal Q elements are corrected.
    verbose : bool
        Print progress.
    radial_method : str
        'fd' (default) for finite-difference or 'spectral' for Laguerre basis.
    n_basis_radial : int
        Number of Laguerre basis functions (spectral method only).
    alpha_radial : float
        Exponential decay parameter for Laguerre basis (spectral method only).
    matrix_method : str
        'quadrature' (default) or 'algebraic' for spectral method.
        When 'algebraic', overlap S and kinetic K use closed-form Laguerre
        recurrence; potential V, P-coupling, Q-coupling stay quadrature.

    Returns
    -------
    result : dict
    """
    import time

    t0 = time.time()

    # --- Step 1: Build algebraic angular solver ---
    solver = AlgebraicAngularSolver(Z, n_basis, l_max)

    # R grid: denser near origin
    R_grid = np.concatenate([
        np.linspace(R_min, 1.0, n_R // 3),
        np.linspace(1.0, 5.0, n_R // 3),
        np.linspace(5.0, R_max, n_R // 3 + 1),
    ])
    R_grid = np.unique(R_grid)

    if verbose:
        print(f"Algebraic coupled-channel: Z={Z}, n_basis={n_basis}, "
              f"l_max={l_max}, {n_channels} channels")
        print(f"Computing adiabatic curves + coupling on {len(R_grid)} R points...")

    # --- Step 2: Compute adiabatic curves, P, and Q matrices ---
    need_exact = (q_mode == 'exact')
    coupling = compute_algebraic_coupling(
        solver, R_grid, n_channels, compute_exact_dPdR=need_exact
    )
    mu_curves = coupling['mu']
    P = coupling['P']
    Q = coupling['Q']
    DBOC_total = coupling['DBOC_total']
    DBOC_internal = coupling['DBOC_internal']
    DBOC_external = coupling['DBOC_external']

    t1 = time.time()
    if verbose:
        print(f"  Angular + coupling solve: {t1 - t0:.2f}s")
        print(f"  DBOC total peak (ch 0):    {np.max(DBOC_total[0]):.6f} Ha")
        print(f"  DBOC internal peak (ch 0): {np.max(DBOC_internal[0]):.6f} Ha")
        print(f"  DBOC external peak (ch 0): {np.max(DBOC_external[0]):.6f} Ha")
        print(f"  Q_00 peak:                 {np.max(np.abs(Q[0, 0])):.6f}")
        if need_exact:
            dPdR = coupling['dPdR']
            Q_ex = coupling['Q_exact']
            print(f"  dP/dR_01 peak:             {np.max(np.abs(dPdR[0, 1])):.6f}")
            print(f"  Q_exact_01 peak:           {np.max(np.abs(Q_ex[0, 1])):.6f}")
            print(f"  Q_closure_01 peak:         {np.max(np.abs(Q[0, 1])):.6f}")

    # --- Step 3: Build splines based on q_mode ---
    use_Q = q_mode in ('full', 'internal', 'exact')
    add_dboc_to_diag = q_mode in ('diagonal',)

    V_eff_splines = []
    for ch in range(n_channels):
        V_eff_ch = effective_potential(R_grid, mu_curves[ch])
        if add_dboc_to_diag:
            V_eff_ch = V_eff_ch + DBOC_total[ch]
        elif q_mode == 'internal':
            # Add external DBOC to diagonal; internal handled by Q
            V_eff_ch = V_eff_ch + DBOC_external[ch]
        V_eff_splines.append(
            CubicSpline(R_grid, V_eff_ch, extrapolate=True)
        )

    P_splines = [[None] * n_channels for _ in range(n_channels)]
    Q_splines_obj = None
    for mu_idx in range(n_channels):
        for nu_idx in range(n_channels):
            P_splines[mu_idx][nu_idx] = CubicSpline(
                R_grid, P[mu_idx, nu_idx], extrapolate=True
            )

    if use_Q:
        Q_splines_obj = [[None] * n_channels for _ in range(n_channels)]
        if q_mode == 'exact':
            Q_data = coupling['Q_exact']
        elif q_mode == 'full':
            Q_data = Q
        else:
            Q_data = np.zeros_like(Q)
        if q_mode == 'internal':
            # Q_internal = P_trunc @ P_trunc at each R point
            for i in range(len(R_grid)):
                P_i = P[:, :, i]
                Q_data[:, :, i] = P_i @ P_i
        for mu_idx in range(n_channels):
            for nu_idx in range(n_channels):
                Q_splines_obj[mu_idx][nu_idx] = CubicSpline(
                    R_grid, Q_data[mu_idx, nu_idx], extrapolate=True
                )

    # --- Step 4: Solve coupled radial equation ---
    if radial_method == 'spectral':
        if verbose:
            print(f"Solving coupled-channel radial equation "
                  f"(spectral, n_basis={n_basis_radial}, {n_channels} ch, "
                  f"q_mode={q_mode})...")
        E_all, F_all, R_grid_rad = solve_coupled_radial_spectral(
            V_eff_splines, P_splines,
            Q_splines=Q_splines_obj,
            n_channels=n_channels,
            n_basis=n_basis_radial, alpha=alpha_radial,
            R_min=R_min,
            n_states=min(5, n_channels * 3),
            include_Q=use_Q,
            matrix_method=matrix_method,
        )
    else:
        if verbose:
            print(f"Solving coupled-channel radial equation "
                  f"(N_R={N_R_radial}, {n_channels} ch, q_mode={q_mode})...")
        E_all, F_all, R_grid_rad = solve_coupled_radial(
            V_eff_splines, P_splines,
            Q_splines=Q_splines_obj,
            n_channels=n_channels,
            R_min=R_min, R_max=R_max, N_R=N_R_radial,
            n_states=min(5, n_channels * 3),
            sigma=sigma,
            include_Q=use_Q,
        )

    t2 = time.time()
    E_exact = -2.903724

    # Channel weights for ground state
    h_rad = R_grid_rad[1] - R_grid_rad[0]
    weights = np.array([
        h_rad * np.sum(F_all[0, ch] ** 2)
        for ch in range(n_channels)
    ])
    if weights.sum() > 0:
        weights /= weights.sum()

    if verbose:
        print(f"  Radial solve: {t2 - t1:.2f}s")
        print(f"\n  Ground state energy: {E_all[0]:.6f} Ha")
        print(f"  Exact He:            {E_exact:.6f} Ha")
        err = abs(E_all[0] - E_exact) / abs(E_exact) * 100
        print(f"  Error:               {err:.4f}%")
        print(f"  Channel weights:     {weights}")
        print(f"  Total time:          {t2 - t0:.2f}s")

    # Also compute single-channel adiabatic energy for comparison
    V_eff_single = effective_potential(R_grid, mu_curves[0])
    V_eff_spline_single = CubicSpline(R_grid, V_eff_single, extrapolate=True)
    if radial_method == 'spectral':
        E_single, _, _ = solve_radial_spectral(
            V_eff_spline_single, n_basis=n_basis_radial,
            alpha=alpha_radial, R_min=R_min, n_states=1,
            matrix_method=matrix_method,
        )
    else:
        E_single, _, _ = solve_radial(
            V_eff_spline_single, R_min, R_max, N_R_radial, n_states=1
        )

    if verbose:
        err_s = abs(E_single[0] - E_exact) / abs(E_exact) * 100
        print(f"\n  Single-channel energy: {E_single[0]:.6f} Ha ({err_s:.4f}%)")
        direction = "UP (correct)" if E_all[0] > E_single[0] else "DOWN"
        print(f"  Coupled vs single: {E_all[0] - E_single[0]:+.6f} Ha ({direction})")

    result_dict = {
        'energy': E_all[0],
        'energy_single_channel': E_single[0],
        'energies': E_all,
        'R_grid_angular': R_grid,
        'mu_curves': mu_curves,
        'V_eff': V_eff_single,
        'DBOC_total': DBOC_total,
        'DBOC_internal': DBOC_internal,
        'DBOC_external': DBOC_external,
        'P_coupling': P,
        'Q_coupling': Q,
        'R_grid_radial': R_grid_rad,
        'wavefunction': F_all[0],
        'channel_weights': weights,
        'solver': solver,
        'Z': Z,
        'l_max': l_max,
        'n_channels': n_channels,
        'n_basis': n_basis,
        'coupled': True,
        'q_mode': q_mode,
    }

    if need_exact:
        result_dict['dPdR'] = coupling['dPdR']
        result_dict['Q_exact'] = coupling['Q_exact']

    return result_dict
