"""
Level 4 natural geometry: σ-only single-channel solver for H₂.

Phase 1 implementation of the two-center, two-electron fiber bundle.
Uses molecule-frame hyperspherical coordinates (R_e, α, θ₁, θ₂) with
the l₁ = l₂ = 0 projection (σ-only, single adiabatic channel).

The angular eigenvalue problem at fixed ρ = R/(2R_e):
    [-1/2 u'' + V_eff(α; ρ) u] = μ u
with Dirichlet BCs u(0) = u(π/2) = 0.

The molecular charge function C_mol^{00}(α; ρ) is computed by averaging
the full two-center potential over θ₁, θ₂ using Gauss-Legendre quadrature.

References:
  - Level 4 design document: papers/core/level4_geometry_design.md
  - Lin, Phys. Rep. 257, 1 (1995) — hyperspherical review
  - Macek, J. Phys. B 1, 831 (1968) — hyperspherical two-electron atoms
"""

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.interpolate import CubicSpline
from typing import Tuple, Optional


def _nuclear_attraction_averaged(
    alpha: np.ndarray,
    rho: float,
    Z: float,
    n_quad: int = 32,
) -> np.ndarray:
    """
    Compute the l₁=l₂=0 (σ-only) projection of the nuclear attraction.

    For each electron i at distance r_i from the molecular midpoint and
    angle θ_i from the internuclear axis:

        -Z/r_iA = -Z / sqrt(r_i² + R²/4 + r_i R cos θ_i)
        -Z/r_iB = -Z / sqrt(r_i² + R²/4 - r_i R cos θ_i)

    With r₁ = R_e cos α, r₂ = R_e sin α, and ρ = R/(2R_e):

        -Z/r_1A = -(Z/R_e) / (cos α · sqrt(1 + (ρ/cos α)² + 2(ρ/cos α) cos θ₁))

    The σ-projection averages over θ₁, θ₂ independently with weight sin θ / 2
    (the l=0 spherical harmonic squared |Y_0^0|² = 1/(4π), normalized).

    Parameters
    ----------
    alpha : ndarray of shape (n_alpha,)
        Correlation angle grid points.
    rho : float
        Dimensionless parameter R / (2 R_e).
    Z : float
        Nuclear charge.
    n_quad : int
        Number of Gauss-Legendre quadrature points for θ integration.

    Returns
    -------
    V_nuc : ndarray of shape (n_alpha,)
        Nuclear attraction contribution to the charge function (divided by R_e).
        Multiply by R_e to get the actual potential.
    """
    # Gauss-Legendre quadrature on cos θ ∈ [-1, 1]
    x_quad, w_quad = np.polynomial.legendre.leggauss(n_quad)

    n_alpha = len(alpha)
    V_nuc = np.zeros(n_alpha)

    cos_a = np.cos(alpha)
    sin_a = np.sin(alpha)

    for ia in range(n_alpha):
        # Electron 1: r₁ = R_e cos α, so r₁/R_e = cos α
        # Distance to nucleus A: r₁A/R_e = sqrt(cos²α + ρ² + 2ρ cos α cos θ₁)
        # Distance to nucleus B: r₁B/R_e = sqrt(cos²α + ρ² - 2ρ cos α cos θ₁)
        ca = cos_a[ia]
        sa = sin_a[ia]

        # Average over θ₁ for electron 1
        # -Z/(R_e · r₁A/R_e) averaged over cos θ₁
        v1 = 0.0
        for iq in range(n_quad):
            ct = x_quad[iq]
            wt = w_quad[iq]
            d1A = np.sqrt(ca**2 + rho**2 + 2 * rho * ca * ct)
            d1B = np.sqrt(ca**2 + rho**2 - 2 * rho * ca * ct)
            # Factor 1/2 from normalizing ∫₋₁¹ (sin θ dθ → d(cos θ)), weight = 1/2
            v1 += wt * (-Z / d1A - Z / d1B) * 0.5

        # Average over θ₂ for electron 2 (same formula with sin α)
        v2 = 0.0
        for iq in range(n_quad):
            ct = x_quad[iq]
            wt = w_quad[iq]
            d2A = np.sqrt(sa**2 + rho**2 + 2 * rho * sa * ct)
            d2B = np.sqrt(sa**2 + rho**2 - 2 * rho * sa * ct)
            v2 += wt * (-Z / d2A - Z / d2B) * 0.5

        V_nuc[ia] = v1 + v2

    return V_nuc


def _ee_repulsion_averaged(
    alpha: np.ndarray,
    n_quad: int = 32,
) -> np.ndarray:
    """
    Compute the l₁=l₂=0 projection of the e-e repulsion.

    1/r₁₂ = 1 / (R_e sqrt(1 - sin 2α cos θ₁₂))

    where cos θ₁₂ = cos θ₁ cos θ₂ + sin θ₁ sin θ₂ cos Φ.

    The σ-projection averages over θ₁, θ₂, Φ with the measure
    (sin θ₁ / 2)(sin θ₂ / 2)(1 / 2π).

    Returns V_ee / R_e (dimensionless charge function contribution).

    Parameters
    ----------
    alpha : ndarray of shape (n_alpha,)
        Correlation angle grid points.
    n_quad : int
        Quadrature points per angular dimension.

    Returns
    -------
    C_ee : ndarray of shape (n_alpha,)
        E-e repulsion charge function (divided by R_e).
    """
    # For l₁=l₂=0, the θ₁₂-dependent 1/r₁₂ can be averaged analytically.
    # Using the multipole expansion:
    #   1/r₁₂ = (1/R_e) Σ_k (r_</r_>)^k P_k(cos θ₁₂) / r_>
    # For k=0 only (l₁=l₂=0 projection):
    #   ⟨1/r₁₂⟩_{l=0} = (1/R_e) · 1/max(r₁, r₂)
    #                   = (1/R_e) · 1/max(cos α, sin α)
    # This is because ∫ P_k(cos θ₁₂) dΩ₁ dΩ₂ = 0 for k > 0
    # when averaged over all angles (Y₀⁰ is isotropic).

    cos_a = np.cos(alpha)
    sin_a = np.sin(alpha)
    max_sc = np.maximum(cos_a, sin_a)

    return 1.0 / max_sc


def compute_Cmol_00(
    alpha_grid: np.ndarray,
    rho: float,
    Z: float = 1.0,
    n_quad: int = 32,
) -> np.ndarray:
    """
    Compute the sigma-only molecular charge function C_mol^{00}(alpha; rho).

    This is the l1=l2=0 projection of the full molecular potential, divided
    by R_e. Both C_ee and V_nuc are returned as charge functions (/ R_e).
    In the angular eigenvalue problem, the potential is:

        V(alpha; R_e, R) = -2 + R_e * C_mol^{00}(alpha; rho)

    Parameters
    ----------
    alpha_grid : ndarray
        Grid of α values in (0, π/2).
    rho : float
        Dimensionless parameter R / (2R_e).
    Z : float
        Nuclear charge (1.0 for H₂).
    n_quad : int
        Quadrature order for angular averaging.

    Returns
    -------
    C_mol : ndarray
        Molecular charge function (dimensionless, multiply by R_e for potential).
    """
    C_ee = _ee_repulsion_averaged(alpha_grid, n_quad)
    V_nuc = _nuclear_attraction_averaged(alpha_grid, rho, Z, n_quad)

    # V_nuc is already the full nuclear potential / R_e
    # C_ee is the e-e charge function / R_e
    # Total charge function: multiply by R_e in the angular equation
    return C_ee + V_nuc


def solve_angular(
    rho: float,
    R_e: float,
    Z: float = 1.0,
    n_alpha: int = 200,
    n_quad: int = 32,
    n_channels: int = 1,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Solve the angular eigenvalue problem at fixed ρ and R_e.

    [-1/2 u'' + V(α; R_e, ρ)] u = μ u

    where V(α) = -2 + R_e · C_mol^{00}(α; ρ) and BCs are u(0) = u(π/2) = 0.

    Parameters
    ----------
    rho : float
        Dimensionless parameter R / (2R_e).
    R_e : float
        Electronic hyperradius (bohr).
    Z : float
        Nuclear charge per nucleus.
    n_alpha : int
        Number of interior FD grid points.
    n_quad : int
        Quadrature order for angular averaging.
    n_channels : int
        Number of eigenvalues to return.

    Returns
    -------
    mu : ndarray of shape (n_channels,)
        Angular eigenvalues.
    vecs : ndarray of shape (n_channels, n_alpha)
        Eigenvectors (u functions on the α grid).
    alpha_grid : ndarray of shape (n_alpha,)
        Interior α grid points.
    """
    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h  # interior points in (0, π/2)

    # Charge function
    C_mol = compute_Cmol_00(alpha, rho, Z, n_quad)

    # Effective potential: Liouville shift + R_e * charge function
    # For l₁ = l₂ = 0: no centrifugal barriers
    V_eff = -2.0 + R_e * C_mol

    # Tridiagonal FD Hamiltonian: -1/2 u'' + V u = μ u
    diag = np.ones(n_alpha) / h**2 + V_eff
    off_diag = -0.5 * np.ones(n_alpha - 1) / h**2

    evals, evecs = eigh_tridiagonal(
        diag, off_diag,
        select='i', select_range=(0, n_channels - 1),
    )

    # Normalize
    for i in range(n_channels):
        norm = np.sqrt(h * np.sum(evecs[:, i]**2))
        if norm > 0:
            evecs[:, i] /= norm

    return evals, evecs.T, alpha


def compute_adiabatic_curve(
    R: float,
    R_e_grid: np.ndarray,
    Z: float = 1.0,
    n_alpha: int = 200,
    n_quad: int = 32,
) -> np.ndarray:
    """
    Compute the effective adiabatic potential U(R_e) at fixed R.

    U(R_e; R) = [μ(R_e; R) + 15/8] / R_e²

    where μ is the angular eigenvalue at ρ = R/(2R_e).

    Nuclear repulsion V_NN = Z²/R is NOT included (added separately).

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    R_e_grid : ndarray
        Electronic hyperradius grid points (bohr).
    Z : float
        Nuclear charge per nucleus.
    n_alpha : int
        FD grid points for angular solve.
    n_quad : int
        Quadrature order.

    Returns
    -------
    U : ndarray
        Effective potential on R_e_grid (Ha).
    """
    n_Re = len(R_e_grid)
    mu_vals = np.zeros(n_Re)

    for i, R_e in enumerate(R_e_grid):
        rho = R / (2.0 * R_e)
        evals, _, _ = solve_angular(rho, R_e, Z, n_alpha, n_quad)
        mu_vals[i] = evals[0]

    U = (mu_vals + 15.0 / 8.0) / R_e_grid**2
    return U


def solve_level4_h2(
    R: float,
    Z: float = 1.0,
    n_alpha: int = 200,
    n_Re: int = 400,
    n_quad: int = 32,
    R_e_min: float = 0.3,
    R_e_max: float = 15.0,
    verbose: bool = True,
) -> dict:
    """
    Full Level 4 σ-only solver for H₂ at fixed internuclear distance R.

    1. Compute adiabatic potential U(R_e) by sweeping R_e
    2. Solve hyperradial equation for E_elec(R)
    3. Add nuclear repulsion: E_total = E_elec + Z²/R

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    Z : float
        Nuclear charge per nucleus (1.0 for H₂).
    n_alpha : int
        FD grid points for angular α equation.
    n_Re : int
        Grid points for hyperradial R_e equation.
    n_quad : int
        Quadrature order for angular averaging.
    R_e_min : float
        Minimum electronic hyperradius (bohr).
    R_e_max : float
        Maximum electronic hyperradius (bohr).
    verbose : bool
        Print diagnostic output.

    Returns
    -------
    result : dict
        Keys: 'E_elec', 'E_total', 'D_e', 'D_e_pct', 'R', 'Z',
        'R_e_grid_angular', 'U_adiabatic', 'R_e_grid_radial', 'wavefunction'.
    """
    import time
    t0 = time.time()

    # --- Step 1: Compute adiabatic curve U(R_e) ---
    # Use non-uniform grid: denser near small R_e where potential varies fast
    R_e_angular = np.concatenate([
        np.linspace(R_e_min, 1.0, 40),
        np.linspace(1.0, 3.0, 40),
        np.linspace(3.0, 6.0, 30),
        np.linspace(6.0, R_e_max, 20),
    ])
    R_e_angular = np.unique(R_e_angular)

    if verbose:
        print(f"Level 4 sigma-only solver for H2 at R = {R:.4f} bohr")
        print(f"  n_alpha={n_alpha}, n_Re={n_Re}, n_quad={n_quad}")
        print(f"  Computing adiabatic curve on {len(R_e_angular)} R_e points...")

    U_angular = compute_adiabatic_curve(R, R_e_angular, Z, n_alpha, n_quad)

    t1 = time.time()
    if verbose:
        print(f"  Angular sweep: {t1 - t0:.2f}s")
        i_min = np.argmin(U_angular)
        print(f"  U_min = {U_angular[i_min]:.6f} Ha at R_e = {R_e_angular[i_min]:.3f}")

    # --- Step 2: Interpolate and solve radial equation ---
    U_spline = CubicSpline(R_e_angular, U_angular, extrapolate=True)

    # Radial grid for FD solve
    h_Re = (R_e_max - R_e_min) / (n_Re + 1)
    R_e_radial = R_e_min + (np.arange(n_Re) + 1) * h_Re

    V_radial = U_spline(R_e_radial)

    # Tridiagonal FD: -1/2 F'' + U(R_e) F = E F
    diag = np.ones(n_Re) / h_Re**2 + V_radial
    off_diag = -0.5 * np.ones(n_Re - 1) / h_Re**2

    evals, evecs = eigh_tridiagonal(
        diag, off_diag,
        select='i', select_range=(0, 0),
    )

    E_elec = evals[0]
    F = evecs[:, 0]

    # Normalize
    norm = np.sqrt(h_Re * np.sum(F**2))
    if norm > 0:
        F /= norm

    t2 = time.time()

    # --- Step 3: Total energy and D_e ---
    V_NN = Z**2 / R
    E_total = E_elec + V_NN

    # Dissociation energy: D_e = E(H) + E(H) - E(H₂)
    # E(H) = -0.5 Ha for each hydrogen atom
    E_atoms = -2 * Z**2 / 2.0  # = -1.0 Ha for H₂
    D_e = E_atoms - E_total  # positive if bound

    # Exact values for comparison
    E_exact = -1.17447  # Ha (Kolos & Wolniewicz)
    D_e_exact = 0.17447  # Ha
    D_e_pct = D_e / D_e_exact * 100 if D_e_exact > 0 else 0.0

    if verbose:
        print(f"  Radial solve: {t2 - t1:.2f}s")
        print(f"\n  === Results at R = {R:.4f} bohr ===")
        print(f"  E_elec     = {E_elec:.6f} Ha")
        print(f"  V_NN       = {V_NN:.6f} Ha")
        print(f"  E_total    = {E_total:.6f} Ha  (exact: {E_exact:.6f})")
        print(f"  E_atoms    = {E_atoms:.6f} Ha")
        print(f"  D_e        = {D_e:.6f} Ha  (exact: {D_e_exact:.6f})")
        print(f"  D_e / D_e_exact = {D_e_pct:.1f}%")
        if D_e_pct > 92.4:
            print(f"  ** IMPROVES on Paper 12 Neumann V_ee (92.4%) **")
        else:
            print(f"  Paper 12 Neumann V_ee achieved 92.4%")
        print(f"  Total time: {t2 - t0:.2f}s")

    return {
        'E_elec': E_elec,
        'E_total': E_total,
        'D_e': D_e,
        'D_e_pct': D_e_pct,
        'R': R,
        'Z': Z,
        'R_e_grid_angular': R_e_angular,
        'U_adiabatic': U_angular,
        'R_e_grid_radial': R_e_radial,
        'wavefunction': F,
    }
