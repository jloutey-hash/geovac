"""Debug: prolate spheroidal beta finder — corrected version.

Key insight: Both angular and radial equations share the same separation
constant A_sep. We use:
  - Angular: matrix Legendre expansion (fast, exact)
  - Radial: self-adjoint FD with SMALL xi_max (captures bound states,
    physical eigenvalues clearly separated from artificial ones)

The Coulomb angular eigenvalue for state (m, n_sph) corresponds to the
(N-1-n_sph)-th eigenvalue of H (counting from bottom), i.e., from the
TOP of the spectrum. This is because the unperturbed (c=0) eigenvalues
are -r(r+1), with the physical ground state at r=m having the LARGEST
(least negative) eigenvalue.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
from scipy.special import obl_cv
from math import factorial
import time


def angular_sep_const(m_abs: int, n_sph: int, c: float,
                      b: float = 0.0, n_basis: int = 40) -> float:
    """Coulomb angular separation constant A_sep.

    This is the eigenvalue of H = L_0 + c^2*eta^2 + b*eta corresponding
    to the state with n_sph angular nodes. Counted from the TOP of the
    spectrum (least negative).

    For b=0: A_sep = -obl_cv(m, m+n_sph, c).
    """
    N = n_basis
    r_vals = np.arange(m_abs, m_abs + N, dtype=float)
    norms = np.array([2.0/(2*r+1) * factorial(int(r+m_abs))/factorial(int(r-m_abs))
                       for r in r_vals])

    nu_mat = np.zeros((N, N))
    for i in range(N - 1):
        r = r_vals[i]
        val = (r - m_abs + 1) * np.sqrt(norms[i+1]) / ((2*r+1) * np.sqrt(norms[i]))
        nu_mat[i+1, i] = val
        nu_mat[i, i+1] = val

    H = np.diag(-r_vals * (r_vals + 1))
    H += c**2 * (nu_mat @ nu_mat)
    H += b * nu_mat

    evals = np.sort(np.linalg.eigvalsh(H))
    # Physical states counted from top: n_sph=0 -> evals[-1], n_sph=1 -> evals[-2]
    return evals[N - 1 - n_sph]


def radial_sep_consts(m_abs: int, c: float, a: float,
                      n_grid: int = 500, xi_max: float = 6.0,
                      n_evals: int = 10) -> np.ndarray:
    """Radial separation constants via self-adjoint FD.

    The radial equation in self-adjoint form:
      d/dxi[(xi^2-1)*dR/dxi] + [-c^2*xi^2 + a*xi - m^2/(xi^2-1)]*R = A_sep*R

    Dirichlet BCs at both ends. Physical eigenvalues are the LARGEST
    (from the top), separated from artificial deep-well states.
    """
    N = n_grid
    xi_min = 1.0 + 0.002
    h = (xi_max - xi_min) / (N + 1)
    xi = xi_min + (np.arange(N) + 1) * h

    # Self-adjoint form: d/dxi[p*dR/dxi] + q*R = A*R
    # p(xi) = xi^2-1, q(xi) = -c^2*xi^2 + a*xi - m^2/(xi^2-1)
    xi2_1 = xi**2 - 1
    q = -c**2 * xi**2 + a * xi - m_abs**2 / xi2_1

    # Midpoint p values
    xi_half_p = xi + h/2
    xi_half_m = xi - h/2
    p_plus = xi_half_p**2 - 1
    p_minus = xi_half_m**2 - 1

    # Tridiagonal matrix
    diag = -(p_plus + p_minus) / h**2 + q
    off_diag = p_plus[:-1] / h**2

    # Get eigenvalues from the TOP (most physical)
    try:
        evals = eigh_tridiagonal(diag, off_diag, eigvals_only=True,
                                  select='i', select_range=(N - n_evals, N - 1))
    except Exception:
        evals = eigh_tridiagonal(diag, off_diag, eigvals_only=True)
        evals = np.sort(evals)

    # Return sorted descending (most positive first = n_rad=0)
    return np.sort(evals)[::-1][:n_evals]


def find_beta(m_abs: int, n_sph: int, n_rad: int,
              p0: float, R: float, Z_A: float, Z_B: float,
              beta_min: float = 0.05, beta_max: float = 5.0,
              n_scan: int = 50, verbose: bool = True) -> list:
    """Find beta where angular and radial separation constants match."""
    c = p0 * R / 2
    Z_total = Z_A + Z_B
    Z_diff = Z_B - Z_A

    def mismatch(beta):
        a = beta * Z_total * R
        b = beta * Z_diff * R
        A_ang = angular_sep_const(m_abs, n_sph, c, b)
        A_rads = radial_sep_consts(m_abs, c, a, n_evals=n_rad + 5)
        if len(A_rads) <= n_rad:
            return float('nan')
        return A_ang - A_rads[n_rad]

    betas = np.linspace(beta_min, beta_max, n_scan)
    vals = [mismatch(bt) for bt in betas]
    vals = np.array(vals)

    if verbose:
        for bt, v in zip(betas, vals):
            if not np.isnan(v):
                print(f"  beta={bt:.3f}: mismatch={v:+.6f}")

    roots = []
    for i in range(len(vals) - 1):
        if not np.isnan(vals[i]) and not np.isnan(vals[i+1]):
            if vals[i] * vals[i+1] < 0:
                try:
                    root = brentq(mismatch, betas[i], betas[i+1], xtol=1e-10)
                    roots.append(root)
                    if verbose:
                        print(f"  ** Root at beta = {root:.8f}")
                except Exception as e:
                    if verbose:
                        print(f"  brentq failed: {e}")
    return roots


if __name__ == "__main__":
    # === Step 1: Verify angular function ===
    print("=== Angular separation constant verification ===")
    for m in [0, 1]:
        for n_sph in range(3):
            n_total = m + n_sph
            A_sep = angular_sep_const(m, n_sph, 2.0, 0.0)
            A_obl = obl_cv(m, n_total, 2.0)
            expected = -A_obl
            diff = abs(A_sep - expected)
            print(f"  m={m}, n_sph={n_sph}: A_sep={A_sep:+.6f}, "
                  f"-obl_cv={expected:+.6f}, diff={diff:.2e}")

    # === Step 2: Radial eigenvalues for H2+ ===
    print("\n=== Radial eigenvalues: H2+ at equilibrium ===")
    E_h2p = -1.1026
    p0_h2p = np.sqrt(-2*E_h2p)
    c_h2p = p0_h2p * 2 / 2  # R=2
    a_h2p = 1.0 * 2.0 * 2   # beta=1, Z_tot=2, R=2
    print(f"  E={E_h2p}, p0={p0_h2p:.4f}, c={c_h2p:.4f}, a={a_h2p:.4f}")

    t0 = time.time()
    A_rads = radial_sep_consts(0, c_h2p, a_h2p, n_evals=5)
    print(f"  Time: {time.time()-t0:.3f}s")
    for i, A in enumerate(A_rads):
        print(f"  A_rad[{i}] (from top) = {A:.6f}")

    A_ang = angular_sep_const(0, 0, c_h2p, 0.0)
    print(f"  A_ang(m=0, n_sph=0) = {A_ang:.6f}")
    if len(A_rads) > 0:
        print(f"  Mismatch(beta=1): {A_ang - A_rads[0]:.6f}")

    # === Step 3: LiH ===
    print("\n=== LiH eigenvalues at beta=1 ===")
    Z_A, Z_B = 3.0, 1.0
    R = 3.015
    p0 = np.sqrt(10.0)
    c = p0 * R / 2
    a_1 = 1.0 * 4.0 * R
    b_1 = 1.0 * (Z_B - Z_A) * R
    print(f"  c={c:.4f}, a={a_1:.4f}, b={b_1:.4f}")

    A_rads = radial_sep_consts(0, c, a_1, n_evals=5)
    for i, A in enumerate(A_rads):
        print(f"  A_rad[{i}] (from top) = {A:.6f}")

    A_ang = angular_sep_const(0, 0, c, b_1)
    print(f"  A_ang(m=0, n_sph=0) = {A_ang:.6f}")
    if len(A_rads) > 0:
        print(f"  Mismatch: {A_ang - A_rads[0]:.6f}")

    # === Step 4: Find betas ===
    print("\n=== Beta search: LiH 1sigma (m=0, n_sph=0, n_rad=0) ===")
    t0 = time.time()
    roots = find_beta(0, 0, 0, p0, R, Z_A, Z_B, n_scan=40, beta_max=3.0)
    dt = time.time() - t0
    print(f"  Search time: {dt:.1f}s")

    if roots:
        for root in roots:
            eps_H = p0**2/2 - root * Z_B * p0
            eps_Li = p0**2/2 - root * Z_A * p0
            print(f"\n  beta = {root:.8f}")
            print(f"  eps(H 1s) = {eps_H:.6f} Ha  (bound = {eps_H < 0})")
            print(f"  eps(Li 1s) = {eps_Li:.6f} Ha")
    else:
        print("\n  No roots found — printing diagnostics")
        for beta_test in [0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]:
            a_p = beta_test * 4.0 * R
            b_p = beta_test * (Z_B - Z_A) * R
            A_a = angular_sep_const(0, 0, c, b_p)
            A_r = radial_sep_consts(0, c, a_p, n_evals=3)
            A_r0 = A_r[0] if len(A_r) > 0 else float('nan')
            print(f"    beta={beta_test:.1f}: A_ang={A_a:+.4f}, "
                  f"A_rad[0]={A_r0:+.4f}, mismatch={A_a - A_r0:+.4f}")
