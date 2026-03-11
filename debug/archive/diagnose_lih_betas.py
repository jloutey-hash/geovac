"""Debug: find LiH molecular Sturmian betas using prolate spheroidal separation.

Now with:
1. Angular eigenvalue verified against obl_cv (exact)
2. Radial eigenvalue with Neumann BC at xi=1 (correct for m=0)
3. Matching condition: A_ang + L_top = 0
"""
import numpy as np
from scipy.optimize import brentq
from scipy.linalg import eigh_tridiagonal
from math import factorial
import time


def angular_sep(m_abs, n_sph, c, b=0.0, n_basis=50):
    """A_sep from angular equation (top of spectrum)."""
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
    H = np.diag(-r_vals * (r_vals + 1)) + c**2 * (nu_mat @ nu_mat) + b * nu_mat
    evals = np.sort(np.linalg.eigvalsh(H))
    return evals[N - 1 - n_sph]


def radial_top_eval(m_abs, c, a, n_grid=1200, xi_max=8.0, n_top=5):
    """Top eigenvalues of L_xi (Neumann BC at left for m=0)."""
    N = n_grid
    xi_min = 1.0 + 0.0005
    h = (xi_max - xi_min) / (N + 1)
    xi = xi_min + (np.arange(N) + 1) * h

    xi2_1 = xi**2 - 1
    q = -c**2 * xi**2 + a * xi - m_abs**2 / xi2_1

    p_plus = (xi + h/2)**2 - 1
    p_minus = (xi - h/2)**2 - 1

    diag = -(p_plus + p_minus) / h**2 + q
    off = p_plus[:-1] / h**2

    if m_abs == 0:
        diag[0] = -p_plus[0] / h**2 + q[0]  # Neumann BC

    evals = eigh_tridiagonal(diag, off, eigvals_only=True)
    return np.sort(evals)[::-1][:n_top]


def find_beta(m_abs, n_sph, n_rad, p0, R, Z_A, Z_B,
              beta_min=0.05, beta_max=5.0, n_scan=60, verbose=True):
    """Find beta where A_ang + L_top[n_rad] = 0."""
    c = p0 * R / 2
    Z_total = Z_A + Z_B
    Z_diff = Z_B - Z_A

    def mismatch(beta):
        a = beta * Z_total * R
        b = beta * Z_diff * R
        A_a = angular_sep(m_abs, n_sph, c, b)
        L_ev = radial_top_eval(m_abs, c, a, n_top=n_rad + 3)
        if len(L_ev) <= n_rad:
            return float('nan')
        return A_a + L_ev[n_rad]  # matching: A_ang = -L_top

    betas = np.linspace(beta_min, beta_max, n_scan)
    vals = [mismatch(bt) for bt in betas]
    vals = np.array(vals)

    if verbose:
        for bt, v in zip(betas, vals):
            if not np.isnan(v):
                print(f"  beta={bt:.3f}: A_ang + L_top = {v:+.4f}")

    roots = []
    for i in range(len(vals) - 1):
        if not np.isnan(vals[i]) and not np.isnan(vals[i+1]):
            if vals[i] * vals[i+1] < 0:
                try:
                    root = brentq(mismatch, betas[i], betas[i+1], xtol=1e-8)
                    roots.append(root)
                    if verbose:
                        print(f"  ** ROOT at beta = {root:.8f}")
                except Exception as e:
                    if verbose:
                        print(f"  brentq failed: {e}")
    return roots


if __name__ == "__main__":
    Z_A, Z_B = 3.0, 1.0
    R = 3.015
    p0 = np.sqrt(10.0)
    c = p0 * R / 2

    print("=" * 60)
    print(f"LiH betas: Z_A={Z_A}, Z_B={Z_B}, R={R}, p0={p0:.4f}")
    print(f"c = p0*R/2 = {c:.4f}")
    print("=" * 60)

    # Find beta for 1sigma ground state (m=0, n_sph=0, n_rad=0)
    print("\n--- 1sigma: m=0, n_sph=0, n_rad=0 ---")
    t0 = time.time()
    roots = find_beta(0, 0, 0, p0, R, Z_A, Z_B, beta_max=3.0, n_scan=50)
    dt = time.time() - t0
    print(f"  Time: {dt:.1f}s")

    if roots:
        beta = roots[0]
        eps_H = p0**2/2 - beta * Z_B * p0
        eps_Li = p0**2/2 - beta * Z_A * p0
        print(f"\n  beta = {beta:.8f}")
        print(f"  Sturmian diag (H, Z=1):  eps = {eps_H:.6f} Ha "
              f"({'BOUND' if eps_H < 0 else 'UNBOUND'})")
        print(f"  Sturmian diag (Li, Z=3): eps = {eps_Li:.6f} Ha")
    else:
        print("\n  No roots found for 1sigma!")

    # Find betas for other states
    for label, m, ns, nr in [
        ("2sigma", 0, 1, 0),
        ("3sigma", 0, 0, 1),
        ("1pi", 1, 0, 0),
    ]:
        print(f"\n--- {label}: m={m}, n_sph={ns}, n_rad={nr} ---")
        t0 = time.time()
        roots = find_beta(m, ns, nr, p0, R, Z_A, Z_B, beta_max=3.0, n_scan=50)
        dt = time.time() - t0
        print(f"  Time: {dt:.1f}s")
        if roots:
            beta = roots[0]
            eps_H = p0**2/2 - beta * Z_B * p0
            print(f"  beta = {beta:.8f}, eps(H) = {eps_H:.6f} Ha")
