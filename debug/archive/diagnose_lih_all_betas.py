"""Compute ALL LiH molecular betas up to nmax=3 and write results.

Molecular orbital quantum numbers (m, n_sph, n_rad):
  m: azimuthal quantum number
  n_sph: angular nodes in eta direction
  n_rad: radial nodes in xi direction
  n_total = m + n_sph (total angular, used in basis indexing)
  n_principal = 1 + n_total + n_rad (total principal quantum number)

For nmax=3: n_principal ≤ 3, so:
  n=1: (0,0,0) -> 1 state
  n=2: (0,1,0), (0,0,1), (1,0,0) -> 3 states
  n=3: (0,2,0), (0,1,1), (0,0,2), (1,1,0), (1,0,1), (2,0,0) -> 6 states
Total: 10 molecular orbitals (matching n^2 = 1+4+9... wait, that's 10 not 14)
Actually for n=1: 1, n=2: 4, n=3: 9 -> but these count m and -m separately.

Let me be more careful: for each (m, n_sph, n_rad), the degeneracy is:
  m=0: 1 (just m=0)
  m>0: 2 (m and -m)

For nmax=3:
  n=1: (0,0,0) -> 1 orbital -> 1 state
  n=2: (0,1,0), (0,0,1), (1,0,0)[x2] -> 4 states
  n=3: (0,2,0), (0,1,1), (0,0,2), (1,1,0)[x2], (1,0,1)[x2], (2,0,0)[x2]
       -> 3 + 3*2 = 9 states
Total: 14 states (= 1+4+9)
"""
import numpy as np
from scipy.optimize import brentq
from scipy.linalg import eigh_tridiagonal
from math import factorial
import time


def angular_sep(m_abs, n_sph, c, b=0.0, n_basis=50):
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
        diag[0] = -p_plus[0] / h**2 + q[0]
    evals = eigh_tridiagonal(diag, off, eigvals_only=True)
    return np.sort(evals)[::-1][:n_top]


def find_beta(m_abs, n_sph, n_rad, p0, R, Z_A, Z_B,
              beta_min=0.05, beta_max=5.0, n_scan=60):
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
        return A_a + L_ev[n_rad]

    betas = np.linspace(beta_min, beta_max, n_scan)
    vals = [mismatch(bt) for bt in betas]
    vals = np.array(vals)

    roots = []
    for i in range(len(vals) - 1):
        if not np.isnan(vals[i]) and not np.isnan(vals[i+1]):
            if vals[i] * vals[i+1] < 0:
                try:
                    root = brentq(mismatch, betas[i], betas[i+1], xtol=1e-8)
                    roots.append(root)
                except Exception:
                    pass
    return roots


if __name__ == "__main__":
    Z_A, Z_B = 3.0, 1.0
    R = 3.015
    p0 = np.sqrt(10.0)
    c = p0 * R / 2

    print("=" * 70)
    print(f"LiH molecular Sturmian betas (prolate spheroidal)")
    print(f"Z_A={Z_A}, Z_B={Z_B}, R={R}, p0={p0:.4f}, c={c:.4f}")
    print("=" * 70)

    # Enumerate all molecular orbitals up to n_principal = 3
    # n_principal = 1 + (m + n_sph) + n_rad
    orbitals = []
    for n_princ in range(1, 4):
        for m in range(n_princ):
            for n_sph in range(n_princ - m):
                n_rad = n_princ - 1 - m - n_sph
                if n_rad >= 0:
                    orbitals.append((m, n_sph, n_rad, n_princ))

    print(f"\n{'n':>3s} {'(m,ns,nr)':>12s} {'beta':>10s} {'eps_Li':>10s} "
          f"{'eps_H':>10s} {'H_bound':>8s} {'time':>6s}")
    print("-" * 65)

    results = []
    for m, n_sph, n_rad, n_princ in orbitals:
        t0 = time.time()
        roots = find_beta(m, n_sph, n_rad, p0, R, Z_A, Z_B, beta_max=5.0)
        dt = time.time() - t0

        if roots:
            beta = roots[0]
            eps_Li = p0**2/2 - beta * Z_A * p0
            eps_H = p0**2/2 - beta * Z_B * p0
            bound = "YES" if eps_H < 0 else "NO"
            degen = 1 if m == 0 else 2
            results.append((n_princ, m, n_sph, n_rad, beta, eps_Li, eps_H, degen))
            print(f"{n_princ:3d} ({m},{n_sph},{n_rad}){'':>6s} {beta:10.6f} "
                  f"{eps_Li:10.4f} {eps_H:10.4f} {bound:>8s} {dt:5.1f}s")
        else:
            print(f"{n_princ:3d} ({m},{n_sph},{n_rad}){'':>6s} {'NO ROOT':>10s}"
                  f" {'---':>10s} {'---':>10s} {'---':>8s} {dt:5.1f}s")

    # Summary
    print("\n" + "=" * 70)
    print("Summary:")
    print(f"  Total molecular orbitals found: {len(results)}")
    print(f"  H-bound orbitals: {sum(1 for r in results if r[6] < 0)}")
    print(f"  Beta range: {min(r[4] for r in results):.4f} to "
          f"{max(r[4] for r in results):.4f}")
    print(f"  H 1s binding threshold: beta > {p0/2:.4f}")

    # Atomic limit check
    print("\n  Atomic limit check (beta vs 1/n):")
    for n_p, m, ns, nr, beta, eLi, eH, deg in results:
        atomic_beta = 1.0 / n_p
        print(f"    n={n_p} (m={m},ns={ns},nr={nr}): "
              f"beta={beta:.4f}, 1/n={atomic_beta:.4f}, "
              f"ratio={beta/atomic_beta:.4f}")

    # Check 2 gate
    beta_1s = results[0][4] if results else 0
    eps_H_1s = p0**2/2 - beta_1s * Z_B * p0
    print(f"\n  CHECK 2 (H 1s bound): eps = {eps_H_1s:.4f} Ha -> "
          f"{'PASS' if eps_H_1s < 0 else 'FAIL'}")
    if eps_H_1s > 0:
        beta_needed = p0 / (2 * Z_B)
        gap = beta_needed - beta_1s
        print(f"    beta(1sigma) = {beta_1s:.4f}, "
              f"needed = {beta_needed:.4f}, gap = {gap:.4f}")
        print(f"    Root cause: 1sigma is the Li 1s core MO (beta ≈ 1/n = 1)")
        print(f"    The H atom enters via the 2sigma MO (beta = "
              f"{results[1][4]:.4f} if available)")
