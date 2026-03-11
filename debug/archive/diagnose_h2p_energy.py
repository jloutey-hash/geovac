"""Debug: H2+ energy from prolate spheroidal separation — Neumann BC fix.

For m=0: R(xi=1) is FINITE (regular solution). The correct BC is
  dR/dxi(xi=1) = 0 (Neumann), NOT R(1)=0 (Dirichlet).

Using the Neumann BC at the left boundary fixes the ground state eigenvalue.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal, eigh
from scipy.optimize import brentq
from scipy.special import obl_cv
from math import factorial
import time


def angular_sep(m_abs, n_sph, c, b=0.0, n_basis=40):
    """A_sep from angular equation (eigenvalue of H_eta from top)."""
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


def radial_evals_neumann(m_abs, c, a, n_grid=600, xi_max=8.0, n_evals=5):
    """Radial eigenvalues of L_xi with NEUMANN BC at left for m=0.

    L_xi = d/dxi[(xi^2-1)d/dxi] - c^2*xi^2 + a*xi - m^2/(xi^2-1)
    """
    N = n_grid
    xi_min = 1.0 + 0.001  # very close to 1 but avoid singularity
    h = (xi_max - xi_min) / (N + 1)
    xi = xi_min + (np.arange(N) + 1) * h  # N interior points

    xi2_1 = xi**2 - 1
    q = -c**2 * xi**2 + a * xi - m_abs**2 / xi2_1

    p_plus = (xi + h/2)**2 - 1
    p_minus = (xi - h/2)**2 - 1

    # Standard tridiagonal: Dirichlet at both ends
    diag = -(p_plus + p_minus) / h**2 + q
    off = p_plus[:-1] / h**2

    if m_abs == 0:
        # Neumann BC at left: dR/dxi = 0 at xi = xi_min
        # Ghost point: f_{-1} = f_0
        # Modified first row: L_{0,0} = -p_plus[0]/h^2 + q[0]
        # (the p_minus contribution vanishes because f_{-1} = f_0)
        diag[0] = -p_plus[0] / h**2 + q[0]
        # L_{0,1} = p_plus[0]/h^2 (unchanged)
    # else: Dirichlet R(xi_min)=0 is correct for m>0

    # Right Dirichlet R(xi_max)=0 is correct for all m

    evals = eigh_tridiagonal(diag, off, eigvals_only=True)
    evals = np.sort(evals)
    return evals[::-1][:n_evals]  # descending (most positive first)


# H2+ at R=2 bohr
Z_A, Z_B = 1.0, 1.0
R = 2.0

print("=" * 60)
print("H2+ energy scan with Neumann BC (m=0)")
print("=" * 60)

print(f"\n{'E':>8s} {'c':>8s} {'A_ang':>10s} {'L_top':>10s} "
      f"{'A+L':>10s}")
print("-" * 55)

prev_sum = None
E_cross = None
for E in np.linspace(-0.3, -2.0, 80):
    p0 = np.sqrt(-2*E)
    c = p0 * R / 2
    a = (Z_A + Z_B) * R
    b = (Z_B - Z_A) * R

    A_ang = angular_sep(0, 0, c, b)
    L_ev = radial_evals_neumann(0, c, a)
    L_top = L_ev[0]
    s = A_ang + L_top  # matching: A_ang = -L_top, so s should be 0

    if prev_sum is not None and prev_sum * s < 0 and E_cross is None:
        E_cross = E

    print(f"{E:8.4f} {c:8.4f} {A_ang:10.4f} {L_top:10.4f} {s:+10.4f}")
    prev_sum = s

if E_cross:
    print(f"\n  Zero crossing near E = {E_cross:.4f} Ha")
    print(f"  (Expected: E = -1.1026 Ha for H2+ at R=2)")

# Also try convergence: vary n_grid
print("\n=== Radial eigenvalue convergence test (E=-1.1, c=1.483) ===")
c_test = np.sqrt(2*1.1) * 2 / 2
a_test = 4.0
for ng in [100, 200, 400, 800, 1600]:
    ev = radial_evals_neumann(0, c_test, a_test, n_grid=ng)
    print(f"  n_grid={ng:5d}: L_top = {ev[0]:.8f}")
