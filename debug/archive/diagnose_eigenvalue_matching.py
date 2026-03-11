"""Debug: fix angular + radial eigenvalue solvers for prolate spheroidal coords.

Angular spheroidal equation (oblate convention, c^2 > 0):
  d/dnu[(1-nu^2)dS/dnu] + [A + c^2*nu^2 - m^2/(1-nu^2)]S = 0

where A is the oblate eigenvalue (scipy obl_cv returns this).
Our Coulomb separation constant lambda = -A.

The matrix method expands S in associated Legendre functions P_r^m.
The KEY recurrence for nu is:
  (2r+1)*nu*P_r^m = (r-m+1)*P_{r+1}^m + (r+m)*P_{r-1}^m
"""
import numpy as np
from scipy.special import obl_cv
from math import factorial


def angular_eigenvalue_oblate(m_abs, n_total, c):
    """Compute oblate angular eigenvalue A using matrix method.

    n_total = m + n_sph (the total angular quantum number in scipy convention).
    Returns A such that obl_cv(m, n_total, c) should equal A.
    """
    N = 40  # basis size
    r_vals = np.arange(m_abs, m_abs + N, dtype=float)

    # Normalization: int_{-1}^1 [P_r^m(x)]^2 dx = 2/(2r+1) * (r+m)!/(r-m)!
    norms = np.array([2.0/(2*r+1) * factorial(int(r+m_abs))/factorial(int(r-m_abs))
                       for r in r_vals])

    # nu matrix in orthonormal basis: Y_r = P_r^m / sqrt(norm_r)
    # <Y_s|nu|Y_r> is tridiagonal from the recurrence
    nu_mat = np.zeros((N, N))
    for i in range(N - 1):
        r = r_vals[i]
        # From (2r+1)*nu*P_r = (r-m+1)*P_{r+1} + (r+m)*P_{r-1}
        # <P_{r+1}|nu|P_r> = (r-m+1)/(2r+1) * norm_{r+1}  ... NO
        # Multiply both sides by P_{r+1} and integrate:
        # (2r+1)*<P_{r+1}|nu|P_r> = (r-m+1)*<P_{r+1}|P_{r+1}> = (r-m+1)*norm_{r+1}
        # So <P_{r+1}|nu|P_r> = (r-m+1)*norm_{r+1}/(2r+1)
        # In orthonormal basis:
        # <Y_{r+1}|nu|Y_r> = <P_{r+1}|nu|P_r> / sqrt(norm_{r+1}*norm_r)
        #                   = (r-m+1)*sqrt(norm_{r+1}) / ((2r+1)*sqrt(norm_r))
        val = (r - m_abs + 1) * np.sqrt(norms[i+1]) / ((2*r+1) * np.sqrt(norms[i]))
        nu_mat[i+1, i] = val
        nu_mat[i, i+1] = val

    # The angular equation in the P_r^m basis:
    # L_0 P_r^m = -r(r+1) P_r^m  (where L_0 = d/dnu[(1-nu^2)d/dnu] - m^2/(1-nu^2))
    # Full equation: L_0 S + (A + c^2*nu^2) S = 0
    # => -r(r+1) a_r + A*a_r + c^2 * sum_s <r|nu^2|s> a_s = 0
    # => H a = -A a  where H_{rs} = -r(r+1)*delta_{rs} + c^2*<r|nu^2|s>
    # Wait: L_0 S + A S + c^2*nu^2 S = 0
    # => L_0 S + c^2*nu^2 S = -A S
    # In matrix form: (L_0 + c^2*nu^2) a = -A * a
    # H = L_0 + c^2*nu^2, eigenvalues of H are -A

    nu2_mat = nu_mat @ nu_mat

    H = np.diag(-r_vals * (r_vals + 1))  # L_0 diagonal
    H += c**2 * nu2_mat

    evals = np.sort(np.linalg.eigvalsh(H))
    # evals are -A in ascending order (most negative first)
    # A = -evals, so A values are descending
    # We want A sorted ascending (smallest first)
    A_vals = -evals[::-1]  # ascending

    # Index: n_total - m_abs gives the index within this m block
    idx = n_total - m_abs
    return A_vals[idx]


if __name__ == "__main__":
    print("=== Angular eigenvalue verification ===\n")
    for m in range(3):
        for n_total in range(m, m + 4):
            for c_test in [0.5, 2.0, 5.0]:
                A_our = angular_eigenvalue_oblate(m, n_total, c_test)
                A_scipy = obl_cv(m, n_total, c_test)
                diff = abs(A_our - A_scipy)
                status = "OK" if diff < 1e-4 else "FAIL"
                print(f"  m={m}, n={n_total}, c={c_test}: "
                      f"ours={A_our:.6f}, scipy={A_scipy:.6f}, "
                      f"diff={diff:.2e} [{status}]")
        print()
