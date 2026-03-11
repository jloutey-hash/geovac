"""Debug: fix radial eigenvalue solver for prolate spheroidal coords.

Radial equation (mu in [1, inf)):
  d/dmu[(mu^2-1)dR/dmu] + [-A + c^2*mu^2 + a*mu - m^2/(mu^2-1)]R = 0

where A is the separation constant. For bound states R -> 0 as mu -> inf.
The n_rad-th eigenvalue has n_rad nodes in (1, inf).

Note: the SIGN convention here uses -A (not +lambda). The angular equation
gives A from obl_cv. Both equations share the SAME A.

Matching condition: A_angular(c, b) = A_radial(c, a, n_rad)
"""
import numpy as np
from scipy.special import obl_cv
from math import factorial


def radial_eigenvalue(m_abs, n_rad, c, a, n_steps=10000, mu_max=50.0):
    """Find the separation constant A for the n_rad-th radial bound state.

    The radial eq with separation const A:
    (mu^2-1)*R'' + 2*mu*R' + [-A + c^2*mu^2 + a*mu - m^2/(mu^2-1)]*R = 0

    Node count INCREASES with A (more negative -A = deeper well = more nodes).
    We bisect to find A where node_count transitions from n_rad-1 to n_rad.
    """
    h = (mu_max - 1.0) / n_steps

    def count_nodes(A_val):
        mu = 1.0 + 1e-6
        if m_abs == 0:
            R_val = 1.0
            R_prime = 0.01
        else:
            R_val = (mu - 1.0)**(m_abs / 2.0)
            R_prime = (m_abs / 2.0) * (mu - 1.0)**(m_abs / 2.0 - 1.0)

        nodes = 0
        for step in range(n_steps):
            mu = 1.0 + 1e-6 + step * h
            mu2_1 = mu**2 - 1.0
            if mu2_1 < 1e-12:
                mu2_1 = 1e-12

            # -A + c^2*mu^2 + a*mu - m^2/(mu^2-1)
            pot = -A_val + c**2 * mu**2 + a * mu - m_abs**2 / mu2_1
            R_pp = (-2 * mu * R_prime + pot * R_val) / mu2_1

            R_new = R_val + h * R_prime
            R_prime_new = R_prime + h * R_pp

            if R_new * R_val < 0 and step > 10:
                nodes += 1

            R_val = R_new
            R_prime = R_prime_new

            if abs(R_val) > 1e40:
                R_val /= 1e40
                R_prime /= 1e40

        return nodes

    # Determine node count direction vs A
    # Test at A = 0 and A = 100
    n0 = count_nodes(0.0)
    n100 = count_nodes(100.0)
    n_neg = count_nodes(-50.0)
    print(f"  Node counts: A=-50 -> {n_neg}, A=0 -> {n0}, A=100 -> {n100}")

    # Scan to understand behavior
    A_scan = np.linspace(-50, 300, 50)
    print(f"  A scan (m={m_abs}, n_rad={n_rad}, c={c:.2f}, a={a:.2f}):")
    prev_nn = None
    bracket = None
    for A_val in A_scan:
        nn = count_nodes(A_val)
        if nn <= 5:  # Only print interesting range
            print(f"    A={A_val:8.2f}: nodes={nn}")
        if prev_nn is not None and prev_nn < n_rad + 1 and nn >= n_rad + 1:
            # Node count increased past target
            bracket = (A_val - (A_scan[1] - A_scan[0]), A_val)
        if prev_nn is not None and prev_nn >= n_rad + 1 and nn < n_rad + 1:
            bracket = (A_val - (A_scan[1] - A_scan[0]), A_val)
        prev_nn = nn

    if bracket:
        print(f"  Bracket found: A in [{bracket[0]:.2f}, {bracket[1]:.2f}]")

    return bracket


if __name__ == "__main__":
    print("=== Radial eigenvalue investigation ===\n")

    # Simple test case
    print("--- Test: c=2, a=3, m=0, n_rad=0 ---")
    radial_eigenvalue(0, 0, 2.0, 3.0)

    print("\n--- Test: c=5, a=6, m=0, n_rad=0 ---")
    radial_eigenvalue(0, 0, 5.0, 6.0)

    # LiH parameters at beta=1
    print("\n--- LiH at beta=1: c=4.77, a=6.03, m=0, n_rad=0 ---")
    p0 = np.sqrt(10.0)
    R = 3.015
    beta = 1.0
    c = beta * p0 * R / 2
    a = beta * 4.0 * R / 2  # Z_tot = 4
    radial_eigenvalue(0, 0, c, a)
