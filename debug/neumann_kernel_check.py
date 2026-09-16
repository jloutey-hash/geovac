"""Does the general-m Neumann kernel actually reproduce 1/r12?

The mu = 0 validation in `prolate_ci_general_m.validate_mu0` compares against
`geovac.neumann_vee`, which only ever evaluates the m = 0 term -- so it cannot
see an error in the m-dependent prefactor.  This checks the kernel itself,
pointwise, against 1/|r1 - r2| computed in Cartesian coordinates.

Also checks Paper 12's printed Eq. (neumann_full), which differs from the
standard form: it carries (2 - delta_m0) (l-m)!/(l+m)! with no (2l+1) and no
(-1)^m, while the standard expansion is

  1/r12 = (2/R) sum_l sum_m (-1)^m (2l+1) [(l-|m|)!/(l+|m|)!]^2
              P_l^|m|(xi_<) Q_l^|m|(xi_>) P_l^|m|(eta1) P_l^|m|(eta2) e^{i m dphi}

Run:  python debug/neumann_kernel_check.py
"""

from __future__ import annotations

import math
import numpy as np
from scipy.special import lpmv


def cart(R, xi, eta, phi):
    rho = math.sqrt((xi * xi - 1.0) * (1.0 - eta * eta))
    return np.array([(R / 2) * rho * math.cos(phi),
                     (R / 2) * rho * math.sin(phi),
                     (R / 2) * xi * eta])


def P_lm_xi(l, m, xi):
    """P_l^m(xi) for xi > 1 (hyperbolic branch): (xi^2-1)^{m/2} d^m P_l/dxi^m."""
    c = np.zeros(l + 1); c[l] = 1.0
    d = np.polynomial.legendre.legder(c, m) if m > 0 else c
    return (xi * xi - 1.0) ** (m / 2.0) * np.polynomial.legendre.legval(xi, d)


def Q_lm_xi(l, m, xi):
    """Q_l^m(xi) for xi > 1: (xi^2-1)^{m/2} d^m Q_l/dxi^m."""
    from numpy.polynomial import polynomial as P
    from numpy.polynomial import legendre as L
    d = [0.5 * math.log((xi + 1.0) / (xi - 1.0))]
    for k in range(1, m + 1):
        km = k - 1
        d.append(-((-1.0) ** km * math.factorial(km) * 0.5
                   * (1.0 / (xi - 1.0) ** k - 1.0 / (xi + 1.0) ** k)))
    wc = np.zeros(1)
    for k in range(1, l + 1):
        a = L.leg2poly(np.eye(l + 1)[k - 1][:k])
        b = L.leg2poly(np.eye(l + 1)[l - k][: l - k + 1])
        wc = P.polyadd(wc, P.polymul(a, b) / k)
    pl = L.leg2poly(np.eye(l + 1)[l])
    out = 0.0
    for a in range(m + 1):
        pla = P.polyder(pl, a) if a > 0 else pl
        out += math.comb(m, a) * P.polyval(xi, pla) * d[m - a]
    if l >= 1:
        wm = P.polyder(wc, m) if m > 0 else wc
        out -= P.polyval(xi, wm)
    return (xi * xi - 1.0) ** (m / 2.0) * out


def neumann_standard(R, p1, p2, l_max):
    """(2/R) sum (-1)^m (2l+1) [(l-m)!/(l+m)!]^2 P Q P P e^{i m dphi}, real form."""
    xi1, eta1, phi1 = p1
    xi2, eta2, phi2 = p2
    lo, hi = (xi1, xi2) if xi1 < xi2 else (xi2, xi1)
    tot = 0.0
    for l in range(l_max + 1):
        for m in range(0, l + 1):
            r = math.factorial(l - m) / math.factorial(l + m)
            pref = (2 - (1 if m == 0 else 0)) * (-1.0) ** m * (2 * l + 1) * r * r
            tot += (pref * P_lm_xi(l, m, lo) * Q_lm_xi(l, m, hi)
                    * lpmv(m, l, eta1) * lpmv(m, l, eta2)
                    * math.cos(m * (phi1 - phi2)))
    return (2.0 / R) * tot


def neumann_paper12(R, p1, p2, l_max):
    """Eq. (neumann_full) exactly as Paper 12 prints it."""
    xi1, eta1, phi1 = p1
    xi2, eta2, phi2 = p2
    lo, hi = (xi1, xi2) if xi1 < xi2 else (xi2, xi1)
    tot = 0.0
    for l in range(l_max + 1):
        for m in range(0, l + 1):
            pref = ((2 - (1 if m == 0 else 0))
                    * math.factorial(l - m) / math.factorial(l + m))
            tot += (pref * P_lm_xi(l, m, lo) * Q_lm_xi(l, m, hi)
                    * lpmv(m, l, eta1) * lpmv(m, l, eta2)
                    * math.cos(m * (phi1 - phi2)))
    return (2.0 / R) * tot


def main():
    rng = np.random.default_rng(20260914)
    R = 1.4011
    print(f"{'pt':>3} {'1/r12 exact':>14} {'standard':>14} {'rel':>10}"
          f" | {'Paper 12 Eq':>14} {'rel':>10}")
    bad_std = bad_p12 = 0
    for k in range(8):
        xi1 = 1.0 + rng.uniform(0.2, 2.5)
        xi2 = 1.0 + rng.uniform(0.2, 2.5)
        eta1 = rng.uniform(-0.9, 0.9)
        eta2 = rng.uniform(-0.9, 0.9)
        phi1 = rng.uniform(0, 2 * math.pi)
        phi2 = rng.uniform(0, 2 * math.pi)
        p1, p2 = (xi1, eta1, phi1), (xi2, eta2, phi2)
        exact = 1.0 / np.linalg.norm(cart(R, *p1) - cart(R, *p2))
        s = neumann_standard(R, p1, p2, 40)
        q = neumann_paper12(R, p1, p2, 40)
        rs = abs(s - exact) / exact
        rq = abs(q - exact) / exact
        bad_std += rs > 1e-6
        bad_p12 += rq > 1e-6
        print(f"{k:3d} {exact:14.9f} {s:14.9f} {rs:10.2e}"
              f" | {q:14.9f} {rq:10.2e}")
    print()
    print(f"  standard form   : {8-bad_std}/8 points agree to 1e-6")
    print(f"  Paper 12 Eq.    : {8-bad_p12}/8 points agree to 1e-6")

    # r12 in prolate spheroidals: does Paper 12's Eq. (r12_prolate) hold?
    print("\n  Paper 12 Eq. (r12_prolate) omits cos(phi1-phi2)?")
    xi1, eta1, xi2, eta2 = 1.7, 0.3, 2.2, -0.4
    for dphi in (0.0, 1.0, 2.5):
        p1, p2 = (xi1, eta1, 0.0), (xi2, eta2, dphi)
        exact2 = np.sum((cart(R, *p1) - cart(R, *p2)) ** 2)
        r1 = (xi1 ** 2 - 1) * (1 - eta1 ** 2)
        r2 = (xi2 ** 2 - 1) * (1 - eta2 ** 2)
        as_printed = (R ** 2 / 4) * (r1 + r2 - 2 * math.sqrt(r1 * r2)
                                     + (xi1 * eta1 - xi2 * eta2) ** 2)
        with_cos = (R ** 2 / 4) * (r1 + r2 - 2 * math.sqrt(r1 * r2) * math.cos(dphi)
                                   + (xi1 * eta1 - xi2 * eta2) ** 2)
        print(f"    dphi={dphi:4.1f}  exact r12^2={exact2:.8f}  "
              f"as printed={as_printed:.8f}  with cos={with_cos:.8f}")


if __name__ == "__main__":
    main()
