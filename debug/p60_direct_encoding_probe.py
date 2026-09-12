"""Can G be block-encoded DIRECTLY, recovering the second power of n?

v5.11.1 priced the lever: composing P^-1/2 with (I-C) makes alpha inherit
||P^-1/2||^2 ~ n^2, so the total is n^2 where the invariant floor allows n.
The waste is large and named: ||G|| = 0.372 against a composed alpha = 3196.

The candidate.  G's symbol is the RATIO of the two symbols, and the ratio is
bounded -- both numerator and denominator vanish quadratically at chi = pi, so
the quotient tends to (kR)^2/24, and at chi -> 0 it tends to 1/4:

    ratio(s) = (1 - sinc(s)) (s^2 + (kR)^2) / (4 s^2),    s = kR cot(chi/2)

So the Toeplitz-minus-Hankel matrix B built from the ratio's own cosine
coefficients is directly constructible, and a circulant-embedded encoding of it
carries alpha = ||ratio||_inf = O(1) rather than O(n^2).

Whether that HELPS is a different question, and it is the one measured here.
Whitening with X = P^-1/2 B^-1/2 gives X^T (I-C) X = B^-1/2 G B^-1/2, so the
decisive quantity is cond(B^-1/2 G B^-1/2): if it is bounded and small, the
direct object whitens to O(1) residual conditioning and the floor is reached.

Coefficients: subtract the s -> infinity limit 1/4 (whose cosine coefficients
vanish for j >= 1) so the integrand decays like 1/s^3 and the tail is
controlled; add it back into r_0 analytically.
"""
import numpy as np
from geovac.sturmian_sigma_law import sw_cross_block

GL_N = 48
_xg, _wg = np.polynomial.legendre.leggauss(GL_N)
KR = 2.0
M_QUAD = 200_001


def ratio_minus_quarter(s, kR):
    """[(kR)^2 - sinc(s)(s^2+(kR)^2)] / (4 s^2), the ratio symbol less its 1/4 limit."""
    with np.errstate(divide="ignore", invalid="ignore"):
        v = (kR**2 - np.sinc(s / np.pi) * (s**2 + kR**2)) / (4 * s**2)
    # s -> 0 limit: ((kR)^2/6 - 1)/4
    return np.where(s < 1e-8, (kR**2 / 6.0 - 1.0) / 4.0, v)


def ratio_coeff(j, kR, Smax=3000.0):
    """(1/pi) int_0^pi cos(j chi) ratio(chi) dchi, via the s substitution."""
    edges, s = [0.0], 0.0
    while s < Smax:
        s += min(np.pi / 2, np.pi / max(2 * j * kR / (kR * kR + s * s), 1e-300))
        edges.append(min(s, Smax))
    e = np.asarray(edges)
    mid, half = 0.5 * (e[:-1] + e[1:]), 0.5 * (e[1:] - e[:-1])
    sv = (mid[:, None] + half[:, None] * _xg[None, :]).ravel()
    wv = (half[:, None] * _wg[None, :]).ravel()
    integ = (np.cos(2 * j * np.arctan2(kR, sv)) * ratio_minus_quarter(sv, kR)
             * 2 * kR / (kR * kR + sv * sv))
    val = float(np.dot(wv, integ) / np.pi)
    return val + 0.25 if j == 0 else val


def tmh(coeffs, n):
    g = lambda j: coeffs.get(j, 0.0)
    return np.array([[g(abs(a - b)) - g(a + b) for b in range(1, n + 1)]
                     for a in range(1, n + 1)])


def tridiag(n):
    return (np.diag(2.0 * np.ones(n)) + np.diag(np.ones(n - 1), 1)
            + np.diag(np.ones(n - 1), -1))


def inv_sqrt(M):
    ev, U = np.linalg.eigh(M)
    assert ev.min() > 0, f"not PD, lam_min={ev.min():.3e}"
    return U @ np.diag(ev ** -0.5) @ U.T


# sup of the ratio symbol -- the alpha a direct encoding would carry
chi = np.linspace(1e-9, np.pi - 1e-12, 2_000_00)
s_of_chi = KR / np.tan(chi / 2)
sup_ratio = float(np.max(np.abs(ratio_minus_quarter(s_of_chi, KR) + 0.25)))
print(f"||ratio||_inf = {sup_ratio:.4f}   <- the alpha a DIRECT encoding carries")
print(f"(composed alpha at n=160 was 3196; ||G|| = 0.372)\n")

print(f"  {'n':>5} {'cond(B)':>9} {'||G-B||/||G||':>14} {'cond(B^-1/2 G B^-1/2)':>23}"
      f" {'alpha_X/floor':>14}")
for n in (20, 40, 80, 160):
    coeffs = {j: ratio_coeff(j, KR) for j in range(0, 2 * n + 2)}
    B = tmh(coeffs, n)
    C = sw_cross_block(KR, n, M=M_QUAD)
    A = np.eye(n) - C
    P_is = inv_sqrt(tridiag(n))
    G = P_is @ A @ P_is

    Bis = inv_sqrt(B)
    resid = Bis @ G @ Bis
    X = P_is @ Bis
    floor = np.linalg.norm(inv_sqrt(A), 2)
    print(f"  {n:5d} {np.linalg.cond(B):9.4f} {np.linalg.norm(G-B,2)/np.linalg.norm(G,2):14.4f}"
          f" {np.linalg.cond(resid):23.4f} {np.linalg.norm(X,2)/floor:14.4f}")

print("\nReading: cond(B^-1/2 G B^-1/2) bounded => the DIRECT object whitens to O(1)")
print("residual conditioning, and alpha_X sits at the invariant floor => total O(n).")
