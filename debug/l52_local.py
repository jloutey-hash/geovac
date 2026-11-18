"""Record L52 (docs/predictions/2026-11-17_aa_bc_contact_amplitudes.md), shared local-pinch formulas.

Contact x* = (D1^2 + D2^2)/(2 D1 D2), where |A - C|^2 = 2 D1 D2 (x* - x) vanishes (B origin, C = (0,0,D1),
A = D2 (sin g, 0, cos g), x = cos g; NA/AB convention).

Pinch formula (Feynman parameters, 3-D): for an analytic weight g,
    int |s|^p |s - a|^q g(A + s) d^3 s  has non-analytic part
    C_pq |a|^(p+q+3) [ int_0^1 u^((p+1)/2) (1-u)^((q+1)/2) g(A + u a) du  +  O(|a|^2) ],
    C_pq = pi^(3/2) Gamma(-(p+q+3)/2) / (Gamma(-p/2) Gamma(-q/2)),
s = r - A, a = C - A; p = power of the A-side singular term, q = power of the C-side cusp term.
At x = x*, |A + u(C - A) - B|^2 = (1-u) D2^2 + u D1^2 (real), so the leading amplitude is a 1-D real integral.

Amplitude of (x* - x)^gamma, gamma = (p+q+3)/2:
    Phi = v_p * w_q * C_pq * (2 D1 D2)^gamma * int_0^1 u^((p+1)/2) (1-u)^((q+1)/2) e^{-zB R_B(u)} du,
v_p = coefficient of |s|^p in V_A, w_q = N_BC (-zC)^q / q!  (cusp coefficient of rho_BC at C, times e^{-zB r_B}).

Legendre coefficients of (x* - x)^gamma (Rodrigues, L-fold integration by parts; positive integrand):
    T_L = (2L+1)/2 * (-1)^L * [Gamma(gamma+1)/Gamma(gamma+1-L)] / (2^L L!) * int_{-1}^{1} (1-x^2)^L (x*-x)^(gamma-L) dx.
"""
from __future__ import annotations

import mpmath as mp


def C_pq(p, q):
    return mp.pi ** 1.5 * mp.gamma(-mp.mpf(p + q + 3) / 2) / (mp.gamma(-mp.mpf(p) / 2) * mp.gamma(-mp.mpf(q) / 2))


def xstar(D1, D2):
    return (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)


def seg(p, q, D1, D2, zB):
    f = lambda u: u ** (mp.mpf(p + 1) / 2) * (1 - u) ** (mp.mpf(q + 1) / 2) * mp.exp(-zB * mp.sqrt((1 - u) * D2 ** 2 + u * D1 ** 2))
    return mp.quad(f, [0, 0.5, 1])


def amplitude(p, q, vp, D1, D2, zB, zC):
    """coefficient of (x* - x)^((p+q+3)/2) in int V rho_BC, V's |s|^p coefficient vp, normalized 1s B, C."""
    N = mp.sqrt(zB ** 3 * zC ** 3) / mp.pi
    wq = N * (-zC) ** q / mp.factorial(q)
    gam = mp.mpf(p + q + 3) / 2
    return vp * wq * C_pq(p, q) * (2 * D1 * D2) ** gam * seg(p, q, D1, D2, zB)


def T(gam, xs, L):
    gam = mp.mpf(gam)
    pref = (2 * L + 1) / mp.mpf(2) * (-1) ** L * mp.gammaprod([gam + 1], [gam + 1 - L]) / (mp.mpf(2) ** L * mp.factorial(L))
    x0 = 1 / (xs + mp.sqrt(xs * xs - 1)) if xs > 1 else mp.mpf(1)
    w = 1 / mp.sqrt(L + 1)
    pts = sorted(set([mp.mpf(-1), max(-1, x0 - 6 * w), x0, min(1, x0 + 6 * w), mp.mpf(1)]))
    f = lambda x: (1 - x * x) ** L * (xs - x) ** (gam - L)
    return pref * mp.quad(f, pts)


def ratio_extrap(Ls, R):
    """fit R_L = r0 + c1/L + c2/L^2 + c3/L^3 by least squares on the given L; return r0, c1."""
    A = mp.matrix([[1, mp.mpf(1) / L, mp.mpf(1) / L ** 2, mp.mpf(1) / L ** 3] for L in Ls])
    b = mp.matrix(R)
    sol = mp.lu_solve(A.T * A, A.T * b)
    return sol[0], sol[1]
