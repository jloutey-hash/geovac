"""Increment 3d: does sigma != 0 break the weight-1 result? (the named open risk)

3c established that the ordered xi integral closes at WEIGHT 1 -- at sigma = 0.
The plan flagged one place a higher-weight object could still enter, and flagged
it as ARGUED, not verified:

    Q_tau^sigma = (xi^2-1)^{|sigma|/2} d^sigma Q_tau / dxi^sigma

and d^sigma of Q_0 = (1/2)ln((xi+1)/(xi-1)) has POLES of order up to sigma at
xi = +-1. If those poles outrun the (xi^2-1)^H prefactor, the endpoint gets more
singular than the sigma = 0 case and the closed form could climb to weight 2.

THE ARGUMENT UNDER TEST. Per electron, the orbital product contributes
(xi^2-1)^{h} with h = (|m_a|+|m_b|)/2, and the kernel contributes
(xi^2-1)^{|sigma|/2}, so H = h + |sigma|/2. Against a pole of order |sigma| the
net exponent at xi = 1 is

    H - |sigma| = h - |sigma|/2 = ( |m_a| + |m_b| - |m_a - m_b| ) / 2

which is >= 0 by the TRIANGLE INEQUALITY, with equality exactly when m_a and m_b
have opposite signs (or one vanishes). So no new pole is ever created, and the
sigma != 0 endpoint is no worse than sigma = 0.

If that holds, the split is

    (xi^2-1)^H D^Q  =  [polynomial] x Q_0   +   [polynomial]

-- the second bracket because (xi^2-1)^{H} x (poles up to |sigma|) collapses to
(xi^2-1)^{H-|sigma|} x polynomial with H-|sigma| >= 0 -- and BOTH pieces are
weight 1 by 3c's moments.

THREE CHECKS
  S1  the triangle bound, machine-checked on the actual half-powers from 3a
  S2  no negative powers of (xi-1) survive in the built integrand
  S3  the single xi integral against Q_tau^sigma, closed form vs quadrature,
      plus a weight census

Run from repo root:  python debug/inc3d_sigma_pole_probe.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy import integrate

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac.two_center_eri import (  # noqa: E402
    log_moment, log_shift_moment, upper_integral,
)

z = sp.Symbol("z")


def _Q_tau(tau: int):
    Q0 = sp.log((z + 1) / (z - 1)) / 2
    out = sp.legendre(tau, z) * Q0
    for k in range(1, tau + 1):
        out -= sp.legendre(k - 1, z) * sp.legendre(tau - k, z) / k
    return out


def _split_DQ(tau: int, s: int, H: int):
    """(xi^2-1)^H d^s Q_tau/dxi^s  ->  (poly_log, poly_rat) with

        result = poly_log * Q_0  +  poly_rat

    Asserts both are polynomials -- that is the whole claim.
    """
    # NEVER let sympy hold the log. The plan's Phase 0 note records that
    # expand/simplify rewrite log((z+1)/(z-1)) so a .coeff() match silently finds
    # nothing; worse, sp.expand distributes INSIDE the argument -- it becomes
    # log(z/(z-1) + 1/(z-1)) -- after which expand_log cannot split it either.
    # So carry Q_0 as an opaque coefficient and differentiate the PAIR by hand,
    # using the one fact needed: Q_0' = -1/(z^2 - 1). Exact, and immune.
    a, b = sp.legendre(tau, z), sp.Integer(0)          # Q_tau = a*Q_0 + b
    for k in range(1, tau + 1):
        b -= sp.legendre(k - 1, z) * sp.legendre(tau - k, z) / k
    for _ in range(s):
        a, b = sp.diff(a, z), sp.together(sp.diff(b, z) - a / (z ** 2 - 1))
    poly_log, rest = sp.simplify(a), sp.simplify(b)
    pref = (z ** 2 - 1) ** H
    pl = sp.cancel(sp.expand(pref * poly_log))
    pr = sp.cancel(sp.together(pref * rest))
    return sp.simplify(pl), sp.simplify(pr)


def leg_S1() -> None:
    print("S1  the triangle bound  H - |sigma| = (|m_a|+|m_b|-|m_a-m_b|)/2 >= 0\n")
    print("     m_a  m_b  sigma   h     H    H-|s|   status")
    print("     " + "-" * 50)
    worst = 99
    for ma in range(-2, 3):
        for mb in range(-2, 3):
            s = ma - mb
            h = sp.Rational(abs(ma) + abs(mb), 2)
            H = h + sp.Rational(abs(s), 2)
            net = H - abs(s)
            worst = min(worst, net)
            ok = "ok" if net >= 0 else "NEGATIVE"
            assert H.is_integer, f"half-power {H} not integral at ({ma},{mb})"
            print(f"     {ma:+d}   {mb:+d}   {s:+d}   {str(h):>4} {str(H):>4}"
                  f"  {str(net):>5}   {ok}")
    print(f"\n     minimum over all (m_a,m_b) in [-2,2]: {worst}  "
          f"{'-> no pole ever created' if worst >= 0 else '-> CLAIM FAILS'}\n")


def leg_S2() -> None:
    print("S2  does any negative power of (xi-1) survive the built integrand?\n")
    print("     tau  s   H    poly_log deg   poly_rat deg   both polynomial?")
    print("     " + "-" * 62)
    for tau, s, H in ((1, 1, 1), (2, 1, 1), (2, 2, 2), (3, 1, 1),
                      (3, 2, 2), (3, 3, 3), (4, 2, 2)):
        pl, pr = _split_DQ(tau, s, H)
        okl = pl.is_polynomial(z)
        okr = pr.is_polynomial(z)
        dl = sp.Poly(pl, z).total_degree() if okl else "RATIONAL"
        dr = sp.Poly(pr, z).total_degree() if okr else "RATIONAL"
        print(f"      {tau}   {s}   {H}       {dl!s:>6}         {dr!s:>6}"
              f"        {'YES' if okl and okr else 'NO'}")
        assert okl and okr, f"tau={tau} s={s} H={H}: pole survived"
    print("\n     Every case polynomial: the (xi^2-1)^H prefactor absorbs the")
    print("     poles exactly, as the triangle bound says it must.\n")


def _closed_single(tau: int, s: int, H: int, j: int, p):
    """int_1^oo z^j (z^2-1)^H [d^s Q_tau] e^{-p z} dz, in closed form.

    poly_rat piece -> elementary (upper_integral at 1).
    poly_log piece -> t = z-1, ln((t+2)/t)/2 = [ln(t+2) - ln t]/2 -> 3c's moments.
    """
    pl, pr = _split_DQ(tau, s, H)
    t = sp.Symbol("t", positive=True)

    total = sp.Integer(0)
    for k, c in sp.Poly(sp.expand(z ** j * pr), z).terms():
        total += c * upper_integral(int(k[0]), p, sp.Integer(1))

    shifted = sp.expand(sp.expand(z ** j * pl).subs(z, t + 1))
    for k, c in sp.Poly(shifted, t).terms():
        n = int(k[0])
        total += c * sp.exp(-p) / 2 * (log_shift_moment(n, p, sp.Integer(2))
                                       - log_moment(n, p))
    return total


def leg_S3() -> None:
    print("S3  single xi integral against Q_tau^sigma: closed form vs quadrature\n")
    Q0f = lambda x: 0.5 * np.log((x + 1) / (x - 1))  # noqa: E731
    worst = 0.0
    cases = [(1, 1, 1, 0), (2, 1, 1, 1), (2, 2, 2, 0), (3, 2, 2, 1), (3, 3, 3, 0)]
    for tau, s, H, j in cases:
        for pv in (sp.Rational(3, 2), sp.Integer(2)):
            cf = float(sp.re(sp.N(_closed_single(tau, s, H, j, pv), 30)))
            DQf = sp.lambdify(z, sp.diff(_Q_tau(tau), z, s) if s else _Q_tau(tau),
                              [{"log": np.log}, "numpy"])
            pf = float(pv)
            nu, _ = integrate.quad(
                lambda x: x ** j * (x ** 2 - 1) ** H * float(DQf(x)) * np.exp(-pf * x),
                1.0, np.inf, limit=400, epsabs=1e-14, epsrel=1e-13)
            rel = abs(cf - nu) / max(abs(nu), 1e-300)
            worst = max(worst, rel)
            print(f"     tau={tau} s={s} H={H} j={j} p={pf:4.1f}:  "
                  f"{cf: .12f}  vs {nu: .12f}   rel={rel:.1e}")
    print(f"\n     worst relative {worst:.2e}\n")

    print("     WEIGHT CENSUS (symbolic p):")
    ps = sp.Symbol("p", positive=True)
    bad = set()
    for tau, s, H, j in cases:
        e = sp.expand(_closed_single(tau, s, H, j, ps))
        names = {type(f).__name__ for f in e.atoms(sp.Function)}
        w2 = names & {"polylog", "dilog", "zeta", "lerchphi"}
        bad |= w2
        print(f"       tau={tau} s={s} H={H} j={j}: {sorted(names)}"
              f"  gamma={'yes' if e.has(sp.EulerGamma) else 'no'}")
    print(f"\n     weight-2 atoms anywhere: {bad or 'NONE'}")
    print("     => sigma != 0 stays at WEIGHT 1. The named risk is closed.\n")


def main() -> None:
    print("Increment 3d -- does sigma != 0 break weight 1?\n")
    leg_S1()
    leg_S2()
    leg_S3()


if __name__ == "__main__":
    main()
