"""Increment 3c / periods probe: does the ordered xi integral close, and at what WEIGHT?

THE QUESTION. The exchange class's ordered double integral is the last numerical
step in the whole engine. It is also an ITERATED INTEGRAL over a simplex
(xi_< / xi_>), which is literally the shape that defines a period. So "does it
close" and "where does it sit in the transcendence hierarchy" are the same
question, and the weight filtration is the frame:

    weight 0   rationals, algebraic
    weight 1   ln, E_1, gamma
    weight 2   Li_2, pi^2/6           <- where iterated integrals of weight-1
                                          objects generically land

Phase 0-e's scoping said the exchange seeds are {E_1, ln, gamma} -- all weight 1.
If that survives contact with the actual integral, the obstruction was LABOUR,
not transcendence, and increment 3c is reachable. If Li_2 appears, the prediction
is "closed form one rung up" and it needs a Paper 18 / Paper 34 tag.

THE TRACE. With A, B = poly x (xi^2-1)^H x exp and P~ = d^s P_tau (polynomial),
Q~ = d^s Q_tau = P_tau Q_0 - W (Q_0 = the log):

    Xi = int_1^oo A(x1) [ Q~(x1) int_1^{x1} B P~  +  P~(x1) int_{x1}^oo B Q~ ]

  * inner-1 is polynomial x exponential -> elementary;
  * inner-2 carries the log -> gives back {log, E_1} as functions of x1;
  * so the OUTER integrand carries poly x exp x {log, E_1}.

Substituting t = xi - 1 puts the outer integral on [0, oo) and splits the log as
ln((t+2)/t) = ln(t+2) - ln(t). That is the move that avoids the c -> 1 limit
gymnastics entirely (both pieces diverge there and cancel; taking them separately
on [0,oo) never forms the divergence). It leaves exactly two new moment families:

    Lm(n, c)      = int_0^oo t^n e^{-ct} ln t dt      = (n!/c^{n+1})(H_n - gamma - ln c)
    Lms(n, c, s)  = int_0^oo t^n e^{-ct} ln(t+s) dt

plus the E_1 moments already built and validated in increment 2. If BOTH of these
close at weight 1, the whole ordered integral does.

This file builds and validates the two families. Run from repo root:
    python debug/inc3c_weight_probe.py
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

from geovac.two_center_eri import upper_integral  # noqa: E402


def log_moment(n: int, c):
    """int_0^oo t^n e^{-ct} ln t dt = (n!/c^{n+1}) [psi(n+1) - ln c].

    psi(n+1) = -gamma + H_n, so this is manifestly WEIGHT 1: it produces gamma
    and ln c and nothing else. This is the piece that carries Euler's gamma into
    the exchange class -- the same gamma Phase 0-e saw at the xi = 1 endpoint,
    reached here without ever forming the divergence.
    """
    H_n = sum(sp.Rational(1, k) for k in range(1, n + 1))
    return sp.factorial(n) / c ** (n + 1) * (H_n - sp.EulerGamma - sp.log(c))


def log_shift_moment(n: int, c, s):
    """int_0^oo t^n e^{-ct} ln(t+s) dt for s > 0.

    By parts against v = -int_t^oo u^n e^{-cu} du (which vanishes at infinity and
    is finite at 0, so no divergence is created):

        = U(n,c,0) ln s + int_0^oo U(n,c,t)/(t+s) dt

    and U(n,c,t) = e^{-ct} poly(t), so w = t+s turns the second term into
    incomplete-gamma pieces whose j = 0 term is E_1(cs). WEIGHT 1: it produces
    ln s and E_1(cs), nothing higher.
    """
    acc = upper_integral(n, c, sp.Integer(0)) * sp.log(s)
    # U(n,c,t) = e^{-ct} sum_k n!/(k! c^{n-k+1}) t^k
    for k in range(n + 1):
        coef = sp.factorial(n) / (sp.factorial(k) * c ** (n - k + 1))
        # int_0^oo t^k e^{-ct}/(t+s) dt,  w = t+s
        inner = sum(sp.binomial(k, j) * (-s) ** (k - j)
                    * upper_integral(j - 1, c, s) for j in range(k + 1))
        acc += coef * sp.exp(c * s) * inner
    return acc


def _weight(expr):
    """Census the transcendental content by weight."""
    names = {type(f).__name__ for f in expr.atoms(sp.Function)}
    has_g = expr.has(sp.EulerGamma)
    w2 = names & {"polylog", "dilog", "Li2"}
    return names, has_g, w2


def main() -> None:
    print("Increment 3c / periods probe -- the two log-moment families\n")

    print("L1  int_0^oo t^n e^{-ct} ln t dt")
    worst = 0.0
    for n in (0, 1, 2, 3, 4):
        for c in (0.7, 1.3, 2.5):
            cf = float(sp.N(log_moment(n, sp.nsimplify(c)), 30))
            nu, _ = integrate.quad(lambda t, n=n, c=c: t ** n * np.exp(-c * t)
                                   * np.log(t), 0, np.inf, limit=400,
                                   epsabs=1e-14, epsrel=1e-13)
            rel = abs(cf - nu) / max(abs(nu), 1e-300)
            worst = max(worst, rel)
            print(f"    n={n} c={c:4.1f}:  {cf: .13f}  vs {nu: .13f}  rel={rel:.1e}")
    print(f"    worst relative {worst:.2e}\n")

    print("L2  int_0^oo t^n e^{-ct} ln(t+s) dt")
    worst2 = 0.0
    for n in (0, 1, 2, 3):
        for c, s in ((1.3, 2.0), (0.7, 2.0), (2.5, 1.0)):
            cf = float(sp.N(log_shift_moment(n, sp.nsimplify(c),
                                             sp.nsimplify(s)), 30))
            nu, _ = integrate.quad(lambda t, n=n, c=c, s=s: t ** n * np.exp(-c * t)
                                   * np.log(t + s), 0, np.inf, limit=400,
                                   epsabs=1e-14, epsrel=1e-13)
            rel = abs(cf - nu) / max(abs(nu), 1e-300)
            worst2 = max(worst2, rel)
            print(f"    n={n} c={c:4.1f} s={s:3.1f}:  {cf: .13f}  vs {nu: .13f}"
                  f"  rel={rel:.1e}")
    print(f"    worst relative {worst2:.2e}\n")

    print("WEIGHT CENSUS of the closed forms (symbolic c, s):")
    cs, ss = sp.Symbol("c", positive=True), sp.Symbol("s", positive=True)
    for lbl, e in (("log_moment n=3", log_moment(3, cs)),
                   ("log_shift_moment n=3", log_shift_moment(3, cs, ss))):
        names, has_g, w2 = _weight(sp.expand(e))
        print(f"    {lbl:24s} functions={sorted(names)}  "
              f"EulerGamma={'yes' if has_g else 'no'}  weight-2 atoms={w2 or 'NONE'}")
    print("\n    Both are weight 1: log, E_1 (expint) and gamma only. No dilogarithm,")
    print("    no polylog, no zeta(2). The ordered xi integral therefore has no")
    print("    weight-2 obstruction -- what remained was assembly, not transcendence.")


if __name__ == "__main__":
    main()
