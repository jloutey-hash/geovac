"""Rung 3b -- principled weight-graded probe of the integrated collinear T2 against the
classical period/quasi-period ring of the central CM fiber (tau=i, rho=1/2), using the
high-precision value from the fast evaluator (debug/routeC_fast_evaluator.py).

v = T2_collinear = 0.3953557659017139... (dps34, N=40 vs 52 agree 5.1e-18 => ~17 digits).

Weight-graded basis from the tau=i fiber's period varpi=K(1/2)=Gamma(1/4)^2/(4 sqrt pi)
and quasi-period E2=E(1/2), plus pi (weights: 1=w0; pi,varpi,E2=w1; the products=w2).
A small-integer, precision-stable relation => T2 is a classical single-fiber
period/quasi-period combination (a genuine partial closure!).  NO low-height relation
=> T2 needs a genuinely-new weight-2 elliptic constant (elliptic dilogarithm / iterated
Eisenstein integral for Gamma(2)) = the frontier (obstacle i).

Honest precision caveat: ~17 digits + small maxcoeff => this is a BOUNDED exclusion of
low-height period-ring closures, not a proof of transcendental independence.
"""
from __future__ import annotations
import mpmath as mp
mp.mp.dps = 30

# fast-evaluator value (17 stable digits); more digits extensible via the evaluator.
V = mp.mpf('0.39535576590171392')

varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(mp.pi))   # = K(1/2)
E2 = mp.ellipe(mp.mpf('0.5'))                                  # E(m=1/2)
pi = mp.pi


def try_pslq(name, basis, names, maxcoeff, tol_exp=-14):
    rel = mp.pslq([V] + basis, tol=mp.mpf(10) ** tol_exp, maxcoeff=maxcoeff, maxsteps=10 ** 6)
    # A genuine CLOSURE of V requires the V-coefficient (rel[0]) nonzero AND small height.
    # rel[0]==0 is a basis-internal identity (dependent basis), NOT a closure of V.
    if rel is None:
        tag = "(none)"
    elif rel[0] == 0:
        tag = "(V-coeff 0 => basis-internal identity, NOT a closure of V)"
    elif max(abs(c) for c in rel) <= 40:
        tag = "<< SMALL, V-coeff nonzero => CANDIDATE CLOSURE"
    else:
        tag = "(high-height => not a closure)"
    print(f"  {name:28s}: {rel}\n      {tag}")
    return rel


def main():
    print("Rung 3b -- weight-graded period-ring probe of the integrated collinear T2\n")
    print(f"  V = {mp.nstr(V, 17)}  (fast evaluator, ~17 digits)")
    print(f"  varpi=K(1/2)={mp.nstr(varpi,12)}  E(1/2)={mp.nstr(E2,12)}\n")

    # (1) period ring {1, pi, varpi} up to weight 2
    try_pslq("period ring w<=2",
             [mp.mpf(1), pi, varpi, pi ** 2, varpi ** 2, varpi * pi],
             ["1", "pi", "varpi", "pi^2", "varpi^2", "varpi*pi"], maxcoeff=10 ** 4)

    # (2) + quasi-period E2 (full period+quasiperiod weight<=2 ring of the tau=i fiber)
    try_pslq("period+quasiperiod w<=2",
             [mp.mpf(1), pi, varpi, E2, pi ** 2, varpi ** 2, varpi * pi,
              E2 ** 2, E2 * varpi, E2 * pi],
             ["...+E2..."], maxcoeff=10 ** 3)

    # (3) classical polylog weight<=2 control (Catalan, ln2, Li2(1/2))
    try_pslq("classical polylog w<=2",
             [mp.mpf(1), pi, pi ** 2, mp.log(2), mp.log(2) ** 2, mp.catalan,
              mp.polylog(2, mp.mpf('0.5'))],
             ["classical"], maxcoeff=10 ** 4)

    print("\nRead: a SMALL-coefficient, precision-stable relation would be a partial closure;")
    print("high-height/none across all three (at ~17 digits, maxcoeff<=1e4) is a BOUNDED")
    print("exclusion => the integrated value is not a low-height classical / single-fiber")
    print("period-ring combination, consistent with a genuine Gamma(2) elliptic MMV whose")
    print("explicit form needs the iterated-Eisenstein basis (obstacle i = the frontier).")


if __name__ == '__main__':
    main()
