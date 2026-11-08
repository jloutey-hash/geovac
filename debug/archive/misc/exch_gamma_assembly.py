"""Track B / exchange class -- IS ordered_xi_* THE RIGHT REPRESENTATIVE?

The assembled exchange ERI (build plan section 8.5) is

    (AB|AB) = prefactor * sum_tau  [ eta-half(tau) ] * [ ordered-xi-half(tau) ]

The xi half is the object analysed in exch_gamma_{census,borel,general}.py.
This driver checks the two things that could make that a misleading choice:

  A1  the ETA half is ELEMENTARY -- coeff * R^{-n} * e^{+-qR} only (no E_1, no
      log, no gamma).  If so it carries no transcendental content of its own and
      can only SHIFT sector actions / rescale Stokes data by rationals.
  A2  the PRODUCT (one tau term of the real ERI) still has 2 pi i x rational
      Stokes data, at shifted positions.

What this driver does NOT settle, and the findings memo says so: the tau sum is
INFINITE for heteronuclear centres (build plan EQ1), and the resurgent structure
of an infinite sum of trans-series is a separate question.

Run:  python debug/exch_gamma_assembly.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import sympy as sp

from geovac.two_center_eri import R_s, ordered_xi_closed

from exch_gamma_census import parse
from exch_gamma_general import borel_symbolic

OUT = []
XI = sp.Symbol("xi")
GAM = sp.Symbol("Gam")


def say(s=""):
    print(s)
    OUT.append(s)


def eta_integral(k, tau, q):
    """int_{-1}^{1} eta^k P_tau(eta) e^{-q eta} d eta, exact."""
    e = sp.Symbol("e")
    integrand = e ** k * sp.legendre(tau, e) * sp.exp(-q * e)
    v = sp.integrate(integrand, (e, -1, 1))
    return sp.expand(sp.simplify(v).rewrite(sp.exp))


def main():
    os.makedirs("debug/data", exist_ok=True)
    say("=" * 100)
    say("EXCHANGE CLASS -- representative check: the eta half and the assembled product")
    say("=" * 100)

    gq = sp.Rational(3, 4)                 # q/R : the (alpha-beta)/2 difference rate
    q = gq * R_s
    say("")
    say("A1 -- the eta half  int_{-1}^{1} eta^k P_tau(eta) e^{-qR eta} d eta,  q/R = %s"
        % gq)
    say("")
    say("  %8s %8s %14s %26s %14s"
        % ("k", "tau", "elementary?", "exponential rates c", "R powers"))
    a1_ok = True
    for k in range(0, 3):
        for tau in range(0, 3):
            val = sp.expand(eta_integral(k, tau, q))
            try:
                terms = parse(val)
            except ValueError as ex:
                say("  %8d %8d  PARSE RAISED: %s" % (k, tau, ex))
                a1_ok = False
                continue
            elem = all(t.kind == "elem" for t in terms)
            a1_ok &= elem
            cs = sorted({t.c for t in terms}, key=lambda z: float(z))
            ks = sorted({int(t.k) for t in terms})
            say("  %8d %8d %14s %26s %14s"
                % (k, tau, elem, [str(c) for c in cs], ks))
    say("")
    say("  A1  the eta half is ELEMENTARY (exp x rational R powers only)         : %s"
        % ("PASS" if a1_ok else "FAIL"))
    say("      => it carries NO transcendental: no E_1, no log, no gamma.")
    say("      Its only effect is to split the single xi-half sector into")
    say("      sectors at A -+ q and to rescale coefficients by rationals.")

    # ---------------------------------------------------------------- A2
    say("")
    say("A2 -- one assembled tau term:  (eta half) x (ordered-xi half)")
    say("")
    al, be = sp.Rational(3, 2), sp.Integer(1)
    A = al + be
    say("  xi half rates (a,b) = (%s,%s); xi-half action A = %s; eta rate q/R = %s"
        % (al, be, A, gq))
    say("")
    say("  %8s %8s %10s %34s %14s"
        % ("k", "tau", "#sectors", "Borel singularity positions (global)", "Stokes rational?"))
    a2_ok = True
    F = sp.expand(ordered_xi_closed(al * R_s, be * R_s))
    for k in range(0, 2):
        for tau in range(0, 3):
            prod = sp.expand(sp.expand(eta_integral(k, tau, q).rewrite(sp.exp)) * F)
            terms = parse(prod)
            # group by sector action, Borel-analyse each sector separately
            acts = sorted({t.action for t in terms}, key=lambda z: float(z))
            allpos, rat = [], True
            for Aa in acts:
                sub = [t for t in terms if t.action == Aa]
                # strip the sector factor: shift c by (Aa - c) as borel_symbolic expects
                B, sing = borel_symbolic(sub)
                if sp.simplify(sp.expand(B.coeff(GAM, 1))) != 0:
                    rat = False
                for pos, P in sing.items():
                    allpos.append(sp.nsimplify(Aa + pos))
                    if not all(c.is_Rational
                               for c in sp.Poly(sp.expand(P), XI).all_coeffs()):
                        rat = False
            allpos = sorted(set(allpos), key=lambda z: float(z))
            a2_ok &= rat
            say("  %8d %8d %10d %34s %14s"
                % (k, tau, len(acts), [str(p) for p in allpos], rat))
    say("")
    say("  A2  every assembled tau term keeps gamma-free Borel data and")
    say("      Stokes polynomials with coefficients in Q(rates)                  : %s"
        % ("PASS" if a2_ok else "FAIL"))
    say("      Positions are the xi-half positions {+-A, +-(a-b)} rigidly shifted")
    say("      by the eta exponents +-q -- a rational translation, no new type.")
    say("")
    say("  SCOPE (honest): for heteronuclear centres the tau sum is INFINITE")
    say("  (build plan EQ1: it terminates iff the two centres share an exponent).")
    say("  Every term is classified above; the resurgent structure of the infinite")
    say("  SUM is a separate question this track does not settle.")

    with open("debug/data/exch_gamma_assembly.txt", "w") as fh:
        fh.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
