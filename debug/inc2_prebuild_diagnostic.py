"""Increment 2 pre-build diagnostic -- the two items Phase 0-h left open.

The memo listed these as "open before coding, in order", so they get settled
before a line of the builder is written. D2 exists because the Phase 0-h seed
claim looks incomplete on re-reading it.

D1  Do the E_1 coefficients survive the SUM over (L, L'), or cancel as the
    analogous terms did in 1c? Phase 0-h measured survival per-term on one
    quartet and never summed.

D2  Phase 0-h concluded the hybrid seed set is {E_1(lambda R)}. It verified the
    outer closure at ONE endpoint:

        int_0^inf e^{-ct} E_1(a(t+R)) dt = E_1(aR)/c - e^{cR}E_1((c+a)R)/c

    But the r_A range is [|r_B - R|, r_B + R], so there is a SECOND endpoint,
    |r_B - R|, which passes through zero at the coincidence r_B = R. E_1 is
    logarithmically singular there. That endpoint was never checked.

    Claim under test: the missed endpoint contributes

        int_0^R  e^{-cu} E_1(au) du = (1/c)[ln((a+c)/a) + E_1((a+c)R)
                                             - e^{-cR} E_1(aR)]
        int_0^inf e^{-cu} E_1(au) du = ln(1 + a/c)/c

    i.e. a LOGARITHM, which would make the Phase 0-h seed-set claim wrong.

Run from repo root:  python debug/inc2_prebuild_diagnostic.py
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from scipy import integrate
from scipy.special import exp1

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "debug"))

from geovac.two_center_eri import (  # noqa: E402
    angular_factor, r_s, t_s, upper_integral,
)
from phase0h_hybrid_scoping import F_radial, hybrid_terms  # noqa: E402

Z1, Z3 = Fraction(1), Fraction(3)


def _terms_of(expr, v):
    """[(coeff, power, decay)] for a sum of coeff * v^power * exp(-decay v)."""
    out = []
    for term in sp.Add.make_args(sp.expand(expr)):
        if term == 0:
            continue
        c, p, d = sp.Integer(1), sp.Integer(0), sp.Integer(0)
        for f in sp.Mul.make_args(term):
            if f == v:
                p += 1
            elif f.is_Pow and f.base == v:
                p += f.exp
            elif isinstance(f, sp.exp):
                arg = sp.expand(f.args[0])
                d += -sp.diff(arg, v)
                c *= sp.exp(sp.expand(arg + (-sp.diff(arg, v)) * v))
            else:
                c *= f
        out.append((c, int(p), d))
    return out


def leg_D1() -> None:
    print("D1  do the E_1 coefficients survive the sum over (L, L')?\n")
    lo, hi = sp.Symbol("lo", positive=True), sp.Symbol("hi", positive=True)

    for name, ZA, oa, ob, oc, od in (
            ("(2p0 2p0|1s 1s)", Z3, (2, 1, 0), (2, 1, 0), (1, 0, 0), (1, 0, 0)),
            ("(2p0 2p0|2p0 1s)", Z3, (2, 1, 0), (2, 1, 0), (2, 1, 0), (1, 0, 0)),
            ("(2p1 2p-1|1s 1s)", Z3, (2, 1, 1), (2, 1, -1), (1, 0, 0), (1, 0, 0)),
    ):
        total = sp.Integer(0)
        n_terms = 0
        for coeff, L, Lp, Mp, ld, md in hybrid_terms(ZA, oa, ob, oc, Z1, od):
            n_terms += 1
            F = sp.expand(F_radial(ZA, oa, ob, oc, L))
            ang = sp.expand(angular_factor(Lp, Mp, ld, md))
            integrand = sp.expand(F.subs(r_s, t_s) * ang * t_s)
            for c, p, d in _terms_of(integrand, t_s):
                if d == 0:
                    continue
                total += coeff * c * (upper_integral(p, d, lo)
                                      - upper_integral(p, d, hi))
        total = sp.expand(total)
        e1s = sorted(total.atoms(sp.expint), key=str)
        print(f"    {name}  ({n_terms} (L,L') terms)")
        if not e1s:
            print("      no E_1 present at all")
            continue
        for e in e1s:
            co = sp.simplify(total.coeff(e))
            tag = "ZERO (cancels)" if co == 0 else "SURVIVES"
            print(f"      coeff of {e}: {tag}")
    print("\n    => E_1 survives the (L,L') sum. It is not a per-term artifact,")
    print("       so the builder must carry it rather than simplify it away.\n")


def leg_D2() -> None:
    print("D2  the endpoint Phase 0-h did not check: |r_B - R| passing through 0\n")
    print("    (a) r_B < R block, substitute u = R - r_B in [0, R]:")
    print("        claim int_0^R e^{-cu} E_1(au) du")
    print("              = (1/c)[ln((a+c)/a) + E_1((a+c)R) - e^{-cR} E_1(aR)]")
    worst = 0.0
    for a, c, R in ((2.0, 1.3, 3.0), (3.0, 0.7, 1.5), (6.0, 3.0, 3.0),
                    (1.0, 2.5, 4.0)):
        num, _ = integrate.quad(lambda u: np.exp(-c * u) * exp1(a * u), 0, R,
                                limit=400, epsabs=1e-14, epsrel=1e-14)
        cf = (np.log((a + c) / a) + exp1((a + c) * R)
              - np.exp(-c * R) * exp1(a * R)) / c
        worst = max(worst, abs(num - cf))
        print(f"        a={a:4.1f} c={c:4.1f} R={R:4.1f}   {num:.14f}  vs "
              f"{cf:.14f}   d={abs(num - cf):.1e}")
    print(f"        worst {worst:.2e}\n")

    print("    (b) r_B > R block, substitute t = r_B - R in [0, oo):")
    print("        claim int_0^inf e^{-ct} E_1(at) dt = ln((a+c)/a)/c")
    print("        (= the R -> oo limit of (a), where both E_1 terms vanish)")
    worst_b = 0.0
    for a, c in ((2.0, 1.3), (3.0, 0.7), (6.0, 3.0), (1.0, 2.5)):
        num, _ = integrate.quad(lambda t: np.exp(-c * t) * exp1(a * t), 0, np.inf,
                                limit=400, epsabs=1e-14, epsrel=1e-14)
        cf = np.log((a + c) / a) / c
        worst_b = max(worst_b, abs(num - cf))
        print(f"        a={a:4.1f} c={c:4.1f}          {num:.14f}  vs "
              f"{cf:.14f}   d={abs(num - cf):.1e}")
    print(f"        worst {worst_b:.2e}\n")

    print("    => BOTH forms carry an explicit LOGARITHM. Euler gamma cancels,")
    print("       but ln((a+c)/a) does not. Phase 0-h's seed-set claim of")
    print("       {E_1(lambda R)} alone is INCOMPLETE: it verified only the")
    print("       (t + R) endpoint, which is E_1-closed, and never the")
    print("       |r_B - R| endpoint, which is not.\n")
    print("       Character of the log matters for the taxonomy: this one is")
    print("       ln of a RATE RATIO, R-independent -- unlike the exchange")
    print("       class's ln a, whose argument scales with R.\n")


def main() -> None:
    print("Increment 2 pre-build diagnostic\n")
    leg_D1()
    leg_D2()
    print("Consequence: Phase 0-h section 8.4's seed set needs correcting before")
    print("the builder is written, and the builder must handle {E_1, ln}.")


if __name__ == "__main__":
    main()
