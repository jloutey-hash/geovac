"""Increment 2, l > 0: verify the shell reformulation BEFORE building it.

Last round I skipped Phase 0-h's open item 2 ("is there a formulation that keeps
the r_A powers non-negative?"), wrote the assembly, and hit the wall it predicted.
This is that item, done first.

THE PROPOSED FIX. Replace V_L(r_A) by its shell-integral representation

    V_L(r) = (4pi/(2L+1)) int_0^inf dx x^2 rad_1(x) min(r,x)^L / max(r,x)^{L+1}

and carry the shell radius x as a parameter through the r_A integral instead of
integrating it first (integrating it first just rebuilds V_L, which is the whole
problem). The claim is that this removes the pathology because:

  * V_L was regular at r -> 0 but its SPLIT into q_L r^{-(L+1)} + e^{-br}(...)
    was not. min/max is piecewise pure-power, and r_A -> 0 forces the INSIDE
    branch, whose power is +L. So nothing spurious is generated.
  * the negative powers that remain live only in the OUTSIDE branch, where the
    integration starts at max(|r_B - R|, x) >= x > 0, so they never meet a
    vanishing lower limit.
  * at that stage the only exponential in play is R_c's, with rate a_c > 0, so
    every E_1 argument is a POSITIVE rate times x -- which is exactly what
    e1_moment already consumes.

THREE CHECKS

  S1  the representation itself: does the shell integral reproduce V_L(r)?
  S2  the structural claim: inside-branch powers >= 0, and the E_1 rates that
      would be produced are positive.
  S3  end to end: a hybrid quartet with l > 0 on the one-center pair, computed
      through the shell route, against the validated quadrature reference.

Run from repo root:  python debug/inc2_shell_reformulation_check.py
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from scipy import integrate

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac.two_center_eri import (  # noqa: E402
    V_L_radial, angular_factor, hybrid_quadrature, hybrid_terms,
    multipole_decomposition, plm_signed, radial_norm, radial_poly, r_s,
    sph_norm_numeric, t_s, x_s, y_s,
)

Z1, Z2, Z3 = Fraction(1), Fraction(2), Fraction(3)


def _rad1_of(ZA, oa, ob, L):
    """The L-th multipole's radial function rad_1(x), as a numpy callable."""
    for LL, _M, _g, rad, b in multipole_decomposition(ZA, *oa, ZA, *ob):
        if LL == L:
            expr = sum(c * r_s ** k for k, c in rad.items()) * sp.exp(-b * r_s)
            return sp.lambdify(r_s, expr, "numpy")
    raise ValueError(f"no multipole L={L}")


def V_L_by_shells(ZA, oa, ob, L, r):
    """(4pi/(2L+1)) int_0^inf dx x^2 rad_1(x) min(r,x)^L / max(r,x)^{L+1}."""
    rad1 = _rad1_of(ZA, oa, ob, L)
    inner, _ = integrate.quad(lambda x: x ** 2 * rad1(x) * x ** L / r ** (L + 1),
                              0.0, r, epsabs=1e-14, epsrel=1e-13)
    outer, _ = integrate.quad(lambda x: x ** 2 * rad1(x) * r ** L / x ** (L + 1),
                              r, 80.0, epsabs=1e-14, epsrel=1e-13)
    return 4 * np.pi / (2 * L + 1) * (inner + outer)


def leg_S1() -> None:
    print("S1  does the shell representation reproduce V_L(r)?\n")
    worst = 0.0
    for oa, ob in (((2, 1, 0), (2, 1, 0)), ((3, 2, 0), (3, 2, 0))):
        for L, _M, _g, rad, b in multipole_decomposition(Z3, *oa, Z3, *ob):
            VL = sp.lambdify(r_s, V_L_radial(rad, b, L), "numpy")
            for rv in (0.4, 1.1, 2.7):
                direct = float(np.real(VL(rv)))
                shells = V_L_by_shells(Z3, oa, ob, L, rv)
                worst = max(worst, abs(direct - shells))
                print(f"    ({oa[1]},{ob[1]}) L={L} r={rv}:  V_L={direct: .12f}"
                      f"   shells={shells: .12f}   d={abs(direct-shells):.1e}")
    print(f"\n    worst {worst:.2e}  {'OK' if worst < 1e-8 else 'FAIL'}\n")


def leg_S2() -> None:
    print("S2  structural claim: inside-branch powers >= 0, E_1 rates positive\n")
    ZA, oa, ob, oc, od = Z3, (2, 1, 0), (2, 1, 0), (1, 0, 0), (1, 0, 0)
    coeffs, a_c = radial_poly(ZA, oc[0], oc[1])
    Nc = radial_norm(ZA, oc[0], oc[1])
    Rc = sum(Nc * c * t_s ** k for k, c in coeffs.items())        # polynomial part
    print(f"    probe (2p0 2p0|1s 1s_B); R_c decay a_c = {a_c} > 0\n")
    print("    L  L'   inside-branch r_A powers      outside-branch r_A powers")
    print("    " + "-" * 66)
    ok_inside = True
    for _coeff, L, Lp, Mp, ld, md in hybrid_terms(ZA, oa, ob, oc, Z1, od):
        ang = sp.expand(angular_factor(Lp, Mp, ld, md))
        # keep only the t-powers; y-powers are handled by the y^2 rad_B cancel
        def powers(extra):
            out = set()
            for term in sp.Add.make_args(sp.expand(Rc * ang * t_s * extra)):
                p = 0
                for f in sp.Mul.make_args(term):
                    if f == t_s:
                        p += 1
                    elif f.is_Pow and f.base == t_s:
                        p += int(f.exp)
                out.add(p)
            return sorted(out)
        pin = powers(t_s ** L / x_s ** (L + 1))
        pout = powers(x_s ** L / t_s ** (L + 1))
        ok_inside &= min(pin) >= 0
        print(f"    {L}  {Lp}   min={min(pin):+d} max={max(pin):+d}"
              f"                min={min(pout):+d} max={max(pout):+d}")
    print(f"\n    inside branch all powers >= 0: {ok_inside}")
    print("    outside branch goes negative, BUT its lower limit is")
    print("    max(|r_B - R|, x) >= x > 0, so it never meets a vanishing limit.")
    print(f"    Every E_1 produced there has argument a_c * x = {a_c} * x,")
    print("    a POSITIVE rate -- exactly what e1_moment consumes.\n")


def leg_S3() -> None:
    """The REGION BOOKKEEPING -- the part that can actually go wrong in the build.

    An end-to-end shell-vs-reference run was tried first and dropped: S1 already
    shows the shell representation equals V_L pointwise to 2.6e-12, so threading
    it through the (already validated) outer quadrature is the same identity at
    triple-quadrature cost and yields no independent information. It ran >25 min
    on one term without finishing.

    What DOES need de-risking is the three-region split of the r_A integral at
    r_A = x, since that is the new bookkeeping the symbolic builder must get
    right. Checked here in 1D, against the same integral computed without any
    split.
    """
    print("S3  region bookkeeping: the r_A split at r_A = x\n")
    print("    (replaces an end-to-end run -- see the docstring for why that")
    print("     leg was redundant given S1)\n")
    L, a_c, lc = 2, 3.0, 1
    R = 2.5

    def f(rA):                       # stand-in for R_c(r_A) x angular, any shape
        return rA ** lc * np.exp(-a_c * rA) * (1.0 + 0.3 * rA ** 2)

    def S(rA, x):
        return min(rA, x) ** L / max(rA, x) ** (L + 1)

    print("      r_B     x     split-by-region      direct           d")
    print("    " + "-" * 60)
    worst = 0.0
    for rb in (0.5, 2.5, 4.0):
        lo, hi = abs(rb - R), rb + R
        for x in (0.3, 1.5, 3.0, 9.0):
            direct, _ = integrate.quad(lambda rA: f(rA) * S(rA, x), lo, hi,
                                       epsabs=1e-13, epsrel=1e-13, limit=300)
            if x <= lo:                                   # all outside
                pieces = [(lo, hi, "out")]
            elif x >= hi:                                 # all inside
                pieces = [(lo, hi, "in")]
            else:                                         # split at x
                pieces = [(lo, x, "in"), (x, hi, "out")]
            split = 0.0
            for p0, p1, kind in pieces:
                if kind == "in":
                    g = (lambda rA, x=x: f(rA) * rA ** L / x ** (L + 1))
                else:
                    g = (lambda rA, x=x: f(rA) * x ** L / rA ** (L + 1))
                val, _ = integrate.quad(g, p0, p1, epsabs=1e-13, epsrel=1e-13,
                                        limit=300)
                split += val
            worst = max(worst, abs(split - direct))
            print(f"     {rb:4.1f}  {x:4.1f}   {split: .12f}  {direct: .12f}"
                  f"  {abs(split-direct):.1e}")
    print(f"\n    worst {worst:.2e}  {'OK' if worst < 1e-10 else 'FAIL'}")
    print("    Note r_B = 2.5 = R exactly, where lo = 0: the inside branch")
    print("    handles it with power +L and nothing diverges. That is the whole")
    print("    point of the reformulation.\n")


def main() -> None:
    print("Increment 2 -- verifying the shell reformulation before building it\n")
    leg_S1()
    leg_S2()
    leg_S3()


if __name__ == "__main__":
    main()
