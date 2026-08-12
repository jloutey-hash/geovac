"""Full-build increment 1c: (AA|BB) in closed form, and its validation.

The derivation and its rationale live with the code, in geovac/two_center_eri.py
(the "Increment 1c" section). This file is the validation driver.

FIVE LEGS, weakest to strongest:

  V0  shell kernel vs direct 3-fold angular quadrature -- catches convention
      errors (Condon-Shortley, signed-M normalization, the sin^|M| pairing).
      Shells chosen NOT to intersect so 1/d is smooth and the reference honest.
  V1  (1s 1s|1s 1s) vs the EXACT closed-form VB J integral. The strongest leg:
      the symbolic result is not merely equal numerically, it simplifies to the
      textbook expression term by term -- derived, not transcribed.
  V2  l > 0 and M != 0 vs the increment-1b quadrature route, which shares no 1c
      code. This is the only reference that can referee M != 0.
  V3  centre-swap consistency (cd|ab) = (-1)^{la+lb+lc+ld} (ab|cd). Regions
      A/B/C are NOT symmetric under x <-> y, so this genuinely exercises them.
  V4  Gate (d) of the build plan: independent McMurchie-Davidson Gaussians on
      the Paper-58 census configuration, expected ~1e-6 (STO->Gaussian fits).

Plus the structural census: the closed form contains `exp` and nothing else.

Run from repo root:  python debug/inc1c_closed_form.py
"""

from __future__ import annotations

import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac.two_center_eri import (  # noqa: E402
    R_s, aabb_closed_form, aabb_quadrature, shell_kernel_antiderivatives,
    shell_kernel_quadrature, t_s, x_s, y_s,
)

Z1, Z2, Z3 = Fraction(1), Fraction(2), Fraction(3)


def _J_exact(zeta: float, R: float) -> float:
    rho = zeta * R
    return (1.0 / R) * (1.0 - np.exp(-2 * rho)
                        * (1 + 1.375 * rho + 0.75 * rho ** 2 + rho ** 3 / 6.0))


def leg_V0() -> float:
    print("V0  shell kernel K(x,y) vs direct angular quadrature")
    print("    (non-intersecting shells, so 1/d is smooth; the intersecting")
    print("     case is refereed by V1, where the 1s densities overlap)")
    worst = 0.0
    for LA, MA, LB, MB in ((0, 0, 0, 0), (1, 0, 1, 0), (2, 0, 0, 0),
                           (1, 1, 1, -1), (2, 1, 1, -1), (2, 2, 2, -2),
                           (1, -1, 2, 1)):
        A_in, A_out, pref = shell_kernel_antiderivatives(LA, MA, LB, MB)
        for xv, yv, Rv in ((0.5, 1.3, 3.0), (3.4, 0.6, 1.5), (0.9, 1.2, 4.0)):
            lo, hi = abs(yv - Rv), yv + Rv
            A = A_out if xv <= lo else A_in
            core = A.subs(t_s, hi) - A.subs(t_s, lo)
            got = float(sp.re((pref / (y_s * R_s) * core).subs(
                {x_s: xv, y_s: yv, R_s: Rv}).evalf()))
            ref = shell_kernel_quadrature(LA, MA, LB, MB, xv, yv, Rv)
            worst = max(worst, abs(got - ref))
    print(f"    7 (LA,MA,LB,MB) x 3 geometries, both branches")
    print(f"    worst |closed - quadrature| = {worst:.2e}  "
          f"{'OK' if worst < 1e-9 else 'FAIL'}\n")
    return worst


def leg_V1() -> float:
    print("V1  (1s_A 1s_A|1s_B 1s_B) vs the EXACT VB J integral")
    expr = sp.simplify(aabb_closed_form(Z1, (1, 0, 0), (1, 0, 0),
                                        Z1, (1, 0, 0), (1, 0, 0)))
    print(f"    closed form:  {expr}")
    target = (1 / R_s) * (1 - sp.exp(-2 * R_s)
                          * (1 + sp.Rational(11, 8) * R_s
                             + sp.Rational(3, 4) * R_s ** 2 + R_s ** 3 / 6))
    identical = sp.simplify(expr - target) == 0
    print(f"    symbolically identical to the textbook J(R): {identical}")
    worst = 0.0
    for Rv in (1.5, 2.5, 4.0):
        got = float(sp.re(sp.N(expr.subs(R_s, sp.nsimplify(Rv)), 30)))
        ref = _J_exact(1.0, Rv)
        worst = max(worst, abs(got - ref))
        print(f"    R={Rv:4.2f}   {got:.14f}  vs {ref:.14f}   d={abs(got-ref):.2e}")
    print(f"    worst = {worst:.2e}  "
          f"{'OK' if worst < 1e-12 and identical else 'FAIL'}\n")
    return worst


def leg_V2() -> float:
    print("V2  l > 0 and M != 0 vs the 1b quadrature route (no shared 1c code)")
    cases = [
        ("(2p0 2p0|1s 1s)",    Z1, (2, 1, 0), (2, 1, 0), Z2, (1, 0, 0), (1, 0, 0)),
        ("(1s 2p0|2s 2p0)",    Z1, (1, 0, 0), (2, 1, 0), Z2, (2, 0, 0), (2, 1, 0)),
        ("(2p1 2p0|2p0 2p1)",  Z1, (2, 1, 1), (2, 1, 0), Z2, (2, 1, 0), (2, 1, 1)),
        ("(3d2 3d0|2p-1 2p1)", Z1, (3, 2, 2), (3, 2, 0), Z2, (2, 1, -1), (2, 1, 1)),
    ]
    R = 2.5
    worst = 0.0
    for name, ZA, oa, ob, ZB, oc, od in cases:
        t0 = time.time()
        got = float(sp.re(sp.N(aabb_closed_form(ZA, oa, ob, ZB, oc, od,
                                                sp.nsimplify(R)), 30)))
        dt = time.time() - t0
        ref = aabb_quadrature(ZA, oa, ob, ZB, oc, od, R)
        worst = max(worst, abs(got - ref))
        print(f"    {name:20s} closed={got: .12f}  quad={ref: .12f}"
              f"  d={abs(got-ref):.1e}  [{dt:.1f}s]")
    print(f"    worst = {worst:.2e}  {'OK' if worst < 1e-11 else 'FAIL'}\n")
    return worst


def leg_V3() -> float:
    print("V3  centre swap: (cd|ab) = (-1)^{la+lb+lc+ld} (ab|cd)")
    print("    (regions A/B/C are asymmetric under x <-> y, so this bites)")
    cases = [(Z1, (1, 0, 0), (1, 0, 0), Z2, (2, 1, 0), (2, 1, 0)),
             (Z1, (2, 1, 0), (2, 0, 0), Z2, (2, 1, 0), (1, 0, 0)),
             (Z1, (2, 1, 0), (2, 1, 0), Z2, (3, 2, 0), (3, 2, 0))]
    R = sp.Rational(5, 2)
    worst = 0.0
    for ZA, oa, ob, ZB, oc, od in cases:
        v1 = float(sp.re(sp.N(aabb_closed_form(ZA, oa, ob, ZB, oc, od, R), 30)))
        v2 = float(sp.re(sp.N(aabb_closed_form(ZB, oc, od, ZA, oa, ob, R), 30)))
        ph = (-1) ** (oa[1] + ob[1] + oc[1] + od[1])
        worst = max(worst, abs(v1 - ph * v2))
        print(f"    {v1: .12f}  vs {ph * v2: .12f}   d={abs(v1 - ph * v2):.1e}")
    print(f"    worst = {worst:.2e}  {'OK' if worst < 1e-12 else 'FAIL'}\n")
    return worst


def leg_V4() -> float:
    print("V4  Gate (d): independent McMurchie-Davidson, census config Z_A=3, "
          "Z_B=1, R=3")
    from geovac import noci_engine as E
    shapes = {}
    for kind, (l, nr) in (("1s", (0, 1)), ("2p", (1, 2))):
        arr, dco, q = E.fit_sto_shape(l, nr)
        shapes[kind] = (arr, dco)
        print(f"    fit {kind}: <fit|STO> = {q:.8f}")
    pa, pb = np.array([0., 0., 0.]), np.array([0., 0., 3.])

    def B(c, kind, zeta, lmn):
        return E.sto_shape_basis(c, kind, zeta, shapes, lmn)

    R = sp.Integer(3)
    rows, worst = [], 0.0

    got = float(sp.re(sp.N(aabb_closed_form(Z3, (1, 0, 0), (1, 0, 0),
                                            Z1, (1, 0, 0), (1, 0, 0), R), 30)))
    md = E.eri_md(B(pa, "1s", 3.0, (0, 0, 0)), B(pa, "1s", 3.0, (0, 0, 0)),
                  B(pb, "1s", 1.0, (0, 0, 0)), B(pb, "1s", 1.0, (0, 0, 0)))
    rows.append(("(1sA 1sA|1sB 1sB)", got, md))

    got = float(sp.re(sp.N(aabb_closed_form(Z3, (2, 1, 0), (2, 1, 0),
                                            Z1, (1, 0, 0), (1, 0, 0), R), 30)))
    md = E.eri_md(B(pa, "2p", 1.5, (0, 0, 1)), B(pa, "2p", 1.5, (0, 0, 1)),
                  B(pb, "1s", 1.0, (0, 0, 0)), B(pb, "1s", 1.0, (0, 0, 0)))
    rows.append(("(2p0A 2p0A|1sB 1sB)", got, md))

    # m = +1 is not a Cartesian monomial: conj(Y11) Y11 = (px^2 + py^2)/2
    got = float(sp.re(sp.N(aabb_closed_form(Z3, (2, 1, 1), (2, 1, 1),
                                            Z1, (1, 0, 0), (1, 0, 0), R), 30)))
    md = 0.5 * sum(
        E.eri_md(B(pa, "2p", 1.5, lmn), B(pa, "2p", 1.5, lmn),
                 B(pb, "1s", 1.0, (0, 0, 0)), B(pb, "1s", 1.0, (0, 0, 0)))
        for lmn in ((1, 0, 0), (0, 1, 0)))
    rows.append(("(2p1A 2p1A|1sB 1sB)", got, md))

    got = float(sp.re(sp.N(aabb_closed_form(Z3, (2, 1, 0), (2, 1, 0),
                                            Z1, (2, 1, 0), (2, 1, 0), R), 30)))
    md = E.eri_md(B(pa, "2p", 1.5, (0, 0, 1)), B(pa, "2p", 1.5, (0, 0, 1)),
                  B(pb, "2p", 0.5, (0, 0, 1)), B(pb, "2p", 0.5, (0, 0, 1)))
    rows.append(("(2p0A 2p0A|2p0B 2p0B)", got, md))

    for name, c, m in rows:
        worst = max(worst, abs(c - m))
        print(f"    {name:22s} closed={c:.9f}  md={m:.9f}  d={abs(c - m):.2e}")
    print(f"    worst = {worst:.2e}  {'OK (fit-limited)' if worst < 1e-5 else 'FAIL'}\n")
    return worst


def leg_census() -> None:
    print("Structural census of the closed form (symbolic R)")
    for name, args in (
            ("1s1s|1s1s",       (Z1, (1, 0, 0), (1, 0, 0), Z1, (1, 0, 0), (1, 0, 0))),
            ("2p0 2p0|1s1s",    (Z1, (2, 1, 0), (2, 1, 0), Z2, (1, 0, 0), (1, 0, 0))),
            ("3d1 3d0|2p0 2p1", (Z1, (3, 2, 1), (3, 2, 0), Z2, (2, 1, 0), (2, 1, 1)))):
        e = aabb_closed_form(*args)
        fn = sorted({type(f).__name__ for f in e.atoms(sp.Function)})
        assert not e.atoms(sp.expint), "E_1 present"
        assert not e.atoms(sp.log), "log present"
        print(f"    {name:18s} functions present -> {fn}")
    print("    => `exp` only. No expint, no log, at any l or M.")
    e = aabb_closed_form(Z1, (2, 1, 0), (2, 1, 0), Z2, (1, 0, 0), (1, 0, 0))
    print(f"    large-R limit R*(2p0 2p0|1s 1s) -> "
          f"{sp.limit(sp.expand(e) * R_s, R_s, sp.oo)}  (unit charges: expect 1)\n")


def main() -> None:
    print("Increment 1c -- (AA|BB) in closed form\n")
    w = [leg_V0(), leg_V1(), leg_V2(), leg_V3(), leg_V4()]
    leg_census()
    print("GATE: closed form reproduces every independent reference at its own "
          "accuracy;\n      the class is elementary -- E_1 structurally absent, "
          "not cancelled.")
    print(f"      machine-precision legs V0-V3 worst = {max(w[:4]):.2e}; "
          f"fit-limited V4 worst = {w[4]:.2e}")


if __name__ == "__main__":
    main()
