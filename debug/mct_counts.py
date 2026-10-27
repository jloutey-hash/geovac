"""Phase 2 of docs/predictions/2026-10-27_master_count_tabulation.md (MC0 C2-C5, MC1, MC2, MC2b, MC4, MC6).

Critical-point (Euler-characteristic) counts.
  Route X : omega = d(-Dhat x) + sum nu_i dlog g_i  on the x-line (exponential twist, D-module count).
  Route M : omega = sum nu_i dlog g_i on C*, generic exponent on x (Mellin dual, AFST Euler integral).
  Route F : conic quotient of the two-offset fibre, omega = d(-D1 x(t) - D2 w(t)) + sum_branch nu_i dlog.
  Route LP: Lee-Pomeransky critical points of sum nu_i log x_i - s log G (controls C2, C3), Groebner mod p,
            counting function copied from lit_scan_euler_kite.py.
Count = number of distinct zeros of omega (finite chart) that are not removed points; exact over Q and mod p0.
usage: python debug/mct_counts.py [controls|slice|fibre|all]
"""
from __future__ import annotations

import itertools
import random
import sys
import time
from fractions import Fraction as Fr

import sympy as sp

P0 = 2**61 - 1
NU = [sp.Rational(3, 7), sp.Rational(5, 11), sp.Rational(2, 13), sp.Rational(7, 17), sp.Rational(4, 19)]
S_LP = sp.Rational(9, 23)
x, tt = sp.symbols("x t")


def count_1d(omega: sp.Expr, var: sp.Symbol, removed: list) -> dict:
    """Distinct zeros of the rational 1-form coefficient omega, excluding removed points (polys in var)."""
    num, den = sp.fraction(sp.cancel(sp.together(omega)))
    out = {}
    for label, mod in (("Q", None), ("p0", P0)):
        kw = {} if mod is None else {"modulus": mod}
        N = sp.Poly(sp.Poly(num, var).clear_denoms()[1].as_expr(), var, **kw)
        Nsq = sp.Poly(sp.quo(N, sp.gcd(N, N.diff(var))), var, **kw)
        for g in removed:
            G = sp.Poly(sp.Poly(g, var).clear_denoms()[1].as_expr(), var, **kw)
            if G.degree() <= 0:
                continue
            while True:
                h = sp.gcd(Nsq, G)
                if h.degree() <= 0:
                    break
                Nsq = sp.Poly(sp.quo(Nsq, h), var, **kw)
        out[label] = Nsq.degree()
        out[label + "_degN"] = N.degree()
    return out


def omega_X(dhat, facs):
    """facs: list of (poly, exponent); exponent 0 entries are dropped (point not removed)."""
    facs = [(g, e) for g, e in facs if e != 0]
    om = -dhat + sum(e * sp.diff(g, x) / g for g, e in facs)
    return om, [g for g, _ in facs]


def omega_M(facs):
    facs = [(g, e) for g, e in facs if e != 0]
    om = sum(e * sp.diff(g, x) / g for g, e in facs)
    return om, [g for g, _ in facs]


def quad(sig):
    return sig * x**2 + 1 - sig


def slice_X(dhat, sig, n1, n2, generic: bool):
    ex1 = 1 - n1
    if generic:
        facs = [(x - 1, NU[0]), (x + 1, NU[1]), (quad(sig), NU[2]), (x, NU[3] if ex1 != 0 else 0)]
    else:
        facs = [(x - 1, sp.Rational(-1, 2)), (x + 1, sp.Rational(-1, 2)), (quad(sig), sp.Rational(-n2, 2)),
                (x, ex1)]
    if sig == 0:
        facs = [f for f in facs if f[0] != quad(sig)]
    om, rem = omega_X(dhat, facs)
    return count_1d(om, x, rem)


def slice_M(sig, n1, n2, generic: bool):
    ex_x = NU[4] + 1 - n1  # generic Mellin exponent -s + 1 - n1 with s generic
    if generic:
        facs = [(x - 1, NU[0]), (x + 1, NU[1]), (quad(sig), NU[2]), (x, ex_x)]
    else:
        facs = [(x - 1, sp.Rational(-1, 2)), (x + 1, sp.Rational(-1, 2)), (quad(sig), sp.Rational(-n2, 2)),
                (x, ex_x)]
    om, rem = omega_M(facs)
    return count_1d(om, x, rem)


# ---------------- Route LP (controls) ----------------
def lp_count(G, xs, p=32003):
    n = len(xs)
    t = sp.Symbol("tsat")
    eqs = [sp.expand(NU[i] * G - S_LP * xs[i] * sp.diff(G, xs[i])) for i in range(n)]
    eqs.append(t * sp.Mul(*xs) * G - 1)
    eqs = [sp.together(e).as_numer_denom()[0] for e in eqs]
    gb = sp.groebner(eqs, t, *xs, order="grevlex", modulus=p)
    lead = [sp.Poly(e, t, *xs, modulus=p).monoms(order="grevlex")[0] for e in gb.exprs]
    nv = n + 1
    for i in range(nv):
        if not any(all(m[j] == 0 for j in range(nv) if j != i) for m in lead):
            return "not zero-dimensional"
    bound = [min(m[i] for m in lead if all(m[j] == 0 for j in range(nv) if j != i)) for i in range(nv)]
    cnt = 0
    for mono in itertools.product(*[range(bb) for bb in bound]):
        if not any(all(mono[j] >= m[j] for j in range(nv)) for m in lead):
            cnt += 1
    return cnt


def controls():
    print("== MC0 controls (C1 kite run separately via lit_scan_euler_kite.py)", flush=True)
    res = {}
    # C2 bubble: Euclidean, unequal masses, generic P^2
    x1, x2, x3 = sp.symbols("x1 x2 x3")
    m1, m2, m3, P2 = sp.Rational(1), sp.Rational(9, 4), sp.Rational(49, 9), sp.Rational(37, 101)
    U = x1 + x2
    G = sp.expand(U + x1 * x2 * P2 + U * (m1 * x1 + m2 * x2))
    t0 = time.process_time()
    res["C2"] = lp_count(G, [x1, x2])
    print("C2 bubble:", res["C2"], "(expect 3)", f"{time.process_time()-t0:.1f}s", flush=True)
    U = x1 * x2 + x1 * x3 + x2 * x3
    G = sp.expand(U + x1 * x2 * x3 * P2 + U * (m1 * x1 + m2 * x2 + m3 * x3))
    t0 = time.process_time()
    res["C3"] = lp_count(G, [x1, x2, x3])
    print("C3 sunrise:", res["C3"], "(expect 7)", f"{time.process_time()-t0:.1f}s", flush=True)
    D = sp.Rational(13, 7)
    a, b = NU[0], NU[1]
    c4 = {
        "C4 Gamma": count_1d(-D + NU[0] / x, x, [x]),
        "C4 Bessel generic": count_1d(-D + a / (x - 1) + b / (x + 1), x, [x - 1, x + 1]),
        "C4 Bessel -1/2": count_1d(-D - sp.Rational(1, 2) / (x - 1) - sp.Rational(1, 2) / (x + 1), x,
                                   [x - 1, x + 1]),
        "C4 Kummer": count_1d(-D + (NU[0] - 1) / x + (NU[1] - NU[0] - 1) / (1 + x), x, [x, 1 + x]),
    }
    zz = sp.Rational(5, 17)
    c4["C5 Gauss 2F1"] = count_1d(NU[0] / x - NU[1] / (1 - x) - NU[2] * zz / (1 - zz * x), x,
                                 [x, 1 - x, 1 - zz * x])
    expect = {"C4 Gamma": 1, "C4 Bessel generic": 2, "C4 Bessel -1/2": 2, "C4 Kummer": 2, "C5 Gauss 2F1": 2}
    for k, v in c4.items():
        print(k, v, "(expect", expect[k], ")", flush=True)
        res[k] = v
    return res


def rand_pts(rng, n):
    pts = []
    while len(pts) < n:
        dh = Fr(rng.randint(1, 97), rng.randint(1, 97))
        sg = Fr(rng.randint(1, 97), rng.randint(1, 97))
        if sg in (0, 1):
            continue
        pts.append((sp.Rational(dh.numerator, dh.denominator), sp.Rational(sg.numerator, sg.denominator)))
    return pts


def slice_items():
    rng = random.Random(20261027)
    pts5 = rand_pts(rng, 5)
    pts3 = rand_pts(rng, 3)
    print("random points MC1:", pts5, "\nrandom points MC2:", pts3, flush=True)
    print("== MC1 (n1=n2=1)", flush=True)
    for dh, sg in pts5:
        print(" pt", dh, sg, "X generic", slice_X(dh, sg, 1, 1, True), "X actual", slice_X(dh, sg, 1, 1, False),
              "M generic", slice_M(sg, 1, 1, True), "M actual", slice_M(sg, 1, 1, False), flush=True)
    print("== MC2 (n1,n2 in {1,2,3})", flush=True)
    for dh, sg in pts3:
        for n1, n2 in itertools.product((1, 2, 3), repeat=2):
            print(" pt", dh, sg, "n", (n1, n2), "X gen", slice_X(dh, sg, n1, n2, True)["Q"],
                  "X act", slice_X(dh, sg, n1, n2, False), "M gen", slice_M(sg, n1, n2, True)["Q"],
                  "M act", slice_M(sg, n1, n2, False)["Q"], flush=True)
    print("== MC2b (n1 in 1..5, n2 in {1,3,5})", flush=True)
    for dh, sg in pts3:
        for n1, n2 in itertools.product((1, 2, 3, 4, 5), (1, 3, 5)):
            print(" pt", dh, sg, "n", (n1, n2), "X act", slice_X(dh, sg, n1, n2, False),
                  "X gen", slice_X(dh, sg, n1, n2, True)["Q"], flush=True)
    print("== MC4 (drop locus, n1=n2=1, Route X actual exponents and generic)", flush=True)
    dfix, sfix = pts5[0]
    sig_list = [s for _, s in rand_pts(rng, 20)]
    d_list = [d for d, _ in rand_pts(rng, 20)]
    vals = [(("sigma", s), slice_X(dfix, s, 1, 1, False)["Q"], slice_X(dfix, s, 1, 1, True)["Q"]) for s in sig_list]
    vals += [(("Dhat", d), slice_X(d, sfix, 1, 1, False)["Q"], slice_X(d, sfix, 1, 1, True)["Q"]) for d in d_list]
    for v in vals:
        print(" ", v, flush=True)
    print(" all random = 4:", all(v[1] == 4 and v[2] == 4 for v in vals), flush=True)
    print(" sigma=1:", slice_X(dfix, sp.Integer(1), 1, 1, False), slice_X(dfix, sp.Integer(1), 1, 1, True), flush=True)
    print(" sigma=0:", slice_X(dfix, sp.Integer(0), 1, 1, False), slice_X(dfix, sp.Integer(0), 1, 1, True), flush=True)
    print(" Dhat=0:", slice_X(sp.Integer(0), sfix, 1, 1, False), slice_X(sp.Integer(0), sfix, 1, 1, True), flush=True)


def fibre_items():
    rng = random.Random(20261028)
    print("== MC6 (two-offset fibre, conic quotient, Route F)", flush=True)
    for i in range(3):
        c1 = sp.Rational(rng.randint(2, 40), rng.randint(41, 97))
        r = sp.Rational(rng.randint(1, 30), rng.randint(31, 60))  # c2/c1 = r^2, rational square
        c2 = c1 * r**2
        D1 = sp.Rational(rng.randint(1, 50), rng.randint(1, 50))
        D2 = sp.Rational(rng.randint(1, 50), rng.randint(1, 50))
        X = (2 * c1 * tt - c1 * tt**2 - c2) / (c2 - c1 * tt**2)
        W = 1 + tt * (X - 1)
        assert sp.simplify(c2 * (X**2 - 1) - c1 * (W**2 - 1)) == 0
        br = sp.factor(sp.numer(sp.together(X**2 - 1)))
        branch = [tt, tt - 1, tt - c2 / c1]
        for D2v, lab in ((D2, "D2!=0"), (sp.Integer(0), "D2=0 (reported)")):
            f = -D1 * X - D2v * W
            om_g = sp.diff(f, tt) + sum(NU[j] / branch[j] for j in range(3))
            om_a = sp.diff(f, tt) - sum(sp.Rational(1, 2) / b for b in branch)
            rem = branch + [c2 - c1 * tt**2]
            print(f" pt{i} c1={c1} c2={c2} D1={D1} D2={D2v} [{lab}] numer(x^2-1)={br}", flush=True)
            print("   generic:", count_1d(om_g, tt, rem), " actual -1/2:", count_1d(om_a, tt, rem), flush=True)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    if what in ("controls", "all"):
        controls()
    if what in ("slice", "all"):
        slice_items()
    if what in ("fibre", "all"):
        fibre_items()
    print(f"wall {time.time()-t0:.1f}s", flush=True)
