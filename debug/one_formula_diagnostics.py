"""Diagnostics on the one function L(a, b; c) (not registered; exploratory, after the GA registration 5e40da9a).
 (1) Role permutations at gamma = 3/2: which of the three inverse lengths sits in the momentum slot.
     j(x) for (masses {t,1}, momentum gamma) = J1; ({t,gamma}, 1); ({1,gamma}, t).  Compared with J1 and the
     EC record's isogenous J2 by traces of Frobenius: a_p^2 at several good primes for curves of these j at
     rational x0 (isogenous up to twist <=> all a_p^2 agree).
 (2) gamma = 1: j_C(x) = 16 (x+12)^3/(x-4)^2 (the Legendre family at gamma = 1) against the elliptic curve under
     F_3's Clausen square, j_K(x) = 16 (16-x)^3/x^2 (derived: 3F2 = 2F1(1/4,1/4;1;x/4)^2 and
     2F1(1/4,1/4;1;4y(1-y)) = 2F1(1/2,1/2;1;y), so y(1-y) = x/16).  Same a_p^2 test.
 (3) The dispersion-root form L(a,b;c) = (pi/4) int_{-oo}^{oo} dk / (Delta_a(k) Delta_b(k-c)), at real c and at
     c = i gamma, with each root as the product of principal square roots of its two linear factors.
Usage: python debug/one_formula_diagnostics.py"""
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

x = sp.Symbol("x")
DATA = Path(__file__).parent / "data"
OUT = {}


def j_from_z2(z2):
    return sp.cancel(64 * (z2 + 3) ** 3 / (1 - z2) ** 2)


def ap_squared(j, p):
    """a_p^2 for a curve with invariant j mod p (any twist); None if bad reduction or j in {0, 1728}."""
    jp = sp.Rational(j).p * pow(sp.Rational(j).q, -1, p) % p if sp.Rational(j).q % p else None
    if jp is None or jp % p == 0 or (jp - 1728) % p == 0:
        return None
    A = 3 * jp * (1728 - jp) % p
    B = 2 * jp * (1728 - jp) ** 2 % p
    if (4 * A**3 + 27 * B**2) % p == 0:
        return None
    sq = {}
    for y in range(p):
        sq[y * y % p] = sq.get(y * y % p, 0) + 1
    n = sum(sq.get((X**3 + A * X + B) % p, 0) for X in range(p)) + 1
    return (p + 1 - n) ** 2


def compare(j1, j2, x0s, primes):
    agree = total = 0
    for x0 in x0s:
        for p in primes:
            a1, a2 = ap_squared(j1.subs(x, x0), p), ap_squared(j2.subs(x, x0), p)
            if a1 is None or a2 is None:
                continue
            total += 1
            agree += a1 == a2
    return agree, total


if __name__ == "__main__":
    g = sp.Rational(3, 2)
    P = 16 * x**2 + 152 * x + 25
    J1 = P**3 / (x * (4 * x - 1) ** 2 * (4 * x - 25) ** 2)
    d2 = json.load(open(DATA / "c2b_ident_candidate_J2.json"))          # EC stage C form: J = P^3 / (lam D)
    w = json.loads(d2["widths"]) if isinstance(d2["widths"], str) else d2["widths"]
    J2 = sp.sympify(d2["P"]) ** 3 / (sp.Rational(d2["lam"]) * x ** w[0] * (4 * x - 1) ** w[1] * (4 * x - 25) ** w[2])   # stage C's D_of form
    assert sp.factor(sp.numer(sp.together(J2 - 1728))) == sp.factor((4 * x - 5) ** 2 * (16 * x**2 - 112 * x + 25) ** 2)
    assert sp.cancel(J2 - J1) != 0
    fam = {
        "masses {t,1}, momentum gamma": j_from_z2((x + 1 - g**2) ** 2 / (4 * x)),
        "masses {t,gamma}, momentum 1": j_from_z2((x + g**2 - 1) ** 2 / (4 * x * g**2)),
        "masses {1,gamma}, momentum t": j_from_z2((1 + g**2 - x) ** 2 / (4 * g**2)),
    }
    assert sp.cancel(fam["masses {t,1}, momentum gamma"] - J1) == 0
    primes = [p for p in sp.primerange(5, 80)]
    x0s = [sp.Rational(k) for k in (2, 3, 7, 11)]
    OUT["role_permutations"] = {}
    for name, j in fam.items():
        row = {"j": str(sp.factor(j)), "vs_J1": compare(j, J1, x0s, primes)}
        if J2 is not None:
            row["vs_J2"] = compare(j, J2, x0s, primes)
        OUT["role_permutations"][name] = row
        print(name, "| j =", sp.factor(j), "| a_p^2 agree with J1:", row["vs_J1"], "| with J2:", row.get("vs_J2"))
    if J2 is not None:
        OUT["J1_vs_J2_control"] = compare(J1, J2, x0s, primes)
        print("control J1 vs J2 (isogenous, must agree):", OUT["J1_vs_J2_control"])
    jC = 16 * (x + 12) ** 3 / (x - 4) ** 2
    assert sp.cancel(j_from_z2(x / 4) - jC) == 0
    jK = 16 * (16 - x) ** 3 / x**2
    y = sp.Symbol("y")
    assert sp.cancel((256 * (1 - y + y**2) ** 3 / (y**2 * (1 - y) ** 2)).subs(y * (1 - y), x / 16).subs(y, sp.Rational(1, 2) - sp.sqrt(1 - x / 4) / 2) - jK) == 0
    OUT["gamma1_vs_clausen"] = compare(jC, jK, [sp.Rational(k) for k in (2, 3, 5, 7, 11)], primes)
    print("gamma = 1 Legendre family vs F_3's Clausen curve, a_p^2 agree:", OUT["gamma1_vs_clausen"])
    OUT["gamma1_self_control"] = compare(jC, jC, [sp.Rational(2)], primes)
    l = sp.Symbol("l")
    OUT["legendre_2isogeny_control"] = compare((256 * (l**2 - l + 1) ** 3 / (l**2 * (l - 1) ** 2)).subs(l, x), (16 * (l**2 + 14 * l + 1) ** 3 / (l * (1 - l) ** 4)).subs(l, x), [sp.Rational(k) for k in (2, 3, 5, 7)], primes)
    print("Legendre 2-isogeny control (must agree):", OUT["legendre_2isogeny_control"])
    # (3) dispersion-root form
    mp.mp.dps = 25
    def Lcos(a, b, c):
        return mp.quad(lambda yv: mp.besselk(0, a * yv) * mp.besselk(0, b * yv) * mp.cos(c * yv), [0, 0.5, 1, 2, 4, 8, 16, mp.inf])
    def Lcosh(a, b, gm):
        return mp.quad(lambda yv: mp.besselk(0, a * yv) * mp.besselk(0, b * yv) * mp.cosh(gm * yv), [0, 0.5, 1, 2, 4, 8, 16, 32, mp.inf])
    def kform(a, b, c):
        f = lambda k: 1 / (mp.sqrt(k + 1j * a) * mp.sqrt(k - 1j * a) * mp.sqrt(k - c + 1j * b) * mp.sqrt(k - c - 1j * b))
        return mp.pi / 4 * mp.quad(f, [-mp.inf, -8, -2, -0.5, 0, 0.5, 2, 8, mp.inf])
    a, b, c = mp.mpf("0.7"), mp.mpf("1.3"), mp.mpf(2)
    r1 = abs(Lcos(a, b, c) - kform(a, b, c)) / abs(Lcos(a, b, c))
    a, b, gm = mp.mpf("1.7"), mp.mpf(2), mp.mpf(3) / 2       # gamma < b: no branch point crosses the contour
    r2 = abs(Lcosh(a, b, gm) - mp.re(kform(a, b, 1j * gm))) / abs(Lcosh(a, b, gm))
    a, b, gm = mp.mpf("1.7"), mp.mpf(1), mp.mpf(3) / 2       # gamma > b: the branch point i(gamma - b) has crossed
    v3, k3 = Lcosh(a, b, gm), kform(a, b, 1j * gm)
    OUT["kform"] = {"real_c": float(r1), "imag_c_gamma_below_b": float(r2),
                    "imag_c_gamma_above_b": {"cosh": str(v3), "kform": str(k3), "ratio": str(k3 / v3)}}
    print("k-form vs cosine form, real c:", float(r1))
    print("k-form vs cosh form, gamma < b:", float(r2))
    print("gamma > b: cosh =", v3, " k-form =", k3)
    json.dump(OUT, open(DATA / "one_formula_diagnostics.json", "w"), indent=1, default=str)
