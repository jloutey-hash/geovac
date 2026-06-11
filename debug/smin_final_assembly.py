"""Final S_min closed form: substitute the identified t-value reductions
(round-1 dossier + s3_pslq stage-A 220-dps identifications) into the round-1
decomposition, collect exactly in sympy, and verify numerically at 130 dps
against the direct Hurwitz-tail anchor.

Inputs: 8-t-value decomposition of debug/smin_dossier_round1_memo.md Q2;
t(2,1), t(2,3) from round 1; t(4,1), t(4,3) from s3_decomp_pslq stage A.
Output: debug/data/smin_final_assembly.json
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath
import sympy as sp

mpmath.mp.dps = 130

pi, ln2, z3, z5, z7 = sp.symbols("pi ln2 z3 z5 z7")

ZETA_EVEN = {2: pi**2 / 6, 4: pi**4 / 90, 6: pi**6 / 945, 8: pi**8 / 9450}
ZETA_ODD = {3: z3, 5: z5, 7: z7}


def lam(s: int):
    base = ZETA_EVEN[s] if s % 2 == 0 else ZETA_ODD[s]
    return (1 - sp.Rational(1, 2**s)) * base


# identified double t-values (t(b,c) = sum_{o2>o1>=1 odd} o2^-b o1^-c)
T = {
    (2, 1): sp.Rational(-7, 16) * z3 + sp.Rational(1, 8) * pi**2 * ln2,
    (4, 1): (sp.Rational(1, 96) * pi**4 * ln2
             - sp.Rational(1, 64) * pi**2 * z3 - sp.Rational(31, 64) * z5),
    (2, 3): sp.Rational(1, 16) * pi**2 * z3 - sp.Rational(31, 64) * z5,
    (4, 3): (sp.Rational(1, 128) * pi**4 * z3
             - sp.Rational(5, 128) * pi**2 * z5
             - sp.Rational(127, 256) * z7),
    (2, 2): (lam(2)**2 - lam(4)) / 2,
    (4, 4): (lam(4)**2 - lam(8)) / 2,
}
PAIR24 = lam(2) * lam(4) - lam(6)  # t(2,4) + t(4,2)  [stuffle]

EC = {1: 1, 2: -3, 3: -1, 4: 3}
EB = {2: 1, 4: -1}
D = {3: 1, 4: -3, 5: -2, 6: 6, 7: 1, 8: -3}


def boundary(b: int, c: int):
    return ((lam(b) - 1)
            + sp.Rational(1, 3**c) * (lam(b) - 1 - sp.Rational(1, 3**b)))


S = sp.Integer(0)
for c, e1 in EC.items():
    for b, e2 in EB.items():
        coef = 64 * e1 * e2
        if (b, c) in ((2, 4), (4, 2)):
            S += coef * (-boundary(b, c))
        else:
            S += coef * (T[(b, c)] - boundary(b, c))
S += 64 * 3 * PAIR24
S += 32 * sum(ds * (lam(s) - 1 - sp.Rational(1, 3**s)) for s, ds in D.items())
S = sp.expand(S)

print("S_min closed form:")
print(" ", S)

f = sp.lambdify((pi, ln2, z3, z5, z7), S, "mpmath")
val = f(mpmath.pi, mpmath.log(2), mpmath.zeta(3), mpmath.zeta(5),
        mpmath.zeta(7))

anchor = mpmath.nsum(
    lambda k: (2 * mpmath.zeta(2, int(k) + mpmath.mpf(3) / 2)
               - mpmath.mpf(1) / 2
               * mpmath.zeta(4, int(k) + mpmath.mpf(3) / 2)) ** 2,
    [1, mpmath.inf], method="levin")

res = abs(val - anchor)
print("closed   =", mpmath.nstr(val, 60))
print("anchor   =", mpmath.nstr(anchor, 60))
print("residual =", mpmath.nstr(res, 5))

out = {
    "closed_form": str(S),
    "closed_value_100dps": mpmath.nstr(val, 100),
    "anchor_100dps": mpmath.nstr(anchor, 100),
    "residual": mpmath.nstr(res, 5),
    "verified_1e100": bool(res < mpmath.mpf(10) ** -100),
}
Path(__file__).parent.joinpath("data", "smin_final_assembly.json").write_text(
    json.dumps(out, indent=2))
print("verified at 1e-100:", out["verified_1e100"])
