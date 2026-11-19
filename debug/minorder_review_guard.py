"""Reviewer checks on the guard `no_lower_order` (tests/test_paper58_r_operators.py), on Sugiura's formula, on the
asymptotic series used in Lemma GI, and on Theorem A (plants).

G1  memo control K2: f = (ln 6 - ln 2 - ln 3) R e^{-R} + e^{-2R}; true order 1; what does no_lower_order say?
G2  further plants of the same kind that use only constants that DO occur in the panel:
      f = (ln 6 - ln 2 - ln 3) e^{-R} E1(R) + e^{-2R}          (true order 1)
      f = (sqrt(2)^2 - 2) R e^{-R} + e^{-2R}                    (sqrt(2)^2 kept unevaluated; true order 1)
      f = (sqrt(6) - sqrt(2) sqrt(3)) R e^{-R} + e^{-2R}        (true order 1)
G3  the test's upper bound for one-electron integrals is rank AT A POINT of r+1 columns: plant a function of
    true order 3 whose 3 columns are dependent at R = 37/59 only.
G4  constants seen by the guard in each guarded case (is any pair dependent in truth?).
S   1s equal-exponent exchange closed form against Sugiura's formula (independent of geovac), R = 1.4, 0.5, 3.
A   DLMF 6.12: Ei(x) e^{-x} x and E1(x) e^{x} x against sum k!/x^k, sum (-1)^k k!/x^k (optimal truncation).
Writes debug/data/minorder_review_guard.json
"""
from __future__ import annotations

import importlib.util
import json
import sys
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location("_ro_review", ROOT / "tests" / "test_paper58_r_operators.py")
RO = importlib.util.module_from_spec(spec)
spec.loader.exec_module(RO)
R = RO.R
out = {}


def say(k, v):
    out[k] = v
    print(k, "->", v, flush=True)


# ---- G1, G2
plants = {
    "G1_K2": (sp.log(6) - sp.log(2) - sp.log(3)) * R * sp.exp(-R) + sp.exp(-2 * R),
    "G2_E1": (sp.log(6) - sp.log(2) - sp.log(3)) * sp.exp(-R) * sp.expint(1, R) + sp.exp(-2 * R),
    "G2_sqrt6": (sp.sqrt(6) - sp.Mul(sp.sqrt(2), sp.sqrt(3), evaluate=False)) * R * sp.exp(-R) + sp.exp(-2 * R),
    "G2_gamma_cancel": (sp.EulerGamma - sp.EulerGamma) * R * sp.exp(-R) + sp.exp(-2 * R),
}
for k, f in plants.items():
    try:
        _atoms, reps = RO.decompose(f)
        say(k + "_constants_seen", [str(x) for x in reps])
        say(k + "_no_lower_order(2)", bool(RO.no_lower_order(f, 2)))
        say(k + "_no_lower_order(1)", bool(RO.no_lower_order(f, 1)))
        say(k + "_true_value_of_planted_coefficient", str(sp.N(sp.diff(f.subs(sp.exp(-2 * R), 0), R).subs(R, 0), 30)))
    except Exception as e:  # noqa: BLE001
        say(k + "_error", repr(e)[:200])

# ---- G3: pointwise upper bound.  f = e^{-R} + (R - 37/59) e^{-2R} ... build a function with three rates whose
# coordinate columns f, f', f'' are dependent exactly at R0 (true order 3; "order 2 exists" would be accepted?)
x0 = sp.Rational(37, 59)
a, b = sp.symbols("a b")
g = sp.exp(-R) + sp.exp(-2 * R) + (a + b * (R - x0)) * sp.exp(-3 * R)
# columns at R0 on atoms e^{-R}, e^{-2R}, e^{-3R}:  f=(1,1,a), f'=(-1,-2,-3a+b), f''=(1,4,9a-6b)
M = sp.Matrix([[1, 1, a], [-1, -2, -3 * a + b], [1, 4, 9 * a - 6 * b]])
sol = sp.solve(M.det(), b)
say("G3_det_at_R0", str(sp.factor(M.det())))
bb = sol[0].subs(a, 1)
g1 = g.subs({a: 1, b: bb})
say("G3_plant", str(g1))
r_at_point = RO.rank_at(RO.deriv_vectors(g1, 2), x0)
say("G3_rank_at_37/59_of_f_f1_f2", int(r_at_point))
say("G3_rank_at_11/7_of_f_f1_f2", int(RO.rank_at(RO.deriv_vectors(g1, 2), sp.Rational(11, 7))))
say("G3_test_idiom_accepts_order_2_exists", bool(r_at_point == 2))
say("G3_no_lower_order(3)_at_that_point", bool(RO.no_lower_order(g1, 3)))

# ---- G4: constants in guarded cases
cases = {"aabb11": RO.aabb(1, 1), "aabb12": RO.aabb(1, 2), "hyb33": RO.hyb_p(3, 3), "hyb34": RO.hyb_p(3, 4),
         "exch1": RO.exch(1), "exch2": RO.exch(2)}
for k, f in cases.items():
    _a, reps = RO.decompose(f)
    say("G4_constants_" + k, [str(x) for x in reps])

# ---- S: Sugiura
mp.mp.dps = 40
f = RO.exch(1)


def sugiura(Rv):
    Rv = mp.mpf(Rv)
    S = (1 + Rv + Rv ** 2 / 3) * mp.exp(-Rv)
    Sp = (1 - Rv + Rv ** 2 / 3) * mp.exp(Rv)
    return (-mp.exp(-2 * Rv) * (-mp.mpf(25) / 8 + 23 * Rv / 4 + 3 * Rv ** 2 + Rv ** 3 / 3)
            + (6 / Rv) * (S ** 2 * (mp.euler + mp.log(Rv)) + Sp ** 2 * mp.ei(-4 * Rv) - 2 * S * Sp * mp.ei(-2 * Rv))) / 5


for Rv in ("1.4", "0.5", "3"):
    cf = mp.mpf(str(sp.N(f.subs(R, sp.Rational(Rv)), 45)))
    say(f"S_exch1s_R={Rv}", {"geovac": mp.nstr(cf, 25), "sugiura": mp.nstr(sugiura(Rv), 25),
                             "rel": mp.nstr(abs(cf - sugiura(Rv)) / abs(cf), 3)})

# ---- A: asymptotic series
for x in (30, 60):
    x = mp.mpf(x)
    n = int(x)
    s_pos = sum(mp.factorial(k) / x ** k for k in range(n))
    s_neg = sum((-1) ** k * mp.factorial(k) / x ** k for k in range(n))
    say(f"A_Ei_x={int(x)}", mp.nstr(mp.ei(x) * mp.exp(-x) * x - s_pos, 5))
    say(f"A_E1_x={int(x)}", mp.nstr(mp.e1(x) * mp.exp(x) * x - s_neg, 5))
    say(f"A_smallest_term_x={int(x)}", mp.nstr(mp.factorial(n) / x ** n, 5))

(ROOT / "debug" / "data" / "minorder_review_guard.json").write_text(json.dumps(out, indent=1, default=str))
