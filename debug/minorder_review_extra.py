"""Reviewer extras: (X1) redundant-atom plant through the committed guard; (X2) class-C scaling K(R; z) = z K(zR; 1)
for 1s and for 2p-sigma 1s; (X3) (AA|BB) 1s: three rates with nonzero coefficient at further unequal pairs
(Theorem A input); (X4) Theorem A plants: number of rates == order by an exact Wronskian-type rank.
Writes debug/data/minorder_review_extra.json
"""
from __future__ import annotations

import importlib.util
import json
import sys
from fractions import Fraction as F
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location("_ro_review2", ROOT / "tests" / "test_paper58_r_operators.py")
RO = importlib.util.module_from_spec(spec)
spec.loader.exec_module(RO)
from geovac import qfd_core as Q  # noqa: E402
from geovac.two_center_eri import aabb_closed_form  # noqa: E402

R = RO.R
S1 = (1, 0, 0)
out = {}


def say(k, v):
    out[k] = v
    print(k, "->", v, flush=True)


# X1
f3 = sp.exp(-R) * sp.expint(1, 2 * R) + sp.exp(-R) * sp.Ei(-2 * R) + sp.exp(-3 * R)
say("X1_value_of_f3_minus_exp(-3R)_at_R=0.7", str(sp.N((f3 - sp.exp(-3 * R)).subs(R, sp.Rational(7, 10)), 25)))
try:
    say("X1_atoms_seen", [str(a) for a in RO.decompose(f3)[0]])
    say("X1_no_lower_order(f3, 2)", bool(RO.no_lower_order(f3, 2)))
    say("X1_no_lower_order(f3, 3)", bool(RO.no_lower_order(f3, 3)))
except Exception as e:  # noqa: BLE001
    say("X1_error", repr(e)[:200])

# X2
for z in (F(3, 2), 2):
    fz = RO.exch(z)
    f1 = RO.exch(1)
    for Rv in (sp.Rational(7, 10), sp.Rational(9, 5)):
        a = sp.N(fz.subs(R, Rv), 30)
        b = sp.N(sp.Rational(z) * f1.subs(R, sp.Rational(z) * Rv), 30)
        say(f"X2_1s_z={z}_R={Rv}", str(sp.N(abs(a - b) / abs(a), 5)))
P0 = (2, 1, 0)


def ex(ZA, oa, ob, ZB, oc, od):
    t = Q.exchange_closed_form(F(ZA), oa, ob, F(ZB), oc, od, R, tau_max=6, return_terms=True)
    assert all(x == 0 for x in t[-2:])
    return sp.Add(*t)


g1 = ex(2, P0, S1, 1, P0, S1)
g2 = ex(4, P0, S1, 2, P0, S1)
for Rv in (sp.Rational(7, 10), sp.Rational(9, 5)):
    a = sp.N(g2.subs(R, Rv), 30)
    b = sp.N(2 * g1.subs(R, 2 * Rv), 30)
    say(f"X2_2psig_z=2_R={Rv}", str(sp.N(abs(a - b) / abs(a), 5)))
say("X2_2psig_z=2_atoms", [str(a) for a in RO.decompose(g2)[0]])
say("X2_2psig_z=2_no_lower_order_4", bool(RO.no_lower_order(g2, 4)))

# X3
for za, zb in ((1, 3), (2, 5), (F(7, 3), 1), (F(1, 2), F(1, 3)), (1, F(1001, 1000))):
    f = aabb_closed_form(F(za), S1, S1, F(zb), S1, S1)
    atoms, reps = RO.decompose(f)
    say(f"X3_aabb_{za}_{zb}", {"atoms": [str(a) for a in atoms], "constants": [str(k) for k in reps]})

# X4 Theorem A plants: exact rank of the coordinate matrix on exponentials for random rational-coefficient sums
import random
random.seed(5)
x = R
ok = True
for trial in range(40):
    k = random.randint(1, 5)
    rates = random.sample([sp.Rational(p, q) for p in range(-6, 7) for q in (1, 2, 3)], k)
    rates = list(dict.fromkeys(rates))
    coeffs = [sp.Rational(random.randint(1, 9), random.randint(1, 5)) * x ** random.randint(-2, 3)
              + random.randint(0, 3) for _ in rates]
    coeffs = [c if c != 0 else sp.Integer(1) for c in coeffs]
    cols = [list(coeffs)]
    for _ in range(len(rates)):
        cols.append([sp.diff(c, x) + lam * c for c, lam in zip(cols[-1], rates)])
    M = sp.Matrix([[sp.together(c) for c in col] for col in cols]).T
    rk_full = M.rank()
    rk_head = M[:, :len(rates)].rank()
    ok &= (rk_full == len(rates) == rk_head)
say("X4_theoremA_40_random_plants_order_equals_number_of_rates", bool(ok))

(ROOT / "debug" / "data" / "minorder_review_extra.json").write_text(json.dumps(out, indent=1, default=str))
