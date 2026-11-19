"""Controls for minorder_core.analyse: it must be able to say "not proved".

K1  inflated order: the 1s equal-exponent exchange integral with claimed order 5, and (AA|BB) 1s (1,2) with
    claimed order 4  ->  lower bound NOT proved.
K2  constants that are dependent in truth but independent as symbols:
    f = (ln 6 - ln 2 - ln 3) R e^{-R} + e^{-2R}.  Symbol-split rank says order 2; the true order is 1.
    The true-constant minor must refuse order 2 (and accept order 1).
K3  an E1 pair that collapses: f = e^{-R} E1(2R) + e^{-R} Ei(-2R) + e^{-3R} is e^{-3R}: order 1, not 2.
Writes debug/data/minorder_controls.json.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "debug"))

import minorder_core as MC  # noqa: E402
from geovac import qfd_core as Q  # noqa: E402
from geovac.two_center_eri import aabb_closed_form  # noqa: E402

R = MC.R
S1 = (1, 0, 0)
F = Fraction
out = {}

ex = Q.exchange_closed_form(F(1), S1, S1, F(1), S1, S1, R, tau_max=6)
out["K1_exchange_r5"] = MC.analyse(ex, 5, upper=False)["proved_lower"]
out["K1_exchange_r4"] = MC.analyse(ex, 4, upper=False)["proved_lower"]
j = aabb_closed_form(F(1), S1, S1, F(2), S1, S1)
out["K1_aabb_r4"] = MC.analyse(j, 4, upper=False)["proved_lower"]
out["K1_aabb_r3"] = MC.analyse(j, 3, upper=False)["proved_lower"]

f2 = (sp.log(6) - sp.log(2) - sp.log(3)) * R * sp.exp(-R) + sp.exp(-2 * R)
for tag, f in (("log6_log2_log3", f2),):
    try:
        atoms, reps = MC.RO.decompose(f)
        out[f"K2_{tag}_constants"] = [str(k) for k in reps]
        out[f"K2_{tag}_split_says_order_2"] = MC.RO.rank_at(MC.RO.deriv_vectors(f, 1, "split"), sp.Rational(37, 59)) == 2
        out[f"K2_{tag}_true_r2"] = MC.analyse(f, 2, upper=False)["proved_lower"]
        out[f"K2_{tag}_true_r1"] = MC.analyse(f, 1, upper=False)["proved_lower"]
    except Exception as e:   # noqa: BLE001
        out[f"K2_{tag}_error"] = repr(e)

f3 = sp.exp(-R) * sp.expint(1, 2 * R) + sp.exp(-R) * sp.Ei(-2 * R) + sp.exp(-3 * R)
try:
    rec = MC.analyse(f3, 2, upper=False)
    out["K3_r2"] = rec["proved_lower"]
    out["K3_atoms"] = rec["f_atoms"]
except AssertionError as e:
    out["K3_r2"] = "REFUSED: " + str(e)[:80]

for k, v in out.items():
    print(k, v, flush=True)
(ROOT / "debug" / "data" / "minorder_controls.json").write_text(json.dumps(out, indent=1))
