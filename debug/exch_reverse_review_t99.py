"""exch_reverse_review_t99.py -- REVIEW: the t = 99/100 closed-form/hp gap.  Closed form by the memo's own route
(exact Fraction coefficients at t = 99/100, exch_reverse_core.value) at dps 100, and by a FRESH exact solve at
(a, b) = (1, 99/100); compare with the review's own evaluator and with exchange_hp."""
import json
import sys
from fractions import Fraction
import mpmath as mp
sys.path.insert(0, "debug")
import exch_reverse_core as RC
mp.mp.dps = 100
a, b, R = Fraction(1), Fraction(99, 100), Fraction(3, 2)
co = RC.store_coords_Q("2ppi", a, b)
Kst = RC.coords_to_K(a, b, co)
v_st = RC.value(a, b, Kst, R)
print("stored-coords (memo route)", mp.nstr(v_st, 60), flush=True)
Kf, cof, info = RC.fresh_K_Q("2ppi", a, b)
v_fr = RC.value(a, b, Kf, R)
print("fresh solve", mp.nstr(v_fr, 60), info, "fresh==stored", cof == co, flush=True)
mine = mp.mpf(json.load(open("debug/data/exch_reverse_review_diagcf.json"))["2ppi_1_99/100_3/2"])
hp = mp.mpf(json.load(open("debug/data/exch_reverse_review_hp_2ppi_1_99-100_3-2_T14.json"))["value"])
out = {"stored": mp.nstr(v_st, 60), "fresh": mp.nstr(v_fr, 60), "fresh_equals_stored": cof == co,
       "stored_minus_mine": mp.nstr(v_st - mine, 3), "stored_minus_hp": mp.nstr(v_st - hp, 3)}
print(out)
json.dump(out, open("debug/data/exch_reverse_review_t99.json", "w"), indent=1)
