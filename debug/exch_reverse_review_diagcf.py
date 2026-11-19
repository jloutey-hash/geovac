"""exch_reverse_review_diagcf.py -- REVIEW: closed form (own evaluator) at the near-t=1 diagnostic points, dps 100.
-> debug/data/exch_reverse_review_diagcf.json"""
import json
import sys
from fractions import Fraction
import mpmath as mp
sys.path.insert(0, "debug")
import exch_reverse_review_closed as C
mp.mp.dps = 100
q = lambda s: mp.mpf(Fraction(s).numerator) / Fraction(s).denominator
out = {}
for case, a, b, R in (("2ppi", "100/101", "1", "3/2"), ("2ppi", "1", "97/100", "3/2"), ("2ppi", "1", "99/100", "1"),
                      ("2psig", "1", "99/100", "3/2"), ("2ppi", "1", "99/100", "3/2")):
    v = C.K(C.coeff_polys(case), q(a), q(b), q(R))
    out["%s_%s_%s_%s" % (case, a, b, R)] = mp.nstr(v, 70)
    print(case, a, b, R, mp.nstr(v, 55), flush=True)
json.dump(out, open("debug/data/exch_reverse_review_diagcf.json", "w"), indent=1)
