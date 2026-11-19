"""exch_reverse_review_neart1.py -- REVIEW: closed form near t = 1 from below at raised precision (pole cancellation),
vs the exchange_hp JSONs.  -> debug/data/exch_reverse_review_neart1.json"""
import json
import sys
import mpmath as mp
sys.path.insert(0, "debug")
import exch_reverse_review_closed as C
out = []
for case, a, b, R, f in (("2ppi", "1", "99/100", "3/2", "2ppi_1_99-100_3-2_T20"), ("3ddel", "1", "9/10", "7/5", "3ddel_1_9-10_7-5_T20")):
    P = C.coeff_polys(case)
    hp = mp.mpf(json.load(open("debug/data/exch_reverse_review_hp_%s.json" % f))["value"])
    for dps in (60, 90, 130):
        mp.mp.dps = dps
        q = lambda s: mp.mpf(mp.fraction(*[int(x) for x in (s.split("/") + ["1"])[:2]]))
        v = C.K(P, q(a), q(b), q(R))
        d = abs(v - hp)
        r = {"case": case, "t": b, "R": R, "dps": dps, "closed": mp.nstr(v, 50), "abs_diff_vs_hp": mp.nstr(d, 3),
             "digits": round(float(-mp.log10(d / abs(v))), 1)}
        out.append(r); print(r, flush=True)
json.dump(out, open("debug/data/exch_reverse_review_neart1.json", "w"), indent=1)
