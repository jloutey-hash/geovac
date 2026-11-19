"""exch_reverse_review_compare.py -- REVIEW: own closed-form evaluation (exch_reverse_review_closed.py) vs
exchange_hp (exch_reverse_review_hp.py) at the review's new a > b points.  -> debug/data/exch_reverse_review_compare.json"""
import glob
import json

import mpmath as mp

mp.mp.dps = 60
cl = json.load(open("debug/data/exch_reverse_review_closed.json"))
rows = []
for f in sorted(glob.glob("debug/data/exch_reverse_review_hp_*.json")):
    h = json.load(open(f))
    key = "new_%s_%s_%s_%s" % (h["case"], h["pair"][0], h["pair"][1], h["R"])
    if key not in cl or h["tau_max"] < 16:
        continue
    cf, hp = mp.mpf(cl[key]), mp.mpf(h["value"])
    d = abs(cf - hp)
    r = {"case": h["case"], "pair": h["pair"], "t": str(mp.mpf(eval(h["pair"][1] + "/1")) / mp.mpf(eval(h["pair"][0] + "/1"))),
         "R": h["R"], "tau_max": h["tau_max"], "closed": mp.nstr(cf, 50), "hp": h["value"][:52],
         "abs_diff": mp.nstr(d, 3), "digits": round(float(-mp.log10(d / abs(cf))), 1) if d else 99,
         "hp_last_term": h["per_tau"][-1], "hp_time_s": h["time_s"]}
    rows.append(r)
    print(r, flush=True)
json.dump(rows, open("debug/data/exch_reverse_review_compare.json", "w"), indent=1)
