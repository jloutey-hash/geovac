"""exch_reverse_compare.py -- closed form (real branch, a > b) vs exchange_hp, from the JSONs written by
exch_reverse_closed.py and exch_reverse_hp.py.  Run: python debug/exch_reverse_compare.py -> debug/data/exch_reverse_compare.json"""
import glob
import json

import mpmath as mp

mp.mp.dps = 60
rows = []
for f in sorted(glob.glob("debug/data/exch_reverse_hp_*.json")):
    h = json.load(open(f))
    if h["tau_max"] < 16:
        continue
    tg = "%s_%s_%s" % (h["case"], h["pair"][0], h["pair"][1])
    c = json.load(open("debug/data/exch_reverse_closed_%s.json" % tg.replace("/", "-")))
    v = c["values_dps60"].get(h["R"])
    if v is None:
        continue
    cf, hp = mp.mpf(v["fresh"]), mp.mpf(h["value"])
    d = abs(cf - hp)
    rows.append({"case": h["case"], "pair": h["pair"], "R": h["R"], "closed_form": mp.nstr(cf, 50),
                 "exchange_hp": h["value"], "abs_diff": mp.nstr(d, 3),
                 "rel_digits": float(mp.nstr(-mp.log10(d / abs(cf)), 4)) if d else 99.0,
                 "hp_last_terms": h["last_terms"], "hp_tau_max": h["tau_max"], "hp_time_s": h["time_s"]})
    print(rows[-1], flush=True)
json.dump(rows, open("debug/data/exch_reverse_compare.json", "w"), indent=1)
