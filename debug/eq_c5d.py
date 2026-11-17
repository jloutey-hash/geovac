"""c5d: planted (1+R)...(5+R) f_2ppi(1,2) on one prime: minimal operator under DMAX 60 vs 100, coefficient equality."""
import json, os, sys
import numpy as np
sys.path.insert(0, "debug")
import eq_controls as C, eq_gen as G, exch_gen as XG
p = 67108859
ser, _ = G.generate("2ppi", 1, 2, 450, XG.RingP(p))
q = [1]
for root in (1, 2, 3, 4, 5):
    q = list(np.convolve(q, [root, 1]))
qq = np.zeros(451, dtype=np.int64); qq[:len(q)] = [int(x) % p for x in q]
ser = {k: np.mod(np.convolve(v, qq)[:451], p) for k, v in ser.items()}
r60, r100 = C._readings(ser, p, 60), C._readings(ser, p, 100)
m60, m100 = r60["minimal_operator"], r100["minimal_operator"]
out = {"first_hit_60": r60["first_hit"], "first_hit_100": r100["first_hit"],
       "shape_60": [m60["order"], m60["degree"]], "shape_100": [m100["order"], m100["degree"]],
       "coeffs_equal": m60["coeffs_mod_p"] == m100["coeffs_mod_p"]}
json.dump(out, open("debug/data/eq_controls_c5d.json", "w"), indent=1, default=str)
print(out)
