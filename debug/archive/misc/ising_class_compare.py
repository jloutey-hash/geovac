"""Compare Route A (mpmath, own K0) and Route B (series + arb) moments; check closed forms."""
import json
import mpmath as mp
mp.mp.dps = 320
A = json.load(open("debug/data/ising_class_moments_mp.json"))["moments"]
B = json.load(open("debug/data/ising_class_moments_arb.json"))["moments"]
G13 = mp.gamma(mp.mpf(1) / 3)
L3 = (mp.zeta(2, mp.mpf(1)/3) - mp.zeta(2, mp.mpf(2)/3)) / 9   # L_{-3}(2)
closed = {
    "1,0": mp.pi / 2, "1,1": mp.mpf(1), "1,2": mp.pi / 2, "1,3": mp.mpf(4),
    "2,0": mp.pi**2 / 4, "2,1": mp.mpf(1) / 2,
    "3,0": 3 * G13**6 / (32 * mp.pi * mp.cbrt(4)),
    "3,1": mp.mpf(3) / 4 * L3,
    "3,2": G13**6 / (96 * mp.pi * mp.cbrt(4)) - 4 * mp.pi**5 * mp.cbrt(4) / (9 * G13**6),
    "3,3": L3 - mp.mpf(2) / 3,
    "4,1": 7 * mp.zeta(3) / 8, "4,3": 7 * mp.zeta(3) / 32 - mp.mpf(3) / 16,
}
out = {}
for key in A:
    a = mp.mpf(A[key]["value"]); b = mp.mpf(B[key]["value"])
    d = abs(a - b) / abs(a)
    rec = {"rel_diff_AB": mp.nstr(d, 3), "digits_AB": int(-mp.log10(d)) if d > 0 else 290}
    if key in closed:
        c = closed[key]
        rec["rel_diff_A_closed"] = mp.nstr(abs(a - c) / abs(c), 3)
        rec["rel_diff_B_closed"] = mp.nstr(abs(b - c) / abs(c), 3)
    out[key] = rec
    print(key, rec)
# BBBG c_{4,0} = (pi^4/4) 4F3(1/2^4; 1^3; 1), checked at lower precision (slowly convergent series)
with mp.workdps(60):
    h = mp.hyper([0.5] * 4, [1] * 3, 1)
    c40 = mp.mpf(A["4,0"]["value"])
    out["c40_vs_4F3_rel_60dps"] = mp.nstr(abs(mp.pi**4 / 4 * h / c40 - 1), 3)
    print("c40 vs 4F3:", out["c40_vs_4F3_rel_60dps"])
json.dump(out, open("debug/data/ising_class_compare.json", "w"), indent=1)
