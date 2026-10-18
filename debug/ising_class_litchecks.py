"""Numerical checks of the literature statements on c_{n,0} against our moments and L-values.
LTY = Lim-Tu-Yu arXiv:2609.13576 eq. (D_k:example), (D_k:example?); BBBG = arXiv:0801.0891 Table 1."""
import json
import mpmath as mp
mp.mp.dps = 270
MOM = json.load(open("debug/data/ising_class_moments_arb.json"))["moments"]
LV = json.load(open("debug/data/ising_class_lvalues.json"))["values"]
c = lambda n, k: mp.mpf(MOM[f"{n},{k}"]["value"])
L = lambda key: mp.mpf(LV[key])
f48 = "48.3.e.a = 12.3.c.a x chi_-4"; f8 = "8.4.a.a eta(2t)^4eta(4t)^4"; f24 = "24.6.a.a"
checks = {
    "LTY c30 = (pi^2/2) L(48.3.e.a,1)  [proven]": c(3, 0) / (mp.pi**2 / 2 * L(f"{f48}|1")) - 1,
    "LTY c40 = 4 pi^2 L(8.4.a.a,2)  [proven]": c(4, 0) / (4 * mp.pi**2 * L(f"{f8}|2")) - 1,
    "LTY c60 = 108 pi^2 L(24.6.a.a,4)  [numerical, '?']": c(6, 0) / (108 * mp.pi**2 * L(f"{f24}|4")) - 1,
    "RWZ L(8.4.a.a,3) relation check: L(f2,3) = pi^2/4 L(f2,1)": L(f"{f8}|3") / (mp.pi**2 / 4 * L(f"{f8}|1")) - 1,
}
out = {}
for k, v in checks.items():
    out[k] = mp.nstr(abs(v), 3); print(k, out[k])
json.dump(out, open("debug/data/ising_class_litchecks.json", "w"), indent=1)
