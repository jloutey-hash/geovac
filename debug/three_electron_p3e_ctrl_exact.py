import json
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 150
v = {k: mp.mpf(x) for k, x in json.loads(open("debug/data/three_electron_p3e_values_150.json").read()).items()}
# PSLQ reductions at masses (1,2,4), 60 dps: m*k_m + sum k_i B_i = 0
red = {(3, 1, 0, 0): (-1155, {'V': -76, 'D1': 273, 'D2': 348, 'S12': -876}),
       (4, 1, 1, 0): (24255, {'V': -4372, 'D1': 12264, 'D2': 20700, 'S12': -40533}),
       (4, 1, 0, 1): (-97020, {'V': -20179, 'D1': 58695, 'D2': 94356, 'S12': -187350}),
       (5, 1, 1, 1): (-132300, {'V': 95972, 'D1': -260022, 'D2': -452799, 'S12': 824619})}
b, c = 2, 4
weights = {(3, 1, 0, 0): 9 * b * c, (4, 1, 0, 1): 3 * b * c * c, (4, 1, 1, 0): 3 * b * b * c, (5, 1, 1, 1): b * b * c * c}
coef = {n: Fr(0) for n in ("V", "D1", "D2", "S12")}
for key, w in weights.items():
    km, d = red[key]
    for n, k in d.items():
        coef[n] += Fr(-k, km) * w
print("Gctrl =", {n: str(x) for n, x in coef.items()})
pred = sum(mp.mpf(x.numerator) / x.denominator * v[n] for n, x in coef.items())
print("residual at 150 dps:", mp.nstr(pred - v["Gctrl"], 5))
den = 1
for x in coef.values():
    den = den * x.denominator // __import__("math").gcd(den, x.denominator)
print("common denominator:", den, " integer coeffs:", {n: int(x * den) for n, x in coef.items()})
json.dump({n: str(x) for n, x in coef.items()}, open("debug/data/three_electron_p3e_control_reduction.json", "w"), indent=1)
