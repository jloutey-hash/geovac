"""One-formula check: the c = 0 case of the cosine transform reproduces Paper 59's three-centre fibre.
 (1) int_0^oo K0(a x) K0(b x) dx = pi^2/(4 sqrt(ab)) P_(-1/2)((a^2+b^2)/(2ab))     [c = 0 of the cosine form]
 (2) Paper 59 eq:period, int_0^oo dk / sqrt((c1 k^2+1)(c2 k^2+1)) = K(m)/sqrt(cmax), m = 1 - cmin/cmax,
     against the same Legendre value through Parseval: with a_i = 1/sqrt(c_i), the left side of (2) equals
     (2/(pi sqrt(c1 c2))) int K0(a1 x) K0(a2 x) dx (eq:besselmoment at p = 0).
Usage: python debug/cf_one_formula_check.py"""
import json
from pathlib import Path

import mpmath as mp

mp.mp.dps = 30
out = []
for c1, c2 in [(mp.mpf("0.21"), mp.mpf("0.09")), (mp.mpf("0.16"), mp.mpf("0.25")), (mp.mpf("0.05"), mp.mpf("0.2"))]:
    a, b = 1 / mp.sqrt(c1), 1 / mp.sqrt(c2)
    moment = mp.quad(lambda y: mp.besselk(0, a * y) * mp.besselk(0, b * y), [0, 0.5, 1, 2, 4, 8, mp.inf])
    leg = mp.pi**2 / (4 * mp.sqrt(a * b)) * mp.legenp(-0.5, 0, (a * a + b * b) / (2 * a * b), type=3)
    cmax, cmin = max(c1, c2), min(c1, c2)
    period = mp.ellipk(1 - cmin / cmax) / mp.sqrt(cmax)
    via = 2 / (mp.pi * mp.sqrt(c1 * c2)) * leg
    r = {"c": [float(c1), float(c2)], "moment_vs_legendre": float(abs(moment - leg) / leg),
         "period_vs_legendre": float(abs(period - via) / period)}
    out.append(r)
    print(r)
json.dump(out, open(Path(__file__).parent / "data" / "cf_one_formula_check.json", "w"), indent=1)
