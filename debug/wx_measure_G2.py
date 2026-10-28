"""WX phase 0, continued (neighbour cases only, D + 0.05; no scored case):
  - resummation control: sum_l F_l P_l(cos 104.5 deg), l = 0..60, vs direct Boys, for g3n and E3-3n;
  - the single-Gaussian q_l clause: max |q_l - expansion| / (0.5 max(|c5|, 10) / l^5) over l = 40..60, and p_30,
    for g3n and g4n (the WG record's g1 sat at 0.89 of its bound, so this margin is measured before reuse).
Usage: python debug/wx_measure_G2.py   -> debug/data/wx_measure_G2.json"""
import json
from pathlib import Path

import mpmath as mp

import wx_gauss as G

mp.mp.dps = 50
out = {}
x = mp.cos(mp.radians(mp.mpf('104.5')))


def nbr(name):
    k, a, b, D = G.CASES[name]
    return (k, a, b, str(mp.mpf(D) + mp.mpf('0.05')))


for name in ('g3', 'E3-3'):
    spec = nbr(name)
    F = G.F_case(spec, 60)
    s = sum(F[l] * mp.legendre(l, x) for l in range(61))
    out[f"resum_{name}n"] = float(abs(s / G.boys_direct(spec, x) - 1))
for name in ('g3', 'g4'):
    spec = nbr(name)
    F = G.F_case(spec, 61)
    T1 = G.quartets(spec)[0]['T1']
    c = G.q_coeffs(G.single_T(spec))
    worst = 0
    for l in range(40, 61):
        q = F[l + 1] * (2 * l + 3) / (F[l] * T1)
        pred = 1 + sum(cj / mp.mpf(l) ** (j + 2) for j, cj in enumerate(c))
        worst = max(worst, float(abs(q - pred) / (0.5 * max(abs(c[3]), 10) / mp.mpf(l) ** 5)))
    p30 = float(mp.log(F[30] / F[34]) / mp.log(mp.mpf(34) / 30))
    out[f"qclause_{name}n"] = dict(T=float(G.single_T(spec)), c=[float(v) for v in c], worst_over_bound=worst, p30=p30,
                                   positive=all(v > 0 for v in F))
print(json.dumps(out, indent=1))
(Path(__file__).parent / "data" / "wx_measure_G2.json").write_text(json.dumps(out, indent=1))
