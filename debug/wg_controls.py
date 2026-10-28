"""WG phase 2 controls WG-0 (a)-(e) (record docs/predictions/2026-10-28_water_gaussian_contact.md).
Usage: python wg_controls.py a | b <case> | c | d | e      -> debug/data/wg_ctrl_<part>.json"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

import wg_common as C

OUT = Path(__file__).parent / "data"
part = sys.argv[1]
res = {}

if part == 'a':
    sys.path.insert(0, str(Path(__file__).parents[1] / "tests"))
    from test_paper59_water_bond_angle import F_coeffs
    stored = json.load(open(OUT / "wa5_cusp.json"))
    for s, (zX, zH) in {'E1': (1.0, 1.0), 'E2': (2.25, 1.24)}.items():
        S = stored[f"{s} (1.809, 1.809)"]["F"]
        F = F_coeffs(1.809, 1.809, zX, zH, 60, level=2)
        dev = {l: abs(F[l] / S[l] - 1) for l in (10, 20, 30, 40, 50, 60)}
        res[s] = dict(dev=dev, max=max(dev.values()), pass_=max(dev.values()) < 1e-10)
elif part == 'b':
    case = sys.argv[2]
    mp.mp.dps = 50   # first run used 30 digits; the Weber quadrature fell below its floor at g1 l = 30 (see wg_diag_b.py)
    ls = {'g1': [0, 5, 10, 20, 30], 'g2': [0, 5, 10, 20, 30], 'E1-3': [0, 10, 20], 'E1-6': [0, 20]}[case]
    qs = C.quartets(case)
    ser = {l: mp.mpf(0) for l in ls}
    web = {l: mp.mpf(0) for l in ls}
    for Q in qs:
        Fq = C.F_quartet(Q, max(ls))
        for l in ls:
            ser[l] += Fq[l]
            web[l] += C.F_weber_quartet(Q, l)
    dev = {l: float(abs(web[l] / ser[l] - 1)) for l in ls}
    res[case] = dict(dev=dev, ser={l: mp.nstr(ser[l], 20) for l in ls}, web={l: mp.nstr(web[l], 20) for l in ls},
                     max=max(dev.values()), pass_=max(dev.values()) < 1e-10)
    part = f"b_{case}"
elif part == 'c':
    mp.mp.dps = 50
    x = mp.cos(mp.radians(mp.mpf('104.5')))
    for case in ('g1', 'E1-3'):
        F = C.F_case(case, 60)
        s = sum(F[l] * mp.legendre(l, x) for l in range(61))
        d = C.boys_direct(case, x)
        dev = float(abs(s / d - 1))
        res[case] = dict(sum=mp.nstr(s, 25), direct=mp.nstr(d, 25), dev=dev, pass_=dev < 1e-14)
elif part == 'd':
    mp.mp.dps = 50
    T = mp.mpf('1.2')
    F = C.F_quartet(dict(pref=mp.mpf(1), T0=T, T1=T), 41)
    l = 40
    q = F[l + 1] * (2 * l + 3) / (F[l] * T)
    c = [1, 0, -T * (T + 4) / 4, T * (T + 4), T * (2 * T ** 3 + 12 * T ** 2 - 17 * T - 196) / 16,
         -T * (9 * T ** 3 + 48 * T ** 2 + 76 * T - 272) / 8]
    r = float(abs(q - sum(cj / mp.mpf(l) ** j for j, cj in enumerate(c))))
    res = dict(resid=r, pass_=r < 1e-7)
elif part == 'e':
    stored = json.load(open(OUT / "wa5_cusp.json"))
    for s in ('E1', 'E2'):
        S = stored[f"{s} (1.809, 1.809)"]["F"]
        p = {l: float(np.log(S[l] / S[l + 4]) / np.log((l + 4) / l)) for l in (21, 23, 42, 50)}
        res[s] = dict(p=p, pass_=all(7 < v < 9 for v in p.values()))
print(json.dumps(res, indent=1, default=str))
(OUT / f"wg_ctrl_{part}.json").write_text(json.dumps(res, indent=1, default=str))
