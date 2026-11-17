"""Record PP scorer (docs/predictions/2026-10-29_hydrogen_p_functions.md, with Addendum 1). Controls first.
Reads debug/data/pp/*.json written by debug/pp_launch.py; writes data/pp_score.json.
    python debug/pp_score.py
"""
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np

mp.mp.dps = 600
HERE = Path(__file__).parent
DATA = HERE / 'data' / 'pp'
Z = mp.mpf(6) / 5
GEOM = {'U': (mp.mpf(9) / 5, mp.mpf(12) / 5, mp.mpf(3) / 4), 'E': (mp.mpf(19) / 10, mp.mpf(19) / 10, mp.mpf(1))}
PRED = {'PA-U': -4, 'PB-U': -4, 'PC-U': -5, 'PA-E': -10, 'PB-E': -10}
_cache = {}


def q(g, d1, d2, a, b, st=30, lm=150):
    name = f'{g}_{d1}_{d2}_{a}_{b}_{st}' + ('' if lm == 150 else f'_L{lm}')
    if name not in _cache:
        _cache[name] = [mp.mpf(v) for v in json.load(open(DATA / f'{name}.json'))['q']]
    return _cache[name]


def lin(*terms):
    """sum of c * list, elementwise."""
    n = min(len(v) for _, v in terms)
    return [mp.fsum(c * v[i] for c, v in terms) for i in range(n)]


def O1(g, d1, d2, b, st=30, lm=150):
    h = mp.mpf(10) ** -st
    return lin((1 / Z ** 2, q(g, d1, d2, '0', b, st, lm)),
               (-1 / (2 * h * Z), q(g, d1, d2, 'p', b, st, lm)), (1 / (2 * h * Z), q(g, d1, d2, 'm', b, st, lm)))


def O12(g, d1, d2, st=30, lm=150):
    h = mp.mpf(10) ** -st
    return lin((1 / Z ** 2, O1(g, d1, d2, '0', st, lm)),
               (-1 / (2 * h * Z), O1(g, d1, d2, 'p', st, lm)), (1 / (2 * h * Z), O1(g, d1, d2, 'm', st, lm)))


def PA(g, st=30):
    h = mp.mpf(10) ** -st
    return lin((1 / (2 * h), O1(g, 'p', '0', '0', st)), (-1 / (2 * h), O1(g, 'm', '0', '0', st)))


def PB(g):
    h = mp.mpf(10) ** -30
    c = 1 / (4 * h * h)
    return lin((c, O12(g, 'p', 'p')), (-c, O12(g, 'p', 'm')), (-c, O12(g, 'm', 'p')), (c, O12(g, 'm', 'm')))


def PC(g):
    D1, D2, _ = GEOM[g]
    c = O12(g, '0', '0', lm=230)
    return [(2 * L + 1) / (D1 * D2) * mp.fsum(c[L + 1::2]) for L in range(151)]


def beta(F, g):
    """beta from the fit beta_eff(L) = beta + c/L + d/L^2 on L = 100..149 (F known to 150)."""
    rho = GEOM[g][2]
    Ls = np.arange(100, 150)
    be = [float(mp.log(abs(F[L + 1] / F[L]) / rho) / mp.log(1 + mp.mpf(1) / L)) for L in Ls]
    sign_flip = any(mp.sign(F[L + 1]) != mp.sign(F[L]) for L in range(100, 150))
    A = np.vstack([np.ones(len(Ls)), 1 / Ls, 1 / Ls ** 2]).T
    return float(np.linalg.lstsq(A, np.array(be), rcond=None)[0][0]), be[-1], sign_flip


out, lines = {}, []


def say(s):
    print(s, flush=True)
    lines.append(s)


# PP-0(a)
a30, a31 = PA('U', 30), PA('U', 31)
dev_a = max(abs(x / y - 1) for x, y in zip(a31, a30))
ok_a = dev_a < mp.mpf('1e-20')
say(f'PP-0(a) PA-U step 1e-31 vs 1e-30, max rel dev L<=150: {mp.nstr(dev_a, 3)}  {"PASS" if ok_a else "FAIL"}')
# PP-0(b)
ok_b = {}
for g, tr in (('U', -4), ('E', -8)):
    b0, raw, fl = beta(q(g, '0', '0', '0', '0'), g)
    ok_b[g] = abs(b0 - tr) < 0.02
    say(f'PP-0(b) 1s at {g}: beta {b0:.5f} (raw at 149 {raw:.4f}, sign flips {fl}), true {tr}  '
        f'{"PASS" if ok_b[g] else "FAIL"}')
    out[f'1s-{g}'] = b0
# PP1
held = {}
for item, f in (('PA-U', lambda: PA('U')), ('PB-U', lambda: PB('U')), ('PC-U', lambda: PC('U')),
                ('PA-E', lambda: PA('E')), ('PB-E', lambda: PB('E'))):
    b0, raw, fl = beta(f(), item[-1])
    held[item] = abs(b0 - PRED[item]) < 0.02
    out[item] = b0
    say(f'PP1 {item}: beta {b0:.5f} (raw at 149 {raw:.4f}, sign flips {fl}), predicted {PRED[item]}  '
        f'{"HELD" if held[item] else "MISSED"}')
pp1 = all(held.values())
pp2 = out['PA-E'] < -9.5 and out['PB-E'] < -9.5
say(f'PP2 PA-E {out["PA-E"]:.4f}, PB-E {out["PB-E"]:.4f} (< -9.5)  {"HELD" if pp2 else "MISSED"}')
nmiss = sum(not v for v in held.values())
if not ok_a or not all(ok_b.values()) or not pp2 or nmiss >= 2:
    gate = 'STOP'
elif pp1:
    gate = 'GO'
elif nmiss == 1:
    gate = 'BORDERLINE'
else:
    gate = 'STOP'
say(f'GATE {gate}')
json.dump(dict(betas=out, lines=lines, gate=gate), open(HERE / 'data' / 'pp_score.json', 'w'), indent=1)
