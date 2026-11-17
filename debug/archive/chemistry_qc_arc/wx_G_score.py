"""WX part G scorer (record docs/predictions/2026-10-28_water_landau_gaussian_fresh.md).  Controls first; reads the
files written by wx_G_run.py.  -> debug/data/wx_G_score.{json,log}"""
import json
import math
from pathlib import Path

import mpmath as mp

import wx_gauss as G

mp.mp.dps = 50
OUT = Path(__file__).parent / "data"
log = []


def say(s):
    print(s, flush=True)
    log.append(str(s))


def rd(name):
    return json.loads((OUT / name).read_text())


# controls
S, okA = {}, {}
for s in ('E3', 'E4'):
    d = rd(f"wx_G_S_{s}.json")
    dev = max(abs(d['level3'][l] / d['level2'][l] - 1) for l in range(64))
    okA[s] = dev < 2e-5
    S[s] = d['level2']
    say(f"WX-G0(a) {s}: Slater level 2 vs 3 max {dev:.2e} {'PASS' if okA[s] else 'FAIL'}")
b = {k: rd(f"wx_G_ctrl_b_{k}.json") for k in ('g3', 'g4', 'E3-3', 'E4-3')}
okB = {k: v['max'] < 1e-12 for k, v in b.items()}
for k, v in b.items():
    say(f"WX-G0(b) {k}: max {v['max']:.2e} {'PASS' if okB[k] else 'FAIL'}")
c = rd("wx_G_ctrl_c.json")
okC = all(v < 1e-14 for v in c.values())
say(f"WX-G0(c) {c} {'PASS' if okC else 'FAIL'}")
d = rd("wx_G_ctrl_d.json")
okD = d['resid'] < 1e-7
say(f"WX-G0(d) resid {d['resid']:.2e} {'PASS' if okD else 'FAIL'}")


def pw(Fl, l):
    return math.log(Fl[l] / Fl[l + 4]) / math.log((l + 4) / l)


pe = {s: {l: pw(S[s], l) for l in (21, 23, 42, 50)} for s in S}
okE = all(7 < v < 9 for s in pe for v in pe[s].values())
say(f"WX-G0(e) {pe} {'PASS' if okE else 'FAIL'}")

live = {'g3': okB['g3'], 'g4': okB['g4']}
for s in ('E3', 'E4'):
    live[f"{s}-3"] = live[f"{s}-6"] = okA[s] and okB[f"{s}-3"]
say(f"live: {live}")

# comparison
F = {k: [mp.mpf(x) for x in rd(f"wx_G_F_{k}.json")] for k in live if live[k]}
res = {}
Ga = True
for k in ('g3', 'g4'):
    if not live[k]:
        Ga = False
        continue
    Fl = F[k]
    T1 = G.quartets(G.CASES[k])[0]['T1']
    cq = G.q_coeffs(G.single_T(G.CASES[k]))
    pos = all(v > 0 for v in Fl)
    worst = 0
    for l in range(40, 61):
        q = Fl[l + 1] * (2 * l + 3) / (Fl[l] * T1)
        pred = 1 + sum(cj / mp.mpf(l) ** (j + 2) for j, cj in enumerate(cq))
        worst = max(worst, float(abs(q - pred) / (0.5 * max(abs(cq[3]), 10) / mp.mpf(l) ** 5)))
    p30 = float(mp.log(Fl[30] / Fl[34]) / mp.log(mp.mpf(34) / 30))
    held = pos and worst < 1 and p30 > 50
    Ga &= held
    res[k] = dict(positive=pos, worst_over_bound=worst, p30=p30, held=held)
    say(f"WX-Ga {k}: positive {pos}; worst/bound {worst:.3g}; p_30 {p30:.1f} -> {'held' if held else 'MISSED'}")

MODEL = {'E3-3': 13, 'E3-6': 39, 'E4-3': 14, 'E4-6': 43}
lstar = {}
Gb1 = Gb2 = Gb3 = True
for k, ml in MODEL.items():
    if not live[k]:
        say(f"{k}: VOID")
        continue
    Fl, Sl = F[k], S[k[:2]]
    ratio = [float(Fl[l] / Sl[l]) for l in range(64)]
    pos = all(v > 0 for v in Fl)
    ls = next((l for l in range(6, 64) if ratio[l] < 0.5), None)
    lstar[k] = ls
    inwin = ls is not None and ml - 2 <= ls <= ml + 1
    l0 = ml + 10
    p0, p5 = float(mp.log(Fl[l0] / Fl[l0 + 4]) / mp.log(mp.mpf(l0 + 4) / l0)), \
        float(mp.log(Fl[l0 + 5] / Fl[l0 + 9]) / mp.log(mp.mpf(l0 + 9) / (l0 + 5)))
    b3 = p0 > 16 and p5 > p0
    Gb1 &= pos
    Gb2 &= inwin
    Gb3 &= b3
    res[k] = dict(positive=pos, lstar=ls, model=ml, in_window=inwin, l0=l0, p_l0=p0, p_l0p5=p5, b3=b3,
                  ratio={l: ratio[l] for l in range(64)})
    say(f"{k}: positive {pos}; l*_obs {ls} (model {ml}, window {ml - 2}..{ml + 1}) {'in' if inwin else 'OUT'}; "
        f"p_{l0} = {p0:.2f}, p_{l0 + 5} = {p5:.2f} {'held' if b3 else 'MISSED'}")
    say("   ratio: " + ", ".join(f"{l}:{ratio[l]:.4g}" for l in range(64)))
for s in ('E3', 'E4'):
    a, bb = lstar.get(f"{s}-3"), lstar.get(f"{s}-6")
    if a and bb:
        rr = bb / a
        okr = 2.9 <= rr <= 3.5
        Gb2 &= okr
        res[f"{s}_scaling"] = dict(ratio=rr, held=okr)
        say(f"{s}: l*(6)/l*(3) = {rr:.3f} {'in 2.9..3.5' if okr else 'OUT'}")
ncon = sum(live[k] for k in MODEL)
controls = okC and okD and okE and live['g3'] and live['g4'] and ncon >= 3
gate = 'STOP' if (not controls or not (Ga and Gb1 and Gb3)) else ('GO' if Gb2 else 'BORDERLINE')
say(f"Ga {'HELD' if Ga else 'MISSED'}  Gb1 {'HELD' if Gb1 else 'MISSED'}  Gb2 {'HELD' if Gb2 else 'MISSED'}  "
    f"Gb3 {'HELD' if Gb3 else 'MISSED'}  live contracted {ncon}  GATE G {gate}")
(OUT / "wx_G_score.json").write_text(json.dumps(dict(res=res, live=live, gate=gate), indent=1, default=str))
(OUT / "wx_G_score.log").write_text("\n".join(log) + "\n")
