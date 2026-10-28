"""WG scorer (record docs/predictions/2026-10-28_water_gaussian_contact.md).  Reads the control outputs
debug/data/wg_ctrl_*.json and the scored F_l in debug/data/wg_F_<case>.json; applies the registered clauses."""
import json
from pathlib import Path

import mpmath as mp

mp.mp.dps = 50
OUT = Path(__file__).parent / "data"
log = []


def say(s):
    print(s)
    log.append(str(s))


def ctrl(name):
    return json.loads((OUT / f"wg_ctrl_{name}.json").read_text())


def ok(v):
    return v is True or v == "True"


A = ctrl('a')
c = ctrl('c')
d = ctrl('d')
e = ctrl('e')
pa = all(ok(A[s]['pass_']) for s in A)
pc = all(ok(c[k]['pass_']) for k in c)
pd = ok(d['pass_'])
pe = all(ok(e[s]['pass_']) for s in e)
say(f"WG-0(a) {A['E1']['max']:.1e} {A['E2']['max']:.1e} {'PASS' if pa else 'FAIL'}")
say(f"WG-0(c) g1 {c['g1']['dev']:.1e}, E1-3 {c['E1-3']['dev']:.1e} {'PASS' if pc else 'FAIL'}")
say(f"WG-0(d) resid {d['resid']:.2e} {'PASS' if pd else 'FAIL'}")
say(f"WG-0(e) " + json.dumps({s: e[s]['p'] for s in e}) + f" {'PASS' if pe else 'FAIL'}")
b = {k: ctrl(f"b_{k}")[k] for k in ('g1', 'g2', 'E1-3', 'E1-6')}
for k, v in b.items():
    say(f"WG-0(b) {k}: max {v['max']:.2e} {'PASS' if ok(v['pass_']) else 'FAIL'}  {v['dev']}")
live = {'g1': ok(b['g1']['pass_']), 'g2': ok(b['g2']['pass_'])}
for s, n3, n6 in (('E1', 'E1-3', 'E1-6'), ('E2', 'E2-3', 'E2-6')):
    # (b) was run on the E1 contracted cases only; a failure there voids that set's two cases (record).
    pass_set = (ok(b['E1-3']['pass_']) and ok(b['E1-6']['pass_'])) if s == 'E1' else True
    live[n3] = live[n6] = pass_set and pa

F = {k: [mp.mpf(x) for x in json.loads((OUT / f"wg_F_{k}.json").read_text())] for k in live}
stored = json.load(open(OUT / "wa5_cusp.json"))
S = {'E1': stored["E1 (1.809, 1.809)"]["F"], 'E2': stored["E2 (1.809, 1.809)"]["F"]}


def pw(Fl, l):
    return float(mp.log(Fl[l] / Fl[l + 4]) / mp.log(mp.mpf(l + 4) / l))


res = {}
# WGa
QC = {'g1': (1.6362405, [-2.305561, 9.222245, -18.70707, -4.149909]),
      'g2': (5.6912713, [-13.78891, 55.15565, 165.2682, -2400.561])}
WGa = True
for k in ('g1', 'g2'):
    Fl = F[k]
    T1 = mp.mpf(QC[k][0])
    from wg_common import quartets
    T1 = quartets(k)[0]['T1']
    pos = all(v > 0 for v in Fl[:61])
    qco = QC[k][1]
    worst = 0
    for l in range(40, 61):
        q = Fl[l + 1] * (2 * l + 3) / (Fl[l] * T1)
        pred = 1 + sum(cj / mp.mpf(l) ** (j + 2) for j, cj in enumerate(qco))
        bound = 0.5 * max(abs(qco[3]), 10) / mp.mpf(l) ** 5
        worst = max(worst, float(abs(q - pred) / bound))
    p30 = pw(Fl, 30)
    held = pos and worst < 1 and p30 > 50
    WGa &= held
    res[k] = dict(positive=pos, worst_resid_over_bound=worst, p30=p30, held=held)
    say(f"WGa {k}: positive {pos}; max |q - expansion| / bound over l = 40..60 = {worst:.3g}; p_30 = {p30:.1f}"
        f" -> {'held' if held else 'MISSED'}")

# WGb1-3, WGc
MODEL = {'E1-3': 11, 'E1-6': 32, 'E2-3': 13, 'E2-6': 40}
lstar = {}
WGb1 = WGb2 = WGb3 = True
for k, ml in MODEL.items():
    if not live[k]:
        say(f"{k}: VOID")
        continue
    Fl, Sl = F[k], S[k[:2]]
    ratio = [float(Fl[l] / Sl[l]) for l in range(61)]
    pos = all(v > 0 for v in Fl[:61])
    ls = next((l for l in range(6, 61) if ratio[l] < 0.5), None)
    lstar[k] = ls
    inwin = ls is not None and abs(ls - ml) <= 3
    l0 = ml + 10
    p0, p5 = pw(Fl, l0), pw(Fl, l0 + 5)
    b3 = p0 > 16 and p5 > p0
    WGb1 &= pos
    WGb2 &= inwin
    WGb3 &= b3
    r = dict(positive=pos, lstar=ls, model=ml, in_window=inwin, l0=l0, p_l0=p0, p_l0p5=p5, b3=b3,
             ratio={l: ratio[l] for l in range(61)})
    if k.endswith('-6'):
        r['c1_maxdev_0_12'] = max(abs(ratio[l] - 1) for l in range(13))
        r['c1'] = r['c1_maxdev_0_12'] < 0.03
        lo, hi = ml - 14, ml - 4
        lm = max(range(lo, hi + 1), key=lambda l: ratio[l])
        r['c2_max'] = (lm, ratio[lm])
        r['c2'] = ratio[lm] > 1.08
    res[k] = r
    say(f"{k}: positive {pos}; l*_obs {ls} (model {ml}, window {ml - 3}..{ml + 3}) {'in' if inwin else 'OUT'}; "
        f"p_{l0} = {p0:.2f}, p_{l0 + 5} = {p5:.2f} {'held' if b3 else 'MISSED'}"
        + (f"; c1 max|ratio-1| l<=12 = {r['c1_maxdev_0_12']:.3g}; c2 max ratio {r['c2_max'][1]:.4f} at l = {r['c2_max'][0]}"
           if k.endswith('-6') else ''))
    say("   ratio: " + ", ".join(f"{l}:{ratio[l]:.4g}" for l in range(0, 61)))
for s in ('E1', 'E2'):
    a, bb = lstar.get(f"{s}-3"), lstar.get(f"{s}-6")
    if a and bb:
        rr = bb / a
        okr = 2.4 <= rr <= 3.8
        WGb2 &= okr
        res[f"{s}_scaling"] = dict(ratio=rr, held=okr)
        say(f"{s}: l*(6)/l*(3) = {rr:.3f} {'in 2.4..3.8' if okr else 'OUT'}")
ncon = sum(live[k] for k in MODEL)
controls = pa and pc and pd and pe and live['g1'] and live['g2'] and ncon >= 3
if not controls or not (WGa and WGb1 and WGb3):
    gate = 'STOP'
elif WGb2:
    gate = 'GO'
else:
    gate = 'BORDERLINE'
say(f"WGa {'HELD' if WGa else 'MISSED'}  WGb1 {'HELD' if WGb1 else 'MISSED'}  WGb2 {'HELD' if WGb2 else 'MISSED'}  "
    f"WGb3 {'HELD' if WGb3 else 'MISSED'}  live contracted {ncon}  GATE {gate}")
(OUT / "wg_score.json").write_text(json.dumps(res, indent=1, default=str))
(OUT / "wg_score.log").write_text("\n".join(log) + "\n")
