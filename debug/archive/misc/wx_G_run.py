"""WX part G runs (record docs/predictions/2026-10-28_water_landau_gaussian_fresh.md).  Writes data only.
  slater <E3|E4>  : S_l at levels 2 and 3, l = 0..63           -> data/wx_G_S_<set>.json   (controls a, e; S_l)
  kspace <case>   : control (b), 50 digits                      -> data/wx_G_ctrl_b_<case>.json
  resum           : control (c), 50 digits, g3 and E3-3         -> data/wx_G_ctrl_c.json
  toy             : control (d)                                 -> data/wx_G_ctrl_d.json
  F <case>        : scored F_l, l = 0..63, series, 50 digits    -> data/wx_G_F_<case>.json"""
import json
import sys
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
OUT = HERE / "data"
sys.path.insert(0, str(HERE.parent / "tests"))
sys.path.insert(0, str(HERE))
import wg_common as WC        # noqa: E402
import wx_gauss as G          # noqa: E402

mp.mp.dps = 50
mode = sys.argv[1]
arg = sys.argv[2] if len(sys.argv) > 2 else None
if mode == 'slater':
    from test_paper59_water_bond_angle import F_coeffs
    zX, zH, D = (float(v) for v in G.SETS[arg])
    F2 = F_coeffs(D, D, zX, zH, 63, level=2)
    F3 = F_coeffs(D, D, zX, zH, 63, level=3)
    (OUT / f"wx_G_S_{arg}.json").write_text(json.dumps(dict(level2=[float(v) for v in F2], level3=[float(v) for v in F3])))
elif mode == 'kspace':
    spec = G.CASES[arg]
    ls = [0, 5, 10, 20, 30] if spec[0] == 'single' else [0, 10, 20]
    ser = {l: mp.mpf(0) for l in ls}
    web = {l: mp.mpf(0) for l in ls}
    for Q in G.quartets(spec):
        Fq = WC.F_quartet(Q, max(ls))
        for l in ls:
            ser[l] += Fq[l]
            web[l] += WC.F_weber_quartet(Q, l)
    dev = {l: float(abs(web[l] / ser[l] - 1)) for l in ls}
    (OUT / f"wx_G_ctrl_b_{arg}.json").write_text(json.dumps(dict(dev=dev, max=max(dev.values())), indent=1))
elif mode == 'resum':
    x = mp.cos(mp.radians(mp.mpf('104.5')))
    res = {}
    for case in ('g3', 'E3-3'):
        F = G.F_case(G.CASES[case], 60)
        s = sum(F[l] * mp.legendre(l, x) for l in range(61))
        res[case] = float(abs(s / G.boys_direct(G.CASES[case], x) - 1))
    (OUT / "wx_G_ctrl_c.json").write_text(json.dumps(res, indent=1))
elif mode == 'toy':
    T = mp.mpf('1.2')
    F = WC.F_quartet(dict(pref=mp.mpf(1), T0=T, T1=T), 41)
    l = 40
    q = F[l + 1] * (2 * l + 3) / (F[l] * T)
    c = G.q_coeffs(T)
    r = float(abs(q - 1 - sum(cj / mp.mpf(l) ** (j + 2) for j, cj in enumerate(c))))
    (OUT / "wx_G_ctrl_d.json").write_text(json.dumps(dict(resid=r), indent=1))
elif mode == 'F':
    F = G.F_case(G.CASES[arg], 63)
    (OUT / f"wx_G_F_{arg}.json").write_text(json.dumps([mp.nstr(v, 40) for v in F]))
print(mode, arg, "done", flush=True)
