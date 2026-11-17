"""WX phase 0 (before registering; no scored case): measure the controls of the Gaussian part on NEIGHBOUR cases
(each fresh case with D moved by +0.05 bohr) and on the old symmetric Slater sets, so every control bound in the fresh
record is measured, not guessed.
  slater <set>   : Slater S_l, F_coeffs level 2 vs level 3, l = 0..63 (set: E1old, E3n, E4n)
  kspace <case>  : k-space route vs series form at 50 digits (case: g3n, g4n, E3-3n, E4-3n)
Usage: python debug/wx_measure_G.py slater E4n  -> debug/data/wx_measure_G_<mode>_<name>.json"""
import json
import sys
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "tests"))
sys.path.insert(0, str(HERE))
import wg_common as WC        # noqa: E402
import wx_gauss as G          # noqa: E402

mode, name = sys.argv[1], sys.argv[2]
out = {}
if mode == 'slater':
    from test_paper59_water_bond_angle import F_coeffs
    zX, zH, D = {'E1old': (1.0, 1.0, 1.809), 'E3n': (1.5, 1.0, 2.25), 'E4n': (3.0, 1.5, 1.65)}[name]
    F2 = F_coeffs(D, D, zX, zH, 63, level=2)
    F3 = F_coeffs(D, D, zX, zH, 63, level=3)
    dev = {l: float(abs(F3[l] / F2[l] - 1)) for l in range(64)}
    out = dict(dev=dev, max=max(dev.values()), min_F2=float(min(F2)))
else:
    mp.mp.dps = 50
    base = dict(G.CASES)
    nb = {'g3n': 'g3', 'g4n': 'g4', 'E3-3n': 'E3-3', 'E4-3n': 'E4-3'}[name]
    k, a, b, D = base[nb]
    spec = (k, a, b, str(mp.mpf(D) + mp.mpf('0.05')))
    ls = [0, 5, 10, 20, 30] if k == 'single' else [0, 10, 20]
    ser = {l: mp.mpf(0) for l in ls}
    web = {l: mp.mpf(0) for l in ls}
    for Q in G.quartets(spec):
        Fq = WC.F_quartet(Q, max(ls))
        for l in ls:
            ser[l] += Fq[l]
            web[l] += WC.F_weber_quartet(Q, l)
    dev = {l: float(abs(web[l] / ser[l] - 1)) for l in ls}
    out = dict(spec=spec, dev=dev, max=max(dev.values()))
print(name, json.dumps(out), flush=True)
(HERE / "data" / f"wx_measure_G_{mode}_{name}.json").write_text(json.dumps(out, indent=1))
