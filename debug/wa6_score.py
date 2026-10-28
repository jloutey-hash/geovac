"""WA6 scorer (record docs/predictions/2026-10-27_water_cusp_next_term.md)."""
import json
import subprocess
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

mp.mp.dps = 40
DATA = Path(__file__).parent / "data" / "wa5_cusp.json"
ALPHA = {'E1': (mp.mpf('0.0305532186374231'), mp.mpf('0.0285025232663417')),
         'E2': (mp.mpf('0.0110819491640684'), mp.mpf('0.0134090231083288'))}


def leg_cf(a, l):
    """(2l+1)/2 int_{-1}^{1} (1-x)^a P_l dx, closed form."""
    return (2 * l + 1) * mp.mpf(2) ** a * mp.gamma(a + 1) ** 2 / (mp.gamma(a + l + 2) * mp.gamma(a + 1 - l)) * (-1) ** l


def leg_quad(a, l):
    return (2 * l + 1) / mp.mpf(2) * mp.quad(lambda x: (1 - x) ** a * mp.legendre(l, x), mp.linspace(-1, 1, 13))


a7, a9 = mp.mpf(7) / 2, mp.mpf(9) / 2
ok0 = True
for a in (a7, a9):
    for l in (15, 30, 45, 60):
        rel = abs(leg_cf(a, l) / leg_quad(a, l) - 1)
        ok0 &= rel < 1e-12
        print(f"WA6-0(a) a={float(a)} l={l} rel {float(rel):.1e}")
out = subprocess.run([sys.executable, str(Path(__file__).parent / "wa6_local.py")], capture_output=True, text=True).stdout
line = out.splitlines()[0]
ok_b = "16*R**9/14175 + 2*R**7/315" in line
print("WA6-0(b)", line, "PASS" if ok_b else "FAIL")
ok0 &= ok_b
d = json.loads(DATA.read_text())
res = {}
for name in ('E1', 'E2'):
    row = d[f"{name} (1.809, 1.809)"]
    F, fl = np.array(row['F']), np.array(row['floor'])
    noise = max(fl[l] / abs(F[l]) for l in range(30, 61))
    ok0 &= noise < 1e-4
    print(f"WA6-0(c) {name} max noise 30..60 {noise:.1e}")
    A7, A9 = ALPHA[name]
    L0 = {l: float(A7 * leg_cf(a7, l)) for l in range(15, 61)}
    L2 = {l: float(A7 * leg_cf(a7, l) + A9 * leg_cf(a9, l)) for l in range(15, 61)}
    dev2 = {l: abs(F[l] / L2[l] - 1) for l in L2}
    dev0 = {l: abs(F[l] / L0[l] - 1) for l in L0}
    res[name] = dict(dev2={l: dev2[l] for l in (15, 20, 30, 40, 50, 60)}, dev0={l: dev0[l] for l in (15, 30, 60)},
                     max50_60=max(dev2[l] for l in range(50, 61)),
                     ratio30_60=dev2[30] / dev2[60], gain60=dev0[60] / dev2[60],
                     signed={l: F[l] / L2[l] - 1 for l in (30, 40, 50, 60)})
    r = res[name]
    a_ok = dev2[60] < 1e-3 and r['max50_60'] < 2e-3
    b_ok = r['ratio30_60'] > 6 and r['gain60'] >= 10
    print(name, json.dumps(r, default=float, indent=None))
    print(f"   WA6a {'HELD' if a_ok else 'MISSED'}   WA6b {'HELD' if b_ok else 'MISSED'}")
print("WA6-0", "PASS" if ok0 else "FAIL")
(Path(__file__).parent / "data" / "wa6_score.json").write_text(json.dumps(res, default=float, indent=1))
