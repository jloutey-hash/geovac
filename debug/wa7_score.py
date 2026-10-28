"""WA7 scorer (record docs/predictions/2026-10-27_water_zero_contact.md). Controls first.
F_l by the channel formula of tests/test_paper59_water_bond_angle.py, made bilinear in a two-term H function."""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import spherical_jn

sys.path.insert(0, str(Path(__file__).parents[1] / "tests"))
sys.path.insert(0, str(Path(__file__).parent))
from test_paper59_water_bond_angle import _G, _tanh_sinh   # noqa: E402

mp.mp.dps = 40
D = 1.809
FREE = [(1.0, 1.0), (2.0, -2.0 ** -2.5)]     # (exponent, coefficient) on normalized 1s functions
ONE_S = [(1.0, 1.0)]
BETA = {('both', 1.0): (11, '0.003966041796557', '0.004938449424033'),
        ('both', 2.25): (11, '0.000490690948286', '0.000783337174502'),
        ('one', 1.0): (9, '-0.013331310330642', '-0.015077280231292'),
        ('one', 2.25): (9, '-0.001649390915071', '-0.002273103875760')}


def F_combo(zX, H1, H2, lmax, level):
    t, wt = _tanh_sinh(60 * 2 ** level, 0.07 / 2 ** level)
    u, wu = leggauss(400 * 2 ** level)
    u, wu = 0.5 * (u + 1), 0.5 * wu
    k = 3.0 * u / (1 - u)
    wk = wu * 3.0 / (1 - u) ** 2
    G1 = sum(c * _G(k[:, None], t[None, :], D, zX, z) for z, c in H1)
    G2 = sum(c * _G(k[:, None], t[None, :], D, zX, z) for z, c in H2)
    out = []
    for l in range(lmax + 1):
        J = spherical_jn(l, np.outer(k, t) * D)
        out.append(2 * (2 * l + 1) / np.pi * np.sum(wk * ((J * G1) @ wt) * ((J * G2) @ wt)))
    return np.array(out)


def leg_cf(a, l):
    return (2 * l + 1) * mp.mpf(2) ** a * mp.gamma(a + 1) ** 2 / (mp.gamma(a + l + 2) * mp.gamma(a + 1 - l)) * (-1) ** l


def leg_quad(a, l):
    return (2 * l + 1) / mp.mpf(2) * mp.quad(lambda x: (1 - x) ** a * mp.legendre(l, x), mp.linspace(-1, 1, 13))


log = []
def say(s):
    print(s); log.append(s)


ok_a = True
for a in (mp.mpf(9) / 2, mp.mpf(11) / 2, mp.mpf(13) / 2):
    for l in (15, 30):
        rel = abs(leg_cf(a, l) / leg_quad(a, l) - 1)
        ok_a &= rel < 1e-12
        say(f"WA7-0(a) a={float(a)} l={l} rel {float(rel):.1e}")
import subprocess  # noqa: E402
tc = subprocess.run([sys.executable, str(Path(__file__).parents[1] / "debug" / "wa7_twocentre.py")],
                    capture_output=True, text=True).stdout.strip()
ok_b = tc == "8*R**11/93555 + 16*R**9/14175 + 2*R**7/315"
say(f"WA7-0(b) {tc} {'PASS' if ok_b else 'FAIL'}")
# (d) WA6 E1 1s control
F1s = F_combo(1.0, ONE_S, ONE_S, 30, 2)
a7, a9 = mp.mpf(7) / 2, mp.mpf(9) / 2
L2 = float(mp.mpf('0.0305532186374231') * leg_cf(a7, 30) + mp.mpf('0.0285025232663417') * leg_cf(a9, 30))
dev_d = abs(F1s[30] / L2 - 1)
ok_d = abs(dev_d / 1.6e-3 - 1) < 0.1
say(f"WA7-0(d) E1 1s dev_30 {dev_d:.3e} {'PASS' if ok_d else 'FAIL'}")

res = {}
for (case, zX), (p, b1, b2) in BETA.items():
    H1 = FREE if case == 'both' else ONE_S
    Fa = F_combo(zX, H1, FREE, 35, 1)
    Fb = F_combo(zX, H1, FREE, 35, 2)
    ls = range(15, 36)
    noise = max(abs(Fb[l] / Fa[l] - 1) for l in ls)
    void = not noise < 1e-3
    ap, ap2 = mp.mpf(p) / 2, mp.mpf(p + 2) / 2
    L1 = {l: float(mp.mpf(b1) * leg_cf(ap, l)) for l in ls}
    L2 = {l: float(mp.mpf(b1) * leg_cf(ap, l) + mp.mpf(b2) * leg_cf(ap2, l)) for l in ls}
    dev2 = {l: abs(Fb[l] / L2[l] - 1) for l in ls}
    dev1 = {l: abs(Fb[l] / L1[l] - 1) for l in ls}
    sign_ok = all(np.sign(Fb[l]) == np.sign(float(b1)) for l in ls)
    tol = 5e-2 if case == 'both' else 2e-2
    a_ok = dev2[30] < tol and sign_ok
    b_ok = dev2[15] / dev2[30] > 6 and dev1[30] >= 5 * dev2[30]
    r = dict(void=void, noise=noise, dev2={l: dev2[l] for l in (15, 20, 25, 30, 35)},
             dev1={l: dev1[l] for l in (15, 30, 35)}, ratio15_30=dev2[15] / dev2[30], gain30=dev1[30] / dev2[30],
             sign_ok=sign_ok, F30=Fb[30], WA7a=a_ok, WA7b=b_ok, signed30=Fb[30] / L2[30] - 1)
    res[f"{case} zX={zX}"] = r
    say(f"{case} zX={zX}: " + json.dumps(r, default=float))
ok0 = ok_a and ok_b and ok_d
live = [r for r in res.values() if not r['void']]
say(f"WA7-0 (a,b,d) {'PASS' if ok0 else 'FAIL'}; cases surviving (c): {len(live)}")
A = all(r['WA7a'] for r in live); B = all(r['WA7b'] for r in live)
gate = 'STOP' if (not ok0 or len(live) < 2 or not A) else ('GO' if B else 'BORDERLINE')
say(f"WA7a {'HELD' if A else 'MISSED'}  WA7b {'HELD' if B else 'MISSED'}  GATE {gate}")
out = Path(__file__).parent / "data"
(out / "wa7_score.json").write_text(json.dumps(res, default=float, indent=1))
(out / "wa7_score.log").write_text("\n".join(log) + "\n")
