"""WA8 scorer (record docs/predictions/2026-10-27_water_zero_contact_fresh.md). Controls first.
F_l by the channel formula of tests/test_paper59_water_bond_angle.py, bilinear in a two-term H function (as WA7)."""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import spherical_jn

sys.path.insert(0, str(Path(__file__).parents[1] / "tests"))
from test_paper59_water_bond_angle import _G, _tanh_sinh   # noqa: E402

mp.mp.dps = 40


def free(a, b):
    return [(a, 1.0), (b, -(a / b) ** 2.5)]


def F_combo(zX, H1, H2, D, lmax, level):
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


CASES = {'c1': ('both', 2.25, 1.24, (1.24, 3.0), 1.809, 11, '0.008900840359941', '0.019526793502478', 0.090),
         'c2': ('one', 2.25, 1.24, (1.24, 3.0), 1.809, 9, '-0.012027885734131', '-0.022171191146138', 0.035),
         'c3': ('both', 1.0, 1.0, (1.0, 2.0), 2.4, 11, '0.027259672414193', '0.049175250807846', 0.064),
         'c4': ('one', 1.0, 1.0, (1.0, 2.0), 2.4, 9, '-0.052058402179883', '-0.084104743262045', 0.026)}
log = []


def say(s):
    print(s)
    log.append(s)


F7 = F_combo(1.0, free(1.0, 2.0), free(1.0, 2.0), 1.809, 30, 2)[30]
ref = 1.4372867712994588e-14
ok_a = abs(F7 / ref - 1) < 1e-6
say(f"WA8-0(a) WA7 E1 'both' F_30 {F7:.12e} vs {ref:.12e}: {'PASS' if ok_a else 'FAIL'}")
res = {}
for name, (kind, zX, zH, pair, D, p, b1, b2, bound) in CASES.items():
    H2 = free(*pair)
    H1 = H2 if kind == 'both' else [(zH, 1.0)]
    Fa = F_combo(zX, H1, H2, D, 35, 1)
    Fb = F_combo(zX, H1, H2, D, 35, 2)
    ls = range(20, 36)
    noise = max(abs(Fb[l] / Fa[l] - 1) for l in ls)
    ap, ap2 = mp.mpf(p) / 2, mp.mpf(p + 2) / 2
    L1 = {l: float(mp.mpf(b1) * leg_cf(ap, l)) for l in ls}
    L2 = {l: float(mp.mpf(b1) * leg_cf(ap, l) + mp.mpf(b2) * leg_cf(ap2, l)) for l in ls}
    dev2 = {l: abs(Fb[l] / L2[l] - 1) for l in ls}
    dev1_35 = abs(Fb[35] / L1[35] - 1)
    sign_ok = all(np.sign(Fb[l]) == np.sign(L2[l]) for l in ls)
    r = dict(void=not noise < 1e-3, noise=noise, dev2={l: dev2[l] for l in (20, 25, 30, 35)}, dev1_35=dev1_35,
             signed35=Fb[35] / L2[35] - 1, ratio20_35=dev2[20] / dev2[35], gain35=dev1_35 / dev2[35],
             sign_ok=sign_ok, F35=Fb[35], bound=bound)
    r['WA8a'] = bool(dev2[35] < bound and sign_ok)
    r['WA8b'] = bool(r['ratio20_35'] > 3 and r['gain35'] >= 5)
    res[name] = r
    say(f"{name}: " + json.dumps(r, default=float))
live = [r for r in res.values() if not r['void']]
A = all(r['WA8a'] for r in live)
B = all(r['WA8b'] for r in live)
gate = 'STOP' if (not ok_a or len(live) < 3 or not A) else ('GO' if B else 'BORDERLINE')
say(f"cases surviving (b): {len(live)}; WA8a {'HELD' if A else 'MISSED'}  WA8b {'HELD' if B else 'MISSED'}  GATE {gate}")
out = Path(__file__).parent / "data"
(out / "wa8_score.json").write_text(json.dumps(res, default=float, indent=1))
(out / "wa8_score.log").write_text("\n".join(log) + "\n")
