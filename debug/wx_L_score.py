"""WX part L scorer (record docs/predictions/2026-10-28_water_landau_gaussian_fresh.md).  Controls first.
Reads debug/data/wx_L_F_<case>.json (levels 3 and 4).  -> debug/data/wx_L_score.{json,log}"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "tests"))
sys.path.insert(0, str(HERE))
from test_paper59_water_bond_angle import F_coeffs, REFEREE_NG16          # noqa: E402
from numpy.polynomial.legendre import legval                               # noqa: E402
from wa6_local import singular_part, ax, az, wx, wz, Rw                    # noqa: E402
from wl_local import phi, naive                                            # noqa: E402
from wl_extrap import leg                                                  # noqa: E402

mp.mp.dps = 40
LS = list(range(16, 41, 2))
CASES = {'P': (4, 1.0, 1.6, 2.2), 'Q': (3, 1.5, 2.0, 1.5), 'R': (3, 1.24, 1.5, 2.1), 'S': (4, 1.5, 2.1, 1.6)}
log = []


def say(s):
    print(s, flush=True)
    log.append(s)


def basis(xs, nb):
    xs = mp.mpf(xs)
    return {l: [leg(lambda x, a=mp.mpf(7) / 2 + j: (xs - x) ** a, l) for j in range(nb)] for l in LS}


def extract(F, B, nb):
    A = mp.matrix([[B[l][j] / B[l][0] for j in range(nb)] for l in LS])
    y = mp.matrix([mp.mpf(F[l]) / B[l][0] for l in LS])
    return mp.lu_solve(A.T * A, A.T * y)[0]


# L0(a)
c00 = sp.pi ** 2 / 630
one = sp.Integer(1)
mono = [((az, one, wz), sp.Rational(-2, 5)), ((one, az, wz), sp.Rational(2, 5)), ((ax, ax, wx ** 2), sp.Rational(-2, 15)),
        ((az ** 2, one, wz ** 2), sp.Rational(1, 5)), ((az ** 2, az, wz ** 3), sp.Rational(2, 35))]
ok_a = True
for (g1, g2, m), want in mono:
    P = sp.expand(singular_part(g1, g2, maxdeg=4)).coeff(Rw, 7)
    got = sp.simplify(sp.Poly(P, *sp.symbols('wx wy wz', positive=True)).coeff_monomial(m) / c00)
    ok_a &= (got == want)
sym = phi(1, 1, 1.809, 1.809) / naive(1, 1, 1.809, 1.809)
ok_a &= abs(sym - 1) < mp.mpf('1e-20')
say(f"WX-L0(a) moments exact and symmetric Phi/C - 1 = {mp.nstr(sym - 1, 3)}: {'PASS' if ok_a else 'FAIL'}")
# L0(b)
Yv, Zv = np.array([0, 0, 2.0]), np.array([1.5, 0, 0.5])
Fr = F_coeffs(2.0, np.linalg.norm(Zv), 1.0, 1.0, 20)
diff = abs(legval(Yv @ Zv / (2.0 * np.linalg.norm(Zv)), Fr) - REFEREE_NG16)
ok_b = diff < 1e-9
say(f"WX-L0(b) referee diff {diff:.1e}: {'PASS' if ok_b else 'FAIL'}")

# L0(c), (d) for every case before any comparison
data, ctrl = {}, {}
for name, (zX, zH, D1, D2) in CASES.items():
    d = json.loads((HERE / "data" / f"wx_L_F_{name}.json").read_text())
    F3, F4 = d["level3"], d["level4"]
    noise = max(abs(F4[l] / F3[l] - 1) for l in LS)
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    B5 = basis(xs, 5)
    xm = mp.mpf(xs)
    fsyn = lambda x: (xm - x) ** mp.mpf(3.5) * (1 + 5 * (xm - x) - 6 * (xm - x) ** 2 + 3 * (xm - x) ** 4) + mp.mpf('0.3') * mp.exp(2 * x)
    syn = extract({l: leg(fsyn, l) for l in LS}, {l: v[:4] for l, v in B5.items()}, 4)
    ctrl[name] = dict(noise34=noise, void=not noise < 1e-5, synthetic_err=float(syn - 1), ok_d=bool(abs(syn - 1) < 2e-2))
    data[name] = (F3, B5)
    say(f"WX-L0 {name}: (c) level 3 vs 4 = {noise:.2e} {'live' if noise < 1e-5 else 'VOID'}; "
        f"(d) synthetic err {float(syn - 1):.2e} {'PASS' if ctrl[name]['ok_d'] else 'FAIL'}")
live = [k for k in CASES if not ctrl[k]['void']]
ok_d_all = all(ctrl[k]['ok_d'] for k in live)

# comparison
res = {}
for name in live:
    zX, zH, D1, D2 = CASES[name]
    F3, B5 = data[name]
    sc = mp.mpf(2 * D1 * D2) ** mp.mpf(3.5)
    P, N = phi(zX, zH, D1, D2), naive(zX, zH, D1, D2)
    a4 = extract(F3, {l: v[:4] for l, v in B5.items()}, 4) / (P * sc)
    a5 = extract(F3, B5, 5) / (P * sc)
    an = a5 * P / N
    r = dict(ctrl[name], A4=float(a4), A5=float(a5), A_naive=float(an))
    r['La'] = bool(abs(a5 - 1) < 1e-3 and abs(a4 - a5) < 2e-3)
    r['Lb'] = bool(abs(an - 1) > 10 * abs(a5 - 1))
    res[name] = r
    say(f"{name}: A5 = {float(a5):.7f}, A4 = {float(a4):.7f}, A_naive = {float(an):.5f}; "
        f"La {'held' if r['La'] else 'MISSED'}, Lb {'held' if r['Lb'] else 'MISSED'}")
La = all(r['La'] for r in res.values())
Lb = [r['Lb'] for r in res.values()]
gate = ('STOP' if (not (ok_a and ok_b and ok_d_all) or len(live) < 3 or not La)
        else ('GO' if all(Lb) else ('BORDERLINE' if Lb.count(False) == 1 else 'STOP')))
say(f"live {len(live)}; La {'HELD' if La else 'MISSED'}; Lb {'HELD' if all(Lb) else 'MISSED'}; GATE L {gate}")
(HERE / "data" / "wx_L_score.json").write_text(json.dumps(dict(ctrl=ctrl, res=res, gate=gate), indent=1))
(HERE / "data" / "wx_L_score.log").write_text("\n".join(log) + "\n")
