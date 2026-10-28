"""WL scorer (record docs/predictions/2026-10-28_water_landau_unequal.md). Controls first."""
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
CASES = {'A': (4, 1.24, 1.809, 1.4), 'B': (4, 1.24, 1.809, 2.2), 'C': (4, 1.24, 1.809, 2.6),
         'D': (2.25, 1.24, 1.809, 2.6)}
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


# WL-0(a)
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
say(f"WL-0(a) moments exact and symmetric Phi/C - 1 = {mp.nstr(sym - 1, 3)}: {'PASS' if ok_a else 'FAIL'}")
# WL-0(b)
Yv, Zv = np.array([0, 0, 2.0]), np.array([1.5, 0, 0.5])
Fr = F_coeffs(2.0, np.linalg.norm(Zv), 1.0, 1.0, 20)
diff = abs(legval(Yv @ Zv / (2.0 * np.linalg.norm(Zv)), Fr) - REFEREE_NG16)
ok_b = diff < 1e-9
say(f"WL-0(b) referee diff {diff:.1e}: {'PASS' if ok_b else 'FAIL'}")

res = {}
for name, (zX, zH, D1, D2) in CASES.items():
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    F1 = F_coeffs(D1, D2, zX, zH, 40, level=1)
    F2 = F_coeffs(D1, D2, zX, zH, 40, level=2)
    noise = max(abs(F2[l] / F1[l] - 1) for l in LS)
    void = not noise < 1e-7
    B5 = basis(xs, 5)
    B4 = {l: B5[l][:4] for l in LS}
    # (d) synthetic control at this x*
    xm = mp.mpf(xs)
    fsyn = lambda x: (xm - x) ** mp.mpf(3.5) * (1 + 5 * (xm - x) - 6 * (xm - x) ** 2 + 3 * (xm - x) ** 4) + mp.mpf('0.3') * mp.exp(2 * x)
    Fsyn = {l: leg(fsyn, l) for l in LS}
    syn = extract(Fsyn, B4, 4)
    ok_d = abs(syn - 1) < 2e-2
    scale = mp.mpf(2 * D1 * D2) ** mp.mpf(3.5)
    P, N = phi(zX, zH, D1, D2), naive(zX, zH, D1, D2)
    a4 = extract(F2, B4, 4) / (P * scale)
    a5 = extract(F2, B5, 5) / (P * scale)
    an = a4 * P / N
    r = dict(void=void, noise=noise, synthetic_err=float(syn - 1), ok_d=bool(ok_d), A4=float(a4), A5=float(a5),
             A_naive=float(an), raw_r40=float(F2[40] / (P * scale * B5[40][0])))
    r['WLa'] = bool(abs(a4 - 1) < 2e-2 and abs(a4 - a5) < 1e-2)
    r['WLb'] = bool(abs(an - 1) > 3 * abs(a4 - 1))
    res[name] = r
    say(f"{name}: " + json.dumps(r))
live = [r for r in res.values() if not r['void']]
ok_d_all = all(r['ok_d'] for r in live)
A = all(r['WLa'] for r in live)
Bc = [r['WLb'] for r in live]
gate = ('STOP' if (not (ok_a and ok_b and ok_d_all) or len(live) < 3 or not A)
        else ('GO' if all(Bc) else ('BORDERLINE' if Bc.count(False) == 1 else 'STOP')))
say(f"WL-0(d) {'PASS' if ok_d_all else 'FAIL'}; live cases {len(live)}; WLa {'HELD' if A else 'MISSED'}; "
    f"WLb {'HELD' if all(Bc) else 'MISSED'}; GATE {gate}")
(HERE / "data" / "wl_score.json").write_text(json.dumps(res, indent=1))
(HERE / "data" / "wl_score.log").write_text("\n".join(log) + "\n")
