"""Record PP controls, measured on NON-scored cases before registering.
  C1  pp_exact.coefficients2 with z1 = z2 reproduces ao_exact.coefficients term by term (L = 0..12).
  C2  z1 != z2 against the independent float quadrature F_combo of debug/wa7_score.py (normalized 1s, D1 = D2 = 1.809).
  C3  cost and finite-difference stability of one D1 derivative at a control geometry, steps h and h/10.
    python debug/pp_controls.py
"""
import sys
import time
from fractions import Fraction

import mpmath as mp

sys.argv = sys.argv[:1]
import ao_exact
import na_exact
import pp_exact

K = na_exact.Field(None)


def num(Q, dps):
    mp.mp.dps = dps
    cv, out = {}, []
    for q in Q:
        s = mp.mpf(0)
        for k, v in q.items():
            if k not in cv:
                cv[k] = ao_exact.value(k, mp)
            s += mp.mpf(int(v.p)) / int(v.q) * cv[k]
        out.append(s)
    return out


# C1
a = ao_exact.coefficients('1.7', '2.5', '1.4', '1.1', 12, K)
b = pp_exact.coefficients2('1.7', '2.5', '1.4', '1.1', '1.1', 12, K)
print('C1 identical term by term:', a == b, flush=True)

# C2
from pathlib import Path   # noqa: E402
import numpy as np   # noqa: E402
from numpy.polynomial.legendre import leggauss   # noqa: E402
from scipy.special import spherical_jn   # noqa: E402
sys.path.insert(0, str(Path(__file__).parents[1] / "tests"))
from test_paper59_water_bond_angle import _G, _tanh_sinh   # noqa: E402


def F_combo(zX, H1, H2, lmax, level, D=1.809):
    """debug/wa7_score.py F_combo, unchanged except D as an argument (that file runs its scorer on import)."""
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


zX, z1, z2 = 1.4, 1.1, 1.6
Fq = F_combo(zX, [(z1, 1.0)], [(z2, 1.0)], 8, 1)
Q = pp_exact.coefficients2('1.809', '1.809', '1.4', '1.1', '1.6', 8, K)
Fe = [v * 4 * zX ** 3 * (z1 * z2) ** 1.5 for v in num(Q, 120)]
for L in range(9):
    print('C2 L', L, 'exact', mp.nstr(Fe[L], 12), 'quad', '%.12g' % Fq[L], 'rel', '%.1e' % float(Fq[L] / Fe[L] - 1),
          flush=True)

# C3: d/dD1 at the control geometry (1.7, 2.5, 1.4, 1.1, 1.1), L up to 150, steps 1e-30 and 1e-31
LM = 150
dps = 6 * LM + 260
for h in (Fraction(1, 10 ** 30), Fraction(1, 10 ** 31)):
    t0 = time.time()
    p = num(pp_exact.coefficients2(Fraction(17, 10) + h, '2.5', '1.4', '1.1', '1.1', LM, K), dps)
    m = num(pp_exact.coefficients2(Fraction(17, 10) - h, '2.5', '1.4', '1.1', '1.1', LM, K), dps)
    d = [(x - y) / (2 * h.numerator / mp.mpf(h.denominator)) for x, y in zip(p, m)]
    print('C3 h', h, 'two runs', round(time.time() - t0), 's', flush=True)
    if h == Fraction(1, 10 ** 30):
        d30 = d
    else:
        for L in (10, 50, 100, 150):
            print('C3 L', L, 'dF/dD1', mp.nstr(d[L], 12), 'step agreement', mp.nstr(abs(d[L] / d30[L] - 1), 3),
                  flush=True)
