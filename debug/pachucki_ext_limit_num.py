"""Numerical check of the a -> b limit (debug/data/pachucki_ext_limit.json) and of Sugiura's formula against
geovac.qfd_core.exchange_hp at a = b = 1, R = 1.4."""
import json, sys
import mpmath as mp
import sympy as sp
sys.path.insert(0, ".")
from geovac import qfd_core as Q
mp.mp.dps = 40
d = json.load(open("debug/data/pachucki_ext_limit.json"))
R, e2, Lg, X2, X4 = sp.symbols("R e2 Lg X2 X4")
Rv = mp.mpf("1.4")
vals = {R: Rv, e2: mp.exp(-2 * Rv), Lg: mp.euler + mp.log(Rv), X2: mp.e1(2 * Rv), X4: mp.exp(2 * Rv) * mp.e1(4 * Rv)}
lim = sp.sympify(d["limit"])
S = 1 + R + R ** 2 / 3; Sp = 1 - R + R ** 2 / 3
sug = sp.Rational(1, 5) * (-e2 * (-sp.Rational(25, 8) + sp.Rational(23, 4) * R + 3 * R ** 2 + R ** 3 / 3)
                           + 6 / R * (e2 * S ** 2 * Lg - Sp ** 2 * X4 + 2 * S * Sp * X2))
f = lambda e: sp.lambdify([R, e2, Lg, X2, X4], e, "mpmath")(*[vals[s] for s in (R, e2, Lg, X2, X4)])
vh, per = Q.exchange_hp(sp.Integer(1), (1, 0, 0), (1, 0, 0), sp.Integer(1), (1, 0, 0), (1, 0, 0), sp.Rational(7, 5), tau_max=24, dps=40)
print("limit  :", mp.nstr(f(lim), 35))
print("Sugiura:", mp.nstr(f(sug), 35))
print("hp     :", mp.nstr(vh, 35), " last terms", [mp.nstr(x, 3) for x in per[-3:]])
