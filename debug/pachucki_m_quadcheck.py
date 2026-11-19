"""Independent quadrature check of the r12-weighted masters (no Pachucki equation, no Neumann expansion).

Prolate spheroidal coordinates, phi integral in closed form (complete elliptic E), 4-D tensor Gauss quadrature:
    f_2(R)            = R/(16 pi^2) int int psi psi r12 / (r1A r1B r2A r2B)
    (4pi)^2/R f2_1111 = int int psi(1) psi(2) r12 d^3r1 d^3r2,       psi = e^{-a r_A - b r_B}
    f_4 analogues with r12^3 (phi integral by Gauss-Legendre).
Accuracy is limited by the cone point r12 = 0 (expect 4-6 digits).  Run: python debug/pachucki_m_quadcheck.py A B R [N]
"""
import json
import sys
from fractions import Fraction

import mpmath as mp
import numpy as np
from scipy.special import ellipe

sys.path.insert(0, "debug")
sys.path.insert(0, "tests")
import pachucki_ext_core as E  # noqa: E402
import pachucki_m_core as M  # noqa: E402
import _pachucki_exchange as PX  # noqa: E402

a, b, Rq = Fraction(sys.argv[1]), Fraction(sys.argv[2]), Fraction(sys.argv[3])
N = int(sys.argv[4]) if len(sys.argv) > 4 else 56
af, bf, R = float(a), float(b), float(Rq)

xl, wl = np.polynomial.laguerre.laggauss(N)
sc = (R / 2) * (af + bf)
xi = 1 + xl / sc
wxi = wl / sc * np.exp(xl)                    # weight for int_1^oo dxi g(xi)  (g includes its own exponential)
et, wet = np.polynomial.legendre.leggauss(N)
ph, wph = np.polynomial.legendre.leggauss(40)
ph = (ph + 1) * np.pi / 2
wph = wph * np.pi / 2 * 2                     # int_0^{2pi} = 2 int_0^pi

X1, E1 = np.meshgrid(xi, et, indexing="ij")
W1 = np.outer(wxi, wet)
rho = (R / 2) * np.sqrt((X1**2 - 1) * (1 - E1**2))
z = (R / 2) * X1 * E1
ex = np.exp(-(R / 2) * ((af + bf) * X1 + (af - bf) * E1))
jac = (R / 2) ** 3 * (X1**2 - E1**2)
r1 = rho.ravel(); z1 = z.ravel()
g_m = (W1 * ex * (R / 2)).ravel()             # measure d^3r/(rA rB) (phi integral separate)
g_v = (W1 * ex * jac).ravel()                 # measure d^3r
res = {"f2": 0.0, "I2": 0.0, "f4": 0.0, "I4": 0.0}
for i in range(len(r1)):
    A_ = r1[i] ** 2 + r1**2 + (z1[i] - z1) ** 2
    B_ = 2 * r1[i] * r1
    S = A_ + B_
    S[S == 0] = 1e-300
    Phi1 = 4 * np.sqrt(S) * ellipe(np.clip(2 * B_ / S, 0, 1))
    Phi3 = (np.clip(A_[:, None] - B_[:, None] * np.cos(ph)[None, :], 0, None) ** 1.5) @ wph
    res["f2"] += g_m[i] * np.dot(g_m, Phi1)
    res["I2"] += g_v[i] * np.dot(g_v, Phi1)
    res["f4"] += g_m[i] * np.dot(g_m, Phi3)
    res["I4"] += g_v[i] * np.dot(g_v, Phi3)
res["f2"] *= R / (16 * np.pi**2) * 2 * np.pi
res["f4"] *= R / (16 * np.pi**2) * 2 * np.pi
res["I2"] *= 2 * np.pi
res["I4"] *= 2 * np.pi

# closed forms
ctx = E.Ctx("Q", a, b)
D = E.downset([(1, 1, 1, 1)])
calc0, f0j, Fs, P2, nfree = E.solve_all(ctx, D, check=True)
calc = M.MCalc(ctx)
RJ = M.RJets(ctx, D, calc)
f2j = M.solve_fn(calc, RJ, D, 1, f0j)
f4j = M.solve_fn(calc, RJ, D, 3, f2j, f0j)
mp.mp.dps = 30
alg = PX.Algebra(a, b)
x = mp.mpf(Rq.numerator) / Rq.denominator
f0n = PX.f0_numeric(alg, x)
Z, O = (0, 0, 0, 0), (1, 1, 1, 1)
cf = {"f2": PX.evaluate(alg, f2j[Z], x, f0n), "I2": 16 * mp.pi**2 / x * PX.evaluate(alg, f2j[O], x, f0n),
      "f4": PX.evaluate(alg, f4j[Z], x, f0n), "I4": 16 * mp.pi**2 / x * PX.evaluate(alg, f4j[O], x, f0n)}
out = {"pair": [str(a), str(b)], "R": str(Rq), "N": N}
for k in ("f2", "I2", "f4", "I4"):
    rel = abs(res[k] - float(cf[k])) / abs(float(cf[k]))
    out[k] = {"quadrature": res[k], "closed_form": mp.nstr(cf[k], 20), "rel_diff": rel}
    print(k, "quadrature %.10g  closed form %s  rel diff %.2e" % (res[k], mp.nstr(cf[k], 12), rel), flush=True)
tg = f"{a}_{b}_{Rq}".replace("/", "-")
json.dump(out, open(f"debug/data/pachucki_m_quadcheck_{tg}.json", "w"), indent=1)
