"""l > 0 exchange integrals whose orbital products are polynomials in the distances
(m = 0 families): same Pachucki-ODE derivation as pachucki_bound_1s.py, with higher
multi-index exponent derivatives.

  z_A = (rA^2 + R^2 - rB^2)/(2R),   z_B = z_A - R,   r^2 = rA^2 or rB^2.
  r1A^i r1B^j r2A^k r2B^l e^{...} / r12  ->  (+-) d_u3^{i+1} d_u2^{j+1} d_w2^{k+1} d_w3^{l+1} f   (sign (-1)^{i+j+k+l})

Cases: 2psA  = (2p_sigma_A 1s_B | 2p_sigma_A 1s_B)
       3d0A  = (3d_0_A 1s_B | 3d_0_A 1s_B)
       2ps2ps = (2p_sigma_A 2p_sigma_B | 2p_sigma_A 2p_sigma_B)
Run: python debug/pachucki_bound_lgt0.py CASE A B   (A, B = orbital exponents, A < B)
"""
from __future__ import annotations

import itertools
import json
import sys
import time

import mpmath as mp
import sympy as sp

sys.path.insert(0, "debug")
import pachucki_bound_core as C  # noqa: E402
from pachucki_bound_core import r, LC  # noqa: E402

case = sys.argv[1]
a = sp.Rational(sys.argv[2])
b = sp.Rational(sys.argv[3])
PMIN, PMAX = -14, 14
t0 = time.time()
rA, rB = sp.symbols("rA rB", positive=True)
zA = (rA**2 + r**2 - rB**2) / (2 * r)
zB = zA - r
if case == "2psA":
    P, N2, hp = zA, a**5 / sp.pi * b**3 / sp.pi, ((2 * a, (2, 1, 0)), (b, (1, 0, 0)))
elif case == "3d0A":
    P, N2, hp = 3 * zA**2 - rA**2, a**7 / (18 * sp.pi) * b**3 / sp.pi, ((3 * a, (3, 2, 0)), (b, (1, 0, 0)))
elif case == "2ps2ps":
    P, N2, hp = zA * zB, a**5 / sp.pi * b**5 / sp.pi, ((2 * a, (2, 1, 0)), (2 * b, (2, 1, 0)))
else:
    raise SystemExit("unknown case")
# integrand polynomial: P(r1) P(r2);  K = N2 (4 pi)^2 / R * sum coeff * (+-) d^alpha f
P1 = sp.Poly(sp.expand(P * r**4), rA, rB)          # clear 1/R powers: multiply by r^4, divide later
terms = []
for (i, j), c1 in P1.terms():
    for (k, l), c2 in P1.terms():
        al = (i + 1, j + 1, k + 1, l + 1)
        sign = (-1) ** (i + j + k + l)
        terms.append((al, sp.expand(sign * c1 * c2)))      # coefficient carries r^8 (two factors of r^4)
need = set()
for al, _ in terms:
    need |= set(C.sub_multi(al))
need = sorted(need, key=lambda x: (sum(x), x))
print(case, "multi-indices needed:", len(need))

F, p2, lnabs = C.pachucki_F_and_p2()
sub = C.theta0(a, b)
lnval = sp.simplify(lnabs.subs(sub))
lnc = (lnval,)
u, w = a, b
Lw = sp.log(abs((w - u) / (w + u)))
Fm1 = (1 / (2 * (u**2 - w**2))) * (
    (sp.exp(2 * r * w) - sp.exp(2 * r * u)) * sp.Ei(-2 * r * (u + w))
    - sp.exp(-2 * r * w) * (sp.Ei(2 * r * (w - u)) - 2 * sp.Ei(r * (w - u)) + Lw)
    + sp.exp(-2 * r * u) * (sp.Ei(2 * r * (u - w)) - 2 * sp.Ei(r * (u - w)) + Lw))
C.FM1.update(C.to_func(sp.expand(Fm1).subs(Lw, -lnval), lnc))
assert C.fequal(C.fdiff(C.FM1), C.to_func(F.subs(sub), lnc))

D = C.multi_derivs(F, p2, need, sub)
Fs = {al: C.to_func(D[al][0], lnc) for al in need}
P2 = {al: D[al][1] for al in need}
print("derivatives done %.0fs" % (time.time() - t0))
keys = set()
for al in need:
    keys |= C.keys_closure(Fs[al])
keys = sorted(k for k in keys if not (k[0] == "E" and k[1] == 0))
ans_keys = keys + [("f",)]
fS = {(0, 0, 0, 0): {("f",): {0: sp.Integer(1)}}}
nfree = 0
for al in need:
    if sum(al) == 0:
        continue
    src = C.fscale(Fs[al], -1)
    for be in C.sub_multi(al):
        U = tuple(x - y for x, y in zip(al, be))
        if sum(U) == 0:
            continue
        cc = C.binom_multi(al, U) * P2[U]
        if cc != 0:
            src = C.fadd(src, C.fscale(fS[be], cc, 1))
    sol, free = C.solve_L0(src, ans_keys, PMIN, PMAX)
    fS[al] = sol
    nfree += len(free)
print("all f_alpha solved in the 11-function span; free ansatz columns:", nfree, " %.0fs" % (time.time() - t0))

Kfun: dict = {}
for al, cf in terms:
    for pw, cc in sp.Poly(cf, r).terms():
        Kfun = C.fadd(Kfun, C.fscale(fS[al], 16 * sp.Rational(N2 * sp.pi**2) * cc, pw[0] - 8 - 1))
closure = sorted(C.keys_closure(Kfun), key=str)
rows, cur = [], Kfun
for i in range(len(closure) + 1):
    rows.append([sum(c * r**p for p, c in cur.get(k, {}).items()) for k in closure])
    cur = C.fdiff(cur)
from sympy.polys.matrices import DomainMatrix  # noqa: E402
rank = DomainMatrix.from_Matrix(sp.Matrix(rows)).convert_to(sp.QQ.frac_field(r, LC)).rank()
print("closure keys:", len(closure), " rank (operator order bound):", rank)

# numeric check against geovac.qfd_core.exchange_hp
mp.mp.dps = 40
lcnum = mp.log(mp.mpf(sp.Rational(a + b, b - a).p) / sp.Rational(a + b, b - a).q)


def ev(A, x, f0):
    s = mp.mpf(0)
    for k, d in A.items():
        if k[0] == "f":
            bk = f0
        elif k[0] == "E":
            bk = mp.exp(mp.mpf(k[1].p) / k[1].q * x)
        else:
            bk = mp.exp(mp.mpf(k[1].p) / k[1].q * x) * mp.ei(mp.mpf(k[2].p) / k[2].q * x)
        for p, c in d.items():
            cc = sp.Poly(c, LC)
            cv = sum(mp.mpf(sp.Rational(co).p) / sp.Rational(co).q * lcnum**m[0]
                     for m, co in zip(cc.monoms(), cc.coeffs()))
            s += cv * x**p * bk
    return s


x = mp.mpf("1.4")
f0 = mp.quad(lambda s: ev(C.FM1, s, None) / s, [x, x + 2, x + 10, mp.inf])
val = ev(Kfun, x, f0)
sys.path.insert(0, ".")
from geovac import qfd_core as Q  # noqa: E402
(ZA, oA), (ZB, oB) = hp
refs = {}
for T in (16, 24):
    with mp.workdps(40):
        vh, per = Q.exchange_hp(sp.Rational(ZA), oA, oB, sp.Rational(ZB), oA, oB, sp.Rational("1.4"),
                                tau_max=T, dps=40)
    refs[T] = vh
print("Pachucki route K(1.4) =", mp.nstr(val, 35))
print("exchange_hp tau<=16   =", mp.nstr(refs[16], 35))
print("exchange_hp tau<=24   =", mp.nstr(refs[24], 35))
out = {"case": case, "pair": [str(a), str(b)], "multi_indices": len(need), "free_ansatz_columns": nfree,
       "closure_keys": [str(k) for k in closure], "rank": int(rank),
       "K_1.4_pachucki": mp.nstr(val, 35), "K_1.4_exchange_hp_tau16": mp.nstr(refs[16], 35),
       "K_1.4_exchange_hp_tau24": mp.nstr(refs[24], 35), "time_s": round(time.time() - t0, 1)}
tg = f"{a}_{b}".replace("/", "-")
json.dump(out, open(f"debug/data/pachucki_bound_{case}_{tg}.json", "w"), indent=1)
print("done %.0fs" % (time.time() - t0))
