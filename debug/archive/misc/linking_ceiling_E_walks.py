"""Pattern check: the MUM holomorphic solution of F_n's operator (v = 1/t^2, y = v^{1/2} g) has
coefficients C(2m,m)^2 * W_{n-1}(2m) / scale^m, W_j(2m) = sum over j-part compositions of
(m!/(a1!...aj!))^2 (even moments of the j-step planar random walk): W_1 = 1, W_2 = C(2m,m),
W_3 = sum C(m,k)^2 C(2k,k), W_4 = Domb numbers."""
import sympy as sp, json
from fractions import Fraction as Fr
from itertools import product
u, th = sp.symbols("u th")
src = open("debug/linking_ceiling_B_ode.py").read().split("import os")[0]
ns = {}; exec(src, ns)
def W(j, m):
    tot = 0
    def rec(parts, left, k):
        nonlocal tot
        if k == j - 1:
            ps = parts + [left]; c = sp.factorial(m)
            for p in ps: c //= sp.factorial(p)
            tot += c**2; return
        for a in range(left + 1): rec(parts + [a], left - a, k + 1)
    rec([], m, 0); return tot
for n in (2, 3, 4, 5):
    Q = ns["recurrence_to_t_ode"](ns["sym_power_ode"](n - 1))
    K = len(Q) - 1
    R = {K - i: sp.Poly(sp.expand(Q[i].subs(u, -2*th).subs(th, th + sp.Rational(1, 2))), th) for i in range(K + 1)}
    G = [Fr(1)]
    for m in range(1, 9):
        s = sum(Fr(int(sp.Rational(R[k].as_expr().subs(th, m - k)).p), int(sp.Rational(R[k].as_expr().subs(th, m - k)).q)) * G[m - k]
                for k in range(1, K + 1) if m - k >= 0)
        r0 = sp.Rational(R[0].as_expr().subs(th, m))
        G.append(-s / Fr(int(r0.p), int(r0.q)))
    target = [sp.binomial(2*m, m)**2 * W(n - 1, m) for m in range(9)]
    ratio = [G[m] / Fr(int(target[m])) for m in range(9)]
    scale = ratio[1]
    ok = all(ratio[m] == scale**m for m in range(9))
    print(f"n={n}: G = C(2m,m)^2 W_{n-1}(2m) * ({scale})^m for m<9: {ok}")
