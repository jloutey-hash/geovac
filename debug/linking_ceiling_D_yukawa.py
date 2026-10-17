"""PL2 step D: is F_4's operator a genuine Calabi-Yau (Sp4) operator, not a symmetric cube or a
reducible one? Step B's operator in w = t^2 has exponents 1/2 (x4) at w = inf. With v = 1/w,
y = v^{1/2} g, x = v/4 it becomes (check 1, symbolic) the MUM operator
  L = th^4 - x(2th+1)^2(10th^2+10th+3) + 9x^2(2th+1)^2(2th+3)^2.
(2) holomorphic coefficients A_n integral; (3) Yukawa coupling K(q) = 1 + sum n_d d^3 q^d/(1-q^d):
a symmetric cube of an order-2 operator gives K = 1 exactly (all n_d = 0); a CY operator gives
nonzero n_d. Exact arithmetic: power series in x with Fraction coefficients, eps-expansion by hand."""
import sympy as sp, json
from fractions import Fraction as Fr
u, th = sp.symbols("u th")
Q = [sp.sympify(q) for q in json.load(open("debug/data/linking_ceiling_B_ode.json"))["4"]["Q"]]
Lw = [sp.expand(q.subs(u, 2*th)) for q in Q]
Lv = {2 - i: sp.expand(c.subs(th, -th)) for i, c in enumerate(Lw)}
Lx = {k: sp.expand(c.subs(th, th + sp.Rational(1, 2)) * 4**k) for k, c in Lv.items()}
lead = Lx[0].as_poly(th).LC()
Lx = {k: sp.expand(c / lead) for k, c in Lx.items()}
target = {0: th**4, 1: -(2*th+1)**2*(10*th**2+10*th+3), 2: 9*(2*th+1)**2*(2*th+3)**2}
ok1 = all(sp.expand(Lx[k] - target[k]) == 0 for k in target)
print("(1) transformed step-B operator equals L:", ok1, flush=True)

N, E = 30, 4                                   # x-order, eps-order
def emul(a, b): return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(E)]
def einv(a):
    r = [Fr(1) / a[0]] + [Fr(0)] * (E - 1)
    for k in range(1, E): r[k] = -sum(a[i]*r[k-i] for i in range(1, k+1)) / a[0]
    return r
def epoly(f, n):  # f polynomial in s, evaluated at s = n + e, as eps-series
    s = sp.Symbol("s"); p = sp.Poly(sp.expand(f(s + n)), s)
    return [Fr(int(p.coeff_monomial(s**k).p), int(p.coeff_monomial(s**k).q)) if k <= p.degree() else Fr(0) for k in range(E)]
P0 = lambda s: s**4
P1 = lambda s: (2*s+1)**2*(10*s**2+10*s+3)
P2 = lambda s: 9*(2*s+1)**2*(2*s+3)**2
A = [[Fr(1)] + [Fr(0)] * (E - 1)]
for n in range(1, N):
    rhs = emul(epoly(P1, n - 1), A[n-1])
    if n >= 2: rhs = [r - q for r, q in zip(rhs, emul(epoly(P2, n - 2), A[n-2]))]
    A.append(emul(rhs, einv(epoly(P0, n))))
f = [[A[n][k] for n in range(N)] for k in range(E)]       # f_k(x) coefficients
print("(2) A_n integral:", all(c.denominator == 1 for c in f[0]), [int(c) for c in f[0][:8]], flush=True)

def mul(a, b): return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(N)]
def inv(a):
    r = [Fr(1)/a[0]] + [Fr(0)]*(N-1)
    for k in range(1, N): r[k] = -sum(a[i]*r[k-i] for i in range(1, k+1)) / a[0]
    return r
def theta(a): return [k * a[k] for k in range(N)]
def exp0(a):  # exp of a series with a[0] = 0
    r = [Fr(1)] + [Fr(0)]*(N-1)
    for k in range(1, N): r[k] = sum(j * a[j] * r[k-j] for j in range(1, k+1)) / k
    return r
i0 = inv(f[0]); g1 = mul(f[1], i0); g2 = mul(f[2], i0)
H = [a - b/2 for a, b in zip(g2, mul(g1, g1))]
dtd = inv([Fr(1) + c if k == 0 else c for k, c in enumerate(theta(g1))])
Dt = lambda p: mul(dtd, theta(p))
K = Dt(Dt(H)); K[0] += 1
# mirror map q = x exp(g1); invert to x(q), substitute into K
qx = [Fr(0)] + exp0(g1)[:N-1]                                # q as series in x
def compose(a, b):  # a(b(q)), b[0] = 0
    r = [Fr(0)]*N; pw = [Fr(1)] + [Fr(0)]*(N-1)
    for k in range(N):
        if a[k]: r = [x + a[k]*y for x, y in zip(r, pw)]
        pw = mul(pw, b)
    return r
xq = [Fr(0), Fr(1)] + [Fr(0)]*(N-2)
for _ in range(N):
    c = compose(qx, xq); xq = [x - (y - (1 if k == 1 else 0)) for k, (x, y) in enumerate(zip(xq, c))]
Kq = compose(K, xq)
M = 12
print("(3) K(q) =", [str(c) for c in Kq[:M]], flush=True)
rem = list(Kq[:M]); nd = {}
for d in range(1, M):
    nd[d] = rem[d] / d**3
    for m in range(2*d, M, d): rem[m] -= nd[d] * d**3
print("    n_d =", [str(nd[d]) for d in range(1, M)], flush=True)
json.dump({"operator_ok": ok1, "A": [str(c) for c in f[0][:12]], "K": [str(c) for c in Kq[:M]],
           "n_d": {d: str(nd[d]) for d in nd}}, open("debug/data/linking_ceiling_D_yukawa.json", "w"), indent=1)
