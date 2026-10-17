"""Step D follow-up: (a) local exponents of L at the finite singular points (apparent vs conifold);
(b) Yukawa / instanton numbers in the integral variable z = x/4 (4^n A_n integral)."""
import sympy as sp
from fractions import Fraction as Fr
from sympy.functions.combinatorial.numbers import stirling
x, th = sp.symbols("x th")
L = th**4 - x*(2*th+1)**2*(10*th**2+10*th+3) + 9*x**2*(2*th+1)**2*(2*th+3)**2
a = [0]*5
for (k, j), cf in sp.Poly(sp.expand(L), x, th).terms():
    for m in range(j+1): a[m] += cf * x**k * stirling(j, m, kind=2) * x**m
a = [sp.expand(c) for c in a]
for x0 in [sp.Rational(1, 36), sp.Rational(1, 4)]:
    rho = 3 - a[3].subs(x, x0) / sp.diff(a[4], x).subs(x, x0)
    print(f"x0 = {x0}: exponents 0, 1, 2, {rho}")
# (b) integral variable z = x/4
P0 = lambda s: s**4
P1 = lambda s: 4*(2*s+1)**2*(10*s**2+10*s+3)
P2 = lambda s: 144*(2*s+1)**2*(2*s+3)**2
N, E = 30, 4
def emul(p, q): return [sum(p[i]*q[k-i] for i in range(k+1)) for k in range(E)]
def einv(p):
    r_ = [Fr(1)/p[0]] + [Fr(0)]*(E-1)
    for k in range(1, E): r_[k] = -sum(p[i]*r_[k-i] for i in range(1, k+1)) / p[0]
    return r_
def ep(f, n):
    s = sp.Symbol("s"); P = sp.Poly(sp.expand(f(s+n)), s)
    return [Fr(int(P.coeff_monomial(s**k))) for k in range(E)]
A = [[Fr(1)] + [Fr(0)]*(E-1)]
for n in range(1, N):
    rhs = emul(ep(P1, n-1), A[n-1])
    if n >= 2: rhs = [p - q for p, q in zip(rhs, emul(ep(P2, n-2), A[n-2]))]
    A.append(emul(rhs, einv(ep(P0, n))))
f = [[A[n][k] for n in range(N)] for k in range(E)]
print("B_n = 4^n A_n integral:", all(c.denominator == 1 for c in f[0]), [int(c) for c in f[0][:7]])
def mul(p, q): return [sum(p[i]*q[k-i] for i in range(k+1)) for k in range(N)]
def inv(p):
    r_ = [Fr(1)/p[0]] + [Fr(0)]*(N-1)
    for k in range(1, N): r_[k] = -sum(p[i]*r_[k-i] for i in range(1, k+1)) / p[0]
    return r_
def theta(p): return [k*p[k] for k in range(N)]
def exp0(p):
    r_ = [Fr(1)] + [Fr(0)]*(N-1)
    for k in range(1, N): r_[k] = sum(j*p[j]*r_[k-j] for j in range(1, k+1)) / k
    return r_
def compose(p, q):
    r_ = [Fr(0)]*N; pw = [Fr(1)] + [Fr(0)]*(N-1)
    for k in range(N):
        if p[k]: r_ = [u + p[k]*v for u, v in zip(r_, pw)]
        pw = mul(pw, q)
    return r_
i0 = inv(f[0]); g1 = mul(f[1], i0); g2 = mul(f[2], i0)
H = [p - q/2 for p, q in zip(g2, mul(g1, g1))]
dtd = inv([Fr(1) + c if k == 0 else c for k, c in enumerate(theta(g1))])
Dt = lambda p: mul(dtd, theta(p))
K = Dt(Dt(H)); K[0] += 1
qz = [Fr(0)] + exp0(g1)[:N-1]
print("mirror map q(z) integral:", all(c.denominator == 1 for c in qz[:N]))
zq = [Fr(0), Fr(1)] + [Fr(0)]*(N-2)
for _ in range(N):
    c = compose(qz, zq); zq = [u - (v - (1 if k == 1 else 0)) for k, (u, v) in enumerate(zip(zq, c))]
Kq = compose(K, zq)
M = 10
print("K(q) =", [str(c) for c in Kq[:M]])
rem = list(Kq[:M]); nd = {}
for d in range(1, M):
    nd[d] = rem[d] / d**3
    for m in range(2*d, M, d): rem[m] -= nd[d]*d**3
print("n_d =", [str(nd[d]) for d in range(1, M)])
