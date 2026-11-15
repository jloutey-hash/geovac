"""Record AE2 control (written by the fresh reviewer of the AE2 draft, read and run by the PM, 2026-11-15; reproduces
the PM's closed-form AE amplitudes at E1 to 12 digits).  Independent local expansion of the non-analytic part of F(a) = int V_A(u-a) rho(u) du about C,
to t^{11/2} (t = 1-x), all terms. rho terms w^k R^m (w=u_z, R=|u|, m odd); V terms v_n s^n (n odd)."""
import sympy as sp
import mpmath as mp
from fractions import Fraction as Fr
from collections import defaultdict
mp.mp.dps = 40

def c(q):
    return 2 ** (q + 3) * mp.pi ** 1.5 * mp.gamma(mp.mpf(q + 3) / 2) / mp.gamma(-mp.mpf(q) / 2)
def K(m, n):
    return c(m) * c(n) / c(m + n + 3)

def dz(poly):  # poly: {(a,b): coeff} for z^a r^b ; d/dz
    out = defaultdict(lambda: 0)
    for (a, b), v in poly.items():
        if a: out[(a - 1, b)] += v * a
        if b: out[(a + 1, b - 2)] += v * b
    return dict(out)
def dzj(j, q):
    p = {(0, q): 1}
    for _ in range(j): p = dz(p)
    return p
def decompose(k, m):
    """z^k r^m = sum_j b_j dz^j r^{q_j}; returns list (j, q, b)."""
    target = {(k, m): Fr(1)}
    res = []
    while target:
        target = {key: v for key, v in target.items() if v != 0}
        if not target: break
        a, b = max(target, key=lambda x: x[0])
        coef = target[(a, b)]
        q = b + 2 * a
        p = dzj(a, q)
        top = p[(a, b)]
        f = Fr(coef) / Fr(top)
        res.append((a, q, f))
        for key, v in p.items():
            target[key] = target.get(key, 0) - f * v
        target = {key: v for key, v in target.items() if v != 0}
    return res

def expansion(kind, D, zB, zC, zA, maxdeg=5, tmax=Fr(11, 2)):
    w, R, e = sp.symbols('w R e', positive=True)
    rB = sp.sqrt(D ** 2 + 2 * D * w * e + R ** 2 * e ** 2)
    taylor = sp.series(sp.exp(-zB * rB), e, 0, maxdeg + 1).removeO()
    oddC = sum((-zC) ** kk * (R * e) ** kk / sp.factorial(kk) for kk in range(1, maxdeg + 1, 2))
    ang = {'S': 1, 'PC': w * e, 'PD': (3 * w ** 2 - R ** 2) * e ** 2 / 2}[kind]
    rho = sp.expand(taylor * oddC * ang)
    rho = sum(rho.coeff(e, d) for d in range(0, maxdeg + 1))
    rho = sp.expand(rho)
    s = sp.symbols('s')
    V = sp.series(1 / s - sp.exp(-2 * zA * s) * (zA + 1 / s), s, 0, 9).removeO()
    vn = {n: V.coeff(s, n) for n in (1, 3, 5, 7)}
    out = defaultdict(lambda: mp.mpf(0))
    for term in sp.Add.make_args(rho):
        cw = sp.Poly(term, w, R)
        for (k, m), cf in zip(cw.monoms(), cw.coeffs()):
            for n, v in vn.items():
                if v == 0: continue
                for (j, q, b) in decompose(k, m):
                    Q = q + n + 3
                    kk = K(q, n)
                    for (a, bb), vv in dzj(j, Q).items():
                        pw = Fr(a) + Fr(bb, 2)
                        if pw > tmax: continue
                        val = mp.mpf(sp.N(cf * v, 50)) * float(b) * kk * vv * (-D) ** a * (2 * D * D) ** (mp.mpf(bb) / 2)
                        out[pw] += val
    return dict(out)

if __name__ == '__main__':
    import sys
    cases = {'E1': ('2', '13/10', '1', '9/10'), 'E2': ('12/5', '1', '6/5', '4/5'),
             'E3': ('3', '3', '1', '9/10'), 'E4': ('5/2', '7/10', '13/10', '11/10')}
    for cs in sys.argv[1:]:
        D, zB, zC, zA = (sp.Rational(x) for x in cases[cs])
        D_ = mp.mpf(sp.N(D, 40))
        NB = mp.sqrt(mp.mpf(sp.N(zB, 40)) ** 3 / mp.pi); z = mp.mpf(sp.N(zC, 40))
        N = {'S': mp.sqrt(z ** 3 / mp.pi), 'PC': mp.sqrt(z ** 5 / mp.pi),
             'PD': mp.sqrt(5 * (2 * z) ** 7 / (4 * mp.pi * 720))}
        f0 = 1  # exp(-zB D) already inside taylor
        for kind in ('S', 'PC', 'PD'):
            ex = expansion(kind, D, zB, zC, zA)
            print(cs, kind, {str(k): mp.nstr(v * NB * N[kind], 12) for k, v in sorted(ex.items()) if abs(v) > 1e-30})
