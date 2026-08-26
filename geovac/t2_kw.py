"""The (KW) representation of the Paper 59 collinear observable T2.

    T2 = (8/pi) int_0^inf dk int_0^1 dw cos(kw) R(k,w)^2,
    R(k,w) = int_0^1 cos(kw(s-1/2)) P(s,k) ds,
    P(x,k) = c e^{-D} (D^-3 + 3 D^-4 + 3 D^-5),  c = x(1-x),  D = sqrt(c k^2+1)

(Paper 59 eq:kw).  The exact factorization dissolves the (s,t)-frame outer wall;
the truncation tail is analytic (Watson expansion, TAIL below), so certified
digits are a matter of quadrature depth.  Promoted from debug/beta2_t2_kw_core.py
+ beta2_t2_tail.py + _fastgl.py (v4.106.x QA remediation) so the certification is
executable from the tracked tree.  API:  int_F(a, b, npan, nnode, Ns, p)
for the k-integral on [a,b], plus TAIL(K, F_buckets(MX)) for the analytic
remainder; T2 = (8/pi) * (int_F(0,K,...) + TAIL(K, F_buckets(MX))).
Production certification record: debug/beta2_track_a_findings.md (66 digits,
six parameter-disjoint runs; u1/u2 extend cross-validation).  Backing test:
tests/test_paper59_t2_value.py (mpmath witness leg).
"""
from __future__ import annotations

import mpmath as mp

# ------------------------- Gauss-Legendre nodes (cached)
import numpy as np

def _legendre_and_deriv(N, x):
    # P_N(x), P_N'(x) via recurrence (all mpf)
    p0 = mp.mpf(1); p1 = x
    if N == 0: return p0, mp.mpf(0)
    for k in range(2, N+1):
        p0, p1 = p1, ((2*k-1)*x*p1 - (k-1)*p0)/k
    pN = p1
    dP = N*(x*pN - p0)/(x*x - 1)
    return pN, dP

_CACHE = {}
def fast_gl(N):
    key = (N, mp.mp.dps)
    if key in _CACHE: return _CACHE[key]
    x0, w0 = np.polynomial.legendre.leggauss(N)   # double precision seed
    xs = []; ws = []
    for xd in x0:
        x = mp.mpf(float(xd))
        for _ in range(5):
            pN, dP = _legendre_and_deriv(N, x)
            dx = pN/dP; x -= dx
            if abs(dx) < mp.mpf(10)**(-(mp.mp.dps+6)): break
        pN, dP = _legendre_and_deriv(N, x)
        xs.append(x); ws.append(2/((1-x*x)*dP*dP))
    _CACHE[key] = (xs, ws)
    return xs, ws

# ------------------------- the (KW) integrand
# ---------------------------------------------------------------- basic pieces
def P(x, k):
    c = x * (1 - x)
    D = mp.sqrt(c * k * k + 1)
    return c * mp.e ** (-D) * (1 / D ** 3 + 3 / D ** 4 + 3 / D ** 5)


_SNODES = {}
def s_nodes(Ns, p):
    """GL nodes/weights for  2 int_0^{1/2} f(s) ds  under s = (1/2) y^p, y in [0,1]."""
    key = (Ns, p, mp.mp.dps)
    if key in _SNODES:
        return _SNODES[key]
    xs, ws = fast_gl(Ns)
    S, W = [], []
    for xg, wg in zip(xs, ws):
        y = (xg + 1) / 2
        s = y ** p / 2
        jac = p * y ** (p - 1) / 2          # ds/dy
        S.append(s)
        W.append(2 * (wg / 2) * jac)        # the outer factor 2 of R folded in
    _SNODES[key] = (S, W)
    return S, W


def R_setup(k, Ns, p):
    """Precompute (1/2 - s_j) and W_j*P(s_j,k) for fixed k."""
    S, W = s_nodes(Ns, p)
    a = []; c = []
    for s, w in zip(S, W):
        a.append((mp.mpf(1) / 2 - s) * k)
        c.append(w * P(s, k))
    return a, c


def R_of_w(a, c, w):
    tot = mp.mpf(0)
    for aj, cj in zip(a, c):
        tot += cj * mp.cos(aj * w)
    return tot


def F_of_k(k, Ns, Nw, p):
    """F(k) = int_0^1 cos(kw) R(k,w)^2 dw."""
    a, c = R_setup(k, Ns, p)
    xs, ws = fast_gl(Nw)
    tot = mp.mpf(0)
    for xg, wg in zip(xs, ws):
        w = (xg + 1) / 2
        Rv = R_of_w(a, c, w)
        tot += (wg / 2) * mp.cos(k * w) * Rv * Rv
    return tot


def nw_for(k, base=40, slope=mp.mpf('0.62')):
    return int(slope * float(k)) + base


def int_F(a, b, npan, nnode, Ns, p=6, nw=None, prog=None):
    """Panelled Gauss-Legendre for int_a^b F(k) dk."""
    a = mp.mpf(a); b = mp.mpf(b)
    xs, ws = fast_gl(nnode)
    h = (b - a) / npan
    tot = mp.mpf(0)
    for ip in range(npan):
        lo = a + ip * h
        for xg, wg in zip(xs, ws):
            k = lo + h * (xg + 1) / 2
            NW = nw(k) if nw else nw_for(k)
            NS = Ns(k) if callable(Ns) else Ns
            tot += (h * wg / 2) * F_of_k(k, NS, NW, p)
        if prog: prog(ip, npan, tot)
    return tot

# ------------------------- the analytic Watson tail
# --------------------------------------------------------------- moments mu_n
def mu_moments(N):
    """mu_n = int_1^inf e^{-v} (v^2-1)^{n+1} (v^-2 + 3 v^-3 + 3 v^-4) dv, n = 0..N."""
    old = mp.mp.dps
    mp.mp.dps = old + 60
    G = {}
    def gam(m):                     # int_1^inf e^{-v} v^m dv = Gamma(m+1, 1)
        if m not in G:
            G[m] = mp.gammainc(m + 1, 1)
        return G[m]
    out = []
    for n in range(N + 1):
        n1 = n + 1
        tot = mp.mpf(0)
        for i in range(n1 + 1):
            cb = mp.binomial(n1, i) * (-1) ** (n1 - i)
            tot += cb * (gam(2 * i - 2) + 3 * gam(2 * i - 3) + 3 * gam(2 * i - 4))
        out.append(tot)
    mp.mp.dps = old
    return [+v for v in out]


# ------------------------------------------------------- truncated 3-var series
# key (mx, mw, my) -> mpf coefficient ;  truncate mx <= MX
def smul(a, b, MX):
    out = {}
    for (i1, j1, l1), c1 in a.items():
        if c1 == 0: continue
        for (i2, j2, l2), c2 in b.items():
            i = i1 + i2
            if i > MX: continue
            key = (i, j1 + j2, l1 + l2)
            out[key] = out.get(key, mp.mpf(0)) + c1 * c2
    return out


def sadd(a, b, sc=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, mp.mpf(0)) + sc * v
    return out


def build_AB(MX):
    """A(w,x), B(w,x) as dicts {(mx,mw): mpf} with y contracted against mu_n."""
    MU = mu_moments(MX + 2)
    # Catalan numbers
    Cat = [mp.mpf(1)]
    for m in range(1, MX + 3):
        Cat.append(Cat[-1] * (2 * (2 * m - 1)) / (m + 1))
    # phi = w * x * g = sum_m Cat_m w y^{m+1} x^{2m+1}
    phi = {}
    m = 0
    while 2 * m + 1 <= MX:
        phi[(2 * m + 1, 1, m + 1)] = Cat[m]
        m += 1
    # h = sum_m binom(2m,m) y^m x^{2m}
    h = {}
    m = 0
    while 2 * m <= MX:
        h[(2 * m, 0, m)] = mp.binomial(2 * m, m)
        m += 1
    # cos(phi), sin(phi)
    one = {(0, 0, 0): mp.mpf(1)}
    cosp = dict(one); sinp = dict(phi)
    pw = dict(phi)                       # phi^1
    r = 1
    fact = mp.mpf(1)
    while True:
        pw = smul(pw, phi, MX)           # phi^{r+1}
        r += 1
        if not pw: break
        fact = mp.factorial(r)
        term = {k: v / fact for k, v in pw.items()}
        if r % 2 == 0:
            cosp = sadd(cosp, term, (-1) ** (r // 2))
        else:
            sinp = sadd(sinp, term, (-1) ** ((r - 1) // 2))
        if r > 2 * MX + 4: break
    A3 = smul(h, cosp, MX)
    B3 = smul(h, sinp, MX)
    def contract(S):
        out = {}
        for (i, j, l), c in S.items():
            key = (i, j)
            out[key] = out.get(key, mp.mpf(0)) + c * MU[l]
        return out
    return contract(A3), contract(B3)


def pmul(a, b, MX):
    """{(mx,mw)} product truncated at mx <= MX."""
    out = {}
    for (i1, j1), c1 in a.items():
        for (i2, j2), c2 in b.items():
            i = i1 + i2
            if i > MX: continue
            key = (i, j1 + j2)
            out[key] = out.get(key, mp.mpf(0)) + c1 * c2
    return out


def R_asym(k, w, A, B):
    """R(k,w) from the Watson series (validation helper)."""
    x = 1 / mp.mpf(k)
    sa = mp.mpf(0); sb = mp.mpf(0)
    for (i, j), c in A.items(): sa += c * x ** i * w ** j
    for (i, j), c in B.items(): sb += c * x ** i * w ** j
    return 4 * x ** 4 * (mp.cos(k * w / 2) * sa + mp.sin(k * w / 2) * sb)


# ------------------------------------------------- w-integrals -> (k-power, trig) buckets
# kinds: '1', 'cos1', 'sin1', 'cos2', 'sin2'
_RE_E = {0: ('cos', 1), 1: ('sin', 1), 2: ('cos', -1), 3: ('sin', -1)}
_IM_E = {0: ('sin', 1), 1: ('cos', -1), 2: ('sin', -1), 3: ('cos', 1)}
_RE_C = {0: 1, 1: 0, 2: -1, 3: 0}
_IM_C = {0: 0, 1: -1, 2: 0, 3: 1}


def w_int_buckets(n, kind, coef, acc):
    """Add  coef * int_0^1 w^n {1|cos(kw)|cos(2kw)|sin(2kw)} dw  into acc[(p, trig)] (k^-p * trig)."""
    if kind == 'dc':
        acc[(0, '1')] = acc.get((0, '1'), mp.mpf(0)) + coef / (n + 1)
        return
    a = 1 if kind == 'c1' else 2
    tag = '1' if a == 1 else '2'
    take_re = (kind != 's2')
    fact = mp.mpf(1)                       # n!/(n-j)!  for j=0
    for j in range(n + 1):
        if j > 0:
            fact *= (n - j + 1)
        m = (j + 1) % 4
        sgn = (-1) ** j
        base = coef * sgn * fact / mp.mpf(a) ** (j + 1)
        if take_re:
            trig, s = _RE_E[m]
        else:
            trig, s = _IM_E[m]
        key = (j + 1, trig + tag)
        acc[key] = acc.get(key, mp.mpf(0)) + base * s
    # boundary term  -(-1)^n n! (i a k)^{-(n+1)}
    m = (n + 1) % 4
    c = _RE_C[m] if take_re else _IM_C[m]
    if c:
        base = -coef * (-1) ** n * fact * (n) and None   # placeholder (fact already = n!)
    fact_n = mp.factorial(n)
    if c:
        acc[(n + 1, '1')] = acc.get((n + 1, '1'), mp.mpf(0)) \
            - coef * (-1) ** n * fact_n * c / mp.mpf(a) ** (n + 1)


def F_buckets(MX):
    """F(k) = sum_{(p,trig)} C * k^-p * trig   (asymptotic, order x^MX)."""
    A, B = build_AB(MX)
    A2 = pmul(A, A, MX); B2 = pmul(B, B, MX); AB = pmul(A, B, MX)
    U = dict(A2)
    for k_, v in B2.items(): U[k_] = U.get(k_, mp.mpf(0)) + v
    V = dict(A2)
    for k_, v in B2.items(): V[k_] = V.get(k_, mp.mpf(0)) - v
    acc = {}
    half = mp.mpf(1) / 2; quart = mp.mpf(1) / 4
    for src, fac, kind in ((U, half, 'c1'), (V, quart, 'dc'), (V, quart, 'c2'), (AB, half, 's2')):
        for (mx, mw), cv in src.items():
            if cv == 0: continue
            w_int_buckets(mw, kind, 16 * fac * cv, acc_shift(acc, 8 + mx))
    return acc


class acc_shift:
    """Wrapper adding a fixed shift to the k-power index when writing into acc."""
    def __init__(self, acc, shift):
        self.acc = acc; self.shift = shift
    def get(self, key, default):
        return self.acc.get((key[0] + self.shift, key[1]), default)
    def __setitem__(self, key, val):
        self.acc[(key[0] + self.shift, key[1])] = val


def F_asym(k, acc):
    k = mp.mpf(k); tot = mp.mpf(0)
    c1 = mp.cos(k); s1 = mp.sin(k); c2 = mp.cos(2 * k); s2 = mp.sin(2 * k)
    TR = {'1': mp.mpf(1), 'cos1': c1, 'sin1': s1, 'cos2': c2, 'sin2': s2}
    for (p, trig), c in acc.items():
        tot += c * TR[trig] * k ** (-p)
    return tot


_EP = {}
def _E(P, a, K):
    """int_K^inf k^-P e^{i a k} dk."""
    key = (P, a, str(K), mp.mp.dps)
    if key not in _EP:
        z = mp.mpc(0, -a)
        _EP[key] = z ** (P - 1) * mp.gammainc(1 - P, z * K)
    return _EP[key]


def TAIL(K, acc):
    K = mp.mpf(K); tot = mp.mpf(0)
    for (p, trig), c in acc.items():
        if trig == '1':
            tot += c * K ** (1 - p) / (p - 1)
        else:
            a = 1 if trig.endswith('1') else 2
            E = _E(p, a, K)
            tot += c * (E.real if trig.startswith('cos') else E.imag)
    return tot
