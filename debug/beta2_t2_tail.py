"""ANALYTIC large-k tail for the (KW) frame:   TAIL(K) = int_K^inf F(k) dk.

Derivation (exact, then asymptotic in x = 1/k):
  R(k,w) = 2 int_0^{1/2} cos(k w (1/2 - s)) P(s,k) ds.
  Substitute the EXACT map  c = s(1-s) = (v^2-1) x^2   (v = D = sqrt(c k^2+1)), so
      s = sum_m Cat_m c^{m+1},   ds = 2 v x^2 dv / sqrt(1-4c),   P = c e^{-v}(v^-3+3v^-4+3v^-5).
  =>  P ds = 2 x^4 Psi(v) h dv,   Psi(v) = e^{-v} (v^2-1)(v^-2 + 3 v^-3 + 3 v^-4),
      h = (1-4 y x^2)^{-1/2},  y := v^2-1,   and  k w (1/2 - s) = k w/2 - w x g,
      g = sum_m Cat_m y^{m+1} x^{2m}   (= k^2 s).
  => R = 4 x^4 [ cos(kw/2) A(w,x) + sin(kw/2) B(w,x) ] + O(e^{-k/2}),
     A = <h cos(w x g)>,  B = <h sin(w x g)>,  <.> = int_1^inf Psi(v) (.) dv  (y^n -> mu_n).
  A is even in x, B odd.  mu_n = int_1^inf Psi(v) y^n dv = Gamma-values at 1 (grows like (2n)!;
  the series is Gevrey with optimal-truncation floor e^{-k/2} -- exactly the neglected v>V region).

  F(k) = int_0^1 cos(kw) R^2 dw
       = 16 x^8 int_0^1 dw [ (A^2+B^2)/2 cos(kw) + (A^2-B^2)/4 + (A^2-B^2)/4 cos(2kw)
                              + A B /2 sin(2kw) ]
  Every w-integral of w^n against {1, cos kw, cos 2kw, sin 2kw} is ELEMENTARY, giving a finite
  sum of  c * k^{-p} * {1, cos k, sin k, cos 2k, sin 2k}.  Then
      int_K^inf k^{-p} dk           = K^{1-p}/(p-1)
      int_K^inf k^{-p} e^{i a k} dk = (-i a)^{p-1} Gamma(1-p, -i a K).
"""
from __future__ import annotations
import mpmath as mp


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
