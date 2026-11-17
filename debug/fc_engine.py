"""High-precision engine for the Legendre coefficients F_L of the four-centre 1s Slater integral (XY|ZW) on a
one-parameter slice (engine for a record to be pre-registered by the PM; this file computes no scored quantity).

Geometry: X at the origin, Y = D1 z-hat; Z = a n-hat, W = b n-hat on one ray from X (0 < a <= b); x = n-hat . z-hat.
Densities are UNNORMALISED Slater products:
    rho1(r) = e^{-zX |r - X|} e^{-zY |r - Y|},    rho2(r) = e^{-zZ |r - Z|} e^{-zW |r - W|},
    I(x) = int int rho1(r1) rho2(r2) / |r1 - r2| = sum_L F_L P_L(x).
Multipoles are LEGENDRE COEFFICIENTS about X: rho1(r) = sum_L rho1_L(r) P_L(r-hat . z-hat), likewise rho2 about n-hat,
so that (addition theorem, both angular integrals)
    F_L = (4 pi / (2L+1))^2 int int rho1_L(r1) rho2_L(r2) r<^L / r>^(L+1) r1^2 r2^2 dr1 dr2.
(With the solid-angle normalisation rho~_L = (4 pi/(2L+1)) rho_L the prefactor is the brief's single 4 pi/(2L+1).)

rho1_L is exact (na_exact's G_L, i_0 = sinh x / x, k_0 = e^{-x}/x):
    rho1_L(r) = -(2L+1) e^{-zX r} [ (2L+1) i_L(x) k_L(y) + x i_{L+1}(x) k_L(y) - y i_L(x) k_{L+1}(y) ],
    x = zY min(r, D1), y = zY max(r, D1);  i_L by downward recurrence from arb Bessel I, k_L by upward recurrence.
rho2_L is numeric: rho2_L(r) = (2L+1)/2 int_0^pi rho2(r, theta) P_L(cos theta) sin theta d theta, on ONE theta grid
shared by every r (uniform Gauss panels for the P_L oscillation, plus panels graded geometrically toward theta = 0,
where the Z and W cusps put complex singularities at theta ~ i |r - a| / sqrt(r a)), so the whole stage is one
arb_mat product P^T V.  |r - Z|^2 = (r - a)^2 + 4 r a sin^2(theta/2) (no cancellation near theta = 0).
The angular stage cancels about L log10(1/rho) digits at large L (rho2_L is exponentially small away from a, b), which
is why it runs at ~100 digits; the radial stage has no cancellation (every rho_L >= 0: the densities are absolutely
monotone in the cosine, so all Legendre coefficients are positive).

Radial stage: r grid = Gauss panels with breakpoints 0, r0, D1, a, b, r_hi and geometric ratio <= q; the inner r1
integral is carried node to node,
    A^(r) = int_0^r rho1_L (s/r)^{L+1} s ds,  B^(r) = int_r^inf rho1_L (r/s)^L s ds,  psi_L = A^ + B^,
(all terms positive, scaled recursions), and F_L = (4 pi/(2L+1))^2 sum_j W_j r_j^2 rho2_L(r_j) psi_L(r_j).
Truncation at r_hi: the densities there are ~ e^{-2 r_hi}, below 1e-30 of F_0 and far below F_L for L >= 1.

Usage:  coefficients(case, Lmax, level) -> (list of arb F_L, info dict);  CASES, LEVELS, landau(case).
"""
from __future__ import annotations

import time
from fractions import Fraction

import flint
from flint import arb, arb_mat, ctx

CASES = {
    'G1': dict(a='1', b='9/4', D1='8/5', zX='13/10', zY='1', zZ='6/5', zW='1'),
    'G2': dict(a='1', b='4', D1='6/5', zX='13/10', zY='1', zZ='6/5', zW='7/10'),
    # K4 reduction: AB case A1 mapped (X = B, Y = C, Z = W = A at D2, zZ + zW = mu = 2 zeta_A = 9/5)
    'K4A1': dict(a='13/5', b='13/5', D1='9/5', zX='13/10', zY='1', zZ='1', zW='4/5'),
}

LEVELS = {
    'A': dict(prec=340, th_w=8, n_u=100, n_g=64, q=Fraction(5, 4), wmax=2, n_r=40, m_sub=12, r_hi=40, r0='1/20'),
    'B': dict(prec=440, th_w=10, n_u=112, n_g=80, q=Fraction(6, 5), wmax=Fraction(3, 2), n_r=48, m_sub=16,
              r_hi=44, r0='1/25'),
}


def params(case):
    if isinstance(case, dict):
        return {k: Fraction(v) for k, v in case.items()}
    return {k: Fraction(v) for k, v in CASES[case].items()}


def landau(case):
    """x_Z*, x_W*, rho_Z, rho_W (exact Fractions) and the dominant point (larger rho)."""
    p = params(case)
    D1, a, b = p['D1'], p['a'], p['b']
    xZ = (D1 ** 2 + a ** 2) / (2 * D1 * a)
    xW = (D1 ** 2 + b ** 2) / (2 * D1 * b)
    rZ = min(D1, a) / max(D1, a)
    rW = min(D1, b) / max(D1, b)
    return dict(xZ=xZ, xW=xW, rhoZ=rZ, rhoW=rW, dominant='Z' if rZ > rW else ('W' if rW > rZ else 'tie'))


def A(q):
    q = Fraction(q)
    return arb(q.numerator) / q.denominator


_GL = {}


def gl(n):
    key = (n, ctx.prec)
    if key not in _GL:
        out = []
        for k in range(n):
            x, w = arb.legendre_p_root(n, k, weight=True)
            out.append((arb(x.mid()), arb(w.mid())))
        _GL[key] = out
    return _GL[key]


def panel(lo, hi, n):
    lo, hi = (v if isinstance(v, arb) else A(v) for v in (lo, hi))
    m, h = (lo + hi) / 2, (hi - lo) / 2
    return [(m + h * x, h * w) for x, w in gl(n)]


# ---------------------------------------------------------------- radial grid
def r_breaks(p, lev):
    q, wmax, r_hi = lev['q'], Fraction(lev['wmax']), Fraction(lev['r_hi'])
    keys = sorted({Fraction(lev['r0']), p['D1'], p['a'], p['b'], r_hi})
    br = [Fraction(0), keys[0]]
    for lo, hi in zip(keys[:-1], keys[1:]):
        cur = lo
        while True:
            nxt = min(cur * q, cur + wmax)
            if nxt >= hi * Fraction(999, 1000):
                br.append(hi)
                break
            br.append(nxt)
            cur = nxt
    return br


def r_grid(p, lev):
    br = r_breaks(p, lev)
    nodes = []
    for lo, hi in zip(br[:-1], br[1:]):
        nodes.extend(panel(lo, hi, lev['n_r']))
    return nodes, br


# ---------------------------------------------------------------- rho1_L (exact)
class Rho1:
    def __init__(self, p, Lmax):
        self.p = p
        self.Lmax = Lmax
        self.zX, self.zY, self.D1 = A(p['zX']), A(p['zY']), A(p['D1'])
        c = self.zY * self.D1
        self.i_c = self.i_all(c)     # i_L(zY D1), used for r > D1
        self.k_c = self.k_all(c)     # k_L(zY D1), used for r < D1

    def i_all(self, x):
        top = self.Lmax + 2
        pre = (arb.pi() / (2 * x)).sqrt()
        i = [None] * (top + 1)
        i[top] = pre * x.bessel_i(arb(top) + arb(1) / 2)
        i[top - 1] = pre * x.bessel_i(arb(top - 1) + arb(1) / 2)
        ix = 1 / x
        for L in range(top - 1, 0, -1):
            i[L - 1] = i[L + 1] + (2 * L + 1) * ix * i[L]
        return i

    def k_all(self, y):
        top = self.Lmax + 2
        e = (-y).exp()
        iy = 1 / y
        k = [e * iy, e * (iy + iy * iy)]
        for L in range(1, top):
            k.append(k[L - 1] + (2 * L + 1) * iy * k[L])
        return k

    def __call__(self, r, inner):
        """rho1_L(r), L = 0..Lmax; inner = (r < D1)."""
        if inner:
            x, y = self.zY * r, self.zY * self.D1
            i, k = self.i_all(x), self.k_c
        else:
            x, y = self.zY * self.D1, self.zY * r
            i, k = self.i_c, self.k_all(y)
        e = (-self.zX * r).exp()
        out = []
        for L in range(self.Lmax + 1):
            c = 2 * L + 1
            out.append(-c * e * ((c * i[L] + x * i[L + 1]) * k[L] - y * i[L] * k[L + 1]))
        return out


# ---------------------------------------------------------------- rho2_L (angular quadrature)
def theta_grid(p, r_nodes, lev):
    pi = arb.pi()
    w = pi / lev['th_w']
    pts = []
    for j in range(1, lev['th_w']):
        pts.extend(panel(w * j, w * (j + 1), lev['n_u']))
    a, b = float(p['a']), float(p['b'])
    eps = min(min(abs(float(r) - a) / (float(r) * a) ** 0.5, abs(float(r) - b) / (float(r) * b) ** 0.5)
              for r, _ in r_nodes)
    th_min = max(eps / 2, 1e-30)
    hi = w
    while True:
        lo = hi / 2
        pts.extend(panel(lo, hi, lev['n_g']))
        hi = lo
        if float(hi) < th_min:
            break
    pts.extend(panel(arb(0), hi, lev['n_g']))
    return pts


def legendre_matrix(th, Lmax):
    n = len(th)
    P = arb_mat(Lmax + 1, n)
    for k, (t, _) in enumerate(th):
        mu = t.cos()
        p0, p1 = arb(1), mu
        P[0, k] = p0
        if Lmax >= 1:
            P[1, k] = p1
        for L in range(1, Lmax):
            p0, p1 = p1, ((2 * L + 1) * mu * p1 - L * p0) / (L + 1)
            P[L + 1, k] = p1
    return P


def rho2_matrix(p, r_list, Lmax, th, P=None, chunk=150):
    """rho2_L(r_j) for L = 0..Lmax as a list of lists [j][L] (arb)."""
    if P is None:
        P = legendre_matrix(th, Lmax)
    a, b, zZ, zW = A(p['a']), A(p['b']), A(p['zZ']), A(p['zW'])
    S = [((t / 2).sin()) ** 2 for t, _ in th]
    wt = [w * t.sin() for t, w in th]
    n = len(th)
    out = []
    for c0 in range(0, len(r_list), chunk):
        rs = r_list[c0:c0 + chunk]
        V = arb_mat(n, len(rs))
        for jj, r in enumerate(rs):
            da2, db2 = (r - a) ** 2, (r - b) ** 2
            fa, fb = 4 * r * a, 4 * r * b
            for k in range(n):
                s = (da2 + fa * S[k]).sqrt()
                t = (db2 + fb * S[k]).sqrt()
                V[k, jj] = wt[k] * (-(zZ * s + zW * t)).exp()
        M = P * V
        for jj in range(len(rs)):
            out.append([arb((M[L, jj] * (2 * L + 1) / 2).mid()) for L in range(Lmax + 1)])
    return out


# ---------------------------------------------------------------- inner radial potential psi_L
def psi_at_nodes(p, nodes, br, Lmax, lev, rho1=None):
    rho1 = rho1 or Rho1(p, Lmax)
    D1 = p['D1']
    # merged sorted points: 0, nodes, breakpoints (> 0)
    pts = [(0.0, arb(0), 'b')] + [(float(r), r, j) for j, (r, _) in enumerate(nodes)] + \
          [(float(v), A(v), 'b') for v in br[1:]]
    pts.sort(key=lambda t: t[0])
    m = lev['m_sub']
    L1 = Lmax + 1
    IA, IB = [], []
    for (fp, rp, _), (fq, rq, _) in zip(pts[:-1], pts[1:]):
        ia = [arb(0)] * L1
        ib = [arb(0)] * L1
        if fq > fp:
            inner = (fp + fq) / 2 < float(D1)
            for s, w in panel(rp, rq, m):
                vals = rho1(s, inner)
                ws = w * s
                u, v = s / rq, (rp / s if fp > 0 else arb(0))
                pu, pv = u, arb(1)
                for L in range(L1):
                    t = ws * vals[L]
                    ia[L] += t * pu
                    ib[L] += t * pv
                    pu *= u
                    pv *= v
        IA.append(ia)
        IB.append(ib)
    n = len(pts)
    Ah = [None] * n
    Bh = [None] * n
    Ah[0] = [arb(0)] * L1
    for i in range(1, n):
        fp, rp = pts[i - 1][0], pts[i - 1][1]
        rq = pts[i][1]
        rat = rp / rq if fp > 0 else arb(0)
        pw = rat
        cur = []
        for L in range(L1):
            cur.append(Ah[i - 1][L] * pw + IA[i - 1][L])
            pw *= rat
        Ah[i] = cur
    Bh[n - 1] = [arb(0)] * L1
    for i in range(n - 2, -1, -1):
        fp, rp = pts[i][0], pts[i][1]
        rq = pts[i + 1][1]
        rat = rp / rq if fp > 0 else arb(0)
        pw = arb(1)
        cur = []
        for L in range(L1):
            cur.append(Bh[i + 1][L] * pw + IB[i][L])
            pw *= rat
        Bh[i] = cur
    psi = [None] * len(nodes)
    for i, (_, _, tag) in enumerate(pts):
        if tag != 'b':
            psi[tag] = [Ah[i][L] + Bh[i][L] for L in range(L1)]
    return psi


# ---------------------------------------------------------------- driver
def coefficients(case, Lmax, level='A', verbose=False):
    lev = LEVELS[level]
    old = ctx.prec
    ctx.prec = lev['prec']
    try:
        p = params(case)
        t0 = time.time()
        nodes, br = r_grid(p, lev)
        th = theta_grid(p, nodes, lev)
        P = legendre_matrix(th, Lmax)
        t1 = time.time()
        R2 = rho2_matrix(p, [r for r, _ in nodes], Lmax, th, P)
        t2 = time.time()
        psi = psi_at_nodes(p, nodes, br, Lmax, lev)
        t3 = time.time()
        pi4 = 4 * arb.pi()
        F = []
        for L in range(Lmax + 1):
            s = arb(0)
            for j, (r, w) in enumerate(nodes):
                s += w * r * r * R2[j][L] * psi[j][L]
            F.append((pi4 / (2 * L + 1)) ** 2 * s)
        t4 = time.time()
        info = dict(level=level, n_r=len(nodes), n_theta=len(th), n_panels_r=len(br) - 1,
                    t_setup=t1 - t0, t_rho2=t2 - t1, t_psi=t3 - t2, t_sum=t4 - t3, t_total=t4 - t0)
        if verbose:
            print(case, info, flush=True)
        return F, info
    finally:
        ctx.prec = old


def rho2_point(case, r, Lmax, level='A'):
    """rho2_L at a single radius (for controls), on the theta grid the engine would use for that r alone."""
    lev = LEVELS[level]
    p = params(case)
    rr = A(r) if not isinstance(r, arb) else r
    th = theta_grid(p, [(rr, None)], lev)
    return rho2_matrix(p, [rr], Lmax, th)[0]


def rho1_point(case, r, Lmax):
    p = params(case)
    rr = A(r) if not isinstance(r, arb) else r
    return Rho1(p, Lmax)(rr, Fraction(r) < p['D1'] if not isinstance(r, arb) else float(r) < float(p['D1']))


def legendre_sum(F, x):
    """sum_L F_L P_L(x) (arb)."""
    x = A(x) if not isinstance(x, arb) else x
    p0, p1 = arb(1), x
    s = F[0] * p0 + (F[1] * p1 if len(F) > 1 else 0)
    for L in range(1, len(F) - 1):
        p0, p1 = p1, ((2 * L + 1) * x * p1 - L * p0) / (L + 1)
        s += F[L + 1] * p1
    return s


if __name__ == '__main__':
    # smoke: low L only (no scored quantity)
    for c in ('G1', 'G2'):
        F, info = coefficients(c, 12, 'A', verbose=True)
        print(c, [F[L].mid().str(12, radius=False) for L in (0, 1, 4, 12)], landau(c))
