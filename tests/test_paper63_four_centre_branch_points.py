"""Guard for record FC (docs/predictions/2026-11-16_four_centre_branch_points.md, scored GO): the four-centre 1s
Slater integral (XY|ZW) on the one-ray slice, as a function of x = cos(angle between the ray XZW and the bond XY), has
a branch point (x* - x)^{7/2} at each moving contact (Y meets Z, Y meets W).  Paper 63 sec:four_centre, the paragraph
"[MEASURED] As a function of one angle, the four-center integral has one branch point at each moving contact".
Written 2026-11-16 as a separate step after the scoring (CLAUDE.md sec. 9).

Geometry: X at the origin, Y = D1 z-hat; Z = a n-hat, W = b n-hat (0 < a <= b); x = n-hat . z-hat.  Unnormalised
densities rho1 = e^{-zX|r-X|} e^{-zY|r-Y|}, rho2 = e^{-zZ|r-Z|} e^{-zW|r-W|}; I(x) = sum_L F_L P_L(x), with
    F_L = (4 pi/(2L+1))^2 int int rho1_L rho2_L r<^L / r>^(L+1) r1^2 r2^2,
rho1_L, rho2_L the Legendre coefficients of the densities about X.  The record's engine (debug/fc_engine.py, frozen
sha256 9fef5cf6...) is ported below unchanged in substance; only a reduced precision level 'T' is added for the fast
tier (the record's level A is used by the slow test).  Nothing is imported from debug/.

Claims guarded, each with the wrong answer it rejects (details in each docstring):
(a) the contact points x_Z* = (D1^2 + a^2)/(2 D1 a), x_W* = (D1^2 + b^2)/(2 D1 b) and rates rho = min/max, exact,
    for G1 and G2, and which one dominates.  Rejects a and b swapped, a wrong rho, a wrong dominant point.
(b) F_L > 0 for L <= 40 (forced by absolute monotonicity; a negative coefficient is a numerical bug).
(c) the three-centre reduction b = a (K4A1 = record AB's case A1) against AB's exact coefficients, L <= 12.
(d) the summed series at G1 and G2, x = 1/2, against the R2 Gaussian-fit referee (fc_controls.json K2, n = 26).
(e) beta_eff at L = 100 lies between the midpoints of in-test toys (x* - x)^nu, nu = 5/2, 7/2, 9/2, both cases;
    rejects 5/2 and 9/2.  The toys themselves are checked against the record's K6 values at L = 60.
(slow) beta_eff at L = 200, 220, 240 (record level A, which the record's FC0 shows agrees with level B to 2e-20)
    within the record's FC1/FC2 bands and within 0.02 of the 7/2 toy.

Fire-test plants and outcomes: the comment block at the end of this file.
"""
from __future__ import annotations

from fractions import Fraction

import pytest
from flint import arb, arb_mat, ctx

# --------------------------------------------------------------------------------------- ported engine (fc_engine)
CASES = {
    'G1': dict(a='1', b='9/4', D1='8/5', zX='13/10', zY='1', zZ='6/5', zW='1'),
    'G2': dict(a='1', b='4', D1='6/5', zX='13/10', zY='1', zZ='6/5', zW='7/10'),
    # K4 reduction: AB case A1 mapped (X = B, Y = C, Z = W = A at D2, zZ + zW = mu = 2 zeta_A = 9/5)
    'K4A1': dict(a='13/5', b='13/5', D1='9/5', zX='13/10', zY='1', zZ='1', zW='4/5'),
}

LEVELS = {
    'A': dict(prec=340, th_w=8, n_u=100, n_g=64, q=Fraction(5, 4), wmax=2, n_r=40, m_sub=12, r_hi=40, r0='1/20'),
    # reduced level for the fast tier (this guard only)
    'T': dict(prec=256, th_w=6, n_u=64, n_g=40, q=Fraction(3, 2), wmax=2, n_r=24, m_sub=10, r_hi=36, r0='1/20'),
}

PRED_NU = Fraction(7, 2)


def params(case):
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


class Rho1:
    """rho1_L exact: -(2L+1) e^{-zX r} [(2L+1) i_L(x) k_L(y) + x i_{L+1}(x) k_L(y) - y i_L(x) k_{L+1}(y)]."""

    def __init__(self, p, Lmax):
        self.p = p
        self.Lmax = Lmax
        self.zX, self.zY, self.D1 = A(p['zX']), A(p['zY']), A(p['D1'])
        c = self.zY * self.D1
        self.i_c = self.i_all(c)
        self.k_c = self.k_all(c)

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


def psi_at_nodes(p, nodes, br, Lmax, lev):
    rho1 = Rho1(p, Lmax)
    D1 = p['D1']
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


def coefficients(case, Lmax, level):
    lev = LEVELS[level]
    old = ctx.prec
    ctx.prec = lev['prec']
    try:
        p = params(case)
        nodes, br = r_grid(p, lev)
        th = theta_grid(p, nodes, lev)
        P = legendre_matrix(th, Lmax)
        R2 = rho2_matrix(p, [r for r, _ in nodes], Lmax, th, P)
        psi = psi_at_nodes(p, nodes, br, Lmax, lev)
        pi4 = 4 * arb.pi()
        F = []
        for L in range(Lmax + 1):
            s = arb(0)
            for j, (r, w) in enumerate(nodes):
                s += w * r * r * R2[j][L] * psi[j][L]
            F.append((pi4 / (2 * L + 1)) ** 2 * s)
        return F
    finally:
        ctx.prec = old


def legendre_sum(F, x):
    x = A(x)
    p0, p1 = arb(1), x
    s = F[0] * p0 + F[1] * p1
    for L in range(1, len(F) - 1):
        p0, p1 = p1, ((2 * L + 1) * x * p1 - L * p0) / (L + 1)
        s += F[L + 1] * p1
    return s


# --------------------------------------------------------------------------------------- toys and beta_eff
def toy_coeffs(x_star, nu, Lmax, n, prec):
    """Legendre coefficients c_L of (x* - x)^nu, L = 0..Lmax, by n-point Gauss-Legendre (record K6)."""
    old = ctx.prec
    ctx.prec = prec
    try:
        xa, nua = A(x_star), A(nu)
        c = [arb(0)] * (Lmax + 1)
        for x, w in gl(n):
            fv = w * ((xa - x).log() * nua).exp()
            p0, p1 = arb(1), x
            c[0] += fv
            c[1] += fv * x
            for L in range(1, Lmax):
                p0, p1 = p1, ((2 * L + 1) * x * p1 - L * p0) / (L + 1)
                c[L + 1] += fv * p1
        return [c[L] * (2 * L + 1) / 2 for L in range(Lmax + 1)]
    finally:
        ctx.prec = old


def dominant(case):
    lan = landau(case)
    if lan['dominant'] == 'Z':
        return lan['xZ'], lan['rhoZ']
    return lan['xW'], lan['rhoW']


def beta_eff(F, rho, L):
    """log|G_{L+1}/G_{L-1}| / log((L+1)/(L-1)), G_L = F_L / rho^L."""
    old = ctx.prec
    ctx.prec = 256
    try:
        r = A(rho)
        g = abs(F[L + 1] / F[L - 1]) / r ** 2
        return float((g.log() / (arb(L + 1) / (L - 1)).log()).mid())
    finally:
        ctx.prec = old


def toy_betas(case, L, n, prec):
    x_star, rho = dominant(case)
    out = {}
    for nu in (PRED_NU - 1, PRED_NU, PRED_NU + 1):
        out[nu] = beta_eff(toy_coeffs(x_star, nu, L + 1, n, prec), rho, L)
    return out


# --------------------------------------------------------------------------------------- shared fast series
L_FAST = 101
_FAST = {}


def fast(case):
    if case not in _FAST:
        _FAST[case] = coefficients(case, L_FAST, 'T')
    return _FAST[case]


def flo(x):
    return float(x.mid())


# --------------------------------------------------------------------------------------- (a)
def test_landau_points_exact():
    """Rejects: a and b swapped in x_Z*, x_W* or rho; rho = D1/a instead of min/max; the wrong dominant point.
    Values are the record's case table: G1 x_Z* 89/80, rho_Z 5/8, x_W* 3049/2880, rho_W 32/45 (W dominates);
    G2 x_Z* 61/60, rho_Z 5/6, x_W* 109/60, rho_W 3/10 (Z dominates).  Also checks x* = (1/rho + rho)/2, the
    identity |Y - C|^2 = 2 D1 c (x* - x) that ties each contact to its rate."""
    want = {'G1': (Fraction(89, 80), Fraction(3049, 2880), Fraction(5, 8), Fraction(32, 45), 'W'),
            'G2': (Fraction(61, 60), Fraction(109, 60), Fraction(5, 6), Fraction(3, 10), 'Z')}
    for case, (xz, xw, rz, rw, dom) in want.items():
        lan = landau(case)
        assert (lan['xZ'], lan['xW'], lan['rhoZ'], lan['rhoW'], lan['dominant']) == (xz, xw, rz, rw, dom), lan
        for xs, r in ((lan['xZ'], lan['rhoZ']), (lan['xW'], lan['rhoW'])):
            assert xs == (1 / r + r) / 2


# --------------------------------------------------------------------------------------- (b)
@pytest.mark.parametrize('case', ['G1', 'G2'])
def test_positivity(case):
    """F_L > 0 for L <= 40 (forced: both densities are absolutely monotone in the angle cosine and the radial kernel
    is positive).  Rejects a sign error in a multipole or kernel that leaves some F_L negative."""
    F = fast(case)
    bad = [L for L in range(41) if not F[L] > 0]
    assert not bad, bad


# --------------------------------------------------------------------------------------- (c)
# Record AB's exact coefficients for case A1 (debug/aabc_exact.py: pieces('A1', 12, Field(None)), combine with
# zA = 9/10, constants evaluated by na_exact.value at 500 digits), multiplied by the engine's convention factor
# 2 pi^2 / zA^3 (record FC control K4: engine / (sum q c) = 2 pi^2 / zA^3).  Computed once, 2026-11-16.
AB_A1 = [
    '1.784622221920514515562977330690395639463',
    '0.3576581824093901743389945608179418988148',
    '0.07853129974705347101569654491914415867192',
    '0.01878453879306009227471141596868396261614',
    '0.004930401130712857549253065333979139887916',
    '0.001421642046166370948657193074119677743983',
    '0.0004481586050347342931632883555121610893503',
    '0.0001531937066429101486219221142859261206234',
    '0.00005624962928729676211238542032527996552614',
    '0.00002197988566210178889421019101881107700954',
    '0.000009062951163248631689784980822715590756306',
    '0.000003914064901991471434405539892511744660374',
    '0.000001759330290957584915750390150352648150129',
]


def test_three_centre_reduction():
    """b = a, Z = W (K4A1) reproduces record AB's exact (AA|BC) coefficients for L <= 12 to 1e-20 relative (record
    K4: 1.1e-28 at level A).  Rejects any error in the radial kernel, the multipole normalisation (4 pi/(2L+1))^2,
    rho1_L, or the 2L+1 in rho2_L that survives the three-centre limit."""
    F = coefficients('K4A1', 12, 'T')
    old = ctx.prec
    ctx.prec = 256
    try:
        devs = [abs(flo(abs(F[L] / arb(AB_A1[L]) - 1))) for L in range(13)]
    finally:
        ctx.prec = old
    assert max(devs) < 1e-20, devs


# --------------------------------------------------------------------------------------- (d)
# R2 referee values, debug/data/fc_controls.json K2 (n = 26 even-tempered Gaussian fit; n = 22 vs 26 spread 3.2e-12
# at G1, 1.7e-9 at G2).  The engine sum (level A, L <= 260) differs from these by -1.07e-11 (G1) and -6.24e-11 (G2).
R2_REF = {'G1': 1.220722671476179, 'G2': 1.042916834956128}
R2_TOL = {'G1': 3e-11, 'G2': 2e-10}


@pytest.mark.parametrize('case', ['G1', 'G2'])
def test_sum_vs_referee(case):
    """sum_{L <= 101} F_L P_L(1/2) against R2 within 3e-11 (G1) and 2e-10 (G2); the tail beyond L = 101 is below
    1e-15.  Rejects a perturbed referee (1e-9), any L-independent factor, a dropped radial breakpoint, and every
    engine error that moves the low-L coefficients by more than about 1e-10."""
    s = flo(legendre_sum(fast(case), '1/2'))
    assert abs(s - R2_REF[case]) < R2_TOL[case], (s, s - R2_REF[case])


# --------------------------------------------------------------------------------------- (e)
# Record K6 toy beta_eff at L = 60 (debug/data/fc_controls.json), for checking the in-test toy instrument.
K6_L60 = {'G1': {Fraction(5, 2): -3.26489030459136, Fraction(7, 2): -4.450250589822434,
                 Fraction(9, 2): -5.683958542488235},
          'G2': {Fraction(5, 2): -3.47327456068501, Fraction(7, 2): -4.789954090652098,
                 Fraction(9, 2): -6.18003844695905}}


@pytest.mark.parametrize('case', ['G1', 'G2'])
def test_toys_reproduce_K6(case):
    """The in-test toys reproduce the record's K6 beta_eff at L = 60 for nu = 5/2, 7/2, 9/2 to 1e-9.  Rejects a toy
    at the wrong x* or rho, or a broken toy recurrence (which would move the bands of the decay test)."""
    x_star, rho = dominant(case)
    for nu, want in K6_L60[case].items():
        got = beta_eff(toy_coeffs(x_star, nu, 61, 220, 256), rho, 60)
        assert abs(got - want) < 1e-9, (nu, got, want)


@pytest.mark.parametrize('case', ['G1', 'G2'])
def test_decay_7_2_at_L100(case):
    """beta_eff(100) of the engine (level T) lies between the toy midpoints around nu = 7/2 (5/2-7/2 and 7/2-9/2) and
    within 0.04 of the 7/2 toy (measured 2026-11-16: -4.2636 vs toy -4.2760 at G1, -4.4984 vs -4.4991 at G2; record
    level B: -4.264, -4.498; the toy spacing is about 1.1).  Rejects a 5/2 or 9/2 branch (planting PRED_NU), a
    wrong rate or contact, and an L-dependent factor such as (2L+1)^{+-1} in the normalisation (shifts beta_eff by
    about +-1)."""
    F = fast(case)
    _, rho = dominant(case)
    b = beta_eff(F, rho, 100)
    t = toy_betas(case, 100, 260, 256)
    lo = (t[PRED_NU] + t[PRED_NU + 1]) / 2
    hi = (t[PRED_NU] + t[PRED_NU - 1]) / 2
    assert lo < b < hi, (b, lo, hi, t)
    assert abs(b - t[PRED_NU]) < 0.04, (b, t)


# --------------------------------------------------------------------------------------- slow: L = 220
# Record FC1/FC2 bands (built from the K6 toy midpoints at L = 220).
FC1_BAND = {'G1': (-4.661, -3.602), 'G2': (-4.796, -3.688)}


@pytest.mark.slow
@pytest.mark.parametrize('case', ['G1', 'G2'])
def test_decay_7_2_at_L220_slow(case):
    """Level A to L = 241.  beta_eff at L = 200, 220, 240 lies in the record's FC1/FC2 band and between the in-test
    toy midpoints at L = 220; |beta_eff(220) - toy 7/2| < 0.02 (record: 0.005 at G1, 0.0003 at G2).  Rejects 5/2,
    9/2, a subdominant-rate error of 1% or more, and a slowly growing L-dependent factor invisible at L = 100."""
    F = coefficients(case, 241, 'A')
    _, rho = dominant(case)
    t = toy_betas(case, 220, 700, 640)
    lo = (t[PRED_NU] + t[PRED_NU + 1]) / 2
    hi = (t[PRED_NU] + t[PRED_NU - 1]) / 2
    for L in (200, 220, 240):
        b = beta_eff(F, rho, L)
        assert FC1_BAND[case][0] < b < FC1_BAND[case][1], (L, b)
        assert lo < b < hi, (L, b, lo, hi)
    assert abs(beta_eff(F, rho, 220) - t[PRED_NU]) < 0.02, (beta_eff(F, rho, 220), t)


# TIMINGS (2026-11-16, this machine): fast tests 17 s in all (level T to L = 101: ~6.5 s per case, shared by the
# positivity, referee and decay tests; K4A1 to L = 12 3 s; toys < 0.2 s each).  Slow tests 104 s (--slow; level A to
# L = 241, ~52 s per case, toys at n = 700 included).
#
# FIRE-TEST PLANTS (python debug/qa/fire_test.py tests/test_paper63_four_centre_branch_points.py --plant "OLD=>NEW"
# -k <sel>), 2026-11-16, written and fired as a separate step after the guard; 34 runs, 30 fired; the 4 that did not are explained below.
#     L-dependent normalisation          "(pi4 / (2 * L + 1)) ** 2 * s)" -> "... * s * (L + 1))"   -k "reduction or sum_vs or decay_7_2_at_L100"   FIRED; -k decay_7_2_at_L100 alone   FIRED
#     missing (2L+1) in F_L prefactor    "(pi4 / (2 * L + 1)) ** 2 * s)" -> "(pi4 / (2 * L + 1)) * s)"   -k "reduction or sum_vs or decay_7_2_at_L100"   FIRED
#     missing (2L+1) in rho2_L           "(M[L, jj] * (2 * L + 1) / 2)" -> "(M[L, jj] / 2)"            -k reduction   FIRED; -k decay_7_2_at_L100   FIRED
#     x_Z* with a and b swapped          "xZ = (D1 ** 2 + a ** 2) / (2 * D1 * a)" -> b                -k landau   FIRED; -k decay_7_2_at_L100   FIRED (G2)
#     rho_Z with a -> b                  "rZ = min(D1, a) / max(D1, a)" -> min(D1, b) / max(D1, b)    -k landau   FIRED; -k decay_7_2_at_L100   FIRED (G2)
#     rho_Z inverted                     "rZ = min(D1, a) / max(D1, a)" -> max / min                  -k landau   FIRED
#     9/2 toy                            "PRED_NU = Fraction(7, 2)" -> Fraction(9, 2)                 -k decay_7_2_at_L100   FIRED; -k L220   FIRED
#     5/2 toy                            "PRED_NU = Fraction(7, 2)" -> Fraction(5, 2)                 -k decay_7_2_at_L100   FIRED; -k L220   FIRED
#     kernel sign, rho1_L                "- y * i[L] * k[L + 1]))" -> "+ y ..."                       -k "positivity or reduction or sum_vs"   FIRED
#     kernel sign, rho2 exponent         "(-(zZ * s + zW * t)).exp()" -> "(-(zZ * s - zW * t))"       -k "reduction or sum_vs"   FIRED
#     radial kernel r<^L instead of ^(L+1)  "pu, pv = u, arb(1)" -> "arb(1), arb(1)"                  -k "reduction or sum_vs"   FIRED
#     radial split dropped at a          p['a'] removed from the r_breaks keys                       -k sum_vs   FIRED
#     radial split dropped at D1         p['D1'] removed from the r_breaks keys                      -k "reduction or sum_vs"   FIRED
#     referee G1 perturbed by 1e-9       1.220722671476179 -> 1.220722672476179                      -k sum_vs   FIRED
#     referee G2 perturbed by 1e-9       1.042916834956128 -> 1.042916835956128                      -k sum_vs   FIRED
#     4 pi -> 2 pi                       "pi4 = 4 * arb.pi()" -> "2 * arb.pi()"                       -k "reduction or sum_vs"   FIRED
#     i_L downward recurrence sign       "i[L + 1] + (2 * L + 1) * ix * i[L]" -> "-"                  -k "reduction or sum_vs"   FIRED
#     theta grading removed              "th_min = max(eps / 2, 1e-30)" -> "max(eps / 2, 0.3)"        -k "reduction or sum_vs"   FIRED
#     AB literal L = 5 off by 7e-19 rel  digit 19 of AB_A1[5] changed                                -k reduction   FIRED
#     beta_eff dividing by rho, not rho^2  "/ r ** 2" -> "/ r"                                        -k toys_reproduce   FIRED
#     toy (2L+1)/2 normalisation dropped "c[L] * (2 * L + 1) / 2" -> "c[L] / 2"                       -k toys_reproduce   FIRED; -k decay_7_2_at_L100   FIRED
#     toy P_L recurrence                 "/ (L + 1)" -> "/ (L + 2)" in toy_coeffs                     -k toys_reproduce   FIRED
#     Legendre sum truncated at L = 12   "range(1, len(F) - 1)" -> "range(1, 12)"                     -k sum_vs   FIRED
#   Did not fire, and why (not guard defects):
#     theta grading cut at 1e-3 instead of 1e-30 (-k "reduction or sum_vs"): still graded to 1e-3 plus a 40-node end
#       panel, enough for 1e-20 (K4A1, L <= 12) and the referee's 1e-10; the grading depth is precision margin.
#     beta_eff "/ r" with -k decay_7_2_at_L100: engine and toys share beta_eff, so the band is invariant to it; the
#       same plant fires toys_reproduce (fixed K6 literals).
#     toy "c[1] += fv * x * 2": touches L = 1 only, never read by a beta_eff at L >= 60.
#     Legendre sum truncated at L = 30: the tail beyond L = 30 is below the referee tolerance (3e-11) at x = 1/2.
