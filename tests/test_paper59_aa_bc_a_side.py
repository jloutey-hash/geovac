"""Guard for record AS (docs/predictions/2026-11-15_aa_bc_a_side_ppi_d.md, scored GO): the three-centre integrals

    PY  = (2p_y,A  2p_y,A  | 1s_B 1s_C)    p_pi at A out of the plane of the centres
    PX' = (2p_x',A 2p_x',A | 1s_B 1s_C)    p_pi at A in the plane, perpendicular to B-A
    DA  = (3d_sigma,A 3d_sigma,A | 1s_B 1s_C)   sigma axis along B-A

over Slater functions, as series in x = cos(angle C B A).  Paper 59 sec:scope, the paragraph "[MEASURED] Changing
the pair on A instead".  Written 2026-11-15 as a separate step after the computation and the paper edit (CLAUDE.md
sec. 9).  The record's engine (debug/aaa_exact.py) is ported here; its base layer (Field, Bessel tables, the
exponential integrals, AB's 1s density, the AP kernels and the three-region integrator) is imported from
tests/test_paper59_aa_bc_p.py, which imports tests/test_paper59_aa_bc.py.

PY and PX' are sum_L a_L P_L(x) + sum_L b_L P_L^2(x), P_L^2 = (1 - x^2) P_L'' (no Condon-Shortley phase), (a, b)
the 3-D harmonic split of PY (m = 0 and m = 2 about the B-C axis).  The total as a plain Legendre series is
c_j = a_j - j(j-1) b_j + 2(2j+1) sum_{L > j, L = j mod 2} b_L.  DA = sum_L F_L P_L(x).

Claims guarded, each with the wrong answer it rejects (details in each test's docstring):

(1) per-channel values: a_L, b_L of PY and F_L of DA equal an independent quadrature (the A potential written as
    its Gauss's-law multipoles V_l of the A density, incomplete-gamma form; the lab tensor M_ij for PY; no U
    kernels, no D2-derivatives, no jet algebra).  Rejects a wrong angular factor of the out-of-plane harmonic, a
    wrong trace term, a wrong kernel coefficient, a wrong derivative order.
(2) the summed series at A1, x = 1/2, against the R2 Gaussian-fit referee (n = 34) and the record's 3600-digit
    values (L <= 400), for PY, PX' and DA; the re-expansion of the total c is an exact identity at truncation.
(3) the constant set: mod 2^61 - 1, every coefficient up to L = 20 (a, F from L = 2, b from L = 3) uses exactly
    AB's 21 constants at the same case; no Euler gamma.
(4) decay: beta_eff at L = 40 (fast, A1) and L = 60 and 80 (slow, both cases) sits between the midpoints of the
    in-test toy exponents: PY a, PX' total, DA on the (x* - x)^{7/2} branch; PY b on the 7/2 branch in the P_L^2
    basis; PY total c on the (x* - x)^{9/2} branch.  Rejects "no cancellation" (PY total on 7/2).
(5) P4: R = (a/ta)/(b/tb) positive, and its 1/L extrapolation from two cutoffs lands near x*^2 - 1; the band at the
    tested L is set from the in-test values, not from L = 220.

Fire-test plants and outcomes: the comment block at the end of this file.
"""
from __future__ import annotations

import importlib.util
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

import mpmath as mp
import pytest

_spec = importlib.util.spec_from_file_location('_ap_guard_as', Path(__file__).with_name('test_paper59_aa_bc_p.py'))
AP = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(AP)
AB = AP.AB

Field, params, add, _put, _conv, _addm = AB.Field, AB.params, AB.add, AB._put, AB._conv, AP._addm
PRIME = (1 << 61) - 1


# ---------------------------------------------------------------------------------------------------------------
# Record-AS engine (port of debug/aaa_exact.py; PY, PX', DA only)
# ---------------------------------------------------------------------------------------------------------------


def fam_p(zA):
    """2p pair on A: V_e = (e.grad)^2 U2(|w|) + U0(|w|); 'NA' = 1/s, integer k = d^k Y / dmu^k."""
    z = zA
    U2 = {'NA': Fraction(3, 2) / z ** 2, 0: -Fraction(3, 2) / z ** 2, 1: Fraction(33, 16) / z, 2: Fraction(-9, 8),
          3: z / 4}
    U0 = {'NA': Fraction(1), 0: Fraction(-1), 1: Fraction(5, 4) * z, 2: -z * z / 2}
    return {'U2': U2, 'U0': U0}


def fam_d(zA):
    """3d-sigma pair on A: V_n = (n.grad)^4 U4 + (n.grad)^2 U2 + U0."""
    z = zA
    U4 = {'NA': Fraction(15, 4) / z ** 4, 0: -Fraction(15, 4) / z ** 4, 1: Fraction(2895, 512) / z ** 3,
          2: -Fraction(975, 256) / z ** 2, 3: Fraction(185, 128) / z, 4: Fraction(-5, 16), 5: z / 32}
    U2 = {'NA': 2 / z ** 2, 0: -2 / z ** 2, 1: Fraction(151, 64) / z, 2: Fraction(-23, 32), 3: -Fraction(17, 48) * z,
          4: z ** 2 / 3, 5: -z ** 3 / 12}
    U0 = {'NA': Fraction(1), 0: Fraction(-1), 1: Fraction(137, 96) * z, 2: -Fraction(41, 48) * z ** 2,
          3: Fraction(7, 24) * z ** 3, 4: -z ** 4 / 9, 5: z ** 5 / 18}
    return {'U4': U4, 'U2': U2, 'U0': U0}


def kernel_jet(L, near, coeffs, it, kt, K, mu, D2, nmax):
    """[kernel's j-th D2-derivative, {(e0, b): {n: c}}] for j = 0..nmax; any mu-order (factored, Leibniz)."""
    outs = [dict() for _ in range(nmax + 1)]
    if 'NA' in coeffs:
        na = {k: coeffs['NA'] * v for k, v in AP.kernel_mono(L, near, 'NA', it, kt, K).items()}
        for j in range(nmax + 1):
            if j:
                na = AP.d_D2(na, K)
            for key, v in AP.evaluate(na, mu, D2, K).items():
                for n, c in v.items():
                    _put(outs[j], key[0], key[1], n, c)
    ords = [k for k in coeffs if k != 'NA']
    if not ords:
        return outs
    amax = max(ords)
    G0, H0 = AP._factors(L, near, it, kt, K)
    Gs, Hs = [G0], [H0]
    for _ in range(amax):
        Gs.append(AP._dG(Gs[-1], K, 'mu'))
        Hs.append(AP._dH(Hs[-1], K))
    Hv = []
    for h in Hs:
        v = {}
        for (pa, n, s), c in h.items():
            _addm(v, (s, n), c * K.pw(mu, pa))
        Hv.append(v)
    for j in range(nmax + 1):
        if j:
            Gs = [AP._dG(g, K, 'D2') for g in Gs]
        Gv = []
        for g in Gs:
            v = {}
            for (pa, pd, t0), c in g.items():
                _addm(v, t0, c * K.pw(mu, pa) * K.pw(D2, pd))
            Gv.append(v)
        acc = {}
        for o, a in coeffs.items():
            if o == 'NA':
                continue
            for k in range(o + 1):
                w = a * K(comb(o, k))
                for t0, gv in Gv[k].items():
                    wg = w * gv
                    for (s, n), hv in Hv[o - k].items():
                        _addm(acc, (t0, s, n), wg * hv)
        for (t0, s, n), v in acc.items():
            _put(outs[j], t0 * mu * D2, s * mu, n, v)
    return outs


def boundary_gen(dens_outer, jet_near, jet_far, n, D2, zB, K):
    """Jump terms of the n-th D2-derivative at r = D2: sum_i C(n, i+1) d_r^i [rho (K_near - K_far)^(n-1-i)]."""
    res = {}
    for i in range(n):
        j = n - 1 - i
        w = K(comb(n, i + 1))
        jump = {}
        for src, sg in ((jet_near[j], 1), (jet_far[j], -1)):
            for (e0, b), poly in src.items():
                for m, v in poly.items():
                    _put(jump, e0, b, m, K(sg) * v)
        for (e0d, bd), pd in dens_outer.items():
            for (e0k, bk), pk in jump.items():
                Pc = _conv(pd, pk, K)
                bb = bd + bk - zB
                e0 = e0d + e0k
                for m, v in Pc.items():
                    q = m + 2
                    ff = 1
                    for a in range(i + 1):
                        if a:
                            ff *= (q - a + 1)
                        if ff == 0:
                            break
                        if bb == 0 and i - a > 0:
                            continue
                        c = K(-2) * w * v * K(comb(i, a) * ff) * (K.pw(bb, i - a) if i - a else K(1)) * K.pw(D2, q - a)
                        add(res, (e0 + bb * D2, '1', None), c)
    return res


def _sumd(*ds):
    out = {}
    for d in ds:
        for k, v in d.items():
            add(out, k, v)
    return out


class Engine:
    """G_l^(j)(D2) = d^j/dD2^j int rho_BC U(|r - A|), Legendre channel l, cached per (kernel, l)."""

    def __init__(self, case, K, Lmax):
        self.case, self.K = case, K
        self.D1, self.D2, self.zB, self.zC, self.zA = params(case)
        self.mu = 2 * self.zA
        self.it, self.kt = AB.bessel_tables(K, Lmax + 8)
        self.I = AB.Integrals(K)
        self.dens_cache, self.cache = {}, {}

    def dens(self, L):
        if L not in self.dens_cache:
            self.dens_cache[L] = {inner: AB.density(L, inner, self.D1, self.zC, self.it, self.kt, self.K)
                                  for inner in (True, False)}
        return self.dens_cache[L]

    def G(self, name, fam, L, nmax):
        key = (name, tuple(sorted((str(k), v) for k, v in fam.items())), L, nmax)
        if key in self.cache:
            return self.cache[key]
        K = self.K
        co = {k: K(v) for k, v in fam.items()}
        jets = {near: kernel_jet(L, near, co, self.it, self.kt, K, self.mu, self.D2, nmax) for near in (True, False)}
        d = self.dens(L)
        out = []
        for j in range(nmax + 1):
            bulk = AP.integrate(d, {near: jets[near][j] for near in (True, False)}, self.D1, self.D2, self.zB,
                                self.I, K)
            if j:
                bulk = _sumd(bulk, boundary_gen(d[False], jets[True], jets[False], j, self.D2, self.zB, K))
            out.append({k: v for k, v in bulk.items() if v != 0})
        self.cache[key] = out
        return out


# jets: {(l, j, p): coef} = coef R^{-p} G_l^{(j)}(R)


def _jadd(out, k, v):
    if v:
        out[k] = out.get(k, 0) + v


def jd(jet):
    out = {}
    for (l, j, p), c in jet.items():
        _jadd(out, (l, j, p + 1), -p * c)
        _jadd(out, (l, j + 1, p), c)
    return out


def jscale(jet, a):
    return {k: a * v for k, v in jet.items() if a * v}


def jshift(jet, a):
    return {(l, j, p + 1): a * v for (l, j, p), v in jet.items() if a * v}


def jsum(*jets):
    out = {}
    for jt in jets:
        for k, v in jt.items():
            _jadd(out, k, v)
    return {k: v for k, v in out.items() if v}


def opA(l, jet):
    return jsum(jd(jet), jshift(jet, Fraction(-l)))


def opB(l, jet):
    return jsum(jd(jet), jshift(jet, Fraction(l + 1)))


def apply_dz(field):
    out = {}
    for l, jt in field.items():
        out[l + 1] = jsum(out.get(l + 1, {}), jscale(opA(l, jt), Fraction(l + 1, 2 * l + 1)))
        if l >= 1:
            out[l - 1] = jsum(out.get(l - 1, {}), jscale(opB(l, jt), Fraction(l, 2 * l + 1)))
    return out


def apply_dplus(field, m):
    out = {}
    for l, jt in field.items():
        out[l + 1] = jsum(out.get(l + 1, {}), jscale(opA(l, jt), Fraction(1, 2 * l + 1)))
        if l - 1 >= m + 1:
            out[l - 1] = jsum(out.get(l - 1, {}), jscale(opB(l, jt), Fraction(-1, 2 * l + 1)))
    return out


def apply_lap(field):
    out = {}
    for l, jt in field.items():
        out[l] = jsum(out.get(l, {}), jd(jd(jt)), jshift(jd(jt), Fraction(2)),
                      jshift(jshift(jt, Fraction(1)), Fraction(-l * (l + 1))))
    return out


def py_jets(L):
    srcs = {l: {(l, 0, 0): Fraction(1)} for l in (L - 2, L, L + 2) if l >= 0}
    lap = apply_lap(srcs).get(L, {})
    dzz = apply_dz(apply_dz(srcs)).get(L, {})
    beta = apply_dplus(apply_dplus(srcs, 0), 1).get(L, {}) if L >= 2 else {}
    return {'lap': lap, 'dzz': dzz, 'beta': beta}


def jet_value(jet, Gfun, D2, K):
    out = {}
    for (l, j, p), c in jet.items():
        a = K(c) * K.pw(D2, -p)
        for k, v in Gfun(l)[j].items():
            add(out, k, a * v)
    return {k: v for k, v in out.items() if v != 0}


def coefficients(case, Lmax, K, integral, eng=None):
    """PY / PXp -> (a, b), lists over L of {key: q}; DA -> list over L.  Value = 2 pi N_B N_C sum q c."""
    eng = eng or Engine(case, K, Lmax + 2)
    if integral in ('PY', 'PXp'):
        fam = fam_p(eng.zA)
        W = [eng.G('U0', fam['U0'], L, 0)[0] for L in range(Lmax + 1)]

        def Gf(l):
            return eng.G('U2', fam['U2'], l, 2)
        a, b = [], []
        half = Fraction(1, 2)
        for L in range(Lmax + 1):
            J = py_jets(L)
            mzero = jscale(jsum(J['lap'], jscale(J['dzz'], Fraction(-1))), half)     # (Lap - d_z^2)/2
            if integral == 'PY':
                ja, jb = mzero, jscale(J['beta'], -half)                                # -(d_+^2 + d_-^2)/4
            else:
                ja = jsum(J['lap'], {(L, 2, 0): Fraction(-1)}, jscale(mzero, Fraction(-1)))   # Lap - d_n^2 - d_y^2
                jb = jscale(J['beta'], half)
            a.append({k: v for k, v in _sumd(jet_value(ja, Gf, eng.D2, K), W[L]).items() if v != 0})
            b.append(jet_value(jb, Gf, eng.D2, K))
        return a, b
    if integral == 'DA':
        fam = fam_d(eng.zA)
        out = []
        for L in range(Lmax + 1):
            t = _sumd(eng.G('U4', fam['U4'], L, 4)[4], eng.G('U2', fam['U2'], L, 2)[2],
                      eng.G('U0', fam['U0'], L, 0)[0])
            out.append({k: v for k, v in t.items() if v != 0})
        return out
    raise ValueError(integral)


def _mf(q):
    q = Fraction(q)
    return mp.mpf(q.numerator) / q.denominator


def numeric(Q, case, dps):
    """2 pi N_B N_C sum q c (the A pair's normalization is inside the kernels)."""
    with mp.workdps(dps):
        zB, zC = params(case)[2:4]
        tw = 2 * mp.sqrt(_mf(zB) ** 3 * _mf(zC) ** 3)
        keys = set(k for q in Q for k in q)
        cv = {k: AB.const_value(k) for k in keys}
        return [tw * sum((mp.mpf(int(v.p)) / int(v.q) * cv[k] for k, v in q.items()), mp.mpf(0)) for q in Q]


def total_c(a, b):
    """c_j = a_j - j(j-1) b_j + 2(2j+1) sum_{L > j, L = j mod 2} b_L  (P_L^2 re-expanded in P_L, truncated)."""
    n = len(a)
    tail = [mp.mpf(0), mp.mpf(0)]
    c = [mp.mpf(0)] * n
    for j in range(n - 1, -1, -1):
        c[j] = a[j] - j * (j - 1) * b[j] + 2 * (2 * j + 1) * tail[j % 2]
        tail[j % 2] += b[j]
    return c


def legendre_pair(Lmax, x):
    """P_L(x) and P_L^2(x) = (1 - x^2) P_L''(x) = 2 x P_L' - L(L+1) P_L (Legendre ODE), L = 0..Lmax."""
    P = [mp.mpf(1), x]
    for L in range(1, Lmax):
        P.append(((2 * L + 1) * x * P[L] - L * P[L - 1]) / (L + 1))
    P = P[:Lmax + 1]
    om = 1 - x * x
    P2 = [mp.mpf(0)] * min(2, Lmax + 1) + [2 * x * L * (P[L - 1] - x * P[L]) / om - L * (L + 1) * P[L]
                                          for L in range(2, Lmax + 1)]
    return P, P2


_EXACT = {}


def exact_Q(case, integral, Lmax):
    """Exact coefficient dicts over Q, cached; PY and PXp share one engine per case."""
    key = (case, integral)
    if key not in _EXACT or _EXACT[key][0] < Lmax:
        if integral in ('PY', 'PXp'):
            eng = Engine(case, Field(None), Lmax + 2)
            for it_ in ('PY', 'PXp'):
                _EXACT[(case, it_)] = (Lmax, coefficients(case, Lmax, Field(None), it_, eng=eng))
        else:
            _EXACT[key] = (Lmax, coefficients(case, Lmax, Field(None), integral))
    return _EXACT[key][1]


def series(case, integral, Lmax, dps):
    """{'a', 'b', 'c'} (p pairs) or {'F'} (DA) as mpf lists, L = 0..Lmax."""
    r = exact_Q(case, integral, Lmax)
    if integral == 'DA':
        return {'F': numeric(r[:Lmax + 1], case, dps)}
    a, b = numeric(r[0][:Lmax + 1], case, dps), numeric(r[1][:Lmax + 1], case, dps)
    with mp.workdps(dps):
        return {'a': a, 'b': b, 'c': total_c(a, b)}


# ---------------------------------------------------------------------------------------------------------------
# Independent route (no U kernels, no D2-derivatives, no jets): Gauss's-law multipoles of the A density
# ---------------------------------------------------------------------------------------------------------------


def V_multipole(kind, l, s, z):
    """V_l(s) of the normalized A pair density: p: (z^5/pi) (e.w)^2 e^{-2 z s}, (e.w)^2 = s^2 (1/3 + 2/3 P2);
    d: N_d^2 s^4 P2(u)^2 e^{-2 z s}, P2^2 = 1/5 + 2/7 P2 + 18/35 P4.  f_l(t) = c_l t^k e^{-mu t}, mu = 2 z;
    V_l = 4 pi/(2l+1) [s^{-l-1} int_0^s t^{l+2} f_l + s^l int_s^oo t^{1-l} f_l], lower/upper incomplete gamma."""
    mu = 2 * z
    if kind == 'p':
        k, cl = 2, {0: z ** 5 / mp.pi / 3, 2: 2 * z ** 5 / mp.pi / 3}
    else:
        Nd2 = 5 * (2 * z) ** 7 / (4 * mp.pi * 720)
        k, cl = 4, {0: Nd2 / 5, 2: Nd2 * mp.mpf(2) / 7, 4: Nd2 * mp.mpf(18) / 35}
    ai, ao = l + 3 + k, 2 - l + k
    inner = mp.gammainc(ai, 0, mu * s) / mu ** ai
    outer = mp.gammainc(ao, mu * s, mp.inf) / mu ** ao
    return 4 * mp.pi / (2 * l + 1) * cl[l] * (s ** (-l - 1) * inner + s ** l * outer)


def _legendre_all(lmax, t):
    P = [mp.mpf(1), t]
    for n in range(1, lmax):
        P.append(((2 * n + 1) * t * P[n] - n * P[n - 1]) / (n + 1))
    return P[:lmax + 1]


def _antideriv(c):
    out = {}
    for k, v in c.items():
        out[k + 1] = out.get(k + 1, 0) + v / (2 * k + 1)
        if k >= 1:
            out[k - 1] = out.get(k - 1, 0) - v / (2 * k + 1)
    return out


def _p2_coeffs(c):
    """(1 - x^2) sum c_k P_k = sum_L e_L P_L^2: e = Legendre coefficients of the double antiderivative, L >= 2."""
    h = _antideriv(_antideriv(c))
    return {L: v for L, v in h.items() if L >= 2 and v != 0}


def _times_x(c):
    out = {}
    for k, v in c.items():
        out[k + 1] = out.get(k + 1, 0) + v * Fraction(k + 1, 2 * k + 1)
        if k >= 1:
            out[k - 1] = out.get(k - 1, 0) + v * Fraction(k, 2 * k + 1)
    return out


def direct_quadrature(case, Ls, dps, deg_in=5, deg_out=5):
    """{'a', 'b', 'DA'}: {L: value} by quadrature (port of the record's AS0 check, debug/aaa_ap0.py).
    PY = e_i e_j M_ij(A), M_ij = int rho_BC [Q0(s) delta_ij + q(s) w_i w_j], Q0 = V0 - V2/2, q = 3 V2/(2 s^2),
    w = r - A, A in the xz-plane: m = 0 part (tr M - M_zz)/2, m = 2 part (M_yy - M_xx)/2 (rotation about z), each
    projected through the addition theorem in t = r.A/(rR) with the BC density's tau-moments taken directly.
    DA: V(r - A) = sum_l V_l(s) P_l(u), u = (r t - R)/s, projected as plain Legendre.
    Gauss-Legendre in s (inner) and r on [0, D1, D2, D2+3, D2+10, D2+25, D2+70] (outer)."""
    with mp.workdps(dps):
        D1, D2, zB, zC, zA = (_mf(v) for v in params(case))
        R = D2
        NB, NC = mp.sqrt(zB ** 3 / mp.pi), mp.sqrt(zC ** 3 / mp.pi)
        GL = mp.calculus.quadrature.GaussLegendre(mp.mp)
        n_in = GL.calc_nodes(deg_in, mp.mp.prec)
        n_out = GL.calc_nodes(deg_out, mp.mp.prec)
        lmax = max(Ls) + 2
        names = ('T', 'al', 'be', 'ga', 'X0', 'X1', 'X2', 'DA')
        acc = {k: [mp.mpf(0)] * (lmax + 1) for k in names}
        edges = [mp.mpf(0), D1, D2, D2 + 3, D2 + 10, D2 + 25, D2 + 70]
        for lo, hi in zip(edges[:-1], edges[1:]):
            hh, cc = (hi - lo) / 2, (hi + lo) / 2
            for xo, wo in n_out:
                r = cc + hh * xo
                wr = wo * hh * r * r
                m0 = [mp.mpf(0)] * (lmax + 1)
                m1, m2, x1, x2 = list(m0), list(m0), list(m0), list(m0)
                a_, b_ = abs(r - D1), r + D1
                h2, c2 = (b_ - a_) / 2, (b_ + a_) / 2
                pre = NB * NC * mp.exp(-zB * r)
                for xi, wi in n_in:
                    sc = c2 + h2 * xi
                    tau = (r * r + D1 * D1 - sc * sc) / (2 * r * D1)
                    w = wi * h2 * sc / (r * D1) * pre * mp.exp(-zC * sc)
                    Pl = _legendre_all(lmax, tau)
                    om = 1 - tau * tau
                    for l in range(lmax + 1):
                        p1 = l * ((Pl[l - 1] if l else 0) - tau * Pl[l])          # (1 - tau^2) P_l'
                        p2 = 2 * tau * p1 - l * (l + 1) * om * Pl[l]            # (1 - tau^2)^2 P_l''
                        m0[l] += w * Pl[l]
                        m1[l] += w * tau * Pl[l]
                        m2[l] += w * tau * tau * Pl[l]
                        x1[l] += w * p1
                        x2[l] += w * p2
                kQ0 = [mp.mpf(0)] * (lmax + 1)
                kq, kV0, kD = list(kQ0), list(kQ0), list(kQ0)
                a_, b_ = abs(r - R), r + R
                h2, c2 = (b_ - a_) / 2, (b_ + a_) / 2
                for xi, wi in n_in:
                    s = c2 + h2 * xi
                    t = (r * r + R * R - s * s) / (2 * r * R)
                    w = wi * h2 * s / (r * R)
                    V0, V2 = V_multipole('p', 0, s, zA), V_multipole('p', 2, s, zA)
                    u = (r * t - R) / s
                    P2u = (3 * u * u - 1) / 2
                    P4u = (35 * u ** 4 - 30 * u * u + 3) / 8
                    Q0, q = V0 - V2 / 2, 3 * V2 / (2 * s * s)
                    fd = (V_multipole('d', 0, s, zA) + V_multipole('d', 2, s, zA) * P2u
                          + V_multipole('d', 4, s, zA) * P4u)
                    Pl = _legendre_all(lmax, t)
                    for l in range(lmax + 1):
                        kQ0[l] += w * Q0 * Pl[l]
                        kq[l] += w * q * Pl[l]
                        kV0[l] += w * V0 * Pl[l]
                        kD[l] += w * fd * Pl[l]
                for l in range(lmax + 1):
                    f = mp.mpf(2 * l + 1) / 2
                    rl, t1, t2 = f * m0[l], f * m1[l], f * m2[l]
                    kQ, kqq, kv, kd = f * kQ0[l], f * kq[l], f * kV0[l], f * kD[l]
                    acc['T'][l] += wr * rl * 3 * kv
                    acc['al'][l] += wr * (rl * kQ + r * r * t2 * kqq)
                    acc['be'][l] += wr * (-2 * r * R * t1 * kqq)
                    acc['ga'][l] += wr * (R * R * rl * kqq)
                    acc['X0'][l] += wr * kqq * m0[l]
                    acc['X1'][l] += wr * kqq * r * x1[l]
                    acc['X2'][l] += wr * kqq * r * r * x2[l]
                    acc['DA'][l] += wr * rl * kd
        fp = [4 * mp.pi / (2 * l + 1) for l in range(lmax + 1)]
        for k in ('T', 'al', 'be', 'ga', 'DA'):
            acc[k] = [v * fp[l] for l, v in enumerate(acc[k])]
        Mzz = {}
        for l in range(lmax + 1):
            for part, cvec in (('al', {l: Fraction(1)}), ('be', _times_x({l: Fraction(1)})),
                               ('ga', _times_x(_times_x({l: Fraction(1)})))):
                for L, v in cvec.items():
                    Mzz[L] = Mzz.get(L, 0) + _mf(v) * acc[part][l]
        a = {L: (acc['T'][L] - Mzz[L]) / 2 for L in Ls}
        bf = {}
        for l in range(lmax + 1):
            if l >= 2:
                bf[l] = bf.get(l, 0) - _mf(Fraction(factorial(l - 2), factorial(l + 2))) * acc['X2'][l]
            if l >= 1:
                c = 2 * Fraction(factorial(l - 1), factorial(l + 1))
                dP = {k: Fraction(2 * k + 1) for k in range(l - 1, -1, -2)}          # P_l'
                for L, v in _p2_coeffs(dP).items():
                    bf[L] = bf.get(L, 0) + _mf(c * v) * R * acc['X1'][l]
            for L, v in _p2_coeffs({l: Fraction(1)}).items():
                bf[L] = bf.get(L, 0) - _mf(v) * R * R * acc['X0'][l]
        b = {L: mp.pi * bf.get(L, 0) for L in Ls}
        return dict(a=a, b=b, DA={L: acc['DA'][L] for L in Ls})


# ---------------------------------------------------------------------------------------------------------------
# Toys: Legendre coefficients of (x* - x)^alpha and P_L^2 coefficients of (1 - x^2)(x* - x)^alpha (as in
# debug/aaa_controls.py: Gauss-Legendre projection, P_L^2 part by a double antiderivative)
# ---------------------------------------------------------------------------------------------------------------


def x_star(case):
    D1, D2 = params(case)[:2]
    return (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)


def _anti_list(c):
    g = [mp.mpf(0)] * len(c)
    for L in range(1, len(c) - 1):
        g[L] = c[L - 1] / (2 * L - 1) - c[L + 1] / (2 * L + 3)
    return g


_TOYS = {}


def toy(case, alpha, Lmax, dps=60, deg=7):
    """{'P': Legendre coefficients of (x* - x)^alpha, 'P2': P_L^2 coefficients of (1 - x^2)(x* - x)^alpha}."""
    key = (case, alpha, Lmax, dps, deg)
    if key not in _TOYS:
        with mp.workdps(dps):
            xs = _mf(x_star(case))
            al = _mf(alpha)
            nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(deg, mp.mp.prec)
            n = Lmax + 4
            c = [mp.mpf(0)] * (n + 1)
            for x, w in nodes:
                hv = w * (xs - x) ** al
                p0, p1 = mp.mpf(1), x
                c[0] += hv
                c[1] += hv * x
                for k in range(1, n):
                    p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
                    c[k + 1] += hv * p1
            c = [(2 * k + 1) / mp.mpf(2) * v for k, v in enumerate(c)]
            _TOYS[key] = {'P': c, 'P2': _anti_list(_anti_list(c))}
    return _TOYS[key]


def rho(case):
    xs = _mf(x_star(case))
    return 1 / (xs + mp.sqrt(xs * xs - 1))


def beta_eff(F, case, L):
    """log|G_{L+1}/G_{L-1}| / log((L+1)/(L-1)), G_L = F_L / rho^L."""
    r = rho(case)
    return float(mp.log(abs((F[L + 1] / r ** (L + 1)) / (F[L - 1] / r ** (L - 1)))) / mp.log(mp.mpf(L + 1) / (L - 1)))


# ---------------------------------------------------------------------------------------------------------------
# Shared data for the fast tests: case A1, PY and PX' to L = 81 over Q (one engine, ~40 s), DA to L = 41 (~25 s)
# ---------------------------------------------------------------------------------------------------------------

FAST_L, FAST_L_DA, DPS = 81, 41, 450      # the float conversion loses about 4.5 L digits; dps vs dps + 100 checked


def _fast(case, integral, dps=DPS):
    return series(case, integral, FAST_L_DA if integral == 'DA' else FAST_L, dps)


_DIRECT = {}


def _direct(case, Ls, dps):
    key = (case, Ls, dps)
    if key not in _DIRECT:
        _DIRECT[key] = direct_quadrature(case, Ls, dps)
    return _DIRECT[key]


# ---------------------------------------------------------------------------------------------------------------
# (1) per-channel values against the independent quadrature
# ---------------------------------------------------------------------------------------------------------------


def test_A_potential_multipoles_match_gauss_law():
    """The incomplete-gamma V_l used by the independent route, against Gauss's law integrated directly from the
    density's multipole (p: l = 0, 2; d: l = 0, 2, 4).  Rejects a wrong angular split (1/3, 2/3; 1/5, 2/7, 18/35),
    a wrong radial power or a wrong d normalization, each of which moves V_l by more than 1e-3."""
    with mp.workdps(30):
        z = mp.mpf(9) / 10
        Nd2 = 5 * (2 * z) ** 7 / (4 * mp.pi * 720)
        dens = {('p', 0): lambda t: z ** 5 / mp.pi * t * t * mp.exp(-2 * z * t) / 3,
                ('p', 2): lambda t: z ** 5 / mp.pi * t * t * mp.exp(-2 * z * t) * 2 / 3}
        # d: P2(u)^2 projected on P_l by quadrature in u (not the 1/5, 2/7, 18/35 of the closed form)
        for l in (0, 2, 4):
            w = (2 * l + 1) / mp.mpf(2) * mp.quad(lambda u: ((3 * u * u - 1) / 2) ** 2 * mp.legendre(l, u), [-1, 1])
            dens[('d', l)] = (lambda w: lambda t: Nd2 * w * t ** 4 * mp.exp(-2 * z * t))(w)
        for (kind, l), f in dens.items():
            for s in (mp.mpf('0.3'), mp.mpf('1.7'), mp.mpf('4.0')):
                g = 4 * mp.pi / (2 * l + 1) * (mp.quad(lambda t: t ** (l + 2) * f(t), [0, s]) / s ** (l + 1)
                                               + s ** l * mp.quad(lambda t: t ** (1 - l) * f(t), [s, mp.inf]))
                assert abs(V_multipole(kind, l, s, z) / g - 1) < mp.mpf('1e-25'), (kind, l, s)
        # total charge of the p and d densities is 1 (s V_0 -> 1)
        for kind in ('p', 'd'):
            assert abs(V_multipole(kind, 0, mp.mpf(300), z) * 300 - 1) < mp.mpf('1e-25'), kind


@pytest.mark.parametrize('L', [2, 3])
def test_channels_match_independent_quadrature_A1(L):
    """a_L, b_L of PY and F_L of DA (case A1) against the direct quadrature: relative 1e-20 at 25 digits (measured
    1e-26 to 7e-25).  The quadrature uses the A density's multipole potentials, the lab tensor M_ij and the BC
    density's tau-moments, with no U kernel, no D2-derivative and no jet.  Rejects a wrong angular factor in the
    out-of-plane harmonic (b = -beta/2), a dropped Laplacian or d_z^2 term in a, a wrong U coefficient, and a DA
    with the wrong derivative order (any of these moves a channel by more than 1e-3)."""
    d = _direct('A1', (2, 3), 25)
    py, da = _fast('A1', 'PY'), _fast('A1', 'DA')
    with mp.workdps(40):
        for name, ex in (('a', py['a'][L]), ('b', py['b'][L]), ('DA', da['F'][L])):
            rel = abs(d[name][L] / ex - 1)
            assert rel < mp.mpf('1e-20'), (name, L, mp.nstr(d[name][L], 25), mp.nstr(ex, 25), mp.nstr(rel, 3))


@pytest.mark.slow
def test_channels_match_independent_quadrature_A2():
    """Case A2 (zC != 1, so a misplaced zC power also shows), L = 2, 4, 5: relative 1e-20.  ~60 s."""
    d = direct_quadrature('A2', (2, 4, 5), 25)
    r = exact_Q('A2', 'PY', 5)
    py = {'a': numeric(r[0][:6], 'A2', 80), 'b': numeric(r[1][:6], 'A2', 80)}
    da = numeric(exact_Q('A2', 'DA', 5)[:6], 'A2', 80)
    with mp.workdps(40):
        for L in (2, 4, 5):
            for name, ex in (('a', py['a'][L]), ('b', py['b'][L]), ('DA', da[L])):
                assert abs(d[name][L] / ex - 1) < mp.mpf('1e-20'), (name, L)


# ---------------------------------------------------------------------------------------------------------------
# (2) summed series
# ---------------------------------------------------------------------------------------------------------------

# record AS, debug/data/aaa_p1p3_3600.json: exact split summed to L = 400 at 3600 digits (first 30 digits)
P1_3600 = {('A1', 'PY', '1/2'): '0.170839130315205806095410183122',
           ('A1', 'PXp', '1/2'): '0.174750322791938751287657493742',
           ('A1', 'DA', '1/2'): '0.164452950303524530622586050614',
           ('A1', 'PY', '-1/3'): '0.151564113485113409357883840513',
           ('A1', 'PXp', '-1/3'): '0.153895651805223421335184612953',
           ('A1', 'DA', '-1/3'): '0.160322185058067646851760295738',
           ('A2', 'PY', '0'): '0.187802494860417680960797879255',
           ('A2', 'PXp', '0'): '0.192833622629965559482973089176',
           ('A2', 'DA', '0'): '0.179272588041622512109697848862',
           ('A2', 'PY', '4/5'): '0.212200653562945407138670933421',
           ('A2', 'PXp', '4/5'): '0.215538933793851339925511018231',
           ('A2', 'DA', '4/5'): '0.181466254946353504347909856241'}
# record AS, control R2: AB's even-tempered Gaussian fit, n = 34 (debug/data/aaa_controls.json, 'referee')
R2_34 = {('A1', 'PY', '1/2'): 0.17083913031500453, ('A1', 'PXp', '1/2'): 0.1747503227917389,
         ('A1', 'DA', '1/2'): 0.16445295030332543, ('A1', 'PY', '-1/3'): 0.15156411348497967,
         ('A1', 'PXp', '-1/3'): 0.1538956518050744, ('A1', 'DA', '-1/3'): 0.1603221850583211,
         ('A2', 'PY', '0'): 0.18780249486037548, ('A2', 'PXp', '0'): 0.19283362262992393,
         ('A2', 'DA', '0'): 0.179272588041665, ('A2', 'PY', '4/5'): 0.21220065356289758,
         ('A2', 'PXp', '4/5'): 0.215538933793801, ('A2', 'DA', '4/5'): 0.18146625494643012}


def _summed(S, xq, Lmax):
    x = _mf(Fraction(xq))
    P, P2 = legendre_pair(Lmax, x)
    if 'F' in S:
        return sum(S['F'][L] * P[L] for L in range(Lmax + 1)), None
    v = sum(S['a'][L] * P[L] + S['b'][L] * P2[L] for L in range(Lmax + 1))
    vc = sum(S['c'][L] * P[L] for L in range(Lmax + 1))
    return v, vc


@pytest.mark.parametrize('integral', ['PY', 'PXp'])
def test_summed_series_A1_matches_referee(integral):
    """Sum_{L <= 81} a_L P_L(1/2) + b_L P_L^2(1/2) at A1 against the R2 referee (n = 34) within the record's P1
    tolerance 1e-10 (record: 2.0e-13), and against the record's 3600-digit sum to L = 400 within 1e-22 (the
    truncation at L = 81 leaves ~1e-25).  The total c summed as a P_L series equals the (a, b) sum to 1e-300: the
    re-expansion of P_L^2 in P_L is exact when both run to the truncation c was built at.  Rejects P_L in place of P_L^2 for b, a Condon-
    Shortley sign on P_L^2, a wrong PX' assembly, and a re-expansion with the wrong sign, factor or parity tail."""
    S = _fast('A1', integral)
    with mp.workdps(DPS):
        v, vc = _summed(S, '1/2', FAST_L)
        key = ('A1', integral, '1/2')
        assert abs(v - mp.mpf(R2_34[key])) < 1e-10, mp.nstr(v, 20)
        assert abs(v - mp.mpf(P1_3600[key])) < mp.mpf('1e-22'), mp.nstr(v, 30)
        assert abs(vc - v) < mp.mpf('1e-300'), mp.nstr(abs(vc - v), 3)


def test_summed_series_A1_DA_matches_referee():
    """DA summed to L = 41 at A1, x = 1/2, against R2 (n = 34) within 1e-10 and the 3600-digit value within
    DA_TRUNC = 1e-16 (the truncation at L = 41, measured 8.5e-18).  Rejects a wrong d normalization or a wrong
    U4/U2 coefficient."""
    S = _fast('A1', 'DA')
    with mp.workdps(DPS):
        v, _ = _summed(S, '1/2', FAST_L_DA)
        key = ('A1', 'DA', '1/2')
        assert abs(v - mp.mpf(R2_34[key])) < 1e-10, mp.nstr(v, 20)
        assert abs(v - mp.mpf(P1_3600[key])) < DA_TRUNC, mp.nstr(v, 30)


DA_TRUNC = mp.mpf('1e-16')     # truncation of the DA sum at L = 41, A1, x = 1/2: measured 8.5e-18
SLOW_L = 121


@pytest.mark.slow
def test_summed_series_all_points_slow():
    """All twelve P1 values (both cases, both angles, PY, PX', DA) summed to L = 121 (DA to 81): within 1e-10 of R2
    and 1e-20 of the 3600-digit values; re-expansion exact.  Engine time: the TIMINGS note at the end."""
    for case, xs_ in (('A1', '1/2'), ('A1', '-1/3'), ('A2', '0'), ('A2', '4/5')):
        for integral in ('PY', 'PXp', 'DA'):
            Lm = 81 if integral == 'DA' else SLOW_L
            S = series(case, integral, Lm, 650)
            with mp.workdps(650):
                v, vc = _summed(S, xs_, Lm)
                key = (case, integral, xs_)
                assert abs(v - mp.mpf(R2_34[key])) < 1e-10, (key, mp.nstr(v, 20))
                assert abs(v - mp.mpf(P1_3600[key])) < mp.mpf('1e-20'), (key, mp.nstr(v, 30))
                if vc is not None:
                    assert abs(vc - v) < mp.mpf('1e-400'), key


# ---------------------------------------------------------------------------------------------------------------
# (3) the constant set equals AB's 21
# ---------------------------------------------------------------------------------------------------------------


@pytest.mark.parametrize('case', ['A1', 'A2'])
def test_constant_set_is_ABs_21_no_gamma(case):
    """Mod 2^61 - 1, every L up to 20: a_L of PY and PX' and F_L of DA from L = 2, b_L from L = 3, have exactly the
    constant set of AB at the same case (A1: AB's frozen EXPECTED_A1; A2: AB's own engine, independent of the A
    pair), and no coefficient at any L <= 20 carries Euler's gamma.  Rejects a spurious gamma, a lost or extra
    constant from a wrong kernel or a dropped jump term, and a set that changes with L."""
    K = Field(PRIME)
    ab = AB.EXPECTED_A1 if case == 'A1' else AP._ab_set(case)
    assert len(ab) == 21
    eng = Engine(case, K, 22)
    ay, by = coefficients(case, 20, K, 'PY', eng=eng)
    ax, bx = coefficients(case, 20, K, 'PXp', eng=eng)
    da = coefficients(case, 20, K, 'DA', eng=eng)
    for name, f, L0 in (('PY a', ay, 2), ('PY b', by, 3), ('PXp a', ax, 2), ('PXp b', bx, 3), ('DA', da, 2)):
        for L in range(21):
            assert not [k for k in f[L] if k[1] == 'g'], (name, case, L)
        for L in range(L0, 21):
            S = set(f[L])
            assert S == ab, (name, case, L, sorted(map(repr, S ^ ab)))
    assert by[2] and by[0] == {} and by[1] == {}


# ---------------------------------------------------------------------------------------------------------------
# (4) decay: the 7/2 branch, and 9/2 for the PY total
# ---------------------------------------------------------------------------------------------------------------

ALPHAS = [Fraction(k, 2) for k in (3, 5, 7, 9, 11, 13)]


def window(case, L, alpha, basis):
    """(lo, hi): the midpoints between the toy beta_eff(L) of alpha and those of alpha + 1 and alpha - 1."""
    t = {a: beta_eff(toy(case, a, L + 2)[basis], case, L) for a in ALPHAS}
    i = ALPHAS.index(alpha)
    return (t[ALPHAS[i]] + t[ALPHAS[i + 1]]) / 2, (t[ALPHAS[i]] + t[ALPHAS[i - 1]]) / 2


EXPECT = {('PY', 'a'): (Fraction(7, 2), 'P'), ('PY', 'b'): (Fraction(7, 2), 'P2'), ('PY', 'c'): (Fraction(9, 2), 'P'),
          ('PXp', 'c'): (Fraction(7, 2), 'P'), ('DA', 'F'): (Fraction(7, 2), 'P')}


def _check_decay(case, L, S, which):
    out = {}
    for (integral, n) in which:
        alpha, basis = EXPECT[(integral, n)]
        lo, hi = window(case, L, alpha, basis)
        b = beta_eff(S[integral][n], case, L)
        out[(integral, n)] = (round(b, 3), round(lo, 3), round(hi, 3))
        assert lo < b < hi, (case, L, integral, n, b, (lo, hi))
    return out


def test_decay_branches_A1_L40_and_L60():
    """A1.  beta_eff(40) of PY a, PX' total c' and DA F between the toy midpoints around (x* - x)^{7/2} (toy 7/2:
    -4.616; midpoints -3.99 and -5.28); PY b around the 7/2 toy in the P_L^2 basis (-7.325; midpoints -6.63, -8.05);
    PY total c around the 9/2 toy (-5.937; midpoints -5.28, -6.63).  Measured: a -5.088, c' -5.052, F -4.826,
    b -7.183, c -5.799.  One cutoff past (L = 60, same data): PY a -4.748 and b -6.831 stay in their windows (7/2
    toys -4.420, -6.905).  A PY total on the 7/2 branch (no cancellation; c' and a show what that looks like) sits
    near -5.05 to -5.09, outside the 9/2 window; b dropped from the total makes c = a.  Precision: the a, b and c
    betas at dps and dps + 100 agree to 1e-12."""
    S = {it: _fast('A1', it) for it in ('PY', 'PXp', 'DA')}
    S_hi = _fast('A1', 'PY', DPS + 100)
    for n in ('a', 'b', 'c'):
        assert abs(beta_eff(S['PY'][n], 'A1', 40) - beta_eff(S_hi[n], 'A1', 40)) < 1e-12, n
    _check_decay('A1', 40, S, [('PY', 'a'), ('PY', 'b'), ('PY', 'c'), ('PXp', 'c'), ('DA', 'F')])
    _check_decay('A1', 60, S, [('PY', 'a'), ('PY', 'b')])
    # signs for 20 <= L <= 40 (record, L >= 20: PY a, b negative and c positive; PX' total negative; DA F positive)
    assert all(v < 0 for v in S['PY']['a'][20:41]) and all(v < 0 for v in S['PY']['b'][20:41])
    assert all(v > 0 for v in S['PY']['c'][20:41]) and all(v < 0 for v in S['PXp']['c'][20:41])
    assert all(v > 0 for v in S['DA']['F'][20:42])


@pytest.mark.slow
@pytest.mark.parametrize('case', ['A1', 'A2'])
def test_decay_branches_L60_L80_slow(case):
    """Both cases, L = 60 and 80 (series to L = 121; c needs the b tail beyond L, rho^41 < 3e-7 at A1): every
    scored series in its toy window."""
    S = {it: series(case, it, SLOW_L if it != 'DA' else 81, 650) for it in ('PY', 'PXp', 'DA')}
    for L in (60, 80):
        _check_decay(case, L, S, list(EXPECT))


# ---------------------------------------------------------------------------------------------------------------
# The toys themselves (the windows above and the P4 ratio below rest on them)
# ---------------------------------------------------------------------------------------------------------------

# record AS, control TA (debug/data/aaa_controls.json, 'amp_toys', case A1, 768-node Gauss-Legendre at 120 digits)
TA_220_A1 = {'ta': '5.454419322977413275881174298476923176997e-46',
             'tb': '1.827732382033674341779126726029859156661e-51',
             'tc': '-4.464168549530084167346853526062353167264e-48'}


def test_toys_reproduce_record_amplitude_toys_A1():
    """The in-test toy at L = 220 against the record's control TA (ta: P_L of (x* - x)^{7/2}; tb: P_L^2 of
    (1 - x^2)(x* - x)^{7/2}; tc: P_L of (x* - x)^{9/2}), relative 1e-25; and against the P_L^2 coefficients by
    direct projection with (L-2)!/(L+2)! at L = 5, 9.  Rejects a wrong antiderivative recurrence, a missing
    (2L+1)/2, a Condon-Shortley sign, and a wrong x* (each changes ta or tb at the first digit)."""
    t7 = toy('A1', Fraction(7, 2), 220, dps=100, deg=8)
    t9 = toy('A1', Fraction(9, 2), 220, dps=100, deg=8)
    with mp.workdps(100):
        for name, v in (('ta', t7['P'][220]), ('tb', t7['P2'][220]), ('tc', t9['P'][220])):
            assert abs(v / mp.mpf(TA_220_A1[name]) - 1) < mp.mpf('1e-25'), (name, mp.nstr(v, 20))
    with mp.workdps(40):
        xs = _mf(x_star('A1'))
        for L in (5, 9):
            def f(x):
                P, P2 = legendre_pair(L, x)
                return (1 - x * x) * (xs - x) ** mp.mpf(3.5) * P2[L]
            proj = mp.quad(f, [-1, 1]) * (2 * L + 1) / 2 * mp.factorial(L - 2) / mp.factorial(L + 2)
            assert abs(proj / toy('A1', Fraction(7, 2), 20)['P2'][L] - 1) < mp.mpf('1e-25'), L


# ---------------------------------------------------------------------------------------------------------------
# (5) P4: the leading amplitudes of a and b are in the ratio x*^2 - 1
# ---------------------------------------------------------------------------------------------------------------


def _ratio(case, S, L, tmax):
    """R_L / (x*^2 - 1), R_L = (a_L / ta_L) / (b_L / tb_L), ta / tb the 7/2 toys in the P_L and P_L^2 bases."""
    t = toy(case, Fraction(7, 2), tmax)
    xs = _mf(x_star(case))
    return float((S['a'][L] / t['P'][L]) / (S['b'][L] / t['P2'][L]) / (xs * xs - 1))


def _richardson(Ls, fs):
    """Value at 1/L = 0 of the quadratic in 1/L through three points."""
    h = [mp.mpf(1) / L for L in Ls]
    tot = mp.mpf(0)
    for i in range(3):
        w = mp.mpf(1)
        for j in range(3):
            if j != i:
                w *= (0 - h[j]) / (h[i] - h[j])
        tot += w * fs[i]
    return float(tot)


def test_P4_amplitude_ratio_tends_to_x_star_squared_minus_one_A1():
    """A1.  r_L = R_L / (x*^2 - 1) at L = 40, 60, 80 (measured 1.807, 1.475, 1.335; 1.258 at L = 100).  The leading
    ratio is exactly x*^2 - 1 and the 9/2 terms present in a only add a 1/L excess, so r_L approaches 1 from above:
    asserted (i) R_L > 0 and 1 < r_L < 2.5 at each L (the single-point band, set by these values: it excludes
    R = 0, negative R, and R = 1, which is r = 7.07 here); (ii) r_L decreases with L; (iii) the quadratic-in-1/L
    extrapolation from L = 40, 60, 80 is within 0.05 of 1 (measured 1.017; from 60, 80, 100: 1.006).  Rejects a b
    scaled by 2 (limit 0.5) or by -1, a leading ratio other than x*^2 - 1, and a ratio that does not move with L."""
    S = _fast('A1', 'PY')
    with mp.workdps(DPS):
        r = [_ratio('A1', S, L, FAST_L) for L in (40, 60, 80)]
    assert all(1 < v < 2.5 for v in r), r
    assert r[0] > r[1] > r[2], r
    lim = _richardson((40, 60, 80), r)
    assert abs(lim - 1) < 0.05, (r, lim)


@pytest.mark.slow
@pytest.mark.parametrize('case', ['A1', 'A2'])
def test_P4_amplitude_ratio_slow(case):
    """Both cases, L = 60, 80, 100 (A2 measured 1.373, 1.266, 1.206): r_L in (1, 2), decreasing, and the quadratic
    extrapolation within 0.03 of 1 (measured A1 1.006, A2 1.004; from 40, 60, 80: 1.017, 1.012)."""
    S = series(case, 'PY', SLOW_L, 650)
    with mp.workdps(650):
        r = [_ratio(case, S, L, SLOW_L) for L in (60, 80, 100)]
    assert all(1 < v < 2 for v in r), r
    assert r[0] > r[1] > r[2], r
    assert abs(_richardson((60, 80, 100), r) - 1) < 0.03, r


# TIMINGS (2026-11-15, this machine): fast tests 71 s in all (A1 engine for PY + PX' to L = 81 ~40 s, DA to L = 41
# ~13 s, the independent quadrature ~12 s, constant sets mod p ~7 s, toys ~5 s).  Slow tests 6.8 min in all
# (--slow): the shared engines to L = 121 (PY, PX') and L = 81 (DA) for both cases land in the first slow test
# (~6.2 min); A2 quadrature 13 s; decay 10 s per case; P4 1 s per case.
#
# FIRE-TEST PLANTS (python debug/qa/fire_test.py tests/test_paper59_aa_bc_a_side.py --plant "OLD=>NEW" -k <sel>)
# 2026-11-15, written and fired as a separate step after the guard; 18 of 19 fired (the 19th is a selector note).
#     re-expansion sign          "c[j] = a[j] - j * (j - 1) * b[j]=>c[j] = a[j] + j * (j - 1) * b[j]"      -k "summed_series_A1_matches or decay_branches_A1"  FIRED
#     b dropped from total       "c[j] = a[j] - j * (j - 1) * b[j] + 2 * (2 * j + 1) * tail[j % 2]=>c[j] = a[j]"   -k decay_branches_A1   FIRED (c reads 7/2)
#     parity tail dropped        "+ 2 * (2 * j + 1) * tail[j % 2]=>+ 0 * tail[j % 2]"                      -k summed_series_A1_matches   FIRED
#       ... the same plant with -k decay_branches_A1 alone DID NOT FIRE: without the tail c still decays on the 9/2
#       branch at L = 40 (a different amplitude, the same exponent), so the decay window cannot see it; the exact
#       re-expansion identity in the summed test is the guard for this term.
#     out-of-plane factor        "ja, jb = mzero, jscale(J['beta'], -half)=>ja, jb = mzero, jscale(J['beta'], -half / 2)"   -k channels_match_independent_quadrature_A1; -k P4_amplitude_ratio_tends   FIRED, FIRED
#     out-of-plane sign          "ja, jb = mzero, jscale(J['beta'], -half)=>ja, jb = mzero, jscale(J['beta'], half)"        -k P4_amplitude_ratio_tends   FIRED
#     d_+ lowering sign          "jscale(opB(l, jt), Fraction(-1, 2 * l + 1))=>jscale(opB(l, jt), Fraction(1, 2 * l + 1))"  -k channels_match_independent_quadrature_A1   FIRED
#     Laplacian (trace) dropped  "mzero = jscale(jsum(J['lap'], jscale(J['dzz'], Fraction(-1))), half)=>mzero = jscale(jscale(J['dzz'], Fraction(-1)), half)"   -k channels_match_independent_quadrature_A1   FIRED
#     DA derivative order        "eng.G('U4', fam['U4'], L, 4)[4]=>eng.G('U4', fam['U4'], L, 4)[3]"     -k "channels_match_independent_quadrature_A1 or DA_matches"   FIRED
#     U2 coefficient             "1: Fraction(33, 16) / z=>1: Fraction(33, 15) / z"                       -k channels_match_independent_quadrature_A1   FIRED
#     PX' assembly               "{(L, 2, 0): Fraction(-1)}=>{(L, 2, 0): Fraction(1)}"                    -k summed_series_A1_matches   FIRED
#     P_L for P_L^2              "v = sum(S['a'][L] * P[L] + S['b'][L] * P2[L]=>v = sum(S['a'][L] * P[L] + S['b'][L] * P[L]"   -k summed_series_A1_matches   FIRED
#     spurious Euler gamma       "b.append(jet_value(jb, Gf, eng.D2, K))=>b.append({**jet_value(jb, Gf, eng.D2, K), (Fraction(0), 'g', None): K(1)})"   -k constant_set   FIRED
#     independent-route split    "2: 2 * z ** 5 / mp.pi / 3}=>2: z ** 5 / mp.pi / 3}"                      -k "multipoles or channels_match_independent_quadrature_A1"   FIRED
#     toy antiderivative sign    "g[L] = c[L - 1] / (2 * L - 1) - c[L + 1] / (2 * L + 3)=>g[L] = c[L - 1] / (2 * L - 1) + c[L + 1] / (2 * L + 3)"   -k "toys_reproduce or decay_branches_A1"   FIRED
#     toy x*                     "return (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)=>return (D1 ** 2 + D2 ** 2) / (2 * D1 * D1)"   -k toys_reproduce   FIRED
#     PY total on 7/2            "('PY', 'c'): (Fraction(9, 2), 'P')=>('PY', 'c'): (Fraction(7, 2), 'P')"   -k decay_branches_A1   FIRED
#     referee off 2e-10          "('A1', 'PY', '1/2'): 0.17083913031500453=>('A1', 'PY', '1/2'): 0.17083913051500453"   -k summed_series_A1_matches   FIRED
#   Not a plant: dropping the jump terms (boundary_gen) at r = D2 -- they are identically zero for every derivative
#   used (record AS0, exact at L <= 8), so the plant would be a no-op.
