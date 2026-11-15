"""Guard for record AP (docs/predictions/2026-11-14_aa_bc_p_orbitals.md): the three-centre integrals

    PC = (1s_A 1s_A | 1s_B 2pz_C)      and      PA = (2ps_A 2ps_A | 1s_B 1s_C)

over Slater functions, as Legendre series in x = cos(angle C B A).  Paper 59 sec:scope.

Written 2026-11-14 as a separate step after the computation (CLAUDE.md sec. 9).  The base engine (Field, Bessel
tables, the exponential integrals, AB's 1s density) is imported from tests/test_paper59_aa_bc.py; the AP engine
(debug/aapc_exact.py) is ported here, with both of its kernel paths (symbolic monomials; factored by Leibniz).

Claims guarded (record AP, all clauses held):

(1) values: per Legendre channel, the exact split equals an independent quadrature that never uses the PA
    reduction (PA = d^2/dD2^2 int rho U + int rho W0) or any D2-derivative: PC integrates the multipole of
    (z - D1) e^{-zC r_C} against the 1s potential; PA integrates AB's 1s.1s multipole against the full 2p-sigma
    potential V0(s) + V2(s) P2(n.w/s), whose closed forms are checked against Gauss's law.  And the series summed to
    L = 60 at A1's referee point against the record's 2000-digit values and R2 bands.
(2) the constant set: mod 2^61 - 1, at L = 1, 2, 3, 7, 20 the nonzero constants of F_L are exactly AB's 21 at the
    same case, for PC and PA, A1 and A2; no EulerGamma.
(3) the 7/2 branch: beta_eff(40) of PC and PA in windows that exclude the 9/2 branch (beta near -5.8 at L = 40, the
    W0 piece's own value) and the 5/2 and 3/2 branches; PA's second-derivative piece and W0 piece separately.
(4) the two kernel paths give identical coefficients (L <= 12, Q and mod p): the factored path that scored P2 is
    the AP0-checked symbolic path.

Record AQ (docs/predictions/2026-11-14_aa_bc_ppi.md), PX = (1s_A 1s_A | 1s_B 2px_C), p_x in the ABC plane:
(5) per-channel G_L against a direct-projection quadrature (no Legendre identity); the series in P_L^1 summed at
    two angles of A1 against the record; AB's 21 constants at L = 1..20; beta_eff(40) near -5.8 (L^-5), windows
    that exclude reading the branch as a P_L series and a cancelled 7/2 coefficient.

Record AD (docs/predictions/2026-11-15_aa_bc_d_orbital.md), PD = (1s_A 1s_A | 1s_B 3dsigma_C):
(6) per-channel F_L against a direct projection of h2(r - C) as is; the summed series at two angles of A1 (this is
    what tests the d factor; the decay does not distinguish it); AB's 21 constants; L_max invariance;
    beta_eff(40) windows excluding a cancelled 7/2 coefficient and a 5/2 branch.

Fire-test plants and outcomes: the comment block at the end of this file (AP 17/17; AQ 9/9 plus one plant that is
a no-op at A1 and fires at A2, 2026-11-14; AD 17/17, 2026-11-15; each a separate step).
"""
from __future__ import annotations

import importlib.util
from fractions import Fraction
from fractions import Fraction as F
from math import comb
from pathlib import Path

import mpmath as mp
import pytest

_spec = importlib.util.spec_from_file_location('_ab_guard', Path(__file__).with_name('test_paper59_aa_bc.py'))
AB = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(AB)

Field, params, add, _put, _conv = AB.Field, AB.params, AB.add, AB._put, AB._conv
PRIME = (1 << 61) - 1


# ---------------------------------------------------------------------------------------------------------------
# AP engine (port of debug/aapc_exact.py)
# ---------------------------------------------------------------------------------------------------------------


def _addm(d, key, v):
    d[key] = d[key] + v if key in d else v


def kernel_mono(L, near, which, it, kt, K):
    """{(pa, pd, n, t0, s): c}: c mu^pa D2^pd r^n e^{t0 mu D2} e^{s mu r}, NA or Y kernel."""
    d = {}
    if which == 'NA':
        _addm(d, (0, -L - 1, L, 0, 0) if near else (0, L, -L - 1, 0, 0), K(1))
        return d
    c = K(2 * L + 1)
    if near:
        for i, vk in kt[L].items():
            for s in (1, -1):
                for j, vi in it[L][s].items():
                    _addm(d, (1 - i - j, -i, -j, -1, s), c * vk * vi)
    else:
        for sg in (1, -1):
            for i, vi in it[L][sg].items():
                for j, vk in kt[L].items():
                    _addm(d, (1 - i - j, -i, -j, sg, -1), c * vi * vk)
    return d


def d_mu(m, K):
    out = {}
    for (pa, pd, n, t0, s), c in m.items():
        if pa:
            _addm(out, (pa - 1, pd, n, t0, s), K(pa) * c)
        if t0:
            _addm(out, (pa, pd + 1, n, t0, s), K(t0) * c)
        if s:
            _addm(out, (pa, pd, n + 1, t0, s), K(s) * c)
    return {k: v for k, v in out.items() if v != 0}


def d_D2(m, K):
    out = {}
    for (pa, pd, n, t0, s), c in m.items():
        if pd:
            _addm(out, (pa, pd - 1, n, t0, s), K(pd) * c)
        if t0:
            _addm(out, (pa + 1, pd, n, t0, s), K(t0) * c)
    return {k: v for k, v in out.items() if v != 0}


def evaluate(m, mu, D2, K):
    acc = {}
    for (pa, pd, n, t0, s), c in m.items():
        _addm(acc, (t0, s, n), c * K.pw(mu, pa) * K.pw(D2, pd))
    d = {}
    for (t0, s, n), v in acc.items():
        _put(d, t0 * mu * D2, s * mu, n, v)
    return d


ORDER = {'Y': 0, 'dY': 1, 'd2Y': 2, 'd3Y': 3}


def kernel_symbolic(L, near, coeffs, it, kt, K, mu, D2, nD2):
    tot = {}
    for name, a in coeffs.items():
        if name == 'NA':
            m = kernel_mono(L, near, 'NA', it, kt, K)
        else:
            m = kernel_mono(L, near, 'Y', it, kt, K)
            for _ in range(ORDER[name]):
                m = d_mu(m, K)
        for key, v in m.items():
            _addm(tot, key, a * v)
    for _ in range(nD2):
        tot = d_D2(tot, K)
    return evaluate({k: v for k, v in tot.items() if v != 0}, mu, D2, K)


def _factors(L, near, it, kt, K):
    G, H = {}, {}
    c = K(2 * L + 1)
    if near:
        for i, vk in kt[L].items():
            _addm(G, (1 - i, -i, -1), c * vk)
        for s in (1, -1):
            for j, vi in it[L][s].items():
                _addm(H, (-j, -j, s), vi)
    else:
        for sg in (1, -1):
            for i, vi in it[L][sg].items():
                _addm(G, (1 - i, -i, sg), c * vi)
        for j, vk in kt[L].items():
            _addm(H, (-j, -j, -1), vk)
    return G, H


def _dG(G, K, var):
    out = {}
    for (pa, pd, t0), c in G.items():
        if var == 'mu':
            if pa:
                _addm(out, (pa - 1, pd, t0), K(pa) * c)
            if t0:
                _addm(out, (pa, pd + 1, t0), K(t0) * c)
        else:
            if pd:
                _addm(out, (pa, pd - 1, t0), K(pd) * c)
            if t0:
                _addm(out, (pa + 1, pd, t0), K(t0) * c)
    return out


def _dH(H, K):
    out = {}
    for (pa, n, s), c in H.items():
        if pa:
            _addm(out, (pa - 1, n, s), K(pa) * c)
        if s:
            _addm(out, (pa, n + 1, s), K(s) * c)
    return out


def kernel_factored(L, near, coeffs, it, kt, K, mu, D2, nD2):
    d = {}
    if 'NA' in coeffs:
        na = {k: coeffs['NA'] * v for k, v in kernel_mono(L, near, 'NA', it, kt, K).items()}
        for _ in range(nD2):
            na = d_D2(na, K)
        for (e0, b), poly in evaluate(na, mu, D2, K).items():
            for n, c in poly.items():
                _put(d, e0, b, n, c)
    amax = max(ORDER[k] for k in coeffs if k != 'NA')
    G0, H0 = _factors(L, near, it, kt, K)
    Gs, Hs = [G0], [H0]
    for _ in range(amax):
        Gs.append(_dG(Gs[-1], K, 'mu'))
        Hs.append(_dH(Hs[-1], K))
    for _ in range(nD2):
        Gs = [_dG(g, K, 'D2') for g in Gs]
    Gv, Hv = [], []
    for g in Gs:
        v = {}
        for (pa, pd, t0), c in g.items():
            _addm(v, t0, c * K.pw(mu, pa) * K.pw(D2, pd))
        Gv.append(v)
    for h in Hs:
        v = {}
        for (pa, n, s), c in h.items():
            _addm(v, (s, n), c * K.pw(mu, pa))
        Hv.append(v)
    acc = {}
    for name, a in coeffs.items():
        if name == 'NA':
            continue
        o = ORDER[name]
        for k in range(o + 1):
            w = a * K(comb(o, k))  # Leibniz
            for t0, gv in Gv[k].items():
                for (s, n), hv in Hv[o - k].items():
                    _addm(acc, (t0, s, n), w * gv * hv)
    for (t0, s, n), v in acc.items():
        _put(d, t0 * mu * D2, s * mu, n, v)
    return d


def dens_pz(L, inner, D1, zC, it, kt, K):
    """Multipole L of (z - D1) e^{-zC r_C}: [r (L d_{L-1} + (L+1) d_{L+1}) - (2L+1) D1 d_L] / (2L+1)."""
    d = {}
    inv = K(1) / K(2 * L + 1)
    parts = []
    if L >= 1:
        parts.append((K(L) * inv, 1, AB.density(L - 1, inner, D1, zC, it, kt, K)))
    parts.append((K(L + 1) * inv, 1, AB.density(L + 1, inner, D1, zC, it, kt, K)))
    parts.append((-K(Fraction(D1)), 0, AB.density(L, inner, D1, zC, it, kt, K)))  # node term
    for w, sh, dd in parts:
        for (e0, b), poly in dd.items():
            for n, v in poly.items():
                _put(d, e0, b, n + sh, w * v)
    return d


def dens_px(L, inner, D1, zC, it, kt, K):
    """Record AQ: m = 1 multipole of x e^{-zC r_C}, r (dens_{L-1} - dens_{L+1}) / (2L+1); none at L = 0."""
    d = {}
    if L == 0:
        return d
    inv = K(1) / K(2 * L + 1)
    for w, dd in ((inv, AB.density(L - 1, inner, D1, zC, it, kt, K)),
                  (-inv, AB.density(L + 1, inner, D1, zC, it, kt, K))):  # m = 1 sign
        for (e0, b), poly in dd.items():
            for n, v in poly.items():
                _put(d, e0, b, n + 1, w * v)
    return d


def _dacc(out, dd, w, shift):
    for (e0, b), poly in dd.items():
        for n, v in poly.items():
            _put(out, e0, b, n + shift, w * v)


def _tmul(E, L, K):
    """(T E)_L = [L E_{L-1} + (L+1) E_{L+1}] / (2L+1): multiply by t in AB's normalization."""
    out = {}
    inv = K(1) / K(2 * L + 1)
    if L >= 1:
        _dacc(out, E(L - 1), K(L) * inv, 0)
    _dacc(out, E(L + 1), K(L + 1) * inv, 0)
    return out


def dens_dz2(L, inner, D1, zC, it, kt, K):
    """Record AD: multipole of h2(r - C) e^{-zC r_C}, h2 = (3 r^2 t^2 - r^2 - 4 D1 r t + 2 D1^2)/2 about B."""
    cache = {}

    def E(k):
        if k not in cache:
            cache[k] = AB.density(k, inner, D1, zC, it, kt, K)
        return cache[k]

    def TE(k):
        return _tmul(E, k, K)
    out = {}
    D1q = K(Fraction(D1))
    half = K(Fraction(1, 2))
    _dacc(out, _tmul(TE, L, K), K(3) * half, 2)   # 3 r^2 T^2
    _dacc(out, E(L), -half, 2)                     # - r^2
    _dacc(out, TE(L), K(-2) * D1q, 1)              # - 4 D1 r T
    _dacc(out, E(L), D1q * D1q, 0)                 # + 2 D1^2
    return out


def integrate(dens, kers, D1, D2, zB, I, K):
    regions = [(Fraction(0), D1, True, True), (D1, D2, False, True), (D2, None, False, False)]
    res = {}
    for (A, B, inner, near) in regions:
        for (e0d, bd), pd in dens[inner].items():
            for (e0k, bk), pk in kers[near].items():
                P = _conv(pd, pk, K)
                if not P:
                    continue
                bb = bd + bk - zB
                tab = I.table(bb, A, B, min(min(P) + 2, -1), max(max(P) + 2, -1))
                for n, v in P.items():
                    for (e1, kind, arg), u in tab[n + 2].items():
                        add(res, (e0d + e0k + e1, kind, arg), K(-2) * v * u)
    return res


def boundary_D2(dens_outer, dK_near, dK_far, D2, zB, K):
    """Jump term of d^2/dD2^2 at r = D2."""
    jump = {}
    for src, sg in ((dK_near, 1), (dK_far, -1)):
        for (e0, b), poly in src.items():
            for n, v in poly.items():
                _put(jump, e0, b, n, K(sg) * v)
    res = {}
    for (e0d, bd), pd in dens_outer.items():
        for (e0k, bk), pk in jump.items():
            P = _conv(pd, pk, K)
            bb = bd + bk - zB
            for n, v in P.items():
                add(res, (e0d + e0k + bb * D2, '1', None), K(-2) * v * K.pw(D2, n + 2))
    return res


def coefficients(case, Lmax, K, integral, path='factored'):
    """[{key: q}] per L, F_L = 2 sqrt(zB^3 zC^m) sum q c; for PA also the pieces [(U2 + B), W0]."""
    D1, D2, zB, zC, zA = params(case)
    mu = 2 * zA
    it, kt = AB.bessel_tables(K, Lmax + 3)
    I = AB.Integrals(K)
    za = K(Fraction(zA))
    kev = kernel_factored if path == 'factored' else kernel_symbolic
    out, pieces = [], []
    for L in range(Lmax + 1):
        if integral in ('PC', 'PX', 'PD'):
            dfun = {'PC': dens_pz, 'PX': dens_px, 'PD': dens_dz2}[integral]
            dens = {inner: dfun(L, inner, D1, zC, it, kt, K) for inner in (True, False)}
            co = {'NA': K(1), 'Y': K(-1), 'dY': za}
            res = integrate(dens, {nr: kev(L, nr, co, it, kt, K, mu, D2, 0) for nr in (True, False)},
                            D1, D2, zB, I, K)
        else:
            dens = {inner: AB.density(L, inner, D1, zC, it, kt, K) for inner in (True, False)}
            cU = {'NA': K(3) / (K(2) * za * za), 'Y': -K(3) / (K(2) * za * za), 'dY': K(33) / (K(16) * za),
                  'd2Y': K(Fraction(-9, 8)), 'd3Y': za / K(4)}
            cW = {'NA': K(1), 'Y': K(-1), 'dY': K(5) * za / K(4), 'd2Y': -za * za / K(2)}
            kU2 = {nr: kev(L, nr, cU, it, kt, K, mu, D2, 2) for nr in (True, False)}
            kU1 = {nr: kev(L, nr, cU, it, kt, K, mu, D2, 1) for nr in (True, False)}
            kW = {nr: kev(L, nr, cW, it, kt, K, mu, D2, 0) for nr in (True, False)}
            pU = integrate(dens, kU2, D1, D2, zB, I, K)
            for k, v in boundary_D2(dens[False], kU1[True], kU1[False], D2, zB, K).items():
                add(pU, k, v)
            pW = integrate(dens, kW, D1, D2, zB, I, K)
            res = dict(pU)
            for k, v in pW.items():
                add(res, k, v)
            pieces.append(({k: v for k, v in pU.items() if v != 0}, {k: v for k, v in pW.items() if v != 0}))
        out.append({k: v for k, v in res.items() if v != 0})
    return (out, pieces) if integral == 'PA' else out


def _mf(q):
    return mp.mpf(q.numerator) / q.denominator


def numeric(Q, case, integral, dps):
    with mp.workdps(dps):
        zB, zC = params(case)[2:4]
        tw = 2 * mp.sqrt(_mf(zB) ** 3 * _mf(zC) ** (3 if integral == 'PA' else 5))  # PX normalization
        if integral == 'PD':   # 2 pi N_B N_d
            tw = 2 * mp.pi * mp.sqrt(_mf(zB) ** 3 / mp.pi) * mp.sqrt(5 * (2 * _mf(zC)) ** 7 / (4 * mp.pi * 720))
        keys = set(k for q in Q for k in q)
        cv = {k: AB.const_value(k) for k in keys}
        return [tw * sum((mp.mpf(int(v.p)) / int(v.q) * cv[k] for k, v in q.items()), mp.mpf(0)) for q in Q]


_CACHE = {}


def exact(case, integral, Lmax):
    key = (case, integral)
    if key not in _CACHE or len(_CACHE[key][0]) <= Lmax:
        r = coefficients(case, max(Lmax, 61), Field(None), integral)
        _CACHE[key] = r if integral == 'PA' else (r, None)
    f, pc = _CACHE[key]
    return f[:Lmax + 1], (pc[:Lmax + 1] if pc else None)


# ---------------------------------------------------------------------------------------------------------------
# Independent route: quadrature (no reduction, no D2-derivatives, no Bessel kernels)
# ---------------------------------------------------------------------------------------------------------------


def V0_2ps(s, z):
    """Monopole potential of the normalized 2p-sigma density (z^5/pi) r^2 cos^2 e^{-2 z r}."""
    return 1 / s - (1 / s + 3 * z / 2 + z * z * s + z ** 3 * s * s / 3) * mp.exp(-2 * z * s)


def V2_2ps(s, z):
    """Quadrupole radial potential (coefficient of P2(cos)); evaluated with guard digits (cancels near s = 0)."""
    with mp.workdps(mp.mp.dps + 150):
        s = mp.mpf(s)
        v = (3 / (z * z * s ** 3) - (2 * z ** 3 * s * s / 3 + 2 * z * z * s + 4 * z + 6 / s + 6 / (z * s * s)
                                    + 3 / (z * z * s ** 3)) * mp.exp(-2 * z * s))
    return +v


def V_l_gauss(l, s, z):
    """Multipole l of the 2p-sigma potential from the density by Gauss's law (independent of the closed forms)."""
    w = mp.mpf(1) / 3 if l == 0 else mp.mpf(2) / 3
    rho = lambda r: z ** 5 / mp.pi * w * r * r * mp.exp(-2 * z * r)
    return 4 * mp.pi / (2 * l + 1) * (mp.quad(lambda r: rho(r) * r ** (l + 2), [0, s]) / s ** (l + 1)
                                      + s ** l * mp.quad(lambda r: rho(r) * r ** (1 - l), [s, mp.inf]))


def _multipole(func, r, D, L, nodes):
    """(2L+1)/2 int P_L(t) func(s, t) dt, s = |r - D e_z| (integrated in s; func(s, t) s smooth)."""
    a, b = abs(r - D), r + D
    h, c = (b - a) / 2, (b + a) / 2
    tot = mp.mpf(0)
    for u, w in nodes:
        s = c + h * u
        t = (r * r + D * D - s * s) / (2 * r * D)
        tot += w * func(s, t) * s * AB._legendre(L, t)
    return (2 * L + 1) / mp.mpf(2) * tot * h / (r * D)


def quadrature_FL(case, integral, L, dps):
    with mp.workdps(dps):
        D1, D2, zB, zC, zA = (_mf(v) for v in params(case))
        m = 5 if integral == 'PC' else 3
        N = mp.sqrt(zB ** 3 * zC ** m) / mp.pi
        nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(5, mp.mp.prec)
        if integral == 'PC':
            dens = lambda s, t, r: (r * t - D1) * mp.exp(-zC * s)
            pot = lambda s, t, r: AB.potential_1s(s, zA)
        else:
            dens = lambda s, t, r: mp.exp(-zC * s)

            def pot(s, t, r):
                c = (r * t - D2) / s            # cosine of w = r - A against n = A/|A|
                return V0_2ps(s, zA) + V2_2ps(s, zA) * (3 * c * c - 1) / 2

        def integrand(r):
            rhoL = N * mp.exp(-zB * r) * _multipole(lambda s, t: dens(s, t, r), r, D1, L, nodes)
            vL = _multipole(lambda s, t: pot(s, t, r), r, D2, L, nodes)
            return 4 * mp.pi / (2 * L + 1) * r * r * rhoL * vL
        return mp.quad(integrand, [0, D1, D2, D2 + 4, mp.inf])


# ---------------------------------------------------------------------------------------------------------------
# (1) values
# ---------------------------------------------------------------------------------------------------------------


def test_2ps_potential_closed_forms_match_gauss_law():
    """V0 and V2 of the 2p-sigma density against Gauss's law; rejects a wrong polynomial term or the 1/3 vs 2/3
    angular split (each moves V by 1e-3 or more)."""
    with mp.workdps(30):
        for z in (mp.mpf(9) / 10, mp.mpf(4) / 5):
            for s in (mp.mpf('0.4'), mp.mpf('1.7'), mp.mpf('4.1')):
                assert abs(V0_2ps(s, z) - V_l_gauss(0, s, z)) < mp.mpf('1e-25')
                assert abs(V2_2ps(s, z) - V_l_gauss(2, s, z)) < mp.mpf('1e-25')


@pytest.mark.parametrize('integral,case,L', [('PC', 'A1', 3), ('PC', 'A2', 1), ('PA', 'A1', 2), ('PA', 'A2', 5)])
def test_exact_split_per_channel_matches_independent_quadrature(integral, case, L):
    """F_L of the exact split equals the quadrature of the density multipole against the potential's multipole
    (PA: the full 2p-sigma potential, no reduction), relative 1e-18 at 22 digits."""
    f, _ = exact(case, integral, L)
    ex = numeric(f, case, integral, 80)[L]
    q = quadrature_FL(case, integral, L, 22)
    with mp.workdps(40):
        rel = abs(q / ex - 1)
        assert rel < mp.mpf('1e-18'), (integral, case, L, mp.nstr(q, 22), mp.nstr(ex, 22), mp.nstr(rel, 3))


@pytest.mark.slow
@pytest.mark.parametrize('integral', ['PC', 'PA'])
def test_exact_split_per_channel_high_precision(integral):
    """Both cases, L = 0, 4, 8, at 30 digits, relative 1e-26."""
    for case in ('A1', 'A2'):
        f, _ = exact(case, integral, 8)
        ex = numeric(f, case, integral, 80)
        for L in (0, 4, 8):
            q = quadrature_FL(case, integral, L, 30)
            with mp.workdps(50):
                assert abs(q / ex[L] - 1) < mp.mpf('1e-26'), (integral, case, L)


P1_VALUE = {'PC': '-0.208962996502063071758676195484', 'PA': '0.203231256960263779073627583616'}  # record AP, A1
R2_REFEREE = {'PC': ('-0.20896299649', 1.5e-10), 'PA': ('0.20323125696', 1e-10)}               # record AP, K2


@pytest.mark.parametrize('integral', ['PC', 'PA'])
def test_summed_series_matches_record_A1(integral):
    """Sum_{L <= 60} F_L P_L(1/2), case A1, against the record's 2000-digit value (summed to L = 240) to 1e-14 and
    the R2 Gaussian-fit referee within its band.  Split evaluated at 300 digits."""
    f, _ = exact('A1', integral, 60)
    Fv = numeric(f, 'A1', integral, 300)
    with mp.workdps(60):
        x = mp.mpf(1) / 2
        val = sum(Fv[L] * AB._legendre(L, x) for L in range(61))
        assert abs(Fv[60]) < mp.mpf('1e-16'), mp.nstr(Fv[60], 3)
        assert abs(val - mp.mpf(P1_VALUE[integral])) < mp.mpf('1e-14'), mp.nstr(val, 25)
        ref, band = R2_REFEREE[integral]
        assert abs(val - mp.mpf(ref)) < band, mp.nstr(val, 15)


# ---------------------------------------------------------------------------------------------------------------
# (2) the constant set equals AB's 21
# ---------------------------------------------------------------------------------------------------------------


def _ab_set(case):
    K = Field(PRIME)
    pc = AB.pieces(case, 2, K)
    return set(AB.combine(pc, params(case)[4], K)[2])


@pytest.mark.parametrize('case', ['A1', 'A2'])
@pytest.mark.parametrize('integral', ['PC', 'PA'])
def test_constant_set_is_ABs_21_no_gamma(integral, case):
    """Mod 2^61 - 1, L = 1, 2, 3, 7, 20: no EulerGamma; the set equals AB's 21 at the same case (A1: the frozen
    EXPECTED_A1 of AB's guard).  Rejects a dropped node term (removes keys at some L), a spurious gamma, a lost
    boundary term or a wrong U coefficient that leaves a new key or loses one."""
    K = Field(PRIME)
    r = coefficients(case, 20, K, integral)
    f = r[0] if integral == 'PA' else r
    ab = AB.EXPECTED_A1 if case == 'A1' else _ab_set(case)
    assert len(ab) == 21
    for L in (1, 2, 3, 7, 20):
        S = set(f[L])
        assert not [k for k in S if k[1] == 'g'], (integral, case, L)
        assert S == ab, (integral, case, L, sorted(map(repr, S ^ ab)))


# ---------------------------------------------------------------------------------------------------------------
# (3) the 7/2 branch
# ---------------------------------------------------------------------------------------------------------------


def _beta(Fs, case, L):
    return AB._beta(Fs, case, L)


def test_decay_is_the_seven_halves_branch_A1():
    """Case A1, beta_eff(40) (record: PC -4.636, PA -4.386; PA pieces: second-derivative -3.831, W0 -5.798).
    Windows: PC [-4.85, -4.4], PA [-4.6, -4.15], second-derivative piece [-4.05, -3.6], W0 [-6.05, -5.55].  A
    surviving 3/2 or 5/2 branch puts beta_eff(40) near -2 or -3; the 9/2 branch alone sits where W0 does.  Also one
    sign for 5 <= L <= 41.  Precision self-checked (dps vs dps + 100 agree to 1e-30)."""
    case, L = 'A1', 40
    out = {}
    for integral in ('PC', 'PA'):
        f, pc = exact(case, integral, L + 1)
        series = [('f', f)]
        if pc:
            series += [('U', [p[0] for p in pc]), ('W', [p[1] for p in pc])]
        for name, Q in series:
            lo, hi = numeric(Q, case, integral, 300), numeric(Q, case, integral, 400)
            with mp.workdps(400):
                bl, bh = _beta(lo, case, L), _beta(hi, case, L)
                assert abs(bl - bh) < mp.mpf('1e-30'), (integral, name)
                out[(integral, name)] = float(bh)
            if name == 'f':
                assert all(hi[n] > 0 for n in range(5, L + 2)), (integral, [n for n in range(5, L + 2) if hi[n] <= 0])
    assert -4.85 < out[('PC', 'f')] < -4.4, out
    assert -4.6 < out[('PA', 'f')] < -4.15, out
    assert -4.05 < out[('PA', 'U')] < -3.6, out
    assert -6.05 < out[('PA', 'W')] < -5.55, out


def test_jump_term_vanishes_because_U_has_no_coulomb_tail():
    """The jump term of d^2/dD2^2 at r = D2 is identically zero (L <= 10, both cases, Q): U = 3(1 - e)/(2 zA^2 s)
    + ... has no 1/s singularity, so its kernel's first D2-derivative is continuous at r = D2.  Rejects a U whose
    NA and Y coefficients do not cancel (a Coulomb tail), which makes the jump nonzero."""
    K = Field(None)
    for case in ('A1', 'A2'):
        D1, D2, zB, zC, zA = params(case)
        mu, za = 2 * zA, K(Fraction(zA))
        it, kt = AB.bessel_tables(K, 13)
        cU = {'NA': K(3) / (K(2) * za * za), 'Y': -K(3) / (K(2) * za * za), 'dY': K(33) / (K(16) * za),
              'd2Y': K(Fraction(-9, 8)), 'd3Y': za / K(4)}  # jump-test U
        for L in range(11):
            dens = AB.density(L, False, D1, zC, it, kt, K)
            k1 = {nr: kernel_factored(L, nr, cU, it, kt, K, mu, D2, 1) for nr in (True, False)}
            B = boundary_D2(dens, k1[True], k1[False], D2, zB, K)
            assert not [v for v in B.values() if v != 0], (case, L)


# ---------------------------------------------------------------------------------------------------------------
# (4) the two kernel paths agree
# ---------------------------------------------------------------------------------------------------------------


@pytest.mark.parametrize('p', [None, PRIME])
def test_factored_kernel_path_equals_symbolic_path(p):
    """L <= 12, both cases, PC and PA: identical coefficient dicts.  Rejects a wrong Leibniz weight or a D2-
    derivative applied to the r factor."""
    K = Field(p)
    for case in ('A1', 'A2'):
        for integral in ('PC', 'PA'):
            a = coefficients(case, 12, K, integral, 'factored')
            b = coefficients(case, 12, K, integral, 'symbolic')
            assert a == b, (case, integral)


# ---------------------------------------------------------------------------------------------------------------
# Record AQ: (1s_A 1s_A | 1s_B 2px_C), p_x in the ABC plane, PX = sum_{L >= 1} G_L P_L^1(x)
# ---------------------------------------------------------------------------------------------------------------


def _PL_dPL(L, t):
    p0, p1 = mp.mpf(1), t
    if L == 0:
        return p0, mp.mpf(0)
    for n in range(1, L):
        p0, p1 = p1, ((2 * n + 1) * t * p1 - n * p0) / (n + 1)
    return p1, L * (t * p1 - p0) / (t * t - 1)


def _PL1_list(Lmax, x):
    """P_L^1(x) = sqrt(1 - x^2) P_L'(x), no Condon-Shortley phase."""
    return [mp.mpf(0)] + [mp.sqrt(1 - x * x) * _PL_dPL(L, x)[1] for L in range(1, Lmax + 1)]


def quadrature_GL_px(case, L, dps):
    """G_L by direct projection: rho^1_L = (2L+1)/(2L(L+1)) int N e^{-zB r} r (1-t^2) P_L'(t) e^{-zC s} dt (no
    Legendre identity), against v_L of the exact 1s potential, 4 pi/(2L+1) int r^2 rho^1_L v_L dr."""
    with mp.workdps(dps):
        D1, D2, zB, zC, zA = (_mf(v) for v in params(case))
        N = mp.sqrt(zB ** 3 * zC ** 5) / mp.pi
        nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(5, mp.mp.prec)
        fac = mp.mpf(2 * L + 1) / (2 * L * (L + 1))

        def ang(func, r, D):
            a, b = abs(r - D), r + D
            h, c = (b - a) / 2, (b + a) / 2
            tot = mp.mpf(0)
            for u, w in nodes:
                s_ = c + h * u
                t = (r * r + D * D - s_ * s_) / (2 * r * D)
                tot += w * func(s_, t) * s_
            return tot * h / (r * D)

        def integrand(r):
            rho1 = fac * N * mp.exp(-zB * r) * r * ang(
                lambda s_, t: mp.exp(-zC * s_) * (1 - t * t) * _PL_dPL(L, t)[1], r, D1)
            vL = (2 * L + 1) / mp.mpf(2) * ang(lambda s_, t: AB.potential_1s(s_, zA) * _PL_dPL(L, t)[0], r, D2)
            return 4 * mp.pi / (2 * L + 1) * r * r * rho1 * vL
        return mp.quad(integrand, [0, D1, D2, D2 + 4, mp.inf])


@pytest.mark.parametrize('case,L', [('A1', 3), ('A2', 1)])
def test_px_channel_matches_direct_projection(case, L):
    """G_L of the exact split equals the direct-projection quadrature, relative 1e-18 at 22 digits (record AQ0:
    <= 6e-23).  Rejects a lost m = 1 sign (d_{L-1} + d_{L+1}), a shifted recurrence, the 2px normalization."""
    f, _ = exact(case, 'PX', L)
    ex = numeric(f, case, 'PX', 80)[L]
    q = quadrature_GL_px(case, L, 22)
    with mp.workdps(40):
        assert abs(q / ex - 1) < mp.mpf('1e-18'), (case, L, mp.nstr(q, 22), mp.nstr(ex, 22))


@pytest.mark.slow
def test_px_channel_high_precision():
    for case in ('A1', 'A2'):
        f, _ = exact(case, 'PX', 7)
        ex = numeric(f, case, 'PX', 80)
        for L in (2, 7):
            q = quadrature_GL_px(case, L, 30)
            with mp.workdps(50):
                assert abs(q / ex[L] - 1) < mp.mpf('1e-26'), (case, L)


# record AQ: 2000-digit values (L <= 240) and the R2 Gaussian-fit referee
PX_VALUES = {F(1, 2): ('0.0459315847412031954161933793045', '0.04593158474'),
             F(0): ('0.0413929924096152832160027409986', '0.04139299240')}


@pytest.mark.parametrize('xq', [F(1, 2), F(0)])
def test_px_summed_series_at_two_angles_A1(xq):
    """Sum_{L <= 60} G_L P_L^1(x) at A1, x = 1/2 and 0 (two angles catch a reduction error smooth in L), against the
    record's 2000-digit values to 1e-14 and the R2 referee within 2e-10.  Rejects P_L in place of P_L^1."""
    f, _ = exact('A1', 'PX', 60)
    G = numeric(f, 'A1', 'PX', 300)
    with mp.workdps(60):
        x = _mf(xq)
        P1 = _PL1_list(60, x)
        val = sum(G[L] * P1[L] for L in range(1, 61))  # series in P_L^1
        ex, r2 = PX_VALUES[xq]
        assert abs(val - mp.mpf(ex)) < mp.mpf('1e-14'), mp.nstr(val, 25)
        assert abs(val - mp.mpf(r2)) < 2e-10


@pytest.mark.parametrize('case', ['A1', 'A2'])
def test_px_constant_set_is_ABs_21(case):
    """Mod 2^61 - 1, L = 1, 2, 3, 7, 20: no gamma; exactly AB's 21 at the same case; nothing at L = 0."""
    K = Field(PRIME)
    f = coefficients(case, 20, K, 'PX')
    ab = AB.EXPECTED_A1 if case == 'A1' else _ab_set(case)
    assert f[0] == {}
    for L in (1, 2, 3, 7, 20):
        S = set(f[L])
        assert not [k for k in S if k[1] == 'g'], (case, L)
        assert S == ab, (case, L, sorted(map(repr, S ^ ab)))


def test_px_decay_is_L_minus_five():
    """beta_eff(40) of G_L (record: A1 -5.810, A2 -5.658; toy with the bare 7/2 branch: -5.937, -5.766).  Windows
    A1 [-6.0, -5.6], A2 [-5.85, -5.45].  The same branch read as a P_L series would sit near -4.6; a 9/2 branch
    (the 7/2 coefficient cancelled) near -7.  One sign for 5 <= L <= 41; dps vs dps + 100 agree to 1e-30."""
    L = 40
    for case, (lo_w, hi_w) in (('A1', (-6.0, -5.6)), ('A2', (-5.85, -5.45))):
        f, _ = exact(case, 'PX', L + 1)
        lo, hi = numeric(f, case, 'PX', 300), numeric(f, case, 'PX', 400)
        with mp.workdps(400):
            bl, bh = _beta(lo, case, L), _beta(hi, case, L)
            assert abs(bl - bh) < mp.mpf('1e-30')
        assert lo_w < float(bh) < hi_w, (case, float(bh))
        assert all(hi[n] > 0 for n in range(5, L + 2)), case


# ---------------------------------------------------------------------------------------------------------------
# Record AD: (1s_A 1s_A | 1s_B 3dsigma_C), 3ds = N_d h2(r - C) e^{-zC rC}, PD = sum_L F_L P_L(x)
# ---------------------------------------------------------------------------------------------------------------


def quadrature_FL_d(case, L, dps):
    """F_L of PD with the density multipole projected directly: h2(r - C) = (3 (r t - D1)^2 - s^2)/2 as is (no T
    recurrence, no polynomial rearrangement), against v_L of the exact 1s potential."""
    with mp.workdps(dps):
        D1, D2, zB, zC, zA = (_mf(v) for v in params(case))
        NB = mp.sqrt(zB ** 3 / mp.pi)
        Nd = mp.sqrt(5 * (2 * zC) ** 7 / (4 * mp.pi * 720))
        nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(5, mp.mp.prec)

        def ang(func, r, D):
            a, b = abs(r - D), r + D
            h, c = (b - a) / 2, (b + a) / 2
            tot = mp.mpf(0)
            for u, w in nodes:
                s_ = c + h * u
                t = (r * r + D * D - s_ * s_) / (2 * r * D)
                tot += w * func(s_, t) * s_
            return tot * h / (r * D)

        def integrand(r):
            rho = (2 * L + 1) / mp.mpf(2) * NB * Nd * mp.exp(-zB * r) * ang(
                lambda s_, t: mp.exp(-zC * s_) * (3 * (r * t - D1) ** 2 - s_ * s_) / 2 * _PL_dPL(L, t)[0], r, D1)
            vL = (2 * L + 1) / mp.mpf(2) * ang(lambda s_, t: AB.potential_1s(s_, zA) * _PL_dPL(L, t)[0], r, D2)
            return 4 * mp.pi / (2 * L + 1) * r * r * rho * vL
        return mp.quad(integrand, [0, D1, D2, D2 + 4, mp.inf])


@pytest.mark.parametrize('case,L', [('A1', 0), ('A1', 4), ('A2', 2)])
def test_d_channel_matches_direct_projection(case, L):
    """F_L of the exact split equals the direct-projection quadrature, relative 1e-18 at 22 digits (record AD0:
    <= 1.1e-22).  Rejects a wrong h2 coefficient, a single T where T^2 belongs, the 3d normalization."""
    f, _ = exact(case, 'PD', L)
    ex = numeric(f, case, 'PD', 80)[L]
    q = quadrature_FL_d(case, L, 22)
    with mp.workdps(40):
        assert abs(q / ex - 1) < mp.mpf('1e-18'), (case, L, mp.nstr(q, 22), mp.nstr(ex, 22))


PD_VALUES = {F(1, 2): ('0.104757560074952207333361729138', '0.10475756007475293'),   # record AD: 2000 digits,
             F(0): ('0.119232388104143576937469926637', '0.1192323881040522')}     # L <= 240; R2 n = 34


@pytest.mark.parametrize('xq', [F(1, 2), F(0)])
def test_d_summed_series_at_two_angles_A1(xq):
    """Sum_{L <= 60} F_L P_L(x) at A1, x = 1/2 and 0, against the record's 2000-digit values to 1e-14 and the n = 34
    Gaussian referee within 1e-10.  The values are what tests the d factor (the decay does not)."""
    f, _ = exact('A1', 'PD', 60)
    Fv = numeric(f, 'A1', 'PD', 300)
    with mp.workdps(60):
        x = _mf(xq)
        val = sum(Fv[L] * AB._legendre(L, x) for L in range(61))
        ex, r2 = PD_VALUES[xq]
        assert abs(val - mp.mpf(ex)) < mp.mpf('1e-14'), mp.nstr(val, 25)
        assert abs(val - mp.mpf(r2)) < 1e-10


@pytest.mark.parametrize('case', ['A1', 'A2'])
def test_d_constant_set_is_ABs_21_and_Lmax_invariant(case):
    """Mod 2^61 - 1: no gamma; exactly AB's 21 at L = 1, 2, 3, 7, 20; and F_L for L <= 25 unchanged when the internal
    L_max goes 30 -> 40 (a truncated top-of-range T^2 recurrence would break it)."""
    K = Field(PRIME)
    f = coefficients(case, 20, K, 'PD')
    ab = AB.EXPECTED_A1 if case == 'A1' else _ab_set(case)
    for L in (1, 2, 3, 7, 20):
        S = set(f[L])
        assert not [k for k in S if k[1] == 'g'], (case, L)
        assert S == ab, (case, L, sorted(map(repr, S ^ ab)))
    assert coefficients(case, 30, K, 'PD')[:26] == coefficients(case, 40, K, 'PD')[:26]


def test_d_decay_seven_halves_branch():
    """beta_eff(40) of F_L (record: A1 -4.931, A2 -4.775; toy with the bare 7/2 branch -4.616, -4.500).  Windows A1
    [-5.1, -4.75], A2 [-4.95, -4.6]: a cancelled 7/2 coefficient (9/2 branch) reads near -5.8 at L = 40, a 5/2
    branch near -3.6.  One sign for 5 <= L <= 41; dps vs dps + 100 agree to 1e-30."""
    L = 40
    for case, (lo_w, hi_w) in (('A1', (-5.1, -4.75)), ('A2', (-4.95, -4.6))):
        f, _ = exact(case, 'PD', L + 1)
        lo, hi = numeric(f, case, 'PD', 300), numeric(f, case, 'PD', 400)
        with mp.workdps(400):
            bl, bh = _beta(lo, case, L), _beta(hi, case, L)
            assert abs(bl - bh) < mp.mpf('1e-30')
        assert lo_w < float(bh) < hi_w, (case, float(bh))
        assert all(hi[n] > 0 for n in range(5, L + 2)), case


# FIRE-TEST PLANTS (python debug/qa/fire_test.py tests/test_paper59_aa_bc_p.py --plant "OLD=>NEW" -k <test>)
#     node term dropped      "(-K(Fraction(D1)), 0,=>(-K(Fraction(0)), 0,"
#     recurrence weight      "parts.append((K(L + 1) * inv, 1,=>parts.append((K(L) * inv, 1,"
#     Coulomb tail in U      "za / K(4)}  # jump-test U=>za / K(4), 'Y': -K(1)}  # jump-test U"
#   Not plants: dropping the boundary term changes nothing (it is identically zero, guarded above); dropping the
#   node term changes coefficients but not the constant set (caught by the quadrature test, not the set test).
#     U coefficient          "K(33) / (K(16) * za)=>K(33) / (K(15) * za)"
#     W0 sign                "'d2Y': -za * za / K(2)}=>'d2Y': za * za / K(2)}"
#     one D2 derivative      "kU2 = {nr: kev(L, nr, cU, it, kt, K, mu, D2, 2)=>kU2 = {nr: kev(L, nr, cU, it, kt, K, mu, D2, 1)"
#     Leibniz weight         "w = a * K(comb(o, k))  # Leibniz=>w = a * K(1)  # Leibniz"
#     spurious gamma         "out.append({k: v for k, v in res.items() if v != 0})=>res[(Fraction(0), 'g', None)] = K(1); out.append({k: v for k, v in res.items() if v != 0})"
#     V2 closed form         "+ 6 / (z * s * s)=>+ 5 / (z * s * s)"
#   record AQ (PX):
#     m = 1 sign lost        "(-inv, AB.density(L + 1, inner, D1, zC, it, kt, K))):  # m = 1 sign=>(inv, AB.density(L + 1, inner, D1, zC, it, kt, K))):  # m = 1 sign"
#     (the normalization plant is a no-op at A1, where zC = 1; it fires at A2)
#     PX normalization       "** (3 if integral == 'PA' else 5))  # PX normalization=>** (3 if integral in ('PA', 'PX') else 5))  # PX normalization"
#     P_L for P_L^1          "val = sum(G[L] * P1[L] for L in range(1, 61))  # series in P_L^1=>val = sum(G[L] * AB._legendre(L, x) for L in range(1, 61))  # series in P_L^1"
#     shifted recurrence     "(inv, AB.density(L - 1, inner, D1, zC, it, kt, K)),=>(inv, AB.density(L, inner, D1, zC, it, kt, K)),"
#   record AD (PD):
#     h2 coefficient         "_dacc(out, _tmul(TE, L, K), K(3) * half, 2)   # 3 r^2 T^2=>_dacc(out, _tmul(TE, L, K), K(2) * half, 2)   # 3 r^2 T^2"
#     T once for T^2         "_dacc(out, _tmul(TE, L, K), K(3) * half, 2)   # 3 r^2 T^2=>_dacc(out, TE(L), K(3) * half, 2)   # 3 r^2 T^2"
#     linear-term sign       "_dacc(out, TE(L), K(-2) * D1q, 1)              # - 4 D1 r T=>_dacc(out, TE(L), K(2) * D1q, 1)              # - 4 D1 r T"
#     constant term dropped  "_dacc(out, E(L), D1q * D1q, 0)                 # + 2 D1^2=>_dacc(out, E(L), 0 * D1q, 0)                 # + 2 D1^2"
#     3d normalization       "(4 * mp.pi * 720))=>(4 * mp.pi * 600))"
#     L_max-dependent truncation (two plants, applied together; real newlines):
#                            "    mu = 2 * zA\n    it, kt = AB.bessel_tables(K, Lmax + 3)=>    mu = 2 * zA\n    globals()['_TOP'] = Lmax - 8\n    it, kt = AB.bessel_tables(K, Lmax + 3)"
#                            "    _dacc(out, E(L + 1), K(L + 1) * inv, 0)\n    return out=>    if L + 1 <= _TOP:\n        _dacc(out, E(L + 1), K(L + 1) * inv, 0)\n    return out"
#                            (a fixed cut such as `if L < 35` is not a plant here: it never triggers at L <= 32)
