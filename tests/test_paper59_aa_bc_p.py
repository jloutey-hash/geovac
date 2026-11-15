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

Fire-test plants and outcomes: the comment block at the end of this file.
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
        if integral == 'PC':
            dens = {inner: dens_pz(L, inner, D1, zC, it, kt, K) for inner in (True, False)}
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
        tw = 2 * mp.sqrt(_mf(zB) ** 3 * _mf(zC) ** (5 if integral == 'PC' else 3))
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
