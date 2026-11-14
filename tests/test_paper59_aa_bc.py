"""Guard for record AB (docs/predictions/2026-11-14_aa_bc_three_centre.md): the three-centre Coulomb integral
(1s_A 1s_A | 1s_B 1s_C) over Slater 1s functions, as a Legendre series in x = cos(angle C B A).

Written 2026-11-14 as a separate step after the computation (CLAUDE.md sec. 9).  Self-contained: the exact
constant-split engine is a port of debug/na_exact.py + debug/aabc_exact.py (no debug/ import).

Claims guarded (the HELD clauses; P3a missed and P4 undecided are NOT guarded):

(1) the reduction (AA|BC) = NA - Y_mu - zeta_A W_mu, mu = 2 zeta_A, per Legendre channel, against an independent
    route: radial quadrature of the density multipole times the Legendre multipole of the 1s potential
    V_A(s) = 1/s - (zeta_A + 1/s) e^{-2 zeta_A s} (itself checked against a radial integral of the density).  The
    quadrature never splits the potential into NA, Y and W.  And the series summed to L = 60 at the referee point
    of case A1 against the record's 2000-digit value and its Gaussian (R2) referee.
(2) the fixed constant set: at L = 2, 3, 7, 20 (mod 2^61 - 1) the nonzero-coefficient constants of F_L are the
    same 21 (case A1; A2 too), inside the record's candidate set C, no EulerGamma, exactly four logs.
(3) the 3/2 cancellation: beta_eff(40) of the combined coefficients is near -4.5 (record: -4.51, -> -4 slowly),
    while NA alone is near -2.14; both in windows that exclude the -2 and -3 alternatives.

Fire-test plants and outcomes: the comment block at the end of this file (kept there so that the first
occurrence of each anchor is the code it mutates).
Results are in docs/claim_test_matrix.md (Paper 59 row for record AB).
"""
from __future__ import annotations

from fractions import Fraction
from fractions import Fraction as F

import flint
import mpmath as mp
import pytest

# ---------------------------------------------------------------------------------------------------------------
# Engine (port of debug/na_exact.py and debug/aabc_exact.py).  F_L / (2 pi N) = sum_c q_{L,c} c over constants
# c = e^{e0} * {1, E1(z), ln v, gamma}, keyed (e0, kind, arg).
# ---------------------------------------------------------------------------------------------------------------


class Field:
    def __init__(self, p=None):
        self.p = p

    def __call__(self, q):
        q = Fraction(q)
        if self.p is None:
            return flint.fmpq(q.numerator, q.denominator)
        return flint.nmod(q.numerator, self.p) / flint.nmod(q.denominator, self.p)

    def pw(self, q, k):
        q = Fraction(q)
        if k == 0:
            return self(1)
        c = self.__dict__.setdefault('_pw', {})
        if q not in c:
            base = self(q)
            c[q] = ([self(1)], [self(1)], base, (self(1) / base) if q != 0 else None)
        pos, neg, base, inv = c[q]
        lst, mult = (pos, base) if k > 0 else (neg, inv)
        while len(lst) <= abs(k):
            lst.append(lst[-1] * mult)
        return lst[abs(k)]

    def ifact(self, n):
        c = self.__dict__.setdefault('_if', [self(1)])
        while len(c) <= n:
            c.append(c[-1] / self(len(c)))
        return c[n]


def bessel_tables(K, Lmax):
    """i[L] = {sigma: {j: coef}} for e^{sigma x} x^{-j};  k[L] = {j: coef} for e^{-x} x^{-j}."""
    h = K(Fraction(1, 2))
    i = [{1: {1: h}, -1: {1: -h}}, {1: {1: h, 2: -h}, -1: {1: h, 2: h}}]
    k = [{1: K(1)}, {1: K(1), 2: K(1)}]
    for L in range(1, Lmax):
        c = K(2 * L + 1)
        new = {}
        for s in (1, -1):
            d = dict(i[L - 1][s])
            for j, v in i[L][s].items():
                d[j + 1] = d.get(j + 1, K(0)) - c * v
            new[s] = d
        i.append(new)
        d = dict(k[L - 1])
        for j, v in k[L].items():
            d[j + 1] = d.get(j + 1, K(0)) + c * v
        k.append(d)
    return i, k


def add(d, key, v):
    d[key] = d[key] + v if key in d else v


class Integrals:
    """J_n(b; A, B) = int_A^B r^n e^{b r} dr; A = 0 finite part, B = None infinity."""

    def __init__(self, K):
        self.K = K
        self.cache = {}

    def endpoint(self, b, R, m):
        K = self.K
        if R is None:
            return {}
        if R == 0:
            if m > 0:
                return {}
            return {(Fraction(0), '1', None): K.pw(b, -m) * K.ifact(-m)}
        return {(Fraction(b) * R, '1', None): K.pw(R, m)}

    def table(self, b, A, B, nmin, nmax):
        key = (b, A, B)
        tab = self.cache.get(key, {})
        K = self.K
        if -1 not in tab:
            if b == 0:
                if A == 0:
                    J = {(Fraction(0), 'ln', Fraction(B)): K(1)}
                else:
                    J = {(Fraction(0), 'ln', Fraction(B)): K(1), (Fraction(0), 'ln', Fraction(A)): K(-1)}
            else:
                if A == 0:
                    J = {(Fraction(0), 'g', None): K(-1), (Fraction(0), 'ln', abs(Fraction(b))): K(-1)}
                else:
                    J = {(Fraction(0), 'E1', -Fraction(b) * A): K(1)}
                if B is not None:
                    add(J, (Fraction(0), 'E1', -Fraction(b) * B), K(-1))
            tab[-1] = J
        n = max(tab) + 1
        while n <= nmax:
            J = {}
            if b == 0:
                for kk, v in self.endpoint(0, B, n + 1).items():
                    add(J, kk, v / K(n + 1))
                for kk, v in self.endpoint(0, A, n + 1).items():
                    add(J, kk, -v / K(n + 1))
            else:
                ib = K(1) / K(Fraction(b))
                for kk, v in self.endpoint(b, B, n).items():
                    add(J, kk, v * ib)
                for kk, v in self.endpoint(b, A, n).items():
                    add(J, kk, -v * ib)
                for kk, v in tab[n - 1].items():
                    add(J, kk, -K(n) * v * ib)
            tab[n] = J
            n += 1
        n = min(tab) - 1
        while n >= nmin:
            J = {}
            inv = K(1) / K(n + 1)
            for kk, v in self.endpoint(b, B, n + 1).items():
                add(J, kk, v * inv)
            for kk, v in self.endpoint(b, A, n + 1).items():
                add(J, kk, -v * inv)
            if b != 0:
                for kk, v in tab[n + 1].items():
                    add(J, kk, -K(Fraction(b)) * v * inv)
            tab[n] = J
            n -= 1
        self.cache[key] = tab
        return tab


def const_value(key):
    e0, kind, arg = key
    f = mp.exp(mp.mpf(e0.numerator) / e0.denominator)
    if kind == '1':
        return f
    if kind == 'g':
        return f * mp.euler
    a = mp.mpf(arg.numerator) / arg.denominator
    if kind == 'ln':
        return f * mp.log(a)
    return f * (mp.e1(a) if a > 0 else -mp.ei(-a))


CASES = {'A1': dict(D1='9/5', D2='13/5', zB='13/10', zC='1', zA='9/10'),
         'A2': dict(D1='3/2', D2='12/5', zB='1', zC='6/5', zA='4/5')}


def params(case):
    g = {k: Fraction(v) for k, v in CASES[case].items()}
    return g['D1'], g['D2'], g['zB'], g['zC'], g['zA']


def _put(d, e0, b, n, v):
    dd = d.setdefault((Fraction(e0), Fraction(b)), {})
    dd[n] = dd[n] + v if n in dd else v


def _ksum(K, tab, arg):
    return sum((K.pw(arg, -j) * v for j, v in tab.items()), K(0))


def density(L, inner, D1, zH, it, kt, K):
    d = {}
    if inner:
        y = zH * D1
        kL, kL1 = _ksum(K, kt[L], y), _ksum(K, kt[L + 1], y)
        for s in (1, -1):
            for j, v in it[L][s].items():
                _put(d, -y, s * zH, -j, K.pw(zH, -j) * v * (K(2 * L + 1) * kL - K(y) * kL1))
            for j, v in it[L + 1][s].items():
                _put(d, -y, s * zH, 1 - j, K.pw(zH, 1 - j) * v * kL)
    else:
        x = zH * D1
        for s in (1, -1):
            iL = sum((K.pw(x, -j) * v for j, v in it[L][s].items()), K(0))
            iL1 = sum((K.pw(x, -j) * v for j, v in it[L + 1][s].items()), K(0))
            for j, v in kt[L].items():
                _put(d, s * x, -zH, -j, K.pw(zH, -j) * v * (K(2 * L + 1) * iL + K(x) * iL1))
            for j, v in kt[L + 1].items():
                _put(d, s * x, -zH, 1 - j, -K.pw(zH, 1 - j) * v * iL)
    return d


def kernel(L, near, D2, mu, which, it, kt, K):
    """NA: r<^L / r>^(L+1);  Y: mu (2L+1) i_L(mu r<) k_L(mu r>);  dY = dY/dmu."""
    d = {}
    c = K(2 * L + 1)
    if which == 'NA':
        if near:
            _put(d, 0, 0, L, K.pw(D2, -L - 1))
        else:
            _put(d, 0, 0, -L - 1, K.pw(D2, L))
        return d
    if near:
        y = mu * D2
        kL, kL1 = _ksum(K, kt[L], y), _ksum(K, kt[L + 1], y)
        for s in (1, -1):
            for j, v in it[L][s].items():
                if which == 'Y':
                    _put(d, -y, s * mu, -j, K(mu) * c * kL * v * K.pw(mu, -j))
                else:
                    _put(d, -y, s * mu, -j, c * (c * kL - K(y) * kL1) * v * K.pw(mu, -j))
            if which == 'dY':
                for j, v in it[L + 1][s].items():
                    _put(d, -y, s * mu, 1 - j, c * kL * v * K.pw(mu, 1 - j))
    else:
        x = mu * D2
        for s in (1, -1):
            iL = sum((K.pw(x, -j) * v for j, v in it[L][s].items()), K(0))
            iL1 = sum((K.pw(x, -j) * v for j, v in it[L + 1][s].items()), K(0))
            for j, v in kt[L].items():
                if which == 'Y':
                    _put(d, s * x, -mu, -j, K(mu) * c * iL * v * K.pw(mu, -j))
                else:
                    _put(d, s * x, -mu, -j, c * (c * iL + K(x) * iL1) * v * K.pw(mu, -j))
            if which == 'dY':
                for j, v in kt[L + 1].items():
                    _put(d, s * x, -mu, 1 - j, -c * iL * v * K.pw(mu, 1 - j))
    return d


def _conv(a, b, K):
    na, nb = max(a), max(b)
    la = [K(0)] * (na - min(a) + 1)
    lb = [K(0)] * (nb - min(b) + 1)
    for n, v in a.items():
        la[na - n] = v
    for n, v in b.items():
        lb[nb - n] = v
    if K.p is None:
        P = flint.fmpq_poly(la) * flint.fmpq_poly(lb)
        co = [P[i] for i in range(P.degree() + 1)]
    else:
        P = flint.nmod_poly([int(v) for v in la], K.p) * flint.nmod_poly([int(v) for v in lb], K.p)
        co = [flint.nmod(int(P[i]), K.p) for i in range(P.degree() + 1)]
    return {na + nb - i: v for i, v in enumerate(co) if v != 0}


def pieces(case, Lmax, K, which=('NA', 'Y', 'dY')):
    D1, D2, zX, zH, zA = params(case)
    mu = 2 * zA  # engine
    it, kt = bessel_tables(K, Lmax + 2)
    I = Integrals(K)
    regions = [(Fraction(0), D1, True, True), (D1, D2, False, True), (D2, None, False, False)]
    out = {w: [] for w in which}
    for L in range(Lmax + 1):
        dens = {inner: density(L, inner, D1, zH, it, kt, K) for inner in (True, False)}
        for w in which:
            res = {}
            for (A, B, inner, near) in regions:
                ker = kernel(L, near, D2, mu, w, it, kt, K)
                for (e0d, bd), pd in dens[inner].items():
                    for (e0k, bk), pk in ker.items():
                        P = _conv(pd, pk, K)
                        bb = bd + bk - zX
                        e0 = e0d + e0k
                        nlo, nhi = min(P) + 2, max(P) + 2
                        tab = I.table(bb, A, B, min(nlo, -1), max(nhi, -1))
                        m2 = K(-2)
                        for n, v in P.items():
                            vv = m2 * v
                            for (e1, kind, arg), u in tab[n + 2].items():
                                add(res, (e0 + e1, kind, arg), vv * u)
            out[w].append({k: v for k, v in res.items() if v != 0})
    return out


def combine(pc, zA, K):
    """f = NA - Y_mu - zA W_mu, with W_mu = -dY/dmu, i.e. f = NA - Y + zA dY, per L."""
    out = []
    za = K(Fraction(zA))
    for qN, qY, qd in zip(pc['NA'], pc['Y'], pc['dY']):
        r = dict(qN)
        for k, v in qY.items():
            add(r, k, -v)  # Y
        for k, v in qd.items():
            add(r, k, za * v)  # W
        out.append({k: v for k, v in r.items() if v != 0})
    return out


def _mpq(v):
    return mp.mpf(int(v.p)) / int(v.q)


def _mf(q):
    return mp.mpf(q.numerator) / q.denominator


def numeric_coeffs(Q, case, dps):
    """F_L (full normalization: 2 sqrt(zB^3 zC^3) times sum q c), evaluated at dps digits."""
    with mp.workdps(dps):
        zB, zC = params(case)[2:4]
        tw = 2 * mp.sqrt(_mf(zB) ** 3 * _mf(zC) ** 3)
        keys = set(k for q in Q for k in q)
        cv = {k: const_value(k) for k in keys}
        return [tw * sum((_mpq(v) * cv[k] for k, v in q.items()), mp.mpf(0)) for q in Q]


_CACHE = {}


def exact_split(case, Lmax):
    """Exact split over Q, L = 0..Lmax (cached; a larger cached split is sliced)."""
    if case not in _CACHE or len(_CACHE[case][1]) <= Lmax:
        K = Field(None)
        pc = pieces(case, max(Lmax, 61), K)
        _CACHE[case] = (pc, combine(pc, params(case)[4], K))
    pc, f = _CACHE[case]
    return {k: v[:Lmax + 1] for k, v in pc.items()}, f[:Lmax + 1]


# ---------------------------------------------------------------------------------------------------------------
# Independent route (no NA / Y / W split, no Bessel kernels, no mu): quadrature
# ---------------------------------------------------------------------------------------------------------------


def potential_1s(s, za):
    """Potential of the normalized density (za^3/pi) e^{-2 za r} at distance s (closed form under test)."""
    return 1 / s - (za + 1 / s) * mp.exp(-2 * za * s)


def potential_1s_quadrature(s, za):
    """Same potential from Gauss's law: (1/s) int_0^s 4 pi r^2 rho + int_s^inf 4 pi r rho."""
    rho = lambda r: za ** 3 / mp.pi * mp.exp(-2 * za * r)
    return (mp.quad(lambda r: 4 * mp.pi * r * r * rho(r), [0, s]) / s
            + mp.quad(lambda r: 4 * mp.pi * r * rho(r), [s, mp.inf]))


def _legendre(L, t):
    p0, p1 = mp.mpf(1), t
    if L == 0:
        return p0
    for l in range(1, L):
        p0, p1 = p1, ((2 * l + 1) * t * p1 - l * p0) / (l + 1)
    return p1


def _multipole(func, r, D, L, nodes):
    """(2L+1)/2 int_{-1}^{1} func(s) P_L(t) dt, s = |r - D e_z|, in the variable s (t = (r^2+D^2-s^2)/(2rD)),
    where func(s) * s is analytic: Gauss-Legendre in s."""
    a, b = abs(r - D), r + D
    h, c = (b - a) / 2, (b + a) / 2
    tot = mp.mpf(0)
    for u, w in nodes:
        s = c + h * u
        t = (r * r + D * D - s * s) / (2 * r * D)
        tot += w * func(s) * s * _legendre(L, t)
    return (2 * L + 1) / mp.mpf(2) * tot * h / (r * D)


def quadrature_FL(case, L, dps=30):
    """F_L of (AA|BC) = 4 pi/(2L+1) int r^2 rho_L(r) v_L(r) dr, rho_L the multipole of phi_B phi_C about B
    (axis BC), v_L the multipole of V_A(|r - A|) about B (axis BA)."""
    with mp.workdps(dps):
        D1, D2, zB, zC, zA = (_mf(v) for v in params(case))
        N = mp.sqrt(zB ** 3 * zC ** 3) / mp.pi
        nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(5, mp.mp.prec)  # 48 nodes
        dens = lambda s: mp.exp(-zC * s)
        pot = lambda s: potential_1s(s, zA)

        def integrand(r):
            rhoL = N * mp.exp(-zB * r) * _multipole(dens, r, D1, L, nodes)
            vL = _multipole(pot, r, D2, L, nodes)
            return 4 * mp.pi / (2 * L + 1) * r * r * rhoL * vL
        return mp.quad(integrand, [0, D1, D2, D2 + 4, mp.inf])


# ---------------------------------------------------------------------------------------------------------------
# (1) the reduction
# ---------------------------------------------------------------------------------------------------------------

REF_X = {'A1': F(1, 2), 'A2': F(-1, 3)}
P1_VALUE = {'A1': '0.21286449876124374093484119906'}    # record AB, P1: exact split, 2000 digits, L <= 240
R2_REFEREE = {'A1': ('0.2128644988', 3e-11)}            # record AB, K2: Gaussian-fit referee, stated band


def test_potential_of_1s_density_matches_gauss_law():
    """Rejects a wrong screened potential (e.g. e^{-za s} or a sign flip on za) used by the quadrature route."""
    with mp.workdps(30):
        for za in (mp.mpf(9) / 10, mp.mpf(4) / 5):
            for s in (mp.mpf('0.3'), mp.mpf('1.7'), mp.mpf('4.1')):
                assert abs(potential_1s(s, za) - potential_1s_quadrature(s, za)) < mp.mpf('1e-25')


def _check_channels(case, Ls, dps, tol):
    _, f = exact_split(case, max(Ls))
    ex = numeric_coeffs(f, case, 80)
    for L in Ls:
        q = quadrature_FL(case, L, dps)
        with mp.workdps(dps + 20):                      # compare above the quadrature's precision
            rel = abs(q / ex[L] - 1)
            assert rel < mp.mpf(tol), (case, L, mp.nstr(q, 25), mp.nstr(ex[L], 25), mp.nstr(rel, 3))


@pytest.mark.parametrize('case,L', [('A1', 7), ('A2', 2)])
def test_reduction_per_channel_matches_independent_quadrature(case, L):
    """F_L of NA - Y_mu - zA W_mu (mu = 2 zA, exact split) equals F_L by quadrature of the density multipole against
    V_A's multipole, relative 1e-20 at 22 digits (measured 3e-23 / 6e-24).  Rejects mu = zA, a flipped W sign, a
    dropped Y (each moves F_L at the 1e-2 level or more)."""
    _check_channels(case, (L,), 22, '1e-20')


@pytest.mark.slow
@pytest.mark.parametrize('case', ['A1', 'A2'])
def test_reduction_per_channel_high_precision(case):
    """Same at L = 0, 3, 7, 30 digits, relative 1e-28 (measured <= 2e-31)."""
    _check_channels(case, (0, 3, 7), 30, '1e-28')


def test_summed_series_matches_referee_A1():
    """Sum_{L <= 60} F_L P_L(1/2) at case A1 vs the record's 2000-digit P1 value (summed to L = 240), tolerance
    1e-15 (F_60 = 2.2e-18, rho = 9/13, so the dropped tail is below 1e-17), and vs the independent R2 Gaussian-fit
    referee 0.2128644988 within its stated +-3e-11 band plus 1e-11.  The split is evaluated at 300 digits (it
    cancels about 100 digits by L = 60; 200 digits leave 10 correct)."""
    Lmax = 60
    _, f = exact_split('A1', Lmax)
    Fv = numeric_coeffs(f, 'A1', 300)
    with mp.workdps(60):
        x = _mf(REF_X['A1'])
        val = sum(Fv[L] * _legendre(L, x) for L in range(Lmax + 1))
        assert abs(Fv[Lmax]) < mp.mpf('1e-17'), mp.nstr(Fv[Lmax], 3)
        assert abs(val - mp.mpf(P1_VALUE['A1'])) < mp.mpf('1e-15'), mp.nstr(val, 25)
        ref, band = R2_REFEREE['A1']
        assert abs(val - mp.mpf(ref)) < band + 1e-11, mp.nstr(val, 15)


# ---------------------------------------------------------------------------------------------------------------
# (2) the fixed constant set
# ---------------------------------------------------------------------------------------------------------------

EXPECTED_A1 = {
    (F(-117, 50), '1', None), (F(-162, 25), '1', None), (F(-189, 50), '1', None), (F(-209, 50), '1', None),
    (F(-389, 50), '1', None), (F(-513, 50), '1', None), (F(-9, 5), '1', None),
    (F(-162, 25), 'E1', F(-27, 10)), (F(-162, 25), 'E1', F(13, 10)), (F(-162, 25), 'E1', F(189, 50)),
    (F(-72, 25), 'E1', F(13, 10)), (F(-72, 25), 'E1', F(369, 50)), (F(-72, 25), 'E1', F(9, 10)),
    (F(-9, 5), 'E1', F(299, 50)), (F(9, 5), 'E1', F(299, 50)),
    (F(162, 25), 'E1', F(533, 50)), (F(72, 25), 'E1', F(533, 50)),
    (F(-162, 25), 'ln', F(1, 2)), (F(-162, 25), 'ln', F(3, 2)), (F(-162, 25), 'ln', F(21, 10)),
    (F(-162, 25), 'ln', F(41, 10)),
}
EXPECTED_A2_LOGS = {F(3, 5), F(7, 5), F(9, 5), F(19, 5)}
PRIME = (1 << 61) - 1


def candidate_set(case, na_keys):
    """The record's candidate set C (P2, negative-power rule), as keys."""
    D1, D2, zB, zC, zA = params(case)
    mu = 2 * zA
    C = set(na_keys)
    P_in = -zC * D1 - mu * D2
    inner = [-zB + s * zC + t * mu for s in (1, -1) for t in (1, -1)]
    shared = [-zB - zC + t * mu for t in (1, -1)]
    for b in inner:
        C |= {(P_in, 'E1', -b * D1), (P_in, 'ln', abs(b)), (P_in, '1', None), (P_in + b * D1, '1', None)}
    for b in shared:
        for s in (1, -1):
            C.add((s * zC * D1 - mu * D2, 'E1', -b * D1))
    for s in (1, -1):
        C.add((s * zC * D1 - mu * D2, 'E1', (zB + zC - mu) * D2))
        for t in (1, -1):
            C.add((s * zC * D1 + t * mu * D2, 'E1', (zB + zC + mu) * D2))
    for s in (1, -1):
        P = s * zC * D1 - mu * D2
        for b in shared:
            C.add((P + b * D1, '1', None))
            C.add((P + b * D2, '1', None))
        for t in (1, -1):
            C.add((s * zC * D1 + t * mu * D2 + (-zB - zC - mu) * D2, '1', None))
    return C


@pytest.mark.parametrize('case', ['A1', 'A2'])
def test_constant_set_is_fixed_no_gamma_four_logs(case):
    """Mod 2^61 - 1, L = 2, 3, 7, 20: no EulerGamma; exactly four logs (the four inner rates |b|); the same 21
    constants at every L, equal to the record's set (A1) and inside the candidate set C.  Rejects a spurious gamma,
    a dropped log, a set that changes with L, mu = zA (moves the mu-carrying arguments)."""
    K = Field(PRIME)
    pc = pieces(case, 20, K)
    f = combine(pc, params(case)[4], K)
    na_keys = set(k for q in pc['NA'] for k in q)
    C = candidate_set(case, na_keys)
    D1, D2, zB, zC, zA = params(case)
    logs_expected = ({F(1, 2), F(3, 2), F(21, 10), F(41, 10)} if case == 'A1' else EXPECTED_A2_LOGS)
    S2 = set(f[2])
    for L in (2, 3, 7, 20):
        S = set(f[L])
        assert not [k for k in S if k[1] == 'g'], (case, L, 'EulerGamma present')
        assert {k[2] for k in S if k[1] == 'ln'} == logs_expected, (case, L)
        assert S == S2, (case, L, S ^ S2)
        assert len(S) == 21, (case, L, len(S))
        assert S <= C, (case, L, S - C)
    if case == 'A1':
        assert S2 == EXPECTED_A1, (S2 ^ EXPECTED_A1)
    # P2d: the mu-free part is NA's own: Y and dY carry no weight on NA's two E1 constants
    na_e1 = [k for k in na_keys if k[1] == 'E1']
    assert len(na_e1) == 2
    for L in range(21):
        assert not any(k in pc['Y'][L] or k in pc['dY'][L] for k in na_e1), (case, L)


# ---------------------------------------------------------------------------------------------------------------
# (3) the 3/2 cancellation: decay L^-4 for f against L^-2 for NA
# ---------------------------------------------------------------------------------------------------------------


def _beta(Fs, case, L):
    D1, D2 = (_mf(v) for v in params(case)[:2])
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    rho = 1 / (xs + mp.sqrt(xs * xs - 1))
    G = lambda n: Fs[n] / rho ** n
    return mp.log(abs(G(L + 1) / G(L - 1))) / mp.log(mp.mpf(L + 1) / (L - 1))


def _decay(case, Lmax, dps):
    """{piece: (F_L at dps digits, F_L at dps + 100 digits)} for f and NA, L = 0..Lmax + 1."""
    pc, f = exact_split(case, Lmax + 1)
    return {name: (numeric_coeffs(Q, case, dps), numeric_coeffs(Q, case, dps + 100))
            for name, Q in (('f', f), ('NA', pc['NA']))}


def test_combination_decays_as_L_minus_four_while_NA_decays_as_L_minus_two():
    """Case A1, exact split, beta_eff(L) = log(G_{L+1}/G_{L-1}) / log((L+1)/(L-1)), G_L = F_L / rho^L, at L = 40
    (record: f -4.51, NA -2.14).  Windows: f in [-4.8, -4.25], NA in [-2.35, -1.95]; a surviving (x*-x)^{5/2}
    branch (L^-3) or (x*-x)^{3/2} branch (L^-2) gives beta_eff(40) far above -4.25 (fire-tested: mu = zA and the
    flipped W sign both leave an L^-3 piece; dropping Y leaves L^-2).  Also F_L > 0 for 5 <= L <= 41.  Precision is
    self-checked: beta at dps and dps + 100 agree to 1e-30."""
    case, L = 'A1', 40
    d = _decay(case, L, 300)
    with mp.workdps(400):
        b = {}
        for name, (lo, hi) in d.items():
            bl, bh = _beta(lo, case, L), _beta(hi, case, L)
            assert abs(bl - bh) < mp.mpf('1e-30'), (name, mp.nstr(bl, 10), mp.nstr(bh, 10))
            b[name] = float(bh)
        assert -4.8 < b['f'] < -4.25, b
        assert -2.35 < b['NA'] < -1.95, b
        fh = d['f'][1]
        assert all(fh[n] > 0 for n in range(5, L + 2)), [n for n in range(5, L + 2) if fh[n] <= 0]


@pytest.mark.slow
@pytest.mark.parametrize('case,window', [('A1', (-4.40, -4.20)), ('A2', (-4.33, -4.13))])
def test_decay_reaches_minus_four_at_L_100(case, window):
    """Both cases at L = 100 (record: A1 -4.226, A2 -4.177; NA -2.062 / -2.044): f within 0.1 of the record and
    below -4.1, NA in [-2.15, -1.95]."""
    L = 100
    d = _decay(case, L, 700)
    with mp.workdps(800):
        bf = float(_beta(d['f'][1], case, L))
        bn = float(_beta(d['NA'][1], case, L))
        assert abs(float(_beta(d['f'][0], case, L)) - bf) < 1e-20
        assert window[0] < bf < window[1], bf
        assert -2.15 < bn < -1.95, bn


# FIRE-TEST PLANTS (python debug/qa/fire_test.py tests/test_paper59_aa_bc.py --plant "OLD=>NEW" -k <test>)
#     mu = zeta_A           "mu = 2 * zA  # engine=>mu = 1 * zA  # engine"
#     W sign flipped        "add(r, k, za * v)  # W=>add(r, k, -za * v)  # W"  (f = NA - Y - zA dY)
#     Y dropped             "add(r, k, -v)  # Y=>add(r, k, 0 * v)  # Y"
#     spurious gamma        "r = dict(qN)=>r = dict(qN); r[(Fraction(0), 'g', None)] = K(1)"
#     log constant dropped  "(Fraction(0), 'ln', abs(Fraction(b))): K(-1)}=>}"  (engine finite part)
#     one expected log gone "(F(-162, 25), 'ln', F(41, 10)),=>"
