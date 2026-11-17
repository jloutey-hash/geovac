"""Exact Legendre coefficients of the four-centre 1s integral (XY|ZW) on FC's slice (record FX engine; computes no
scored quantity).  Geometry, densities and normalisation as debug/fc_engine.py:

    F_L = (4 pi/(2L+1))^2 int r^2 rho2_L(r) psi_L(r) dr,   psi_L(r) = int s^2 rho1_L(s) r<^L / r>^(L+1) ds,

and this module returns  G_L = F_L / pi^2 = sum_c q_{L,c} c  with q exact (fmpq) or mod p (nmod), c constant keys
of debug/fx_reducer.py.

rho1_L (about X): record NA's exact density (aabc_exact.density, Bessel tables of na_exact).
rho2_L (about n):
  * matched slice (zW = sqrt(a/b) zZ, debug/fx_engine.py): the u-representation, built from the bivariate polynomial
    Q_l(u, r) = (2ra)^l P_l(mu(u, r)) by the Legendre three-term rule (one step per L, so L <= Lmax costs Lmax steps):
        rho2_L = (2L+1)/(2a)^(L+1) r^-(L+1) [ sum_{u = u+, -u-}  +- ( e^{-zZ u}(S + U)(u)  -  C E1(zZ u) ) ],
    S, U, C the antiderivative data of e^{-zZ u} h(u) (h = Q_L (u - g^2 u^-3)/4) from the same ladders as the reducer;
  * a = b (two centres at one point; FC's K4 reduction to record AB): rho2 = e^{-mu |r - a n|}, mu = zZ + zW, whose
    Legendre coefficient is -(2L+1) x [record NA's density with zX = 0, zH = mu, D1 = a] (record AB's dY kernel).
psi_L is assembled per panel from exact antiderivatives Phi = int s^L g, Gamma = int s^-(L+1) g (g = s^2 rho1_L):
e^{bs} x Laurent + C E1(-b s), constants at D1 and the finite part at 0.  Every panel product r^2 rho2_L psi_L is
handed to the reducer; the r-constant pieces of psi are integrated once per panel under a tag atom and multiplied out.

API:  CASES;  Engine(case, K).run(Lmax) -> generator of (L, {CK: q});  generic sequences via generic_sequence().
"""
from __future__ import annotations

import sys
from fractions import Fraction as Fr
from pathlib import Path

import flint

sys.path.insert(0, str(Path(__file__).parent))
import na_exact as NE          # noqa: E402
import aabc_exact as AB        # noqa: E402
import fx_reducer as X         # noqa: E402

CASES = {
    'M1': dict(a='1', b='4', D1='6/5', zX='13/10', zY='1', zZ='6/5', zW='3/5'),
    # replacement for FX design M2 (memo section 5): a = 4/9, b = 1, D1 = 6/5, kappa = 2/3; zZ, zW kept from M2
    'M2r': dict(a='4/9', b='1', D1='6/5', zX='13/10', zY='1', zZ='3/2', zW='1'),
    'M2': dict(a='4/5', b='9/5', D1='6/5', zX='13/10', zY='1', zZ='3/2', zW='1'),
    # a = b: record AB case A1 (X = B, Y = C, Z = W = A at D2 = 13/5, zZ + zW = 2 zA = 9/5), as FC's K4A1
    'K4A1': dict(a='13/5', b='13/5', D1='9/5', zX='13/10', zY='1', zZ='1', zW='4/5'),
}
TAG_LO, TAG_HI = ('~tag', 'lo'), ('~tag', 'hi')
ONE = (Fr(0), ())


def params(case):
    d = CASES[case] if isinstance(case, str) else case
    p = {k: Fr(v) for k, v in d.items()}
    if p['a'] != p['b']:
        k2 = p['a'] / p['b']
        from math import isqrt
        n, m = isqrt(k2.numerator), isqrt(k2.denominator)
        assert n * n == k2.numerator and m * m == k2.denominator, 'sqrt(a/b) must be rational'
        p['kappa'] = Fr(n, m)
        assert p['zW'] == p['kappa'] * p['zZ'], 'not on the matched line'
    return p


def _laurent_from(R, d):
    """{n: coef} -> (offset, poly) meaning r^offset * poly(r)."""
    d = {n: v for n, v in d.items() if X.nz(v)}
    if not d:
        return 0, R.poly([])
    lo, hi = min(d), max(d)
    return lo, R.poly([d.get(n, R.K(0)) for n in range(lo, hi + 1)])


class Engine:
    def __init__(self, case, K):
        self.p = params(case)
        self.K, self.R = K, X.Ring(K)
        p = self.p
        self.matched = p['a'] != p['b']
        br = sorted({Fr(0), p['D1'], p['a'], p['b']})
        self.panels = list(zip(br, br[1:] + [None]))
        self.stats = {}

    # ------------------------------------------------------------------ rho1 side: psi_L on a panel
    def _anti(self, lam, lau):
        """Indefinite antiderivative of e^{lam s} r^off P(s): (e-Laurent dict, C) with  F = e^{lam s} Lau(s) - C E1(-lam s)."""
        K, R = self.K, self.R
        off, P = lau
        cs = R.coeffs(P)
        pos = {off + i: v for i, v in enumerate(cs) if off + i >= 0}
        neg = {-(off + i): v for i, v in enumerate(cs) if off + i < 0}
        red = X.Reducer(R, 0, 0)
        out = {}
        if pos:
            S = red._poly_anti(lam, R.poly([pos.get(n, K(0)) for n in range(max(pos) + 1)]))
            for n, v in enumerate(R.coeffs(S)):
                X.add(out, n, v)
        C = K(0)
        if neg:
            c = [K(0)] + [neg.get(k, K(0)) for k in range(1, max(neg) + 1)]
            U, C = red._pp_anti(lam, c)
            for k in range(1, len(U)):
                X.add(out, -k, U[k])
        return out, C

    def _rho1_parts(self, L, inner, it, kt):
        """g s^L and g s^-(L+1) (g = s^2 rho1_L) as {(e0, rate): {n: coef}}."""
        K, p = self.K, self.p
        dens = AB.density(L, inner, p['D1'], p['zY'], it, kt, K)
        m = K(-(2 * L + 1))
        hi, lo = {}, {}
        for (e0, b), d in dens.items():
            key = (Fr(e0), Fr(b) - p['zX'])
            hi[key] = {n + 2 + L: m * v for n, v in d.items()}
            lo[key] = {n + 1 - L: m * v for n, v in d.items()}
        return hi, lo

    def _antis(self, parts):
        """[(e0, rate, laurent dict, C)] for an antiderivative of a parts dict."""
        out = []
        for (e0, b), d in parts.items():
            assert b != 0
            lau, C = self._anti(b, _laurent_from(self.R, d))
            out.append((e0, b, lau, C))
        return out

    def _at(self, anti, R0):
        """Antiderivative value at R0 > 0 (constants) or finite part at 0; {CK: coef}."""
        K = self.K
        res = {}
        for e0, b, lau, C in anti:
            if R0 == 0:
                v, f = K(0), Fr(1)
                for n, c in lau.items():
                    if n <= 0:
                        k = -n
                        f = Fr(b) ** k / _fact(k)
                        v += c * K(f)
                if X.nz(v):
                    X.add(res, (e0, ()), v)
                if X.nz(C):                             # -C E1(-b s) -> -C(-gamma - ln(-b))
                    X.add(res, (e0, (('g',),)), C)
                    X.add(res, (e0, (('ln', -b),)), C)
            else:
                x = K(R0)
                v = sum((c * (x ** n if n >= 0 else (K(1) / x) ** (-n)) for n, c in lau.items()), K(0))
                if X.nz(v):
                    X.add(res, (e0 + b * R0, ()), v)
                if X.nz(C):
                    X.add(res, (e0, (('E1', -b * R0),)), -C)
        return res

    def psi(self, L, panel, it, kt):
        """psi_L on a panel: (function terms [(CK, rate, (off, poly), E)], const_lo {CK: q}, const_hi {CK: q})."""
        K, R, p = self.K, self.R, self.p
        A, B = panel
        inner = B is not None and B <= p['D1']
        hi_in, lo_in = self._rho1_parts(L, True, it, kt)
        PhiI, GamI = self._antis(hi_in), self._antis(lo_in)
        hi_out, lo_out = self._rho1_parts(L, False, it, kt)
        PhiO, GamO = self._antis(hi_out), self._antis(lo_out)
        neg = lambda d: {k: -v for k, v in d.items()}    # noqa: E731

        def merge(*ds):
            o = {}
            for d in ds:
                for k, v in d.items():
                    X.add(o, k, v)
            return o
        if inner:
            phi_f, phi_s, gam_f, gam_s = PhiI, 1, GamI, -1
            c_lo = neg(self._at(PhiI, Fr(0)))
            c_hi = merge(self._at(GamI, p['D1']), neg(self._at(GamO, p['D1'])))
        else:
            phi_f, phi_s, gam_f, gam_s = PhiO, 1, GamO, -1
            c_lo = merge(self._at(PhiI, p['D1']), neg(self._at(PhiI, Fr(0))), neg(self._at(PhiO, p['D1'])))
            c_hi = {}
        terms = []
        for anti, sgn, shift in ((phi_f, phi_s, -L - 1), (gam_f, gam_s, L)):
            for e0, b, lau, C in anti:
                d = {n + shift: K(sgn) * v for n, v in lau.items()}
                terms.append(((e0, ()), b, _laurent_from(R, d), ()))
                if X.nz(C):                            # sgn * (-C) E1(-b r) r^shift
                    terms.append(((e0, ()), Fr(0), _laurent_from(R, {shift: -K(sgn) * C}), ((-b, Fr(0)),)))
        terms.append(((Fr(0), (TAG_LO,)), Fr(0), _laurent_from(R, {-L - 1: K(1)}), ()))
        terms.append(((Fr(0), (TAG_HI,)), Fr(0), _laurent_from(R, {L: K(1)}), ()))
        return terms, {k: v for k, v in c_lo.items() if X.nz(v)}, {k: v for k, v in c_hi.items() if X.nz(v)}

    # ------------------------------------------------------------------ rho2 side: f = r^2 rho2_L on a panel
    def _Q_steps(self, Lmax):
        """Yields Q_L (dict u-power -> r-poly) for L = 0..Lmax."""
        K, R, p = self.K, self.R, self.p
        a, b = p['a'], p['b']
        g = R.poly([K(-(b - a) * a), K(0), K((b - a) / b)])     # (b - a)(r^2 - ab)/b
        nu = {2: R.poly([K(Fr(-1, 4))]),
              0: R.poly([K(a * a), K(0), K(1)]) - g * K(Fr(1, 2)),
              -2: -(g * g) * K(Fr(1, 4))}
        r2a2 = R.poly([K(0), K(0), K(4 * a * a)])
        Qm, Q = None, {0: R.poly([K(1)])}
        yield Q
        Qm, Q = Q, dict(nu)
        yield Q
        for l in range(1, Lmax):
            new = {}
            c1, c2 = K(Fr(2 * l + 1, l + 1)), K(Fr(l, l + 1))
            for n, P in Q.items():
                for s, v in nu.items():
                    X.add(new, n + s, P * v * c1)
            for n, P in Qm.items():
                X.add(new, n, -(P * r2a2) * c2)
            Qm, Q = Q, {n: P for n, P in new.items() if P.degree() >= 0}
            yield Q

    def f_terms_matched(self, L, Q, panel):
        """[(e0, lam, N, I, pole, Kp, E)]: f = r^2 rho2_L = sum e^{e0 + lam r} N(r)/(r^I (r - pole)^Kp) E(r)."""
        K, R, p = self.K, self.R, self.p
        a, b, kap, z = p['a'], p['b'], p['kappa'], p['zZ']
        g = R.poly([K(-(b - a) * a), K(0), K((b - a) / b)])
        g2 = g * g
        h = {}
        for n, P in Q.items():
            X.add(h, n + 1, P * K(Fr(1, 4)))
            X.add(h, n - 3, -(P * g2) * K(Fr(1, 4)))
        h = {n: P for n, P in h.items() if P.degree() >= 0}
        lam = K(-z)
        il = K(1) / lam
        Np = max(h)
        S = {Np: h[Np] * il}
        for k in range(Np - 1, -1, -1):
            S[k] = (h.get(k, R.poly([])) - S[k + 1] * K(k + 1)) * il
        M = -min(h)
        c = {m: h.get(-m, R.poly([])) for m in range(1, M + 1)}
        U = {}
        if M >= 2:
            U[M - 1] = -c[M] * (K(1) / K(M - 1))
            for j in range(M - 1, 1, -1):
                U[j - 1] = (U[j] * lam - c[j]) * (K(1) / K(j - 1))
            C = c[1] - U[1] * lam
        else:
            C = c[1]
        pref = K(Fr(2 * L + 1) / (2 * a) ** (L + 1))
        I = L - 1                                    # r^2 * r^-(L+1)
        A_, B_ = panel
        mid = A_ + 1 if B_ is None else (A_ + B_) / 2
        sa, sb = (1 if mid > a else -1), (1 if mid > b else -1)
        branches = [(Fr(1), 1 + kap, a + kap * b), (Fr(-1), sa + kap * sb, -sa * a - kap * sb * b)]
        out = []
        for sgn, al, be in branches:
            ul = R.lin(be, al)
            acc = S[Np]
            for k in range(Np - 1, -1, -1):
                acc = acc * ul + S[k]
            e0, lr = -z * be, -z * al
            out.append((e0, lr, acc * (pref * K(sgn)), I, None, 0, ()))
            if M >= 2:
                Kp = M - 1
                accU = R.poly([])
                for k in range(1, Kp + 1):
                    accU = accU * ul + U[k]
                pole = -be / al
                out.append((e0, lr, accU * (pref * K(sgn) / K(al) ** Kp), I, pole, Kp, ()))
            out.append((Fr(0), Fr(0), C * (-pref * K(sgn)), I, None, 0, ((z * al, z * be),)))
        return out

    def f_terms_equal(self, L, panel, it, kt):
        K, R, p = self.K, self.R, self.p
        mu = p['zZ'] + p['zW']
        A_, B_ = panel
        inner = B_ is not None and B_ <= p['a']
        dens = AB.density(L, inner, p['a'], mu, it, kt, K)
        m = K(-(2 * L + 1))
        out = []
        for (e0, bb), d in dens.items():
            off, P = _laurent_from(R, {n + 2: m * v for n, v in d.items()})
            out.append((Fr(e0), Fr(bb), P, -off, None, 0, ()))
        return out

    # ------------------------------------------------------------------ assembly
    def run(self, Lmax):
        K, R = self.K, self.R
        it, kt = NE.bessel_tables(K, Lmax + 2)
        Qgen = self._Q_steps(Lmax) if self.matched else None
        for L in range(Lmax + 1):
            Q = next(Qgen) if Qgen else None
            total = {}
            for panel in self.panels:
                fts = self.f_terms_matched(L, Q, panel) if self.matched else self.f_terms_equal(L, panel, it, kt)
                pts, c_lo, c_hi = self.psi(L, panel, it, kt)
                red = X.Reducer(R, panel[0], panel[1])
                for (e0, lam, N, I, pole, Kp, E) in fts:
                    if N.degree() < 0:
                        continue
                    for (ck, b, (off, P), E2) in pts:
                        if P.degree() < 0:
                            continue
                        red.put_rational(lam + b, E + E2, X.ck_shift(ck, e0), N * P, I - off, pole, Kp)
                res = red.reduce()
                for k, v in res.items():
                    e0, atoms = k
                    if TAG_LO in atoms or TAG_HI in atoms:
                        cdict = c_lo if TAG_LO in atoms else c_hi
                        rest = (e0, tuple(a for a in atoms if a not in (TAG_LO, TAG_HI)))
                        for ck2, w in cdict.items():
                            X.add(total, X.ck_mul(rest, ck2), v * w)
                    else:
                        X.add(total, k, v)
            norm = K(Fr(16, (2 * L + 1) ** 2))
            yield L, X.canon({k: v * norm for k, v in total.items()}, K)


_FACT = [1]


def _fact(k):
    while len(_FACT) <= k:
        _FACT.append(_FACT[-1] * len(_FACT))
    return _FACT[k]


def coefficients(case, Lmax, K):
    return [q for _, q in Engine(case, K).run(Lmax)]


def generic_sequence(Qs, p, seed=2026):
    """Random generic combination of the coefficient sequences (as aabc_na_posctrl): list of nmod."""
    import random
    keys = sorted(set(k for q in Qs for k in q), key=repr)
    rnd = random.Random(seed)
    w = {k: flint.nmod(rnd.randrange(1, p), p) for k in keys}
    zero = flint.nmod(0, p)
    return [sum((w[k] * q.get(k, zero) for k in keys), zero) for q in Qs], keys
