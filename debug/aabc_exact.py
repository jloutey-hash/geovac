"""Record AB (docs/predictions/2026-11-14_aa_bc_three_centre.md), Phase-2 engine: exact constant split of the
Legendre coefficients of (1s_A 1s_A | 1s_B 1s_C) in x = cos(angle C B A).

    (AA|BC) = NA - Y_mu - zeta_A W_mu,   W_mu = -dY_mu/dmu,   mu = 2 zeta_A,

and each piece's coefficient F_L / (2 pi N) = sum_c q_{L,c} c over constants c = e^{e0} {1, E1(z), ln v, gamma}, as in
na_exact (whose Field, Bessel tables, radial integrals J_n and constant evaluation are reused unchanged).

Kernels replacing NA's w_L = r<^L / r>^(L+1) (r<, r> of (r, D2)):
    Y : mu (2L+1) i_L(mu r<) k_L(mu r>)
    dY = dY/dmu : (2L+1) [ (2L+1) i_L(x) k_L(y) + x i_{L+1}(x) k_L(y) - y i_L(x) k_{L+1}(y) ],  x = mu r<, y = mu r>
(i_L' = i_{L+1} + L i_L / x, k_L' = -k_{L+1} + L k_L / x).  f's kernel is w_L - Y + zeta_A dY.
The density factor (Legendre multipole of e^{-zB rB - zC rC} about B) is na_exact's, copied.
Products of density and kernel series are polynomial convolutions in 1/r (flint), then integrated term by term with
na_exact.Integrals (finite part at r = 0).

Usage (functions; drivers below):
    pieces(case, Lmax, K, which=('NA', 'Y', 'dY'))   -> {piece: [ {key: q} per L ]}
"""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import flint

sys.path.insert(0, str(Path(__file__).parent))
import na_exact as E

CASES = {'A1': dict(D1='9/5', D2='13/5', zB='13/10', zC='1', zA='9/10'),
         'A2': dict(D1='3/2', D2='12/5', zB='1', zC='6/5', zA='4/5')}
REF_X = {'A1': '1/2', 'A2': '-1/3'}


def params(case):
    g = {k: Fraction(v) for k, v in CASES[case].items()}
    return g['D1'], g['D2'], g['zB'], g['zC'], g['zA']


def _put(d, e0, b, n, v):
    dd = d.setdefault((Fraction(e0), Fraction(b)), {})
    dd[n] = dd[n] + v if n in dd else v


def _ksum(K, tab, arg):
    return sum((K.pw(arg, -j) * v for j, v in tab.items()), K(0))


def density(L, inner, D1, zH, it, kt, K):
    """na_exact's density factor (G_L / -2 without e^{-zX r}) as {(e0, b): {n: coef}}: e^{e0} r^n e^{b r}."""
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
    """Kernel in the region r < D2 (near) or r > D2 as {(e0, b): {n: coef}}."""
    d = {}
    c = K(2 * L + 1)
    if which == 'NA':
        if near:
            _put(d, 0, 0, L, K.pw(D2, -L - 1))
        else:
            _put(d, 0, 0, -L - 1, K.pw(D2, L))
        return d
    if near:                                          # x = mu r, y = mu D2
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
    else:                                             # x = mu D2, y = mu r
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
    """product of two {n: coef} as {n: coef} via flint polynomials in 1/r."""
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


def pieces(case, Lmax, K, which=('NA', 'Y', 'dY'), params_override=None):
    D1, D2, zX, zH, zA = params_override or params(case)
    mu = 2 * zA
    assert D1 < D2
    it, kt = E.bessel_tables(K, Lmax + 2)
    I = E.Integrals(K)
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
                                E.add(res, (e0 + e1, kind, arg), vv * u)
            out[w].append({k: v for k, v in res.items() if v != 0})
    return out


def combine(pc, zA, K):
    """f = NA - Y + zA dY  (W = -dY), per L."""
    out = []
    za = K(Fraction(zA))
    for qN, qY, qd in zip(pc['NA'], pc['Y'], pc['dY']):
        r = dict(qN)
        for k, v in qY.items():
            E.add(r, k, -v)
        for k, v in qd.items():
            E.add(r, k, za * v)
        out.append({k: v for k, v in r.items() if v != 0})
    return out


if __name__ == '__main__':
    # self-check: the engine with NA's kernel reproduces na_exact exactly (both fields)
    import time
    for case in CASES:
        D1, D2, zB, zC, zA = params(case)
        for K in (E.Field(None), E.Field((1 << 61) - 1)):
            t = time.time()
            mine = pieces(case, 30, K, which=('NA',))['NA']
            ref = E.coefficients(str(D1), str(D2), str(zB), str(zC), 30, K)
            ok = all({k: v for k, v in a.items()} == {k: v for k, v in b.items()} for a, b in zip(mine, ref))
            print(case, 'p' if K.p else 'Q', 'NA-kernel == na_exact:', ok, round(time.time() - t, 1), 's')
