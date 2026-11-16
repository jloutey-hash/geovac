"""Record AP (docs/predictions/2026-11-14_aa_bc_p_orbitals.md), Phase-2 engine: exact constant split of the Legendre
coefficients (in x = cos(angle C B A)) of

    PC = (1s_A 1s_A | 1s_B 2pz_C)   and   PA = (2ps_A 2ps_A | 1s_B 1s_C),

extending debug/aabc_exact.py (record AB).  na_exact's Field, Bessel tables, Integrals and value are unchanged.

PC: AB's kernel f = NA - Y + zA dY (V_A of a 1s density), density multipoles from AB's by z - D1 = r cos - D1:
    dens'_L = [ r (L dens_{L-1} + (L+1) dens_{L+1}) - (2L+1) D1 dens_L ] / (2L+1).
PA: PA = d^2/dD2^2 int rho U + int rho W0, the kernels written in AB's family
    U  <->  3/(2 zA^2) NA - 3/(2 zA^2) Y + 33/(16 zA) dY - 9/8 d2Y + zA/4 d3Y
    W0 <->  NA - Y + 5 zA/4 dY - zA^2/2 d2Y
(Y <-> e^{-mu s}/s, d^kY = d^k/dmu^k Y; dY <-> -e^{-mu s}, d2Y <-> s e^{-mu s}, d3Y <-> -s^2 e^{-mu s}).

Kernels are built as symbolic monomials  c mu^pa D2^pd r^n e^{t0 mu D2} e^{s mu r}  (`kernel_mono`), so mu- and
D2-derivatives are exact (`d_mu`, `d_D2`); `evaluate` turns them into aabc_exact's {(e0, b): {n: coef}} form.
`boundary_D2` is the jump term of d^2/dD2^2 at r = D2:  [-2 r^2 e^{-zB r} dens(r) (dK_near - dK_far)]_{r = D2}.
"""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import aabc_exact as X  # noqa: E402
import na_exact as E  # noqa: E402

CASES, REF_X, params = X.CASES, X.REF_X, X.params
FAST = True   # factored kernel path; False = the AP0-checked symbolic path (identical output, checked)


# ---------------------------------------------------------------- symbolic kernel monomials


def _addm(d, key, v):
    d[key] = d[key] + v if key in d else v


def kernel_mono(L, near, which, it, kt, K):
    """{(pa, pd, n, t0, s): c} for the NA or Y kernel in the region r < D2 (near) or r > D2."""
    d = {}
    if which == 'NA':
        if near:
            _addm(d, (0, -L - 1, L, 0, 0), K(1))
        else:
            _addm(d, (0, L, -L - 1, 0, 0), K(1))
        return d
    assert which == 'Y'
    c = K(2 * L + 1)
    if near:     # mu (2L+1) i_L(mu r) k_L(mu D2):  k_L(y) = e^{-y} sum kt_i y^-i ; i_L(x) = sum_s e^{s x} it_j x^-j
        for i, vk in kt[L].items():
            for s in (1, -1):
                for j, vi in it[L][s].items():
                    _addm(d, (1 - i - j, -i, -j, -1, s), c * vk * vi)
    else:        # mu (2L+1) i_L(mu D2) k_L(mu r)
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
    """-> {(e0, b): {n: coef}} at rational mu, D2."""
    acc, pmu, pD2 = {}, {}, {}
    for (pa, pd, n, t0, s), c in m.items():
        if pa not in pmu:
            pmu[pa] = K.pw(mu, pa)
        if pd not in pD2:
            pD2[pd] = K.pw(D2, pd)
        _addm(acc, (t0, s, n), c * pmu[pa] * pD2[pd])
    d = {}
    for (t0, s, n), v in acc.items():
        X._put(d, t0 * mu * D2, s * mu, n, v)
    return d


def kernel_family(L, near, name, it, kt, K):
    """Monomials of 'NA', 'Y', 'dY', 'd2Y', 'd3Y' (mu-derivatives of Y)."""
    if name == 'NA':
        return kernel_mono(L, near, 'NA', it, kt, K)
    m = kernel_mono(L, near, 'Y', it, kt, K)
    for _ in range({'Y': 0, 'dY': 1, 'd2Y': 2, 'd3Y': 3}[name]):
        m = d_mu(m, K)
    return m


def combo(L, near, coeffs, it, kt, K, nD2=0):
    """sum_k coeffs[k] * kernel_k, then nD2 derivatives in D2."""
    tot = {}
    for name, a in coeffs.items():
        for key, v in kernel_family(L, near, name, it, kt, K).items():
            _addm(tot, key, a * v)
    for _ in range(nD2):
        tot = d_D2(tot, K)
    return {k: v for k, v in tot.items() if v != 0}


# ---------------------------------------------------------------- factored fast path (same output as combo+evaluate)
# The Y kernel is mu (2L+1) g(mu D2) h(mu r): one Bessel factor in D2, one in r.  Keeping the two factors apart
# and applying mu-derivatives by Leibniz makes each kernel O(L) terms instead of O(L^2).


def _factors(L, near, it, kt, K):
    """G = {(pa, pd, t0): c} (mu^pa D2^pd e^{t0 mu D2}, includes mu (2L+1)); H = {(pa, n, s): c}."""
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


def _dG_mu(G, K):
    out = {}
    for (pa, pd, t0), c in G.items():
        if pa:
            _addm(out, (pa - 1, pd, t0), K(pa) * c)
        if t0:
            _addm(out, (pa, pd + 1, t0), K(t0) * c)
    return out


def _dG_D2(G, K):
    out = {}
    for (pa, pd, t0), c in G.items():
        if pd:
            _addm(out, (pa, pd - 1, t0), K(pd) * c)
        if t0:
            _addm(out, (pa + 1, pd, t0), K(t0) * c)
    return out


def _dH_mu(H, K):
    out = {}
    for (pa, n, s), c in H.items():
        if pa:
            _addm(out, (pa - 1, n, s), K(pa) * c)
        if s:
            _addm(out, (pa, n + 1, s), K(s) * c)
    return out


def kernel_eval(L, near, coeffs, it, kt, K, mu, D2, nD2=0):
    """== evaluate(combo(L, near, coeffs, it, kt, K, nD2), mu, D2, K), computed through the factored form."""
    from math import comb
    d = {}
    if 'NA' in coeffs:
        na = {k: coeffs['NA'] * v for k, v in kernel_mono(L, near, 'NA', it, kt, K).items()}
        for _ in range(nD2):
            na = d_D2(na, K)
        for k, v in evaluate(na, mu, D2, K).items():
            for n, c in v.items():
                X._put(d, k[0], k[1], n, c)
    order = {'Y': 0, 'dY': 1, 'd2Y': 2, 'd3Y': 3}
    amax = max([order[k] for k in coeffs if k != 'NA'], default=-1)
    if amax < 0:
        return d
    G0, H0 = _factors(L, near, it, kt, K)
    Gs, Hs = [G0], [H0]
    for _ in range(amax):
        Gs.append(_dG_mu(Gs[-1], K))
        Hs.append(_dH_mu(Hs[-1], K))
    for _ in range(nD2):
        Gs = [_dG_D2(g, K) for g in Gs]
    pmu, pD2 = {}, {}

    def pw(cache, q, k):
        if k not in cache:
            cache[k] = K.pw(q, k)
        return cache[k]
    Gv = []          # Gv[k] = {t0: value}
    for g in Gs:
        v = {}
        for (pa, pd, t0), c in g.items():
            _addm(v, t0, c * pw(pmu, mu, pa) * pw(pD2, D2, pd))
        Gv.append(v)
    Hv = []          # Hv[k] = {(s, n): value}
    for h in Hs:
        v = {}
        for (pa, n, s), c in h.items():
            _addm(v, (s, n), c * pw(pmu, mu, pa))
        Hv.append(v)
    acc = {}
    for name, a in coeffs.items():
        if name == 'NA':
            continue
        o = order[name]
        for k in range(o + 1):
            w = a * K(comb(o, k))
            for t0, gv in Gv[k].items():
                wg = w * gv
                for (s, n), hv in Hv[o - k].items():
                    _addm(acc, (t0, s, n), wg * hv)
    for (t0, s, n), v in acc.items():
        X._put(d, t0 * mu * D2, s * mu, n, v)
    return d


# ---------------------------------------------------------------- densities


def dens_1s(L, inner, D1, zC, it, kt, K):
    return X.density(L, inner, D1, zC, it, kt, K)


def dens_pz(L, inner, D1, zC, it, kt, K):
    """Multipole L of (z - D1) times AB's 1s_C density factor (powers of r shift by one for the r cos term)."""
    d = {}
    inv = K(1) / K(2 * L + 1)
    parts = []
    if L >= 1:
        parts.append((K(L) * inv, 1, X.density(L - 1, inner, D1, zC, it, kt, K)))
    parts.append((K(L + 1) * inv, 1, X.density(L + 1, inner, D1, zC, it, kt, K)))
    parts.append((-K(Fraction(D1)), 0, X.density(L, inner, D1, zC, it, kt, K)))
    for w, sh, dd in parts:
        for (e0, b), poly in dd.items():
            for n, v in poly.items():
                X._put(d, e0, b, n + sh, w * v)
    return d


def dens_px(L, inner, D1, zC, it, kt, K):
    """Record AQ: m = 1 multipole of x e^{-zC r_C} in AB's normalization, r (dens_{L-1} - dens_{L+1}) / (2L+1),
    from sqrt(1 - t^2) P_L^1 = L(L+1)/(2L+1) (P_{L-1} - P_{L+1}); the series is in P_L^1(x).  L = 0 is absent."""
    d = {}
    if L == 0:
        return d
    inv = K(1) / K(2 * L + 1)
    for w, dd in ((inv, X.density(L - 1, inner, D1, zC, it, kt, K)),
                  (-inv, X.density(L + 1, inner, D1, zC, it, kt, K))):
        for (e0, b), poly in dd.items():
            for n, v in poly.items():
                X._put(d, e0, b, n + 1, w * v)
    return d


def _dacc(out, dd, w, shift, K):
    for (e0, b), poly in dd.items():
        for n, v in poly.items():
            X._put(out, e0, b, n + shift, w * v)


def _tmul(E, L, K):
    """(T E)_L = [L E_{L-1} + (L+1) E_{L+1}] / (2L+1): multiply by t = cos(theta) in AB's normalization."""
    out = {}
    inv = K(1) / K(2 * L + 1)
    if L >= 1:
        _dacc(out, E(L - 1), K(L) * inv, 0, K)
    _dacc(out, E(L + 1), K(L + 1) * inv, 0, K)
    return out


def dens_dz2(L, inner, D1, zC, it, kt, K):
    """Record AD: multipole L of h2(r - C) e^{-zC r_C}, h2 = (3 r^2 t^2 - r^2 - 4 D1 r t + 2 D1^2)/2 about B."""
    cache = {}

    def E(k):
        if k not in cache:
            cache[k] = X.density(k, inner, D1, zC, it, kt, K)
        return cache[k]

    def TE(k):
        return _tmul(E, k, K)
    out = {}
    D1q = K(Fraction(D1))
    half = K(Fraction(1, 2))
    _dacc(out, _tmul(TE, L, K), K(3) * half, 2, K)       # 3 r^2 T^2
    _dacc(out, E(L), -half, 2, K)                         # - r^2
    _dacc(out, TE(L), K(-2) * D1q, 1, K)                  # - 4 D1 r T  (x 1/2)
    _dacc(out, E(L), D1q * D1q, 0, K)                     # + 2 D1^2   (x 1/2)
    return out


def _ecache(inner, D1, zC, it, kt, K):
    cache = {}

    def E(k):
        if k not in cache:
            cache[k] = X.density(k, inner, D1, zC, it, kt, K)
        return cache[k]
    return E


def dens_dxz(L, inner, D1, zC, it, kt, K):
    """Record AX: m = 1 multipole of x (z - D1) e^{-zC r_C} about B in AB's normalization (the P_L^1 expansion
    coefficient divided by 2L + 1, as dens_px).  (z - D1) e^{-zC s} has AB-normalized Legendre coefficients
    g_k = r (T e)_k - D1 e_k; then dens_px's step r (g_{L-1} - g_{L+1}) / (2L+1).  L = 0 is absent."""
    out = {}
    if L == 0:
        return out
    E = _ecache(inner, D1, zC, it, kt, K)
    D1q = K(Fraction(D1))

    def g(k):
        d = {}
        _dacc(d, _tmul(E, k, K), K(1), 1, K)
        _dacc(d, E(k), -D1q, 0, K)
        return d
    inv = K(1) / K(2 * L + 1)
    _dacc(out, g(L - 1), inv, 1, K)
    _dacc(out, g(L + 1), -inv, 1, K)
    return out


def dens_dx2y2(L, inner, D1, zC, it, kt, K):
    """Record AX: m = 2 multipole of (x^2 - y^2) e^{-zC r_C} = r^2 (1 - t^2) cos(2 phi) e^{-zC s} about B, in AB's
    normalization (P_L^2 expansion coefficient / (2L+1); P_L^2 = (1 - t^2) P_L'', no Condon-Shortley phase).
    sqrt(1-t^2) P_l^m = (P_{l+1}^{m+1} - P_{l-1}^{m+1}) / (2l+1) twice, on c_k = (2k+1) e_k:
        a_L = e_{L-2}/(2L-1) - 2(2L+1) e_L/((2L-1)(2L+3)) + e_{L+2}/(2L+3),   dens = r^2 a_L / (2L+1).
    L = 0, 1 are absent."""
    out = {}
    if L < 2:
        return out
    E = _ecache(inner, D1, zC, it, kt, K)
    _dacc(out, E(L - 2), K(1) / K((2 * L - 1) * (2 * L + 1)), 2, K)
    _dacc(out, E(L), K(-2) / K((2 * L - 1) * (2 * L + 3)), 2, K)
    _dacc(out, E(L + 2), K(1) / K((2 * L + 1) * (2 * L + 3)), 2, K)
    return out


# ---------------------------------------------------------------- radial integration


def integrate(dens, kers, D1, D2, zB, I, K):
    """sum over AB's three regions of -2 int r^2 e^{-zB r} dens(r) K(r) dr -> {key: coef}.
    dens: {inner: {(e0,b):{n:c}}}, kers: {near: {(e0,b):{n:c}}}."""
    regions = [(Fraction(0), D1, True, True), (D1, D2, False, True), (D2, None, False, False)]
    if D1 == D2:   # record AE: equal arms, the middle region is empty
        regions = [regions[0], regions[2]]
    res = {}
    m2 = K(-2)
    for (A, B, inner, near) in regions:
        for (e0d, bd), pd in dens[inner].items():
            for (e0k, bk), pk in kers[near].items():
                P = X._conv(pd, pk, K)
                if not P:
                    continue
                bb = bd + bk - zB
                e0 = e0d + e0k
                tab = I.table(bb, A, B, min(min(P) + 2, -1), max(max(P) + 2, -1))
                for n, v in P.items():
                    vv = m2 * v
                    for (e1, kind, arg), u in tab[n + 2].items():
                        E.add(res, (e0 + e1, kind, arg), vv * u)
    return res


def boundary_D2(dens_outer, dK_near, dK_far, D2, zB, K):
    """[-2 r^2 e^{-zB r} dens(r) (dK_near - dK_far)]_{r = D2} as {key: coef} (exponentials only)."""
    jump = {}
    for src, sg in ((dK_near, 1), (dK_far, -1)):
        for (e0, b), poly in src.items():
            for n, v in poly.items():
                X._put(jump, e0, b, n, K(sg) * v)
    res = {}
    for (e0d, bd), pd in dens_outer.items():
        for (e0k, bk), pk in jump.items():
            P = X._conv(pd, pk, K)
            bb = bd + bk - zB
            for n, v in P.items():
                E.add(res, (e0d + e0k + bb * D2, '1', None), K(-2) * v * K.pw(D2, n + 2))
    return res


# ---------------------------------------------------------------- the two integrals


def coefficients(case, Lmax, K, integral, params_override=None, pieces_out=False):
    """List over L of {key: q} with F_L = 2 pi N sum q c  (N = N_B N_C of the BC pair; PA's A pair is in V)."""
    D1, D2, zB, zC, zA = params_override or params(case)
    mu = 2 * zA
    assert D1 <= D2
    it, kt = E.bessel_tables(K, Lmax + 3)
    I = E.Integrals(K)
    za = K(Fraction(zA))
    out, outp = [], []

    def kev(L, near, co, nD2):
        if FAST:
            return kernel_eval(L, near, co, it, kt, K, mu, D2, nD2)
        return evaluate(combo(L, near, co, it, kt, K, nD2=nD2), mu, D2, K)
    for L in range(Lmax + 1):
        dfun = {'PC': dens_pz, 'PX': dens_px, 'PD': dens_dz2, 'DX': dens_dxz,
                'DD': dens_dx2y2}.get(integral, dens_1s)
        dens = {inner: dfun(L, inner, D1, zC, it, kt, K) for inner in (True, False)}
        if integral in ('PC', 'AB', 'PX', 'PD', 'DX', 'DD'):
            co = {'NA': K(1), 'Y': K(-1), 'dY': za}
            kers = {near: kev(L, near, co, 0) for near in (True, False)}
            res = integrate(dens, kers, D1, D2, zB, I, K)
            parts = None
        else:
            cU = {'NA': K(3) / (K(2) * za * za), 'Y': -K(3) / (K(2) * za * za), 'dY': K(33) / (K(16) * za),
                  'd2Y': K(Fraction(-9, 8)), 'd3Y': za / K(4)}
            cW = {'NA': K(1), 'Y': K(-1), 'dY': K(5) * za / K(4), 'd2Y': -za * za / K(2)}
            kU2 = {near: kev(L, near, cU, 2) for near in (True, False)}
            kU1 = {near: kev(L, near, cU, 1) for near in (True, False)}
            kW = {near: kev(L, near, cW, 0) for near in (True, False)}
            pU = integrate(dens, kU2, D1, D2, zB, I, K)
            pB = boundary_D2(dens[False], kU1[True], kU1[False], D2, zB, K)
            pW = integrate(dens, kW, D1, D2, zB, I, K)
            res = {}
            for p in (pU, pB, pW):
                for k, v in p.items():
                    E.add(res, k, v)
            parts = {'U2': pU, 'B': pB, 'W0': pW}
        out.append({k: v for k, v in res.items() if v != 0})
        if pieces_out:
            outp.append(parts)
    return (out, outp) if pieces_out else out


def twopiN(case, integral):
    """2 pi N as an mpmath number (current precision)."""
    import mpmath as mp
    D1, D2, zB, zC, zA = params(case)
    q = lambda v: mp.mpf(v.numerator) / v.denominator  # noqa: E731
    if integral == 'PD':   # 2 pi N_B N_d, N_d = sqrt(5 (2 zC)^7 / (4 pi 720))
        return 2 * mp.pi * mp.sqrt(q(zB) ** 3 / mp.pi) * mp.sqrt(5 * (2 * q(zC)) ** 7 / (4 * mp.pi * 720))
    if integral in ('DX', 'DD'):   # record AX: 2 pi N_B N_Q, N_xz = sqrt(15/(4 pi)) R_n, N_x2y2 = sqrt(15/(16 pi)) R_n
        Rn = mp.sqrt((2 * q(zC)) ** 7 / 720)
        ang = mp.sqrt(15 / (4 * mp.pi)) if integral == 'DX' else mp.sqrt(15 / (16 * mp.pi))
        return 2 * mp.pi * mp.sqrt(q(zB) ** 3 / mp.pi) * ang * Rn
    return 2 * mp.sqrt(q(zB) ** 3 * q(zC) ** (5 if integral in ('PC', 'PX') else 3))


# ---------------------------------------------------------------- AP0 self-checks


def selfcheck():
    import time
    import mpmath as mp
    ok = True
    for case in CASES:
        D1, D2, zB, zC, zA = params(case)
        mu = 2 * zA
        for K in (E.Field(None), E.Field((1 << 61) - 1)):
            t = time.time()
            # (i) PC machinery with 1s density = AB
            mine = coefficients(case, 30, K, 'AB')
            pc = X.pieces(case, 30, K)
            ref = X.combine(pc, zA, K)
            a = all(m == r for m, r in zip(mine, ref))
            # (ii) mu-jet of Y == aabc_exact's dY kernel
            it, kt = E.bessel_tables(K, 33)
            b = True
            for L in range(31):
                for near in (True, False):
                    mj = evaluate(kernel_family(L, near, 'dY', it, kt, K), mu, D2, K)
                    rf = X.kernel(L, near, D2, mu, 'dY', it, kt, K)
                    clean = lambda d: {k: {n: v for n, v in p.items() if v != 0} for k, p in d.items()}  # noqa
                    m1 = {k: v for k, v in clean(mj).items() if v}
                    r1 = {k: v for k, v in clean(rf).items() if v}
                    b &= m1 == r1
            print(case, 'p' if K.p else 'Q', '(i) AB reproduced:', a, ' (ii) mu-jet == dY:', b,
                  round(time.time() - t, 1), 's', flush=True)
            ok &= a and b
    # (iii) D2-jet with jump term vs 2000-digit central difference of the exact F_L of int rho U
    K = E.Field(None)
    mp.mp.dps = 2000
    h = mp.mpf(10) ** -300
    for case in CASES:
        D1, D2, zB, zC, zA = params(case)
        Ls = (0, 1, 5, 20)
        _, pp = coefficients(case, max(Ls), K, 'PA', pieces_out=True)

        def Uval(d2, L):
            """exact int rho U at D2 = d2 (rational), by the same engine with nD2 = 0."""
            it, kt = E.bessel_tables(K, L + 3)
            I = E.Integrals(K)
            za = K(Fraction(zA))
            cU = {'NA': K(3) / (K(2) * za * za), 'Y': -K(3) / (K(2) * za * za), 'dY': K(33) / (K(16) * za),
                  'd2Y': K(Fraction(-9, 8)), 'd3Y': za / K(4)}
            dens = {inner: dens_1s(L, inner, D1, zC, it, kt, K) for inner in (True, False)}
            kU = {near: evaluate(combo(L, near, cU, it, kt, K), 2 * zA, d2, K) for near in (True, False)}
            q = integrate(dens, kU, D1, d2, zB, I, K)
            return sum(mp.mpf(int(v.p)) / int(v.q) * E.value(k, mp) for k, v in q.items())

        hq = Fraction(1, 10 ** 300)
        for L in Ls:
            fd = (Uval(D2 + hq, L) - 2 * Uval(D2, L) + Uval(D2 - hq, L)) / h ** 2
            ex = sum(mp.mpf(int(v.p)) / int(v.q) * E.value(k, mp)
                     for part in ('U2', 'B') for k, v in pp[L][part].items())
            err = abs(fd - ex)
            print(case, 'L', L, '(iii) |D2-jet - central diff| =', mp.nstr(err, 3), ' value', mp.nstr(ex, 12),
                  flush=True)
            ok &= err < mp.mpf(10) ** -500
    print('AP0 self-checks pass:', ok)
    return ok


if __name__ == '__main__':
    selfcheck()
