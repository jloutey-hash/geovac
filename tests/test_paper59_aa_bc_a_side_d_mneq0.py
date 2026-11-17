"""Guard for record AY (docs/predictions/2026-11-16_aa_bc_a_side_d_mneq0.md, scored GO): the three-centre integrals

    DNX = (dnx dnx | 1s_B 1s_C),  w M w = sqrt3 (n.w)(e'.w)
    DNY = (dny dny | 1s_B 1s_C),  w M w = sqrt3 (n.w)(e_y.w)
    DXY = (dxy dxy | 1s_B 1s_C),  w M w = sqrt3 (e'.w)(e_y.w)
    DXX = (dxx dxx | 1s_B 1s_C),  w M w = (sqrt3/2)((e'.w)^2 - (e_y.w)^2)

(the four m != 0 real Slater 3d pairs on A; n = A/|A|, e' in the plane of the centres perpendicular to B-A, e_y out
of it) over Slater functions, as series in x = cos(angle C B A).  Paper 59 sec:scope, the paragraph "[MEASURED] The
rest of the 3d shell on A follows the same rule".  Written 2026-11-16 as a separate step after the computation and
the paper edit (CLAUDE.md sec. 9).  The record's engine (debug/aay_exact.py: universal d kernels U4, U2, U0, the
A-frame operators and the CS-free ladder) is ported here; its base layer (aaa_exact's Engine: the AB density, the
kernel jets and the three-region integrator) is the port already in record AS's guard,
tests/test_paper59_aa_bc_a_side.py, imported with importlib (which in turn imports tests/test_paper59_aa_bc_p.py and
tests/test_paper59_aa_bc.py).

I(x) = sum_L [f0_L P_L(x) + f2_L P_L^2(x) + f4_L P_L^4(x)], P_L^m = (1 - x^2)^{m/2} d^m P_L / dx^m (no Condon-Shortley
phase), m the azimuthal order about n (A's frame).  DNX, DNY carry m in {0, 2}; DXY, DXX m in {0, 4}.  The total as a
plain Legendre series is c = f0 + leg[(1 - x^2) D^2 f2] + leg[(1 - x^2)^2 D^4 f4], exact at truncation.

Claims guarded, each with the wrong answer it rejects (details in each test's docstring):

(0) the exact symmetry identities f0(DNY) = f0(DNX), f2(DNY) = -f2(DNX), f0(DXY) = f0(DXX), f4(DXY) = -f4(DXX), and
    the sum rule f0 of DA + DNX + DNY + DXY + DXX = 5 S (S the spherically averaged pair), key for key over Q.
(1) per-channel values f{m}_L against an independent quadrature (Gauss's-law multipoles of the A density, the
    A-frame harmonic content of (w M w)^2 found by quadrature on the sphere, P_L^m projection; no U kernels, no
    D2-derivatives, no ladder).
(2) the summed series at A1, x = 1/2, against the R2 Gaussian-fit referee (n = 34) and the record's 3600-digit
    values (L <= 400); the re-expansion of the total c equals the split sum.
(3) the constant set: mod 2^61 - 1, every f0 (from L = 2), f2 (from L = 2), f4 (from L = 4) up to L = 12 has exactly
    AB's 21 constants at the same case; no Euler gamma.
(4) decay: beta_eff between the midpoints of in-test toy exponents.  DNY, DXY totals on (x* - x)^{9/2}; DNX, DXX
    totals on 7/2; f2, f4 on the 7/2 branch in their P_L^m bases.  Rejects "no cancellation".
(5) P4: R_2 = (f0/ta)/(f2/tb2) of DNY positive and tending to x*^2 - 1; R_4 of DXY negative and tending to
    -(x*^2 - 1)^2, by a 1/L extrapolation; bands set from the in-test values.

Fire-test plants and outcomes: the comment block at the end of this file.
"""
from __future__ import annotations

import importlib.util
from fractions import Fraction
from math import factorial
from pathlib import Path

import mpmath as mp
import pytest

_spec = importlib.util.spec_from_file_location('_as_guard_ay', Path(__file__).with_name('test_paper59_aa_bc_a_side.py'))
AS = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(AS)
AP, AB = AS.AP, AS.AB

Field, params, add = AB.Field, AB.params, AB.add
PRIME = (1 << 61) - 1
NEW = ('DNX', 'DNY', 'DXY', 'DXX')
MSET = {'DNX': (0, 2), 'DNY': (0, 2), 'DXY': (0, 4), 'DXX': (0, 4), 'DA': (0,)}


# ---------------------------------------------------------------------------------------------------------------
# Record-AY engine (port of debug/aay_exact.py)
# ---------------------------------------------------------------------------------------------------------------


def fam_u(zA):
    """Universal d-pair kernels (N_d^2 inside): V_M = (d M d)^2 U4 + (d M^2 d) U2 + tr(M^2) U0, tr M^2 = 3/2.
    Keys: 'NA' = 1/s, integer k = d^k Y / dmu^k."""
    z = Fraction(zA)
    U4 = {'NA': Fraction(5, 3) / z ** 4, 0: -Fraction(5, 3) / z ** 4, 1: Fraction(965, 384) / z ** 3,
          2: -Fraction(325, 192) / z ** 2, 3: Fraction(185, 288) / z, 4: Fraction(-5, 36), 5: z / 72}
    U2 = {'NA': Fraction(8, 3) / z ** 2, 0: -Fraction(8, 3) / z ** 2, 1: Fraction(31, 8) / z, 2: Fraction(-29, 12),
          3: Fraction(7, 9) * z, 4: -z ** 2 / 9}
    U0 = {'NA': Fraction(2, 3), 0: Fraction(-2, 3), 1: Fraction(11, 12) * z, 2: -z ** 2 / 2, 3: z ** 3 / 9}
    return {'U4': U4, 'U2': U2, 'U0': U0}


# a polynomial in (d_+, d_-, d_z) is {(a, b, c): Fraction};  d_x = (d_+ + d_-)/2,  d_y^2 = -(d_+ - d_-)^2/4


def _pmul(p, q):
    out = {}
    for (a, b, c), u in p.items():
        for (d, e, f), v in q.items():
            k = (a + d, b + e, c + f)
            out[k] = out.get(k, 0) + u * v
    return {k: v for k, v in out.items() if v}


def _padd(*ps, w=None):
    out = {}
    for i, p in enumerate(ps):
        a = 1 if w is None else w[i]
        for k, v in p.items():
            out[k] = out.get(k, 0) + a * v
    return {k: v for k, v in out.items() if v}


H = Fraction(1, 2)
DX = {(1, 0, 0): H, (0, 1, 0): H}
DX2 = _pmul(DX, DX)
DY2 = {(2, 0, 0): -Fraction(1, 4), (1, 1, 0): Fraction(1, 2), (0, 2, 0): -Fraction(1, 4)}
DZ2 = {(0, 0, 2): Fraction(1)}
LAP = _padd(DX2, DY2, DZ2)


def operators(name):
    """{'U4': (d M d)^2, 'U2': d M^2 d, 'U0': tr M^2} in the A frame (x' = e', y' = e_y, z' = n)."""
    q = Fraction(3, 4)
    if name == 'DA':
        dMd = _padd(DZ2, LAP, w=[Fraction(3, 2), -H])
        return {'U4': _pmul(dMd, dMd), 'U2': _padd(DX2, DY2, DZ2, w=[Fraction(1, 4), Fraction(1, 4), 1]),
                'U0': {(0, 0, 0): Fraction(3, 2)}}
    if name == 'DNX':
        return {'U4': _pmul(_pmul(DX2, DZ2), {(0, 0, 0): Fraction(3)}), 'U2': _padd(DX2, DZ2, w=[q, q]),
                'U0': {(0, 0, 0): Fraction(3, 2)}}
    if name == 'DNY':
        return {'U4': _pmul(_pmul(DY2, DZ2), {(0, 0, 0): Fraction(3)}), 'U2': _padd(DY2, DZ2, w=[q, q]),
                'U0': {(0, 0, 0): Fraction(3, 2)}}
    if name == 'DXY':
        return {'U4': _pmul(_pmul(DX2, DY2), {(0, 0, 0): Fraction(3)}), 'U2': _padd(DX2, DY2, w=[q, q]),
                'U0': {(0, 0, 0): Fraction(3, 2)}}
    if name == 'DXX':
        d = _padd(DX2, DY2, w=[1, -1])
        return {'U4': _pmul(_pmul(d, d), {(0, 0, 0): q}), 'U2': _padd(DX2, DY2, w=[q, q]),
                'U0': {(0, 0, 0): Fraction(3, 2)}}
    raise ValueError(name)


# ladder on single-source jets: a jet is {(j, p): Fraction} = coef R^{-p} G^{(j)}(R); a field is {(l, m): jet}


def _jd(jet):
    out = {}
    for (j, p), c in jet.items():
        if p:
            out[(j, p + 1)] = out.get((j, p + 1), 0) - p * c
        out[(j + 1, p)] = out.get((j + 1, p), 0) + c
    return out


def _jlin(*terms):
    out = {}
    for a, jet, sh in terms:
        if not a:
            continue
        for (j, p), c in jet.items():
            k = (j, p + sh)
            out[k] = out.get(k, 0) + a * c
    return {k: v for k, v in out.items() if v}


def _Aop(l, jet):
    return _jlin((1, _jd(jet), 0), (Fraction(-l), jet, 1))


def _Bop(l, jet):
    return _jlin((1, _jd(jet), 0), (Fraction(l + 1), jet, 1))


def _fadd(field, l, m, jet, a):
    if l < abs(m) or l < 0 or not a:
        return
    cur = field.setdefault((l, m), {})
    for k, v in jet.items():
        cur[k] = cur.get(k, 0) + a * v


def step(field, op):
    """One of d_+, d_-, d_z on f P_l^{|m|} e^{i m phi} (CS-free ladder)."""
    out = {}
    for (l, m), jet in field.items():
        am = abs(m)
        d = Fraction(1, 2 * l + 1)
        Aj, Bj = _Aop(l, jet), _Bop(l, jet)
        if op == 'z':
            _fadd(out, l + 1, m, Aj, d * (l - am + 1))
            _fadd(out, l - 1, m, Bj, d * (l + am))
            continue
        up = (op == '+' and m >= 0) or (op == '-' and m <= 0)
        nm = m + (1 if op == '+' else -1)
        if up:
            _fadd(out, l + 1, nm, Aj, d)
            _fadd(out, l - 1, nm, Bj, -d)
        else:
            _fadd(out, l + 1, nm, Aj, -d * (l - am + 1) * (l - am + 2))
            _fadd(out, l - 1, nm, Bj, d * (l + am) * (l + am - 1))
    return {k: {kk: vv for kk, vv in v.items() if vv} for k, v in out.items()}


def pole(L, m, poly):
    """[poly(d_+, d_-, d_z) (G_L(R) P_L^{|m|} e^{i m phi})] at the pole, as a jet."""
    out = {}
    for (a, b, c), w in poly.items():
        if m + a - b != 0:
            continue
        f = {(L, m): {(0, 0): Fraction(1)}}
        for _ in range(a):
            f = step(f, '+')
        for _ in range(b):
            f = step(f, '-')
        for _ in range(c):
            f = step(f, 'z')
        for (l, mm), jet in f.items():
            if mm == 0:
                for k, v in jet.items():
                    out[k] = out.get(k, 0) + w * v
    return {k: v for k, v in out.items() if v}


def channel_jets(name, L):
    """{m: {kernel: jet}}: f{m}_L = (L-m)!/(L+m)! sum_{sigma} [O_M G_L P_L^m e^{i sigma m phi}]_pole (m > 0)."""
    ops = operators(name)
    out = {}
    for m in (0, 1, 2, 3, 4):
        if m > L:
            continue
        rec = {}
        for kern, poly in ops.items():
            if m == 0:
                jet = pole(L, 0, poly)
            else:
                jp, jm = pole(L, m, poly), pole(L, -m, poly)
                jet = {k: jp.get(k, 0) + jm.get(k, 0) for k in set(jp) | set(jm)}
                c = Fraction(factorial(L - m), factorial(L + m))
                jet = {k: c * v for k, v in jet.items() if v}
            if jet:
                rec[kern] = jet
        if rec:
            out[m] = rec
    return out


NMAX = {'U4': 4, 'U2': 2, 'U0': 0}


def coefficients(case, Lmax, K, integral, eng):
    """{m: list over L of {key: q}}; value = 2 pi N_B N_C sum q c."""
    D2 = eng.D2
    fam = fam_u(eng.zA)
    out = {}
    for L in range(Lmax + 1):
        for m, rec in channel_jets(integral, L).items():
            acc = {}
            for kern, jet in rec.items():
                Gs = eng.G(kern + '_d', fam[kern], L, NMAX[kern])
                for (j, p), c in jet.items():
                    a = K(c) * K.pw(D2, -p)
                    for k, v in Gs[j].items():
                        add(acc, k, a * v)
            out.setdefault(m, [dict() for _ in range(Lmax + 1)])[L] = {k: v for k, v in acc.items() if v != 0}
    return out


def _leg_deriv(f):
    """Legendre coefficients of f': g_k = (2k+1) sum_{j > k, j - k odd} f_j."""
    n = len(f)
    g = [0 * f[0]] * n
    tail = [0 * f[0], 0 * f[0]]
    for k in range(n - 1, -1, -1):
        g[k] = (2 * k + 1) * tail[(k + 1) % 2]
        tail[k % 2] = tail[k % 2] + f[k]
    return g


def _leg_omx2(f):
    """Legendre coefficients of (1 - x^2) f (length + 2)."""
    n = len(f)
    z = 0 * f[0]
    xf = [z] * (n + 1)
    for k, v in enumerate(f):
        xf[k + 1] = xf[k + 1] + v * (k + 1) / (2 * k + 1)
        if k:
            xf[k - 1] = xf[k - 1] + v * k / (2 * k + 1)
    xxf = [z] * (n + 2)
    for k, v in enumerate(xf):
        xxf[k + 1] = xxf[k + 1] + v * (k + 1) / (2 * k + 1)
        if k:
            xxf[k - 1] = xxf[k - 1] + v * k / (2 * k + 1)
    return [(f[k] if k < n else z) - xxf[k] for k in range(n + 2)]


def total_coeffs(f):
    """Plain Legendre coefficients c_0..c_{Lmax+4} of sum_L sum_m f{m}_L P_L^m(x)."""
    n = len(next(iter(f.values())))
    z = 0 * next(iter(f.values()))[0]
    c = [z] * (n + 4)
    for m, v in f.items():
        g = list(v)
        for _ in range(m):
            g = _leg_deriv(g)
        for _ in range(m // 2):
            g = _leg_omx2(g)
        for k, w in enumerate(g):
            c[k] = c[k] + w
    return c


_EXACT = {}


def exact_Q(case, Lmax):
    """{integral: {m: [ {key: q} ]}} over Q for DNX, DNY, DXY, DXX, DA, one engine per case, cached."""
    if case not in _EXACT or _EXACT[case][0] < Lmax:
        K = Field(None)
        eng = AS.Engine(case, K, Lmax + 2)
        _EXACT[case] = (Lmax, {it: coefficients(case, Lmax, K, it, eng) for it in NEW + ('DA',)})
    return _EXACT[case][1]


_SER = {}


def series(case, Lmax, dps):
    """{integral: {0: [...], 2 or 4: [...], 'c': [...]}} as mpf lists, L = 0..Lmax (c to Lmax: c_j uses f_{>= j})."""
    key = (case, Lmax, dps)
    if key not in _SER:
        r = exact_Q(case, Lmax)
        out = {}
        for it, ch in r.items():
            s = {m: AS.numeric(v[:Lmax + 1], case, dps) for m, v in ch.items()}
            with mp.workdps(dps):
                s['cfull'] = total_coeffs({m: s[m] for m in ch})
                s['c'] = s['cfull'][:Lmax + 1]
            out[it] = s
        _SER[key] = out
    return _SER[key]


# ---------------------------------------------------------------------------------------------------------------
# Independent route (no U kernels, no D2-derivatives, no ladder): Gauss's-law multipoles of the A density
# ---------------------------------------------------------------------------------------------------------------

CL = {0: Fraction(1, 5), 2: Fraction(2, 7), 4: Fraction(18, 35)}    # AS's d multipole split of P2(u)^2


def g_unit(l, s, z):
    """Potential of N_d^2 t^4 e^{-2 z t} P_l (unit harmonic weight): AS's incomplete-gamma V_l divided by its c_l."""
    return AS.V_multipole('d', l, s, z) / AS._mf(CL[l])


def quad_form(name, wx, wy, wz):
    """w M w in the A frame (x' = e', y' = e_y, z' = n), straight from the record's table."""
    s3 = mp.sqrt(3)
    return {'DA': (3 * wz * wz - (wx * wx + wy * wy + wz * wz)) / 2, 'DNX': s3 * wz * wx, 'DNY': s3 * wz * wy,
            'DXY': s3 * wx * wy, 'DXX': s3 / 2 * (wx * wx - wy * wy)}[name]


def plm_poly(l, m):
    """Ascending Fraction coefficients of P_l^m(t) = (1 - t^2)^{m/2} d^m P_l/dt^m, m even (no CS phase)."""
    P = [[Fraction(1)], [Fraction(0), Fraction(1)]]
    for n in range(1, l):
        a = [Fraction(0)] + [c * Fraction(2 * n + 1, n + 1) for c in P[n]]
        b = [c * Fraction(n, n + 1) for c in P[n - 1]] + [Fraction(0)] * 2
        P.append([x - y for x, y in zip(a, b)])
    p = P[l][:l + 1]
    for _ in range(m):
        p = [k * p[k] for k in range(1, len(p))] or [Fraction(0)]
    for _ in range(m // 2):
        q = [Fraction(0)] * (len(p) + 2)
        for k, c in enumerate(p):
            q[k] += c
            q[k + 2] -= c
        p = q
    return p


def _pev(c, t):
    v = mp.mpf(0)
    for a in reversed(c):
        v = v * t + AS._mf(a)
    return v


def harmonic_content(name, dps=40):
    """{(l, m): h} with (w M w)^2 = sum h P_l^m(cos th) cos(m ph) on the unit sphere, by Gauss-Legendre in cos th
    (8 nodes) and a 16-point rule in ph: exact for this degree-4 polynomial, no closed form used."""
    with mp.workdps(dps):
        nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(2, mp.mp.prec)     # 6 nodes, exact to 11
        out = {}
        for l in (0, 2, 4):
            for m in range(0, l + 1, 2):
                pl = plm_poly(l, m)
                num = den = mp.mpf(0)
                for t, w in nodes:
                    st = mp.sqrt(1 - t * t)
                    Y0 = _pev(pl, t)
                    for k in range(16):
                        ph = 2 * mp.pi * k / 16
                        Y = Y0 * mp.cos(m * ph)
                        f = quad_form(name, st * mp.cos(ph), st * mp.sin(ph), t) ** 2
                        num += w * f * Y
                        den += w * Y * Y
                if abs(num / den) > mp.mpf(10) ** (-dps + 5):
                    out[(l, m)] = num / den
        return out


def direct_channels(case, Ls, names, dps, deg_in=5, deg_out=5):
    """{(name, L, m): f{m}_L} by quadrature (port of the record's AY0 'direct' check, debug/aay_ap0.py):
    f{m}_L = 2 pi (L-m)!/(L+m)! sum_l h_{lm} int r^2 rho_L(r) int dt P_L^m(t) g_l(s) P_l^m((r t - R)/s), with rho_L
    the BC density's Legendre channel about B taken by quadrature in the B-C angle."""
    H = {nm: harmonic_content(nm) for nm in names}
    with mp.workdps(dps):
        D1, D2, zB, zC, zA = (AS._mf(v) for v in params(case))
        R = D2
        NB, NC = mp.sqrt(zB ** 3 / mp.pi), mp.sqrt(zC ** 3 / mp.pi)
        GL = mp.calculus.quadrature.GaussLegendre(mp.mp)
        n_in, n_out = GL.calc_nodes(deg_in, mp.mp.prec), GL.calc_nodes(deg_out, mp.mp.prec)
        Lm = [(L, m) for L in Ls for m in (0, 2, 4) if m <= L]
        lm = [(l, m) for l in (0, 2, 4) for m in (0, 2, 4) if m <= l]
        polys = {k: plm_poly(*k) for k in set(Lm) | set(lm) | {(L, 0) for L in Ls}}
        acc = {(L, m, l): mp.mpf(0) for (L, m) in Lm for l in (0, 2, 4) if m <= l}
        edges = [mp.mpf(0), D1, D2, D2 + 3, D2 + 10, D2 + 25, D2 + 70]
        for lo, hi in zip(edges[:-1], edges[1:]):
            hh, cc = (hi - lo) / 2, (hi + lo) / 2
            for xo, wo in n_out:
                r = cc + hh * xo
                wr = wo * hh * r * r
                m0 = {L: mp.mpf(0) for L in Ls}
                a_, b_ = abs(r - D1), r + D1
                h2, c2 = (b_ - a_) / 2, (b_ + a_) / 2
                pre = NB * NC * mp.exp(-zB * r)
                for xi, wi in n_in:
                    sc = c2 + h2 * xi
                    tau = (r * r + D1 * D1 - sc * sc) / (2 * r * D1)
                    w = wi * h2 * sc / (r * D1) * pre * mp.exp(-zC * sc)
                    for L in Ls:
                        m0[L] += w * _pev(polys[(L, 0)], tau)
                rhoL = {L: (2 * L + 1) * m0[L] / 2 for L in Ls}
                J = {k: mp.mpf(0) for k in acc}
                a_, b_ = abs(r - R), r + R
                h2, c2 = (b_ - a_) / 2, (b_ + a_) / 2
                for xi, wi in n_in:
                    s = c2 + h2 * xi
                    t = (r * r + R * R - s * s) / (2 * r * R)
                    w = wi * h2 * s / (r * R)
                    cw = (r * t - R) / s
                    gv = {l: g_unit(l, s, zA) for l in (0, 2, 4)}
                    pl = {k: _pev(polys[k], cw) for k in lm}
                    pL = {k: _pev(polys[k], t) for k in Lm}
                    for (L, m, l) in J:
                        J[(L, m, l)] += w * pL[(L, m)] * gv[l] * pl[(l, m)]
                for k in acc:
                    acc[k] += wr * rhoL[k[0]] * J[k]
        out = {}
        for name in names:
            for (L, m) in Lm:
                v = sum(H[name].get((l, m), 0) * acc[(L, m, l)] for l in (0, 2, 4) if (L, m, l) in acc)
                out[(name, L, m)] = 2 * mp.pi * AS._mf(Fraction(factorial(L - m), factorial(L + m))) * v
        return out


# ---------------------------------------------------------------------------------------------------------------
# Toys: coefficients of (x* - x)^alpha in the P_L, P_L^2, P_L^4 bases ((1 - x^2)^{m/2} (x* - x)^alpha in P_L^m:
# an m-fold Legendre antiderivative of the P_L coefficients)
# ---------------------------------------------------------------------------------------------------------------

_TOY4 = {}


def toy(case, alpha, Lmax, dps=60, deg=7):
    """{'P', 'P2', 'P4'}; valid to L = Lmax (AS's toy run 2 further, P4 = two more antiderivatives of P2)."""
    key = (case, alpha, Lmax, dps, deg)
    if key not in _TOY4:
        t = AS.toy(case, alpha, Lmax + 2, dps, deg)
        with mp.workdps(dps):
            _TOY4[key] = {'P': t['P'], 'P2': t['P2'], 'P4': AS._anti_list(AS._anti_list(t['P2']))}
    return _TOY4[key]


x_star, beta_eff, _mf = AS.x_star, AS.beta_eff, AS._mf


def plm_values(Lmax, m, x):
    """P_L^m(x), L = 0..Lmax, CS-free: P_m^m = (2m-1)!! (1-x^2)^{m/2}, (L-m+1) P_{L+1}^m = (2L+1) x P_L^m - (L+m) P_{L-1}^m."""
    out = [mp.mpf(0)] * (Lmax + 1)
    if m > Lmax:
        return out
    df = mp.mpf(1)
    for k in range(1, 2 * m, 2):
        df *= k
    out[m] = df * (1 - x * x) ** (m // 2)
    if m + 1 <= Lmax:
        out[m + 1] = (2 * m + 1) * x * out[m]
    for L in range(m + 1, Lmax):
        out[L + 1] = ((2 * L + 1) * x * out[L] - (L + m) * out[L - 1]) / (L - m + 1)
    return out


# ---------------------------------------------------------------------------------------------------------------
# Shared data for the fast tests: case A1, all five integrals to L = FAST_L over Q (one engine)
# ---------------------------------------------------------------------------------------------------------------

FAST_L, DPS = 51, 450       # c_j uses f_{>= j}: read c at L <= FAST_L - 21 only
SLOW_L, SLOW_DPS = 121, 650


def _fast(dps=DPS):
    return series('A1', FAST_L, dps)


_DIRECT = {}


def _direct(case, Ls, dps):
    key = (case, Ls, dps)
    if key not in _DIRECT:
        _DIRECT[key] = direct_channels(case, Ls, NEW + ('DA',), dps)
    return _DIRECT[key]


# ---------------------------------------------------------------------------------------------------------------
# (0) the exact symmetry identities
# ---------------------------------------------------------------------------------------------------------------


def test_exact_pair_symmetries_over_Q():
    """Key for key over Q, every L <= FAST_L at A1: f0(DNY) = f0(DNX), f2(DNY) = -f2(DNX), f0(DXY) = f0(DXX),
    f4(DXY) = -f4(DXX); DNX, DNY carry exactly m in {0, 2} and DXY, DXX exactly {0, 4} (DA {0}); every f2 (L >= 2) and
    f4 (L >= 4) is nonzero, so the identities are not 0 = 0.  Rejects an out-of-plane factor d_y^2 with the wrong sign
    (f2(DNY) = +f2(DNX)), an x'/y' mix-up that puts m = 2 into DXY, and a DXX quartic (d_x^2 - d_y^2)^2 written as
    d_x^4 - d_y^4 (which breaks f0(DXY) = f0(DXX))."""
    r = exact_Q('A1', FAST_L)
    for nm in NEW + ('DA',):
        assert tuple(sorted(r[nm])) == MSET[nm], (nm, sorted(r[nm]))
    for a, b, m in (('DNX', 'DNY', 2), ('DXX', 'DXY', 4)):
        for L in range(FAST_L + 1):
            assert r[a][0][L] == r[b][0][L], (a, b, 'f0', L)
            assert r[b][m][L] == {k: -v for k, v in r[a][m][L].items()}, (a, b, m, L)
            if L >= m:
                assert r[a][m][L], (a, m, L)


# ---------------------------------------------------------------------------------------------------------------
# (1) per-channel values against the independent quadrature
# ---------------------------------------------------------------------------------------------------------------


def test_harmonic_content_matches_the_record_table():
    """The sphere quadrature of (w M w)^2 that feeds the independent route, against the record's rational table
    (AY0 log): DNX 1/5, 1/7, +1/14, -12/35, +1/35; DXY 1/5, -2/7, 3/35, -1/280; DNY, DXX the m != 0 signs flipped;
    DA 1/5, 2/7, 18/35.  Rejects a wrong w M w in quad_form (a dropped sqrt3, x' and y' swapped, a 1/2 lost on
    DXX), which would otherwise leave the quadrature agreeing with itself."""
    F = Fraction
    want = {'DA': {(0, 0): F(1, 5), (2, 0): F(2, 7), (4, 0): F(18, 35)},
            'DNX': {(0, 0): F(1, 5), (2, 0): F(1, 7), (2, 2): F(1, 14), (4, 0): F(-12, 35), (4, 2): F(1, 35)},
            'DNY': {(0, 0): F(1, 5), (2, 0): F(1, 7), (2, 2): F(-1, 14), (4, 0): F(-12, 35), (4, 2): F(-1, 35)},
            'DXY': {(0, 0): F(1, 5), (2, 0): F(-2, 7), (4, 0): F(3, 35), (4, 4): F(-1, 280)},
            'DXX': {(0, 0): F(1, 5), (2, 0): F(-2, 7), (4, 0): F(3, 35), (4, 4): F(1, 280)}}
    with mp.workdps(40):
        for nm, w in want.items():
            h = harmonic_content(nm)
            assert set(h) == set(w), (nm, sorted(h))
            for k, v in w.items():
                assert abs(h[k] - _mf(v)) < mp.mpf('1e-30'), (nm, k)


@pytest.mark.parametrize('integral', NEW + ('DA',))
def test_channels_match_independent_quadrature_A1(integral):
    """f{m}_L (case A1, L = 2..5, every m) against the direct quadrature: relative 1e-20 at 25 digits (measured
    3e-26 to 9e-24); a channel the engine has as zero is zero in the quadrature (below 1e-30).  The quadrature uses
    the A density's multipole potentials, the A-frame harmonic content of (w M w)^2 and P_L^m projection, with no U
    kernel, no D2-derivative and no ladder.  Rejects a wrong U4/U2/U0 coefficient, a wrong operator ((d M d)^2 or
    d M^2 d), a wrong ladder coefficient, a wrong (L-m)!/(L+m)! and a wrong D2-derivative order."""
    d = _direct('A1', (2, 3, 4, 5), 25)
    S = series('A1', 5, 80)[integral]
    with mp.workdps(40):
        for (nm, L, m), v in d.items():
            if nm != integral:
                continue
            ex = S[m][L] if m in S else mp.mpf(0)
            if ex == 0:
                assert abs(v) < mp.mpf('1e-30'), (nm, L, m, mp.nstr(v, 5))
                continue
            rel = abs(v / ex - 1)
            assert rel < mp.mpf('1e-20'), (nm, L, m, mp.nstr(v, 25), mp.nstr(ex, 25), mp.nstr(rel, 3))


@pytest.mark.slow
def test_channels_match_independent_quadrature_A2():
    """Case A2 (zC != 1), L = 2, 4, 5, all five integrals: relative 1e-20."""
    d = direct_channels('A2', (2, 4, 5), NEW + ('DA',), 25)
    S = series('A2', 5, 80)
    with mp.workdps(40):
        for (nm, L, m), v in d.items():
            ex = S[nm][m][L] if m in S[nm] else mp.mpf(0)
            if ex == 0:
                assert abs(v) < mp.mpf('1e-30'), (nm, L, m)
            else:
                assert abs(v / ex - 1) < mp.mpf('1e-20'), (nm, L, m)


# ---------------------------------------------------------------------------------------------------------------
# (3) the constant set equals AB's 21
# ---------------------------------------------------------------------------------------------------------------

CONST_L = 12


@pytest.mark.parametrize('case', ['A1', 'A2'])
def test_constant_set_is_ABs_21_no_gamma(case):
    """Mod 2^61 - 1, every L up to 12: for each of DNX, DNY, DXY, DXX, f0 (from L = 2) and f2 (from L = 2) or f4
    (from L = 4) each have exactly the constant set of AB at the same case (A1: AB's frozen EXPECTED_A1; A2: AB's
    own engine, independent of the A pair), so their union does too, and no coefficient at any L <= 12 carries
    Euler's gamma.  Rejects a spurious gamma, a constant lost or added by a wrong kernel family or a dropped jump
    term, and a set that changes with L."""
    K = Field(PRIME)
    ab = AB.EXPECTED_A1 if case == 'A1' else AP._ab_set(case)
    assert len(ab) == 21
    eng = AS.Engine(case, K, CONST_L + 2)
    for nm in NEW:
        ch = coefficients(case, CONST_L, K, nm, eng)
        for m, f in ch.items():
            for L in range(CONST_L + 1):
                assert not [k for k in f[L] if k[1] == 'g'], (nm, m, case, L)
            for L in range(max(2, m), CONST_L + 1):
                S = set(f[L])
                assert S == ab, (nm, m, case, L, sorted(map(repr, S ^ ab)))


# ---------------------------------------------------------------------------------------------------------------
# (2) summed series
# ---------------------------------------------------------------------------------------------------------------

# record AY, debug/data/aay_p1p3_<case>.json: the split summed to L = 400 at 3600 digits (first 30 digits)
P1_3600 = {('A1', 'DNX', '1/2'): '0.156032509629553532370024000273',
           ('A1', 'DNY', '1/2'): '0.152091673380763215389996102846',
           ('A1', 'DXY', '1/2'): '0.138205955286622360926082120198',
           ('A1', 'DXX', '1/2'): '0.138361276079340676701841001869',
           ('A1', 'DNX', '-1/3'): '0.147298002172441018699205204328',
           ('A1', 'DNY', '-1/3'): '0.143600436032060860159288817068',
           ('A1', 'DXY', '-1/3'): '0.126702316721690692959268008546',
           ('A1', 'DXX', '-1/3'): '0.126783484396503709323918689674',
           ('A2', 'DNX', '0'): '0.17134199494804779660253362513',
           ('A2', 'DNY', '0'): '0.165789209076208061534090233485',
           ('A2', 'DXY', '0'): '0.1508949659434039139527381987',
           ('A2', 'DXX', '0'): '0.151069948002653869148242774891',
           ('A2', 'DNX', '4/5'): '0.175463108217065647088013719724',
           ('A2', 'DNY', '4/5'): '0.173250334089233638858065421401',
           ('A2', 'DXY', '4/5'): '0.162015767371004248771265114279',
           ('A2', 'DXX', '4/5'): '0.162066393117948367330228779313'}
# record AY, control R2: AB's even-tempered Gaussian fit, n = 34 (debug/data/aay_controls.json, 'referee')
R2_34 = {('A1', 'DNX', '1/2'): 0.15603250962948195, ('A1', 'DNY', '1/2'): 0.15209167338069288,
         ('A1', 'DXY', '1/2'): 0.13820595528656945, ('A1', 'DXX', '1/2'): 0.13836127607928503,
         ('A1', 'DNX', '-1/3'): 0.1472980021727219, ('A1', 'DNY', '-1/3'): 0.14360043603206465,
         ('A1', 'DXY', '-1/3'): 0.12670231672164406, ('A1', 'DXX', '-1/3'): 0.12678348439647485,
         ('A2', 'DNX', '0'): 0.17134199494807917, ('A2', 'DNY', '0'): 0.16578920907627664,
         ('A2', 'DXY', '0'): 0.1508949659434819, ('A2', 'DXX', '0'): 0.1510699480027287,
         ('A2', 'DNX', '4/5'): 0.17546310821713246, ('A2', 'DNY', '4/5'): 0.1732503340893014,
         ('A2', 'DXY', '4/5'): 0.16201576737107562, ('A2', 'DXX', '4/5'): 0.16206639311801938}


def _summed(case, S, xq, Lmax):
    """(split sum sum_L sum_m f{m}_L P_L^m(x), plain-Legendre sum of the full re-expanded total c)."""
    x = _mf(Fraction(xq))
    Pm = {m: plm_values(Lmax + 4, m, x) for m in (0, 2, 4)}
    v = sum(S[m][L] * Pm[m][L] for m in S if m != 'c' and m != 'cfull' for L in range(Lmax + 1))
    vc = sum(c * Pm[0][j] for j, c in enumerate(S['cfull']))
    return v, vc


FAST_TRUNC = mp.mpf('1e-18')     # truncation of the sums at L = 51, A1, x = 1/2: measured 2e-22 to 1.2e-20
SLOW_TRUNC = mp.mpf('1e-29')     # at L = 121 (all 16 points: measured <= 5e-31, the 30 stored digits)


@pytest.mark.parametrize('integral', NEW)
def test_summed_series_A1_matches_referee(integral):
    """Sum_{L <= 51} sum_m f{m}_L P_L^m(1/2) at A1 against the R2 referee (n = 34) within the record's P1 tolerance
    1e-10 (record: <= 2.8e-13), and against the record's 3600-digit sum to L = 400 within FAST_TRUNC.  The full
    re-expanded total c summed as a P_L series equals the split sum to 1e-400 (an exact identity at truncation).
    Rejects an m-part with the wrong sign or dropped, a wrong normalization or kernel, P_L^m evaluated wrongly, and
    a re-expansion with a wrong derivative or (1 - x^2) step."""
    S = _fast()[integral]
    with mp.workdps(DPS):
        v, vc = _summed('A1', S, '1/2', FAST_L)
        key = ('A1', integral, '1/2')
        assert abs(v - mp.mpf(R2_34[key])) < 1e-10, mp.nstr(v, 20)
        assert abs(v - mp.mpf(P1_3600[key])) < FAST_TRUNC, mp.nstr(v, 30)
        assert abs(vc - v) < mp.mpf('1e-400'), mp.nstr(abs(vc - v), 3)


@pytest.mark.slow
def test_summed_series_all_points_slow():
    """All sixteen P1 values (both cases, both angles, four integrals) summed to L = 121: within 1e-10 of R2 and
    SLOW_TRUNC of the 3600-digit values; re-expansion exact.  Engine time: the TIMINGS note at the end."""
    for case, xs_ in (('A1', '1/2'), ('A1', '-1/3'), ('A2', '0'), ('A2', '4/5')):
        Sc = series(case, SLOW_L, SLOW_DPS)
        for integral in NEW:
            with mp.workdps(SLOW_DPS):
                v, vc = _summed(case, Sc[integral], xs_, SLOW_L)
                key = (case, integral, xs_)
                assert abs(v - mp.mpf(R2_34[key])) < 1e-10, (key, mp.nstr(v, 20))
                assert abs(v - mp.mpf(P1_3600[key])) < SLOW_TRUNC, (key, mp.nstr(v, 30))
                assert abs(vc - v) < mp.mpf('1e-600'), key


# ---------------------------------------------------------------------------------------------------------------
# The toys (the windows and the P4 ratios rest on them)
# ---------------------------------------------------------------------------------------------------------------

# record AY, control TA (debug/data/aay_toys.json, case A1, L = 220)
TA_220_A1 = {'ta': '5.454419322977413275881174298476923176997e-46',
             'tb2': '1.827732382033674341779126726029859156661e-51',
             'tb4': '6.449390605723611065720946013566367385285e-57',
             'tc': '-4.464168549530084167346853526062353167264e-48'}


def test_toys_reproduce_record_amplitude_toys_A1():
    """The in-test toys at L = 220 against the record's control TA (ta: P_L of (x* - x)^{7/2}; tb2, tb4: P_L^2,
    P_L^4 coefficients of (1 - x^2)^{m/2} (x* - x)^{7/2}; tc: P_L of (x* - x)^{9/2}), relative 1e-25; and the P_L^4
    coefficients against a direct projection with (L-4)!/(L+4)! at L = 6, 9.  Rejects a wrong number of
    antiderivatives for P4, a wrong antiderivative recurrence, and a wrong x*."""
    t7 = toy('A1', Fraction(7, 2), 220, dps=100, deg=8)
    t9 = toy('A1', Fraction(9, 2), 220, dps=100, deg=8)
    with mp.workdps(100):
        for name, v in (('ta', t7['P'][220]), ('tb2', t7['P2'][220]), ('tb4', t7['P4'][220]), ('tc', t9['P'][220])):
            assert abs(v / mp.mpf(TA_220_A1[name]) - 1) < mp.mpf('1e-25'), (name, mp.nstr(v, 20))
    with mp.workdps(40):
        xs = _mf(x_star('A1'))
        t = toy('A1', Fraction(7, 2), 20, dps=40)
        for L in (6, 9):
            pl = plm_poly(L, 4)
            proj = mp.quad(lambda x: (1 - x * x) ** 2 * (xs - x) ** mp.mpf(3.5) * _pev(pl, x), [-1, 1])
            proj *= (2 * L + 1) / mp.mpf(2) * mp.factorial(L - 4) / mp.factorial(L + 4)
            assert abs(proj / t['P4'][L] - 1) < mp.mpf('1e-25'), L


# ---------------------------------------------------------------------------------------------------------------
# (4) decay: 7/2 for the in-plane totals and the parts, 9/2 for the out-of-plane totals
# ---------------------------------------------------------------------------------------------------------------

ALPHAS = [Fraction(k, 2) for k in (3, 5, 7, 9, 11, 13)]
BASIS = {0: 'P', 2: 'P2', 4: 'P4', 'c': 'P'}


def window(case, L, alpha, basis):
    """(lo, hi): the midpoints between the toy beta_eff(L) of alpha and those of alpha + 1 and alpha - 1."""
    t = {a: beta_eff(toy(case, a, L + 2)[basis], case, L) for a in ALPHAS}
    i = ALPHAS.index(alpha)
    return (t[ALPHAS[i]] + t[ALPHAS[i + 1]]) / 2, (t[ALPHAS[i]] + t[ALPHAS[i - 1]]) / 2


EXPECT = {('DNX', 0): Fraction(7, 2), ('DNX', 2): Fraction(7, 2), ('DNX', 'c'): Fraction(7, 2),
          ('DNY', 'c'): Fraction(9, 2),
          ('DXX', 0): Fraction(7, 2), ('DXX', 4): Fraction(7, 2), ('DXX', 'c'): Fraction(7, 2),
          ('DXY', 'c'): Fraction(9, 2)}
FAST_DECAY = [('DNX', 0), ('DNX', 2), ('DNX', 'c'), ('DNY', 'c'), ('DXX', 4), ('DXY', 'c')]


def _check_decay(case, L, S, which):
    for (integral, n) in which:
        lo, hi = window(case, L, EXPECT[(integral, n)], BASIS[n])
        b = beta_eff(S[integral][n], case, L)
        assert lo < b < hi, (case, L, integral, n, b, (lo, hi))


def _gap(S, case, L):
    """beta_eff(L) of each out-of-plane total minus that of its in-plane partner."""
    return (beta_eff(S['DNY']['c'], case, L) - beta_eff(S['DNX']['c'], case, L),
            beta_eff(S['DXY']['c'], case, L) - beta_eff(S['DXX']['c'], case, L))


def test_decay_branches_A1_L30():
    """A1, L = 30 (c read 21 below the truncation; equal at L <= 51 and L <= 121 to 1e-4).  Between the toy
    midpoints: DNX f0 and total c on (x* - x)^{7/2}; DNX f2 on 7/2 in the P_L^2 basis; DXX f4 on 7/2 in the P_L^4
    basis; DNY and DXY totals on the 9/2 branch.  Measured: DNX f0 -4.965, c -5.006, f2 -7.102, DXX f4 -10.569,
    DNY c -5.760, DXY c -6.875.  (DXX's f0 and c at L = 30, -5.955 and -6.170, still carry the large 9/2 term of
    the m = 4 split and enter the 7/2 window only at L ~ 60-80: the slow test.)  The cancellation is asserted
    directly too: each out-of-plane total decays faster than its in-plane partner by more than half a power of L,
    beta(DNY c) - beta(DNX c) < -1/2 and beta(DXY c) - beta(DXX c) < -1/2 (measured -0.754, -0.705; -1
    asymptotically, while no cancellation gives 0).  Signs for 20 <= L <= 30 (record, L >= 20): DNX f0 -, f2 +,
    c -; DNY f0 -, f2 -, c +; DXY f0 +, f4 -, c -; DXX f0 +, f4 +, c +.  Precision: every beta at dps and dps + 100
    agrees to 1e-12.  Rejects "no cancellation" (f_m of DNY or DXY with the sign of DNX or DXX), a dropped m-part
    in the total, and a 7/2 branch for an out-of-plane total."""
    S, S_hi = _fast(), _fast(DPS + 100)
    for it, ch in S.items():
        for n in ch:
            if n != 'cfull':
                assert abs(beta_eff(ch[n], 'A1', 30) - beta_eff(S_hi[it][n], 'A1', 30)) < 1e-12, (it, n)
    _check_decay('A1', 30, S, FAST_DECAY)
    gN, gX = _gap(S, 'A1', 30)
    assert gN < -0.5 and gX < -0.5, (gN, gX)
    sg = {('DNX', 0): -1, ('DNX', 2): 1, ('DNX', 'c'): -1, ('DNY', 0): -1, ('DNY', 2): -1, ('DNY', 'c'): 1,
          ('DXY', 0): 1, ('DXY', 4): -1, ('DXY', 'c'): -1, ('DXX', 0): 1, ('DXX', 4): 1, ('DXX', 'c'): 1}
    for (it, n), s in sg.items():
        assert all(v * s > 0 for v in S[it][n][20:31]), (it, n)


@pytest.mark.slow
@pytest.mark.parametrize('case', ['A1', 'A2'])
def test_decay_branches_L80_L100_slow(case):
    """Both cases, L = 80 and 100 (series to L = 121): every part and total in its toy window, DXX's f0 and c
    included (A1 L = 80: DXX c -4.852 against the window's lower edge -4.902, the closest margin), and both
    out-of-plane gaps below -1/2 (measured -0.87 to -0.97)."""
    S = series(case, SLOW_L, SLOW_DPS)
    for L in (80, 100):
        _check_decay(case, L, S, list(EXPECT))
        gN, gX = _gap(S, case, L)
        assert gN < -0.5 and gX < -0.5, (case, L, gN, gX)


# ---------------------------------------------------------------------------------------------------------------
# (5) P4: R_2 = (f0/ta)/(f2/tb2) of DNY tends to +(x*^2 - 1), R_4 = (f0/ta)/(f4/tb4) of DXY to -(x*^2 - 1)^2
# ---------------------------------------------------------------------------------------------------------------


def _ratios(case, S, L, tmax):
    """(r2, r4, R_2(DNX) / R_2(DNY) + 1, R_4(DXX) / R_4(DXY) + 1): r2 = R_2(DNY) / (x*^2 - 1),
    r4 = R_4(DXY) / (-(x*^2 - 1)^2); ta, tb the 7/2 toys in the P_L and P_L^m bases."""
    t = toy(case, Fraction(7, 2), tmax)
    xs = _mf(x_star(case))
    e = xs * xs - 1

    def R(it, m):
        return (S[it][0][L] / t['P'][L]) / (S[it][m][L] / t['P' + str(m)][L])
    return (float(R('DNY', 2) / e), float(R('DXY', 4) / (-e * e)),
            float(R('DNX', 2) / R('DNY', 2) + 1), float(R('DXX', 4) / R('DXY', 4) + 1))


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


def test_P4_amplitude_ratios_A1():
    """A1, L = 30, 40, 50 (f0 and f_m are local in L, so no truncation enters).  r2_L = R_2 / (x*^2 - 1) of DNY:
    measured 2.150, 1.767, 1.571 (1.20 at L = 120, 1.104 at 220); r4_L = R_4 / (-(x*^2 - 1)^2) of DXY: 3.988,
    2.715, 2.178 (1.35 at 120, 1.174 at 220).  The leading ratios are exactly +(x*^2 - 1) and -(x*^2 - 1)^2 and an
    f0-only 9/2 term adds a 1/L excess of relative size m x*/(x*^2 - 1), so both approach 1 from above, R_4 more
    slowly.  Asserted: (i) r2, r4 > 1 at each L (so R_2 > 0 and R_4 < 0) with r2 < 2.5 and r4 < 5; (ii) both
    decrease with L; (iii) the quadratic-in-1/L extrapolation of r2 is within 0.1 of 1 (measured 1.045) and that of
    r4 lies in (1, 2) (measured 1.727; 1.024 from L = 80, 100, 120, the slow test); (iv) R(DNX) = -R(DNY) and
    R(DXX) = -R(DXY) to 1e-12.  Rejects R_4 tending to -(x*^2 - 1) (r4 about 7x larger), R_2 or R_4 with the wrong
    sign, an f_m scaled by 2 or 1/2 (limits 0.5 or 2), and a ratio that does not move with L."""
    S = _fast()
    with mp.workdps(DPS):
        rs = [_ratios('A1', S, L, FAST_L) for L in (30, 40, 50)]
    r2, r4 = [r[0] for r in rs], [r[1] for r in rs]
    assert all(1 < v < 2.5 for v in r2) and all(1 < v < 5 for v in r4), (r2, r4)
    assert r2[0] > r2[1] > r2[2] and r4[0] > r4[1] > r4[2], (r2, r4)
    assert abs(_richardson((30, 40, 50), r2) - 1) < 0.1, r2
    assert 1 < _richardson((30, 40, 50), r4) < 2, r4
    for r in rs:
        assert abs(r[2]) < 1e-12 and abs(r[3]) < 1e-12, r


@pytest.mark.slow
@pytest.mark.parametrize('case', ['A1', 'A2'])
def test_P4_amplitude_ratios_slow(case):
    """Both cases, L = 80, 100, 120: r2 in (1, 1.5) and r4 in (1, 2), decreasing; quadratic extrapolation of r2
    within 0.01 of 1 (measured A1 1.0025, A2 1.0017) and of r4 within 0.05 (A1 1.0235, A2 1.0123)."""
    S = series(case, SLOW_L, SLOW_DPS)
    with mp.workdps(SLOW_DPS):
        rs = [_ratios(case, S, L, SLOW_L) for L in (80, 100, 120)]
    r2, r4 = [r[0] for r in rs], [r[1] for r in rs]
    assert all(1 < v < 1.5 for v in r2) and all(1 < v < 2 for v in r4), (r2, r4)
    assert r2[0] > r2[1] > r2[2] and r4[0] > r4[1] > r4[2], (r2, r4)
    assert abs(_richardson((80, 100, 120), r2) - 1) < 0.01, r2
    assert abs(_richardson((80, 100, 120), r4) - 1) < 0.05, r4


# TIMINGS (2026-11-16, this machine): fast tests 48 s in all (A1 engine for all five integrals to L = 51 ~25 s, the
# independent quadrature ~12 s, toys ~5 s, constant sets mod p ~2.5 s).  Slow tests 8.8 min in all (--slow, run
# alongside the fire tests): the shared engines to L = 121 for both cases land in the first slow test (A1 ~280 s,
# A2 ~210 s alone); A2 quadrature 11 s; decay 9-10 s per case; P4 under 1 s per case.
#
# FIRE-TEST PLANTS (python debug/qa/fire_test.py tests/test_paper59_aa_bc_a_side_d_mneq0.py --plant "OLD=>NEW" -k <sel>)
# 2026-11-16, written and fired as a separate step after the guard; 24 of 24 fired.
#     e_y factor sign (d_y^2 -> d_x^2)   DY2 = {... -Fraction(1, 4) ...} -> all + signs              -k "exact_pair_symmetries or decay_branches_A1"   FIRED
#     DNY quartic = DNX quartic          DNY 'U4': _pmul(_pmul(DY2, DZ2), ...) -> DX2               -k decay_branches_A1   FIRED (no cancellation)
#     DXX quartic as d_x^4 - d_y^4       "return {'U4': _pmul(_pmul(d, d), {(0, 0, 0): q})" -> _padd(DX2^2, DY2^2, w=[q, -q])   -k exact_pair_symmetries   FIRED
#     m = 4 part dropped from total      total_coeffs loop skips m == 4                              -k decay_branches_A1   FIRED
#     m = 2 part dropped from total      total_coeffs loop skips m == 2                              -k decay_branches_A1   FIRED
#     U4 kernel coefficient              "5: z / 72}=>5: z / 71}"                                    -k quadrature_A1   FIRED
#     U0 kernel coefficient              "3: z ** 3 / 9}=>3: z ** 3 / 8}"                            -k quadrature_A1   FIRED
#     ladder lowering sign               "Aj, -d * (l - am + 1) * (l - am + 2)" -> "+d"              -k quadrature_A1   FIRED
#     ladder raising sign                "_fadd(out, l - 1, nm, Bj, -d)" -> "d"                      -k quadrature_A1   FIRED
#     (L-m)!/(L+m)! wrong                "factorial(L + m))" -> "factorial(L + m - 1))"              -k quadrature_A1   FIRED
#     wrong R power in jet assembly      "K.pw(D2, -p)" -> "K.pw(D2, 1 - p)"                         -k quadrature_A1   FIRED
#     f_m scaled by 2                    "Fraction(factorial(L - m), ..." -> "Fraction(2 * factorial(L - m), ..."   -k P4_amplitude_ratios_A1   FIRED
#     spurious Euler gamma               coefficients adds (Fraction(0), 'g', None): K(1)            -k constant_set   FIRED
#     DNY total on the 7/2 toy           "('DNY', 'c'): Fraction(9, 2)," -> Fraction(7, 2)           -k decay_branches_A1   FIRED
#     DXY total on the 7/2 toy           "('DXY', 'c'): Fraction(9, 2)}" -> Fraction(7, 2)           -k decay_branches_A1   FIRED
#     R_4 limit -(x*^2 - 1)              "R('DXY', 4) / (-e * e)" -> "/ (-e)"                        -k P4_amplitude_ratios_A1   FIRED
#     R_2 limit with the wrong sign      "R('DNY', 2) / e" -> "/ -e"                                 -k P4_amplitude_ratios_A1   FIRED
#     P4 toy with two antiderivatives    "AS._anti_list(AS._anti_list(t['P2']))" -> one              -k "toys_reproduce or P4_amplitude_ratios_A1"   FIRED
#     re-expansion (1 - x^2) sign        "- xxf[k] for k" -> "+ xxf[k] for k"                        -k summed_series_A1   FIRED
#     P_L^m recurrence start             "(2 * m + 1) * x * out[m]" -> "(2 * m + 2)"                 -k summed_series_A1   FIRED
#     independent route: DXY form        "'DXY': s3 * wx * wy," -> "s3 * wx * wx,"                   -k "harmonic or quadrature_A1"   FIRED
#     independent route: multipole split "4: Fraction(18, 35)}" -> "Fraction(18, 36)}" in CL         -k quadrature_A1   FIRED
#     referee off by 2e-10               "('A1', 'DXY', '1/2'): 0.13820595528656945" -> ...548...    -k summed_series_A1   FIRED
#     DNY f2 sign in the sign table      "('DNY', 2): -1, ('DNY', 'c'): 1," -> "('DNY', 2): 1, ..."  -k decay_branches_A1   FIRED
#   Not a plant: dropping the jump terms at r = D2 -- identically zero for every derivative used (record AY0, L <= 8);
#   the no-cancellation plants act on the operator, since the f_m sign of DNY/DXY is fixed by e_y alone.
