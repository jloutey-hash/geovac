"""Exact constant split of the bond-angle series of (AA|BC) with an m != 0 3d pair on A (record AY; inventory row 8,
last A-side item).  Geometry of AB/AS: B at the origin, C = (0, 0, D1), A = D2 (sin g, 0, cos g), x = cos g; the BC
pair is 1s_B 1s_C.  A-frame: n = A/|A| (z'), e' = (cos g, 0, -sin g) (x', in plane, perpendicular to B-A),
e_y = (0, 1, 0) (y').  Real normalized Slater 3d on A, N_d (w M w) e^{-zA s}, tr M^2 = 3/2 (N_d as record AS's DA):

    DNX = (dnx dnx | 1s_B 1s_C),  w M w = sqrt3 (n.w)(e'.w)
    DNY = (dny dny | 1s_B 1s_C),  w M w = sqrt3 (n.w)(e_y.w)
    DXY = (dxy dxy | 1s_B 1s_C),  w M w = sqrt3 (e'.w)(e_y.w)
    DXX = (dxx dxx | 1s_B 1s_C),  w M w = (sqrt3/2) ((e'.w)^2 - (e_y.w)^2)
    DA  = record AS's 3d-sigma pair, w M w = (3 (n.w)^2 - w^2)/2     (recomputed here as a control)

KERNELS (Fourier route, exact; derivation re-run by aay_ap0.check_family).  For ANY traceless M with tr M^2 = 3/2,
the potential of N_d^2 (w M w)^2 e^{-2 zA s} is
    V_M(w) = (d M d)^2 U4(s) + (d M^2 d) U2(s) + tr(M^2) U0(s),     d = grad_w,
with universal radial U4, U2, U0 in the engine family {1/s} + {d^k Y/dmu^k, k <= 5}  (U4 = (4/9) x AS's DA U4).

SPLIT (the engine's convention).  Expand the BC density about B, rho_BC = sum_L rho_L(r) P_L(cos theta_r).  Channel L
contributes a function of A that is the single lab harmonic G_L(|A|) P_L(z.A_hat) for each radial kernel U.  Rewrite
P_L(z.A_hat) in the A-frame (addition theorem, lab z at polar angle g, azimuth pi in that frame), apply the constant
A-frame operator O_M(d_A) and evaluate at the pole A_hat = z'.  Every squared density here is even in x' and in y'
separately, so only even m in {0, 2, 4} survive and

    I(x) = sum_L [ f0_L P_L(x) + f2_L P_L^2(x) + f4_L P_L^4(x) ],    P_L^m = (1 - x^2)^{m/2} d^m P_L/dx^m (no CS phase),

    f{m}_L = (L-m)!/(L+m)! * sum_{sigma = +-1} [ O_M  (G_L(R) P_L^m(cos theta') e^{i sigma m phi'}) ]_{pole}   (m > 0)
    f0_L   =                 [ O_M (G_L(R) P_L(cos theta')) ]_{pole}.

Each f{m}_L is a finite combination of G_L^{(j)}(D2) D2^{-p}, j <= 4, of the SAME L (local in L).  m = 0 is the A-frame
azimuthally averaged part; m = 2 occurs for DNX, DNY (f2(DNY) = -f2(DNX), f0 equal); m = 4 for DXY, DXX
(f4(DXY) = -f4(DXX), f0 equal); DA has m = 0 only.  The total as a plain Legendre series: total_coeffs() below
(c = f0 + leg[(1-x^2) D^2 f2] + leg[(1-x^2)^2 D^4 f4], D = d/dx on Legendre vectors); exact at truncation.

Ladder rules (CS-free, verified by sympy; A_l f = f' - l f/R, B_l f = f' + (l+1) f/R), on f P_l^{|m|} e^{i m phi}:
    d_z   -> [(l-|m|+1) A_l f P_{l+1}^{|m|} + (l+|m|) B_l f P_{l-1}^{|m|}] / (2l+1)
    raise -> [A_l f P_{l+1}^{|m|+1} - B_l f P_{l-1}^{|m|+1}] / (2l+1)          (d_+ for m >= 0, d_- for m <= 0)
    lower -> [-(l-|m|+1)(l-|m|+2) A_l f P_{l+1}^{|m|-1} + (l+|m|)(l+|m|-1) B_l f P_{l-1}^{|m|-1}] / (2l+1)
Output convention as aaa_exact: list over L of {key: q}, value = 2 pi N_B N_C sum q c  (twopiN).
"""
from __future__ import annotations

import sys
from fractions import Fraction
from math import factorial
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import aaa_exact as A  # noqa: E402
import na_exact as E  # noqa: E402

CASES, params = A.CASES, A.params
NEW = ('DNX', 'DNY', 'DXY', 'DXX')


def fam_u(zA: Fraction) -> dict:
    """universal d-pair kernels (N_d^2 inside), keys as aaa_exact: 'NA' = 1/s, k = d^k Y/dmu^k."""
    z = Fraction(zA)
    U4 = {'NA': Fraction(5, 3) / z ** 4, 0: -Fraction(5, 3) / z ** 4, 1: Fraction(965, 384) / z ** 3,
          2: -Fraction(325, 192) / z ** 2, 3: Fraction(185, 288) / z, 4: Fraction(-5, 36), 5: z / 72}
    U2 = {'NA': Fraction(8, 3) / z ** 2, 0: -Fraction(8, 3) / z ** 2, 1: Fraction(31, 8) / z, 2: Fraction(-29, 12),
          3: Fraction(7, 9) * z, 4: -z ** 2 / 9}
    U0 = {'NA': Fraction(2, 3), 0: Fraction(-2, 3), 1: Fraction(11, 12) * z, 2: -z ** 2 / 2, 3: z ** 3 / 9}
    return {'U4': U4, 'U2': U2, 'U0': U0}


# ---------------------------------------------------------------- operators in (d_+, d_-, d_z)
# a polynomial is {(a, b, c): Fraction} for d_+^a d_-^b d_z^c.  d_x = (d_+ + d_-)/2, d_y^2 = -(d_+ - d_-)^2/4.


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


def operators(name: str) -> dict:
    """{'U4': quartic (d M d)^2, 'U2': quadratic (d M^2 d), 'U0': tr M^2}, A-frame (x' = e', y' = e_y, z' = n).
    'PY', 'PXp', 'PA' (p pairs, kernels of aaa_exact.fam_p: U2 with (e.d)^2, U0) are the record-AS positive controls."""
    q = Fraction(3, 4)
    if name == 'DA':      # M = diag(-1/2, -1/2, 1):  d M d = (3 d_z^2 - Lap)/2,  M^2 = diag(1/4, 1/4, 1)
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
    if name in ('PY', 'PXp', 'PA'):
        return {'U2': {'PY': DY2, 'PXp': DX2, 'PA': DZ2}[name], 'U0': {(0, 0, 0): Fraction(1)}}
    raise ValueError(name)


# ---------------------------------------------------------------- ladder on single-source jets
# a jet is {(j, p): Fraction} = coef * R^{-p} G^{(j)}(R); a field is {(l, m): jet} for jet * P_l^{|m|} e^{i m phi}


def _jd(jet):
    out = {}
    for (j, p), c in jet.items():
        if p:
            out[(j, p + 1)] = out.get((j, p + 1), 0) - p * c
        out[(j + 1, p)] = out.get((j + 1, p), 0) + c
    return out


def _jlin(*terms):
    """sum of a * jet (a Fraction), plus shifts: terms are (a, jet, extra_p)."""
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


def step(field: dict, op: str) -> dict:
    """op in {'+', '-', 'z'}."""
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


def pole(L: int, m: int, poly: dict) -> dict:
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


def channel_jets(name: str, L: int) -> dict:
    """{m: {kernel: jet}} for m in {0, 2, 4} (m <= L), the f{m}_L of the module docstring."""
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


# ---------------------------------------------------------------- engine


class Engine(A.Engine):
    """aaa_exact.Engine (same densities, kernel jets, jump terms) with the universal d kernels."""

    def G_kernel(self, kern: str, L: int, nmax: int, which: str = 'd') -> list:
        fam = (fam_u(self.zA) if which == 'd' else A.fam_p(self.zA))[kern]
        return self.G(kern + '_' + which, fam, L, nmax)


NMAX = {'U4': 4, 'U2': 2, 'U0': 0}


def coefficients(case: str, Lmax: int, K, integral: str, eng: Engine = None, params_override=None) -> dict:
    """{m: list over L of {key: q}} for m in the integral's set ({0, 2} DNX/DNY, {0, 4} DXY/DXX, {0} DA; p-pair
    controls PY/PXp {0, 2}, PA {0})."""
    eng = eng or Engine(case, K, Lmax, params_override)
    which = 'p' if integral in ('PY', 'PXp', 'PA') else 'd'
    D2 = eng.D2
    out = {}
    for L in range(Lmax + 1):
        cj = channel_jets(integral, L)
        for m, rec in cj.items():
            acc = {}
            for kern, jet in rec.items():
                nm = max(j for j, _ in jet)
                Gs = eng.G_kernel(kern, L, NMAX[kern] if which == 'd' else {'U2': 2, 'U0': 0}[kern], which)
                assert nm <= len(Gs) - 1
                for (j, p), c in jet.items():
                    a = K(c) * K.pw(D2, -p)
                    for k, v in Gs[j].items():
                        E.add(acc, k, a * v)
            out.setdefault(m, [dict() for _ in range(Lmax + 1)])[L] = {k: v for k, v in acc.items() if v != 0}
    return out


def twopiN(case: str):
    return A.twopiN(case)


# ---------------------------------------------------------------- total as a plain Legendre series


def _leg_deriv(f: list) -> list:
    """Legendre coefficients of f' : g_k = (2k+1) sum_{j > k, j - k odd} f_j."""
    n = len(f)
    g = [0 * f[0]] * n
    tail = [0 * f[0], 0 * f[0]]
    for k in range(n - 1, -1, -1):
        g[k] = (2 * k + 1) * tail[(k + 1) % 2]
        tail[k % 2] = tail[k % 2] + f[k]
    return g


def _leg_omx2(f: list) -> list:
    """Legendre coefficients of (1 - x^2) f (length grows by 2)."""
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


def total_coeffs(f: dict) -> list:
    """f = {m: list of numbers (mp or Fraction) over L = 0..Lmax}, m in {0, 2, 4}; returns the plain Legendre
    coefficients c_0..c_{Lmax+4} of sum_L sum_m f{m}_L P_L^m(x) (exact polynomial identity at truncation)."""
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


if __name__ == '__main__':
    import time
    K = E.Field(None)
    t0 = time.time()
    for name in NEW + ('DA',):
        print(name, {L: {m: sorted(r) for m, r in channel_jets(name, L).items()} for L in (0, 2, 5)})
    print(round(time.time() - t0, 2), 's')
