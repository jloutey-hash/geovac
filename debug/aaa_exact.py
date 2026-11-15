"""Exact constant split of the bond-angle series of two (AA|BC) integrals with a non-s pair on A
(B at the origin, C = (0, 0, D1), A = D2 (sin g, 0, cos g), x = cos g; the BC pair stays 1s_B 1s_C):

    PY  = (2p_y,A 2p_y,A | 1s_B 1s_C)       p_pi at A perpendicular to the ABC plane
    PXp = (2p_x',A 2p_x',A | 1s_B 1s_C)     p_pi at A in the plane, x' = (cos g, 0, -sin g)  (comes free)
    DA  = (3d-sigma_A 3d-sigma_A | 1s_B 1s_C)  sigma axis n = A/|A|, density N_d^2 h2(w)^2 e^{-2 zA s}

Imports aapc_exact / aabc_exact / na_exact unchanged.

KERNELS (debug/aaa_symbolic.py, Fourier route, exact):  V_e(w) = sum_j (e.grad)^{2j} U_{2j}(|w|) with
    p:  U2 = aapc_exact's U, U0 = aapc_exact's W0 (reproduced key for key),
    d:  U4, U2, U0 in the family {NA} + {d^k Y, k <= 5}.
So with I_U(A) = int rho_BC U(|r - A|) = sum_L G_L(|A|) P_L(cos theta_A)  (G_L = the engine's integral at D2 = |A|):

    DA_L  = G4_L'''' + G2_L'' + W_L                         (directional derivatives along n = A/|A| are d/dD2)
    PY    = d_y^2 I_U2 + I_W,   d_y^2 = (Lap - d_z^2)/2 - (d_+^2 + d_-^2)/4,
            PY(x) = sum_L [a_L P_L(x) + b_L P_L^2(x)],  P_L^2 = (1 - x^2) P_L''  (no Condon-Shortley phase),
            a_L = [(Lap - d_z^2) Phi]_L / 2 + W_L,   b_L = -beta_L / 2,   d_+^2 Phi = e^{2 i phi} sum beta_L P_L^2.
    PXp   = (Lap - d_n^2 - d_y^2) I_U2 + I_W:  a'_L = [Lap Phi]_L - G_L'' - [(Lap - d_z^2) Phi]_L/2 + W_L,  b'_L = +beta_L/2.
The (a_L, b_L) split is the 3D harmonic one of PY(A) as a function on the A-sphere (m = 0 and m = +-2 about z, the
p_y direction held fixed in the lab); a function of x alone does not fix it.  Formulas used (checked by sympy):
    Lap (f P_l) = (f'' + 2 f'/R - l(l+1) f/R^2) P_l
    d_z (f P_l) = (l+1)/(2l+1) (A_l f) P_{l+1} + l/(2l+1) (B_l f) P_{l-1}
    d_+ (f P_l^m e^{i m phi}) = e^{i(m+1)phi}/(2l+1) [(A_l f) P_{l+1}^{m+1} - (B_l f) P_{l-1}^{m+1}]
    A_l f = f' - l f/R,  B_l f = f' + (l+1) f/R.
Each a_L, b_L is a finite combination of G_{L'}^{(j)}, L' in {L-2, L, L+2}, j <= 2, with coefficients in Q[1/D2].

D2-DERIVATIVES WITH JUMP TERMS (general n):  F(D2) = int_0^{D2} a(r, D2) dr + int_{D2}^oo b(r, D2) dr gives
    F^(n) = int d^n_D (a | b) + sum_{i=0}^{n-1} C(n, i+1) d_r^i [rho(r) (d_D^{n-1-i} K_near - d_D^{n-1-i} K_far)]_{r=D2},
(rho = -2 r^2 e^{-zB r} dens).  n = 2 reduces to aapc_exact.boundary_D2 (checked key for key).

Output convention as aapc_exact: list over L of {key: q}, value = 2 pi N_B N_C sum q c.
"""
from __future__ import annotations

import sys
from fractions import Fraction
from math import comb
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import aabc_exact as X  # noqa: E402
import aapc_exact as P  # noqa: E402
import na_exact as E  # noqa: E402

CASES, params = X.CASES, X.params

# family coefficients as functions of zA (from aaa_symbolic; hard-coded so the engine does not need sympy).
# key 'NA' = 1/s; integer k = d^k Y / dmu^k = (-1)^k s^{k-1} e^{-mu s}.


def fam_p(zA: Fraction) -> dict:
    z = zA
    U2 = {'NA': Fraction(3, 2) / z ** 2, 0: -Fraction(3, 2) / z ** 2, 1: Fraction(33, 16) / z, 2: Fraction(-9, 8),
          3: z / 4}
    U0 = {'NA': Fraction(1), 0: Fraction(-1), 1: Fraction(5, 4) * z, 2: -z * z / 2}
    V0 = {'NA': Fraction(1), 0: Fraction(-1), 1: Fraction(3, 2) * z, 2: -z * z, 3: z ** 3 / 3}   # spherical average
    return {'U2': U2, 'U0': U0, 'V0': V0}


def fam_d(zA: Fraction) -> dict:
    z = zA
    U4 = {'NA': Fraction(15, 4) / z ** 4, 0: -Fraction(15, 4) / z ** 4, 1: Fraction(2895, 512) / z ** 3,
          2: -Fraction(975, 256) / z ** 2, 3: Fraction(185, 128) / z, 4: Fraction(-5, 16), 5: z / 32}
    U2 = {'NA': 2 / z ** 2, 0: -2 / z ** 2, 1: Fraction(151, 64) / z, 2: Fraction(-23, 32), 3: -Fraction(17, 48) * z,
          4: z ** 2 / 3, 5: -z ** 3 / 12}
    U0 = {'NA': Fraction(1), 0: Fraction(-1), 1: Fraction(137, 96) * z, 2: -Fraction(41, 48) * z ** 2,
          3: Fraction(7, 24) * z ** 3, 4: -z ** 4 / 9, 5: z ** 5 / 18}
    return {'U4': U4, 'U2': U2, 'U0': U0}


# ---------------------------------------------------------------- kernels (general mu-order; nD2 jet)


def kernel_jet(L: int, near: bool, coeffs: dict, it, kt, K, mu: Fraction, D2: Fraction, nmax: int) -> list:
    """[kernel D2-derivative of order j evaluated, as {(e0, b): {n: c}}] for j = 0..nmax; coeffs {'NA'|k: field}.
    Same factored construction as aapc_exact.kernel_eval, any mu-order k."""
    outs = [dict() for _ in range(nmax + 1)]
    if 'NA' in coeffs:
        na = {k: coeffs['NA'] * v for k, v in P.kernel_mono(L, near, 'NA', it, kt, K).items()}
        for j in range(nmax + 1):
            if j:
                na = P.d_D2(na, K)
            for key, v in P.evaluate(na, mu, D2, K).items():
                for n, c in v.items():
                    X._put(outs[j], key[0], key[1], n, c)
    ords = [k for k in coeffs if k != 'NA']
    if not ords:
        return outs
    amax = max(ords)
    G0, H0 = P._factors(L, near, it, kt, K)
    Gs, Hs = [G0], [H0]
    for _ in range(amax):
        Gs.append(P._dG_mu(Gs[-1], K))
        Hs.append(P._dH_mu(Hs[-1], K))
    pmu, pD2 = {}, {}

    def pw(cache, q, k):
        if k not in cache:
            cache[k] = K.pw(q, k)
        return cache[k]
    Hv = []
    for h in Hs:
        v = {}
        for (pa, n, s), c in h.items():
            P._addm(v, (s, n), c * pw(pmu, mu, pa))
        Hv.append(v)
    for j in range(nmax + 1):
        if j:
            Gs = [P._dG_D2(g, K) for g in Gs]
        Gv = []
        for g in Gs:
            v = {}
            for (pa, pd, t0), c in g.items():
                P._addm(v, t0, c * pw(pmu, mu, pa) * pw(pD2, D2, pd))
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
                        P._addm(acc, (t0, s, n), wg * hv)
        for (t0, s, n), v in acc.items():
            X._put(outs[j], t0 * mu * D2, s * mu, n, v)
    return outs


def boundary_gen(dens_outer: dict, jet_near: list, jet_far: list, n: int, D2: Fraction, zB: Fraction, K) -> dict:
    """sum_{i=0}^{n-1} C(n, i+1) d_r^i [rho (K_near^(n-1-i) - K_far^(n-1-i))]_{r = D2}, rho = -2 r^2 e^{-zB r} dens."""
    res = {}
    for i in range(n):
        j = n - 1 - i
        w = K(comb(n, i + 1))
        jump = {}
        for src, sg in ((jet_near[j], 1), (jet_far[j], -1)):
            for (e0, b), poly in src.items():
                for m, v in poly.items():
                    X._put(jump, e0, b, m, K(sg) * v)
        for (e0d, bd), pd in dens_outer.items():
            for (e0k, bk), pk in jump.items():
                Pc = X._conv(pd, pk, K)
                bb = bd + bk - zB
                e0 = e0d + e0k
                for m, v in Pc.items():
                    q = m + 2                       # r^q e^{bb r}, times -2
                    # d_r^i (r^q e^{bb r}) = sum_a C(i, a) q(q-1)..(q-a+1) r^{q-a} bb^{i-a} e^{bb r}
                    ff = 1
                    for a in range(i + 1):
                        if a:
                            ff *= (q - a + 1)
                        if ff == 0:
                            break
                        if bb == 0 and i - a > 0:
                            continue
                        c = K(-2) * w * v * K(comb(i, a) * ff) * (K.pw(bb, i - a) if i - a else K(1)) * K.pw(D2, q - a)
                        E.add(res, (e0 + bb * D2, '1', None), c)
    return res


def _sumd(*ds, w=None):
    out = {}
    for i, d in enumerate(ds):
        a = None if w is None else w[i]
        for k, v in d.items():
            E.add(out, k, v if a is None else a * v)
    return out


class Engine:
    """Per (case, field): G_l^(j) for a kernel family, cached; densities of AB (1s_B 1s_C)."""

    def __init__(self, case: str, K, Lmax: int, params_override=None):
        self.case, self.K = case, K
        self.D1, self.D2, self.zB, self.zC, self.zA = params_override or params(case)
        assert self.D1 < self.D2
        self.mu = 2 * self.zA
        self.it, self.kt = E.bessel_tables(K, Lmax + 8)
        self.I = E.Integrals(K)
        self.dens_cache = {}
        self.cache = {}

    def dens(self, L):
        if L not in self.dens_cache:
            self.dens_cache[L] = {inner: X.density(L, inner, self.D1, self.zC, self.it, self.kt, self.K)
                                  for inner in (True, False)}
        return self.dens_cache[L]

    def G(self, name: str, fam: dict, L: int, nmax: int) -> list:
        """[d^j/dD2^j int rho U]_{L}, j = 0..nmax, as {key: q} (bulk + jump terms)."""
        key = (name, tuple(sorted((str(k), v) for k, v in fam.items())), L, nmax)   # family in the key: the p and
        if key in self.cache:                                                    # d pairs share kernel names
            return self.cache[key]
        K = self.K
        co = {k: K(v) for k, v in fam.items()}
        jets = {near: kernel_jet(L, near, co, self.it, self.kt, K, self.mu, self.D2, nmax) for near in (True, False)}
        d = self.dens(L)
        out = []
        for j in range(nmax + 1):
            bulk = P.integrate(d, {near: jets[near][j] for near in (True, False)}, self.D1, self.D2, self.zB, self.I, K)
            if j:
                bnd = boundary_gen(d[False], jets[True], jets[False], j, self.D2, self.zB, K)
                bulk = _sumd(bulk, bnd)
            out.append({k: v for k, v in bulk.items() if v != 0})
        self.cache[key] = out
        return out

    def jump_parts(self, name: str, fam: dict, L: int, nmax: int) -> list:
        """the jump terms alone, j = 1..nmax (diagnostic)."""
        K = self.K
        co = {k: K(v) for k, v in fam.items()}
        jets = {near: kernel_jet(L, near, co, self.it, self.kt, K, self.mu, self.D2, nmax) for near in (True, False)}
        d = self.dens(L)
        return [{k: v for k, v in boundary_gen(d[False], jets[True], jets[False], j, self.D2, self.zB, K).items()
                 if v != 0} for j in range(1, nmax + 1)]


# ---------------------------------------------------------------- jets: radial operators on G_l(R)
# a jet is {(l, j, p): Fraction} meaning coef * R^{-p} * G_l^{(j)}(R)


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


def jshift(jet, a):   # a * jet / R
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


def apply_dz(field: dict) -> dict:
    """field {l: jet} (m = 0 harmonics) -> d_z field."""
    out = {}
    for l, jt in field.items():
        out[l + 1] = jsum(out.get(l + 1, {}), jscale(opA(l, jt), Fraction(l + 1, 2 * l + 1)))
        if l >= 1:
            out[l - 1] = jsum(out.get(l - 1, {}), jscale(opB(l, jt), Fraction(l, 2 * l + 1)))
    return out


def apply_dplus(field: dict, m: int) -> dict:
    """field {l: jet} of P_l^m e^{i m phi} harmonics -> d_+ field with m + 1 (terms with l' < m + 1 vanish)."""
    out = {}
    for l, jt in field.items():
        out[l + 1] = jsum(out.get(l + 1, {}), jscale(opA(l, jt), Fraction(1, 2 * l + 1)))
        if l - 1 >= m + 1:
            out[l - 1] = jsum(out.get(l - 1, {}), jscale(opB(l, jt), Fraction(-1, 2 * l + 1)))
    return out


def apply_lap(field: dict) -> dict:
    out = {}
    for l, jt in field.items():
        out[l] = jsum(out.get(l, {}), jd(jd(jt)), jshift(jd(jt), Fraction(2)),
                      jshift(jshift(jt, Fraction(1)), Fraction(-l * (l + 1))))
    return out


def py_jets(L: int) -> dict:
    """{'lap': jet, 'dzz': jet, 'beta': jet} = target-L components of Lap Phi, d_z^2 Phi, beta (d_+^2 Phi), built
    from sources l in {L-2, L, L+2}."""
    srcs = {l: {(l, 0, 0): Fraction(1)} for l in (L - 2, L, L + 2) if l >= 0}
    lap = apply_lap(srcs).get(L, {})
    dzz = apply_dz(apply_dz(srcs)).get(L, {})
    beta = apply_dplus(apply_dplus(srcs, 0), 1).get(L, {}) if L >= 2 else {}
    return {'lap': lap, 'dzz': dzz, 'beta': beta}


def jet_value(jet: dict, Gfun, D2: Fraction, K) -> dict:
    out = {}
    for (l, j, p), c in jet.items():
        a = K(c) * K.pw(D2, -p)
        for k, v in Gfun(l)[j].items():
            E.add(out, k, a * v)
    return {k: v for k, v in out.items() if v != 0}


# ---------------------------------------------------------------- the integrals


def coefficients(case: str, Lmax: int, K, integral: str, params_override=None, eng: Engine = None):
    """integral in {'PY', 'PXp', 'PA', 'DA', 'S'}.
    PY / PXp -> (a, b): lists over L of {key: q}  (series sum_L a_L P_L(x) + b_L P_L^2(x));
    PA, DA, S -> list over L of {key: q}  (S = the spherically averaged 2p pair, kernel V0 of the p density)."""
    eng = eng or Engine(case, K, Lmax + 2, params_override)
    zA = eng.zA
    if integral in ('PY', 'PXp', 'PA', 'S'):
        fam = fam_p(zA)
        if integral == 'S':
            return [eng.G('V0', fam['V0'], L, 0)[0] for L in range(Lmax + 1)]
        W = [eng.G('U0', fam['U0'], L, 0)[0] for L in range(Lmax + 1)]
        if integral == 'PA':
            return [_sumd(eng.G('U2', fam['U2'], L, 2)[2], W[L]) for L in range(Lmax + 1)]

        def Gf(l):
            return eng.G('U2', fam['U2'], l, 2)
        a, b = [], []
        half = Fraction(1, 2)
        for L in range(Lmax + 1):
            J = py_jets(L)
            mzero = jscale(jsum(J['lap'], jscale(J['dzz'], Fraction(-1))), half)
            if integral == 'PY':
                ja, jb = mzero, jscale(J['beta'], -half)
            else:
                ja = jsum(J['lap'], {(L, 2, 0): Fraction(-1)}, jscale(mzero, Fraction(-1)))
                jb = jscale(J['beta'], half)
            a.append(_sumd(jet_value(ja, Gf, eng.D2, K), W[L]))
            b.append(jet_value(jb, Gf, eng.D2, K))
            a[-1] = {k: v for k, v in a[-1].items() if v != 0}
        return a, b
    if integral == 'DA':
        fam = fam_d(zA)
        out = []
        for L in range(Lmax + 1):
            t = _sumd(eng.G('U4', fam['U4'], L, 4)[4], eng.G('U2', fam['U2'], L, 2)[2], eng.G('U0', fam['U0'], L, 0)[0])
            out.append({k: v for k, v in t.items() if v != 0})
        return out
    raise ValueError(integral)


def twopiN(case: str) -> 'mp.mpf':
    """2 pi N_B N_C (the A pair's normalization is inside the kernels)."""
    return P.twopiN(case, 'AB')


# ---------------------------------------------------------------- equivalence checks against aapc_exact


def check_against_aapc(Lcheck: int = 30) -> bool:
    """(1) family U2, U0 at nD2 = 0, 1, 2 == aapc_exact.kernel_eval with cU, cW;  (2) boundary_gen(n = 2) ==
    aapc_exact.boundary_D2;  (3) this module's PA == aapc_exact's PA, key for key, L <= Lcheck; Q and mod 2^61-1."""
    ok = True
    for case in CASES:
        for K in (E.Field(None), E.Field((1 << 61) - 1)):
            D1, D2, zB, zC, zA = params(case)
            eng = Engine(case, K, Lcheck + 2)
            za = K(zA)
            cU = {'NA': K(3) / (K(2) * za * za), 'Y': -K(3) / (K(2) * za * za), 'dY': K(33) / (K(16) * za),
                  'd2Y': K(Fraction(-9, 8)), 'd3Y': za / K(4)}
            fam = fam_p(zA)
            co = {k: K(v) for k, v in fam['U2'].items()}
            k1 = k2 = True
            for L in (0, 1, 5, Lcheck):
                for near in (True, False):
                    jet = kernel_jet(L, near, co, eng.it, eng.kt, K, eng.mu, D2, 2)
                    for j in range(3):
                        ref = P.kernel_eval(L, near, cU, eng.it, eng.kt, K, eng.mu, D2, j)
                        cl = lambda d: {k: {n: v for n, v in p.items() if v != 0} for k, p in d.items()}  # noqa
                        a = {k: v for k, v in cl(jet[j]).items() if v}
                        r = {k: v for k, v in cl(ref).items() if v}
                        k1 &= a == r
                jn = kernel_jet(L, True, co, eng.it, eng.kt, K, eng.mu, D2, 1)
                jf = kernel_jet(L, False, co, eng.it, eng.kt, K, eng.mu, D2, 1)
                dd = eng.dens(L)
                mine = {k: v for k, v in boundary_gen(dd[False], jn, jf, 2, D2, zB, K).items() if v != 0}
                ref = {k: v for k, v in P.boundary_D2(dd[False], jn[1], jf[1], D2, zB, K).items() if v != 0}
                k2 &= mine == ref
            mine = coefficients(case, Lcheck, K, 'PA', eng=eng)
            ref = P.coefficients(case, Lcheck, K, 'PA')
            k3 = all(a == {k: v for k, v in r.items() if v != 0} for a, r in zip(mine, ref))
            print(case, 'Q' if K.p is None else 'p', 'kernel jets == aapc kernel_eval:', k1,
                  ' jump n=2 == boundary_D2:', k2, ' PA == aapc PA (L <= %d):' % Lcheck, k3, flush=True)
            ok &= k1 and k2 and k3
    return ok


if __name__ == '__main__':
    print('aaa_exact vs aapc_exact equivalence:', check_against_aapc())
