"""Self-checks for debug/aay_exact.py (record AY: m != 0 3d pairs on A against 1s_B 1s_C).  Small L only; nothing
here computes a scored quantity (no summed series at referee points, no decay, no constants beyond L = 8).

    python debug/aay_ap0.py family   : the hard-coded universal kernels == a fresh sympy Fourier derivation (exact);
                                       the kernel identity (d M d)^2 U4 + (d M^2 d) U2 + tr M^2 U0 == the Gauss-law
                                       multipole potential, all five M, random points; Gauss potentials vs quadrature
                                       and Poisson
    python debug/aay_ap0.py exact    : DA, PA key for key == aaa_exact (L <= 12); sum rule DA + DNX + DNY + DXY + DXX
                                       == 5 S (S: spherically averaged pair, Gauss V_0 kernel) per m, key for key,
                                       L <= 30, both cases; PY / PXp totals (this split) == record AS totals (j <= 8,
                                       L <= 60); total_coeffs vs direct evaluation of the split (random vectors)
    python debug/aay_ap0.py direct   : exact channels f{m}_L, L in (0, 1, 2, 3, 8), vs direct quadrature (Gauss-law
                                       potentials, A-frame harmonic projection, BC density tau-moments; no U kernels)
    python debug/aay_ap0.py deriv    : jump terms at r = D2 (U4 to n = 4, U2 to 2), L <= 8; 4th / 2nd D2-derivatives
                                       vs 5-point / 3-point differences at 400 digits
    python debug/aay_ap0.py smoke    : constant set at L <= 8 (mod 2^61 - 1) vs AB's 21 (already known: allowed)
Results merged into debug/data/aay_controls.json under 'AY0'.
"""
from __future__ import annotations

import json
import random
import sys
import time
from fractions import Fraction
from math import factorial
from pathlib import Path

import mpmath as mp
import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import aaa_exact as A  # noqa: E402
import aaa_symbolic as S  # noqa: E402
import aay_exact as Y  # noqa: E402
import na_exact as E  # noqa: E402

OUT = HERE / 'data' / 'aay_controls.json'
LS = (0, 1, 2, 3, 8)
SQ3 = sp.sqrt(3)
xs_, ys_, zs_ = sp.symbols('wx wy wz', real=True)
QUAD = {   # A-frame (x' = e', y' = e_y, z' = n): w M w
    'DA': (3 * zs_ ** 2 - (xs_ ** 2 + ys_ ** 2 + zs_ ** 2)) / 2,
    'DNX': SQ3 * zs_ * xs_, 'DNY': SQ3 * zs_ * ys_, 'DXY': SQ3 * xs_ * ys_,
    'DXX': SQ3 / 2 * (xs_ ** 2 - ys_ ** 2)}
CL = {0: Fraction(1, 5), 2: Fraction(2, 7), 4: Fraction(18, 35)}   # aaa's V_l = N_d^2 c_l g_l


def load() -> dict:
    return json.loads(OUT.read_text()) if OUT.exists() else {}


def save_section(key: str, val) -> None:
    """one part file per section (parallel runs); aay_controls.py merge folds them into aay_controls.json."""
    (HERE / 'data' / f'aay_part_AY0_{key}.json').write_text(json.dumps(val, indent=1))


def mf(q) -> mp.mpf:
    q = Fraction(q)
    return mp.mpf(q.numerator) / q.denominator


# ---------------------------------------------------------------- family + Gauss law


def family():
    K_, M_ = sp.symbols('K M', positive=True)
    mu = sp.sqrt(M_)
    F = 8 * sp.pi * mu / (K_ + M_) ** 2
    Nd2 = (5 * (2 * S.zA) ** 7 / (4 * sp.pi * 720)).subs(S.zA, mu / 2)
    hats = {'U4': 4 * sp.pi * Nd2 * 16 * sp.diff(F, K_, 4) / K_, 'U2': -4 * sp.pi * Nd2 * 32 * sp.diff(F, K_, 3) / K_,
            'U0': 4 * sp.pi * Nd2 * 8 * sp.diff(F, K_, 2) / K_}
    Us, ok = {}, True
    for k, h in hats.items():
        U = S._inv_ft(sp.together(h), K_, M_)
        Us[k] = sp.expand(sp.simplify(U.subs(M_, S.MU ** 2)))
        fam = S.to_family(Us[k])
        for zv in (Fraction(9, 10), Fraction(4, 5)):
            got = S.family_at(fam, zv)
            hard = {kk: v for kk, v in Y.fam_u(zv)[k].items()}
            ok &= got == hard
    print('family: hard-coded universal kernels == sympy Fourier derivation (zA = 9/10, 4/5):', ok, flush=True)
    # identity vs Gauss-law potentials, all five M, random points
    x, y, z = sp.symbols('x y z', real=True)
    R = sp.sqrt(x ** 2 + y ** 2 + z ** 2)
    V, _ = S.multipole_potentials('d')
    g = {l: V[l] / sp.Rational(CL[l].numerator, CL[l].denominator) for l in (0, 2, 4)}   # N_d^2 g_l
    worst = {}
    random.seed(11)
    pts = [{x: sp.Rational(random.randint(-30, 30), 17), y: sp.Rational(random.randint(-30, 30), 19),
            z: sp.Rational(random.randint(-30, 30), 13), S.zA: sp.Rational(random.choice((9, 4, 7)), 10)}
           for _ in range(3)]
    for name, q in QUAD.items():
        Mm = sp.hessian(q, (xs_, ys_, zs_)) / 2
        lhs = 0
        U4, U2, U0 = (Us[k].subs(S.s, R) for k in ('U4', 'U2', 'U0'))
        dvars = (x, y, z)

        def dMd(f, Mx):
            return sum(Mx[i, j] * sp.diff(f, dvars[i], dvars[j]) for i in range(3) for j in range(3) if Mx[i, j] != 0)
        lhs = dMd(dMd(U4, Mm), Mm) + dMd(U2, Mm * Mm) + (Mm * Mm).trace() * U0
        # Gauss side: harmonic parts of (w M w)^2 = H4 + s^2 H2 + s^4 H0 (H0 = 1/5 here)
        qq = sp.expand(q.subs({xs_: x, ys_: y, zs_: z}) ** 2)
        lap = lambda f: sp.expand(sum(sp.diff(f, v, 2) for v in dvars))  # noqa
        H0 = lap(lap(qq)) / 120
        H2 = sp.expand((lap(qq) - 20 * R ** 2 * H0) / 14)
        H4 = sp.expand(qq - R ** 2 * H2 - R ** 4 * H0)
        rhs = g[0].subs(S.s, R) * H0 + g[2].subs(S.s, R) * H2 / R ** 2 + g[4].subs(S.s, R) * H4 / R ** 4
        w = 0.0
        for pt in pts:
            d = sp.N((lhs - rhs).subs(pt), 40)
            mag = sp.N(rhs.subs(pt), 40)
            w = max(w, abs(float(d / mag)))
        worst[name] = w
        print(f'family: kernel identity vs Gauss-law potential, {name}: max rel {w:.2e}', flush=True)
    # Gauss potentials vs quadrature + Poisson
    t = sp.symbols('t', positive=True)
    gw = 0
    for l in (0, 2, 4):
        Vl = sp.lambdify((S.s, S.zA), g[l], 'mpmath')
        f = sp.lambdify((t, S.zA), 5 * (2 * S.zA) ** 7 / (4 * sp.pi * 720) * t ** 4 * sp.exp(-2 * S.zA * t), 'mpmath')
        with mp.workdps(40):
            for zv in (mp.mpf('0.9'), mp.mpf('0.8')):
                for s0 in (mp.mpf('0.3'), mp.mpf('1.7'), mp.mpf('4.0')):
                    inner = mp.quad(lambda u: u ** (l + 2) * f(u, zv), [0, s0])
                    outer = mp.quad(lambda u: u ** (1 - l) * f(u, zv), [s0, mp.inf])
                    gg = 4 * mp.pi / (2 * l + 1) * (s0 ** (-l - 1) * inner + s0 ** l * outer)
                    gw = max(gw, abs(gg / Vl(s0, zv) - 1))
                    d1 = mp.diff(lambda u: Vl(u, zv), s0, 1)
                    d2 = mp.diff(lambda u: Vl(u, zv), s0, 2)
                    lhs = d2 + 2 * d1 / s0 - l * (l + 1) * Vl(s0, zv) / s0 ** 2
                    gw = max(gw, abs(lhs / (-4 * mp.pi * f(s0, zv)) - 1))
    print(f'family: Gauss potentials g_l (l = 0, 2, 4) vs quadrature + Poisson, max rel {mp.nstr(gw, 3)}', flush=True)
    save_section('family', dict(hardcoded_equals_derivation=ok, identity_vs_gauss_max_rel=worst,
                                gauss_vs_quadrature_poisson_max_rel=mp.nstr(gw, 3)))


# ---------------------------------------------------------------- exact checks


def spherical_fam(zA: Fraction) -> dict:
    V, _ = S.multipole_potentials('d')
    return S.family_at(S.to_family(V[0]), zA)     # potential of N_d^2 s^4 e^{-mu s} / 5


def exact():
    res = {}
    for case in A.CASES:
        K = E.Field(None)
        Lc = 30
        eng = Y.Engine(case, K, Lc)
        t0 = time.time()
        ch = {nm: Y.coefficients(case, Lc, K, nm, eng=eng) for nm in Y.NEW + ('DA',)}
        tex = round(time.time() - t0, 1)
        Sd = [eng.G('V0_d', spherical_fam(eng.zA), L, 0)[0] for L in range(Lc + 1)]
        ok0 = ok2 = ok4 = sym = True
        for L in range(Lc + 1):
            tot = A._sumd(*[ch[nm][0][L] for nm in ch])
            tot = {k: v for k, v in tot.items() if v != 0}
            ok0 &= tot == {k: 5 * v for k, v in Sd[L].items() if v != 0}
            if L >= 2:
                s2 = {k: v for k, v in A._sumd(ch['DNX'][2][L], ch['DNY'][2][L]).items() if v != 0}
                ok2 &= not s2
                sym &= ch['DNX'][0][L] == ch['DNY'][0][L]
            if L >= 4:
                s4 = {k: v for k, v in A._sumd(ch['DXY'][4][L], ch['DXX'][4][L]).items() if v != 0}
                ok4 &= not s4
                sym &= ch['DXY'][0][L] == ch['DXX'][0][L]
        ms = {nm: sorted(c.keys()) for nm, c in ch.items()}
        print(f'{case} sum rule (L <= {Lc}, Q, key for key): f0 DA+DNX+DNY+DXY+DXX == 5 S: {ok0};  f2 DNX+DNY == 0: '
              f'{ok2};  f4 DXY+DXX == 0: {ok4};  f0(DNX)=f0(DNY), f0(DXY)=f0(DXX): {sym};  m sets {ms}  ({tex} s)',
              flush=True)
        # DA, PA vs aaa_exact
        Lk = 12
        e2 = Y.Engine(case, K, Lk)
        da = Y.coefficients(case, Lk, K, 'DA', eng=e2)[0]
        pa = Y.coefficients(case, Lk, K, 'PA', eng=e2)[0]
        rda = A.coefficients(case, Lk, K, 'DA')
        rpa = A.coefficients(case, Lk, K, 'PA')
        kda = all(a == {k: v for k, v in r.items() if v != 0} for a, r in zip(da, rda))
        kpa = all(a == {k: v for k, v in r.items() if v != 0} for a, r in zip(pa, rpa))
        print(f'{case} DA == aaa_exact DA, PA == aaa_exact PA (key for key, L <= {Lk}): {kda}, {kpa}', flush=True)
        # PY, PXp totals vs record AS's totals (both engines to L = Lp; compare c_j, j <= 8, at 1200 digits)
        ptot_all = {}
        for Lp in (60, 120):
          if True:
            e3 = Y.Engine(case, K, Lp)
            ea = A.Engine(case, K, Lp + 2)
            ptot = {}
            with mp.workdps(1200):
                tw = A.twopiN(case)
                ev = lambda q: tw * sum(mp.mpf(int(v.p)) / int(v.q) * E.value(k, mp) for k, v in q.items())  # noqa
                for nm in ('PY', 'PXp'):
                    mine = Y.coefficients(case, Lp, K, nm, eng=e3)
                    cm = Y.total_coeffs({m: [ev(q) for q in v] for m, v in mine.items()})
                    a, b = A.coefficients(case, Lp, K, nm, eng=ea)
                    aa, bb = [ev(q) for q in a], [ev(q) for q in b]
                    tail = [mp.mpf(0), mp.mpf(0)]
                    cc = [mp.mpf(0)] * (Lp + 1)
                    for j in range(Lp, -1, -1):
                        cc[j] = aa[j] - j * (j - 1) * bb[j] + 2 * (2 * j + 1) * tail[j % 2]
                        tail[j % 2] += bb[j]
                    ptot[nm] = mp.nstr(max(abs(cm[j] / cc[j] - 1) for j in range(9)), 3)
            ptot_all[Lp] = ptot
            print(f'{case} PY, PXp totals (this split, L <= {Lp}) vs record AS totals, j <= 8, max rel: {ptot}', flush=True)
        res[case] = dict(sum_rule_f0=ok0, sum_rule_f2=ok2, sum_rule_f4=ok4, pair_symmetries=sym, m_sets=ms,
                         DA_keyeq_aaa=kda, PA_keyeq_aaa=kpa, PY_PXp_total_vs_AS_rel_j_le_8=ptot_all)
    # total_coeffs vs direct evaluation, random vectors
    random.seed(3)
    with mp.workdps(80):
        f = {m: [mp.mpf(random.uniform(-1, 1)) for _ in range(41)] for m in (0, 2, 4)}
        f[2][0] = f[2][1] = 0
        for i in range(4):
            f[4][i] = 0
        c = Y.total_coeffs(f)
        worst = mp.mpf(0)
        for xv in (mp.mpf('0.37'), mp.mpf('-0.81')):
            P = [mp.legendre(L, xv) for L in range(len(c))]
            split = sum(f[m][L] * (1 - xv * xv) ** (m // 2) * mp.diff(lambda t: mp.legendre(L, t), xv, m)
                        for m in f for L in range(41))
            worst = max(worst, abs(sum(ci * Pi for ci, Pi in zip(c, P)) / split - 1))
    print('total_coeffs vs direct split evaluation (random vectors, L <= 40):', mp.nstr(worst, 3), flush=True)
    res['total_coeffs_random_check'] = mp.nstr(worst, 3)
    save_section('exact', res)


# ---------------------------------------------------------------- direct quadrature


def plm_poly(l: int, m: int) -> list:
    """coefficients (ascending) of P_l^m for even m as a polynomial in t (CS-free)."""
    t = sp.symbols('t')
    e = sp.expand((1 - t ** 2) ** (m // 2) * sp.diff(sp.legendre(l, t), t, m))
    return [Fraction(int(c.p), int(c.q)) for c in sp.Poly(e, t).all_coeffs()[::-1]] if e != 0 else [Fraction(0)]


def hcoef(q) -> dict:
    """(w M w)^2 on the unit sphere = sum h_{l m} P_l^m(cos th) cos(m ph), l, m even."""
    th, ph = sp.symbols('th ph', real=True)
    sub = {xs_: sp.sin(th) * sp.cos(ph), ys_: sp.sin(th) * sp.sin(ph), zs_: sp.cos(th)}
    f = sp.expand(q.subs(sub) ** 2)
    t = sp.symbols('t')
    out = {}
    for l in (0, 2, 4):
        for m in range(0, l + 1, 2):
            Y_ = sp.expand((1 - t ** 2) ** (m // 2) * sp.diff(sp.legendre(l, t), t, m)).subs(t, sp.cos(th)) * sp.cos(m * ph)
            num = sp.integrate(sp.integrate(sp.expand(f * Y_ * sp.sin(th)), (ph, 0, 2 * sp.pi)), (th, 0, sp.pi))
            den = sp.integrate(sp.integrate(sp.expand(Y_ * Y_ * sp.sin(th)), (ph, 0, 2 * sp.pi)), (th, 0, sp.pi))
            v = sp.nsimplify(sp.simplify(num / den))
            if v != 0:
                out[(l, m)] = Fraction(int(v.p), int(v.q))
    return out


def direct_channels(case: str, deg_out: int, deg_in: int, dps: int, H: dict) -> dict:
    """f{m}_L for every integral in H ({name: {(lam, m): h}}), L in LS, by quadrature."""
    V, _ = S.multipole_potentials('d')
    gf = {l: sp.lambdify((S.s, S.zA), V[l] / sp.Rational(CL[l].numerator, CL[l].denominator), 'mpmath')
          for l in (0, 2, 4)}
    polys = {}
    with mp.workdps(dps):
        D1, D2, zB, zC, zA = (mf(v) for v in A.params(case))
        R = D2
        NB, NC = mp.sqrt(zB ** 3 / mp.pi), mp.sqrt(zC ** 3 / mp.pi)
        GL = mp.calculus.quadrature.GaussLegendre(mp.mp)
        n_in = GL.calc_nodes(deg_in, mp.mp.prec)
        n_out = GL.calc_nodes(deg_out, mp.mp.prec)
        Lm = [(L, m) for L in LS for m in (0, 2, 4) if m <= L]
        lm = [(l, m) for l in (0, 2, 4) for m in (0, 2, 4) if m <= l]
        for (l, m) in set(Lm) | set(lm):
            polys[(l, m)] = [mf(c) for c in plm_poly(l, m)]

        def pev(key, t):
            c = polys[key]
            v = mp.mpf(0)
            for a in reversed(c):
                v = v * t + a
            return v
        acc = {(L, m, l): mp.mpf(0) for (L, m) in Lm for l in (0, 2, 4) if m <= l}
        edges = [mp.mpf(0), D1, D2, D2 + 3, D2 + 10, D2 + 25, D2 + 70]
        for lo, hi in zip(edges[:-1], edges[1:]):
            hh, cc = (hi - lo) / 2, (hi + lo) / 2
            for xo, wo in n_out:
                r = cc + hh * xo
                wr = wo * hh * r * r
                # BC density Legendre channels rho_L(r) = (2L+1)/2 int rho P_L(tau) dtau
                m0 = {L: mp.mpf(0) for L in LS}
                a_, b_ = abs(r - D1), r + D1
                h2, c2 = (b_ - a_) / 2, (b_ + a_) / 2
                pre = NB * NC * mp.exp(-zB * r)
                for xi, wi in n_in:
                    sc = c2 + h2 * xi
                    tau = (r * r + D1 * D1 - sc * sc) / (2 * r * D1)
                    w = wi * h2 * sc / (r * D1) * pre * mp.exp(-zC * sc)
                    for L in LS:
                        m0[L] += w * pev((L, 0), tau)
                rhoL = {L: (2 * L + 1) * m0[L] / 2 for L in LS}
                # A-frame kernels J_{L m l}(r) = int dt P_L^m(t) g_l(s) P_l^m(cos theta_w)
                J = {k: mp.mpf(0) for k in acc}
                a_, b_ = abs(r - R), r + R
                h2, c2 = (b_ - a_) / 2, (b_ + a_) / 2
                for xi, wi in n_in:
                    s = c2 + h2 * xi
                    t = (r * r + R * R - s * s) / (2 * r * R)
                    w = wi * h2 * s / (r * R)
                    cw = (r * t - R) / s
                    gv = {l: gf[l](s, zA) for l in (0, 2, 4)}
                    pl = {k: pev(k, cw) for k in lm}
                    pL = {k: pev(k, t) for k in Lm}
                    for (L, m, l) in J:
                        J[(L, m, l)] += w * pL[(L, m)] * gv[l] * pl[(l, m)]
                for k in acc:
                    acc[k] += wr * rhoL[k[0]] * J[k]
        out = {}
        for name, h in H.items():
            for (L, m) in Lm:
                v = sum(mf(h.get((l, m), 0)) * acc[(L, m, l)] for l in (0, 2, 4) if (L, m, l) in acc)
                out[(name, L, m)] = 2 * mp.pi * mf(Fraction(factorial(L - m), factorial(L + m))) * v
        return out


def direct():
    H = {nm: hcoef(q) for nm, q in QUAD.items()}
    print('A-frame harmonic content h_{l m} of (w M w)^2:', {k: {str(kk): str(v) for kk, v in h.items()} for k, h in H.items()},
          flush=True)
    res = {}
    dps = 60
    for case in A.CASES:
        t0 = time.time()
        K = E.Field(None)
        eng = Y.Engine(case, K, max(LS))
        ch = {nm: Y.coefficients(case, max(LS), K, nm, eng=eng) for nm in H}
        d1 = direct_channels(case, 6, 6, dps, H)
        d2 = direct_channels(case, 7, 6, dps, H)
        rec = {}
        with mp.workdps(300):
            tw = A.twopiN(case)
            ev = lambda q: tw * sum(mp.mpf(int(v.p)) / int(v.q) * E.value(k, mp) for k, v in q.items())  # noqa
            for (nm, L, m), _ in sorted(d2.items()):
                q = ch[nm].get(m, [{}] * (max(LS) + 1))[L]
                e = ev(q)
                if e == 0:
                    rec[f'{nm}_L{L}_m{m}'] = dict(exact_zero=True, direct=mp.nstr(d2[(nm, L, m)], 5))
                    continue
                r1, r2 = abs(d1[(nm, L, m)] / e - 1), abs(d2[(nm, L, m)] / e - 1)
                rec[f'{nm}_L{L}_m{m}'] = dict(rel_96=mp.nstr(r1, 3), rel_192=mp.nstr(r2, 3))
        worst = max(mp.mpf(v['rel_192']) for v in rec.values() if 'rel_192' in v)
        zeros = [k for k, v in rec.items() if 'exact_zero' in v]
        print(f'{case} direct quadrature vs exact channels (outer GL 96 / 192 per piece, inner 96, {dps} digits): '
              f'worst rel (192) {mp.nstr(worst, 3)}; exact-zero channels {zeros}; {round(time.time() - t0)} s', flush=True)
        for k, v in rec.items():
            print('   ', k, v, flush=True)
        res[case] = dict(worst_rel_192=mp.nstr(worst, 3), channels=rec)
    save_section('direct', res)


# ---------------------------------------------------------------- D2-derivatives


def deriv():
    K = E.Field(None)
    res = {}
    for case in A.CASES:
        eng = Y.Engine(case, K, 10)
        nz = []
        for kern, n in (('U4', 4), ('U2', 2)):
            fam = Y.fam_u(eng.zA)[kern]
            for L in range(9):
                for j, part in enumerate(eng.jump_parts(kern + '_d', fam, L, n), start=1):
                    if part:
                        nz.append((kern, L, j))
        D1, D2, zB, zC, zA = A.params(case)
        fd = {}
        with mp.workdps(400):
            hq = Fraction(1, 10 ** 40)
            h = mp.mpf(1) / 10 ** 40
            ev = lambda q: sum(mp.mpf(int(v.p)) / int(v.q) * E.value(k, mp) for k, v in q.items())  # noqa
            for kern, n in (('U4', 4), ('U2', 2)):
                fam = Y.fam_u(zA)[kern]
                for L in (0, 3, 8):
                    g = {k: ev(Y.Engine(case, K, 10, params_override=(D1, D2 + k * hq, zB, zC, zA))
                               .G(kern + '_d', fam, L, 0)[0]) for k in range(-2, 3)}
                    if n == 4:
                        fdv = (g[2] - 4 * g[1] + 6 * g[0] - 4 * g[-1] + g[-2]) / h ** 4
                    else:
                        fdv = (g[1] - 2 * g[0] + g[-1]) / h ** 2
                    ex = ev(eng.G(kern + '_d', fam, L, n)[n])
                    fd[f'{kern}_n{n}_L{L}'] = mp.nstr(abs(fdv / ex - 1), 3)
        print(f'{case} jump terms at r = D2 (U4 to n = 4, U2 to n = 2, L <= 8): nonzero in {nz or "none"}; '
              f'derivatives vs differences (h = 1e-40, 400 digits): {fd}', flush=True)
        res[case] = dict(jump_nonzero=[list(map(str, t)) for t in nz], fd_rel=fd)
    save_section('deriv', res)


def smoke():
    p = (1 << 61) - 1
    K = E.Field(p)
    ab = json.loads((HERE / 'data' / f'aabc_p2_{p}.json').read_text())
    res = {}
    for case in A.CASES:
        ab21 = set(ab[case]['constants'])
        eng = Y.Engine(case, K, 8)
        rec = {}
        for nm in Y.NEW:
            ch = Y.coefficients(case, 8, K, nm, eng=eng)
            for m, f in ch.items():
                L0 = 2 if m == 0 else m
                sets = [set(map(repr, f[L])) for L in range(L0, 9)]
                rec[f'{nm}_m{m}'] = dict(all_equal_AB21=all(s == ab21 for s in sets), L_range=[L0, 8])
        print(case, 'smoke (L <= 8, mod 2^61-1): constant set == AB 21:', rec, flush=True)
        res[case] = rec
    save_section('smoke_L_le_8', res)


if __name__ == '__main__':
    {'family': family, 'exact': exact, 'direct': direct, 'deriv': deriv, 'smoke': smoke}[sys.argv[1]]()
