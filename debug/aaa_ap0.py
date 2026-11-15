"""Self-checks for debug/aaa_exact.py (PY, PXp, DA: non-s pairs on A against 1s_B 1s_C).  Small L only.

(i)  Per channel, the exact coefficients at L = 0, 1, 2, 3, 8 (both cases) against an INDEPENDENT quadrature that
     projects the A potential as written about B (Gauss's-law multipole potentials V_l of the A density, no U kernels,
     no derivative identity, no jet algebra):
       DA, PA: V(r - A) = sum_l V_l(s) P_l(u), u = (r t - R)/s, is a function of (r, R, t = r.A/rR): plain Legendre.
       PY: the lab tensor M_ij(A) = int rho_BC [Q0(s) delta_ij + q(s) w_i w_j], Q0 = V0 - V2/2, q = 3 V2 / (2 s^2),
           at A in the xz-plane; rotation about z gives PY's m = 0 part (tr M - M_zz)/2 and m = 2 part
           (M_yy - M_xx)/2; their channels follow from the addition theorem in t, with the BC density's
           tau-moments (rho, tau rho, tau^2 rho, (1 - tau^2) P_l', (1 - tau^2)^2 P_l'') projected directly.
     Gauss-Legendre in s (inner) and in r on [0, D1, D2, D2+3, D2+10, D2+25, D2+70] (outer), all integrands analytic
     on each piece; two node counts reported.
(ii) p_x'^2 + p_y^2 + p_sigma^2 = r^2: PY + PXp + PA = 3 S (S = spherically averaged 2p pair, its own kernel V0):
     exactly, per L, L <= 30 (key for key in Q), and as summed series (L <= 30) at two angles.
(iii) the V_l closed forms against Gauss's law by quadrature, the Poisson equation (numerical derivative), and the
     total charge.
Log: data/aaa_ap0.log (tee'd by the caller).
"""
from __future__ import annotations

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
import na_exact as E  # noqa: E402

LS = (0, 1, 2, 3, 8)


def mf(q) -> mp.mpf:
    q = Fraction(q)
    return mp.mpf(q.numerator) / q.denominator


def potentials():
    """mpmath callables V_l(s, zA) for the p and d densities (Gauss's formula, sympy)."""
    out = {}
    for kind in ('p', 'd'):
        V, _ = S.multipole_potentials(kind)
        out[kind] = {l: sp.lambdify((S.s, S.zA), e, 'mpmath') for l, e in V.items()}
    return out


def legendre_all(lmax: int, t):
    P = [mp.mpf(1), t]
    for n in range(1, lmax):
        P.append(((2 * n + 1) * t * P[n] - n * P[n - 1]) / (n + 1))
    return P[:lmax + 1]


# ---------------------------------------------------------------- Legendre algebra (exact, Fractions)


def antideriv(c: dict) -> dict:
    """Legendre coefficients of an antiderivative of sum c_k P_k (constant dropped): int P_k = (P_{k+1} - P_{k-1})/(2k+1)."""
    out = {}
    for k, v in c.items():
        out[k + 1] = out.get(k + 1, 0) + v / (2 * k + 1)
        if k >= 1:
            out[k - 1] = out.get(k - 1, 0) - v / (2 * k + 1)
    return out


def p2_coeffs(c: dict) -> dict:
    """(1 - x^2) f(x) = sum_L e_L P_L^2(x) for f = sum c_k P_k:  e_L = Legendre coefficients of f's double
    antiderivative (L >= 2)."""
    h = antideriv(antideriv(c))
    return {L: v for L, v in h.items() if L >= 2 and v != 0}


def deriv_legendre(l: int) -> dict:
    """P_l' = sum_{k = l-1, l-3, ...} (2k+1) P_k."""
    return {k: Fraction(2 * k + 1) for k in range(l - 1, -1, -2)}


def times_x(c: dict) -> dict:
    out = {}
    for k, v in c.items():
        out[k + 1] = out.get(k + 1, 0) + v * Fraction(k + 1, 2 * k + 1)
        if k >= 1:
            out[k - 1] = out.get(k - 1, 0) + v * Fraction(k, 2 * k + 1)
    return out


# ---------------------------------------------------------------- direct quadrature


def direct(case: str, deg_out: int, deg_in: int, dps: int, Vf: dict):
    """returns dict of per-L mp values: a (PY m=0), b (PY m=2), PA, DA; L in LS."""
    with mp.workdps(dps):
        D1, D2, zB, zC, zA = (mf(v) for v in A.params(case))
        R = D2
        NB, NC = mp.sqrt(zB ** 3 / mp.pi), mp.sqrt(zC ** 3 / mp.pi)
        GL = mp.calculus.quadrature.GaussLegendre(mp.mp)
        n_in = GL.calc_nodes(deg_in, mp.mp.prec)
        n_out = GL.calc_nodes(deg_out, mp.mp.prec)
        lmax = max(LS) + 2
        names = ('T', 'al', 'be', 'ga', 'X0', 'X1', 'X2', 'PA', 'DA')
        acc = {k: [mp.mpf(0)] * (lmax + 1) for k in names}
        Vp, Vd = Vf['p'], Vf['d']
        edges = [mp.mpf(0), D1, D2, D2 + 3, D2 + 10, D2 + 25, D2 + 70]
        for lo, hi in zip(edges[:-1], edges[1:]):
            hh, cc = (hi - lo) / 2, (hi + lo) / 2
            for xo, wo in n_out:
                r = cc + hh * xo
                wr = wo * hh * r * r
                # BC density tau-moments
                m0 = [mp.mpf(0)] * (lmax + 1)
                m1, m2, x1, x2 = list(m0), list(m0), list(m0), list(m0)
                a_, b_ = abs(r - D1), r + D1
                h2, c2 = (b_ - a_) / 2, (b_ + a_) / 2
                pre = NB * NC * mp.exp(-zB * r)
                for xi, wi in n_in:
                    sc = c2 + h2 * xi
                    tau = (r * r + D1 * D1 - sc * sc) / (2 * r * D1)
                    w = wi * h2 * sc / (r * D1) * pre * mp.exp(-zC * sc)
                    Pl = legendre_all(lmax, tau)
                    om = 1 - tau * tau
                    for l in range(lmax + 1):
                        p1 = l * ((Pl[l - 1] if l else 0) - tau * Pl[l])          # (1 - tau^2) P_l'
                        p2 = 2 * tau * p1 - l * (l + 1) * om * Pl[l]            # (1 - tau^2)^2 P_l''
                        m0[l] += w * Pl[l]
                        m1[l] += w * tau * Pl[l]
                        m2[l] += w * tau * tau * Pl[l]
                        x1[l] += w * p1
                        x2[l] += w * p2
                # A kernels in t = r.A / (r R)
                kQ0 = [mp.mpf(0)] * (lmax + 1)
                kq, kV0, kS, kD = list(kQ0), list(kQ0), list(kQ0), list(kQ0)
                a_, b_ = abs(r - R), r + R
                h2, c2 = (b_ - a_) / 2, (b_ + a_) / 2
                for xi, wi in n_in:
                    s = c2 + h2 * xi
                    t = (r * r + R * R - s * s) / (2 * r * R)
                    w = wi * h2 * s / (r * R)
                    V0, V2 = Vp[0](s, zA), Vp[2](s, zA)
                    u = (r * t - R) / s
                    P2u = (3 * u * u - 1) / 2
                    P4u = (35 * u ** 4 - 30 * u * u + 3) / 8
                    Q0, q = V0 - V2 / 2, 3 * V2 / (2 * s * s)
                    fs = V0 + V2 * P2u
                    fd = Vd[0](s, zA) + Vd[2](s, zA) * P2u + Vd[4](s, zA) * P4u
                    Pl = legendre_all(lmax, t)
                    for l in range(lmax + 1):
                        kQ0[l] += w * Q0 * Pl[l]
                        kq[l] += w * q * Pl[l]
                        kV0[l] += w * V0 * Pl[l]
                        kS[l] += w * fs * Pl[l]
                        kD[l] += w * fd * Pl[l]
                for l in range(lmax + 1):
                    f = mp.mpf(2 * l + 1) / 2
                    rl, t1, t2 = f * m0[l], f * m1[l], f * m2[l]
                    kQ, kqq, kv, ks, kd = f * kQ0[l], f * kq[l], f * kV0[l], f * kS[l], f * kD[l]
                    acc['T'][l] += wr * rl * 3 * kv
                    acc['al'][l] += wr * (rl * kQ + r * r * t2 * kqq)
                    acc['be'][l] += wr * (-2 * r * R * t1 * kqq)
                    acc['ga'][l] += wr * (R * R * rl * kqq)
                    acc['X0'][l] += wr * kqq * m0[l]
                    acc['X1'][l] += wr * kqq * r * x1[l]
                    acc['X2'][l] += wr * kqq * r * r * x2[l]
                    acc['PA'][l] += wr * rl * ks
                    acc['DA'][l] += wr * rl * kd
        fp = [4 * mp.pi / (2 * l + 1) for l in range(lmax + 1)]
        for k in ('T', 'al', 'be', 'ga', 'PA', 'DA'):
            acc[k] = [v * fp[l] for l, v in enumerate(acc[k])]
        # M_zz Legendre coefficients
        Mzz = {}
        for l in range(lmax + 1):
            for part, cvec in (('al', {l: Fraction(1)}), ('be', times_x({l: Fraction(1)})),
                               ('ga', times_x(times_x({l: Fraction(1)})))):
                for L, v in cvec.items():
                    Mzz[L] = Mzz.get(L, 0) + mf(v) * acc[part][l]
        a = {L: (acc['T'][L] - Mzz[L]) / 2 for L in LS}
        # (M_yy - M_xx)/2 in P_L^2
        bf = {}
        for l in range(lmax + 1):
            if l >= 2:
                c = -Fraction(factorial(l - 2), factorial(l + 2))
                bf[l] = bf.get(l, 0) + mf(c) * acc['X2'][l]
            if l >= 1:
                c = 2 * Fraction(factorial(l - 1), factorial(l + 1))
                for L, v in p2_coeffs(deriv_legendre(l)).items():
                    bf[L] = bf.get(L, 0) + mf(c * v) * R * acc['X1'][l]
            for L, v in p2_coeffs({l: Fraction(1)}).items():
                bf[L] = bf.get(L, 0) - mf(v) * R * R * acc['X0'][l]
        b = {L: mp.pi * bf.get(L, 0) for L in LS}
        return dict(a=a, b=b, PA={L: acc['PA'][L] for L in LS}, DA={L: acc['DA'][L] for L in LS})


def exact_values(case: str, dps: int):
    K = E.Field(None)
    eng = A.Engine(case, K, max(LS) + 2)
    a, b = A.coefficients(case, max(LS), K, 'PY', eng=eng)
    pa = A.coefficients(case, max(LS), K, 'PA', eng=eng)
    da = A.coefficients(case, max(LS), K, 'DA', eng=eng)
    with mp.workdps(dps):
        tw = A.twopiN(case)

        def val(q):
            return tw * sum(mp.mpf(int(v.p)) / int(v.q) * E.value(k, mp) for k, v in q.items())
        return dict(a={L: val(a[L]) for L in LS}, b={L: val(b[L]) for L in LS}, PA={L: val(pa[L]) for L in LS},
                    DA={L: val(da[L]) for L in LS}), eng


def check_i(Vf) -> bool:
    ok = True
    dps = 60
    for case in A.CASES:
        t0 = time.time()
        ex, _ = exact_values(case, 300)
        d1 = direct(case, 6, 6, dps, Vf)
        d2 = direct(case, 7, 6, dps, Vf)
        print(f'{case} (i) exact vs direct quadrature (outer GL 96 / 192 nodes per piece, inner 96), '
              f'{round(time.time() - t0)} s', flush=True)
        with mp.workdps(dps):
            for ch in ('a', 'b', 'PA', 'DA'):
                for L in LS:
                    e = ex[ch][L]
                    if ch == 'b' and L < 2:
                        print(f'   {ch:2s} L={L}: exact {mp.nstr(e, 5)} (P_L^2 vanishes for L < 2)', flush=True)
                        ok &= e == 0
                        continue
                    r1, r2 = abs(d1[ch][L] / e - 1), abs(d2[ch][L] / e - 1)
                    good = r2 < mp.mpf('1e-18')
                    ok &= bool(good)
                    print(f'   {ch:2s} L={L}: exact {mp.nstr(e, 26)}  rel diff {mp.nstr(r1, 3)} / {mp.nstr(r2, 3)}'
                          f'  {"ok" if good else "FAIL"}', flush=True)
    return ok


def check_ii() -> bool:
    ok = True
    Lc = 30
    for case in A.CASES:
        K = E.Field(None)
        eng = A.Engine(case, K, Lc + 2)
        ay, by = A.coefficients(case, Lc, K, 'PY', eng=eng)
        ax, bx = A.coefficients(case, Lc, K, 'PXp', eng=eng)
        pa = A.coefficients(case, Lc, K, 'PA', eng=eng)
        sS = A.coefficients(case, Lc, K, 'S', eng=eng)
        exact = True
        for L in range(Lc + 1):
            lhs = A._sumd(ay[L], ax[L], pa[L])
            lhs = {k: v for k, v in lhs.items() if v != 0}
            rhs = {k: 3 * v for k, v in sS[L].items() if v != 0}
            exact &= lhs == rhs
            bsum = {k: v for k, v in A._sumd(by[L], bx[L]).items() if v != 0}
            exact &= not bsum
        print(f'{case} (ii) PY + PXp + PA == 3 S per L, exactly (Q, key for key, L <= {Lc}), b parts cancel: {exact}',
              flush=True)
        ok &= exact
        with mp.workdps(600):
            tw = A.twopiN(case)
            ev = lambda q: tw * sum(mp.mpf(int(v.p)) / int(v.q) * E.value(k, mp) for k, v in q.items())  # noqa
            for xr in (mp.mpf('0.3'), mp.mpf('-0.6')):
                Pl = legendre_all(Lc + 1, xr)
                P2 = []
                for L in range(Lc + 1):     # P_L^2 = (1 - x^2) P_L'' via the Legendre ODE
                    p1 = L * ((Pl[L - 1] if L else 0) - xr * Pl[L])
                    P2.append((2 * xr * p1 - L * (L + 1) * (1 - xr * xr) * Pl[L]) / (1 - xr * xr))
                tot = sum(ev(ay[L]) * Pl[L] + ev(by[L]) * P2[L] + ev(ax[L]) * Pl[L] + ev(bx[L]) * P2[L]
                          + ev(pa[L]) * Pl[L] for L in range(Lc + 1))
                ref = 3 * sum(ev(sS[L]) * Pl[L] for L in range(Lc + 1))
                rel = abs(tot / ref - 1)
                print(f'   summed (L <= {Lc}) at x = {mp.nstr(xr, 3)} (600 digits): |(PY + PXp + PA)/(3 S) - 1| = '
                      f'{mp.nstr(rel, 3)}', flush=True)
                ok &= rel < mp.mpf('1e-40')
    return ok


def check_iii(Vf) -> bool:
    ok = True
    t = sp.symbols('t', positive=True)
    for kind in ('p', 'd'):
        V, fl = S.multipole_potentials(kind)
        for zval in (mp.mpf('0.9'), mp.mpf('0.8')):
            with mp.workdps(40):
                for l, Vl in Vf[kind].items():
                    f = sp.lambdify((t, S.zA), fl[l], 'mpmath')
                    worst = mp.mpf(0)
                    for s0 in (mp.mpf('0.3'), mp.mpf('1.7'), mp.mpf('4.0')):
                        inner = mp.quad(lambda y: y ** (l + 2) * f(y, zval), [0, s0])
                        outer = mp.quad(lambda y: y ** (1 - l) * f(y, zval), [s0, mp.inf])
                        g = 4 * mp.pi / (2 * l + 1) * (s0 ** (-l - 1) * inner + s0 ** l * outer)
                        worst = max(worst, abs(g / Vl(s0, zval) - 1))
                        # Poisson: V'' + 2 V'/s - l(l+1) V / s^2 = -4 pi f_l
                        d1 = mp.diff(lambda y: Vl(y, zval), s0, 1)
                        d2 = mp.diff(lambda y: Vl(y, zval), s0, 2)
                        lhs = d2 + 2 * d1 / s0 - l * (l + 1) * Vl(s0, zval) / s0 ** 2
                        worst = max(worst, abs(lhs / (-4 * mp.pi * f(s0, zval)) - 1))
                    good = worst < mp.mpf('1e-25')
                    ok &= bool(good)
                    print(f'(iii) {kind} V_{l} zA={mp.nstr(zval, 2)}: max rel dev from Gauss law + Poisson at '
                          f's = 0.3, 1.7, 4.0: {mp.nstr(worst, 3)} {"ok" if good else "FAIL"}', flush=True)
                charge = Vf[kind][0](mp.mpf(200), zval) * 200
                ok &= abs(charge - 1) < mp.mpf('1e-30')
                print(f'(iii) {kind} total charge s V_0(s) at s = 200: {mp.nstr(charge, 30)}', flush=True)
    # the kernel decomposition sum_j (e.grad)^{2j} U_2j == sum_l V_l P_l (exact sympy, random points)
    for kind in ('p', 'd'):
        r = S.check_identity(kind)
        print(f'(iii) {kind}: sum_j (e.grad)^(2j) U_2j - sum_l V_l P_l(u), relative, 4 random points: {r}', flush=True)
        ok &= r == 0
    return ok


def check_jumps() -> None:
    """Reported: the jump terms at r = D2 of every D2-derivative used (U2_p to n = 2; U2_d to 2, U4_d to 4) at
    L <= 8, exactly in Q; and the 4th D2-derivative of I_U4 against a central difference at D2 +- k 10^-40."""
    K = E.Field(None)
    for case in A.CASES:
        eng = A.Engine(case, K, 10)
        zA = eng.zA
        fp, fd = A.fam_p(zA), A.fam_d(zA)
        nz = []
        for (nm, fam, n) in (('U2_p', fp['U2'], 2), ('U2_d', fd['U2'], 2), ('U4_d', fd['U4'], 4)):
            for L in range(9):
                for j, part in enumerate(eng.jump_parts(nm, fam, L, n), start=1):
                    if part:
                        nz.append((nm, L, j))
        print(f'{case} (iv) jump terms at r = D2, all n used, L <= 8: nonzero in {nz if nz else "none (all identically 0)"}',
              flush=True)
        D1, D2, zB, zC, _ = A.params(case)
        with mp.workdps(400):
            h = Fraction(1, 10 ** 40)
            ev = lambda q: sum(mp.mpf(int(v.p)) / int(v.q) * E.value(k, mp) for k, v in q.items())  # noqa
            for L in (0, 3, 8):
                g = {k: ev(A.Engine(case, K, 10, params_override=(D1, D2 + k * h, zB, zC, zA)).G('U4_d', fd['U4'], L, 0)[0])
                     for k in range(-2, 3)}
                fdv = (g[2] - 4 * g[1] + 6 * g[0] - 4 * g[-1] + g[-2]) / (mp.mpf(1) / 10 ** 40) ** 4
                ex = ev(eng.G('U4_d', fd['U4'], L, 4)[4])
                print(f'   L={L}: d^4/dD2^4 int rho U4 exact vs 5-point difference: rel {mp.nstr(abs(fdv / ex - 1), 3)}',
                      flush=True)


def main():
    if sys.argv[1:] == ['jumps']:
        check_jumps()
        return
    t0 = time.time()
    Vf = potentials()
    o3 = check_iii(Vf)
    o2 = check_ii()
    o1 = check_i(Vf)
    print('AAA0 self-checks: (i)', o1, ' (ii)', o2, ' (iii)', o3, ' all:', o1 and o2 and o3,
          round(time.time() - t0), 's', flush=True)


if __name__ == '__main__':
    main()
