"""Record L52 controls (measured before registering).
    python debug/l52_controls.py sym     : odd coefficients of V_A (Coulomb, Yukawa line) -> 3/2, 5/2 coefficients
    python debug/l52_controls.py const   : C_pq * Beta  ==  corpus K(p, q) of debug/aaeq2_local.py (same pinch constant)
    python debug/l52_controls.py equal   : equal-arm 1s 7/2 amplitude vs record AE2b's closed form (aaeq2_amp.log, E3),
                                           and T_L vs its pred_F
    python debug/l52_controls.py na [dps]: NA's 3/2 amplitude at N1, N2 from na_exact's exact F_L, L <= 240
    python debug/l52_controls.py toy     : instrument control: ratio extrapolation on exact toy coefficients
    python debug/l52_controls.py toy92   : the same for 9/2 toys with large nuisances (A1, A3 x*)
Output: debug/data/l52_controls_<mode>.json
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import l52_local as LL

DATA = HERE / 'data'


def sym():
    s, z, k = sp.symbols('s zeta kappa', positive=True)
    mu = 2 * z
    VC = 1 / s - (z + 1 / s) * sp.exp(-mu * s)
    al = 1 / (mu ** 2 - k ** 2) ** 2
    VY = mu ** 4 * al * (sp.exp(-k * s) - sp.exp(-mu * s)) / s - mu ** 3 / (2 * (mu ** 2 - k ** 2)) * sp.exp(-mu * s)
    out = {}
    for name, V in (('coulomb', VC), ('yukawa', VY)):
        ser = sp.series(V, s, 0, 8).removeO()
        out[name] = {str(n): str(sp.factor(sp.simplify(ser.coeff(s, n)))) for n in (-1, 1, 3, 5, 7)}
    # Yukawa -> Coulomb at kappa -> 0, and the kappa-independence of v3
    out['yukawa_kappa0_minus_coulomb'] = str(sp.simplify(sp.limit(VY, k, 0) - VC))
    # Generic bounded spherical density rho = sum c_n s^n: Poisson (Helmholtz) gives v_{-1} = v_1 = 0 identically.
    # V_kappa solves (-lap + kappa^2) V = 4 pi rho_A, rho_A = (zeta^3/pi) e^{-mu s} (radial Laplacian), for s > 0
    lap = sp.diff(s ** 2 * sp.diff(VY, s), s) / s ** 2
    out['yukawa_helmholtz_residual'] = str(sp.simplify(-lap + k ** 2 * VY - 4 * sp.pi * z ** 3 / sp.pi * sp.exp(-mu * s)))
    out['note'] = ('(-lap + kappa^2) V = 4 pi rho, rho bounded: s^-1 needs a point charge, s^1 needs rho ~ 1/s; '
                   'v3 = -4 pi rho_1 / 12 for every kappa')
    return out


def const():
    sys.path.insert(0, str(HERE))
    import aaeq2_local as A
    mp.mp.dps = 30
    res = {}
    for (p, q) in ((-1, 1), (3, 1), (5, 1), (3, 3), (-1, 3)):
        mine = LL.C_pq(p, q) * mp.beta(mp.mpf(p + 3) / 2, mp.mpf(q + 3) / 2)
        res[f'{p},{q}'] = dict(mine=mp.nstr(mine, 20), corpus_K=mp.nstr(A.K(q, p), 20))
    return res


def equal():
    """E3: D = 3, zB = 3, zC = 1, zA = 9/10 (record AE2)."""
    mp.mp.dps = 40
    D, zB, zC, zA = mp.mpf(3), mp.mpf(3), mp.mpf(1), mp.mpf(9) / 10
    v3 = (2 * zA) ** 4 / 24
    amp = LL.amplitude(3, 1, v3, D, D, zB, zC)
    ae = mp.mpf('1.3356449576164494005e-6')
    predF = {220: '5.7665712651334e-18', 400: '4.86608260103526e-20', 580: '2.4976942938652e-21'}
    # AE2b's amplitude multiplies |a|^7 = (2D^2)^{7/2} (1-x)^{7/2}; ours multiplies (x*-x)^{7/2}
    mine = amp / (2 * D * D) ** mp.mpf(3.5)
    tchk = {str(L): mp.nstr(amp * LL.T(3.5, mp.mpf(1), L) / mp.mpf(v) - 1, 5)
            for L, v in predF.items()}
    return dict(mine=mp.nstr(mine, 20), AE2b=mp.nstr(ae, 20), rel=mp.nstr(mine / ae - 1, 5), T_vs_predF_rel=tchk)


def na(dps=2000, Lmax=240):
    import na_exact as E
    out = {}
    for case, (D1, D2, zB, zC) in {'N1': ('9/5', '13/5', '13/10', '1'), 'N2': ('3/2', '12/5', '1', '6/5')}.items():
        t0 = time.time()
        Q = E.coefficients(D1, D2, zB, zC, Lmax, E.Field(None))
        mp.mp.dps = dps
        q = [Fraction(v) for v in (D1, D2, zB, zC)]
        d1, d2, zb, zc = (mp.mpf(v.numerator) / v.denominator for v in q)
        tw = 2 * mp.sqrt(zb ** 3 * zc ** 3)
        F = [tw * sum(mp.mpf(int(v.p)) / int(v.q) * E.value(k, mp) for k, v in qq.items()) for qq in Q]
        mp.mp.dps = 60
        xs = LL.xstar(d1, d2)
        amp = LL.amplitude(-1, 1, 1, d1, d2, zb, zc)
        Ls = list(range(120, Lmax + 1, 20))
        R = [mp.mpf(F[L]) / (amp * LL.T(1.5, xs, L)) for L in Ls]
        r0, c1 = LL.ratio_extrap(Ls, R)
        out[case] = dict(amp=mp.nstr(amp, 15), R={str(L): mp.nstr(r, 10) for L, r in zip(Ls, R)},
                         R_extrap=mp.nstr(r0, 10), c1=mp.nstr(c1, 6), secs=round(time.time() - t0))
        print(case, out[case], flush=True)
    return out


def toy92():
    """9/2 toys with nuisances (5, -10) and (-8, 20) in (x*-x)^{11/2}, ^{13/2}, at A1's and A3's x*."""
    mp.mp.dps = 60
    out = {}
    for nm, xs in (('A1', LL.xstar(mp.mpf(9) / 5, mp.mpf(13) / 5)), ('A3', LL.xstar(mp.mpf(2), mp.mpf(3)))):
        Ls = list(range(120, 241, 20))
        for c, d in ((5, -10), (-8, 20)):
            R = [(LL.T(4.5, xs, L) + c * LL.T(5.5, xs, L) + d * LL.T(6.5, xs, L)) / LL.T(4.5, xs, L) for L in Ls]
            r0, c1 = LL.ratio_extrap(Ls, R)
            out[f'{nm} c={c} d={d}'] = dict(R240=float(R[-1]), R_extrap=float(r0))
    return out


def toy():
    """f = (x*-x)^{7/2} + 3 (x*-x)^{9/2} - 5 (x*-x)^{11/2}: exact T_L; does the extrapolation recover 1 at L=120..240?"""
    mp.mp.dps = 60
    xs = LL.xstar(mp.mpf(9) / 5, mp.mpf(13) / 5)
    Ls = list(range(120, 241, 20))
    R = [(LL.T(3.5, xs, L) + 3 * LL.T(4.5, xs, L) - 5 * LL.T(5.5, xs, L)) / LL.T(3.5, xs, L) for L in Ls]
    r0, c1 = LL.ratio_extrap(Ls, R)
    return dict(R={str(L): mp.nstr(r, 10) for L, r in zip(Ls, R)}, R_extrap=mp.nstr(r0, 10), c1=mp.nstr(c1, 6),
                R240_minus_1=mp.nstr(R[-1] - 1, 5))


if __name__ == '__main__':
    mode = sys.argv[1]
    t = time.time()
    res = {'sym': sym, 'const': const, 'equal': equal, 'toy': toy, 'toy92': toy92}[mode]() if mode != 'na' else \
        na(int(sys.argv[2]) if len(sys.argv) > 2 else 2000)
    res['secs'] = round(time.time() - t)
    (DATA / f'l52_controls_{mode}.json').write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
