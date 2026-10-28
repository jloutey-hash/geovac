"""WX Gaussian part: wg_common.py and wg_derive.py's model with the bond length D as a parameter (both used a fixed
D = 1.809).  Formulas unchanged; see docs/predictions/2026-10-28_water_gaussian_contact.md for their derivation.

Cases are (kind, zX, zH, D): kind 'single' takes (alpha_X, alpha_H); kind 3 or 6 is STO-3G / STO-6G on all centres."""
import mpmath as mp

import wg_common as WC

CASES = {
    'g3': ('single', '0.6', '1.5', '2.2'),
    'g4': ('single', '0.5', '1.2', '1.6'),
    'E3-3': (3, '1.5', '1.0', '2.2'), 'E3-6': (6, '1.5', '1.0', '2.2'),
    'E4-3': (3, '3.0', '1.5', '1.6'), 'E4-6': (6, '3.0', '1.5', '1.6'),
}
SETS = {'E3': ('1.5', '1.0', '2.2'), 'E4': ('3.0', '1.5', '1.6')}


def orbitals(spec):
    kind, zX, zH, _ = spec
    if kind == 'single':
        return WC.single(zX), WC.single(zH)
    return WC.contraction(kind, mp.mpf(zX)), WC.contraction(kind, mp.mpf(zH))


def quartets(spec):
    D = mp.mpf(spec[3])
    OX, OH = orbitals(spec)
    out = []
    for aa, ca in OX:
        for aj, cj in OH:
            p = aa + aj
            M1 = ca * cj * (mp.pi / p) ** 1.5 * mp.exp(-aa * aj * D ** 2 / p)
            s1 = aj / p
            for ab, cb in OX:
                for ak, ck in OH:
                    q = ab + ak
                    M2 = cb * ck * (mp.pi / q) ** 1.5 * mp.exp(-ab * ak * D ** 2 / q)
                    s2 = ak / q
                    mu = p * q / (p + q)
                    out.append(dict(M=M1 * M2, mu=mu, pref=M1 * M2 * 2 * mp.sqrt(mu / mp.pi),
                                    T0=mu * D ** 2 * (s1 ** 2 + s2 ** 2), T1=2 * mu * D ** 2 * s1 * s2,
                                    sP=s1 * D, sQ=s2 * D))
    return out


def F_case(spec, lmax):
    tot = [mp.mpf(0)] * (lmax + 1)
    for Q in quartets(spec):
        for l, v in enumerate(WC.F_quartet(Q, lmax)):
            tot[l] += v
    return tot


def boys_direct(spec, x):
    tot = mp.mpf(0)
    for Q in quartets(spec):
        T = Q['T0'] - Q['T1'] * x
        F0 = mp.sqrt(mp.pi / T) / 2 * mp.erf(mp.sqrt(T)) if T > 0 else mp.mpf(1)
        tot += Q['pref'] * F0
    return tot


def q_coeffs(T):
    """q_l = 1 + sum_j c_j / l^j, j = 2..5, on the symmetric line T0 = T1 = T (record WG)."""
    T = mp.mpf(T)
    return [-T * (T + 4) / 4, T * (T + 4), T * (2 * T ** 3 + 12 * T ** 2 - 17 * T - 196) / 16,
            -T * (9 * T ** 3 + 48 * T ** 2 + 76 * T - 272) / 8]


def single_T(spec):
    _, aX, aH, D = spec
    aX, aH, D = mp.mpf(aX), mp.mpf(aH), mp.mpf(D)
    return aH ** 2 * D ** 2 / (aX + aH)


# ---- crossover model (wg_derive.py, D as a parameter)
def V_slater(w, zH):
    N = mp.sqrt(zH ** 3 / mp.pi)
    Q = N * 8 * mp.pi / zH ** 3
    if w < mp.mpf('1e-6'):
        return Q ** 2 * (5 * zH / 16)
    J = 1 / w - mp.exp(-zH * w) * (1 / w + 11 * zH / 16 + 3 * zH ** 2 * w / 16 + zH ** 3 * w ** 2 / 48)
    return Q ** 2 * J


def leg_slater(zH, D, l):
    f = lambda x: V_slater(D * mp.sqrt(2 * (1 - x)), zH)
    return (2 * l + 1) / mp.mpf(2) * mp.quad(lambda x: f(x) * mp.legendre(l, x), mp.linspace(-1, 1, 9))


def leg_gauss(n, zH, D, l):
    tot = 0
    for ai, ci in WC.contraction(n, zH):
        for aj, cj in WC.contraction(n, zH):
            qi, qj = ci * (mp.pi / ai) ** 1.5, cj * (mp.pi / aj) ** 1.5
            m = ai * aj / (ai + aj)
            T = 2 * m * D ** 2
            Fl = WC.F_quartet(dict(pref=mp.mpf(1), T0=T, T1=T), l)[l]
            tot += qi * qj * 2 * mp.sqrt(m / mp.pi) * Fl
    return tot


def chi_slater(z, r):
    return mp.sqrt(z ** 3 / mp.pi) * mp.exp(-z * r)


def chi_gauss(n, z, r):
    return sum(c * mp.exp(-a * r ** 2) for a, c in WC.contraction(n, z))


def model(sname, n, lmax=60):
    zX, zH, D = (mp.mpf(v) for v in SETS[sname])
    lsv = {l: leg_slater(zH, D, l) for l in range(4, lmax + 1)}
    cx = (chi_gauss(n, zX, D) / chi_slater(zX, D)) ** 2
    lg = {l: leg_gauss(n, zH, D, l) for l in range(4, lmax + 1)}
    rat = {l: float(cx * lg[l] / lsv[l]) for l in range(4, lmax + 1)}
    lstar = next((l for l in range(6, lmax + 1) if rat[l] < 0.5), None)
    amax = max(a for a, _ in WC.contraction(n, zH))
    return dict(ratio=rat, lstar=lstar, alpha_max_H=float(amax), D_sqrt_amax=float(D * mp.sqrt(amax)))
