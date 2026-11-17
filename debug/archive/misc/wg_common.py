"""WG phase 2 shared code (record docs/predictions/2026-10-28_water_gaussian_contact.md).

Quartets of s-Gaussian primitives for (XY|XZ), symmetric geometry, and their Legendre coefficients by the series form
of the record:  F_l = pref (2l+1) sum_k T1^(l+2k) / (2^k k! (2l+2k+1)!!) e^(-T0) M(1, l+2k+3/2, T0) / (2l+4k+1).
M(1, b, T0) on b = m + 3/2 by the downward recurrence M_b = 1 + (T0/b) M_(b+1), seeded by mpmath hyp1f1.
"""
import mpmath as mp

D = mp.mpf('1.809')
Z2H = mp.mpf('1.24') ** 2
STO = {
    3: ([mp.mpf('3.42525091'), mp.mpf('0.62391373'), mp.mpf('0.16885540')],
        [mp.mpf('0.15432897'), mp.mpf('0.53532814'), mp.mpf('0.44463454')]),
    6: ([mp.mpf('0.3552322122E+02'), mp.mpf('0.6513143725E+01'), mp.mpf('0.1822142904E+01'),
         mp.mpf('0.6259552659E+00'), mp.mpf('0.2430767471E+00'), mp.mpf('0.1001124280E+00')],
        [mp.mpf('0.9163596281E-02'), mp.mpf('0.4936149294E-01'), mp.mpf('0.1685383049E+00'),
         mp.mpf('0.3705627997E+00'), mp.mpf('0.4164915298E+00'), mp.mpf('0.1303340841E+00')]),
}


def contraction(n, zeta):
    al, d = STO[n]
    al = [a / Z2H for a in al]
    terms = [(zeta ** 2 * a, di * (2 * zeta ** 2 * a / mp.pi) ** mp.mpf(0.75)) for a, di in zip(al, d)]
    S = sum(ci * cj * (mp.pi / (ai + aj)) ** 1.5 for ai, ci in terms for aj, cj in terms)
    return [(a, c / mp.sqrt(S)) for a, c in terms]


def single(alpha):
    return [(mp.mpf(alpha), (2 * mp.mpf(alpha) / mp.pi) ** mp.mpf(0.75))]


CASES = {'g1': ('single', 1, 1), 'g2': ('single', '0.3', 2),
         'E1-3': (3, 1, 1), 'E1-6': (6, 1, 1), 'E2-3': (3, '2.25', '1.24'), 'E2-6': (6, '2.25', '1.24')}


def orbitals(case):
    kind, zX, zH = CASES[case]
    if kind == 'single':
        return single(zX), single(zH)
    return contraction(kind, mp.mpf(zX)), contraction(kind, mp.mpf(zH))


def quartets(case):
    """List of dicts: pref, T0, T1, M1M2, mu, sP, sQ (|P|, |Q|)."""
    OX, OH = orbitals(case)
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


def kummer_table(T0, mtop):
    """M(1, m + 3/2, T0) for m = 0..mtop, downward recurrence from a hyp1f1 seed."""
    M = [mp.mpf(0)] * (mtop + 1)
    M[mtop] = mp.hyp1f1(1, mtop + mp.mpf(3) / 2, T0)
    for m in range(mtop - 1, -1, -1):
        M[m] = 1 + T0 / (m + mp.mpf(3) / 2) * M[m + 1]
    return M


def F_quartet(Q, lmax, tol_exp=None):
    """Legendre coefficients l = 0..lmax of one quartet, series form."""
    tol = mp.mpf(10) ** (-(tol_exp or (mp.mp.dps - 5)))
    T0, T1 = Q['T0'], Q['T1']
    K = int(T1) + 80
    mtop = lmax + 2 * K + 4
    M = kummer_table(T0, mtop)
    eT = mp.exp(-T0)
    out = []
    t0 = eT * M[0]                     # l = 0, k = 0 term: T1^0/1!! e^-T0 M_0 / 1
    for l in range(lmax + 1):
        if l > 0:
            t0 = t0 * T1 / (2 * l + 1) * M[l] / M[l - 1] * (2 * l - 1) / (2 * l + 1)
        tot, t, k = t0, t0, 0
        while True:
            m = l + 2 * k
            t = (t * T1 ** 2 / (2 * (k + 1) * (2 * l + 2 * k + 3)) * M[m + 2] / M[m]
                 * (2 * l + 4 * k + 1) / (2 * l + 4 * k + 5))
            k += 1
            tot += t
            if k > T1 and t < tot * tol:
                break
            if l + 2 * k + 2 > mtop:
                raise RuntimeError("Kummer table too short")
        out.append(Q['pref'] * (2 * l + 1) * tot)
    return out


def F_case(case, lmax):
    tot = [mp.mpf(0)] * (lmax + 1)
    for Q in quartets(case):
        for l, v in enumerate(F_quartet(Q, lmax)):
            tot[l] += v
    return tot


def F_weber_quartet(Q, l):
    """Independent route: M1 M2 (2(2l+1)/pi) int_0^inf e^(-k^2/(4 mu)) j_l(k |P|) j_l(k |Q|) dk."""
    a = 1 / (4 * Q['mu'])
    sP, sQ = Q['sP'], Q['sQ']
    kmax = mp.sqrt((mp.mp.dps * mp.log(10) + 20) / a)
    npc = int(kmax * max(sP, sQ) / mp.pi) + 8
    pts = mp.linspace(0, kmax, npc)

    def f(k):
        if k == 0:
            return mp.mpf(1) if l == 0 else mp.mpf(0)
        return (mp.exp(-a * k ** 2) * mp.pi / (2 * k * mp.sqrt(sP * sQ))
                * mp.besselj(l + 0.5, k * sP) * mp.besselj(l + 0.5, k * sQ))
    return Q['M'] * 2 * (2 * l + 1) / mp.pi * mp.quad(f, pts)


def boys_direct(case, x):
    tot = mp.mpf(0)
    for Q in quartets(case):
        T = Q['T0'] - Q['T1'] * x
        F0 = mp.sqrt(mp.pi / T) / 2 * mp.erf(mp.sqrt(T)) if T > 0 else mp.mpf(1)
        tot += Q['pref'] * F0
    return tot
