"""Record PP (hydrogen p functions in the bond-angle series): debug/ao_exact.coefficients with a separate exponent on
each hydrogen, z1 on Y (distance D1) and z2 on Z (distance D2).  Returns q_L with F_L / (4 zX^3 (z1 z2)^(3/2)) =
sum_c q_{L,c} c, i.e. the raw integral of unnormalized exponentials up to a constant independent of every parameter.
Nothing else in ao_exact is changed; with z1 = z2 it must reproduce ao_exact.coefficients term by term.
"""
from fractions import Fraction

from na_exact import Integrals, bessel_tables
from ao_exact import KInt, Mul, add, dotJ, dotK, g_of, kmul, phi_at, phi_sum


def coefficients2(D1, D2, zX, z1, z2, Lmax, K, Lmin=0):
    D1, D2, zX, z1, z2 = (Fraction(v) for v in (D1, D2, zX, z1, z2))
    it, kt = bessel_tables(K, Lmax + 2)
    J = Integrals(K)
    Kt = KInt(K, J)
    M = Mul(K)
    lo, hi = sorted([D1, D2])
    pieces = [(Fraction(0), lo), (lo, hi), (hi, None)] if lo != hi else [(Fraction(0), lo), (lo, None)]
    out = []
    for L in range(Lmin, Lmax + 1):
        res = {}
        info = []
        for (A, B) in pieces:
            mid = A + 1 if B is None else (A + B) / 2
            g1 = g_of(L, mid < D1, D1, zX, z1, it, kt, K, 0)
            g2C = g_of(L, mid < D2, D2, zX, z2, it, kt, K, L)
            g2T = g_of(L, mid < D2, D2, zX, z2, it, kt, K, -L - 1)
            info.append((A, B, g1, g2C, g2T))
        fullC, fullT = [], []
        for (A, B, g1, g2C, g2T) in info:
            dC, dT = {}, {}
            for (e0, b), poly in g2C.items():
                dotJ(J, poly, b, A, B, dC, ((e0, ()), K(1)), K)
            for (e0, b), poly in g2T.items():
                dotJ(J, poly, b, A, B, dT, ((e0, ()), K(1)), K)
            fullC.append(dC)
            fullT.append(dT)
        for k, (A, B, g1, g2C, g2T) in enumerate(info):
            constC, constT = {}, {}
            for j in range(k):
                for kk, v in fullC[j].items():
                    add(constC, kk, v)
            for j in range(k + 1, len(info)):
                for kk, v in fullT[j].items():
                    add(constT, kk, v)
            varC, varT = [], []
            for (e0, b), poly in g2C.items():
                P, kap = phi_sum(poly, b, K)
                varC.append((e0, b, P, kap))
                for kk, v in phi_at(P, kap, b, A, K).items():
                    add(constC, kmul((e0, ()), kk), -v)
            for (e0, b), poly in g2T.items():
                P, kap = phi_sum(poly, b, K)
                varT.append((e0, b, P, kap))
                for kk, v in phi_at(P, kap, b, B, K).items():
                    add(constT, kmul((e0, ()), kk), v)
            for (e1, b1), p1 in g1.items():
                pC = {n - L - 1: v for n, v in p1.items()}
                pT = {n + L: v for n, v in p1.items()}
                for kk, v in constC.items():
                    dotJ(J, pC, b1, A, B, res, (kmul((e1, ()), kk), v), K)
                for kk, v in constT.items():
                    dotJ(J, pT, b1, A, B, res, (kmul((e1, ()), kk), v), K)
                for (e0, b, P, kap), sgn, pp in [(x, 1, pC) for x in varC] + [(x, -1, pT) for x in varT]:
                    prod = M.mul(pp, P)
                    dotJ(J, prod, b1 + b, A, B, res, ((e1 + e0, ()), K(sgn)), K)
                    if kap != 0:
                        dotK(Kt, pp, b1, -b, A, B, res, ((e1 + e0, ()), K(sgn) * kap), K)
        out.append({k: v for k, v in res.items() if v != 0})
    return out
