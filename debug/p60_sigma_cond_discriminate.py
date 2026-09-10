"""Is the cond(M) blow-up BASIS overcompleteness, or specific to the V_0 metric?

cond(M) rose 15 -> 1.2e6 across l_max = 0..2 while |beta-1| fell 540x.  M is the
V_0 (Shibuya-Wulfman) matrix.  Two-centre high-l shared-scale Sturmians are a
classically overcomplete set, so the L2 overlap S should blow up too if the
cause is the BASIS.  If S stays tame and only M degrades, the statement is about
the METHOD's metric instead.  Same machinery, no 1/r kernel.
"""
import os, sys
from fractions import Fraction
import numpy as np
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from geovac.qfd_core import I2c, _gam
from geovac.two_center_eri import radial_poly, radial_norm
from debug.p60_sigma_onebody import orbital_rArB
from debug.p60_mol_sturmian_sigma import w_element


def s_element(oi, oj, R):
    """<i|j>, both m=0 -- w_element with no 1/r kernel."""
    if oi[0] == oj[0]:
        if oi[3] != oj[3]:
            return sp.Integer(0)
        ci, ai = radial_poly(Fraction(oi[1]), oi[2], oi[3])
        cj, aj = radial_poly(Fraction(oj[1]), oj[2], oj[3])
        Ni = radial_norm(Fraction(oi[1]), oi[2], oi[3])
        Nj = radial_norm(Fraction(oj[1]), oj[2], oj[3])
        return sp.expand(Ni * Nj * sum(v1 * v2 * _gam(k1 + k2 + 2, ai + aj)
                                       for k1, v1 in ci.items() for k2, v2 in cj.items()))
    pi_, ai = orbital_rArB(oi[0], Fraction(oi[1]), oi[2], oi[3], R=R)
    pj_, aj = orbital_rArB(oj[0], Fraction(oj[1]), oj[2], oj[3], R=R)
    aA = (ai if oi[0] == "A" else 0) + (aj if oj[0] == "A" else 0)
    aB = (ai if oi[0] == "B" else 0) + (aj if oj[0] == "B" else 0)
    ang = sp.sqrt((2 * oi[3] + 1) * (2 * oj[3] + 1))
    tot = sp.Integer(0)
    for (i1, j1), v1 in pi_.items():
        for (i2, j2), v2 in pj_.items():
            tot += v1 * v2 * I2c(i1 + i2, j1 + j2, sp.nsimplify(aA), sp.nsimplify(aB), R)
    return sp.expand(ang * tot / (4 * sp.pi))


E, R = -1.1026342144, 2.0
k = float(np.sqrt(-2.0 * E)); kf = Fraction(k).limit_denominator(10 ** 12)
Rs = sp.nsimplify(sp.Float(R, 30))
print("H2+ R=2.0.  cond(S) = L2 overlap (basis health);  cond(M) = V_0 metric")
print(" l_max n_max   N    cond(S)        cond(M)      ratio M/S")
for lmax in (0, 1, 2):
    for nmax in (3, 4):
        basis = [(C, kf * n, n, l) for C in ("A", "B")
                 for n in range(1, nmax + 1) for l in range(0, min(n, lmax + 1))]
        N = len(basis)
        S = np.zeros((N, N)); M = np.zeros((N, N))
        for q in range(N):
            for p in range(q, N):
                S[q, p] = S[p, q] = float(sp.N(s_element(basis[q], basis[p], Rs), 30))
                M[q, p] = M[p, q] = float(sp.N(w_element(basis[q], basis[p], "A", Rs)
                                               + w_element(basis[q], basis[p], "B", Rs), 30))
        cS, cM = np.linalg.cond(S), np.linalg.cond(M)
        print("   %d     %d    %2d   %11.1f   %12.1f   %8.2f" % (lmax, nmax, N, cS, cM, cM / cS))
