"""Molecular Sturmian problem in a general-l SIGMA atom-centred basis.

Extends debug/p60_mol_sturmian_atomcentred.py (s-only, l-truncation floor
11.7 mHa) using the closed-form m=0 engine in debug/p60_sigma_onebody.py, whose
l=0 case reproduces qfd_core._inv_r with symbolic difference exactly 0.

Orbitals are (centre, Z, n, l) with Z = n*k so the decay is k for every n --
shared-scale Coulomb Sturmians.  The Sturmian identity then gives

    K[q,p] = n_p k W_own(p)[q,p] ,   K c = beta (Z_A W_A + Z_B W_B) c

VALIDATION: at the exact H2+ electronic energy beta must be 1.  The s-only basis
stalls at |beta-1| = 6.97e-3 (= 11.7 mHa); adding p_sigma, d_sigma must shrink it.
"""
import os, sys
from fractions import Fraction
import numpy as np
import sympy as sp
from scipy.linalg import eig
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from geovac.qfd_core import I2c, _gam
from geovac.two_center_eri import radial_poly, radial_norm
from debug.p60_sigma_onebody import orbital_rArB

EXACT_H2P = {2.0: -1.1026342144, 1.4: -1.2813330}


def w_element(oi, oj, which: str, R):
    """<i| 1/r_C |j>, both m=0.  oi = (centre, Z, n, l)."""
    if oi[0] == oj[0] and oi[0] == which:                     # pure one-centre
        if oi[3] != oj[3]:
            return sp.Integer(0)                              # l != l' -> 0
        ci, ai = radial_poly(Fraction(oi[1]), oi[2], oi[3])
        cj, aj = radial_poly(Fraction(oj[1]), oj[2], oj[3])
        Ni = radial_norm(Fraction(oi[1]), oi[2], oi[3])
        Nj = radial_norm(Fraction(oj[1]), oj[2], oj[3])
        a = ai + aj
        return sp.expand(Ni * Nj * sum(v1 * v2 * _gam(k1 + k2 + 1, a)
                                       for k1, v1 in ci.items() for k2, v2 in cj.items()))
    pi_, ai = orbital_rArB(oi[0], Fraction(oi[1]), oi[2], oi[3], R=R)
    pj_, aj = orbital_rArB(oj[0], Fraction(oj[1]), oj[2], oj[3], R=R)
    aA = (ai if oi[0] == "A" else 0) + (aj if oj[0] == "A" else 0)
    aB = (ai if oi[0] == "B" else 0) + (aj if oj[0] == "B" else 0)
    ang = sp.sqrt((2 * oi[3] + 1) * (2 * oj[3] + 1))
    shift = (-1, 0) if which == "A" else (0, -1)
    tot = sp.Integer(0)
    for (i1, j1), v1 in pi_.items():
        for (i2, j2), v2 in pj_.items():
            i, j = i1 + i2 + shift[0], j1 + j2 + shift[1]
            if i < -1 or j < -1:
                continue
            tot += v1 * v2 * I2c(i, j, sp.nsimplify(aA), sp.nsimplify(aB), R)
    return sp.expand(ang * tot / (4 * sp.pi))


def run(nmax: int, lmax: int, R: float, E: float):
    k = float(np.sqrt(-2.0 * E))
    kf = Fraction(k).limit_denominator(10 ** 12)
    Rs = sp.nsimplify(sp.Float(R, 30))
    basis = [(C, kf * n, n, l) for C in ("A", "B")
             for n in range(1, nmax + 1) for l in range(0, min(n, lmax + 1))]
    N = len(basis)
    WA = np.zeros((N, N)); WB = np.zeros((N, N))
    for q in range(N):
        for p in range(q, N):
            WA[q, p] = WA[p, q] = float(sp.N(w_element(basis[q], basis[p], "A", Rs), 30))
            WB[q, p] = WB[p, q] = float(sp.N(w_element(basis[q], basis[p], "B", Rs), 30))
    K = np.zeros((N, N))
    for p, (C, _z, n, _l) in enumerate(basis):
        K[:, p] = n * k * (WA if C == "A" else WB)[:, p]
    M = WA + WB
    w = eig(K, M, right=False)
    w = np.sort(np.real(w[np.abs(np.imag(w)) < 1e-8]))
    j = int(np.argmin(np.abs(w - 1.0)))
    return N, w[j], np.abs(K - K.T).max() / max(np.abs(K).max(), 1e-30), np.linalg.cond(M)


if __name__ == "__main__":
    R = 2.0; E = EXACT_H2P[R]
    print("H2+ R=2.0, exact E=%.7f.  beta must be 1;  s-only stalls at 6.97e-3 (11.7 mHa)" % E)
    for lmax in (0, 1, 2):
        for nmax in (3, 4):
            N, b, asym, cM = run(nmax, lmax, R, E)
            print("   l_max=%d n_max=%d  N=%2d  beta=%.8f  |b-1|=%.3e  asym=%.1e  cond(M)=%.1f"
                  % (lmax, nmax, N, b, abs(b - 1), asym, cM))
