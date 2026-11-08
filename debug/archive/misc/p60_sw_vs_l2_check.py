"""Does our sigma data support Paper 60's 'SW uniformly better-conditioned than L2'?

Paper 60 eq:sw writes S_SW = (2k^2)^-1 <grad phi|grad phi> + 1/2 <phi|phi>, which
by the Sturmian equation is (1/k^2) K -- our KINETIC-shifted matrix, NOT our
M = -<V_0>.  So cond(M) vs cond(S) is a DIFFERENT comparison from the paper's.
Compare the paper's pair before claiming any tension.
"""
import os, sys
from fractions import Fraction
import numpy as np
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from debug.p60_mol_sturmian_sigma import w_element
from debug.p60_sigma_cond_discriminate import s_element

E, R = -1.1026342144, 2.0
k = float(np.sqrt(-2.0 * E)); kf = Fraction(k).limit_denominator(10 ** 12)
Rs = sp.nsimplify(sp.Float(R, 30))
print("H2+ R=2.0.  cond(S_SW) is the PAPER's object ((1/k^2)K);  cond(S) the L2 overlap")
print(" l_max n_max   N     cond(S_L2)      cond(S_SW)    SW better?")
for lmax in (0, 1, 2):
    for nmax in (3, 4):
        basis = [(C, kf * n, n, l) for C in ("A", "B")
                 for n in range(1, nmax + 1) for l in range(0, min(n, lmax + 1))]
        N = len(basis)
        S = np.zeros((N, N)); WA = np.zeros((N, N)); WB = np.zeros((N, N))
        for q in range(N):
            for p in range(q, N):
                S[q, p] = S[p, q] = float(sp.N(s_element(basis[q], basis[p], Rs), 30))
                WA[q, p] = WA[p, q] = float(sp.N(w_element(basis[q], basis[p], "A", Rs), 30))
                WB[q, p] = WB[p, q] = float(sp.N(w_element(basis[q], basis[p], "B", Rs), 30))
        K = np.zeros((N, N))
        for p, (C, _z, n, _l) in enumerate(basis):
            K[:, p] = n * k * (WA if C == "A" else WB)[:, p]
        K = (K + K.T) / 2
        cS, cSW = np.linalg.cond(S), np.linalg.cond(K)
        print("   %d     %d    %2d   %12.1f   %12.1f    %s"
              % (lmax, nmax, N, cS, cSW, "YES" if cSW < cS else "no"))
