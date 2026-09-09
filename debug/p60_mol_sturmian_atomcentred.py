"""The one-electron molecular Sturmian problem in an ATOM-CENTRED Coulomb-
Sturmian basis -- quadrature-free, and the route Avery's method actually takes.

Why this basis rather than prolate spheroidal:  the exact closed-form engines
(geovac/qfd_core.py one-electron, geovac/two_center_eri.py two-electron) speak
atom-centred functions.  Assembling V C = V_0 B C needs the orbitals in the same
representation as the integrals, so we build them there directly instead of
projecting prolate-spheroidal orbitals across.

Shared-scale Coulomb Sturmians out of a hydrogenic engine:  a hydrogenic R_{n0}
has decay Z/n, so setting Z = n*k gives decay k for EVERY n, and orbital energy
-Z^2/(2n^2) = -k^2/2 independent of n.  That is the Sturmian property.

The construction then needs ONE primitive.  With k = p_0 = sqrt(-2E),

    (-1/2 grad^2 - E) chi_p = (-1/2 grad^2 + k^2/2) chi_p = (n_p k / r_own) chi_p

by the Sturmian identity, so

    K[q,p] = n_p k <chi_q| 1/r_own(p) |chi_p>,     V_0 = -(Z_A W_A + Z_B W_B),

and [-1/2 grad^2 + beta v_0 - E] phi = 0 projects to the generalized eigenproblem

    K c = beta (Z_A W_A + Z_B W_B) c ,      W_C[q,p] = <chi_q| 1/r_C |chi_p>.

Every W_C is a closed form from qfd_core._inv_r.  No quadrature anywhere.

VALIDATION: for H2+ at its exact electronic energy, V_0 IS V, so beta = 1 must
appear in the spectrum.
"""
import os, sys
from fractions import Fraction
import numpy as np
import sympy as sp
from scipy.linalg import eig
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from geovac import qfd_core as Q

EXACT_H2P = {2.0: -1.1026342144, 1.4: -1.2813, 2.5: -1.0021818}


def build(nmax: int, k: float, R: float, ZA=1.0, ZB=1.0):
    """Return (K, M, labels).  M = Z_A W_A + Z_B W_B."""
    kf = Fraction(k).limit_denominator(10**12)
    basis = [(C, kf * n, n) for C in ("A", "B") for n in range(1, nmax + 1)]
    N = len(basis)
    WA = np.zeros((N, N)); WB = np.zeros((N, N))
    for q in range(N):
        for p in range(N):
            WA[q, p] = float(sp.N(Q._inv_r(basis[q], basis[p], "A", R), 30))
            WB[q, p] = float(sp.N(Q._inv_r(basis[q], basis[p], "B", R), 30))
    K = np.zeros((N, N))
    for p, (C, _z, n) in enumerate(basis):
        W_own = WA if C == "A" else WB
        K[:, p] = n * k * W_own[:, p]
    M = ZA * WA + ZB * WB
    return K, M, basis, WA, WB


if __name__ == "__main__":
    R = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    E = EXACT_H2P[R]
    k = np.sqrt(-2.0 * E)
    K, M, basis, WA, WB = build(nmax, k, R)
    print("H2+  R=%.1f  E_exact=%.7f  k=p0=%.6f  basis=%d atom-centred Sturmians"
          % (R, E, k, len(basis)))
    print("bra/ket check on K (exact identity would give K = K^T): max asym = %.2e"
          % (np.abs(K - K.T).max() / max(np.abs(K).max(), 1e-30)))
    print("cond(M) = %.3f" % np.linalg.cond(M))
    w = eig(K, M, right=False)
    w = np.sort(np.real(w[np.abs(np.imag(w)) < 1e-8]))
    print("beta spectrum (real):", np.array2string(w[:6], precision=6))
    j = int(np.argmin(np.abs(w - 1.0)))
    print("beta closest to 1: %.8f   |beta-1| = %.2e   <-- must be ~0" % (w[j], abs(w[j] - 1)))
