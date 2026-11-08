"""Be R12-CI analytic engine (noise-free).  Full matrix elements via block reduction.

All orbitals are s -> every term reduces to radial integrals with monopole kernels
(f0 = monopole of f ; f0sq = monopole of f^2 ; C = L=0 Coulomb 1/r> ; etc.), over the
up-block {1,2} and down-block {3,4} densities and their marginals.  Validated against MC.

Stage 2 (this file): S_11 = INT Phi0^2 F^2.  Then H_01, H_11 added + validated.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

d = np.load("debug/data/be_r12ci_ref.npz")
Z1, Z2, EREF = float(d["z1"]), float(d["z2"]), float(d["Eref"])
Zn = 4.0
N1 = 2.0 * Z1 ** 1.5
OV = N1 * 6.0 / (Z1 + Z2) ** 4
_xg, _wg = leggauss(600); _r0 = 12.5 * (_xg + 1.0); _w0 = 12.5 * _wg
_u = _r0 * np.exp(-Z2 * _r0) - OV * N1 * np.exp(-Z1 * _r0)
N2 = 1.0 / np.sqrt(np.sum(_u * _u * _r0 ** 2 * _w0))
FPI = 4 * np.pi


def a_r(r): return N1 * np.exp(-Z1 * r)
def b_r(r): return N2 * (r * np.exp(-Z2 * r) - OV * N1 * np.exp(-Z1 * r))
def f_g(r): return 1.0 - np.exp(-r)


# grids
NR = 160; Rmax = 22.0
_x, _w = leggauss(NR); rr = 0.5 * Rmax * (_x + 1.0); wr = 0.5 * Rmax * _w
w2 = wr * rr ** 2
NX = 40; xx, wx = leggauss(NX)
A = a_r(rr); B = b_r(rr)
Dp = (np.outer(A, B) - np.outer(B, A)) ** 2          # block density (NR,NR)
rho = FPI * (A ** 2 + B ** 2)                          # electron marginal
Nb = 32 * np.pi ** 2                                   # block norm UB[1]


def _r12(r1, r2, x): return np.sqrt(np.maximum(r1 * r1 + r2 * r2 - 2 * r1 * r2 * x, 1e-30))


def monopole(kernel_of_r):
    """(NR,NR) monopole  1/2 INT_-1^1 kernel(r12) dx  of a function of r12."""
    R1 = rr[:, None, None]; R2 = rr[None, :, None]; X = xx[None, None, :]
    K = kernel_of_r(_r12(R1, R2, X))
    return 0.5 * np.tensordot(K, wx, axes=([2], [0]))


F0 = monopole(f_g)                                     # monopole f
F0sq = monopole(lambda r: f_g(r) ** 2)                 # monopole f^2
Cmono = monopole(lambda r: 1.0 / r)                    # monopole 1/r12 == 1/r> (L=0 Coulomb)


def UB(Ofun):
    """8 pi^2 INT r1^2 r2^2 Dp [INT_-1^1 O dx] : up-block moment of O(r1,r2,x)."""
    R1 = rr[:, None, None]; R2 = rr[None, :, None]; X = xx[None, None, :]
    ang = np.tensordot(Ofun(R1, R2, X), wx, axes=([2], [0]))
    return 8 * np.pi ** 2 * np.einsum("i,j,ij,ij->", w2, w2, Dp, ang)


def CROSS1(P1, P3, K):
    """(4pi)^2 INT P1(r1) P3(r3) K(r1,r3) r1^2 r3^2 : single cross leg between marginals."""
    return FPI ** 2 * (P1 * w2) @ K @ (P3 * w2)


# m_f(r1) = INT |D_up(1,2)|^2 f12 d3r2 = 4pi INT r2^2 Dp(r1,r2) f0(r1,r2) dr2   (f-dressed marginal)
m_f = FPI * (Dp * F0) @ w2


def SU(P1weight, K):
    """(4pi)^3 INT r1^2 P1weight(r1) [INT r3^2 r4^2 Ddn(r3,r4)^2 K(r1,r3)K(r1,r4)] :
       one up electron (weight P1weight, e.g. rho) bridges to BOTH down electrons via K."""
    # inner(r1) = (K[r1,:]*w2) @ Dp @ (K[r1,:]*w2)
    KD = K * w2[None, :]                                # (NR,NR) rows indexed by r1
    inner = np.einsum("ij,jk,ik->i", KD, Dp, KD)
    return FPI ** 3 * np.sum(w2 * P1weight * inner)


def DISJ(K13, K24):
    """(4pi)^4 INT Dup(1,2)^2 Ddn(3,4)^2 K13(r1,r3) K24(r2,r4) : disjoint double bridge."""
    M = (Dp * w2[None, :]) @ K24 @ (w2[:, None] * Dp)   # M[i1,i3]? build carefully
    # M[i1,i3] = sum_{i2,i4} Dp[i1,i2] w2_i2 K24[i2,i4] w2_i4 Dp[i4,i3]  (Ddn symmetric)
    return FPI ** 4 * (w2) @ (M * K13) @ (w2)


def S11():
    t_intra = 2 * UB(lambda r1, r2, x: f_g(_r12(r1, r2, x)) ** 2) * Nb
    t_f12f34 = 2 * UB(lambda r1, r2, x: f_g(_r12(r1, r2, x))) ** 2
    t_self = 4 * CROSS1(rho, rho, F0sq)
    t_su = 8 * SU(rho, F0)
    t_disj = 4 * DISJ(F0, F0)
    t_f12cf = 16 * CROSS1(m_f, rho, F0)
    tot = t_intra + t_f12f34 + t_self + t_su + t_disj + t_f12cf
    return tot, dict(intra=t_intra, f12f34=t_f12f34, self=t_self, su=t_su, disj=t_disj, f12cf=t_f12cf)


if __name__ == "__main__":
    S00 = Nb ** 2
    S11v, parts = S11()
    print("S_00 =", S00)
    print("S_11 =", S11v, "  S_11/S_00 =", S11v / S00)
    for k, v in parts.items():
        print(f"   {k:8s}: {v/S00:+.5f}")
