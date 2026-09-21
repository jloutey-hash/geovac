"""Be R12-CI, analytic engine (trustworthy, noise-free).  Stage 1: block foundation.

Phi_0 = D_up(1,2) D_dn(3,4), D_up = a(r1)b(r2)-b(r1)a(r2), a=1s,b=2s (radial R, INT R^2 r^2 dr=1).
The measure factorizes into up-block {1,2} and down-block {3,4}.  Any multiplicative
expectation reduces to block moments.  Angular reduction: for an s-orbital pair, an operator
depending on (r1,r2,cos th12) integrates as
    INT |D_up|^2 O d3r1 d3r2 = 8 pi^2 INT r1^2 r2^2 D_up(r1,r2)^2 [INT_-1^1 O(r1,r2,x) dx] dr1 dr2.

Stage 1 validates: S_00, the gate H_00/S_00 = E0, S_01 (vs the monopole Slater-Condon form).
Later stages add the singular V_ne/V_ee pieces (block-reduced) and kinetic (smooth MC).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

d = np.load("debug/data/be_r12ci_ref.npz")
Z1, Z2, EREF = float(d["z1"]), float(d["z2"]), float(d["Eref"])
Zn = 4.0
N1 = 2.0 * Z1 ** 1.5
OV = N1 * 6.0 / (Z1 + Z2) ** 4
_xg, _wg = leggauss(600); _rr0 = 12.5 * (_xg + 1.0); _wr0 = 12.5 * _wg
_u = _rr0 * np.exp(-Z2 * _rr0) - OV * N1 * np.exp(-Z1 * _rr0)
N2 = 1.0 / np.sqrt(np.sum(_u * _u * _rr0 ** 2 * _wr0))


def a_rad(r):
    return N1 * np.exp(-Z1 * r)


def b_rad(r):
    return N2 * (r * np.exp(-Z2 * r) - OV * N1 * np.exp(-Z1 * r))


def f_gem(r):
    return 1.0 - np.exp(-r)


# ------- radial + angular grids for block quadrature ------------------------ #
NR = 200
Rmax = 22.0
_x, _w = leggauss(NR)
rr = 0.5 * Rmax * (_x + 1.0)
wr = 0.5 * Rmax * _w
NX = 48
xx, wx = leggauss(NX)              # cos(theta12) in [-1,1]

# precompute orbitals on the grid
A = a_rad(rr); B = b_rad(rr)
# radial weight r^2 dr
w2 = wr * rr ** 2

# D_up(r1,r2)^2 on the (r1,r2) grid
Dup2 = (np.outer(A, B) - np.outer(B, A)) ** 2      # (NR,NR)


def r12_of(r1, r2, x):
    return np.sqrt(np.maximum(r1 * r1 + r2 * r2 - 2 * r1 * r2 * x, 1e-30))


def UB(Ofun):
    """up-block moment  8 pi^2 INT r1^2 r2^2 D_up^2 [INT_-1^1 O dx] dr1 dr2.
       Ofun(r1grid[:,None,None], r2grid[None,:,None], x[None,None,:]) -> (NR,NR,NX)."""
    R1 = rr[:, None, None]; R2 = rr[None, :, None]; X = xx[None, None, :]
    O = Ofun(R1, R2, X)                              # (NR,NR,NX)
    ang = np.tensordot(O, wx, axes=([2], [0]))       # INT_-1^1 O dx  -> (NR,NR)
    return 8 * np.pi ** 2 * np.einsum("i,j,ij,ij->", w2, w2, Dup2, ang)


# one-electron-marginal dressings for cross (up-down) terms:
#   rho_up(r1) = INT |D_up(1,2)|^2 d3r2  (electron-1 marginal of the up block)
def rho_up_marginal():
    # INT d3r2 D_up^2 = 4pi INT r2^2 D_up(r1,r2)^2 dr2 ; D_up^2 angle-indep in r2
    m = 4 * np.pi * (Dup2 * w2[None, :]).sum(axis=1)  # (NR,)  function of r1
    return m


if __name__ == "__main__":
    print("=" * 66)
    print("Stage 1: block foundation + gate")
    print("=" * 66)
    S00_up = UB(lambda r1, r2, x: np.ones_like(r1 + r2 + x))
    S00 = S00_up ** 2                                 # up-block x down-block (identical)
    print(f"UB[1] = {S00_up:.4f}   (expect 32 pi^2 = {32*np.pi**2:.4f})")
    print(f"S_00 = UB[1]^2 = {S00:.4f}")

    # monopole f-integral cross check: S_01 via block reduction
    # F = f12 + f34 + (f13+f14+f23+f24)
    # <f12> block = UB[f12]*UB[1] ; <cross f13> = INT rho_up(r1) rho_dn(r3) f0(r1,r3)
    f12_up = UB(lambda r1, r2, x: f_gem(r12_of(r1, r2, x)))     # UB[f12]
    S01_intra = 2 * f12_up * S00_up                            # f12 (up) + f34 (dn), each *other block norm
    # cross: rho marginals + monopole f
    rho = rho_up_marginal()                                    # (NR,)
    def f_monopole(r1grid, r3grid):
        a = r1grid[:, None]; b = r3grid[None, :]
        out = np.zeros((a.shape[0], b.shape[1]))
        for x, w in zip(xx, wx):
            out += 0.5 * w * f_gem(r12_of(a, b, x))
        return out
    f0 = f_monopole(rr, rr)
    # <f_13> = INT rho_up(r1) rho_dn(r3) f_13 d3r1 d3r3 ; INT dOmega1 dOmega3 f = (4pi)^2 f0
    cross_f = (4 * np.pi) ** 2 * (rho * w2) @ f0 @ (rho * w2)
    S01_cross = 4 * cross_f                                    # 4 up-down pairs
    S01 = S01_intra + S01_cross
    print(f"\nS_01 = {S01:.4f}   (intra {S01_intra:.4f} + cross {S01_cross:.4f})")
    print(f"S_01/S_00 = {S01/S00:.6f}")

    # closed-form monopole Slater-Condon check: J^f0_aa + J^f0_bb + 4 J^f0_ab - 2 K^f0_ab
    def J0(P, Q):
        return (P * w2) @ f0 @ (Q * w2)
    Jaa = J0(A * A, A * A); Jbb = J0(B * B, B * B); Jab = J0(A * A, B * B); Kab = J0(A * B, A * B)
    sc = Jaa + Jbb + 4 * Jab - 2 * Kab
    print(f"Slater-Condon S_01/S_00 (monopole) = {sc:.6f}   diff {S01/S00 - sc:+.2e}")
    np.savez("debug/data/be_r12ci_analytic_stage1.npz", S00=S00, S01=S01)
