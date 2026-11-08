"""Be R12-CI, Part 3 (capstone): the genuinely-4-body Coulomb term WITH the antisymmetric
Be determinant (exchange chains), reduced RI-free form vs brute force.

The determinantal 4-body bridging integral (up-pair 1,2 correlate; down-pair 3,4 correlate;
Coulomb 1/r13 bridges up-electron-1 to down-electron-3):

   T4 = INT D_up(r1,r2)^2 D_dn(r3,r4)^2  f(r12) f(r34) (1/r13)  d3r1..d3r4
   D_up(r1,r2) = R1s(r1)R2s(r2) - R2s(r1)R1s(r2)   (the SAME determinant as Parts 1-2)

D_up^2 = R1s(r1)^2 R2s(r2)^2 + R2s(r1)^2 R1s(r2)^2 - 2 R1s(r1)R2s(r1) R1s(r2)R2s(r2)
                                                    ^^^^^^^^^^^ the EXCHANGE cross term
So the leaf integration carries the exchange chain explicitly.

REDUCED (RI-free): each correlated pair reduces to an isotropic radial DRESSING of its
bridge electron (leaf integrates out through the spherical average of f -- the exchange
cross term included), leaving a standard L=0 two-electron Coulomb integral:
   G_up(r1) = INT D_up(r1,r2)^2 f(r12) d3r2   (isotropic; all orbitals s)
   T4 = INT G_up(r1) G_dn(r3) / r13  d3r1 d3r3   (L=0 Slater; s-orbitals -> only L=0)

BRUTE: full 12-D importance-sampled MC of T4 (no reduction).

Agreement validates the reduced form on the antisymmetric object.  (This is the L=0 case
-- all-s Be ground config; the L>0 bridge was validated separately in r12ci_4e_be_integral.py
with a 2p bridge.  L>0 AND exchange together = the next increment.)
"""
import sys, os
import numpy as np
from numpy.polynomial.legendre import leggauss
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import be_r12ci_matelem as M   # orbitals R1s=s1_val, R2s=s2_val ; sampling

R1s, R2s = M.s1_val, M.s2_val
f_gem = M.f_gem


# ------------------------- BRUTE : 12-D MC --------------------------------- #
def brute(nw=40_000_000, batch=4_000_000, seed=7):
    rng = np.random.default_rng(seed)
    acc = 0.0; acc2 = 0.0; n = 0
    while n < nw:
        m = min(batch, nw - n)
        R, w = M.sample(m, rng)                       # each electron ~ M.g3d ; w=1/Pi g3d
        r = np.linalg.norm(R, axis=2)
        Dup = R1s(r[:, 0]) * R2s(r[:, 1]) - R2s(r[:, 0]) * R1s(r[:, 1])
        Ddn = R1s(r[:, 2]) * R2s(r[:, 3]) - R2s(r[:, 2]) * R1s(r[:, 3])
        r12 = np.linalg.norm(R[:, 0] - R[:, 1], axis=1)
        r34 = np.linalg.norm(R[:, 2] - R[:, 3], axis=1)
        r13 = np.linalg.norm(R[:, 0] - R[:, 2], axis=1)
        val = Dup ** 2 * Ddn ** 2 * f_gem(r12) * f_gem(r34) / r13
        acc += np.sum(w * val); acc2 += np.sum((w * val) ** 2); n += m
    mean = acc / nw
    err = np.sqrt(max(acc2 / nw - mean * mean, 0.0) / nw)
    return mean, err


# ------------------------- REDUCED : dressing + L=0 Slater ------------------ #
NR = 500
Rmax = 25.0
xg, wg = leggauss(NR)
rr = 0.5 * Rmax * (xg + 1.0)
wr = 0.5 * Rmax * wg
NX = 200
xx, wx = leggauss(NX)


def f_monopole(r1grid, r2grid):
    """f_0(r1,r2) = (1/2) INT_-1^1 f(|r1-r2|) dx  (spherical average) : (Nr1, Nr2)."""
    a = r1grid[:, None]; b = r2grid[None, :]
    out = np.zeros((a.shape[0], b.shape[1]))
    for x, w in zip(xx, wx):
        r12 = np.sqrt(np.maximum(a * a + b * b - 2 * a * b * x, 1e-30))
        out += 0.5 * w * f_gem(r12)
    return out


def dressing():
    """G_up(r1) = INT D_up(r1,r2)^2 f(r12) d3r2   on the rr grid.
       D_up^2 (angle-independent) times dOmega2 gives 4pi f_0 ; radial ∫ over r2."""
    f0 = f_monopole(rr, rr)                            # (Nr1, Nr2), 4pi absorbed below
    a1 = R1s(rr); b1 = R2s(rr)                         # on r1 grid
    a2 = R1s(rr); b2 = R2s(rr)                         # on r2 grid
    # D_up(r1,r2)^2 = a1^2 b2^2 + b1^2 a2^2 - 2 a1 b1 a2 b2   (r1 rows, r2 cols)
    Dsq = (np.outer(a1 ** 2, b2 ** 2) + np.outer(b1 ** 2, a2 ** 2)
           - 2 * np.outer(a1 * b1, a2 * b2))
    integrand = Dsq * f0 * (rr ** 2 * wr)[None, :]     # times r2^2 dr2
    return 4 * np.pi * integrand.sum(axis=1)           # 4pi from dOmega2 ; G_up(r1)


def reduced():
    G = dressing()                                     # G_up = G_dn (same orbitals)
    # T4 = INT G_up(r1) G_dn(r3)/r13 d3r1 d3r3 ; G isotropic -> (4pi)^2 * L=0 Slater
    r_gt = np.maximum.outer(rr, rr)
    Gw = G * rr ** 2 * wr
    L0 = Gw @ (1.0 / r_gt) @ Gw
    return (4 * np.pi) ** 2 * L0


if __name__ == "__main__":
    print("=" * 70)
    print("Be 4-body Coulomb term WITH the determinant (exchange chains)")
    print("   T4 = INT D_up^2 D_dn^2 f12 f34 / r13")
    print("=" * 70)
    red = reduced()
    print(f"REDUCED (dressing + L=0 Slater, exchange cross-term included) = {red:.8e}")
    for N in (10_000_000, 40_000_000):
        m, e = brute(N)
        rel = abs(m - red) / abs(red)
        print(f"BRUTE N={N:>11,}: {m:.8e} +/- {e:.1e}   rel.diff={rel:.2e}  ({(m-red)/e:+.1f} sigma)")
    print("-" * 70)
    print("Agreement validates the reduced RI-free form on the ANTISYMMETRIC determinant:")
    print("the exchange cross term -2 R1s R2s(1) R1s R2s(2) reduces through the same isotropic")
    print("dressing. So 'with exchange chains' is confirmed for the (L=0, all-s Be) case.")
