"""Two-center LiH R12-CI energy -- assembling the Be {Phi0, F Phi0} 2x2 (be_r12ci_full.py)
in the prolate two-center geometry, with the validated two-center 4-body reduction
(lih_r12_4body_integral.py sigma + lih_r12_4body_pi_channel.py pi/delta) living in g.

Ansatz (Be-style, ansatz B of debug/lih_r12_build_plan.md):
    Psi = Phi0 + c (F - Fbar) Phi0 ,   G = (F - Fbar) Phi0 ,   F = sum_{i<j} f(r_ij).
2x2 generalized eigenproblem
    H = [[E0, h],[h, g]] ,  S = [[1, 0],[0, sigma2]] ,
    E0 = <Phi0|H|Phi0>,  Fbar = <Phi0|F|Phi0>,
    sigma2 = <Phi0|(F-Fbar)^2|Phi0>,  h = <Phi0|H(F-Fbar)|Phi0>,  g = <G|H|G>.
The genuinely-4-body content is the chain f12 (1/r13) f34 inside g's V_ee part -- exactly the
integral validated this session.

Minimal reference (the two-center analog of Be's minimal 1s^2 2s^2): a 4-electron determinant
    Phi0 = |1s_A^2 1s_B^2|      (Li core on focus A, H-side pair on focus B)
with 1s_A (exponent za, tight, Li-core-like) and 1s_B (exponent zb, diffuse). A PoC of the
machinery, not a spectroscopic LiH (as the Be R12-CI was a 19%-of-correlation PoC).

STAGE 1 (this file): the f-integral primitives + Fbar (well-conditioned), validated vs MC.
Clean decomposition: every AO 2-body f-integral <p q|f|r s> = INT rho_{pr}(1) rho_{qs}(2) f
reduces to prolate QUADRATURE by dressing the ISOTROPIC member of a pair (Psi_aa, Psi_bb are
radial), EXCEPT (ab|f|ab) where the two-center product rho_ab sits on both sides -> one clean
6-D importance-MC.  The ill-conditioned sigma2 / h (quadrature) and g (the 4-body reduction)
are the next stages -- scaffolded at the end.

Geminal f(r) = exp(-GAM r).  Run from root:  python debug/lih_r12ci_energy.py
"""
import os
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lih_r12_4body_integral import (   # noqa: E402
    R, a, CENTER_A, CENTER_B, f_gem, sample_1s)

rng = np.random.default_rng(20260921)

# ---- minimal two-center model --------------------------------------------- #
Z_A, Z_B = 3.0, 1.0        # Li, H nuclear charges (used later for E0)
ZA, ZB = 2.70, 1.00        # 1s_A (Li-core-like), 1s_B (valence) orbital exponents
GAM = 0.50                 # geminal f(r) = exp(-GAM r)
N_A = np.sqrt(ZA ** 3 / np.pi)   # 1s Slater norm: |1s|^2 = (z^3/pi) e^{-2 z r}, INT=1
N_B = np.sqrt(ZB ** 3 / np.pi)

# ---- prolate grid --------------------------------------------------------- #
NXI, NETA = 240, 96
_xg, _wxg = leggauss(NXI)
_xi_max = 1.0 + 44.0 / (2 * min(ZA, ZB) * a)
XI1D = 1.0 + 0.5 * (_xg + 1.0) * (_xi_max - 1.0)
WXI = 0.5 * (_xi_max - 1.0) * _wxg
_eg, _weg = leggauss(NETA)
ETA1D = _eg.copy(); WETA = _weg.copy()
XI, ETA = np.meshgrid(XI1D, ETA1D, indexing='ij')
RA = a * (XI + ETA); RB = a * (XI - ETA); JAC = (XI ** 2 - ETA ** 2)


def _grid_int(fg):
    """INT fg dtau over the prolate grid (fg on (XI,ETA); phi-independent)."""
    return 2 * np.pi * a ** 3 * np.einsum('i,j,ij,ij->', WXI, WETA, JAC, fg)


rho_A = (ZA ** 3 / np.pi) * np.exp(-2 * ZA * RA)      # |1s_A|^2
rho_B = (ZB ** 3 / np.pi) * np.exp(-2 * ZB * RB)      # |1s_B|^2
orb_A = N_A * np.exp(-ZA * RA)                        # 1s_A
orb_B = N_B * np.exp(-ZB * RB)                        # 1s_B
rho_ab = orb_A * orb_B                                # 2-center transition density
S_AB = _grid_int(rho_ab)                              # <1s_A|1s_B>

# ---- isotropic f-dressing (Be-style spherical average; exact for a 1s) ----- #
_NR = 400; _Rmax = 30.0
_xr, _wr = leggauss(_NR); _rr = 0.5 * _Rmax * (_xr + 1.0); _wrr = 0.5 * _Rmax * _wr
_NX = 160; _xx, _wx = leggauss(_NX)


def f_dress_iso(s_vals, zeta):
    """Psi(s) = INT |1s(zeta)|^2(r2) f(|s-r2|) d3r2  (radial; 1s isotropic about its center)."""
    s = np.asarray(s_vals)[:, None]; r2 = _rr[None, :]
    f0 = np.zeros((s.shape[0], _rr.shape[0]))
    for x, w in zip(_xx, _wx):
        r12 = np.sqrt(np.maximum(s * s + r2 * r2 - 2 * s * r2 * x, 1e-30))
        f0 += 0.5 * w * f_gem(r12)
    radial = 4 * zeta ** 3 * _rr ** 2 * np.exp(-2 * zeta * _rr)   # |1s|^2 * 4pi r2^2
    return f0 @ (radial * _wrr)


Psi_aa = f_dress_iso(RA.ravel(), ZA).reshape(RA.shape)   # f-dressing of |1s_A|^2 (fn of r_A)
Psi_bb = f_dress_iso(RB.ravel(), ZB).reshape(RB.shape)   # f-dressing of |1s_B|^2 (fn of r_B)

# ---- AO 2-body f-integrals ------------------------------------------------- #
# <p q|f|r s> = INT rho_{pr}(1) rho_{qs}(2) f ; pair densities aa,bb,ab. Dress the isotropic one.
V_aaaa = _grid_int(rho_A * Psi_aa)                        # (aa|f|aa)
V_bbbb = _grid_int(rho_B * Psi_bb)                        # (bb|f|bb)
V_aabb = _grid_int(rho_A * Psi_bb)                        # (aa|f|bb) = INT rho_A Psi_bb
V_aaab = _grid_int(rho_ab * Psi_aa)                       # (aa|f|ab) = INT rho_ab Psi_aa
V_bbab = _grid_int(rho_ab * Psi_bb)                       # (bb|f|ab) = INT rho_ab Psi_bb


def mc_ff(sample1, sample2, n=24_000_000, batch=3_000_000):
    """<f> for e1~sample1(n), e2~sample2(n); bounded f -> low variance, batch-means error."""
    means = []; ntot = 0
    while ntot < n:
        r1 = sample1(batch); r2 = sample2(batch)
        means.append(f_gem(np.linalg.norm(r1 - r2, axis=1)).mean()); ntot += batch
    bm = np.array(means); return bm.mean(), bm.std(ddof=1) / np.sqrt(len(bm))


def s_A(n): return sample_1s(n, ZA, CENTER_A)
def s_B(n): return sample_1s(n, ZB, CENTER_B)


def mc_abab(n=24_000_000, batch=3_000_000):
    """(ab|f|ab) = INT rho_ab(1) rho_ab(2) f, by importance from q=|1s_A|^2:
    weight u = rho_ab/q = (N_B/N_A) e^{ZA r_A - ZB r_B} (finite variance since ZA>ZB)."""
    means = []; ntot = 0
    while ntot < n:
        r1 = s_A(batch); r2 = s_A(batch)
        rA1 = np.linalg.norm(r1 - CENTER_A, axis=1); rB1 = np.linalg.norm(r1 - CENTER_B, axis=1)
        rA2 = np.linalg.norm(r2 - CENTER_A, axis=1); rB2 = np.linalg.norm(r2 - CENTER_B, axis=1)
        u1 = (N_B / N_A) * np.exp(ZA * rA1 - ZB * rB1)
        u2 = (N_B / N_A) * np.exp(ZA * rA2 - ZB * rB2)
        d = np.linalg.norm(r1 - r2, axis=1)
        means.append((u1 * u2 * f_gem(d)).mean()); ntot += batch
    bm = np.array(means); return bm.mean(), bm.std(ddof=1) / np.sqrt(len(bm))


if __name__ == "__main__":
    print(f"[model] Z_A={Z_A} Z_B={Z_B} R={R}  za={ZA} zb={ZB}  <a|b>={S_AB:.6f}")
    print(f"[grid ] norm(rho_A)={_grid_int(rho_A):.6f}  norm(rho_B)={_grid_int(rho_B):.6f} (exact 1)")

    print("\n--- AO 2-body f-integrals (f=exp(-GAM r)): quadrature vs MC ---")
    checks = [("(aa|f|aa)", V_aaaa, mc_ff(s_A, s_A)),
              ("(bb|f|bb)", V_bbbb, mc_ff(s_B, s_B)),
              ("(aa|f|bb)", V_aabb, mc_ff(s_A, s_B))]
    for nm, q, (m, e) in checks:
        print(f"  {nm}: quad={q:.6f}  MC={m:.6f}+/-{e:.1e}  rel={abs(q - m) / abs(m):.2e}")
    print(f"  (aa|f|ab): quad={V_aaab:.6f}   (bb|f|ab): quad={V_bbab:.6f}   [1 iso dressing each]")
    Vabab, eab = mc_abab()
    print(f"  (ab|f|ab): importance-MC={Vabab:.6f}+/-{eab:.1e}   [2-center on both sides]")

    # ---- AO f-tensor over pair densities (0=aa,1=bb,2=ab), Loewdin, Fbar -----
    W = np.array([[V_aaaa, V_aabb, V_aaab],
                  [V_aabb, V_bbbb, V_bbab],
                  [V_aaab, V_bbab, Vabab]])
    dens = {('a', 'a'): 0, ('b', 'b'): 1, ('a', 'b'): 2, ('b', 'a'): 2}
    AOs = ['a', 'b']

    def Vf(mu, nu, la, si):     # <mu nu|f|la si> = W[pair(mu,la), pair(nu,si)]
        return W[dens[(mu, la)], dens[(nu, si)]]

    Smat = np.array([[1.0, S_AB], [S_AB, 1.0]])
    sval, svec = np.linalg.eigh(Smat)
    Xlow = svec @ np.diag(1 / np.sqrt(sval)) @ svec.T     # MO_p = sum_mu Xlow[mu,p] AO_mu

    def Vf_MO(p, q, r, s):
        tot = 0.0
        for mi, mu in enumerate(AOs):
            for ni, nu in enumerate(AOs):
                for li, la in enumerate(AOs):
                    for si, sg in enumerate(AOs):
                        tot += (Xlow[mi, p] * Xlow[ni, q] * Xlow[li, r] * Xlow[si, s]
                                * Vf(mu, nu, la, sg))
        return tot

    Jf = lambda p, q: Vf_MO(p, q, p, q)
    Kf = lambda p, q: Vf_MO(p, q, q, p)
    Fbar = Jf(0, 0) + Jf(1, 1) + 4 * Jf(0, 1) - 2 * Kf(0, 1)
    print(f"\n--- Fbar over the Loewdin determinant |m0^2 m1^2| ---")
    print(f"  J^f_00={Jf(0,0):.5f}  J^f_11={Jf(1,1):.5f}  J^f_01={Jf(0,1):.5f}  K^f_01={Kf(0,1):.5f}")
    print(f"  Fbar = <Phi0| sum_ij f_ij |Phi0> = {Fbar:.6f}")

    print("\n" + "-" * 76)
    print("STAGE 1 done: AO f-integral primitives (quad==MC on the 3 isotropic ones) + Fbar.")
    print("NEXT stages: E0 (determinant energy: T + V_ne + Coulomb J/K + V_NN, two-center);")
    print("  sigma2 = <F^2> - Fbar^2 (ill-conditioned -> quadrature, no cancellation);")
    print("  h = <Phi0|HF|Phi0> - Fbar*E0 (ill-conditioned -> quadrature);")
    print("  g = <G|H|G> (well-conditioned) with the 4-body chain f12(1/r13)f34 via the")
    print("  validated reduction (lih_r12_4body_integral / lih_r12_4body_pi_channel).")
