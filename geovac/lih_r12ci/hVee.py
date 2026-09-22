"""LiH R12-CI Stage 4b (part 1): ANALYTIC (RI-free) h_Vee = <Phi0| V_ee (F-Fbar) |Phi0>.

h_Vee = Cov[F, V_ee] under the block density |Phi0|^2, VMC target +0.3082 (lih_r12ci_vmc.py,
geminal exp(-0.5 r)).  It is the SAME separable-density bilinear as sigma^2, but with the
SECOND operator = the Coulomb kernel 1/r instead of a second f:

  <F_A F_B> = 2 alpha1^C + 4 beta1^C                                   [SAME pair, C = A.B]
            + 2 alpha1^A alpha1^B                                       [DISJOINT intra-intra]
            + 2 [ INT S^A Psi^B_rho + INT S^B Psi^A_rho ]               [SHARE-vertex intra x inter]
            + 2 INT rho ( Psi^A_00 Psi^B_11 + Psi^A_11 Psi^B_00 - 2 Psi^A_01 Psi^B_01 ) [SHARE inter-inter]
            + (t1+t2+t3+t4)                                             [DISJOINT inter-inter]
  Cov[F_A,F_B] = <F_A F_B> - (2a1^A+4b1^A)(2a1^B+4b1^B).
  S^K(1) = P00 Psi^K_11 + P11 Psi^K_00 - 2 P01 Psi^K_01   (K-dressed intra slice).
  ti = 1/4 sum_ab ka kb (g/h .W_A. g/h)(h/g .W_B. h/g)  (the 4 crossed disjoint contractions).

For A=B=f this reproduces the sigma^2 code term-by-term (verified below: cov[f,f]==sigma^2).
For h_Vee: A=f, B=Coulomb, C=A.B=Yukawa (e^{-gam r}/r).  Kernel objects needed:
  W^f, Psi^f     (build_kernel, coarse grid)          -- have (sigma^2 machinery)
  W^Y            (Yukawa AO-pair matrix)               -- same-pair C; from hT module
  W^coul, Psi^coul  (aa,bb = closed Hartree; ab = prolate-Neumann POTENTIAL field) -- NEW here.

Run from debug/:  python lih_r12ci_hVee_analytic.py
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import lqn, eval_legendre

from .kernels import (
    Xg, Eg, rA, rB, geo_f, dens, build_kernel, dressings, Wmat, grid_int, X, a, ZA, ZB, GAM,
    XI as XI1D, ETA as ETA1D, WXI, WETA)
from .energy import _hartree_1s, R
from .hT import (yukawa_pot_iso, mc_yukawa_abab, rA_f, rB_f, d_aa, d_ab, d_bb, cvec)
from .basis import sample_block

HVEE_REF = 0.3082      # VMC target (lih_r12ci_vmc.py), geminal exp(-0.5 r)
E2_REF = 3.614436      # <Phi0|V_ee|Phi0> (Stage 2)
SIG2_REF = 0.13699     # sigma^2 self-check target
AOP = ['aa', 'ab', 'bb']

c00, c11, c01 = cvec(0, 0), cvec(1, 1), cvec(0, 1); crho = c00 + c11
P00 = c00[0] * d_aa + c00[1] * d_ab + c00[2] * d_bb
P11 = c11[0] * d_aa + c11[1] * d_ab + c11[2] * d_bb
P01 = c01[0] * d_aa + c01[1] * d_ab + c01[2] * d_bb
rho_g = P00 + P11
# disjoint-comps for gam_dj: (g, h, kappa) with D_p = sum kappa g(x) h(y)
_COMPS = [(c00, c11, 1.0), (c11, c00, 1.0), (c01, c01, -2.0)]

# --------------------------------------------------------------------------- #
# Coulomb dressing field of the transition density rho_ab, via the prolate-Neumann
# POTENTIAL (extracted from the same expansion as _prolate_neumann_coulomb).
#   V_D(xi,eta) = (2/R)(2pi) a^3 sum_l (2l+1) P_l(eta) [ sum_xi' wxi K_l(xi,xi') g_l(xi') ],
#   g_l(xi') = sum_eta' weta JAC(xi',eta') D(xi',eta') P_l(eta'),  K_l = P_l(xi_<) Q_l(xi_>).
# --------------------------------------------------------------------------- #
NXI, NETA = XI1D.size, ETA1D.size                             # coarse-grid axes + weights (imported)
JAC2 = Xg ** 2 - Eg ** 2                                       # (NXI,NETA)


def neumann_potential(Dg2, LMAX=34):
    """Coulomb potential field of density Dg2(xi,eta) [phi-indep] on the coarse grid."""
    Qtab = np.array([lqn(LMAX, x)[0] for x in XI1D])          # (NXI, LMAX+1)
    Pxi = np.array([eval_legendre(l, XI1D) for l in range(LMAX + 1)])   # (LMAX+1, NXI)
    minidx = np.minimum.outer(np.arange(NXI), np.arange(NXI))
    maxidx = np.maximum.outer(np.arange(NXI), np.arange(NXI))
    pref = (2.0 / R) * (2 * np.pi) * a ** 3
    W = JAC2 * Dg2
    V = np.zeros((NXI, NETA))
    for l in range(LMAX + 1):
        Pl_eta = eval_legendre(l, ETA1D)
        g_l = (W * Pl_eta[None, :]) @ WETA                     # (NXI,)
        K = Pxi[l][minidx] * Qtab[:, l][maxidx]                # (NXI,NXI) P_l(xi<)Q_l(xi>)
        radial = K @ (WXI * g_l)                               # (NXI,)
        V += (2 * l + 1) * np.outer(radial, Pl_eta)            # (NXI,NETA)
    return pref * V


# --------------------------------------------------------------------------- #
# generic kernel object: W (3x3, [aa,ab,bb]) + dressing dict {aa,ab,bb} (grid, flat)
# --------------------------------------------------------------------------- #
def make_kernel_f(gam):
    K = build_kernel(gam); W, Psi = Wmat(K)                    # Psi dict flat
    return dict(W=W, Psi=Psi)


def make_kernel_coul(LMAX=34):
    Vaa = _hartree_1s(rA_f, ZA); Vbb = _hartree_1s(rB_f, ZB)   # closed-form isotropic potentials
    dab2 = (d_ab).reshape(NXI, NETA)
    Vab = neumann_potential(dab2, LMAX).reshape(-1)            # transition-density potential
    Psi = {'aa': Vaa, 'ab': Vab, 'bb': Vbb}
    W = np.zeros((3, 3))
    for i, u in enumerate(AOP):
        for j, v in enumerate(AOP):
            W[i, j] = grid_int(dens[u] * Psi[v])
    return dict(W=W, Psi=Psi)


def make_kernel_Y(gam):
    """Yukawa AO-pair matrix (same-pair C=f.coul); no dressing field needed (only SAME terms)."""
    VYA = yukawa_pot_iso(rA_f, ZA, gam); VYB = yukawa_pot_iso(rB_f, ZB, gam)
    IY = {'aaaa': grid_int(d_aa * VYA), 'bbbb': grid_int(d_bb * VYB),
          'aabb': grid_int(d_bb * VYA), 'aaab': grid_int(d_ab * VYA), 'bbab': grid_int(d_ab * VYB)}
    abab, _ = mc_yukawa_abab(gam)
    W = np.array([[IY['aaaa'], IY['aaab'], IY['aabb']],
                  [IY['aaab'], abab, IY['bbab']],
                  [IY['aabb'], IY['bbab'], IY['bbbb']]])         # [aa,ab,bb]
    return dict(W=W)


def _a1b1(W):
    a1 = c00 @ W @ c11 - c01 @ W @ c01
    b1 = 0.25 * (crho @ W @ crho)
    return a1, b1


def _Sfield(Psi):
    """S^K(1) = P00 Psi_11 + P11 Psi_00 - 2 P01 Psi_01  (K-dressed intra slice)."""
    Pf = lambda c: c[0] * Psi['aa'] + c[1] * Psi['ab'] + c[2] * Psi['bb']
    return P00 * Pf(c11) + P11 * Pf(c00) - 2 * P01 * Pf(c01)


def _psi(Psi, c):
    return c[0] * Psi['aa'] + c[1] * Psi['ab'] + c[2] * Psi['bb']


def cov_FA_FB(A, B, C):
    """Cov[F_A, F_B] over |Phi0|^2.  A,B carry W+Psi; C carries W for the product kernel A.B."""
    a1A, b1A = _a1b1(A['W']); a1B, b1B = _a1b1(B['W'])
    a1C, b1C = _a1b1(C['W'])
    same = 2 * a1C + 4 * b1C
    dj_intra = 2 * a1A * a1B
    # SHARE intra x inter
    SA, SB = _Sfield(A['Psi']), _Sfield(B['Psi'])
    PsiA_rho, PsiB_rho = _psi(A['Psi'], crho), _psi(B['Psi'], crho)
    share_delta = 2 * (grid_int(SA * PsiB_rho) + grid_int(SB * PsiA_rho))
    # SHARE inter x inter
    A00, A11, A01 = _psi(A['Psi'], c00), _psi(A['Psi'], c11), _psi(A['Psi'], c01)
    B00, B11, B01 = _psi(B['Psi'], c00), _psi(B['Psi'], c11), _psi(B['Psi'], c01)
    share_gsh = 2 * grid_int(rho_g * (A00 * B11 + A11 * B00 - 2 * A01 * B01))
    # DISJOINT inter x inter (4 crossed contractions)
    WA, WB = A['W'], B['W']
    dj = 0.0
    for ga, ha, ka in _COMPS:
        for gb, hb, kb in _COMPS:
            t1 = (ga @ WA @ gb) * (ha @ WB @ hb)
            t2 = (ha @ WA @ hb) * (ga @ WB @ gb)
            t3 = (ga @ WA @ hb) * (ha @ WB @ gb)
            t4 = (ha @ WA @ gb) * (ga @ WB @ hb)
            dj += ka * kb * (t1 + t2 + t3 + t4)
    dj *= 0.25
    tot = same + dj_intra + share_delta + share_gsh + dj
    FbarA = 2 * a1A + 4 * b1A; FbarB = 2 * a1B + 4 * b1B
    return tot - FbarA * FbarB, dict(same=same, dj_intra=dj_intra, share_delta=share_delta,
                                     share_gsh=share_gsh, dj=dj, FbarA=FbarA, FbarB=FbarB)


def vee_mc(seed_u=21, seed_d=22):
    """MC of <F.V_ee> and <V_ee> over |Phi0|^2 (independent validation of h_Vee)."""
    cfg = dict(nw=10000, burn=3000, nsnap=12, thin=40)
    up, _ = sample_block(seed=seed_u, **cfg); dn, _ = sample_block(seed=seed_d, **cfg)
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    fv, vv = [], []
    for s in range(len(up)):
        pos = np.concatenate([up[s], dn[s]], axis=1)
        d = {ij: np.maximum(np.linalg.norm(pos[:, ij[0]] - pos[:, ij[1]], axis=-1), 1e-12) for ij in pairs}
        F = sum(np.exp(-GAM * d[ij]) for ij in pairs)
        Vee = sum(1.0 / d[ij] for ij in pairs)
        fv.append((F * Vee).mean()); vv.append(Vee.mean())
    fv = np.array(fv); vv = np.array(vv)
    return fv.mean(), fv.std(ddof=1) / np.sqrt(len(fv)), vv.mean(), vv.std(ddof=1) / np.sqrt(len(vv))


if __name__ == "__main__":
    print("=" * 78)
    print("LiH R12-CI Stage 4b (part 1): ANALYTIC h_Vee (RI-free), target +0.3082")
    print("=" * 78)

    Af = make_kernel_f(GAM)                 # f
    Cf2 = make_kernel_f(2 * GAM)            # f^2 (product kernel for the sigma^2 self-check)
    Kcoul = make_kernel_coul()             # Coulomb
    Ky = make_kernel_Y(GAM)                # Yukawa = f.coul (product kernel for h_Vee)

    # ---- GATE 1: Coulomb dressing / W^coul controls ----
    # Neumann potential of the ISOTROPIC rho_aa must reproduce the closed-form Hartree.
    Vaa_neu = neumann_potential(d_aa.reshape(NXI, NETA)).reshape(-1)
    Vaa_cf = _hartree_1s(rA_f, ZA)
    msk = d_aa > 1e-8 * d_aa.max()
    print(f"\n  GATE 1 -- Neumann potential of rho_aa vs closed Hartree: "
          f"weighted rel {np.abs((Vaa_neu-Vaa_cf))[msk].mean()/np.abs(Vaa_cf[msk]).mean():.1e}")
    Wc = Kcoul['W']
    print(f"    W^coul[aa,aa]={Wc[0,0]:.5f} vs 5ZA/8={5*ZA/8:.5f}   "
          f"W^coul[bb,bb]={Wc[2,2]:.5f} vs 5ZB/8={5*ZB/8:.5f}   W^coul[ab,ab]={Wc[1,1]:.5f}")
    a1c, b1c = _a1b1(Wc)
    print(f"    E2 = <V_ee> = 2 a1^coul + 4 b1^coul = {2*a1c+4*b1c:.5f}  vs Stage-2 {E2_REF:.5f}")

    # ---- GATE 2: sigma^2 self-check: cov[f,f] must reproduce 0.137 ----
    sig2, _ = cov_FA_FB(Af, Af, Cf2)
    print(f"\n  GATE 2 -- cov[F_f,F_f] (self-check) = {sig2:.6f}  vs sigma^2 {SIG2_REF:.6f}"
          f"   (rel {abs(sig2-SIG2_REF)/SIG2_REF:.1e})")

    # ---- h_Vee = cov[F_f, F_coul] ----
    hVee, info = cov_FA_FB(Af, Kcoul, Ky)
    print(f"\n  --- h_Vee = Cov[F_f, F_coul] ---")
    for k in ('same', 'dj_intra', 'share_delta', 'share_gsh', 'dj'):
        print(f"    {k:12s} = {info[k]:+.6f}")
    print(f"    Fbar_f={info['FbarA']:.5f}  E2(=Fbar_coul)={info['FbarB']:.5f}")
    print(f"\n  >>> h_Vee = {hVee:+.6f}   VMC target = {HVEE_REF:+.4f}   (dev {abs(hVee-HVEE_REF):.4f})")

    # ---- MC validation ----
    fv, efv, vv, evv = vee_mc()
    Fbar_f = info['FbarA']
    hVee_mc = fv - Fbar_f * vv
    print(f"\n  GATE 3 -- MC: <F.V_ee>={fv:.4f}+/-{efv:.1e}  <V_ee>={vv:.4f}+/-{evv:.1e}")
    print(f"    h_Vee(MC) = <F.V_ee> - Fbar <V_ee> = {hVee_mc:+.5f}   (analytic {hVee:+.5f})")
