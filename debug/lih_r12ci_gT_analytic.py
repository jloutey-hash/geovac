"""LiH R12-CI Stage 4b (part 2c): ANALYTIC (RI-free) g_T = <G|T|G>, G=(F-Fbar)Phi0.

The kinetic part of the diagonal element g.  Via the drift form + a Green/IBP identity
(the SAME identity that collapsed h_T's Part B), g_T reduces to ONE genuinely-new object
(the drift^2 term gT1) plus scalars -- NO vector nabla-f dressing is needed:

  g_T = 1/2 <chi^2 sum_i v_i^2> + <chi sum_i grad_i F . v_i> + 1/2 <sum_i |grad_i F|^2>
      = gT1 + gT2 + gT3,      v_i = grad_i Phi0 / Phi0 .
  IBP:  gT2 + gT3 = -1/2 <chi sum_i grad^2 F> = -gam^2 sigma^2 + 2 gam Cov[F, Y_sum],
    Y_sum = sum_{i<j} e^{-gam r_ij}/r_ij (Yukawa),  grad^2 f = gam^2 f - 2 gam (e^{-gam r}/r).
  =>  g_T = gT1 - gam^2 sigma^2 + 2 gam Cov[F, Y_sum].

Only gT1 = 1/2 <chi^2 sum_i v_i^2> is new.  sum_i |grad_i Phi0|^2 = (1/4)[KD(1,2)D_p(3,4)
+ D_p(1,2)KD(3,4)], KD the gradient pair-density (G_pq = grad m_p . grad m_q).  So
  <F^k sum v^2> = (1/2) sum_{k f-edges} sum_{cU in KD, cD in D_p} wt * contract(f-edges over the
                 blocks with block-U density = KD, block-D = D_p),
a 2-f reduction with NO Coulomb -> NO triangle (all separable).  All of S2=<sum v^2>,
FS=<F sum v^2>, F2S=<F^2 sum v^2> from the SAME reducer -> the ~2% cancellation in
gT1 = 1/2[F2S - 2 Fbar FS + Fbar^2 S2] is grid-consistent.  Cov[F,Y_sum] = cov_FA_FB with
Yukawa kernels (h_Vee structure, Coulomb -> Yukawa(gam), product kernel Yukawa(2gam)).

Run from debug/:  python lih_r12ci_gT_analytic.py
"""
import numpy as np

from lih_r12ci_sigma2_analytic import a, GAM, grid_int, dens, build_kernel, geo_f
from lih_r12ci_hVee_analytic import (
    P00, P11, P01, rho_g, cov_FA_FB, make_kernel_f, make_kernel_Y, SIG2_REF)
from lih_r12ci_hT_analytic import Gpq, cvec, yukawa_pot_iso, rA_f, rB_f, d_aa, d_ab, d_bb, ZA, ZB
from lih_r12ci_gVee_analytic import psi_yuk, dress, Kf, Kf2, PAIRS, COMPS

GT_REF = 1.38640            # VMC g_T (=_g_targets.py)
GT1_REF = 1.16826           # drift^2 target
SIG2 = SIG2_REF             # 0.13699

# gradient one-electron densities + kinetic overlaps
G00, G11, G01 = Gpq(0, 0), Gpq(1, 1), Gpq(0, 1)
T00, T11, T01 = grid_int(G00), grid_int(G11), grid_int(G01)

# KD (gradient pair-density) AO-pair components on block U: (g, h, kappa)
KD_COMPS = [(G00, P11, 1.0), (G11, P00, 1.0), (G01, P01, -2.0),
            (P11, G00, 1.0), (P00, G11, 1.0), (P01, G01, -2.0)]
# D_p components on block D (from hVee): (g, h, kappa)
DP_COMPS = [(P00, P11, 1.0), (P11, P00, 1.0), (P01, P01, -2.0)]


def contract_raw(fedges, base):
    """INT (prod_edges f) * prod_i base_i d(all 4 electrons)  (no kept electron; leaf reduction)."""
    D = dict(base); scalar = 1.0; alive = {1, 2, 3, 4}
    edges = [set(e) for e in fedges]

    def dnb(u):
        s = set()
        for ed in edges:
            if u in ed:
                s |= (ed - {u})
        return s
    while alive:
        u = next((c for c in alive if len(dnb(c)) <= 1), None)
        if u is None:
            raise RuntimeError(f"unexpected 2f topology: {fedges}")
        inc = [ed for ed in edges if u in ed]; nb = dnb(u)
        if not nb:
            scalar *= grid_int(D[u])
        else:
            w = next(iter(nb)); mult = len(inc)
            D[w] = D[w] * dress(Kf if mult == 1 else Kf2, D[u])
        alive.discard(u); edges = [ed for ed in edges if u not in ed]
    return scalar


def FkSv2(k):
    """<F^k sum_i v_i^2> = (1/2) sum_{k ordered f-edges} sum_{cU in KD, cD in Dp} wt * contract."""
    if k == 0:
        edge_sets = [[]]
    elif k == 1:
        edge_sets = [[p] for p in PAIRS]
    else:
        edge_sets = [[p, q] for p in PAIRS for q in PAIRS]
    tot = 0.0
    for edges in edge_sets:
        for gU, hU, kU in KD_COMPS:
            for gD, hD, kD in DP_COMPS:
                base = {1: gU, 2: hU, 3: gD, 4: hD}
                tot += kU * kD * contract_raw(edges, base)
    return 0.5 * tot


def make_kernel_Y_full(gam):
    """Yukawa(gam) kernel with BOTH W (from make_kernel_Y) and the dressing field Psi (psi_yuk)."""
    K = make_kernel_Y(gam)
    Psi = {'aa': psi_yuk(d_aa, gam), 'ab': psi_yuk(d_ab, gam), 'bb': psi_yuk(d_bb, gam)}
    return dict(W=K['W'], Psi=Psi)


if __name__ == "__main__":
    print("=" * 78)
    print("LiH R12-CI g_T (analytic, RI-free) = gT1 - gam^2 sigma^2 + 2 gam Cov[F,Y_sum]")
    print(f"  f=exp(-{GAM} r), gam={GAM}")
    print("=" * 78)

    # ---- Fbar (grid) ----
    Af = make_kernel_f(GAM)
    from lih_r12ci_hVee_analytic import _a1b1, crho, c00, c11, c01
    a1, b1 = _a1b1(Af['W']); Fbar = 2 * a1 + 4 * b1
    print(f"\n  Fbar = {Fbar:.6f}   T00={T00:.5f} T11={T11:.5f}  <T>=T00+T11={T00+T11:.5f}")

    # ---- gT1 = 1/2 [F2S - 2 Fbar FS + Fbar^2 S2], all from the same reducer ----
    S2 = FkSv2(0); FS = FkSv2(1); F2S = FkSv2(2)
    print(f"\n  <sum v^2>      S2  = {S2:.5f}   (=2<T>={2*(T00+T11):.5f})")
    print(f"  <F sum v^2>    FS  = {FS:.5f}")
    print(f"  <F^2 sum v^2>  F2S = {F2S:.5f}")
    gT1 = 0.5 * (F2S - 2 * Fbar * FS + Fbar ** 2 * S2)
    print(f"  gT1 = 1/2[F2S - 2 Fbar FS + Fbar^2 S2] = {gT1:+.6f}   (VMC target {GT1_REF:+.5f})")

    # ---- Cov[F, Y_sum] via cov_FA_FB (Coulomb -> Yukawa(gam); product kernel Yukawa(2gam)) ----
    By = make_kernel_Y_full(GAM)                 # second operator = Yukawa(gam), needs W + Psi
    Cy2 = make_kernel_Y(2 * GAM)                 # product kernel f*Y = Yukawa(2gam), needs W only
    CovFY, cinfo = cov_FA_FB(Af, By, Cy2)
    Ybar = cinfo['FbarB']
    print(f"\n  Cov[F, Y_sum] = {CovFY:+.6f}   (Ybar = <Y_sum> = {Ybar:.5f})")

    # ---- assemble g_T ----
    gT23 = -GAM ** 2 * SIG2 + 2 * GAM * CovFY     # = gT2 + gT3
    g_T = gT1 + gT23
    print(f"\n  gT2+gT3 = -gam^2 sigma^2 + 2 gam Cov[F,Y_sum] = {gT23:+.6f}   (VMC {-0.04750+0.26563:+.5f})")
    print(f"  >>> g_T = gT1 + (gT2+gT3) = {g_T:+.6f}   (VMC target {GT_REF:+.5f})")
