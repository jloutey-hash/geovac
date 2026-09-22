"""LiH R12-CI Stage 4b (part 2a): ANALYTIC (RI-free) g_Vne = <chi^2 V_ne>, chi=F-Fbar.

VMC target (lih_r12ci_vmc.py / _g_targets.py, geminal exp(-0.5 r)): g_Vne = -2.8454.
  g_Vne = <F^2 V_ne> - 2 Fbar <F V_ne> + Fbar^2 <V_ne>   (the <FV>,<V> already validated).
The NEW piece is the trilinear <F^2 V_ne>.  Block decomposition F = S + I, S=f12+f34 (intra),
I=f13+f14+f23+f24 (inter):
  <F^2 V> = <S^2 V> + 2<S I V> + <I^2 V> = A + B + C.
Each reduces to grid integrals of dressed one-electron densities (V_ne one-body -> a v-weight on
one electron; block independence factorizes the rest).  All three A,B,C are MC-validated inline.

Primitives:  P_pq (=m_p m_q), rho=P00+P11 (block marginal, INT=2), v(r)=-Z_A/rA-Z_B/rB (one-body),
  Psi^K_h = dress(K,h) = a^3 K@(geo h)  (K=build_kernel: axially-averaged kernel),
  S^f(1) = P00 Psi^f_11 + P11 Psi^f_00 - 2 P01 Psi^f_01   (f-dressed intra slice),
  rho_vpartner(1) = P00 <vP11> + P11 <vP00> - 2 P01 <vP01>  (v on the marginalized partner),
  G2(1) = 2 Psi^f_00 Psi^f_11 - 2 (Psi^f_01)^2   (two f-legs into the down block).

Run from debug/:  python lih_r12ci_gVne_analytic.py
"""
import numpy as np

from .kernels import (
    geo_f, dens, build_kernel, grid_int, a, ZA, ZB, GAM)
from .energy import Z_A, Z_B
from .hT import cvec, rA_f, rB_f, d_aa, d_ab, d_bb
from .basis import sample_block

GVNE_REF = -2.8454
FBAR_REF = 1.88028
VNE_REF = -20.86474
FVNE_REF = -40.17042

# --- primitives ---
v = -Z_A / rA_f - Z_B / rB_f
c00, c11, c01 = cvec(0, 0), cvec(1, 1), cvec(0, 1)
P00 = c00[0] * d_aa + c00[1] * d_ab + c00[2] * d_bb
P11 = c11[0] * d_aa + c11[1] * d_ab + c11[2] * d_bb
P01 = c01[0] * d_aa + c01[1] * d_ab + c01[2] * d_bb
rho = P00 + P11
Kf = build_kernel(GAM); Kf2 = build_kernel(2 * GAM)


def dress(K, h):
    return a ** 3 * (K @ (geo_f * h))


# f- and f^2-dressings of the standard densities
Pf = {'00': dress(Kf, P00), '11': dress(Kf, P11), '01': dress(Kf, P01), 'rho': dress(Kf, rho)}
Pf2 = {'00': dress(Kf2, P00), '11': dress(Kf2, P11), '01': dress(Kf2, P01), 'rho': dress(Kf2, rho)}
# f-dressings of the v-weighted densities
Pfv = {'00': dress(Kf, v * P00), '11': dress(Kf, v * P11), '01': dress(Kf, v * P01),
       'rho': dress(Kf, v * rho)}

# scalars
vP00i = grid_int(v * P00); vP11i = grid_int(v * P11); vP01i = grid_int(v * P01)
rho_vp = P00 * vP11i + P11 * vP00i - 2 * P01 * vP01i          # rho_vpartner
Pf_rhovp = dress(Kf, rho_vp)
vbar_blk = grid_int(v * rho)                                  # = <V_ne>/2
Ifac = lambda u, key: grid_int(u * Pf[key])                   # I_f[u, P_key]
alpha1 = Ifac(P00, '11') - Ifac(P01, '01')
alpha1_f2 = grid_int(P00 * Pf2['11']) - grid_int(P01 * Pf2['01'])

# f-dressed intra slices
S_f = P00 * Pf['11'] + P11 * Pf['00'] - 2 * P01 * Pf['01']
S_fv = P00 * Pfv['11'] + P11 * Pfv['00'] - 2 * P01 * Pfv['01']   # v on the dressed (partner) electron
# down-block two-leg fields
G2 = 2 * Pf['00'] * Pf['11'] - 2 * Pf['01'] ** 2
G2v = Pfv['00'] * Pf['11'] + Pfv['11'] * Pf['00'] - 2 * Pfv['01'] * Pf['01']   # v on one down leg


def analytic_F2Vne():
    # ---- A = <S^2 V> = 2<f12^2 V> + 2<f12 f34 V> ----
    VII2 = (grid_int(v * P00 * Pf2['11']) + grid_int(v * P11 * Pf2['00'])
            - 2 * grid_int(v * P01 * Pf2['01']))                 # <f12^2 (v1+v2)>
    VI = (grid_int(v * P00 * Pf['11']) + grid_int(v * P11 * Pf['00'])
          - 2 * grid_int(v * P01 * Pf['01']))                    # <f12 (v1+v2)>_up
    f12sqV = VII2 + alpha1_f2 * vbar_blk                          # <f12^2 V>
    f12f34V = 2 * alpha1 * VI                                     # <f12 f34 V>
    A = 2 * f12sqV + 2 * f12f34V

    # ---- B = 2<S I V> = 16 <f12 f13 V> ----
    f12f13_up = 0.25 * (grid_int(v * Pf['rho'] * S_f) + grid_int(Pf['rho'] * S_fv))   # (v1+v2)
    f12f13_dn = 0.25 * grid_int((Pfv['rho'] + Pf_rhovp) * S_f)                          # (v3+v4)
    B = 16 * (f12f13_up + f12f13_dn)

    # ---- C = <I^2 V> = 4<f13^2 V> + 8<f13 f14 V> + 4<f13 f24 V> ----
    f13sqV = 0.5 * (grid_int(v * rho * Pf2['rho']) + grid_int(rho_vp * Pf2['rho']))     # per pair
    C_same = 4 * f13sqV
    f13f14V = 0.25 * (grid_int(v * rho * G2) + grid_int(G2 * rho_vp) + 2 * grid_int(rho * G2v))
    C_shareUD = 8 * f13f14V                                       # share-up + share-down (equal)
    # disjoint <f13 f24 V>: D_p = sum kappa g(x) h(y); comps (g_field,h_field,g_key,h_key,kappa)
    dj = 0.0
    comps2 = [(P00, P11, '00', '11', 1.0), (P11, P00, '11', '00', 1.0), (P01, P01, '01', '01', -2.0)]
    for ga, ha, gk, hk, ka in comps2:
        for gb, hb, gbk, hbk, kb in comps2:
            Ig = grid_int(ga * Pf[gbk]); Ih = grid_int(ha * Pf[hbk])
            Igv = grid_int(v * ga * Pf[gbk]) + grid_int(ga * Pfv[gbk])   # I[vg_a,g_b]+I[g_a,vg_b]
            Ihv = grid_int(v * ha * Pf[hbk]) + grid_int(ha * Pfv[hbk])
            dj += ka * kb * (Igv * Ih + Ig * Ihv)
    f13f24V = 0.25 * dj
    C_dj = 4 * f13f24V
    C = C_same + C_shareUD + C_dj
    return A, B, C, dict(f12sqV=f12sqV, f12f34V=f12f34V, f13sqV=f13sqV, f13f14V=f13f14V,
                         f13f24V=f13f24V, VI=VI, alpha1=alpha1)


def analytic_FVne():
    """<F V_ne> = 2<f12 V> + 4<f13 V>  (cross-check vs -40.170)."""
    VI = (grid_int(v * P00 * Pf['11']) + grid_int(v * P11 * Pf['00']) - 2 * grid_int(v * P01 * Pf['01']))
    f12V = VI + alpha1 * vbar_blk
    f13V = 0.25 * (grid_int(v * rho * Pf['rho']) + grid_int(Pf['rho'] * rho_vp)
                   + grid_int((Pfv['rho'] + Pf_rhovp) * rho))
    return 2 * f12V + 4 * f13V


def mc_breakdown(seed_u=41, seed_d=42):
    cfg = dict(nw=12000, burn=3000, nsnap=16, thin=45)
    up, _ = sample_block(seed=seed_u, **cfg); dn, _ = sample_block(seed=seed_d, **cfg)
    P = [(0, 1), (2, 3)]; INTER = [(0, 2), (0, 3), (1, 2), (1, 3)]
    A_, B_, C_, F2V_, FV_, V_, chi2V_ = [], [], [], [], [], [], []
    Fall = []
    POS = []
    for s in range(len(up)):
        pos = np.concatenate([up[s], dn[s]], axis=1); POS.append(pos)
        d = lambda i, j: np.linalg.norm(pos[:, i] - pos[:, j], axis=-1)
        Fall.append(sum(np.exp(-GAM * d(i, j)) for i, j in P + INTER))
    Fbar = np.concatenate(Fall).mean()
    for pos in POS:
        d = lambda i, j: np.maximum(np.linalg.norm(pos[:, i] - pos[:, j], axis=-1), 1e-12)
        Su = sum(np.exp(-GAM * d(i, j)) for i, j in P)
        Iu = sum(np.exp(-GAM * d(i, j)) for i, j in INTER)
        F = Su + Iu
        rA = np.linalg.norm(pos - np.array([0, 0, -a]), axis=-1)
        rB = np.linalg.norm(pos - np.array([0, 0, +a]), axis=-1)
        Vne = (-Z_A / np.maximum(rA, 1e-12) - Z_B / np.maximum(rB, 1e-12)).sum(axis=1)
        A_.append((Su ** 2 * Vne).mean()); B_.append((2 * Su * Iu * Vne).mean())
        C_.append((Iu ** 2 * Vne).mean()); F2V_.append((F ** 2 * Vne).mean())
        FV_.append((F * Vne).mean()); V_.append(Vne.mean())
        chi2V_.append(((F - Fbar) ** 2 * Vne).mean())
    m = lambda x: np.mean(x); e = lambda x: np.std(x, ddof=1) / np.sqrt(len(x))
    return dict(Fbar=Fbar, A=(m(A_), e(A_)), B=(m(B_), e(B_)), C=(m(C_), e(C_)),
                F2V=(m(F2V_), e(F2V_)), FV=(m(FV_), e(FV_)), V=(m(V_), e(V_)),
                chi2V=(m(chi2V_), e(chi2V_)))


if __name__ == "__main__":
    print("=" * 78)
    print("LiH R12-CI Stage 4b (2a): ANALYTIC g_Vne = <chi^2 V_ne>, target -2.8454")
    print("=" * 78)
    A, B, C, info = analytic_F2Vne()
    F2V = A + B + C
    FVne = analytic_FVne()
    Vne = 2 * vbar_blk
    print(f"\n  scalars: alpha1={alpha1:.5f}  alpha1_f2={alpha1_f2:.5f}  vbar_blk={vbar_blk:.5f}")
    print(f"  <V_ne>  (analytic) = {Vne:.5f}   vs {VNE_REF:.5f}")
    print(f"  <F V_ne>(analytic) = {FVne:.5f}   vs {FVNE_REF:.5f}")
    print(f"\n  <F^2 V_ne> parts:  A(S^2 V)={A:+.4f}  B(2 S I V)={B:+.4f}  C(I^2 V)={C:+.4f}")
    print(f"  <F^2 V_ne> (analytic) = A+B+C = {F2V:+.5f}")

    Fbar = FBAR_REF
    gVne = F2V - 2 * Fbar * FVne + Fbar ** 2 * Vne
    print(f"\n  >>> g_Vne = <F^2 V> - 2 Fbar <F V> + Fbar^2 <V> = {gVne:+.5f}   target {GVNE_REF:+.4f}")

    print("\n  --- MC validation (block sampler) ---")
    mc = mc_breakdown()
    print(f"  Fbar(MC)={mc['Fbar']:.5f}")
    for k in ('A', 'B', 'C', 'F2V', 'FV', 'V', 'chi2V'):
        print(f"    {k:6s} MC = {mc[k][0]:+.4f} +/- {mc[k][1]:.1e}")
    print(f"  analytic A={A:+.4f} B={B:+.4f} C={C:+.4f}  F2V={F2V:+.4f}  FV={FVne:+.4f}  V={Vne:+.4f}")
    print(f"  analytic g_Vne={gVne:+.5f}   MC chi2V={mc['chi2V'][0]:+.5f}+/-{mc['chi2V'][1]:.1e}")
