"""g_Vee groundwork: de-risk the two NEW pieces + get MC targets, before the full assembly.

(1) YUKAWA DRESSING FIELD Psi^Y_h(r) = INT h(r') e^{-gam|r-r'|}/|r-r'| dr', WITHOUT a Yukawa-Neumann:
    Psi^Y_h = Psi^coul_h - Psi^{gs}_h ,  gs(d) = (1 - e^{-gam d})/d  (SMOOTH kernel, -> gam at d=0,
    no grid diagonal singularity).  Psi^coul via neumann_potential (h_Vee) / closed Hartree (iso).
    GATE: Psi^Y_aa (this route) == yukawa_pot_iso(rA,ZA,gam) (closed form).
(2) The 4-BODY BRIDGE inside the energy: <f12 f34 coul_13> = 1/4 I_coul[S^f, S^f], S^f the
    f-dressed intra slice, I_coul via neumann_potential(S^f).  GATE: analytic == block-sampler MC.
(3) MC 6-product breakdown of <F^2 V_ee> = <S^2 C_S>+<S^2 C_I>+2<S I C_S>+2<S I C_I>+<I^2 C_S>+<I^2 C_I>
    (S=f12+f34, I=inter f, C_S=coul_12+coul_34, C_I=inter coul) -> analytic targets.

Run from debug/:  python lih_r12ci_gVee_groundwork.py
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

from lih_r12ci_sigma2_analytic import (
    Xg, Eg, geo_f, dens, build_kernel, grid_int, a, ZA, ZB, GAM,
    rho_cyl_f, zc_f)
from lih_r12ci_energy import Z_A, Z_B, _hartree_1s
from lih_r12ci_hT_analytic import cvec, rA_f, rB_f, d_aa, d_ab, d_bb, yukawa_pot_iso
from lih_r12ci_hVee_analytic import neumann_potential, NXI, NETA
from lih_r12ci_sigma2_mc import sample_block

# --- densities ---
c00, c11, c01 = cvec(0, 0), cvec(1, 1), cvec(0, 1)
P00 = c00[0] * d_aa + c00[1] * d_ab + c00[2] * d_bb
P11 = c11[0] * d_aa + c11[1] * d_ab + c11[2] * d_bb
P01 = c01[0] * d_aa + c01[1] * d_ab + c01[2] * d_bb
rho = P00 + P11
Kf = build_kernel(GAM)


def dress(K, h):
    return a ** 3 * (K @ (geo_f * h))


def build_kernel_smooth(gam):
    """K[i,j] = INT_0^2pi (1 - e^{-gam d})/d dphi'  (d = |r_i - r_j(phi')|); smooth, no singularity."""
    NPHI = 40
    xp, wp = leggauss(NPHI)
    phi = 0.5 * np.pi * (xp + 1.0); wphi = 0.5 * np.pi * wp
    rc = rho_cyl_f; z = zc_f
    NG = rc.size
    K = np.zeros((NG, NG))
    rc2 = rc[:, None] * rc[None, :]
    base = rc[:, None] ** 2 + rc[None, :] ** 2 + (z[:, None] - z[None, :]) ** 2
    for cphi, w in zip(np.cos(phi), wphi):
        d = np.sqrt(np.maximum(base - 2.0 * rc2 * cphi, 1e-24))
        K += 2.0 * w * (1.0 - np.exp(-gam * d)) / d
    return K


Ks = build_kernel_smooth(GAM)


def psi_coul(h):
    """Coulomb potential field of density h (grid flat) via prolate Neumann."""
    return neumann_potential(h.reshape(NXI, NETA)).reshape(-1)


def psi_yuk(h, gam=GAM):
    """Yukawa (screened-Coulomb) potential field: Psi^coul - Psi^{smooth}."""
    return psi_coul(h) - dress(Ks, h)


if __name__ == "__main__":
    print("=" * 78)
    print("g_Vee GROUNDWORK: Yukawa field + 4-body bridge + MC 6-product breakdown")
    print("=" * 78)

    # ---- GATE 1: Yukawa field of the isotropic rho_aa vs closed yukawa_pot_iso ----
    PsiY_aa = psi_yuk(d_aa)
    yk_closed = yukawa_pot_iso(rA_f, ZA, GAM)
    msk = d_aa > 1e-6 * d_aa.max()
    relY = np.abs(PsiY_aa - yk_closed)[msk].mean() / np.abs(yk_closed[msk]).mean()
    print(f"\n  GATE 1 -- Psi^Y_aa (Coul - smooth) vs yukawa_pot_iso closed form: weighted rel {relY:.1e}")
    # also I_Y[ab,ab] via the field vs the h_T MC value 0.00456
    PsiY_ab = psi_yuk(d_ab)
    IYabab = grid_int(d_ab * PsiY_ab)
    print(f"    I_Y[ab,ab] = grid_int(rho_ab Psi^Y_ab) = {IYabab:.5f}   (h_T importance-MC gave 0.00456)")

    # ---- GATE 2: 4-body bridge <f12 f34 coul_13> = 1/4 I_coul[S^f, S^f] ----
    Pf = {'00': dress(Kf, P00), '11': dress(Kf, P11), '01': dress(Kf, P01)}
    S_f = P00 * Pf['11'] + P11 * Pf['00'] - 2 * P01 * Pf['01']     # f-dressed intra slice
    bridge_an = 0.25 * grid_int(S_f * psi_coul(S_f))
    print(f"\n  GATE 2 -- 4-body bridge <f12 f34 coul_13> (analytic) = {bridge_an:.6f}")

    # ---- MC: bridge + 6-product breakdown of <F^2 V_ee> ----
    cfg = dict(nw=12000, burn=3000, nsnap=16, thin=45)
    up, _ = sample_block(seed=51, **cfg); dn, _ = sample_block(seed=52, **cfg)
    Pintra = [(0, 1), (2, 3)]; INTER = [(0, 2), (0, 3), (1, 2), (1, 3)]
    br, SSc, SSi, SIc, SIi, IIc, IIi, F2V, chi2V, FV = [], [], [], [], [], [], [], [], [], []
    Fall = []
    POS = []
    for s in range(len(up)):
        pos = np.concatenate([up[s], dn[s]], axis=1); POS.append(pos)
        dd = lambda i, j: np.linalg.norm(pos[:, i] - pos[:, j], axis=-1)
        Fall.append(sum(np.exp(-GAM * dd(i, j)) for i, j in Pintra + INTER))
    Fbar = np.concatenate(Fall).mean()
    for pos in POS:
        dd = lambda i, j: np.maximum(np.linalg.norm(pos[:, i] - pos[:, j], axis=-1), 1e-12)
        f12 = np.exp(-GAM * dd(0, 1)); f34 = np.exp(-GAM * dd(2, 3))
        S = f12 + f34
        I = sum(np.exp(-GAM * dd(i, j)) for i, j in INTER)
        F = S + I
        c12 = 1.0 / dd(0, 1); c34 = 1.0 / dd(2, 3); C_S = c12 + c34
        C_I = sum(1.0 / dd(i, j) for i, j in INTER)
        Vee = C_S + C_I
        br.append((f12 * f34 / dd(0, 2)).mean())
        SSc.append((S ** 2 * C_S).mean()); SSi.append((S ** 2 * C_I).mean())
        SIc.append((2 * S * I * C_S).mean()); SIi.append((2 * S * I * C_I).mean())
        IIc.append((I ** 2 * C_S).mean()); IIi.append((I ** 2 * C_I).mean())
        F2V.append((F ** 2 * Vee).mean()); FV.append((F * Vee).mean())
        chi2V.append(((F - Fbar) ** 2 * Vee).mean())
    m = lambda x: np.mean(x); e = lambda x: np.std(x, ddof=1) / np.sqrt(len(x))
    print(f"\n  --- MC (block sampler), Fbar={Fbar:.5f} ---")
    print(f"    <f12 f34 coul_13> MC = {m(br):.6f} +/- {e(br):.1e}   (analytic {bridge_an:.6f}, "
          f"dev {abs(bridge_an-m(br))/e(br):.1f} sigma)")
    print(f"\n  6-product breakdown of <F^2 V_ee> (analytic targets):")
    for nm, x in [('<S^2 C_S>', SSc), ('<S^2 C_I>', SSi), ('2<S I C_S>', SIc),
                  ('2<S I C_I>', SIi), ('<I^2 C_S>', IIc), ('<I^2 C_I>', IIi)]:
        print(f"    {nm:12s} = {m(x):+.5f} +/- {e(x):.1e}")
    print(f"    {'<F^2 V_ee>':12s} = {m(F2V):+.5f} +/- {e(F2V):.1e}")
    print(f"    {'<F V_ee>':12s} = {m(FV):+.5f}   {'<chi^2 V_ee>':12s} = {m(chi2V):+.5f} +/- {e(chi2V):.1e}  (=g_Vee, target +0.5361)")
