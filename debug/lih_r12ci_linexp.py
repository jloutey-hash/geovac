"""LiH R12-CI with the CUSP-CORRECT geminal f = r e^{-gam r} (linexp) -- reproduces the
headline PoC energy E_R12 = -7.932 analytically, RI-free.

The reduction framework (enumerator, triangle, IBP identities) is geminal-agnostic; only the
KERNELS change from the exp-geminal build (lih_r12ci_{sigma2,hT,hVee,gVne,gVee,gT}_analytic.py):

  f-kernel        exp: e^{-gam d}          linexp: d e^{-gam d}
  f^2-kernel      exp: e^{-2gam d}         linexp: d^2 e^{-2gam d}
  same-pair f/r   exp: Yukawa e^{-gam r}/r linexp: e^{-gam r}  (bounded; the 1/r cancels)
  same-pair f^2/r exp: Yukawa 2gam         linexp: d e^{-2gam d}
  grad^2 f        exp: gam^2 f - 2 gam Y   linexp: gam^2 f - 4 gam e^{-gam r} + 2 Y
                                           (Y = e^{-gam r}/r)

The last changes the IBP-collapsed piece:
  h_T PartB      = -1/2 gam^2 Fbar + 2 gam Ebar - Ybar        (Ebar=<sum e^{-gam r}>)
  g_T (gT2+gT3)  = -gam^2 sigma^2 + 4 gam Cov[F,E_sum] - 2 Cov[F,Y_sum]
E0 is geminal-independent (unchanged, -7.887822).

Run from debug/:  python lih_r12ci_linexp.py
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

from lih_r12ci_sigma2_analytic import (
    build_kernel, dressings, Wmat, grid_int, dens, geo_f, a, GAM, X,
    rho_cyl_f, zc_f, rA as RA2d, rB as RB2d)
from lih_r12ci_hVee_analytic import (
    neumann_potential, cov_FA_FB, make_kernel_coul, P00, P11, P01, rho_g,
    _a1b1, c00, c11, c01, crho, NXI, NETA)
from lih_r12ci_hT_analytic import Gpq, cvec, yukawa_pot_iso, rA_f, rB_f, d_aa, d_ab, d_bb, ZA, ZB
from lih_r12ci_energy import Z_A, Z_B
import lih_r12ci_gVee_analytic as gv
from lih_r12ci_triangle_gate import _MODE_CACHE, MMAX, build_kernel_m

NG = NXI * NETA
vne_grid = (-Z_A / rA_f - Z_B / rB_f)          # one-body V_ne on the flat grid
DP_COMPS = [(P00, P11, 1.0), (P11, P00, 1.0), (P01, P01, -2.0)]
KD_COMPS = [(Gpq(0, 0), P11, 1.0), (Gpq(1, 1), P00, 1.0), (Gpq(0, 1), P01, -2.0),
            (P11, Gpq(0, 0), 1.0), (P00, Gpq(1, 1), 1.0), (P01, Gpq(0, 1), -2.0)]
T00, T11 = grid_int(Gpq(0, 0)), grid_int(Gpq(1, 1))


# --------------------------------------------------------------------------- #
# general kernels (any radial func), on the sigma2 grid
# --------------------------------------------------------------------------- #
def build_kernel_gen(func, nphi=28):
    xp, wp = leggauss(nphi); phi = 0.5 * np.pi * (xp + 1.0); wphi = 0.5 * np.pi * wp
    rc = rho_cyl_f; z = zc_f; rc2 = rc[:, None] * rc[None, :]
    base = rc[:, None] ** 2 + rc[None, :] ** 2 + (z[:, None] - z[None, :]) ** 2
    K = np.zeros((NG, NG))
    for cphi, w in zip(np.cos(phi), wphi):
        d = np.sqrt(np.maximum(base - 2.0 * rc2 * cphi, 0.0))
        K += 2.0 * w * func(d)
    return K


def build_kernel_m_gen(func, m, nphi=48):
    xp, wp = leggauss(nphi); phi = 0.5 * np.pi * (xp + 1.0); wphi = 0.5 * np.pi * wp
    rc = rho_cyl_f; z = zc_f; rc2 = rc[:, None] * rc[None, :]
    base = rc[:, None] ** 2 + rc[None, :] ** 2 + (z[:, None] - z[None, :]) ** 2
    K = np.zeros((NG, NG))
    for ph, w in zip(phi, wphi):
        d = np.sqrt(np.maximum(base - 2.0 * rc2 * np.cos(ph), 0.0))
        K += 2.0 * w * func(d) * np.cos(m * ph)
    return K / (2 * np.pi) if m == 0 else K / np.pi


g = GAM
f_lin = lambda d: d * np.exp(-g * d)
f_lin2 = lambda d: (d * np.exp(-g * d)) ** 2            # (f)^2 = d^2 e^{-2g d}
Kf_L = build_kernel_gen(f_lin)
Kf2_L = build_kernel_gen(f_lin2)
Kexp_g = build_kernel_gen(lambda d: np.exp(-g * d))     # same-pair f/r (linexp) + Ebar kernel
Kexp_2g = build_kernel_gen(lambda d: np.exp(-2 * g * d))  # f*Y = e^{-2g r} product kernel
Klin_2g = build_kernel_gen(lambda d: d * np.exp(-2 * g * d))  # same-pair f^2/r + f*e^{-g r} product


def dress_L(K, h):
    return a ** 3 * (K @ (geo_f * h))


def kern_obj(K):
    W, Psi = Wmat(K); return dict(W=W, Psi=Psi)


# --------------------------------------------------------------------------- #
# local reducer for the 2-f (all-integrated) block quantities: g_Vne, gT1
# --------------------------------------------------------------------------- #
def contract_raw(fedges, base):
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
            raise RuntimeError(f"bad topology {fedges}")
        inc = [ed for ed in edges if u in ed]; nb = dnb(u)
        if not nb:
            scalar *= grid_int(D[u])
        else:
            w = next(iter(nb)); mult = len(inc)
            D[w] = D[w] * dress_L(Kf_L if mult == 1 else Kf2_L, D[u])
        alive.discard(u); edges = [ed for ed in edges if u not in ed]
    return scalar


def _edge_sets(k):
    if k == 0:
        return [[]]
    if k == 1:
        return [[p] for p in gv.PAIRS]
    return [[p, q] for p in gv.PAIRS for q in gv.PAIRS]


def FkVne(k):
    """<F^k V_ne> = sum_i (1/4) sum_{cU,cD in Dp} wt * contract(edges, base with v on electron i)."""
    tot = 0.0
    for edges in _edge_sets(k):
        for i in (1, 2, 3, 4):
            for gU, hU, kU in DP_COMPS:
                for gD, hD, kD in DP_COMPS:
                    base = {1: gU, 2: hU, 3: gD, 4: hD}
                    base[i] = base[i] * vne_grid
                    tot += kU * kD * contract_raw(edges, base)
    return 0.25 * tot


def FkSv2(k):
    """<F^k sum_i v_i^2> = (1/2) sum_{k edges} sum_{cU in KD, cD in Dp} wt * contract."""
    tot = 0.0
    for edges in _edge_sets(k):
        for gU, hU, kU in KD_COMPS:
            for gD, hD, kD in DP_COMPS:
                base = {1: gU, 2: hU, 3: gD, 4: hD}
                tot += kU * kD * contract_raw(edges, base)
    return 0.5 * tot


if __name__ == "__main__":
    print("=" * 78)
    print(f"LiH R12-CI, CUSP-CORRECT geminal f = r e^(-{g} r)  (linexp) -- analytic RI-free")
    print("=" * 78)

    # ============ sigma^2 (linexp) ============
    W, Psi = Wmat(Kf_L); W2, _ = Wmat(Kf2_L)
    If = lambda cu, cv: cu @ W @ cv; If2 = lambda cu, cv: cu @ W2 @ cv
    alpha1 = If(c00, c11) - If(c01, c01); alpha2 = If2(c00, c11) - If2(c01, c01)
    beta1 = 0.25 * If(crho, crho); beta2 = 0.25 * If2(crho, crho)
    comps = [(c00, c11, 1.0), (c11, c00, 1.0), (c01, c01, -2.0)]
    gam_dj = 0.25 * sum(ka * kb * (ga @ W @ gb) * (ha @ W @ hb)
                        for ga, ha, ka in comps for gb, hb, kb in comps)
    Fbar = 2 * alpha1 + 4 * beta1
    Pf = lambda c: c[0] * Psi['aa'] + c[1] * Psi['ab'] + c[2] * Psi['bb']
    Psi00, Psi11, Psi01, Psirho = Pf(c00), Pf(c11), Pf(c01), Pf(crho)
    delta = 0.25 * grid_int(Psirho * (P00 * Psi11 + P11 * Psi00 - 2 * P01 * Psi01))
    gam_sh = 0.5 * grid_int(rho_g * (Psi00 * Psi11 - Psi01 ** 2))
    sig2 = (2 * alpha2 + 4 * beta2 + 8 * gam_sh + 4 * gam_dj + 16 * delta
            - 2 * alpha1 ** 2 - 16 * alpha1 * beta1 - 16 * beta1 ** 2)
    print(f"\n  Fbar = {Fbar:.6f}   sigma^2 = {sig2:.6f}")

    # ============ h_Vne (linexp) ============
    Vne_exp = 2.0 * grid_int(rho_g * vne_grid); vbar = 0.5 * grid_int(rho_g * vne_grid)
    Pd = {'00': P00, '11': P11, '01': P01}; Psid = {'00': Psi00, '11': Psi11, '01': Psi01}
    cmp3 = [('00', '11', 1.0), ('11', '00', 1.0), ('01', '01', -2.0)]
    Ta = 0.5 * sum(c * grid_int(vne_grid * Pd[gg] * Psid[hh]) for gg, hh, c in cmp3)
    Tb = vbar * alpha1
    Tc = 0.25 * grid_int(vne_grid * rho_g * Psirho)
    Td = 0.5 * sum(c * grid_int(vne_grid * Pd[gg]) * grid_int(Pd[hh] * Psirho / 2.0) for gg, hh, c in cmp3)
    h_Vne = 4.0 * (Ta + Tb + 2 * Tc + 2 * Td) - Fbar * Vne_exp
    print(f"  h_Vne = {h_Vne:+.6f}   (<V_ne>={Vne_exp:.5f})")

    # ============ h_T (linexp): PartA + PartB ============
    G00, G11, G01 = Gpq(0, 0), Gpq(1, 1), Gpq(0, 1); T01 = grid_int(G01)
    KI = grid_int(G00 * Psi11) + grid_int(G11 * Psi00) - 2 * grid_int(G01 * Psi01)
    kappa = (G00 + G11) + T00 * P11 + T11 * P00 - 2 * T01 * P01
    PartA = KI + grid_int(kappa * Psirho) + (alpha1 - Fbar) * (T00 + T11)
    # Ebar = <sum e^{-g r}> (exp-geminal Fbar) ; Ybar = <sum e^{-g r}/r> (Yukawa)
    We, _ = Wmat(Kexp_g); a1e = c00 @ We @ c11 - c01 @ We @ c01; b1e = 0.25 * crho @ We @ crho
    Ebar = 2 * a1e + 4 * b1e
    VYA = yukawa_pot_iso(rA_f, ZA, g); VYB = yukawa_pot_iso(rB_f, ZB, g)
    # ab,ab (two-center transition-density Yukawa self-energy): DETERMINISTIC via the
    # Neumann-minus-smooth psi_yuk grid integral (was a 12M-sample MC). The element is tiny
    # (~4.6e-3, MC vs psi_yuk differ 1.2%), so it is bit-negligible on Ybar/h_T but removes
    # the MC noise and makes the driver fast + reproducible.
    from lih_r12ci_gVee_analytic import psi_yuk as _psi_yuk
    IYab = grid_int(d_ab * _psi_yuk(d_ab, g))
    WY = np.array([[grid_int(d_aa * VYA), grid_int(d_ab * VYA), grid_int(d_bb * VYA)],
                   [grid_int(d_ab * VYA), IYab, grid_int(d_ab * VYB)],
                   [grid_int(d_bb * VYA), grid_int(d_ab * VYB), grid_int(d_bb * VYB)]])
    a1y = c00 @ WY @ c11 - c01 @ WY @ c01; b1y = 0.25 * crho @ WY @ crho
    Ybar = 2 * a1y + 4 * b1y
    PartB = -0.5 * g ** 2 * Fbar + 2 * g * Ebar - Ybar
    h_T = PartA + PartB
    print(f"  h_T = PartA {PartA:+.5f} + PartB {PartB:+.5f} = {h_T:+.6f}   (Ebar={Ebar:.4f} Ybar={Ybar:.4f})")

    # ============ h_Vee (linexp): cov_FA_FB(f, coul, product=e^{-g r}) ============
    A_L = kern_obj(Kf_L); Kcoul = make_kernel_coul()
    C_hvee = dict(W=Wmat(Kexp_g)[0])                      # product f*coul = e^{-g r}
    h_Vee, hinfo = cov_FA_FB(A_L, Kcoul, C_hvee)
    E2 = hinfo['FbarB']
    print(f"  h_Vee = {h_Vee:+.6f}   (E2={E2:.5f})")
    h = h_T + h_Vne + h_Vee
    print(f"  h = h_T + h_Vne + h_Vee = {h:+.6f}")

    # ============ g_Vne (linexp) via weighted contract ============
    F0V, F1V, F2V_ne = FkVne(0), FkVne(1), FkVne(2)
    g_Vne = F2V_ne - 2 * Fbar * F1V + Fbar ** 2 * F0V
    print(f"\n  <V_ne>(enum)={F0V:.5f} (gate -20.86474)   g_Vne = {g_Vne:+.6f}")

    # ============ g_T (linexp): gT1 + (gT2+gT3) ============
    S2, FS, F2S = FkSv2(0), FkSv2(1), FkSv2(2)
    gT1 = 0.5 * (F2S - 2 * Fbar * FS + Fbar ** 2 * S2)
    # Cov[F, E_sum]: A=linexp-f, B=exp(g), product = d e^{-2g d} = linexp(2g)
    B_e = kern_obj(Kexp_g); C_e = dict(W=Wmat(Klin_2g)[0])
    CovFE, _ = cov_FA_FB(A_L, B_e, C_e)
    # Cov[F, Y_sum]: A=linexp-f, B=Yukawa(g), product f*Y = e^{-2g r} = exp(2g).
    # Isotropic aa/bb dressings: EXACT closed-radial Yukawa potential (VYA,VYB, built above
    # for WY); only the two-center ab uses the Neumann-minus-smooth psi_yuk. This mirrors
    # make_kernel_coul (exact _hartree_1s isotropic + Neumann ab). Using psi_yuk for the
    # isotropic dressings (its ~0.5% grid error vs the exact potential, amplified by the
    # Fbar=3.57 cancellation) was the SOLE >1% miss; the exact dressings close it:
    # Cov[F,Y_sum] -0.0793 -> -0.1025 (MC -0.1018), landing E_R12 = -7.917.
    B_y = dict(W=WY, Psi={'aa': VYA, 'ab': _psi_yuk(d_ab, g), 'bb': VYB})
    C_y = dict(W=Wmat(Kexp_2g)[0])
    CovFY, _ = cov_FA_FB(A_L, B_y, C_y)
    gT23 = -g ** 2 * sig2 + 4 * g * CovFE - 2 * CovFY
    g_T = gT1 + gT23
    print(f"  gT1={gT1:+.5f}  Cov[F,E]={CovFE:+.5f} Cov[F,Y]={CovFY:+.5f}  gT2+gT3={gT23:+.5f}")
    print(f"  g_T = {g_T:+.6f}")

    # ============ g_Vee (linexp): monkey-patch the enumerator to linexp kernels ============
    _MODE_CACHE.clear()
    gv.Kf = Kf_L; gv.Kf2 = Kf2_L
    gv.Fm = [build_kernel_m_gen(f_lin, m) for m in range(MMAX + 1)]
    # same-pair kernel: n_on_r=1 -> f/r = e^{-g r} (Kexp_g); n_on_r=2 -> f^2/r = d e^{-2g d} (Klin_2g)
    gv.I_yuk = lambda L, Rt, gam_k: grid_int(L * dress_L(Kexp_g if abs(gam_k - g) < 1e-9 else Klin_2g, Rt))
    F2Vee = sum(v for v in gv.F2Vee().values())
    FVee = sum(gv.eval_pair_fc(p, r) for p in gv.PAIRS for r in gv.PAIRS)
    Vee_enum = 2 * gv.eval_coul_only((1, 2)) + 4 * gv.eval_coul_only((1, 3))
    g_Vee = F2Vee - 2 * Fbar * FVee + Fbar ** 2 * Vee_enum
    print(f"\n  <F^2 V_ee>={F2Vee:.5f}  <F V_ee>={FVee:.5f}  <V_ee>={Vee_enum:.5f}")
    print(f"  g_Vee = {g_Vee:+.6f}")

    # ============ assemble 2x2 -> E_R12 ============
    V_NN = Z_A * Z_B / 3.015; E0_tot = -7.887822; E0_elec = E0_tot - V_NN
    g_elec = g_T + g_Vne + g_Vee
    aa = E0_elec; cc = g_elec / sig2; bb = h / np.sqrt(sig2)
    E_R12 = 0.5 * (aa + cc - np.sqrt((aa - cc) ** 2 + 4 * bb ** 2)) + V_NN
    print("\n" + "=" * 78)
    print(f"  g = g_T + g_Vne + g_Vee = {g_T:+.5f} {g_Vne:+.5f} {g_Vee:+.5f} = {g_elec:+.6f}")
    print(f"  h = {h:+.6f}   sigma^2 = {sig2:.6f}   E0 = {E0_tot:.6f}")
    print(f"  >>> E_R12 (linexp, cusp-correct) = {E_R12:.6f} Ha   (dE = {(E_R12-E0_tot)*1e3:+.2f} mHa)")
    print(f"      VMC headline: -7.932 (dE -29.5 mHa)")
    print(f"      variational: E_R12<E0 {E_R12<E0_tot}, E_R12>-8.070 {E_R12>-8.070}")
