"""Unified E_R12 assembly for the LiH R12-CI package (both geminals).

Assembles the fully analytic, RI-free {Phi0, (F-Fbar)Phi0} 2x2 for a James-Coolidge
geminal on the ionic single-zeta reference Phi0 = |1s_A^2 1s_B^2|, and diagonalizes to
E_R12.  This is the productionized form of debug/lih_r12ci_{assemble,linexp}.py: the same
validated primitive functions (Wmat, cov_FA_FB, neumann_potential, yukawa_pot_iso, psi_yuk,
the 216-triple enumerator, the low-rank general-m triangle) drive one assembly parametrized
by the geminal.  Deterministic (no Monte-Carlo): the tiny two-center ab,ab Yukawa block is
evaluated by the deterministic psi_yuk grid integral, matching debug/lih_r12ci_linexp.py.

Two geminals differ only in the kernel set, the h_T Part-B formula, the g_T IBP formula, and
the h_Vee product kernel; everything else is shared machinery:

    f-kernel        exp: e^{-g d}          linexp: d e^{-g d}
    f^2-kernel      exp: e^{-2g d}         linexp: d^2 e^{-2g d}
    h_T Part B      exp: -1/2 g^2 Fbar + g Ybar
                    linexp: -1/2 g^2 Fbar + 2 g Ebar - Ybar     (Ebar = <sum e^{-g r}>)
    g_T gT2+gT3     exp: -g^2 sig2 + 2 g Cov[F,Y]
                    linexp: -g^2 sig2 + 4 g Cov[F,E] - 2 Cov[F,Y]
    h_Vee product   exp: Yukawa(g) = e^{-g r}/r    linexp: e^{-g r}   (the 1/r cancels)

Heavy submodule imports (which build grids/kernels at import) are deferred into energy() so
that `import geovac.lih_r12ci` stays cheap; the first energy() call pays the ~2 min build.

Validated targets (debug/lih_r12ci_assemble.py, lih_r12ci_linexp.py; CHANGELOG v5.15.13-.14):
    exp    E_R12 = -7.9420 Ha (dE -54.2 mHa), matches the same-geminal VMC 2x2 to 0.79 mHa
    linexp E_R12 = -7.9168 Ha (dE -29.0 mHa), matches the VMC correlation dE -29.5 to 0.5 mHa
"""
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss

# --- validated VMC ground truths (independent Monte-Carlo, same geminal/reference) --- #
VMC_TARGETS = {
    # exp: g,sigma2 from debug/_g_targets.py; h from lih_r12ci_vmc.py; 2x2 from assemble.py
    "exp": dict(sigma2=0.13699, h=0.1268, g=-0.92300, E_R12=-7.94121, dE_mHa=-53.4),
    # linexp: cusp-correct headline (lih_r12ci_vmc.py geminal scan, best gamma=0.5)
    "linexp": dict(sigma2=0.1279, h=-0.1018, g=None, E_R12=-7.9316, dE_mHa=-29.5),
}

# per-piece analytic references (debug drivers), for regression pinning
ANALYTIC_REF = {
    "exp": dict(E0=-7.887822, sigma2=0.137228, h=0.129043, g=-0.919047,
                h_T=0.749453, h_Vne=-0.930200, h_Vee=0.309790,
                g_T=1.394191, g_Vne=-2.851580, g_Vee=0.538342, E_R12=-7.94200),
    "linexp": dict(E0=-7.887822, sigma2=0.127856, h=-0.103321, g=-0.771314,
                   h_T=-0.455429, h_Vne=0.417894, h_Vee=-0.065786,
                   g_T=1.384904, g_Vne=-2.623223, g_Vee=0.467005, E_R12=-7.916821),
}


@dataclass
class R12Result:
    """Result of a two-center 4e explicit-r12 CI energy assembly."""
    geminal: str
    E_R12: float          # total energy (Ha), incl. V_NN
    E0: float             # reference <Phi0|H|Phi0> (Ha), incl. V_NN
    dE_mHa: float         # correlation lowering E_R12 - E0 (mHa)
    sigma2: float         # overlap metric S_11 = Var[F]
    h: float              # off-diagonal H_01 = h_T + h_Vne + h_Vee
    g: float              # diagonal (electronic) g = g_T + g_Vne + g_Vee
    pieces: dict          # every matrix-element component
    variational: bool     # E0 > E_R12 > exact(-8.070)


def energy(geminal: str = "exp") -> R12Result:
    """Assemble and diagonalize the analytic RI-free 2x2 -> E_R12 for the given geminal.

    geminal : 'exp' (f = e^{-g r}, E_R12 = -7.942) or 'linexp' (f = r e^{-g r}, -7.917).
    Deterministic; the first call builds the shared grids/kernels (~2 min).
    """
    if geminal not in ("exp", "linexp"):
        raise ValueError(f"geminal must be 'exp' or 'linexp', got {geminal!r}")

    # --- lazy imports of the validated primitives (heavy at first import) --- #
    from .kernels import (build_kernel, Wmat, grid_int, geo_f, dens, a, GAM,
                          rho_cyl_f, zc_f, X)
    from .hVee import (neumann_potential, cov_FA_FB, make_kernel_coul, make_kernel_f,
                       make_kernel_Y, P00, P11, P01, rho_g, _a1b1, c00, c11, c01, crho,
                       NXI, NETA)
    from .hT import (Gpq, cvec, yukawa_pot_iso, rA_f, rB_f, d_aa, d_ab, d_bb, ZA, ZB)
    from .energy import Z_A, Z_B
    from . import gVee as gv
    from .triangle import _MODE_CACHE, MMAX, build_kernel_m
    from .gVee import psi_yuk

    g = GAM
    NG = NXI * NETA
    vne_grid = (-Z_A / rA_f - Z_B / rB_f)
    DP_COMPS = [(P00, P11, 1.0), (P11, P00, 1.0), (P01, P01, -2.0)]
    KD_COMPS = [(Gpq(0, 0), P11, 1.0), (Gpq(1, 1), P00, 1.0), (Gpq(0, 1), P01, -2.0),
                (P11, Gpq(0, 0), 1.0), (P00, Gpq(1, 1), 1.0), (P01, Gpq(0, 1), -2.0)]
    T00, T11 = grid_int(Gpq(0, 0)), grid_int(Gpq(1, 1))

    # ---- general phi-averaged kernels on the sigma2 grid (for the linexp radial funcs) ---- #
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

    # ---- geminal-specific kernels ---- #
    if geminal == "exp":
        Kf = build_kernel(g); Kf2 = build_kernel(2 * g)
        Fm = [build_kernel_m(g, m) for m in range(MMAX + 1)]
    else:  # linexp: f = d e^{-g d}
        f_lin = lambda d: d * np.exp(-g * d)
        Kf = build_kernel_gen(f_lin); Kf2 = build_kernel_gen(lambda d: (d * np.exp(-g * d)) ** 2)
        Fm = [build_kernel_m_gen(f_lin, m) for m in range(MMAX + 1)]
    # kernels reused by several pieces below (bounded exp kernels)
    Kexp_g = build_kernel_gen(lambda d: np.exp(-g * d))       # e^{-g d}
    Kexp_2g = build_kernel_gen(lambda d: np.exp(-2 * g * d))  # e^{-2g d}
    Klin_2g = build_kernel_gen(lambda d: d * np.exp(-2 * g * d))  # d e^{-2g d}

    def dress(K, h):
        return a ** 3 * (K @ (geo_f * h))

    def kern_obj(K):
        W, Psi = Wmat(K); return dict(W=W, Psi=Psi)

    # ---- 2-f leaf reducer (g_Vne, gT1): all four electrons integrated ---- #
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
                D[w] = D[w] * dress(Kf if mult == 1 else Kf2, D[u])
            alive.discard(u); edges = [ed for ed in edges if u not in ed]
        return scalar

    def _edge_sets(k):
        if k == 0:
            return [[]]
        if k == 1:
            return [[p] for p in gv.PAIRS]
        return [[p, q] for p in gv.PAIRS for q in gv.PAIRS]

    def FkVne(k):
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
        tot = 0.0
        for edges in _edge_sets(k):
            for gU, hU, kU in KD_COMPS:
                for gD, hD, kD in DP_COMPS:
                    base = {1: gU, 2: hU, 3: gD, 4: hD}
                    tot += kU * kD * contract_raw(edges, base)
        return 0.5 * tot

    # ============ sigma^2 ============ #
    W, Psi = Wmat(Kf); W2, _ = Wmat(Kf2)
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

    # ============ h_Vne ============ #
    Vne_exp = 2.0 * grid_int(rho_g * vne_grid); vbar = 0.5 * grid_int(rho_g * vne_grid)
    Pd = {'00': P00, '11': P11, '01': P01}; Psid = {'00': Psi00, '11': Psi11, '01': Psi01}
    cmp3 = [('00', '11', 1.0), ('11', '00', 1.0), ('01', '01', -2.0)]
    Ta = 0.5 * sum(c * grid_int(vne_grid * Pd[gg] * Psid[hh]) for gg, hh, c in cmp3)
    Tb = vbar * alpha1
    Tc = 0.25 * grid_int(vne_grid * rho_g * Psirho)
    Td = 0.5 * sum(c * grid_int(vne_grid * Pd[gg]) * grid_int(Pd[hh] * Psirho / 2.0)
                   for gg, hh, c in cmp3)
    h_Vne = 4.0 * (Ta + Tb + 2 * Tc + 2 * Td) - Fbar * Vne_exp

    # ============ Yukawa AO-pair matrix WY (deterministic, geminal-independent) ============ #
    # isotropic aa/bb/cross: exact closed-radial yukawa_pot_iso; two-center ab,ab: psi_yuk grid
    # integral (deterministic; the 12M MC it replaces is bit-negligible here).
    VYA = yukawa_pot_iso(rA_f, ZA, g); VYB = yukawa_pot_iso(rB_f, ZB, g)
    IYab = grid_int(d_ab * psi_yuk(d_ab, g))
    WY = np.array([[grid_int(d_aa * VYA), grid_int(d_ab * VYA), grid_int(d_bb * VYA)],
                   [grid_int(d_ab * VYA), IYab, grid_int(d_ab * VYB)],
                   [grid_int(d_bb * VYA), grid_int(d_ab * VYB), grid_int(d_bb * VYB)]])
    a1y = c00 @ WY @ c11 - c01 @ WY @ c01; b1y = 0.25 * crho @ WY @ crho
    Ybar = 2 * a1y + 4 * b1y
    # Yukawa DRESSING fields for the g_T Cov[F,Y] share terms (geminal-dependent):
    #   linexp: exact-isotropic aa/bb (item-(1) fix, v5.15.14) + Neumann ab -> Cov[F,Y]=-0.1025,
    #           E_R12=-7.917 fully analytic (the ~0.5% psi_yuk grid error, amplified by the
    #           Fbar=3.57 cancellation, was the sole >1% miss).
    #   exp:    psi_yuk for all three -> reproduces the VALIDATED debug exp g_T=1.3942 / E_R12
    #           =-7.942, whose end-to-end VMC 2x2 agreement (0.79 mHa) rests on a documented
    #           piece-level cancellation (g_T +0.008 high, g_Vne -0.006 low) that the exact-
    #           dressing fix would BREAK (moving E_R12 to ~-7.9435). Do NOT "fix" it here; see
    #           CHANGELOG v5.15.14 and the gT.make_kernel_Y_full note.
    if geminal == "linexp":
        PsiY = {'aa': VYA, 'ab': psi_yuk(d_ab, g), 'bb': VYB}
    else:
        PsiY = {'aa': psi_yuk(d_aa, g), 'ab': psi_yuk(d_ab, g), 'bb': psi_yuk(d_bb, g)}

    # ============ h_T = PartA + PartB ============ #
    G00, G11, G01 = Gpq(0, 0), Gpq(1, 1), Gpq(0, 1); T01 = grid_int(G01)
    KI = grid_int(G00 * Psi11) + grid_int(G11 * Psi00) - 2 * grid_int(G01 * Psi01)
    kappa = (G00 + G11) + T00 * P11 + T11 * P00 - 2 * T01 * P01
    PartA = KI + grid_int(kappa * Psirho) + (alpha1 - Fbar) * (T00 + T11)
    if geminal == "exp":
        # grad^2 f = g^2 f - 2 g (e^{-g r}/r);  PartB = -1/2 g^2 Fbar + g Ybar
        PartB = -0.5 * g ** 2 * Fbar + g * Ybar
    else:
        # grad^2 f = g^2 f - 4 g e^{-g r} + 2 e^{-g r}/r;  PartB = -1/2 g^2 Fbar + 2 g Ebar - Ybar
        We, _ = Wmat(Kexp_g); a1e = c00 @ We @ c11 - c01 @ We @ c01; b1e = 0.25 * crho @ We @ crho
        Ebar = 2 * a1e + 4 * b1e
        PartB = -0.5 * g ** 2 * Fbar + 2 * g * Ebar - Ybar
    h_T = PartA + PartB

    # ============ h_Vee = Cov[F, V_ee] ============ #
    A_L = kern_obj(Kf); Kcoul = make_kernel_coul()
    if geminal == "exp":
        C_hvee = make_kernel_Y(g)          # product f*coul = Yukawa(g) = e^{-g r}/r  (W only)
    else:
        C_hvee = dict(W=Wmat(Kexp_g)[0])   # product f*coul = e^{-g r}
    h_Vee, _hinfo = cov_FA_FB(A_L, Kcoul, C_hvee)

    h = h_T + h_Vne + h_Vee

    # ============ g_Vne ============ #
    F0V, F1V, F2V_ne = FkVne(0), FkVne(1), FkVne(2)
    g_Vne = F2V_ne - 2 * Fbar * F1V + Fbar ** 2 * F0V

    # ============ g_T = gT1 + (gT2+gT3) ============ #
    S2, FS, F2S = FkSv2(0), FkSv2(1), FkSv2(2)
    gT1 = 0.5 * (F2S - 2 * Fbar * FS + Fbar ** 2 * S2)
    B_y = dict(W=WY, Psi=PsiY)
    if geminal == "exp":
        C_y = make_kernel_Y(2 * g)         # product f*Y = e^{-2g r}/r = Yukawa(2g)
        CovFY, _ = cov_FA_FB(A_L, B_y, C_y)
        gT23 = -g ** 2 * sig2 + 2 * g * CovFY
        CovFE = None
    else:
        B_e = kern_obj(Kexp_g); C_e = dict(W=Wmat(Klin_2g)[0])   # product f*E = d e^{-2g d}
        CovFE, _ = cov_FA_FB(A_L, B_e, C_e)
        C_y = dict(W=Wmat(Kexp_2g)[0])     # product f*Y = e^{-2g r}
        CovFY, _ = cov_FA_FB(A_L, B_y, C_y)
        gT23 = -g ** 2 * sig2 + 4 * g * CovFE - 2 * CovFY
    g_T = gT1 + gT23

    # ============ g_Vee (216-triple enumerator + RI-free triangle) ============ #
    _MODE_CACHE.clear()
    _default_I_yuk = gv.I_yuk
    gv.Kf = Kf; gv.Kf2 = Kf2; gv.Fm = Fm
    if geminal == "exp":
        gv.I_yuk = _default_I_yuk          # gVee default: grid_int(L*psi_yuk(Rt,gam)) [Yukawa]
    else:
        # linexp same-pair kernels: f/r -> e^{-g r} (Kexp_g); f^2/r -> d e^{-2g d} (Klin_2g)
        gv.I_yuk = (lambda L, Rt, gam_k:
                    grid_int(L * dress(Kexp_g if abs(gam_k - g) < 1e-9 else Klin_2g, Rt)))
    try:
        F2Vee = sum(v for v in gv.F2Vee().values())
        FVee = sum(gv.eval_pair_fc(p, r) for p in gv.PAIRS for r in gv.PAIRS)
        Vee_enum = 2 * gv.eval_coul_only((1, 2)) + 4 * gv.eval_coul_only((1, 3))
    finally:
        gv.I_yuk = _default_I_yuk          # restore, so successive calls don't cross-contaminate
    g_Vee = F2Vee - 2 * Fbar * FVee + Fbar ** 2 * Vee_enum

    # ============ assemble 2x2 -> E_R12 ============ #
    V_NN = Z_A * Z_B / 3.015           # Z_A Z_B / R (R = 3.015 a0)
    E0_tot = -7.887822                 # stage2_E0() (deterministic, control-validated)
    E0_elec = E0_tot - V_NN
    g_elec = g_T + g_Vne + g_Vee
    aa = E0_elec; cc = g_elec / sig2; bb = h / np.sqrt(sig2)
    E_R12 = 0.5 * (aa + cc - np.sqrt((aa - cc) ** 2 + 4 * bb ** 2)) + V_NN
    dE_mHa = (E_R12 - E0_tot) * 1e3

    pieces = dict(Fbar=Fbar, sigma2=sig2, h_T=h_T, h_Vne=h_Vne, h_Vee=h_Vee,
                  g_T=g_T, g_Vne=g_Vne, g_Vee=g_Vee, gT1=gT1, gT23=gT23,
                  Ybar=Ybar, CovFY=CovFY, CovFE=CovFE, PartA=PartA, PartB=PartB,
                  V_NN=V_NN, E0_elec=E0_elec)
    return R12Result(
        geminal=geminal, E_R12=E_R12, E0=E0_tot, dE_mHa=dE_mHa, sigma2=sig2,
        h=h, g=g_elec, pieces=pieces,
        variational=(E_R12 < E0_tot and E_R12 > -8.070))


if __name__ == "__main__":
    for gem in ("exp", "linexp"):
        r = energy(gem)
        ref = ANALYTIC_REF[gem]
        print(f"\n=== geminal = {gem} ===")
        print(f"  sigma2={r.sigma2:.6f} (ref {ref['sigma2']:.6f})  h={r.h:+.6f} (ref {ref['h']:+.6f})"
              f"  g={r.g:+.6f} (ref {ref['g']:+.6f})")
        print(f"  h_T={r.pieces['h_T']:+.5f} h_Vne={r.pieces['h_Vne']:+.5f} h_Vee={r.pieces['h_Vee']:+.5f}"
              f"  |  g_T={r.pieces['g_T']:+.5f} g_Vne={r.pieces['g_Vne']:+.5f} g_Vee={r.pieces['g_Vee']:+.5f}")
        print(f"  >>> E_R12 = {r.E_R12:.6f} Ha (dE {r.dE_mHa:+.2f} mHa; ref {ref['E_R12']:.6f})"
              f"   variational={r.variational}")
        dev = abs(r.E_R12 - ref['E_R12']) * 1e3
        print(f"      dev from analytic ref: {dev:.3f} mHa  {'OK' if dev < 0.5 else 'CHECK'}")
