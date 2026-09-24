"""LiH R12-CI Stage 3 (part 2) -- ANALYTIC (RI-free) sigma^2 = <Phi0|F^2|Phi0> - Fbar^2.

The structural showpiece: sigma^2 with NO Monte-Carlo and NO resolution-of-identity, to be
compared against the VMC ground truth 0.137 (lih_r12ci_sigma2_mc.py).

Reduction (block density is SEPARABLE):
  D_p(ra,rb)^2 = P00(ra)P11(rb) + P11(ra)P00(rb) - 2 P01(ra)P01(rb),   P_pq = m_p m_q,
so every term of <F^2> collapses to grid integrals of ONE-electron densities P_pq and their
f-DRESSINGS  Psi^f_h(r) = INT h(r') f(|r-r'|) d3r'.  With orthonormal MOs: marginal rho=m0^2+m1^2,
block norm N=2 (exact).  F = sum_{i<j} f(r_ij), 6 pairs; f = exp(-0.5 r) (Stage-1 geminal).

Term inventory (A=f_12 up-intra, B=f_34 dn-intra, C=inter):
  alpha1=<f>_intra, alpha2=<f^2>_intra ; beta1,beta2 = inter rho-rho (f, f^2) moments
  delta  = <f_12 f_13>   (intra x inter, share-vertex 3-body)
  gam_sh = <f_13 f_14>   (inter x inter, share-vertex)
  gam_dj = <f_13 f_24>   (inter x inter, disjoint but block-correlated)
  sigma^2 = 2 alpha2 + 4 beta2 + 8 gam_sh + 4 gam_dj + 16 delta
            - 2 alpha1^2 - 16 alpha1 beta1 - 16 beta1^2
  (verified: f=const -> sigma^2 = 0, coefficients 2+4+8+4+16-2-16-16 cancel.)

alpha/beta/gam_dj = SCALAR f-interactions I_f[P_pq,P_rs] = c_pq . W . c_rs (W = 3x3 AO-pair
f-interaction matrix, = Stage-1 f-tensor).  delta, gam_sh need the DRESSING FIELDS Psi^f_{P_pq}
= combos of Psi^f_{aa,ab,bb}; the two-center Psi^f_{ab} is the one new piece.  P_pq in AO-pair
basis {aa,ab,bb}: c^pq = [X0p X0q, X0p X1q + X1p X0q, X1p X1q].

Run:  python debug/lih_r12ci_sigma2_analytic.py
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

from .energy import ZA, ZB, N_A, N_B, R, a, Z_A, Z_B
# (the Stage-1 f-integrals V_aaaa..V_bbab are imported lazily in __main__ only: energy.py now
#  builds them on first access (~60 s) and nothing else in the package consumes them -- 2026-09-22)
from .basis import X          # Loewdin AO->MO (computed from S_AB)

SIG2_REF = 0.13699          # VMC ground truth (lih_r12ci_sigma2_mc.py)
FBAR_REF = 1.880743         # Stage 1 analytic Fbar
VABAB_REF = 0.004350        # Stage 1 (ab|f|ab) importance-MC (validates the 2-center dressing)
GAM = 0.5                   # geminal f = exp(-GAM r)

# --------------------------------------------------------------------------- #
# 2D prolate (xi,eta) grid  (all densities are phi-independent -> axial symmetry)
# --------------------------------------------------------------------------- #
NXI, NETA, NPHI = 72, 44, 28
xi_max = 1.0 + 42.0 / (2 * min(ZA, ZB) * a)
_xg, _wxg = leggauss(NXI); XI = 1.0 + 0.5 * (_xg + 1.0) * (xi_max - 1.0)
WXI = 0.5 * (xi_max - 1.0) * _wxg
_eg, _weg = leggauss(NETA); ETA = _eg.copy(); WETA = _weg.copy()
Xg, Eg = np.meshgrid(XI, ETA, indexing='ij')                       # (NXI,NETA)
rA = a * (Xg + Eg); rB = a * (Xg - Eg)
JAC = Xg ** 2 - Eg ** 2                                             # d3r = a^3 JAC dxi deta dphi
RHO_CYL = a * np.sqrt(np.maximum((Xg ** 2 - 1.0) * (1.0 - Eg ** 2), 0.0))
ZC = a * Xg * Eg
W2D = np.outer(WXI, WETA)                                           # (NXI,NETA) 2D weights
NG = Xg.size
_flat = lambda M: M.reshape(-1)
rho_cyl_f = _flat(RHO_CYL); zc_f = _flat(ZC)
geo_f = _flat(W2D * JAC)                                            # 2D quadrature weight x Jacobian

# AO-pair densities (amplitudes: orb = N exp(-Z r)); ρ_aa=orb_A^2, ρ_ab=orb_A orb_B, ρ_bb=orb_B^2
orbA = N_A * np.exp(-ZA * rA); orbB = N_B * np.exp(-ZB * rB)
dens = {'aa': _flat(orbA ** 2), 'ab': _flat(orbA * orbB), 'bb': _flat(orbB ** 2)}
AOP = ['aa', 'ab', 'bb']


def build_kernel(gam):
    """K[i,j] = INT_0^2pi exp(-gam |r_i - r_j(phi')|) dphi'  (axially-averaged f-kernel)."""
    xp, wp = leggauss(NPHI)                       # phi' on [0,pi], x2 for [0,2pi] (cos symmetric)
    phi = 0.5 * np.pi * (xp + 1.0); wphi = 0.5 * np.pi * wp
    rc = rho_cyl_f; z = zc_f
    K = np.zeros((NG, NG))
    rc2 = rc[:, None] * rc[None, :]
    base = rc[:, None] ** 2 + rc[None, :] ** 2 + (z[:, None] - z[None, :]) ** 2
    for cphi, w in zip(np.cos(phi), wphi):
        d = np.sqrt(np.maximum(base - 2.0 * rc2 * cphi, 0.0))
        K += 2.0 * w * np.exp(-gam * d)
    return K


def dressings(K):
    """Psi^f_h(r_i) = a^3 sum_j geo_j h_j K[i,j]  for h in {aa,ab,bb}."""
    src = {k: (geo_f * dens[k]) for k in AOP}
    return {k: a ** 3 * (K @ src[k]) for k in AOP}


def Wmat(K):
    """3x3 AO-pair f-interaction matrix  I_f[rho_i,rho_j] = 2pi a^3 sum_k geo_k rho_i,k Psi_j,k."""
    Psi = dressings(K)
    W = np.zeros((3, 3))
    for i, u in enumerate(AOP):
        for j, v in enumerate(AOP):
            W[i, j] = 2 * np.pi * a ** 3 * np.sum(geo_f * dens[u] * Psi[v])
    return W, Psi


def grid_int(field):
    """INT field d3r = 2pi a^3 sum geo * field   (field phi-independent)."""
    return 2 * np.pi * a ** 3 * np.sum(geo_f * field)


# --------------------------------------------------------------------------- #
# Phase 0b (2026-09-22): EXACT ordered-integral prolate-Neumann operator on THIS grid, opt-in.
#   The legacy radial step  K @ (WXI * g)  in hVee.neumann_potential / triangle.coul_mode_potential
#   integrates the kinked kernel P_l^m(xi_<)Q_l^m(xi_>) by a cumulative GL sum over the 72 nodes:
#   O(NXI^-2), +5e-3 relative on the Li-1s self-Coulomb, +8.09 mHa on the Phase-0 <V_ee>.  With
#   USE_EXACT_NEUMANN = True the same consumers (and gVee.psi_coul / psi_yuk through them) use
#   neumann_exact.ExactNeumann instead -- exact relative to the degree-71 interpolant through the
#   GL nodes (5 zeta/8 to 1e-13 on this grid; see the module docstring).  Default False keeps the
#   legacy path BIT-IDENTICAL (tests/test_lih_r12ci.py anchors -7.9168 / -7.9420 on it).
#   Consumers read the flag AT CALL TIME as a module attribute (never by from-import).
# --------------------------------------------------------------------------- #
USE_EXACT_NEUMANN = False
_EXACT_MMAX = 4                     # triangle.MMAX; one operator serves the m=0 and general-m callers
_EXACT_CACHE = {}


def exact_neumann(lmax: int, mmax: int = 0):
    """Cached neumann_exact.ExactNeumann on the kernels.py xi grid, covering (lmax, mmax)."""
    for (lm, mm), op in _EXACT_CACHE.items():
        if lm >= lmax and mm >= mmax:
            return op
    from .neumann_exact import ExactNeumann
    mm = max(mmax, _EXACT_MMAX)
    op = ExactNeumann(XI, 1.0, xi_max, lmax, mm)
    _EXACT_CACHE[(lmax, mm)] = op
    return op


if __name__ == "__main__":
    from .energy import V_aaaa, V_bbbb, V_aabb, V_aaab, V_bbab    # lazy Stage-1 build (~60 s)
    print("=" * 78)
    print("LiH R12-CI Stage 3 (part 2): ANALYTIC sigma^2 (RI-free), target 0.137")
    print(f"  grid {NXI}x{NETA} (xi_max={xi_max:.1f}), NPHI={NPHI}; f=exp(-{GAM} r)")
    print("=" * 78)

    Kf = build_kernel(GAM); Kf2 = build_kernel(2 * GAM)             # f and f^2 kernels
    W, Psi = Wmat(Kf)                                              # f-interaction matrix + fields
    W2, _ = Wmat(Kf2)                                             # f^2-interaction matrix

    # ---- GATE 1: W vs Stage-1 f-tensor (validates the dressing machinery, incl. 2-center ab) ----
    idx = {'aa': 0, 'ab': 1, 'bb': 2}
    checks = [('I_f[aa,aa]', W[0, 0], V_aaaa), ('I_f[bb,bb]', W[2, 2], V_bbbb),
              ('I_f[aa,bb]', W[0, 2], V_aabb), ('I_f[aa,ab]', W[0, 1], V_aaab),
              ('I_f[bb,ab]', W[2, 1], V_bbab), ('I_f[ab,ab]', W[1, 1], VABAB_REF)]
    print("\n  GATE 1 -- W (analytic dressing) vs Stage-1 f-tensor:")
    for nm, got, ref in checks:
        print(f"    {nm} = {got:.6f}  vs {ref:.6f}   (rel {abs(got-ref)/abs(ref):.1e})")

    # ---- MO-density coefficient vectors in the AO-pair basis {aa,ab,bb} ----
    def cvec(p, q):
        return np.array([X[0, p] * X[0, q],
                         X[0, p] * X[1, q] + X[1, p] * X[0, q],
                         X[1, p] * X[1, q]])
    c00, c11, c01 = cvec(0, 0), cvec(1, 1), cvec(0, 1)
    crho = c00 + c11                                               # rho = m0^2 + m1^2

    If = lambda cu, cv: cu @ W @ cv                                # scalar f-interaction
    If2 = lambda cu, cv: cu @ W2 @ cv

    # ---- 2-body / disjoint moments (scalars) ----
    alpha1 = If(c00, c11) - If(c01, c01)
    alpha2 = If2(c00, c11) - If2(c01, c01)
    beta1 = 0.25 * If(crho, crho)
    beta2 = 0.25 * If2(crho, crho)
    # gam_dj = 1/4 sum_ab kappa_a kappa_b (g_a.W.g_b)(h_a.W.h_b); comps (g,h,kappa):
    comps = [(c00, c11, 1.0), (c11, c00, 1.0), (c01, c01, -2.0)]
    gam_dj = 0.0
    for ga, ha, ka in comps:
        for gb, hb, kb in comps:
            gam_dj += ka * kb * (ga @ W @ gb) * (ha @ W @ hb)
    gam_dj *= 0.25

    # ---- GATE 2: reproduce Fbar = 2 alpha1 + 4 beta1 ----
    Fbar = 2 * alpha1 + 4 * beta1
    print(f"\n  GATE 2 -- Fbar = 2 alpha1 + 4 beta1 = {Fbar:.6f}  vs Stage-1 {FBAR_REF:.6f}"
          f"   (rel {abs(Fbar-FBAR_REF)/FBAR_REF:.1e})")

    # ---- 3-body moments (need dressing FIELDS) ----
    Psi_field = lambda c: c[0] * Psi['aa'] + c[1] * Psi['ab'] + c[2] * Psi['bb']   # Psi^f_{P}
    P00 = c00[0] * dens['aa'] + c00[1] * dens['ab'] + c00[2] * dens['bb']
    P11 = c11[0] * dens['aa'] + c11[1] * dens['ab'] + c11[2] * dens['bb']
    P01 = c01[0] * dens['aa'] + c01[1] * dens['ab'] + c01[2] * dens['bb']
    rho_f = P00 + P11
    Psi00, Psi11, Psi01, Psirho = Psi_field(c00), Psi_field(c11), Psi_field(c01), Psi_field(crho)

    # delta = 1/4 INT Psi^f_rho * [P00 Psi11 + P11 Psi00 - 2 P01 Psi01]
    delta = 0.25 * grid_int(Psirho * (P00 * Psi11 + P11 * Psi00 - 2 * P01 * Psi01))
    # gam_sh = 1/2 INT rho * [Psi00 Psi11 - Psi01^2]
    gam_sh = 0.5 * grid_int(rho_f * (Psi00 * Psi11 - Psi01 ** 2))

    # ---- assemble sigma^2 ----
    sig2 = (2 * alpha2 + 4 * beta2 + 8 * gam_sh + 4 * gam_dj + 16 * delta
            - 2 * alpha1 ** 2 - 16 * alpha1 * beta1 - 16 * beta1 ** 2)
    print("\n  --- moments ---")
    print(f"    alpha1={alpha1:.6f} alpha2={alpha2:.6f}  beta1={beta1:.6f} beta2={beta2:.6f}")
    print(f"    delta={delta:.6f}  gam_sh={gam_sh:.6f}  gam_dj={gam_dj:.6f}")
    print(f"\n  >>> sigma^2 (analytic, RI-free) = {sig2:.6f}   vs VMC {SIG2_REF:.6f}"
          f"   (rel {abs(sig2-SIG2_REF)/SIG2_REF:.1e})")

    # ========================================================================= #
    # STAGE 4a (part 1): h_Vne = <Phi0| V_ne (F - Fbar) |Phi0> = <V_ne F> - Fbar <V_ne>
    #   V_ne(r) = -Z_A/r_A - Z_B/r_B is a KNOWN grid multiplier -> no new dressing.
    #   <V_ne F> = 4 <v1 F>,  <v1 F> = <v1 f12> + <v1><f34> + 2<v1 f13> + 2<v1 f23>.
    # ========================================================================= #
    rA_f = _flat(rA); rB_f = _flat(rB)
    vne = -Z_A / rA_f - Z_B / rB_f                          # one-body V_ne on the grid
    Vne_exp = 2.0 * grid_int(rho_f * vne)                   # <sum_i v(i)> = 2 INT rho v
    vbar = 0.5 * grid_int(rho_f * vne)                      # <v>_1electron = INT rho_bar v (=Vne_exp/4)
    Pd = {'00': P00, '11': P11, '01': P01}
    Psid = {'00': Psi00, '11': Psi11, '01': Psi01}
    comps = [('00', '11', 1.0), ('11', '00', 1.0), ('01', '01', -2.0)]     # D = sum c g(x) h(y)
    Ta = 0.5 * sum(c * grid_int(vne * Pd[g] * Psid[h]) for g, h, c in comps)   # <v1 f12>
    Tb = vbar * alpha1                                                         # <v1><f34>
    Tc = 0.25 * grid_int(vne * rho_f * Psirho)                                 # <v1 f13>
    Td = 0.5 * sum(c * grid_int(vne * Pd[g]) * grid_int(Pd[h] * Psirho / 2.0)
                   for g, h, c in comps)                                       # <v1 f23>
    VneF = 4.0 * (Ta + Tb + 2 * Tc + 2 * Td)
    h_Vne = VneF - Fbar * Vne_exp
    print("\n  --- STAGE 4a part 1: h_Vne (analytic, RI-free) ---")
    print(f"    <V_ne> (grid) = {Vne_exp:.5f}   (Stage-2 analytic -20.86474; cross-check)")
    print(f"    Ta={Ta:.5f} Tb={Tb:.5f} Tc={Tc:.5f} Td={Td:.5f}  ->  <V_ne F>={VneF:.5f}")
    print(f"    h_Vne = <V_ne F> - Fbar <V_ne> = {h_Vne:+.6f}   (VMC target: see lih_r12ci_vmc.py)")
    print("  (NEXT: h_T kinetic [gradient densities]; h_Vee + g share the 2-center Coulomb field.)")
