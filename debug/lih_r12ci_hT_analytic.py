"""LiH R12-CI Stage 4a (part 2): ANALYTIC (RI-free) h_T = <Phi0| T (F-Fbar) |Phi0>.

The kinetic block of the off-diagonal element h = H01.  VMC target (lih_r12ci_vmc.py,
geminal f=exp(-0.5 r)): h_T = +0.7485 +/- 0.004.  Built with NO Monte-Carlo (except one
importance-MC cross-check element) and NO resolution-of-identity, on the SAME separable
block density as sigma^2 (lih_r12ci_sigma2_analytic.py).

  h_T = <Phi0|T(F-Fbar)|Phi0>/S00 = PartA + PartB,   T = sum_i (-1/2 grad_i^2),  chi=F-Fbar.
  Green/IBP (exact):  <Phi0|T|chi Phi0> = 1/2 sum_i INT (chi |grad_i Phi0|^2 + Phi0 grad_i F . grad_i Phi0)
    PartA = 1/2 <chi sum_i |grad_i Phi0|^2>   (the drift^2 term),
    PartB = 1/2 <sum_i grad_i F . grad_i Phi0 / Phi0>   (the grad F . grad Phi0 coupling).

PART A -- separable, like sigma^2 but with GRADIENT one-electron densities G_pq = grad m_p . grad m_q:
  sum_i |grad_i Phi0|^2 = KD(1,2) D_p(3,4) + D_p(1,2) KD(3,4),
  KD(a,b) = |grad_a D|^2 + |grad_b D|^2 = G00(a)P11(b)+G11(a)P00(b)-2G01(a)P01(b) + (a<->b),
  =>  PartA = KI + INT kappa Psi^f_rho + (alpha1 - Fbar)(T00+T11),
    KI    = I_f[G00,P11] + I_f[G11,P00] - 2 I_f[G01,P01]   (grid_int(G_pq * Psi^f_{P_rs})),
    kappa = (G00+G11) + T00 P11 + T11 P00 - 2 T01 P01,      T_pq = INT G_pq (kinetic overlap).
  Reuses the sigma^2 dressing fields Psi^f (no new machinery).

PART B -- the coupling term.  Hermiticity/IBP collapses the VECTOR grad-f coupling to a SCALAR:
  1/2 sum_i INT Phi0 grad_i F . grad_i Phi0 = -1/4 sum_i INT Phi0^2 grad_i^2 F,
  sum_i grad_i^2 F = 2 sum_{i<j} (grad^2 f)(r_ij),   grad^2 f = f'' + (2/r) f' = gam^2 f - 2 gam (e^{-gam r}/r).
  =>  PartB = -1/2 <sum_{i<j} grad^2 f(r_ij)> = -1/2 gam^2 Fbar + gam * Ybar,
    Ybar = <sum_{i<j} Y(r_ij)>,  Y(r) = e^{-gam r}/r  (a Yukawa/screened-Coulomb energy).
  Ybar via the 3x3 Yukawa AO-pair matrix W^Y (same Fbar-structure): the isotropic dressings
  I_Y[aa,*], I_Y[bb,*] use the closed radial Yukawa potential of a 1s density (no grid diagonal
  singularity); only I_Y[ab,ab] (both densities 2-center) is an importance-MC, as (ab|ab) was.

Run from debug/:  python lih_r12ci_hT_analytic.py
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

from lih_r12ci_sigma2_analytic import (
    Xg, Eg, rA, rB, geo_f, dens, build_kernel, Wmat, grid_int, X, a, ZA, ZB, GAM)
from lih_r12ci_energy import stage2_E0, _hartree_1s, N_A, N_B
from lih_r12_4body_integral import sample_1s, CENTER_A, CENTER_B
from lih_r12ci_sigma2_mc import sample_block

HT_REF = 0.7485        # VMC target (lih_r12ci_vmc.py), geminal exp(-0.5 r)
rng = np.random.default_rng(20260921)

# --------------------------------------------------------------------------- #
# flat grid fields
# --------------------------------------------------------------------------- #
rA_f = rA.reshape(-1); rB_f = rB.reshape(-1)
d_aa, d_ab, d_bb = dens['aa'], dens['ab'], dens['bb']
# cos(angle between r_A and r_B directions) = (xi^2+eta^2-2)/(xi^2-eta^2), clipped
cosAB = np.clip((Xg ** 2 + Eg ** 2 - 2.0) / np.maximum(Xg ** 2 - Eg ** 2, 1e-30),
                -1.0, 1.0).reshape(-1)


def Gpq(p, q):
    """gradient one-electron density G_pq = grad m_p . grad m_q on the grid.
    grad orb_A = -ZA orb_A rhat_A -> |grad orb_A|^2=ZA^2 rho_A, grad orb_A.grad orb_B=ZA ZB cosAB rho_ab."""
    return (X[0, p] * X[0, q] * ZA ** 2 * d_aa
            + (X[0, p] * X[1, q] + X[1, p] * X[0, q]) * ZA * ZB * cosAB * d_ab
            + X[1, p] * X[1, q] * ZB ** 2 * d_bb)


def cvec(p, q):
    return np.array([X[0, p] * X[0, q], X[0, p] * X[1, q] + X[1, p] * X[0, q], X[1, p] * X[1, q]])


# --------------------------------------------------------------------------- #
# Yukawa potential of a 1s(Z) density: closed radial form (angle-averaged kernel)
#   V^Y(s) = (2 Z^3/(gam s)) INT_0^inf r2 e^{-2Z r2} [e^{-gam|s-r2|} - e^{-gam(s+r2)}] dr2
# angle-average of e^{-gam d}/d over the sphere at radius r2 = (1/(2 s r2 gam))[e^{-gam|s-r2|}-e^{-gam(s+r2)}]
# --------------------------------------------------------------------------- #
_NR = 800; _Rmax = 45.0
_xr, _wr = leggauss(_NR); _r2 = 0.5 * _Rmax * (_xr + 1.0); _w2 = 0.5 * _Rmax * _wr


def yukawa_pot_iso(s, Z, gam):
    s = np.asarray(s, float); sc = np.maximum(s, 1e-12)
    r2 = _r2[None, :]; ss = sc[:, None]
    ker = np.exp(-gam * np.abs(ss - r2)) - np.exp(-gam * (ss + r2))
    integ = (r2 * np.exp(-2.0 * Z * r2) * ker) @ _w2
    return (2.0 * Z ** 3 / (gam * sc)) * integ


def mc_yukawa_abab(gam, n=12_000_000, batch=3_000_000):
    """I_Y[ab,ab] = INT rho_ab(1) rho_ab(2) e^{-gam r12}/r12, importance from |1s_A|^2
    (weight u=rho_ab/rho_A=(N_B/N_A)e^{ZA rA - ZB rB}).  Screened -> milder tail than Coulomb."""
    means = []; ntot = 0
    while ntot < n:
        r1 = sample_1s(batch, ZA, CENTER_A); r2 = sample_1s(batch, ZA, CENTER_A)
        rA1 = np.linalg.norm(r1 - CENTER_A, axis=1); rB1 = np.linalg.norm(r1 - CENTER_B, axis=1)
        rA2 = np.linalg.norm(r2 - CENTER_A, axis=1); rB2 = np.linalg.norm(r2 - CENTER_B, axis=1)
        u1 = (N_B / N_A) * np.exp(ZA * rA1 - ZB * rB1)
        u2 = (N_B / N_A) * np.exp(ZA * rA2 - ZB * rB2)
        d = np.maximum(np.linalg.norm(r1 - r2, axis=1), 1e-12)
        means.append((u1 * u2 * np.exp(-gam * d) / d).mean()); ntot += batch
    bm = np.array(means); return bm.mean(), bm.std(ddof=1) / np.sqrt(len(bm))


def ybar_mc(gam, seed_u=11, seed_d=12):
    """Independent MC of Ybar = <sum_{i<j} e^{-gam r_ij}/r_ij> over the block density |Phi0|^2."""
    cfg = dict(nw=10000, burn=3000, nsnap=12, thin=40)
    up, _ = sample_block(seed=seed_u, **cfg); dn, _ = sample_block(seed=seed_d, **cfg)
    ns = len(up); pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    vals = []
    for s in range(ns):
        pos = np.concatenate([up[s], dn[s]], axis=1)
        Y = sum(np.exp(-gam * np.maximum(np.linalg.norm(pos[:, i] - pos[:, j], axis=-1), 1e-12))
                / np.maximum(np.linalg.norm(pos[:, i] - pos[:, j], axis=-1), 1e-12) for i, j in pairs)
        vals.append(Y.mean())
    v = np.array(vals); return v.mean(), v.std(ddof=1) / np.sqrt(ns)


if __name__ == "__main__":
    print("=" * 78)
    print("LiH R12-CI Stage 4a (part 2): ANALYTIC h_T (RI-free), target +0.7485")
    print(f"  f=exp(-{GAM} r);  grid from sigma2_analytic")
    print("=" * 78)

    # ---- gradient densities + kinetic overlaps ----
    G00, G11, G01 = Gpq(0, 0), Gpq(1, 1), Gpq(0, 1)
    T00, T11, T01 = grid_int(G00), grid_int(G11), grid_int(G01)

    # cross-check: T00+T11 = <Phi0|sum_i T_i|Phi0>  (stage-2 T_exp)
    _, info = stage2_E0(mc_check=False)
    T_exp = info['T_exp']
    print(f"\n  GATE A1 -- kinetic overlaps:  T00={T00:.5f} T11={T11:.5f} T01={T01:.5f}")
    print(f"    T00+T11 = {T00+T11:.5f}   vs stage-2 <T>_Phi0 = {T_exp:.5f}"
          f"   (rel {abs(T00+T11-T_exp)/abs(T_exp):.1e})")

    # ---- f-interaction matrix + dressing fields (same as sigma^2) ----
    Kf = build_kernel(GAM); W, Psi = Wmat(Kf)
    c00, c11, c01 = cvec(0, 0), cvec(1, 1), cvec(0, 1); crho = c00 + c11
    Psi_field = lambda c: c[0] * Psi['aa'] + c[1] * Psi['ab'] + c[2] * Psi['bb']
    Psi00, Psi11, Psi01, Psirho = (Psi_field(c00), Psi_field(c11), Psi_field(c01), Psi_field(crho))
    P00 = c00[0] * d_aa + c00[1] * d_ab + c00[2] * d_bb
    P11 = c11[0] * d_aa + c11[1] * d_ab + c11[2] * d_bb
    P01 = c01[0] * d_aa + c01[1] * d_ab + c01[2] * d_bb

    alpha1 = c00 @ W @ c11 - c01 @ W @ c01
    beta1 = 0.25 * (crho @ W @ crho)
    Fbar = 2 * alpha1 + 4 * beta1

    # ---- PART A ----
    KI = grid_int(G00 * Psi11) + grid_int(G11 * Psi00) - 2 * grid_int(G01 * Psi01)
    kappa = (G00 + G11) + T00 * P11 + T11 * P00 - 2 * T01 * P01
    PartA = KI + grid_int(kappa * Psirho) + (alpha1 - Fbar) * (T00 + T11)
    print(f"\n  --- PART A (drift^2 term) ---")
    print(f"    KI={KI:.5f}  INT kappa*Psirho={grid_int(kappa*Psirho):.5f}  "
          f"(alpha1-Fbar)(T00+T11)={(alpha1-Fbar)*(T00+T11):.5f}")
    print(f"    PartA = {PartA:+.6f}")

    # ---- PART B ----
    # validate yukawa_pot_iso -> Hartree as gam->0
    s_test = np.array([0.3, 1.0, 2.0, 4.0])
    yk = yukawa_pot_iso(s_test, ZA, 1e-3); hh = _hartree_1s(s_test, ZA)
    print(f"\n  GATE B1 -- yukawa_pot_iso(gam->0) vs Hartree(1s_A): "
          f"max rel {np.max(np.abs(yk-hh)/np.abs(hh)):.1e}")

    VYA = yukawa_pot_iso(rA_f, ZA, GAM); VYB = yukawa_pot_iso(rB_f, ZB, GAM)
    IY_aaaa = grid_int(d_aa * VYA); IY_bbbb = grid_int(d_bb * VYB)
    IY_aabb = grid_int(d_bb * VYA); IY_aabb2 = grid_int(d_aa * VYB)   # symmetry check
    IY_aaab = grid_int(d_ab * VYA); IY_bbab = grid_int(d_ab * VYB)
    IY_abab, e_abab = mc_yukawa_abab(GAM)
    print(f"\n  --- PART B (Yukawa) ---")
    print(f"    I_Y: aaaa={IY_aaaa:.5f} bbbb={IY_bbbb:.5f} aabb={IY_aabb:.5f}"
          f"(={IY_aabb2:.5f}, rel {abs(IY_aabb-IY_aabb2)/abs(IY_aabb2):.1e})")
    print(f"         aaab={IY_aaab:.5f} bbab={IY_bbab:.5f}  abab={IY_abab:.5f}+/-{e_abab:.1e} [MC]")
    # AO-pair order MUST match cvec / W^f: [aa, ab, bb]  (not [aa, bb, ab]!)
    WY = np.array([[IY_aaaa, IY_aaab, IY_aabb],
                   [IY_aaab, IY_abab, IY_bbab],
                   [IY_aabb, IY_bbab, IY_bbbb]])
    alpha1Y = c00 @ WY @ c11 - c01 @ WY @ c01
    beta1Y = 0.25 * (crho @ WY @ crho)
    Ybar = 2 * alpha1Y + 4 * beta1Y
    PartB = -0.5 * GAM ** 2 * Fbar + GAM * Ybar
    print(f"    alpha1^Y={alpha1Y:.5f} beta1^Y={beta1Y:.5f}  Ybar=2a+4b={Ybar:.5f}")
    print(f"    PartB = -1/2 gam^2 Fbar + gam Ybar = {-0.5*GAM**2*Fbar:.5f} + {GAM*Ybar:.5f}"
          f" = {PartB:+.6f}")

    ymc, eymc = ybar_mc(GAM)
    print(f"    GATE B2 -- Ybar (analytic) {Ybar:.5f} vs MC {ymc:.5f}+/-{eymc:.1e}"
          f"  (dev {abs(Ybar-ymc)/eymc:.1f} sigma)")

    # ---- TOTAL ----
    h_T = PartA + PartB
    print(f"\n  >>> h_T = PartA + PartB = {PartA:+.6f} + ({PartB:+.6f}) = {h_T:+.6f}")
    print(f"      VMC target = {HT_REF:+.4f}   (dev {abs(h_T-HT_REF):.4f})")
