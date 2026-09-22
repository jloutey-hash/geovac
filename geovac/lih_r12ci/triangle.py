"""GATE: the two-center 3-body TRIANGLE reduction (the owed piece of g_Vee), validated
reduced==MC BEFORE any energy assembly (per debug/lih_r12_build_plan.md Step-3 discipline).

The triangle is the only non-separable term class in <F^2 V_ee>: two f-geminals sharing a
vertex, with the Coulomb on the opposite edge of the 3-cycle. Two block-placements:

  T1 = <f13 f23 coul12>  (coul on the INTRA edge {1,2}; f's share vertex 3)
       -> integrate spectator 4 -> shared vertex 3 carries the MARGINAL rho;
          the two Coulomb electrons 1,2 carry the D_p(1,2) AO-pair components.
  T2 = <f12 f23 coul13>  (coul on an INTER edge {1,3}; f's share vertex 2)
       -> integrate spectator 4 -> shared vertex 2 carries a D_p AO-pair component h,
          Coulomb electron 1 carries partner g, electron 3 carries the marginal rho.

Both reduce to a NON-SEPARABLE 2-point kernel + a per-mode prolate-Neumann Coulomb:
    <..> = (1/4) sum_k kappa_k INT left_k(1) right_k(2) coul(r12) K_mu(1,2) d1 d2,
    K_mu(1,2) = INT mu(3) f(r13) f(r23) d3      (mu = the middle/shared-vertex density),
azimuthal-Fourier in the RELATIVE azimuth:
    K_mu(1,2)   = sum_m K_mu^{(m)}(x1,x2) cos(m(phi1-phi2)),
    K_mu^{(m)}  = d_m a^3 sum_i3 geo_i3 mu_i3 F_m[.,i3] F_m[.,i3],   d_0=2pi, d_{m>=1}=pi,
    F_m[i,j]    = m-th cos-Fourier coeff of exp(-gam|r_i - r_j(phi')|)   (build_kernel_m),
    coul(r12)   = sum_m coul^{(m)}(x1,x2) cos(m(phi1-phi2))            (prolate Neumann),
and the phi1,phi2 integral diagonalizes in m with weight W_0=(2pi)^2, W_{m>=1}=2pi^2.
Low-rank (eigh) K_mu^{(m)} = sum_r lam_r v_r(x1) v_r(x2)  =>  each term a per-mode Coulomb
E_m[left*v_r, right*v_r].  The 1/4 = |Phi0|^2 = D_p D_p / 4 (INT D_p = 2 per block).

Run from debug/:  python lih_r12ci_triangle_gate.py
"""
import warnings

import numpy as np
from numpy.polynomial.legendre import leggauss

warnings.filterwarnings("ignore", category=DeprecationWarning)
from scipy.special import lpmn, lqmn  # noqa: E402

from .kernels import (   # noqa: E402
    Xg, Eg, geo_f, a, GAM, grid_int, rho_cyl_f, zc_f,
    XI as XI1D, ETA as ETA1D, WXI, WETA)
from .hVee import neumann_potential, P00, P11, P01, rho_g  # noqa: E402
from .energy import R  # noqa: E402
from .basis import sample_block  # noqa: E402

NXI, NETA = XI1D.size, ETA1D.size
NG = NXI * NETA
JAC2 = (Xg ** 2 - Eg ** 2)                       # (NXI,NETA)
GW2 = WXI[:, None] * WETA[None, :] * JAC2         # 2D "grid2" weight (no 2pi a^3), (NXI,NETA)
LMAX = 34
MMAX = 4
rng = np.random.default_rng(20260921)

# --- Legendre tables on the sigma2 grid axes (P_l^m(eta), P_l^m(xi), Q_l^m(xi)) --- #
Peta = np.zeros((MMAX + 1, LMAX + 1, NETA))
for i, e in enumerate(ETA1D):
    Peta[:, :, i] = lpmn(MMAX, LMAX, e)[0]
Pxi = np.zeros((MMAX + 1, LMAX + 1, NXI)); Qxi = np.zeros((MMAX + 1, LMAX + 1, NXI))
for i, x in enumerate(XI1D):
    Pxi[:, :, i] = lpmn(MMAX, LMAX, x)[0]
    Qxi[:, :, i] = lqmn(MMAX, LMAX, x)[0]
_minidx = np.minimum.outer(np.arange(NXI), np.arange(NXI))
_maxidx = np.maximum.outer(np.arange(NXI), np.arange(NXI))
import math
_NORM = {(m, l): (math.factorial(l - m) / math.factorial(l + m)) ** 2
         for m in range(MMAX + 1) for l in range(m, LMAX + 1)}


# --------------------------------------------------------------------------- #
# azimuthal-Fourier component F_m[i,j] of the geminal on the sigma2 grid
# --------------------------------------------------------------------------- #
def build_kernel_m(gam, m, nphi=48):
    """F_m[i,j] = (2-delta_{m0})^{-1}(1/pi) INT_0^{2pi} exp(-gam d(phi')) cos(m phi') dphi'
    i.e. the coeff of cos(m*Delta_phi):  F_0=(1/2pi)INT f dphi',  F_{m>=1}=(1/pi)INT f cos(m phi')dphi'.
    (integrand even in phi' -> 2*INT_0^pi)."""
    xp, wp = leggauss(nphi)
    phi = 0.5 * np.pi * (xp + 1.0); wphi = 0.5 * np.pi * wp     # [0,pi]
    rc = rho_cyl_f; z = zc_f
    rc2 = rc[:, None] * rc[None, :]
    base = rc[:, None] ** 2 + rc[None, :] ** 2 + (z[:, None] - z[None, :]) ** 2
    K = np.zeros((NG, NG))
    for ph, w in zip(phi, wphi):
        d = np.sqrt(np.maximum(base - 2.0 * rc2 * np.cos(ph), 0.0))
        K += 2.0 * w * np.exp(-gam * d) * np.cos(m * ph)       # INT_0^{2pi} f cos(m phi')
    return K / (2 * np.pi) if m == 0 else K / np.pi


# --------------------------------------------------------------------------- #
# mode-m prolate-Neumann Coulomb potential of a 2D density-component B (coeff of cos(m Dphi))
#   V^{(m)}_B(xi1,eta1) = (2-d_m0)(2/R)a^6 W_m sum_l (-1)^m (2l+1) norm P_l^m(eta1)
#                          sum_xi2 wxi K_l(xi1,xi2) gB_l(xi2),  gB_l = sum_eta weta JAC B P_l^m(eta)
# and  E_m[A,B] = sum_{xi,eta} wxi weta JAC A * V^{(m)}_B  (grid2 contraction, no extra 2pi a^3).
# GATE: E_0[A,B] == grid_int(A * neumann_potential(B)) = INT A B / r12  (the m=0 machinery).
# --------------------------------------------------------------------------- #
def coul_mode_potential(B2d, m):
    Wm = (2 * np.pi) ** 2 if m == 0 else 2 * np.pi ** 2
    pref = (2.0 - (1.0 if m == 0 else 0.0)) * (2.0 / R) * a ** 6 * Wm
    V = np.zeros((NXI, NETA))
    WB = JAC2 * B2d                                            # (NXI,NETA)
    for l in range(m, LMAX + 1):
        Plm_e = Peta[m, l]
        gB = (WB * Plm_e[None, :]) @ WETA                      # (NXI,)  sum over eta
        Kl = Pxi[m, l][_minidx] * Qxi[m, l][_maxidx]           # (NXI,NXI) P_l^m(xi<)Q_l^m(xi>)
        radial = Kl @ (WXI * gB)                               # (NXI,)
        V += (-1) ** m * (2 * l + 1) * _NORM[(m, l)] * np.outer(radial, Plm_e)
    return pref * V


def grid2_int(F2d):
    return np.sum(GW2 * F2d)


def E_m(A2d, B2d, m):
    return grid2_int(A2d * coul_mode_potential(B2d, m))


# --------------------------------------------------------------------------- #
# low-rank eigen-decomposition of K_mu^{(m)} = d_m a^3 F_m diag(geo*mu) F_m^T
#   (randomized eigh: symmetric, only the top ~rank modes matter; cached by (key,m))
# --------------------------------------------------------------------------- #
_MODE_CACHE = {}


def Kmu_modes(Fm, mu_flat, m, key, rank=30, oversample=16, power=2, tol=1e-11):
    ck = (key, m, rank)
    if ck in _MODE_CACHE:
        return _MODE_CACHE[ck]
    d_m = 2 * np.pi if m == 0 else np.pi
    coef = d_m * a ** 3
    gmu = geo_f * mu_flat

    def matvec(Xc):                                     # K @ Xc  (Xc: NG x cols)
        return coef * (Fm @ (gmu[:, None] * (Fm @ Xc)))
    p = rank + oversample
    Om = rng.standard_normal((NG, p))
    Y = matvec(Om)
    for _ in range(power):
        Y = matvec(Y)
    Q, _ = np.linalg.qr(Y)
    Bmat = Q.T @ matvec(Q); Bmat = 0.5 * (Bmat + Bmat.T)
    w, U = np.linalg.eigh(Bmat)
    idx = np.argsort(np.abs(w))[::-1][:rank]
    w = w[idx]; Vsm = U[:, idx]
    keep = np.abs(w) > tol * np.abs(w).max()
    lam, V = w[keep], (Q @ Vsm)[:, keep]
    _MODE_CACHE[ck] = (lam, V)
    return lam, V


def triangle_raw(Fm_by_m, mu_flat, left_flat, right_flat, mu_key, mmax=MMAX, rank=30):
    """INT left(1) right(2) coul(r12) K_mu(1,2) d1 d2 (RAW, no 1/4, no kappa)."""
    tot = 0.0; per_m = []
    for m in range(mmax + 1):
        lam, V = Kmu_modes(Fm_by_m[m], mu_flat, m, mu_key, rank=rank)
        acc = 0.0
        for r in range(lam.size):
            vr = V[:, r]
            L = (left_flat * vr).reshape(NXI, NETA)
            Rt = (right_flat * vr).reshape(NXI, NETA)
            acc += lam[r] * E_m(L, Rt, m)
        per_m.append(acc); tot += acc
    return tot, per_m


# --------------------------------------------------------------------------- #
# Monte-Carlo ground truth for the two triangle expectations over |Phi0|^2
# --------------------------------------------------------------------------- #
def _mc_triangles(seed_u=71, seed_d=72, nw=14000, burn=3000, nsnap=18, thin=45):
    up, au = sample_block(seed=seed_u, nw=nw, burn=burn, nsnap=nsnap, thin=thin)
    dn, ad = sample_block(seed=seed_d, nw=nw, burn=burn, nsnap=nsnap, thin=thin)
    T1, T2 = [], []
    for s in range(len(up)):
        pos = np.concatenate([up[s], dn[s]], axis=1)     # (nw,4,3): 0,1 up ; 2,3 down
        d = lambda i, j: np.maximum(np.linalg.norm(pos[:, i] - pos[:, j], axis=-1), 1e-12)
        f = lambda i, j: np.exp(-GAM * d(i, j))
        T1.append((f(0, 2) * f(1, 2) / d(0, 1)).mean())  # <f13 f23 coul12>
        T2.append((f(0, 1) * f(1, 2) / d(0, 2)).mean())  # <f12 f23 coul13>
    m = lambda x: np.mean(x); e = lambda x: np.std(x, ddof=1) / np.sqrt(len(x))
    return (m(T1), e(T1)), (m(T2), e(T2)), (au, ad)


if __name__ == "__main__":
    print("=" * 78)
    print("TRIANGLE GATE: two-center 3-body reduction (T1 intra-coul, T2 inter-coul)")
    print(f"  sigma2 grid {NXI}x{NETA}, LMAX={LMAX}, MMAX={MMAX}, f=exp(-{GAM} r), R={R}")
    print("=" * 78)

    # ---- Fourier kernels ----
    Fm = [build_kernel_m(GAM, m) for m in range(MMAX + 1)]

    # ---- GATE 0: E_0 vs the m=0 neumann_potential machinery ----
    A = (P00).reshape(NXI, NETA); Bt = (rho_g).reshape(NXI, NETA)
    e0 = E_m(A, Bt, 0)
    ref = grid_int(P00 * neumann_potential(rho_g.reshape(NXI, NETA)).reshape(-1))
    print(f"\n  GATE 0 -- E_0[P00,rho] = {e0:.6f}  vs grid_int(P00*neumann_pot(rho)) = {ref:.6f}"
          f"   (rel {abs(e0-ref)/abs(ref):.1e})")

    # ---- MC ground truth ----
    (t1mc, e1), (t2mc, e2), (au, ad) = _mc_triangles()
    print(f"\n  MC (block sampler, acc {au:.2f}/{ad:.2f}):")
    print(f"    T1=<f13 f23 coul12> = {t1mc:.6f} +/- {e1:.1e}")
    print(f"    T2=<f12 f23 coul13> = {t2mc:.6f} +/- {e2:.1e}")

    # ---- T1: mu=rho on shared vertex 3; Coulomb electrons 1,2 carry D_p(1,2) comps ----
    # D_p(1,2) = sum kappa g(1) h(2): (P00,P11,+1),(P11,P00,+1),(P01,P01,-2)
    comps = [(P00, P11, 1.0, 'P00', 'P11'), (P11, P00, 1.0, 'P11', 'P00'),
             (P01, P01, -2.0, 'P01', 'P01')]
    T1_raw = 0.0
    for g, h, k, gk, hk in comps:
        v, _ = triangle_raw(Fm, rho_g, g, h, 'rho')
        T1_raw += k * v
    T1_an = 0.25 * T1_raw
    print(f"\n  T1 reduced (analytic, RI-free) = {T1_an:.6f}   vs MC {t1mc:.6f}"
          f"   (dev {abs(T1_an-t1mc)/e1:.1f} sigma, rel {abs(T1_an-t1mc)/abs(t1mc):.1e})")

    # ---- T2: shared vertex 2 carries D_p comp h; electron 1 carries g; electron 3 carries rho ----
    #   <f12 f23 coul13> = 1/4 sum_kappa INT g(1) rho(3) coul(r13) K_h(1,3);  K_h = INT h(2) f12 f23 d2
    T2_raw = 0.0
    for g, h, k, gk, hk in comps:
        v, _ = triangle_raw(Fm, h, g, rho_g, hk)  # mu=h (shared vtx 2), left=g (e1), right=rho (e3)
        T2_raw += k * v
    T2_an = 0.25 * T2_raw
    print(f"  T2 reduced (analytic, RI-free) = {T2_an:.6f}   vs MC {t2mc:.6f}"
          f"   (dev {abs(T2_an-t2mc)/e2:.1f} sigma, rel {abs(T2_an-t2mc)/abs(t2mc):.1e})")

    ok = abs(T1_an - t1mc) < 5 * e1 and abs(T2_an - t2mc) < 5 * e2
    print(f"\n  VERDICT: triangle reduction {'VALIDATED' if ok else 'MISMATCH -- debug'} "
          f"(both within 5 sigma of MC)")
