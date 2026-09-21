"""Increment 6: r12-COUPLED core-Hartree V_H and Huzinaga projector for
frozen-core prolate LiH.

The p=0 add-ons (lih_frozen_core_first._twoelec_VH, lih_projector_test.build)
factorize the 2e matrix element as (1e term) (x) S + S (x) (1e term) because at
p=0 the basis function g_i(1,2) = g_i^(1)(1) g_i^(2)(2) is a product.  With an
r12^p Hylleraas factor (p>=1) the two electrons are correlated and this
factorization FAILS.  This module rebuilds both terms so LiH can be run WITH r12
(p_set={0,1}), apples-to-apples with the HeH+ PoC that converged (-4.5%->-0.5%).

Route (per the track log): 2e prolate quadrature, vectorized over COMBINED powers
(the matrix element depends only on j_i+j_j etc., not on the pair), so it is fast
even at l=3 / p={0,1}.

    phi-kernels (mu=0 sigma; r12 depends on phi only via psi = phi1 - phi2):
      Phi_P[a,b] = INT INT dphi1 dphi2 r12^P   (both electrons' phi integrated)
        P=0 -> (2pi)^2
        P=1 -> 2pi * hR * 4 sqrt(A+B) E(m),  m = 2B/(A+B)   (elliptic E)
        P=2 -> (2pi)^2 hR^2 A[a,b]                          (cos term integrates 0)
      J_p[a,b]   = INT dphi1 r12^p            (single phi; = Phi_p / 2pi)
        J_0 -> 2pi ;  J_1 -> hR * 4 sqrt(A+B) E(m)

    V_H (local, summed over electrons):
      VH_ij = Q1_P[(j_i+j_j,l_i+l_j),(k_i+k_j,m_i+m_j)]        (V_H on electron 1)
            + Q2_P[...]                                        (V_H on electron 2)
      Q1_P[(J,L),(K,M)] = (w V_H xi^J eta^L e^{-2a xi}) . Phi_P . (w xi^K eta^M e^{-2a xi})

    Projector P_1s(1)+P_1s(2), non-local:
      P1_ij = (2pi/N1s) SUM_b w_b [xi_b^{k_i}eta_b^{m_i}e^{-a xi_b}]
                                  [xi_b^{k_j}eta_b^{m_j}e^{-a xi_b}] Btil_i[b] Btil_j[b]
      Btil_i[b] = SUM_a w_a [xi_a^{j_i}eta_a^{l_i}e^{-a xi_a}] 1s[a] J_{p_i}[a,b]
      (P2 = electron-2 projected, (j,l)<->(k,m) swapped).  N1s = <1s|1s>.

VALIDATION (run this file): on a p=0-only basis both must reproduce the existing
p=0 references to grid accuracy, and the quad overlap must reproduce the mpf S.

Run from root:  python debug/lih_r12_coupled.py
"""
import os
import sys
import time
import numpy as np
from scipy.special import ellipe

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from geovac import prolate_recondition as pr                 # noqa: E402
from geovac.hylleraas import build_quadrature_grids          # noqa: E402
import prolate_r12_mpf as m                                  # noqa: E402
from heh_probe import build_basis_full                       # noqa: E402
from prolate_core_hartree import V_H_closed, ZC_LI           # noqa: E402

TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------------------
# Flattened (xi,eta) grid + per-point arrays, and the phi-integrated kernels.
# ---------------------------------------------------------------------------
def grid_arrays(R, alpha, N_xi=32, N_eta=24, xi_max=18.0, zc=ZC_LI):
    """Flattened (xi,eta) quadrature grid with per-point weights and factors.
    Returns a dict of length-ng arrays and the ng x ng phi-kernels Phi_{0,1,2}
    (double-phi) and J_{0,1} (single-phi)."""
    g = build_quadrature_grids(N_xi=N_xi, N_eta=N_eta, N_phi=4, xi_max=xi_max)
    xi, wxi, eta, weta = g['xi'], g['w_xi'], g['eta'], g['w_eta']
    hR = R / 2.0
    XI, ETA = np.meshgrid(xi, eta, indexing='ij')
    WXI, WETA = np.meshgrid(wxi, weta, indexing='ij')
    x = XI.ravel(); e = ETA.ravel()
    w = (WXI * WETA).ravel() * (x**2 - e**2) * (hR**3)      # per-electron weight (no phi)
    rho = np.sqrt(np.clip((x**2 - 1.0) * (1.0 - e**2), 0.0, None))
    r_A = hR * (x + e)
    vH = V_H_closed(np.maximum(r_A, 1e-12), zc)
    ea = np.exp(-alpha * x)                                  # single e^{-a xi}
    e2a = np.exp(-2.0 * alpha * x)                           # e^{-2a xi}
    s1s = np.exp(-zc * r_A)                                  # 1s_A on the grid
    # phi kernels
    xe = x * e
    A = (xe[:, None] - xe[None, :])**2 + rho[:, None]**2 + rho[None, :]**2
    B = 2.0 * rho[:, None] * rho[None, :]
    mparam = np.clip(2.0 * B / (A + B + 1e-300), 0.0, 1.0)
    sqrtApB = np.sqrt(A + B)
    J1 = hR * 4.0 * sqrtApB * ellipe(mparam)                # INT dphi1 r12
    Phi0 = (TWO_PI**2) * np.ones_like(A)
    Phi1 = TWO_PI * J1
    Phi2 = (TWO_PI**2) * (hR**2) * A
    return dict(x=x, e=e, w=w, vH=vH, ea=ea, e2a=e2a, s1s=s1s, hR=hR,
                J0=TWO_PI, J1=J1, Phi0=Phi0, Phi1=Phi1, Phi2=Phi2)


def _phi_of_P(G, P):
    return {0: G['Phi0'], 1: G['Phi1'], 2: G['Phi2']}[P]


# ---------------------------------------------------------------------------
# r12-coupled core-Hartree V_H  ->  <g_i r12^pi | V_H(r1A)+V_H(r2A) | g_j r12^pj>
# ---------------------------------------------------------------------------
def vh_coupled(basis_p, G, alpha):
    """2e core-Hartree matrix on a mixed p={0,1} basis, r12-coupled."""
    x, e, w, vH, e2a = G['x'], G['e'], G['w'], G['vH'], G['e2a']
    n = len(basis_p)
    # combined-power moment tables per P (P = p_i+p_j in {0,1,2})
    # need elec-1 (J,L) from (j,l) and elec-2 (K,M) from (k,m); powers up to 2*max
    jmax = max(max(b.j, b.k) for b, _ in basis_p)
    lmax = max(max(b.l, b.m) for b, _ in basis_p)
    Jset = range(0, 2 * jmax + 1)
    Lset = range(0, 2 * lmax + 1)
    JL = [(J, L) for J in Jset for L in Lset]
    idxJL = {jl: i for i, jl in enumerate(JL)}
    # U (V_H-carrying) and V (plain) power-vectors, shape [nJL, ng]
    Uw = np.array([w * vH * (x**J) * (e**L) * e2a for (J, L) in JL])
    Vw = np.array([w * (x**J) * (e**L) * e2a for (J, L) in JL])
    Q1 = {}; Q2 = {}
    for P in (0, 1, 2):
        Phi = _phi_of_P(G, P)
        Q1[P] = Uw @ Phi @ Vw.T          # V_H on electron 1
        Q2[P] = Vw @ Phi @ Uw.T          # V_H on electron 2
    W = np.zeros((n, n))
    for i, (bi, pi) in enumerate(basis_p):
        for jj, (bj, pj) in enumerate(basis_p):
            P = pi + pj
            e1 = idxJL[(bi.j + bj.j, bi.l + bj.l)]     # (J,L) elec-1
            e2 = idxJL[(bi.k + bj.k, bi.m + bj.m)]     # (K,M) elec-2
            W[i, jj] = Q1[P][e1, e2] + Q2[P][e1, e2]
    return W


# ---------------------------------------------------------------------------
# r12-coupled Huzinaga projector  ->  <g_i r12^pi | P_1s(1)+P_1s(2) | g_j r12^pj>
# ---------------------------------------------------------------------------
def projector_coupled(basis_p, G, alpha):
    """Normalized 2e core-orthogonality projector on a mixed p={0,1} basis."""
    x, e, w, ea, s1s = G['x'], G['e'], G['w'], G['ea'], G['s1s']
    J1 = G['J1']
    n = len(basis_p)
    N1s = TWO_PI * np.sum(w * s1s * s1s)                     # <1s|1s>
    wJ1 = (w[:, None] * J1)                                  # for p=1 Btil

    def Btil(jl_powers, p):
        """Btil_i[b] for a single-electron power (j,l) with r12 power p."""
        f1s = (x**jl_powers[0]) * (e**jl_powers[1]) * ea * s1s   # g^(1)*1s on grid
        if p == 0:
            return np.full_like(x, G['J0'] * np.sum(w * f1s))    # const in b
        return f1s @ wJ1                                         # SUM_a w_a f1s[a] J1[a,b]

    # electron-1 projector: project elec-1 (j,l), overlap elec-2 (k,m)
    # W2_i[b] = xi_b^{k_i} eta_b^{m_i} e^{-a xi_b} * Btil_i[b]
    W2_1 = np.empty((n, len(x)))
    W2_2 = np.empty((n, len(x)))
    for i, (bi, pi) in enumerate(basis_p):
        f2_1 = (x**bi.k) * (e**bi.m) * ea                   # elec-2 factor (proj on elec-1)
        W2_1[i] = f2_1 * Btil((bi.j, bi.l), pi)
        f2_2 = (x**bi.j) * (e**bi.l) * ea                   # elec-1 factor (proj on elec-2)
        W2_2[i] = f2_2 * Btil((bi.k, bi.m), pi)
    P1 = (TWO_PI / N1s) * (W2_1 * w) @ W2_1.T
    P2 = (TWO_PI / N1s) * (W2_2 * w) @ W2_2.T
    return P1 + P2


# ===========================================================================
# VALIDATION
# ===========================================================================
def _validate():
    import lih_frozen_core_first as L
    from lih_projector_test import build as build_p0
    R, alpha = 3.015, 1.0

    print("=" * 70)
    print("VALIDATION 1: r12-coupled V_H reduces to p=0 _twoelec_VH on a p=0 basis")
    basis0 = build_basis_full(2, 2, alpha, p_set=(0,))
    G = grid_arrays(R, alpha, N_xi=40, N_eta=30, xi_max=20.0)     # match p0 ref grid
    VH_new = vh_coupled(basis0, G, alpha)
    VH_ref = L._twoelec_VH(basis0, R, alpha)
    rel = np.max(np.abs(VH_new - VH_ref) / np.maximum(np.abs(VH_ref), 1e-10))
    print(f"  n={len(basis0)}  max rel diff = {rel:.3e}   diag[:4] new={np.round(np.diag(VH_new)[:4],5)}")
    print(f"                                     ref={np.round(np.diag(VH_ref)[:4],5)}")
    print("  V_H p=0 REDUCTION", "OK" if rel < 1e-9 else ("GRID-OK" if rel < 5e-3 else "FAIL"))

    print("\nVALIDATION 2: r12-coupled projector reduces to p=0 build() on a p=0 basis")
    _, P_ref = build_p0(basis0, R, alpha, 1.0)                    # lam=1 -> raw P
    P_new = projector_coupled(basis0, G, alpha)
    rel2 = np.max(np.abs(P_new - P_ref) / np.maximum(np.abs(P_ref), 1e-10))
    print(f"  max rel diff = {rel2:.3e}   diag[:4] new={np.round(np.diag(P_new)[:4],5)}")
    print(f"                                ref={np.round(np.diag(P_ref)[:4],5)}")
    print("  PROJECTOR p=0 REDUCTION", "OK" if rel2 < 1e-9 else ("GRID-OK" if rel2 < 5e-3 else "FAIL"))

    print("\nVALIDATION 3: quad overlap of the MIXED p={0,1} basis vs mpf engine S")
    basis01 = build_basis_full(1, 1, alpha, p_set=(0, 1))
    Squad = m.unsym_overlap_quad(basis01, R, alpha)
    S_mpf, _ = m.assemble_hetero(basis01, R, alpha, L.Z_LI, L.Z_H, l_neumann=20, dps=30)
    relS = np.max(np.abs(Squad - S_mpf) / np.maximum(np.abs(S_mpf), 1e-9))
    print(f"  n={len(basis01)}  max rel diff (quad vs mpf) = {relS:.3e}  (grid-limited ~1e-3)")
    print("  BASIS-CONSISTENCY", "OK" if relS < 1e-2 else "CHECK")


if __name__ == '__main__':
    t0 = time.time()
    _validate()
    print(f"\n[{time.time()-t0:.0f}s]")
