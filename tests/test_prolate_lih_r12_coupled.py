"""Guard for the r12-COUPLED core-Hartree V_H and Huzinaga projector
(debug/lih_r12_coupled.py, CHANGELOG v5.15.6) used to answer the frozen-core
prolate-LiH R_eq question.

Two independent things are guarded, each fire-tested against a specific wrong answer:

  (1) p=0 REDUCTION.  On a p=0-only basis both coupled terms must reproduce the
      pre-existing p=0 factorized references (_twoelec_VH, lih_projector_test.build).
      Rejects: a broken assembly / normalization / Phi_0 constant (e.g. Phi_0=2pi
      instead of (2pi)^2 would put V_H off by 2pi; a missing /<1s|1s> would break
      the projector).

  (2) r12 COUPLING (the load-bearing new machinery).  On a mixed p={0,1} basis the
      combined-power / elliptic-kernel V_H must match an INDEPENDENT direct 2e
      prolate quadrature that evaluates r12^P by an explicit phi loop.
      Rejects: a wrong Phi_1 (elliptic-E) kernel, a wrong power-shift, or a wrong
      P=p_i+p_j dispatch -- none of which the p=0 reduction can see.

Fast (small basis, coarse grid, no mpf): ~a few seconds.
"""
import os
import sys

import numpy as np
import pytest

_DEBUG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "debug")
if _DEBUG not in sys.path:
    sys.path.insert(0, _DEBUG)

# Several debug/ drivers parse sys.argv[1] as an int at import time (they are CLI
# scripts); strip pytest's argv so importing them does not raise.
sys.argv = [sys.argv[0]]

from geovac import prolate_recondition as pr                     # noqa: E402
import lih_r12_coupled as C                                       # noqa: E402
from heh_probe import build_basis_full                            # noqa: E402

R, ALPHA = 3.015, 1.0


def test_vh_coupled_reduces_to_p0_reference():
    """(1) r12-coupled V_H == p=0 factorized _twoelec_VH on a p=0 basis."""
    import lih_frozen_core_first as L
    basis0 = build_basis_full(1, 1, ALPHA, p_set=(0,))
    G = C.grid_arrays(R, ALPHA, N_xi=40, N_eta=30, xi_max=20.0)
    new = C.vh_coupled(basis0, G, ALPHA)
    ref = L._twoelec_VH(basis0, R, ALPHA)
    rel = np.max(np.abs(new - ref) / np.maximum(np.abs(ref), 1e-10))
    assert rel < 5e-3, f"V_H p=0 reduction failed: max rel {rel:.2e}"


def test_projector_coupled_reduces_to_p0_reference():
    """(1) r12-coupled projector == p=0 factorized build() on a p=0 basis."""
    from lih_projector_test import build as build_p0
    basis0 = build_basis_full(1, 1, ALPHA, p_set=(0,))
    G = C.grid_arrays(R, ALPHA, N_xi=40, N_eta=30, xi_max=20.0)
    new = C.projector_coupled(basis0, G, ALPHA)
    _, ref = build_p0(basis0, R, ALPHA, 1.0)               # lam=1 -> raw P
    rel = np.max(np.abs(new - ref) / np.maximum(np.abs(ref), 1e-10))
    assert rel < 5e-3, f"projector p=0 reduction failed: max rel {rel:.2e}"


def _vh_direct_quad(basis_p, R, alpha, zc=C.ZC_LI):
    """Independent direct 2e prolate quadrature of
    <g_i r12^pi | V_H(r1A)+V_H(r2A) | g_j r12^pj>, r12^P by explicit phi loop.
    Mirrors prolate_r12_mpf.vne_quad structure with V_H in place of the Coulomb sum."""
    from geovac.hylleraas import build_quadrature_grids
    g = build_quadrature_grids(N_xi=20, N_eta=14, N_phi=20, xi_max=16.0)
    xi, wxi, eta, weta = g['xi'], g['w_xi'], g['eta'], g['w_eta']
    dphi, wphi = g['dphi'], g['w_phi']
    hR = R / 2.0
    n = len(basis_p)
    V = np.zeros((n, n))
    for a in range(len(xi)):
        x1 = xi[a]
        for cc in range(len(xi)):
            x2 = xi[cc]
            ef = np.exp(-2.0 * alpha * (x1 + x2))
            for b in range(len(eta)):
                e1 = eta[b]; Jp1 = x1**2 - e1**2
                vh1 = C.V_H_closed(hR * (x1 + e1), zc)
                for d in range(len(eta)):
                    e2 = eta[d]; Jp2 = x2**2 - e2**2
                    vh2 = C.V_H_closed(hR * (x2 + e2), zc)
                    vh = vh1 + vh2
                    rho1 = np.sqrt(max((x1**2 - 1) * (1 - e1**2), 0.0))
                    rho2 = np.sqrt(max((x2**2 - 1) * (1 - e2**2), 0.0))
                    Aa = (x1 * e1 - x2 * e2)**2 + rho1**2 + rho2**2
                    Bc = 2 * rho1 * rho2
                    r12 = hR * np.sqrt(np.maximum(Aa - Bc * np.cos(dphi), 0.0))
                    wgt = wxi[a] * wxi[cc] * weta[b] * weta[d] * Jp1 * Jp2 * (hR**6) * ef * vh
                    for i in range(n):
                        bi, pi = basis_p[i]
                        gi = x1**bi.j * x2**bi.k * e1**bi.l * e2**bi.m
                        for jjj in range(n):
                            bj, pj = basis_p[jjj]
                            gj = x1**bj.j * x2**bj.k * e1**bj.l * e2**bj.m
                            P = pi + pj
                            rP = r12**P if P != 0 else np.ones_like(dphi)
                            V[i, jjj] += wgt * gi * gj * np.sum(rP * wphi) * 2 * np.pi
    return V


def test_vh_coupling_matches_direct_quadrature():
    """(2) The r12 coupling (Phi_1 elliptic kernel, P=1/2 dispatch) matches an
    INDEPENDENT direct phi-loop quadrature on a mixed p={0,1} basis."""
    specs = [(0, 0, 0, 0), (1, 0, 0, 0)]
    basis = ([(pr.ProductFn(j, l, k, mm, 0, ALPHA), 0) for (j, l, k, mm) in specs]
             + [(pr.ProductFn(j, l, k, mm, 0, ALPHA), 1) for (j, l, k, mm) in specs])
    # both use a matched grid so agreement is to grid accuracy, not method
    G = C.grid_arrays(R, ALPHA, N_xi=20, N_eta=14, xi_max=16.0)
    fast = C.vh_coupled(basis, G, ALPHA)
    ref = _vh_direct_quad(basis, R, ALPHA)
    # the P=1 (odd) block is rows 0-1 x cols 2-3 -- the r12-coupled part
    off = np.s_[:2, 2:]
    rel_off = np.max(np.abs(fast[off] - ref[off]) / np.maximum(np.abs(ref[off]), 1e-8))
    rel_all = np.max(np.abs(fast - ref) / np.maximum(np.abs(ref), 1e-8))
    assert rel_off < 2e-2, f"r12-coupled (P=1) V_H block mismatch: {rel_off:.2e}"
    assert rel_all < 2e-2, f"full mixed V_H mismatch: {rel_all:.2e}"


if __name__ == "__main__":
    test_vh_coupled_reduces_to_p0_reference()
    test_projector_coupled_reduces_to_p0_reference()
    test_vh_coupling_matches_direct_quadrature()
    print("all r12-coupled guards PASS")
