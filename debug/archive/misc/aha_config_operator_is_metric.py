"""The configuration operator IS the composed-basis metric  (2026-08-25, /aha follow-up).

CLAIM (Paper 32 rem:multicenter_composition / rem:config_berry_curvature follow-up).
The v5.1.0 configuration operator F = sum_i P_i -- whose conical intersections and pi Berry
phase over nuclear geometry were characterised as belonging to an "abstract configuration
operator, NOT the electronic Hamiltonian" -- is not abstract at all.  It is SIMILAR to the
composed basis's OVERLAP MATRIX G, whenever the intra-center blocks of G are orthonormal.

PROOF (three lines; the drivers already build both objects).
  G = Gram of the union basis, block (i,j) = <chi_i | chi_j>, diagonal blocks = I.
  Xh := chol(G)^T  (upper triangular, Xh^T Xh = G): columns = the basis vectors expressed in an
  orthonormal ambient frame.  Write W_k = Xh[:, block k].  Then
    (1)  W_k^T W_k = G[k,k] = I           => W_k has ORTHONORMAL columns
    (2)  P_k = W_k pinv(W_k) = W_k W_k^T  (the pinv collapses by (1))
    (3)  F = sum_k W_k W_k^T = Xh Xh^T = Xh (Xh^T Xh) Xh^{-1} = Xh G Xh^{-1}.
  So F ~ G:  SIMILAR, not merely cospectral.  Every spectral statement about F -- the conical
  intersections, the band gaps, the Berry phase of its eigenvectors -- is a statement about the
  composed basis's metric over nuclear geometry.

CONSEQUENCES TESTED HERE
  C1  spec(F) = spec(G) across the whole v5.1.0 arc: linear BeH2, bent BeH2 (incl. the H2O
      angle), and the four-center configuration.
  C2  Step (1) is LOAD-BEARING, not decoration: a deliberately non-orthonormal intra-center
      block breaks the identity.  (Non-tautology control.)
  C3  The four-center loop holonomy is a metric invariant in closed form:
      Tr(P1P2P3P4) = Tr(G12 G23 G34 G41)  -- exactly the S_ij form the v5.1.0 memo wrote.
  C4  The CI is a MID-spectrum metric degeneracy, not a conditioning blow-up: cond(G) is
      monotone through d* (the strong form of the /aha lead is FALSE, recorded here).
  C5  Orthogonalization consequence: around a CI the eigenvectors of the METRIC are
      double-valued (return overlap -1), so eigenvector-based (canonical) orthogonalization
      is path-dependent; symmetric Lowdin G^{-1/2}, a function of the matrix, is single-valued.
"""
from __future__ import annotations
import importlib.util
import json
import sys

import numpy as np
from scipy.linalg import fractional_matrix_power

sys.path.insert(0, "debug")


def _load(name, path):
    s = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


LIN = _load("lin", "debug/beh2_ci_exact_landscape.py")   # linear BeH2, sigma pair per center
BEND = _load("bend", "debug/beh2_bending_ci.py")         # in-plane {s,px,pz} per center
FOUR = _load("four", "debug/four_center_config.py")      # four planar centers

RESULTS = {}


# ---------------------------------------------------------------- shared helpers
def gram_linear(d1, d2):
    S1 = LIN._beh(d1)
    S2 = LIN.PAR @ LIN._beh(d2) @ LIN.PAR
    SH = LIN.PAR @ LIN._hh(d1 + d2) @ LIN.PAR
    I = np.eye(LIN.NS)
    return np.block([[I, S1, S2], [S1.T, I, SH], [S2.T, SH.T, I]])


def F_from_gram(G, block_sizes):
    """Build F = sum_k P_k from G by the drivers' own recipe (Cholesky whitening)."""
    Xh = np.linalg.cholesky(G).T
    F = np.zeros_like(G)
    o = 0
    for b in block_sizes:
        W = Xh[:, o:o + b]
        F += W @ np.linalg.pinv(W)
        o += b
    return F, Xh


def similarity_residual(G, block_sizes):
    """max|spec(F)-spec(G)| and the direct similarity residual |F Xh - Xh G|."""
    F, Xh = F_from_gram(G, block_sizes)
    dspec = float(np.abs(np.linalg.eigvalsh(F) - np.linalg.eigvalsh(G)).max())
    dsim = float(np.abs(F @ Xh - Xh @ G).max())
    return dspec, dsim, F, G


# ---------------------------------------------------------------- C1  spectrum identity
def c1_linear():
    rows = []
    for (d1, d2) in [(2.0, 2.0), (2.3, 2.3), (2.445, 2.445), (2.5, 2.5), (3.0, 3.0),
                     (2.698, 2.266), (4.0, 4.0), (6.0, 6.0)]:
        G = gram_linear(d1, d2)
        if np.linalg.eigvalsh(G).min() < 1e-9:
            continue
        dspec, dsim, F, _ = similarity_residual(G, [2, 2, 2])
        w = np.linalg.eigvalsh(G)
        gaps = np.diff(w)
        k = int(np.argmin(gaps))
        rows.append(dict(d1=d1, d2=d2, dspec=dspec, dsim=dsim,
                         cond=float(w.max() / w.min()), min_eig=float(w.min()),
                         min_gap=float(gaps[k]), cross_band=k, cross_val=float(w[k])))
    return rows


def c1_bending():
    rows = []
    for (d1, d2, theta) in [(2.5, 2.5, 180.0), (2.5, 2.5, 160.0), (2.5, 2.5, 120.0),
                            (2.5, 2.5, 104.5), (2.5, 2.5, 90.0), (3.0, 3.0, 104.5)]:
        a = np.deg2rad((180.0 - theta) / 2.0)
        G = BEND.metric(d1, d2, a)
        if np.linalg.eigvalsh(G).min() < 1e-9:
            continue
        dspec, dsim, F, _ = similarity_residual(G, [3, 3, 3])
        rows.append(dict(d1=d1, d2=d2, theta=theta, dspec=dspec, dsim=dsim,
                         min_eig=float(np.linalg.eigvalsh(G).min())))
    return rows


# ---------------------------------------------------------------- C2  non-tautology control
def c2_control():
    """Break intra-block orthonormality -> the identity must FAIL."""
    G = gram_linear(2.6, 2.6)
    dspec_ok, _, _, _ = similarity_residual(G, [2, 2, 2])
    Gb = G.copy()
    eps = 0.3                      # make center 0's two functions non-orthogonal
    Gb[0, 1] = Gb[1, 0] = eps
    dspec_bad, _, _, _ = similarity_residual(Gb, [2, 2, 2])
    return dict(intra_identity=dspec_ok, intra_broken=dspec_bad, eps=eps)


# ---------------------------------------------------------------- C3  four-center holonomy
def c3_four_center():
    pos = [(0.0, 0.0), (2.6, 0.0), (2.6, 2.6), (0.0, 2.6)]
    G, bs = FOUR.build_G(pos, rank=3)
    Ps = FOUR.projectors(G, bs)
    W_proj = float(np.real(np.trace(Ps[0] @ Ps[1] @ Ps[2] @ Ps[3])))
    o = np.cumsum([0] + list(bs))
    def blk(i, j):
        return G[o[i]:o[i + 1], o[j]:o[j + 1]]
    W_gram = float(np.trace(blk(0, 1) @ blk(1, 2) @ blk(2, 3) @ blk(3, 0)))
    dspec, dsim, _, _ = similarity_residual(G, list(bs))
    return dict(W_from_projectors=W_proj, W_from_gram_blocks=W_gram,
                residual=abs(W_proj - W_gram), dspec=dspec, dsim=dsim)


# ---------------------------------------------------------------- C5  orthogonalization loop
def c5_loop(cx, cy, r, N=400):
    """Transport the near-degenerate METRIC eigenvector once around a loop.

    Returns (eigvec_return, lowdin_return_drift, lowdin_loop_variation).
    Single-valuedness is a RETURN property: `lowdin_return_drift` compares G^{-1/2} at the
    end of the loop to its value at the start (same geometry).  `lowdin_loop_variation` is
    the max excursion *along* the loop -- it is nonzero for any non-constant loop and says
    nothing about multivaluedness; it is reported only to keep the two from being confused.
    """
    v0 = vp = None
    X0 = None
    variation = 0.0
    X_end = None
    for i in range(N + 1):
        th = 2 * np.pi * i / N
        G = gram_linear(cx + r * np.cos(th), cy + r * np.sin(th))
        w, V = np.linalg.eigh(G)
        X = fractional_matrix_power(G, -0.5).real
        if i == 0:
            k = int(np.argmin(np.diff(w)))
            vp = v0 = V[:, k].copy()
            X0 = X
        else:
            ov = V.T @ vp
            j = int(np.argmax(np.abs(ov)))
            vp = V[:, j] * np.sign(ov[j])
            variation = max(variation, float(np.abs(X - X0).max()))
        X_end = X
    return (float(vp @ v0), float(np.abs(X_end - X0).max()), variation)


# ---------------------------------------------------------------- run
if __name__ == "__main__":
    RESULTS["C1_linear"] = c1_linear()
    RESULTS["C1_bending"] = c1_bending()
    RESULTS["C2_control"] = c2_control()
    RESULTS["C3_four_center"] = c3_four_center()
    RESULTS["C5_loops"] = {
        "encircles_central_CI": c5_loop(2.445, 2.445, 0.05),
        "encircles_central_CI_tight": c5_loop(2.445, 2.445, 0.02),
        "CI_free_control": c5_loop(1.60, 1.60, 0.05),
    }

    print("=== C1  linear BeH2:  F ~ G  (similar, not just cospectral) ===")
    print(f"{'d1':>6} {'d2':>6} {'|dspec|':>10} {'|F Xh-Xh G|':>12} {'cond(G)':>9} "
          f"{'min_eig':>8} {'min_gap':>9}  cross@band(val)")
    for r in RESULTS["C1_linear"]:
        print(f"{r['d1']:6.3f} {r['d2']:6.3f} {r['dspec']:10.2e} {r['dsim']:12.2e} "
              f"{r['cond']:9.2f} {r['min_eig']:8.4f} {r['min_gap']:9.2e} "
              f"{r['cross_band']:>10}({r['cross_val']:.3f})")

    print("\n=== C1  bent BeH2 (in-plane 9-dim, incl. the H2O angle) ===")
    for r in RESULTS["C1_bending"]:
        print(f"  d={r['d1']:.2f} theta={r['theta']:6.1f}deg  |dspec|={r['dspec']:.2e} "
              f"|F Xh-Xh G|={r['dsim']:.2e}  min_eig={r['min_eig']:.4f}")

    print("\n=== C2  non-tautology control (break intra-block orthonormality) ===")
    c = RESULTS["C2_control"]
    print(f"   intra-block = I        -> |dspec| = {c['intra_identity']:.2e}   (identity holds)")
    print(f"   intra-block off by {c['eps']} -> |dspec| = {c['intra_broken']:.2e}   (identity FAILS)")

    print("\n=== C3  four centers: loop holonomy is a metric invariant ===")
    f = RESULTS["C3_four_center"]
    print(f"   Tr(P1P2P3P4)          = {f['W_from_projectors']:.12f}")
    print(f"   Tr(G12 G23 G34 G41)   = {f['W_from_gram_blocks']:.12f}   resid {f['residual']:.2e}")
    print(f"   F ~ G on 4 centers: |dspec|={f['dspec']:.2e}  |F Xh-Xh G|={f['dsim']:.2e}")

    print("\n=== C5  orthogonalization around the CI ===")
    print(f"   {'loop':28s} {'eigvec return':>14}  {'Lowdin RETURN':>14}  {'(loop variation)':>17}")
    for k, (sgn, ret, var) in RESULTS["C5_loops"].items():
        print(f"   {k:28s} {sgn:+14.4f}  {ret:14.2e}  {var:17.2e}")

    with open("debug/data/aha_config_operator_is_metric.json", "w") as fh:
        json.dump(RESULTS, fh, indent=2)
    print("\nwrote debug/data/aha_config_operator_is_metric.json")
