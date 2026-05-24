"""SU(4) rate constant analysis -- fourth datapoint for the unified GH-convergence theorem.

Empirical universality test for Conjecture A: c(G) = 4/pi for every compact Lie group.

Prior data (canonical / dual-Coxeter normalization):

  Group  rank |W|  A=4/pi   B=2|W|/pi^r   Extracted c    Sprint
  SU(2)  1    2    1.273    1.273         1.273 (analytical)  Paper 38
  SU(3)  2    6    1.273    1.216         1.243 (2.4%)        Q-rate
  Sp(2)  2    8    1.273    1.621         1.087 (14.6%)       Sp(2)/G_2
  G_2    2    12   1.273    2.432         1.177 (7.6%)        Sp(2)/G_2
  SU(4)  3    24   1.273    1.549         ???                 this sprint

SU(4) is rank 3: the simplest rank-3 group, simply-laced (A_3). The gap between
A (1.273) and B (1.549) is about 17.8% -- distinguishable at prior-sprint
precision (5-15%).

Method (mirrors debug/sp2_g2_rate_constant.py, extended to rank 3):
  1. Generic rank-3 central spectral Fejer kernel K_Lambda(g) on SU(4) = A_3.
  2. Mass-concentration moment gamma_Lambda via Weyl integration on T^3 = [0,2pi]^3.
  3. Panel of Lambda^2 values at adaptive saturated max_dynkin.
  4. Stein-Weiss 2-param fit gamma = (a log L + b)/L.
  5. Verdict A vs B, cross-group ratios.

Conventions
-----------
SU(4) = A_3, simply-laced, all roots length^2 = 2.
Cartan matrix:
    [[ 2, -1,  0],
     [-1,  2, -1],
     [ 0, -1,  2]]
Six positive roots, |W(A_3)| = 24 (= S_4), rho = (1,1,1) in omega basis.

Casimir rescale: C_gen(adjoint=(1,0,1)) = 8 in our short-root-squared-length-2
generic convention; canonical h^v_{SU(4)} = 4. Ratio = 8/4 = 2.
Same rescale as SU(3) (h^v_{SU(N)} = N, and generic Casimir / h^v = 2 for any A_n).
=>  a_can = a_gen / sqrt(2)

Geodesic distance: same dual-omega formalism. Coroot lattice in dual-omega
coords is Z^3 (since coroots of A_n pair with fundamental weights as
alpha_i^vee . omega_j = delta_{ij}, simply-laced).

Cross-references
----------------
- debug/dirac_triangle_extended_verify.py: build_A, SimpleLieAlgebra (rank-agnostic).
- debug/sp2_g2_rate_constant.py: rank-2 driver (reused infrastructure where possible).
- debug/su3_rate_constant.py: SU(3) canonical baseline.
- Paper 40, papers/group1_operator_algebras/paper_40_unified_propinquity_convergence.tex.
- debug/unified_gh_scoping_memo.md: scoping plan.

Output
------
- debug/data/su4_rate_constant.json: full panel + fits + verdict.

Does NOT modify production code in geovac/.
"""

from __future__ import annotations

import json
import math
import time
import sys
from pathlib import Path
from typing import List, Tuple

import numpy as np
import sympy as sp
from sympy import Rational, sqrt, Matrix

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Reuse infrastructure
from dirac_triangle_extended_verify import (
    SimpleLieAlgebra,
    build_A,
)
from sp2_g2_rate_constant import (
    _gl_nodes,
    precompute_weyl_orbit,
    fit_models,
)


# =============================================================================
# Generic rank-3 central spectral Fejer kernel + Weyl integration
# =============================================================================


def enumerate_irreps_under_casimir_a3(
    la: SimpleLieAlgebra, lambda_sq: float, max_dynkin: int = 12
) -> List[Tuple[int, int, int]]:
    """Enumerate all dominant weights (a, b, c) with 0 <= a,b,c <= max_dynkin and
    Casimir <= lambda_sq."""
    out = []
    for a in range(max_dynkin + 1):
        for b in range(max_dynkin + 1):
            for c in range(max_dynkin + 1):
                lam = (a, b, c)
                cas = float(la.casimir(lam))
                if cas <= lambda_sq + 1e-9:
                    out.append(lam)
    return out


def find_saturated_panel_a3(
    la: SimpleLieAlgebra, lambda_sq: float
) -> Tuple[List[Tuple[int, int, int]], int, int]:
    """Find the saturated irrep panel for given Lambda^2 by incrementing
    max_dynkin until the panel size stabilizes."""
    md = 8
    panel_prev = enumerate_irreps_under_casimir_a3(la, lambda_sq, max_dynkin=md)
    while True:
        md_next = md + 4
        panel_next = enumerate_irreps_under_casimir_a3(
            la, lambda_sq, max_dynkin=md_next
        )
        if len(panel_next) == len(panel_prev):
            return panel_prev, md, len(panel_prev)
        panel_prev = panel_next
        md = md_next
        if md > 60:
            return panel_prev, md, len(panel_prev)


def character_grid_a3(
    la: SimpleLieAlgebra,
    lam: Tuple[int, int, int],
    weyl_orbit_a3: list,
    T1: np.ndarray,
    T2: np.ndarray,
    T3: np.ndarray,
) -> np.ndarray:
    """Compute chi_lam(theta_1, theta_2, theta_3) on a 3D grid via the Weyl
    character formula.
    """
    rho = la.rho_omega
    lam_plus_rho = tuple(lam[i] + rho[i] for i in range(la.rank))

    z1 = np.exp(1j * T1)
    z2 = np.exp(1j * T2)
    z3 = np.exp(1j * T3)

    num = np.zeros_like(T1, dtype=np.complex128)
    den = np.zeros_like(T1, dtype=np.complex128)
    for W, sign in weyl_orbit_a3:
        Wlr = (
            int(W[0, 0]) * lam_plus_rho[0] + int(W[0, 1]) * lam_plus_rho[1] + int(W[0, 2]) * lam_plus_rho[2],
            int(W[1, 0]) * lam_plus_rho[0] + int(W[1, 1]) * lam_plus_rho[1] + int(W[1, 2]) * lam_plus_rho[2],
            int(W[2, 0]) * lam_plus_rho[0] + int(W[2, 1]) * lam_plus_rho[1] + int(W[2, 2]) * lam_plus_rho[2],
        )
        Wr = (
            int(W[0, 0]) * rho[0] + int(W[0, 1]) * rho[1] + int(W[0, 2]) * rho[2],
            int(W[1, 0]) * rho[0] + int(W[1, 1]) * rho[1] + int(W[1, 2]) * rho[2],
            int(W[2, 0]) * rho[0] + int(W[2, 1]) * rho[1] + int(W[2, 2]) * rho[2],
        )
        num += sign * z1 ** Wlr[0] * z2 ** Wlr[1] * z3 ** Wlr[2]
        den += sign * z1 ** Wr[0] * z2 ** Wr[1] * z3 ** Wr[2]

    with np.errstate(divide="ignore", invalid="ignore"):
        chi = num / den
    bad = ~np.isfinite(chi) | (np.abs(den) < 1e-12)
    if np.any(bad):
        dim_val = la.dim_weyl(lam)
        chi = np.where(bad, complex(float(dim_val), 0.0), chi)
    return chi


def weyl_denominator_sq_a3(
    la: SimpleLieAlgebra, T1: np.ndarray, T2: np.ndarray, T3: np.ndarray
) -> np.ndarray:
    """|Delta(theta)|^2 = prod_{alpha > 0} 4 sin^2(<alpha, theta>/2)
    in 3D dual-omega coords."""
    out = np.ones_like(T1)
    for alpha in la.positive_roots_omega:
        a1, a2, a3v = int(alpha[0]), int(alpha[1]), int(alpha[2])
        phase = (a1 * T1 + a2 * T2 + a3v * T3) / 2.0
        out = out * 4.0 * np.sin(phase) ** 2
    return out


def geodesic_dist_sq_grid_a3(
    la: SimpleLieAlgebra,
    T1: np.ndarray,
    T2: np.ndarray,
    T3: np.ndarray,
    n_search: int = 1,
) -> np.ndarray:
    """Bi-invariant Killing-form geodesic distance squared from identity to
    theta in dual-omega coords, rank 3.

    Coroot lattice in dual-omega = Z^3. Minimize over (n_1, n_2, n_3) in
    [-n_search, n_search]^3.
    """
    G_inv = la.G_omega.inv()
    # Symmetric 3x3 matrix
    G11 = float(G_inv[0, 0])
    G22 = float(G_inv[1, 1])
    G33 = float(G_inv[2, 2])
    G12 = float(G_inv[0, 1])
    G13 = float(G_inv[0, 2])
    G23 = float(G_inv[1, 2])

    two_pi = 2.0 * math.pi
    best = None
    for n1 in range(-n_search, n_search + 1):
        for n2 in range(-n_search, n_search + 1):
            for n3 in range(-n_search, n_search + 1):
                t1p = T1 - n1 * two_pi
                t2p = T2 - n2 * two_pi
                t3p = T3 - n3 * two_pi
                dsq = (
                    G11 * t1p * t1p + G22 * t2p * t2p + G33 * t3p * t3p
                    + 2.0 * G12 * t1p * t2p
                    + 2.0 * G13 * t1p * t3p
                    + 2.0 * G23 * t2p * t3p
                )
                if best is None:
                    best = dsq
                else:
                    best = np.minimum(best, dsq)
    return best


def verify_haar_normalization_a3(
    la: SimpleLieAlgebra, n_quad: int = 40
) -> Tuple[float, int]:
    """Verify (1/(|W| (2pi)^3)) int_T^3 |Delta|^2 dt = 1 via Gauss-Legendre on
    [0, 2pi]^3.
    """
    weyl_orbit = precompute_weyl_orbit(la)
    W_order = len(weyl_orbit)
    two_pi = 2.0 * math.pi
    nodes, weights = _gl_nodes(0.0, two_pi, n_quad)
    T1, T2, T3 = np.meshgrid(nodes, nodes, nodes, indexing="ij")
    W1, W2, W3 = np.meshgrid(weights, weights, weights, indexing="ij")
    measure = weyl_denominator_sq_a3(la, T1, T2, T3)
    integral = np.sum(measure * W1 * W2 * W3) / (W_order * two_pi ** 3)
    return float(integral), W_order


def gamma_lambda_a3(
    la: SimpleLieAlgebra,
    irreps: List[Tuple[int, int, int]],
    weyl_orbit: list,
    n_quad: int = 40,
) -> float:
    """Compute gamma_Lambda = int_G K_Lambda(g) d(e, g) dg on SU(4) via 3D
    Gauss-Legendre on the maximal torus T^3 = [0, 2pi]^3.

    K_Lambda(theta) = (1/Z) |sum_pi sqrt(dim pi) chi_pi(theta)|^2,
    Z = sum_pi dim pi.

    Optimized for rank 3:
      1. Precompute denominator den(theta) ONCE (depends only on rho).
      2. Precompute powers of z1, z2, z3 up to max_pqr+3 ONCE (avoid repeated
         numpy power calls per irrep).
      3. Loop over Weyl elements outside the irrep loop for the numerator.
    """
    W_order = len(weyl_orbit)
    two_pi = 2.0 * math.pi
    nodes, weights = _gl_nodes(0.0, two_pi, n_quad)
    T1, T2, T3 = np.meshgrid(nodes, nodes, nodes, indexing="ij")
    W1, W2, W3 = np.meshgrid(weights, weights, weights, indexing="ij")

    # Geodesic distance grid (sqrt of d^2_min).
    dsq = geodesic_dist_sq_grid_a3(la, T1, T2, T3)
    d_grid = np.sqrt(np.maximum(dsq, 0.0))

    # Weyl denominator |Delta|^2.
    measure = weyl_denominator_sq_a3(la, T1, T2, T3)

    z1 = np.exp(1j * T1)
    z2 = np.exp(1j * T2)
    z3 = np.exp(1j * T3)

    rho = la.rho_omega
    # Precompute Weyl-shifted rho indices.
    Wr_list = []
    for W, sign in weyl_orbit:
        Wr = (
            int(W[0, 0]) * rho[0] + int(W[0, 1]) * rho[1] + int(W[0, 2]) * rho[2],
            int(W[1, 0]) * rho[0] + int(W[1, 1]) * rho[1] + int(W[1, 2]) * rho[2],
            int(W[2, 0]) * rho[0] + int(W[2, 1]) * rho[1] + int(W[2, 2]) * rho[2],
        )
        Wr_list.append((Wr, sign))

    # Determine max power needed for z1, z2, z3 across all irreps and all
    # Weyl elements. Wlr[i] = sum_j W[i,j] * (lam+rho)[j]. With max |W[i,j]|
    # bounded (in W(A_3) it's 0,1,-1), this is bounded by 3*(max_dynkin+1).
    max_dynkin_plus_rho = max(max(lam) for lam in irreps) + 1
    # max abs entry of W * (lam+rho) bounded by rank * max(W)*max(lam+rho)
    # For A_3 Weyl matrices, entries are in {-1, 0, 1, 2, -2, etc.}; safe bound:
    max_idx = 0
    for W, _ in weyl_orbit:
        for i in range(3):
            row_abs_sum = sum(abs(int(W[i, j])) for j in range(3))
            if row_abs_sum * max_dynkin_plus_rho > max_idx:
                max_idx = row_abs_sum * max_dynkin_plus_rho

    # Precompute z1, z2, z3 to powers -max_idx..max_idx.
    z1_pow = {0: np.ones_like(z1)}
    z2_pow = {0: np.ones_like(z2)}
    z3_pow = {0: np.ones_like(z3)}
    for k in range(1, max_idx + 1):
        z1_pow[k] = z1_pow[k - 1] * z1
        z2_pow[k] = z2_pow[k - 1] * z2
        z3_pow[k] = z3_pow[k - 1] * z3
        z1_pow[-k] = np.conj(z1_pow[k])
        z2_pow[-k] = np.conj(z2_pow[k])
        z3_pow[-k] = np.conj(z3_pow[k])

    # Precompute denominator den(theta) = sum_w sign(w) z^{w rho}.
    den = np.zeros_like(z1)
    for Wr, sign in Wr_list:
        den = den + sign * z1_pow[Wr[0]] * z2_pow[Wr[1]] * z3_pow[Wr[2]]
    # Avoid division by zero at Weyl walls (codim-1, measure zero).
    # We use a small eps replacement at bad points after we form chi.

    # Build kernel sum: sum_pi sqrt(dim_pi) chi_pi(theta) = (1/den) * sum_pi sqrt(dim_pi) N_pi(theta).
    # N_pi(theta) = sum_w sign(w) z^{w(lam+rho)}. Loop over irreps.
    Z = sum(int(la.dim_weyl(lam)) for lam in irreps)
    sum_grid = np.zeros_like(z1)
    for lam in irreps:
        dim_lam = int(la.dim_weyl(lam))
        lam_plus_rho = (lam[0] + rho[0], lam[1] + rho[1], lam[2] + rho[2])
        N = np.zeros_like(z1)
        for W, sign in weyl_orbit:
            Wlr0 = int(W[0, 0]) * lam_plus_rho[0] + int(W[0, 1]) * lam_plus_rho[1] + int(W[0, 2]) * lam_plus_rho[2]
            Wlr1 = int(W[1, 0]) * lam_plus_rho[0] + int(W[1, 1]) * lam_plus_rho[1] + int(W[1, 2]) * lam_plus_rho[2]
            Wlr2 = int(W[2, 0]) * lam_plus_rho[0] + int(W[2, 1]) * lam_plus_rho[1] + int(W[2, 2]) * lam_plus_rho[2]
            N = N + sign * z1_pow[Wlr0] * z2_pow[Wlr1] * z3_pow[Wlr2]
        sum_grid = sum_grid + math.sqrt(dim_lam) * N

    # chi_pi = N_pi / den. So |sum sqrt(dim) chi|^2 = |sum sqrt(dim) N|^2 / |den|^2.
    # But |den|^2 IS the Weyl denominator squared = measure. So:
    #   K(theta) * measure = (|sum sqrt(dim) N|^2 / |den|^2) * |den|^2 = |sum sqrt(dim) N|^2.
    # Big simplification: K * |Delta|^2 = |sum_pi sqrt(dim_pi) N_pi|^2 / Z, with NO division.
    K_times_measure = (np.abs(sum_grid) ** 2) / Z

    integrand = K_times_measure * d_grid * W1 * W2 * W3
    integral = np.sum(integrand) / (W_order * two_pi ** 3)
    return float(integral.real)


# =============================================================================
# Panel runner
# =============================================================================


def run_panel_a3(
    la: SimpleLieAlgebra,
    name: str,
    lambda_sq_panel: List[int],
) -> Tuple[List[dict], dict]:
    """Run gamma_Lambda computation across panel of Lambda^2 values for rank 3."""
    print(f"\n{'=' * 78}")
    print(f"{name}: gamma_Lambda panel (rank 3)")
    print(f"{'=' * 78}")

    weyl_orbit = precompute_weyl_orbit(la)
    W_order = len(weyl_orbit)
    print(f"  |W| = {W_order}")
    print(f"  rank = {la.rank}")
    print(f"  positive roots ({len(la.positive_roots_omega)}): {la.positive_roots_omega}")
    print(f"  rho_omega = {la.rho_omega}")
    print(f"  G_omega = {la.G_omega.tolist()}")
    print(f"  G_omega^-1 = {la.G_omega.inv().tolist()}")

    # Verify Haar normalization
    haar, _ = verify_haar_normalization_a3(la, n_quad=40)
    print(f"  Haar @ n_quad=40: {haar:.10f}  (expect 1.0)")
    if abs(haar - 1.0) > 1e-3:
        print(f"  WARNING: Haar normalization off by {haar - 1.0:.3e}")

    rows = []
    print(
        f"\n  {'Lam^2':>6} {'Lambda':>8} {'max_d':>5} {'n_irr':>5} {'max_pqr':>7} "
        f"{'n_quad':>6} {'gamma':>10} {'c_est':>8} {'t(s)':>7}"
    )
    for lsq in lambda_sq_panel:
        irreps, md, n_irr = find_saturated_panel_a3(la, lsq)
        max_pqr = max(a + b + c for (a, b, c) in irreps) if irreps else 0
        # n_quad: scale with max_pqr. For rank-3 character with phases
        # exp(i * k * theta) up to k ~ max_pqr + a few, need n_quad >=
        # ~2 * (max_pqr + 3). Cap at 80 to bound memory.
        nq = max(40, 2 * max_pqr + 30)
        # Memory bound: n_quad^3 * 16 bytes for complex128 ~= 512 MB at n=80
        nq = min(nq, 80)
        t0 = time.time()
        g = gamma_lambda_a3(la, irreps, weyl_orbit, n_quad=nq)
        dt = time.time() - t0
        L = math.sqrt(lsq)
        c_est = L * g / math.log(L) if L > 1.5 else float("nan")
        rows.append({
            "lambda_sq": lsq,
            "lambda": L,
            "gamma": g,
            "n_irreps": n_irr,
            "max_dynkin_for_saturation": md,
            "max_p_plus_q_plus_r": max_pqr,
            "n_quad": nq,
            "time_seconds": dt,
            "L*gamma/log(L)": c_est,
        })
        print(
            f"  {lsq:>6} {L:>8.3f} {md:>5} {n_irr:>5} {max_pqr:>7} "
            f"{nq:>6} {g:>10.6f} {c_est:>8.4f} {dt:>7.1f}"
        )

    meta = {
        "weyl_group_order": W_order,
        "haar_check": haar,
        "rank": la.rank,
        "n_positive_roots": len(la.positive_roots_omega),
    }
    return rows, meta


def verdict_for_group_a3(
    name: str,
    a_extracted_gen: float,
    weyl_order: int,
    rank: int,
    casimir_rescale: float,
) -> dict:
    """Compare extracted a against A (4/pi) and B (2|W|/pi^r) under
    dual-Coxeter normalization."""
    pi = math.pi
    cA = 4.0 / pi
    cB = 2.0 * weyl_order / (pi ** rank)

    a_canonical = a_extracted_gen / math.sqrt(casimir_rescale)

    err_A = abs(a_canonical - cA) / cA * 100
    err_B = abs(a_canonical - cB) / cB * 100

    gap_AB_pct = abs(cA - cB) / max(cA, cB) * 100

    if err_A < 5.0 and err_B > 15.0:
        verdict = f"A ({name} a_can={a_canonical:.4f} matches 4/pi={cA:.4f} at {err_A:.2f}%)"
    elif err_B < 5.0 and err_A > 15.0:
        verdict = f"B ({name} a_can={a_canonical:.4f} matches 2|W|/pi^r={cB:.4f} at {err_B:.2f}%)"
    elif err_A < 5.0 and err_B < 5.0:
        verdict = f"AMBIGUOUS ({name} both within 5%: A err={err_A:.2f}%, B err={err_B:.2f}%)"
    elif err_A < err_B:
        verdict = f"A-leaning ({name} a_can={a_canonical:.4f}, A err={err_A:.2f}%, B err={err_B:.2f}%)"
    elif err_B < err_A:
        verdict = f"B-leaning ({name} a_can={a_canonical:.4f}, B err={err_B:.2f}%, A err={err_A:.2f}%)"
    else:
        verdict = "FAIL"

    return {
        "name": name,
        "a_gen_extracted": a_extracted_gen,
        "casimir_rescale": casimir_rescale,
        "a_canonical": a_canonical,
        "weyl_order": weyl_order,
        "rank": rank,
        "c_A_4_over_pi": cA,
        "c_B_2W_over_pi_r": cB,
        "err_A_pct": err_A,
        "err_B_pct": err_B,
        "gap_AB_pct": gap_AB_pct,
        "verdict": verdict,
    }


def main():
    out_dir = HERE / "data"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("SU(4) rate constant -- fourth datapoint for the unified GH theorem")
    print("=" * 78)
    print()
    print("  Prior data (dual-Coxeter normalization):")
    print(f"    SU(2): a_can = 4/pi = {4/math.pi:.4f}  (analytical, Paper 38)")
    print(f"    SU(3): a_can = 1.243 (2.4% from A=4/pi, 2.2% from B=12/pi^2)")
    print(f"    Sp(2): a_can = 1.087 (14.6% from A, 32.9% from B; A favored)")
    print(f"    G_2:   a_can = 1.177 (7.6% from A, 51.6% from B; A favored)")
    print()
    print(f"  Conjecture A:           c = 4/pi    = {4/math.pi:.6f}")
    print(f"  Conjecture B (SU(4)):   c = 48/pi^3 = {48/math.pi**3:.6f}")
    gap_pct = abs(4/math.pi - 48/math.pi**3) / max(4/math.pi, 48/math.pi**3) * 100
    print(f"  Gap A vs B (SU(4)):     {gap_pct:.2f}%")

    results = {
        "metadata": {
            "purpose": "Fourth datapoint (rank 3) for universal rate constant test",
            "predictions": {
                "A_4_over_pi": 4.0 / math.pi,
                "B_48_over_pi3": 48.0 / math.pi ** 3,
                "gap_AB_pct": gap_pct,
            },
            "casimir_convention": "generic (short root squared length 2); "
                "dual-Coxeter rescale = C_gen(adj)/h^v = 8/4 = 2",
            "lambda_can_to_gen": "L_can = L_gen / sqrt(2), L^2_can = L^2_gen / 2",
        }
    }

    # =========================================================================
    # SU(4) = A_3 panel
    # =========================================================================
    a3 = build_A(3)

    # Sanity-check Casimir rescale numerically
    adj = (1, 0, 1)
    cas_adj = float(a3.casimir(adj))
    print(f"\n  SU(4) sanity:")
    print(f"    Adjoint (1,0,1) dim = {a3.dim_weyl(adj)} (expect 15)")
    print(f"    Adjoint Casimir (generic) = {cas_adj} (expect 8)")
    print(f"    h^v_{{SU(4)}} = 4")
    print(f"    Rescale = {cas_adj}/4 = {cas_adj/4}")

    # SU(4) generic Lambda^2 panel. Canonical Lambda^2 = generic Lambda^2 / 2.
    # Rank-3 cost grows aggressively with max_pqr; cap at L^2_gen=600 (L_can ~
    # 17.3) which clears the asymptotic threshold L_can >= 7 (L^2_can >= 49,
    # L^2_gen >= 98) with margin. Mirrors what Sp(2)/G_2 used for asymptotic
    # extraction while keeping compute manageable. Per-row cost expected:
    # L^2=100 ~30s, L^2=200 ~3min, L^2=400 ~15min, L^2=600 ~30min => total <2h.
    panel = [8, 12, 16, 20, 28, 36, 48, 60, 80, 100, 140, 200, 300, 400, 600]
    rows_a3, meta_a3 = run_panel_a3(a3, "SU(4)/A_3", panel)
    results["su4"] = {"rows": rows_a3, "meta": meta_a3}

    lams_a3 = [r["lambda"] for r in rows_a3]
    gs_a3 = [r["gamma"] for r in rows_a3]

    # Fits at several exclusion thresholds. SU(4) canonical Lambda = generic / sqrt(2).
    # Asymptotic regime per Q-rate convention: L_can >= 7 -> L_gen >= 7*sqrt(2) ~ 9.9
    fits_all = {}
    # Asymptotic threshold L_can >= 7 -> L_gen >= 9.9. With panel max L_gen ~
    # sqrt(600) ~= 24.5 (L_can ~= 17.3), several exclusion thresholds available.
    for lo in [0.0, 5.0, 9.9, 14.0, 17.3, 20.0]:
        fits_all[f"min_L_{lo}"] = fit_models(
            lams_a3, gs_a3, exclude_below_L=lo, label=f"SU(4) L>={lo}"
        )
    print("\n=== SU(4) fits (L_gen >= X corresponds to L_can >= X/sqrt(2)) ===")
    for key, f in fits_all.items():
        if "two_param" not in f:
            continue
        a = f["two_param"]["a"]
        a_can = a / math.sqrt(2)
        print(f"  {key}: n_data={f['n_data']}, a_gen = {a:.6f}, a_can = {a_can:.6f}")
    results["su4"]["fits"] = fits_all

    # Primary SU(4) estimate: 2-param Stein-Weiss on L_gen >= 9.9 (L_can >= 7,
    # asymptotic regime mirroring su3_rate_constant.py convention)
    primary_key = "min_L_9.9"
    a_a3_gen = fits_all[primary_key]["two_param"]["a"]
    rescale = 2.0  # = C_gen(adjoint)/h^v = 8/4
    a_a3_can = a_a3_gen / math.sqrt(rescale)
    print(f"\n  SU(4) primary estimate (L_gen >= 9.9, L_can >= 7):")
    print(f"    a_gen = {a_a3_gen:.6f}")
    print(f"    a_can = a_gen / sqrt({rescale}) = {a_a3_can:.6f}")
    print(f"    Conjecture A: 4/pi    = {4/math.pi:.6f}")
    print(f"    Conjecture B: 48/pi^3 = {48/math.pi**3:.6f}")

    verd = verdict_for_group_a3(
        "SU(4)", a_a3_gen, weyl_order=24, rank=3, casimir_rescale=rescale
    )
    print(f"    VERDICT: {verd['verdict']}")
    results["su4"]["verdict"] = verd

    # =========================================================================
    # Cross-group ratio table
    # =========================================================================
    # Prior values (a_canonical, dual-Coxeter):
    a_su2_can = 4.0 / math.pi  # analytical
    a_su3_can = 1.243144  # canonical SU(3), from su3_rate_constant.py
    a_sp2_can = 1.0875    # primary Sp(2)
    a_g2_can = 1.1765     # primary G_2
    a_su4_can = a_a3_can

    print("\n" + "=" * 78)
    print("CROSS-GROUP RATIOS (convention-independent)")
    print("=" * 78)
    print(f"  SU(4) / SU(2) = {a_su4_can/a_su2_can:.4f}  (A predicts 1.000)")
    print(f"  SU(4) / SU(3) = {a_su4_can/a_su3_can:.4f}  (A predicts 1.000, "
          f"B predicts (48/pi^3)/(12/pi^2) = 4/pi = {4/math.pi:.4f})")
    print(f"  SU(4) / Sp(2) = {a_su4_can/a_sp2_can:.4f}  (A predicts 1.000)")
    print(f"  SU(4) / G_2   = {a_su4_can/a_g2_can:.4f}  (A predicts 1.000)")

    # B prediction ratios:
    # SU(4) / SU(2) = (48/pi^3) / (4/pi) = 12/pi^2 = 1.216 in B
    # SU(4) / SU(3) = (48/pi^3) / (12/pi^2) = 4/pi = 1.273 in B
    # SU(4) / Sp(2) = (48/pi^3) / (16/pi^2) = 3/pi = 0.955 in B
    # SU(4) / G_2   = (48/pi^3) / (24/pi^2) = 2/pi = 0.637 in B
    print()
    print("  B-prediction ratios (for reference):")
    print(f"    SU(4)/SU(2)_B = 12/pi^2 = {12/math.pi**2:.4f}")
    print(f"    SU(4)/SU(3)_B = 4/pi    = {4/math.pi:.4f}")
    print(f"    SU(4)/Sp(2)_B = 3/pi    = {3/math.pi:.4f}")
    print(f"    SU(4)/G_2_B   = 2/pi    = {2/math.pi:.4f}")

    results["cross_group_ratios"] = {
        "a_canonical_values": {
            "SU(2)": a_su2_can, "SU(3)": a_su3_can, "Sp(2)": a_sp2_can,
            "G_2": a_g2_can, "SU(4)": a_su4_can,
        },
        "ratios_SU4_over": {
            "SU(2)": a_su4_can / a_su2_can,
            "SU(3)": a_su4_can / a_su3_can,
            "Sp(2)": a_su4_can / a_sp2_can,
            "G_2":   a_su4_can / a_g2_can,
        },
        "A_predicts": {"all": 1.0},
        "B_predicts_SU4_over": {
            "SU(2)": 12 / math.pi ** 2,
            "SU(3)": 4 / math.pi,
            "Sp(2)": 3 / math.pi,
            "G_2":   2 / math.pi,
        },
    }

    # =========================================================================
    # Combined verdict
    # =========================================================================
    print("\n" + "=" * 78)
    print("COMBINED VERDICT (5 datapoints)")
    print("=" * 78)
    print(f"  Conjecture A: c(G) = 4/pi = {4/math.pi:.4f}  (rank-invariant)")
    print()
    print(f"  {'Group':>6} {'rank':>4} {'|W|':>4} {'A':>8} {'B':>8} {'extracted':>10} {'A err':>7} {'B err':>7}")
    rows_summary = [
        ("SU(2)", 1, 2, 4/math.pi, 4/math.pi, a_su2_can),
        ("SU(3)", 2, 6, 4/math.pi, 12/math.pi**2, a_su3_can),
        ("Sp(2)", 2, 8, 4/math.pi, 16/math.pi**2, a_sp2_can),
        ("G_2",   2, 12, 4/math.pi, 24/math.pi**2, a_g2_can),
        ("SU(4)", 3, 24, 4/math.pi, 48/math.pi**3, a_su4_can),
    ]
    for name, r, W, cA, cB, a_can in rows_summary:
        eA = abs(a_can - cA) / cA * 100
        eB = abs(a_can - cB) / cB * 100
        print(f"  {name:>6} {r:>4} {W:>4} {cA:>8.4f} {cB:>8.4f} {a_can:>10.4f} {eA:>6.1f}% {eB:>6.1f}%")

    results["combined_summary"] = [
        {"group": name, "rank": r, "weyl_order": W,
         "c_A": cA, "c_B": cB, "a_canonical_extracted": a_can,
         "err_A_pct": abs(a_can - cA) / cA * 100,
         "err_B_pct": abs(a_can - cB) / cB * 100}
        for name, r, W, cA, cB, a_can in rows_summary
    ]

    out_path = out_dir / "su4_rate_constant.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")

    return results


if __name__ == "__main__":
    main()
