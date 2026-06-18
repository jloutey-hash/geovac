"""Sp(2) and G_2 rate constant analysis for the unified GH-convergence theorem.

Discriminator test between two competing conjectures for the leading constant
c(G) in gamma_Lambda(G) ~ c(G) * log(Lambda) / Lambda on a compact Lie group G:

  Conjecture A (rank-invariant): c(G) = 4/pi for every compact Lie group.
  Conjecture B (Weyl-formula):   c(G) = 2 |W(G)| / pi^rank(G).

Both agree at SU(2) (|W|=2, rank=1): 4/pi. They differ at higher rank.
SU(3) (sprint Q-rate, debug/su3_rate_constant.py) gave c ~ 1.243, which sits
between 4/pi ~ 1.273 (A) and 12/pi^2 ~ 1.216 (B) within 2.3% of each --
SU(3) cannot distinguish them. Sp(2) (|W|=8) and G_2 (|W|=12), both rank 2,
give categorically different B predictions:

   Group    rank  |W|   A=4/pi    B=2|W|/pi^r    Gap
   Sp(2)    2     8     1.273     16/pi^2=1.621  27%
   G_2      2     12    1.273     24/pi^2=2.432  91%

A 5% precision extraction at either group cleanly distinguishes A from B.

Method (mirrors debug/su3_rate_constant.py):
  Step 1: Central spectral Fejer kernel K_Lambda(g) on G.
  Step 2: Mass-concentration moment gamma_Lambda = int_G K_Lambda(g) d(e,g)^2 dg
          via Weyl integration on the maximal torus T^2.
  Step 3: Panel of Lambda^2 values [4..1000] at adaptive saturated max_dynkin.
  Step 4: Stein-Weiss 2-param fit gamma = (a log L + b)/L.
  Step 5: Verdict A vs B.

Conventions
-----------
For both groups we use the *Dual-omega torus parameterization*:
  (theta_1, theta_2) in [0, 2 pi]^2,
where theta_i is dual to fundamental weight omega_i, i.e. for any weight
lambda with Dynkin labels (lambda_1, lambda_2) in omega basis:
    <lambda, theta> = lambda_1 theta_1 + lambda_2 theta_2.

A root alpha with omega-basis coords (a_1, a_2) acts as <alpha, theta> =
a_1 theta_1 + a_2 theta_2. Weyl denominator:
    |Delta(theta)|^2 = prod_{alpha > 0} 4 sin^2((a_1 theta_1 + a_2 theta_2)/2).

The Weyl integration formula:
    int_G f dg = (1 / (|W| (2 pi)^2)) int_{T^2} f(theta) |Delta(theta)|^2 dtheta.
This is verified numerically (Haar of constant 1 -> 1).

For the geodesic distance squared, the metric in dual-omega coords is the
inverse Gram matrix G_omega^{-1} on omega basis:
    d^2(0, theta) = theta^T G_omega^{-1} theta.
The lattice of identity elements is the *coroot lattice* (because dual-omega
coords are dual to fundamental weights, and the kernel of exp on the maximal
torus of a simply-connected compact Lie group is 2 pi Q^vee, which in
dual-omega coords is exactly Z^2 -- shifts by 2 pi in each coordinate).
So d^2_min = min over (n_1, n_2) in [-1, 1]^2 of
       (theta - 2 pi (n_1, n_2))^T G_omega^{-1} (theta - 2 pi (n_1, n_2)).

Cross-references
----------------
- debug/dirac_triangle_extended_verify.py: build_C2, build_G2, tensor_product
  infrastructure (Cartan + Freudenthal + Brauer-Klimyk for any simple Lie algebra).
- debug/su3_rate_constant.py: SU(3) blueprint for the rate constant analysis.
- debug/unified_gh_scoping_memo.md: forward plan and verdict structure.
- geovac/central_fejer_su2.py: SU(2) Paper 38 reference.

Output
------
- debug/data/sp2_g2_rate_constant.json: full panel + fits + verdicts.

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
import mpmath
import sympy as sp
from sympy import Rational, sqrt, Matrix

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Reuse the rank-2 group machinery
from dirac_triangle_extended_verify import (
    SimpleLieAlgebra,
    build_A,
    build_C2,
    build_G2,
    tensor_product,
)


# =============================================================================
# Generic rank-2 central spectral Fejer kernel + Weyl integration
# =============================================================================


def enumerate_irreps_under_casimir(
    la: SimpleLieAlgebra, lambda_sq: float, max_dynkin: int = 20
) -> List[Tuple[int, int]]:
    """Enumerate all dominant weights (a, b) with 0 <= a, b <= max_dynkin and
    casimir <= lambda_sq."""
    out = []
    for a in range(max_dynkin + 1):
        for b in range(max_dynkin + 1):
            lam = (a, b)
            c = float(la.casimir(lam))
            if c <= lambda_sq + 1e-9:
                out.append(lam)
    return out


def find_saturated_panel(
    la: SimpleLieAlgebra, lambda_sq: float
) -> Tuple[List[Tuple[int, int]], int, int]:
    """Find the saturated irrep panel for given Lambda^2 by incrementing
    max_dynkin until the panel size stabilizes."""
    md = 20
    panel_prev = enumerate_irreps_under_casimir(la, lambda_sq, max_dynkin=md)
    while True:
        md_next = md + 10
        panel_next = enumerate_irreps_under_casimir(
            la, lambda_sq, max_dynkin=md_next
        )
        if len(panel_next) == len(panel_prev):
            return panel_prev, md, len(panel_prev)
        panel_prev = panel_next
        md = md_next
        if md > 200:
            return panel_prev, md, len(panel_prev)


def _gl_nodes(a: float, b: float, n: int):
    """Gauss-Legendre nodes/weights on [a, b]."""
    nodes, weights = np.polynomial.legendre.leggauss(n)
    return (
        (b - a) / 2 * nodes + (b + a) / 2,
        (b - a) / 2 * weights,
    )


def precompute_weyl_orbit(la: SimpleLieAlgebra) -> List[Tuple[Tuple[int, int], int]]:
    """Enumerate Weyl group as list of (matrix, sign).

    For rank 2, the Weyl group is generated by two simple reflections s_1, s_2.
    In dual-omega coords, s_i: theta -> theta - (alpha_i, theta) * (alpha_i in
    dual basis). For simple weight reflections via Cartan matrix transpose:
    s_i acts on a weight mu in omega basis as
        s_i(mu)_j = mu_j - mu_i A_{ji}
    where A is the Cartan matrix (so column i of A gives the simple root
    alpha_i in omega basis). For action on dual-omega theta coordinates,
    s_i: theta_j -> theta_j - theta_i A_{ji} (same formula since the action
    is by the Cartan matrix transpose on the dual basis).

    Returns list of (W_matrix_as_2x2_tuple, sign).
    For |W(A_2)|=6, |W(C_2)|=8, |W(G_2)|=12.
    """
    rank = la.rank
    cartan = la.cartan

    def s_i_action(i: int):
        """Return 2x2 matrix of s_i acting on dual-omega theta vector."""
        M = sp.eye(rank)
        # s_i: theta_j -> theta_j - theta_i A_{ji}
        # Matrix form: M[j, i] -= A[j, i] for the i-th column adjustment
        for j in range(rank):
            M[j, i] = M[j, i] - cartan[j, i]
        return M

    s_mats = [s_i_action(i) for i in range(rank)]

    # BFS Weyl orbit: start with identity, apply s_i, collect.
    I = sp.eye(rank)
    orbit = {tuple(tuple(row) for row in I.tolist()): (I, 1)}
    frontier = [(I, 1)]
    while frontier:
        new_frontier = []
        for W, sign in frontier:
            for s in s_mats:
                W_new = s * W
                key = tuple(tuple(int(W_new[i, j]) for j in range(rank)) for i in range(rank))
                if key not in orbit:
                    orbit[key] = (W_new, -sign)
                    new_frontier.append((W_new, -sign))
        frontier = new_frontier

    return [(W, sign) for W, sign in orbit.values()]


def character_grid(
    la: SimpleLieAlgebra,
    lam: Tuple[int, int],
    weyl_orbit: list,
    T1: np.ndarray,
    T2: np.ndarray,
) -> np.ndarray:
    """Compute chi_lam(theta_1, theta_2) on a 2D grid via the Weyl character
    formula:

        chi_lam(theta) = N(lam) / N(0)
        N(mu) = sum_w sign(w) exp(i <w(mu+rho), theta>)
              = sum_w sign(w) prod_i z_i^{(w(mu+rho))_i}

    where z_i = exp(i theta_i).

    Returns complex numpy array same shape as T1.
    """
    rho = la.rho_omega
    lam_plus_rho = tuple(lam[i] + rho[i] for i in range(la.rank))

    z1 = np.exp(1j * T1)
    z2 = np.exp(1j * T2)

    num = np.zeros_like(T1, dtype=np.complex128)
    den = np.zeros_like(T1, dtype=np.complex128)
    for W, sign in weyl_orbit:
        # W(lam + rho) in omega basis:
        Wlr = (
            int(W[0, 0]) * lam_plus_rho[0] + int(W[0, 1]) * lam_plus_rho[1],
            int(W[1, 0]) * lam_plus_rho[0] + int(W[1, 1]) * lam_plus_rho[1],
        )
        Wr = (
            int(W[0, 0]) * rho[0] + int(W[0, 1]) * rho[1],
            int(W[1, 0]) * rho[0] + int(W[1, 1]) * rho[1],
        )
        num += sign * z1 ** Wlr[0] * z2 ** Wlr[1]
        den += sign * z1 ** Wr[0] * z2 ** Wr[1]

    # Where denominator is near zero (Weyl walls), the character has a
    # removable singularity equal to the dimension of V_lam. Walls are
    # codim-1 sets that contribute measure zero to the integration; for
    # safety we replace with the dim.
    with np.errstate(divide="ignore", invalid="ignore"):
        chi = num / den
    bad = ~np.isfinite(chi) | (np.abs(den) < 1e-12)
    if np.any(bad):
        dim_val = la.dim_weyl(lam)
        chi = np.where(bad, complex(float(dim_val), 0.0), chi)
    return chi


def weyl_denominator_sq(
    la: SimpleLieAlgebra, T1: np.ndarray, T2: np.ndarray
) -> np.ndarray:
    """|Delta(theta)|^2 = prod_{alpha > 0} 4 sin^2(<alpha, theta>/2)
    where <alpha, theta> = a_1 theta_1 + a_2 theta_2 for alpha = (a_1, a_2) in
    omega basis."""
    out = np.ones_like(T1)
    for alpha in la.positive_roots_omega:
        a1, a2 = int(alpha[0]), int(alpha[1])
        phase = (a1 * T1 + a2 * T2) / 2.0
        out = out * 4.0 * np.sin(phase) ** 2
    return out


def geodesic_dist_sq_grid(
    la: SimpleLieAlgebra,
    T1: np.ndarray,
    T2: np.ndarray,
    n_search: int = 1,
) -> np.ndarray:
    """Bi-invariant Killing-form geodesic distance squared from identity to
    theta in dual-omega coords.

    The metric on the maximal torus in dual-omega coordinates is
        g^{-1} = G_omega^{-1}
    where G_omega is the Gram matrix on the omega basis.

    The lattice of identity elements (kernel of exp) for a simply-connected
    compact Lie group, in dual-omega coords, is 2 pi Z^2 (since dual-omega
    coords are dual to fundamental weights, and the coroot lattice in the
    fundamental-weight basis is precisely the integer lattice -- coroots
    pair with fundamental weights as alpha_i^vee . omega_j = delta_{ij}).

    So we minimize over (n_1, n_2) in [-n_search, n_search]^2.
    """
    G_inv = la.G_omega.inv()
    G11 = float(G_inv[0, 0])
    G12 = float(G_inv[0, 1])  # = G21
    G22 = float(G_inv[1, 1])

    two_pi = 2.0 * math.pi
    best = None
    for n1 in range(-n_search, n_search + 1):
        for n2 in range(-n_search, n_search + 1):
            t1p = T1 - n1 * two_pi
            t2p = T2 - n2 * two_pi
            dsq = G11 * t1p * t1p + 2.0 * G12 * t1p * t2p + G22 * t2p * t2p
            if best is None:
                best = dsq
            else:
                best = np.minimum(best, dsq)
    return best


def verify_haar_normalization(
    la: SimpleLieAlgebra, n_quad: int = 60
) -> float:
    """Verify (1/(|W| (2 pi)^2)) int_T^2 |Delta|^2 dt = 1 for the given group.

    Uses Gauss-Legendre on [0, 2 pi]^2.
    """
    weyl_orbit = precompute_weyl_orbit(la)
    W_order = len(weyl_orbit)

    two_pi = 2.0 * math.pi
    nodes, weights = _gl_nodes(0.0, two_pi, n_quad)
    T1, T2 = np.meshgrid(nodes, nodes, indexing="ij")
    W1, W2 = np.meshgrid(weights, weights, indexing="ij")

    measure = weyl_denominator_sq(la, T1, T2)
    integral = np.sum(measure * W1 * W2)
    integral = integral / (W_order * two_pi * two_pi)
    return float(integral), W_order


def gamma_lambda(
    la: SimpleLieAlgebra,
    irreps: List[Tuple[int, int]],
    weyl_orbit: list,
    n_quad: int = 40,
) -> float:
    """Compute gamma_Lambda = int_G K_Lambda(g) d(e, g) dg via 2D Gauss-Legendre
    on the maximal torus.

    K_Lambda(theta) = (1/Z) |sum_pi sqrt(dim pi) chi_pi(theta)|^2,
    Z = sum_pi dim pi.

    Note: d(e, g) is the geodesic distance, NOT its square. Mirrors Paper 38
    convention: gamma is the *first moment* of the geodesic distance under K.
    (The "mass-concentration moment of d^2" notation in the prompt is a
    slight informalism; the actual quantity is the first moment of d, which
    is the standard Stein-Weiss moment.)
    """
    W_order = len(weyl_orbit)
    two_pi = 2.0 * math.pi
    nodes, weights = _gl_nodes(0.0, two_pi, n_quad)
    T1, T2 = np.meshgrid(nodes, nodes, indexing="ij")
    W1, W2 = np.meshgrid(weights, weights, indexing="ij")

    # Geodesic distance grid (sqrt of d^2_min).
    dsq = geodesic_dist_sq_grid(la, T1, T2)
    d_grid = np.sqrt(np.maximum(dsq, 0.0))

    # Weyl denominator |Delta|^2.
    measure = weyl_denominator_sq(la, T1, T2)

    # Build kernel K(theta) = (1/Z) |sum_pi sqrt(dim_pi) chi_pi(theta)|^2.
    Z = sum(int(la.dim_weyl(lam)) for lam in irreps)
    sum_grid = np.zeros_like(T1, dtype=np.complex128)
    for lam in irreps:
        dim_lam = int(la.dim_weyl(lam))
        chi = character_grid(la, lam, weyl_orbit, T1, T2)
        sum_grid += math.sqrt(dim_lam) * chi
    K_grid = (np.abs(sum_grid) ** 2) / Z

    integrand = K_grid * d_grid * measure * W1 * W2
    integral = np.sum(integrand) / (W_order * two_pi * two_pi)
    return float(integral.real)


# =============================================================================
# Panel runner and fits
# =============================================================================


def run_panel(
    la: SimpleLieAlgebra,
    name: str,
    lambda_sq_panel: List[int],
) -> Tuple[List[dict], dict]:
    """Run gamma_Lambda computation across panel of Lambda^2 values."""
    print(f"\n{'=' * 78}")
    print(f"{name}: gamma_Lambda panel")
    print(f"{'=' * 78}")

    # Precompute Weyl orbit once.
    weyl_orbit = precompute_weyl_orbit(la)
    W_order = len(weyl_orbit)
    print(f"  |W| = {W_order}")
    print(f"  rank = {la.rank}")
    print(f"  positive roots = {la.positive_roots_omega}")

    # Verify Haar normalization at moderate n_quad.
    haar, _ = verify_haar_normalization(la, n_quad=60)
    print(f"  Haar integration check: int_G 1 dg = {haar:.10f}  (expect 1.0)")
    if abs(haar - 1.0) > 1e-3:
        print(f"  WARNING: Haar normalization off by {haar - 1.0:.3e}")
        # Re-run at higher n_quad to confirm
        haar2, _ = verify_haar_normalization(la, n_quad=100)
        print(f"  Haar @ n_quad=100: {haar2:.10f}")

    rows = []
    print(
        f"\n  {'Lam^2':>6} {'Lambda':>8} {'max_d':>6} {'n_irr':>6} {'max_a+b':>8} "
        f"{'n_quad':>7} {'gamma':>10} {'c_est':>8} {'t(s)':>6}"
    )
    for lsq in lambda_sq_panel:
        irreps, md, n_irr = find_saturated_panel(la, lsq)
        max_ab = max(a + b for (a, b) in irreps) if irreps else 0
        # n_quad: at least 2 * max(a+b) + 30 for good quadrature
        nq = max(40, 2 * max_ab + 30)
        nq = min(nq, 250)
        t0 = time.time()
        g = gamma_lambda(la, irreps, weyl_orbit, n_quad=nq)
        dt = time.time() - t0
        L = math.sqrt(lsq)
        c_est = L * g / math.log(L) if L > 1.5 else float("nan")
        rows.append({
            "lambda_sq": lsq,
            "lambda": L,
            "gamma": g,
            "n_irreps": n_irr,
            "max_dynkin_for_saturation": md,
            "max_a_plus_b": max_ab,
            "n_quad": nq,
            "time_seconds": dt,
            "L*gamma/log(L)": c_est,
        })
        print(
            f"  {lsq:>6} {L:>8.3f} {md:>6} {n_irr:>6} {max_ab:>8} "
            f"{nq:>7} {g:>10.6f} {c_est:>8.4f} {dt:>6.1f}"
        )

    meta = {
        "weyl_group_order": W_order,
        "haar_check": haar,
    }
    return rows, meta


def fit_models(
    lams: List[float],
    gammas: List[float],
    exclude_below_L: float = 0.0,
    label: str = "",
) -> dict:
    """Run 2-parameter Stein-Weiss, 3-param single-log, double-log fits."""
    pairs = [(L, g) for L, g in zip(lams, gammas) if L >= exclude_below_L]
    Lp = [p[0] for p in pairs]
    Gp = [p[1] for p in pairs]
    Y = np.array(Gp)
    n_data = len(Y)
    fits = {"label": label, "n_data": n_data, "exclude_below_L": exclude_below_L}

    if n_data < 2:
        return fits

    # 2-parameter Stein-Weiss: gamma = (a log L + b) / L
    X = np.array([[math.log(L) / L, 1.0 / L] for L in Lp])
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    rss = float(np.sum((Y - X @ coef) ** 2))
    fits["two_param"] = {
        "form": "(a log L + b)/L",
        "a": float(coef[0]),
        "b": float(coef[1]),
        "rss": rss,
        "aic": n_data * math.log(rss / n_data) + 4 if rss > 0 else float("-inf"),
        "bic": n_data * math.log(rss / n_data) + 2 * math.log(n_data) if rss > 0 else float("-inf"),
    }

    if n_data < 3:
        return fits

    # 3-parameter single-log: gamma = a1 log L/L + b/L + c/L^2
    X = np.array([[math.log(L) / L, 1.0 / L, 1.0 / (L * L)] for L in Lp])
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    rss = float(np.sum((Y - X @ coef) ** 2))
    fits["single_log"] = {
        "a1": float(coef[0]),
        "b": float(coef[1]),
        "c": float(coef[2]),
        "rss": rss,
        "aic": n_data * math.log(rss / n_data) + 6 if rss > 0 else float("-inf"),
    }

    # 3-parameter double-log: gamma = a2 (log L)^2/L + b2 log L/L + c2/L
    X = np.array([[(math.log(L) ** 2) / L, math.log(L) / L, 1.0 / L] for L in Lp])
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    rss = float(np.sum((Y - X @ coef) ** 2))
    fits["double_log"] = {
        "a2": float(coef[0]),
        "b2": float(coef[1]),
        "c2": float(coef[2]),
        "rss": rss,
        "aic": n_data * math.log(rss / n_data) + 6 if rss > 0 else float("-inf"),
    }

    # Stein-Weiss 3-param: gamma = a log L/L + a b/L + a c/(L log L)
    X = np.array([[math.log(L) / L, 1.0 / L, 1.0 / (L * math.log(L))] for L in Lp])
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    rss = float(np.sum((Y - X @ coef) ** 2))
    fits["stein_weiss_3"] = {
        "a": float(coef[0]),
        "ab": float(coef[1]),
        "ac": float(coef[2]),
        "rss": rss,
        "aic": n_data * math.log(rss / n_data) + 6 if rss > 0 else float("-inf"),
    }

    return fits


def verdict_for_group(
    name: str,
    a_extracted_gen: float,
    weyl_order: int,
    rank: int,
    casimir_rescale: float,
) -> dict:
    """Compare extracted a against A (4/pi) and B (2|W|/pi^r).

    Convention. The conjectures are stated in the *canonical* (dual-Coxeter)
    normalization: C(adjoint) = h^vee. My script uses a "generic" normalization
    where short roots have squared length 2. To compare:

        a_canonical = a_generic / sqrt(casimir_rescale)

    where casimir_rescale = (C_gen of adjoint) / h^vee. For SU(3): 6/3 = 2,
    so sqrt(2). For Sp(2): 12/3 = 4, so sqrt(4) = 2. For G_2: 24/4 = 6, so
    sqrt(6).

    Predictions in canonical convention:
        A (rank-invariant):  c = 4/pi
        B (Weyl-formula):    c = 2 |W| / pi^rank

    Verdict A/B/AMBIGUOUS/HYBRID based on which is closer.
    """
    pi = math.pi
    cA = 4.0 / pi
    cB = 2.0 * weyl_order / (pi ** rank)

    a_canonical = a_extracted_gen / math.sqrt(casimir_rescale)

    err_A = abs(a_canonical - cA) / cA * 100
    err_B = abs(a_canonical - cB) / cB * 100

    gap_AB_pct = abs(cA - cB) / max(cA, cB) * 100

    if err_A < 5.0 and err_B > 15.0:
        verdict = (
            f"A ({name} a_can={a_canonical:.4f} matches 4/pi={cA:.4f} at {err_A:.2f}%)"
        )
    elif err_B < 5.0 and err_A > 15.0:
        verdict = (
            f"B ({name} a_can={a_canonical:.4f} matches 2|W|/pi^r={cB:.4f} at {err_B:.2f}%)"
        )
    elif err_A < 5.0 and err_B < 5.0:
        verdict = (
            f"AMBIGUOUS ({name} both within 5%: A err={err_A:.2f}%, B err={err_B:.2f}%)"
        )
    elif err_A < err_B:
        verdict = (
            f"A-leaning ({name} a_can={a_canonical:.4f}, A err={err_A:.2f}%, B err={err_B:.2f}%)"
        )
    elif err_B < err_A:
        verdict = (
            f"B-leaning ({name} a_can={a_canonical:.4f}, B err={err_B:.2f}%, A err={err_A:.2f}%)"
        )
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
    print("Sp(2) and G_2 rate constant analysis (Conjecture A vs B discriminator)")
    print("=" * 78)

    results = {
        "metadata": {
            "purpose": "Discriminate Conjecture A (c=4/pi rank-invariant) vs B (c=2|W|/pi^r) "
            "on rank-2 non-SU(N) groups",
            "predictions": {
                "A_4_over_pi": 4.0 / math.pi,
                "Sp2_B_16_over_pi2": 16.0 / math.pi ** 2,
                "G2_B_24_over_pi2": 24.0 / math.pi ** 2,
                "Sp2_gap_AB_pct": abs(4.0 / math.pi - 16.0 / math.pi ** 2)
                / max(4.0 / math.pi, 16.0 / math.pi ** 2) * 100,
                "G2_gap_AB_pct": abs(4.0 / math.pi - 24.0 / math.pi ** 2)
                / max(4.0 / math.pi, 24.0 / math.pi ** 2) * 100,
            },
        }
    }

    # =========================================================================
    # SU(3) sanity reproduction (small panel) -- confirms our generic
    # machinery reproduces the SU(3) c value to within ~5%
    # =========================================================================
    print("\n" + "=" * 78)
    print("SU(3) sanity reproduction (small panel) using generic dual-omega code")
    print("=" * 78)
    a2 = build_A(2)
    # SU(3) generic panel matched to canonical SU(3)'s extended range
    # (factor 2 between canonical and generic Casimirs).
    # Canonical [4..1000] -> generic [8..2000]
    su3_panel = [8, 12, 16, 20, 28, 36, 48, 60, 80, 100, 140, 200, 300, 400, 600, 800, 1000, 1400, 2000]
    rows_su3, meta_su3 = run_panel(a2, "SU(3)/A_2", su3_panel)
    lams_su3 = [r["lambda"] for r in rows_su3]
    gs_su3 = [r["gamma"] for r in rows_su3]
    # Canonical L_can >= 7 corresponds to L_gen >= 7*sqrt(2) ~ 9.9
    fits_su3 = fit_models(lams_su3, gs_su3, exclude_below_L=9.9, label="SU(3) sanity L>=7_can")
    a_su3_gen = fits_su3["two_param"]["a"]
    # SU(3) dual-Coxeter rescale: C_gen(adjoint=(1,1)) = 6, h^vee = 3. Ratio: 6/3 = 2.
    su3_rescale = 2.0
    a_su3_can = a_su3_gen / math.sqrt(su3_rescale)
    print(f"\n  SU(3) sanity 2-param a_gen (L>=5): {a_su3_gen:.6f}")
    print(f"  SU(3) canonical-rescaled (factor sqrt(2)): a_can = {a_su3_can:.6f}")
    print(f"  Expected canonical c ~ 1.24 (from su3_rate_constant.py)")
    print(f"  (Note: small panel L<=22.4 is not yet asymptotic; canonical at full panel "
          f"L<=31.6 gave a=1.24, see memo §3)")
    results["su3_sanity"] = {
        "rows": rows_su3,
        "fits": fits_su3,
        "a_gen_two_param": a_su3_gen,
        "casimir_rescale": su3_rescale,
        "a_canonical": a_su3_can,
        "weyl_meta": meta_su3,
    }

    # =========================================================================
    # Sp(2) (C_2) panel
    # =========================================================================
    c2 = build_C2()
    # Sp(2) generic Lambda^2 panel. Sp(2) canonical Lambda^2 = generic Lambda^2 / 4.
    # Canonical SU(3) panel was [4..1000] -> for Sp(2) we use generic [16..4000]
    # so canonical equivalent panel is [4..1000] (matched to SU(3) canonical).
    sp2_panel = [4, 8, 16, 24, 32, 56, 72, 96, 120, 160, 200, 280, 400, 600, 800, 1200, 1600, 2400, 3200, 4000]
    rows_sp2, meta_sp2 = run_panel(c2, "Sp(2)/C_2", sp2_panel)
    results["sp2"] = {"rows": rows_sp2, "meta": meta_sp2}

    lams_sp2 = [r["lambda"] for r in rows_sp2]
    gs_sp2 = [r["gamma"] for r in rows_sp2]
    fits_sp2_all = {}
    # Note: Sp(2) canonical Lambda = generic Lambda / 2.
    # Canonical SU(3) used exclusion L_can >= 7 (asymptotic). For Sp(2),
    # equivalent is L_gen >= 14.
    for lo in [0.0, 8.0, 14.0, 20.0, 30.0]:
        fits_sp2_all[f"min_L_{lo}"] = fit_models(
            lams_sp2, gs_sp2, exclude_below_L=lo, label=f"Sp(2) L>={lo}"
        )
    print("\n=== Sp(2) fits (L_gen >= X corresponds to L_can >= X/2) ===")
    for key, f in fits_sp2_all.items():
        if "two_param" not in f:
            continue
        a = f["two_param"]["a"]
        a_can = a / 2.0
        print(f"  {key}: n_data={f['n_data']}, a_gen = {a:.6f}, a_can = {a_can:.6f}")
    results["sp2"]["fits"] = fits_sp2_all

    # Primary Sp(2) estimate: 2-param Stein-Weiss on L_gen >= 14 (= L_can >= 7,
    # asymptotic regime mirroring su3_rate_constant.py convention)
    a_sp2_gen = fits_sp2_all["min_L_14.0"]["two_param"]["a"]
    # Sp(2) dual-Coxeter rescale: my generic Casimir on adjoint (2,0) is 12,
    # canonical h^vee = 3. Ratio: 12/3 = 4. So a_can = a_gen / 2.
    sp2_rescale = 4.0  # = C_gen(adjoint) / h^vee
    print(f"\n  Sp(2) primary estimate (generic): a_gen = {a_sp2_gen:.6f}")
    print(f"  Casimir rescale factor: {sp2_rescale} (= C_gen(adj=(2,0))/h^vee = 12/3)")
    print(f"  Sp(2) canonical: a_can = {a_sp2_gen/math.sqrt(sp2_rescale):.6f}")
    print(f"  Conjecture A: 4/pi  = {4/math.pi:.6f}")
    print(f"  Conjecture B: 16/pi^2 = {16/math.pi**2:.6f}")
    print(f"  Gap A vs B (Sp(2)): {abs(4/math.pi - 16/math.pi**2)/(16/math.pi**2)*100:.1f}%")
    verd_sp2 = verdict_for_group("Sp(2)", a_sp2_gen, weyl_order=8, rank=2,
                                 casimir_rescale=sp2_rescale)
    print(f"  VERDICT: {verd_sp2['verdict']}")
    results["sp2"]["verdict"] = verd_sp2

    # =========================================================================
    # G_2 panel
    # =========================================================================
    g2 = build_G2()
    # G_2 generic Lambda^2 panel. G_2 canonical Lambda^2 = generic Lambda^2 / 6.
    # To match canonical SU(3) range [4..1000], use generic [24..6000].
    g2_panel = [12, 24, 36, 48, 84, 108, 144, 180, 240, 300, 420, 600, 900, 1200, 1800, 2400, 3000, 4200, 6000]
    rows_g2, meta_g2 = run_panel(g2, "G_2", g2_panel)
    results["g2"] = {"rows": rows_g2, "meta": meta_g2}

    lams_g2 = [r["lambda"] for r in rows_g2]
    gs_g2 = [r["gamma"] for r in rows_g2]
    fits_g2_all = {}
    # G_2 canonical Lambda = generic Lambda / sqrt(6).
    # Canonical L_can >= 7 corresponds to L_gen >= 7*sqrt(6) ~ 17.15.
    for lo in [0.0, 10.0, 17.15, 24.5, 35.0]:
        fits_g2_all[f"min_L_{lo}"] = fit_models(
            lams_g2, gs_g2, exclude_below_L=lo, label=f"G_2 L>={lo}"
        )
    print("\n=== G_2 fits (L_gen >= X corresponds to L_can >= X/sqrt(6)) ===")
    for key, f in fits_g2_all.items():
        if "two_param" not in f:
            continue
        a = f["two_param"]["a"]
        a_can = a / math.sqrt(6)
        print(f"  {key}: n_data={f['n_data']}, a_gen = {a:.6f}, a_can = {a_can:.6f}")
    results["g2"]["fits"] = fits_g2_all

    # Primary G_2 estimate: 2-param Stein-Weiss on L_gen >= 17.15 (= L_can >= 7)
    a_g2_gen = fits_g2_all["min_L_17.15"]["two_param"]["a"]
    # G_2 dual-Coxeter rescale: my generic Casimir on adjoint (0,1) is 24,
    # canonical h^vee = 4. Ratio: 24/4 = 6.
    g2_rescale = 6.0  # = C_gen(adjoint=(0,1))/h^vee = 24/4
    print(f"\n  G_2 primary estimate (generic): a_gen = {a_g2_gen:.6f}")
    print(f"  Casimir rescale factor: {g2_rescale} (= C_gen(adj=(0,1))/h^vee = 24/4)")
    print(f"  G_2 canonical: a_can = {a_g2_gen/math.sqrt(g2_rescale):.6f}")
    print(f"  Conjecture A: 4/pi    = {4/math.pi:.6f}")
    print(f"  Conjecture B: 24/pi^2 = {24/math.pi**2:.6f}")
    print(f"  Gap A vs B (G_2): {abs(4/math.pi - 24/math.pi**2)/(24/math.pi**2)*100:.1f}%")
    verd_g2 = verdict_for_group("G_2", a_g2_gen, weyl_order=12, rank=2,
                                casimir_rescale=g2_rescale)
    print(f"  VERDICT: {verd_g2['verdict']}")
    results["g2"]["verdict"] = verd_g2

    # =========================================================================
    # Combined verdict
    # =========================================================================
    print("\n" + "=" * 78)
    print("COMBINED VERDICT")
    print("=" * 78)
    v_sp2_str = verd_sp2["verdict"]
    v_g2_str = verd_g2["verdict"]
    print(f"  Sp(2): {v_sp2_str}")
    print(f"  G_2:   {v_g2_str}")

    # Determine combined reading
    sp2_lean = "A" if "A" in v_sp2_str.split()[0] else ("B" if "B" in v_sp2_str.split()[0] else "AMBIG")
    g2_lean = "A" if "A" in v_g2_str.split()[0] else ("B" if "B" in v_g2_str.split()[0] else "AMBIG")

    if sp2_lean == "A" and g2_lean == "A":
        combined = "Conjecture A (rank-invariant c = 4/pi)"
    elif sp2_lean == "B" and g2_lean == "B":
        combined = "Conjecture B (Weyl-formula c = 2|W|/pi^r)"
    elif sp2_lean == g2_lean:
        combined = f"Combined {sp2_lean}"
    else:
        combined = f"HYBRID (Sp(2): {sp2_lean}, G_2: {g2_lean})"
    print(f"  Combined: {combined}")

    results["combined_verdict"] = {
        "sp2_lean": sp2_lean,
        "g2_lean": g2_lean,
        "combined": combined,
    }

    out_path = out_dir / "sp2_g2_rate_constant.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")

    return results


if __name__ == "__main__":
    main()
