"""
Hopf-U(1) block decomposition test (Track RH-E).
===================================================

Tests Paper 29 Hypothesis 1 (§5.3) and its Dirac-S^3 analog (Rule B n_max=3):
that the integer-coefficient factor polynomials of the Ihara zeta correspond
to Hashimoto operator restrictions to U(1)-equivariant sectors.

Scalar S^5 Bargmann-Segal N_max=3:
    ζ_G(s)^{-1} = (s-1)^23 (s+1)^23 · P_12(s^2) · P_22(s^2)
    P_12(s^2) = 432 s^12 + ... + 1       (degree 12 in s)
    P_22(s^2) = 829440 s^22 + ... - 1    (degree 22 in s)

Dirac-S^3 Rule B n_max=3:
    ζ_G(s)^{-1} = (s-1)^79 (s+1)^79 (9s^2+1)^4 · P_22(s^2) · P_24(s^2)

Methodology
-----------
By Ihara-Bass, ζ_G(s)^{-1} = (1 - s^2)^{r-c} · det(I - sA + s^2 Q), where
A is the node-level adjacency and Q = diag(d_v - 1). The (1-s^2)^{r-c}
prefactor absorbs all β_1 = r - c "gauge" trivial zeros at s = ±1; the
NON-TRIVIAL factor polynomials live in the V×V node-level determinant.
This reduces the problem from a 2E-dim edge-operator eigendecomposition
to a V-dim NODE eigendecomposition.

U(1) action. The natural U(1) acting on the graph is the phase rotation
of node wavefunctions,  ψ_v → e^{i χ m(v)} ψ_v, where m = m_l (scalar S^5)
or m = m_j (Dirac). The full U(1) has *continuous* character, but on a
real integer adjacency A it reduces to the Z_2 reflection m → -m. (The
continuous U(1) only commutes with A if A is m-translation-invariant,
which it is not for a finite graph with bounded m-range.) The m-reflection
P is an exact symmetry of every m-respecting adjacency:
  - P A = A P      (adjacency commutes with reflection)
  - P Q = Q P      (degree diagonal trivially commutes)
  - P^2 = I        (involution)

Therefore the V-dim space splits orthogonally as V = V_{sym} ⊕ V_{antisym}
(m → -m even vs odd eigenspaces of P), and the matrix (I - sA + s^2 Q)
is block-diagonal in this decomposition. The two block determinants
multiply to the full det, and — by the Paper 29 conjecture — should
correspond to the two factor polynomials P_12, P_22 (scalar) or
P_22, P_24 (Dirac).

Verdict criterion. If HYP 1 is correct:
  (a) det(M_sym) × det(M_antisym) = full det(I - sA + s^2 Q).
  (b) Each block's factored form contains exactly one of the target factor
      polynomials P_d, plus (possibly) small "bookkeeping" factors
      like (s ± 1), (9s^2 + 1)^{#}, which are allowed because the
      prefactor (1-s^2)^{r-c} and the (9s^2+1)^4 factor in Dirac are
      not themselves part of the P_d polynomials.

Author: GeoVac Track RH-E, April 2026.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import sympy as sp

from geovac.nuclear.bargmann_graph import build_bargmann_graph
from geovac.ihara_zeta_dirac import build_dirac_s3_graph


# ---------------------------------------------------------------------------
# Target polynomials
# ---------------------------------------------------------------------------

def target_P12_P22_scalar_S5() -> Tuple[sp.Poly, sp.Poly]:
    """Scalar S^5 Bargmann-Segal N_max=3 factor polynomials."""
    s = sp.symbols("s")
    P12 = sp.Poly(
        432 * s**12 + 666 * s**10 + 374 * s**8
        + 135 * s**6 + 47 * s**4 + 11 * s**2 + 1,
        s,
    )
    P22 = sp.Poly(
        829440 * s**22 + 2453184 * s**20 + 3308104 * s**18
        + 2696682 * s**16 + 1470640 * s**14 + 557227 * s**12
        + 146654 * s**10 + 25709 * s**8 + 2630 * s**6
        + 79 * s**4 - 12 * s**2 - 1,
        s,
    )
    return P12, P22


def target_P22_P24_dirac() -> Tuple[sp.Poly, sp.Poly]:
    """Dirac-S^3 Rule B n_max=3 factor polynomials."""
    s = sp.symbols("s")
    P22 = sp.Poly(
        538876800 * s**22 + 1088750160 * s**20 + 925413984 * s**18
        + 456734400 * s**16 + 148182464 * s**14 + 33353711 * s**12
        + 5288330 * s**10 + 580341 * s**8 + 41332 * s**6
        + 1589 * s**4 + 10 * s**2 - 1,
        s,
    )
    P24 = sp.Poly(
        538876800 * s**24 + 1108028160 * s**22 + 934113744 * s**20
        + 450036324 * s**18 + 145121264 * s**16 + 35278007 * s**14
        + 7106943 * s**12 + 1222735 * s**10 + 170871 * s**8
        + 17793 * s**6 + 1257 * s**4 + 53 * s**2 + 1,
        s,
    )
    return P22, P24


# ---------------------------------------------------------------------------
# Adjacency + node labels
# ---------------------------------------------------------------------------

def bargmann_scalar_data() -> Tuple[np.ndarray, List[Tuple[int, int, int]], Dict]:
    """Return (A, nodes, label_info) for scalar S^5 Bargmann-Segal N_max=3."""
    N_max = 3
    g = build_bargmann_graph(N_max)
    V = g.n_nodes
    A = np.zeros((V, V), dtype=int)
    for (i, j), w in g.adjacency.items():
        A[i, j] = 1
        A[j, i] = 1
    info = {"N_max": N_max, "V": V, "E": int(A.sum()) // 2}
    return A, list(g.nodes), info


def dirac_ruleB_data() -> Tuple[np.ndarray, List, Dict]:
    """Return (A, labels, label_info) for Dirac-S^3 Rule B n_max=3."""
    n_max = 3
    A, labels, deg, desc = build_dirac_s3_graph(n_max, "B")
    V = A.shape[0]
    info = {"n_max": n_max, "V": V, "E": int(A.sum()) // 2, "rule": "B"}
    return A, list(labels), info


# ---------------------------------------------------------------------------
# m-reflection permutation
# ---------------------------------------------------------------------------

def bargmann_m_reflection(nodes: List[Tuple[int, int, int]]) -> np.ndarray:
    """Build the permutation matrix P that maps (N, l, m_l) → (N, l, -m_l)."""
    V = len(nodes)
    idx_map = {nd: i for i, nd in enumerate(nodes)}
    P = np.zeros((V, V), dtype=int)
    for i, (N, l, m) in enumerate(nodes):
        j = idx_map[(N, l, -m)]
        P[i, j] = 1
    return P


def dirac_mj_reflection(labels: List) -> np.ndarray:
    """Build the permutation matrix P mapping (n, κ, 2m_j) → (n, κ, -2m_j)."""
    V = len(labels)
    P = np.zeros((V, V), dtype=int)
    for i, lab in enumerate(labels):
        for j, lab2 in enumerate(labels):
            if (lab2.n_fock == lab.n_fock and lab2.kappa == lab.kappa
                    and lab2.two_m_j == -lab.two_m_j):
                P[i, j] = 1
                break
    return P


# ---------------------------------------------------------------------------
# Sym / antisym subspace construction
# ---------------------------------------------------------------------------

def build_reflection_projectors_bargmann(nodes: List[Tuple[int, int, int]]):
    """Return (U_sym, U_antisym) sympy matrices for the m → -m projection.

    U_sym is V × dim_sym = V × 13 (for N_max=3: 6 fixed m=0 nodes + 7 pair-symm).
    U_antisym is V × 7 (7 pair-antisym).
    """
    V = len(nodes)
    sqrt2 = sp.sqrt(2)
    sym_vectors: List[sp.Matrix] = []
    antisym_vectors: List[sp.Matrix] = []
    idx_map = {nd: i for i, nd in enumerate(nodes)}
    seen = set()
    for i, (N, l, m) in enumerate(nodes):
        if m == 0:
            v = sp.zeros(V, 1)
            v[i] = 1
            sym_vectors.append(v)
        elif m > 0 and (N, l, m) not in seen:
            seen.add((N, l, m))
            j = idx_map[(N, l, -m)]
            vs = sp.zeros(V, 1)
            vs[i] = sp.Rational(1, 2) * sqrt2
            vs[j] = sp.Rational(1, 2) * sqrt2
            va = sp.zeros(V, 1)
            va[i] = sp.Rational(1, 2) * sqrt2
            va[j] = -sp.Rational(1, 2) * sqrt2
            sym_vectors.append(vs)
            antisym_vectors.append(va)
    U_sym = sp.Matrix.hstack(*sym_vectors)
    U_antisym = sp.Matrix.hstack(*antisym_vectors)
    return U_sym, U_antisym


def build_reflection_projectors_dirac(labels: List):
    """Return (U_sym, U_antisym) sympy matrices for the m_j → -m_j projection.

    No fixed points (m_j is half-integer), so U_sym and U_antisym both have
    dim V/2 = 14 (for n_max=3, V=28).
    """
    V = len(labels)
    sqrt2 = sp.sqrt(2)
    sym_vectors: List[sp.Matrix] = []
    antisym_vectors: List[sp.Matrix] = []
    seen = set()
    for i, lab in enumerate(labels):
        if lab.two_m_j > 0 and (lab.n_fock, lab.kappa, lab.two_m_j) not in seen:
            for j, lab2 in enumerate(labels):
                if (lab2.n_fock == lab.n_fock and lab2.kappa == lab.kappa
                        and lab2.two_m_j == -lab.two_m_j):
                    seen.add((lab.n_fock, lab.kappa, lab.two_m_j))
                    vs = sp.zeros(V, 1)
                    vs[i] = sp.Rational(1, 2) * sqrt2
                    vs[j] = sp.Rational(1, 2) * sqrt2
                    va = sp.zeros(V, 1)
                    va[i] = sp.Rational(1, 2) * sqrt2
                    va[j] = -sp.Rational(1, 2) * sqrt2
                    sym_vectors.append(vs)
                    antisym_vectors.append(va)
                    break
    U_sym = sp.Matrix.hstack(*sym_vectors)
    U_antisym = sp.Matrix.hstack(*antisym_vectors)
    return U_sym, U_antisym


# ---------------------------------------------------------------------------
# Block determinant
# ---------------------------------------------------------------------------

def block_det_I_minus_sA_plus_sQ(A_np: np.ndarray, U: sp.Matrix,
                                  Q_np: np.ndarray):
    """Compute det(I - s A_block + s^2 Q_block) for the block defined by U.

    Parameters
    ----------
    A_np : V × V integer adjacency
    U    : V × d sympy projection (orthogonal basis of block)
    Q_np : V × V integer diagonal matrix = diag(d_v - 1)

    Returns
    -------
    sympy polynomial in s (the block determinant, expanded).
    """
    V, d = U.shape
    A_sp = sp.Matrix(A_np.tolist())
    Q_sp = sp.Matrix(Q_np.tolist())
    A_block = sp.simplify(U.T * A_sp * U)
    Q_block = sp.simplify(U.T * Q_sp * U)
    s = sp.symbols("s")
    M = sp.eye(d) - s * A_block + s * s * Q_block
    return sp.expand(M.det())


# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def run_scalar_s5_experiment() -> Dict:
    """Test Hypothesis 1 on scalar S^5 Bargmann-Segal N_max=3."""
    A, nodes, info = bargmann_scalar_data()
    V = info["V"]
    E = info["E"]
    deg = A.sum(axis=1)
    Q = np.diag(deg - 1)

    # Build m-reflection and verify it commutes with A
    P = bargmann_m_reflection(nodes)
    comm = P @ A - A @ P
    assert np.abs(comm).sum() == 0, "A must commute with m-reflection"

    # Fixed points and pair counts
    fixed = [i for i, (N, l, m) in enumerate(nodes) if m == 0]
    non_fixed = V - len(fixed)
    dim_sym_expected = len(fixed) + non_fixed // 2
    dim_antisym_expected = non_fixed // 2

    # Build projectors
    U_sym, U_antisym = build_reflection_projectors_bargmann(nodes)
    assert U_sym.shape == (V, dim_sym_expected)
    assert U_antisym.shape == (V, dim_antisym_expected)

    # Verify orthogonality
    I_sym = sp.simplify(U_sym.T * U_sym)
    I_antisym = sp.simplify(U_antisym.T * U_antisym)
    cross = sp.simplify(U_sym.T * U_antisym)
    assert I_sym == sp.eye(dim_sym_expected), "U_sym not orthonormal"
    assert I_antisym == sp.eye(dim_antisym_expected), "U_antisym not orthonormal"
    assert cross == sp.zeros(dim_sym_expected, dim_antisym_expected), \
        "U_sym, U_antisym not orthogonal"

    # Compute block dets
    print(f"  Computing det(M_sym) (dim {dim_sym_expected}) ...")
    det_sym = block_det_I_minus_sA_plus_sQ(A, U_sym, Q)
    print(f"  Computing det(M_antisym) (dim {dim_antisym_expected}) ...")
    det_antisym = block_det_I_minus_sA_plus_sQ(A, U_antisym, Q)

    s = sp.symbols("s")
    P12, P22 = target_P12_P22_scalar_S5()

    det_sym_factored = sp.factor(det_sym)
    det_antisym_factored = sp.factor(det_antisym)
    det_sym_poly = sp.Poly(det_sym, s)
    det_antisym_poly = sp.Poly(det_antisym, s)

    # Check whether each factor polynomial divides the block det
    sym_contains_P22 = sp.rem(det_sym, P22.as_expr(), s) == 0
    sym_contains_P12 = sp.rem(det_sym, P12.as_expr(), s) == 0
    antisym_contains_P22 = sp.rem(det_antisym, P22.as_expr(), s) == 0
    antisym_contains_P12 = sp.rem(det_antisym, P12.as_expr(), s) == 0

    # Check that product equals full det(I - sA + s²Q)
    full_M = sp.eye(V) - s * sp.Matrix(A.tolist()) + s * s * sp.Matrix(Q.tolist())
    full_det = sp.expand(full_M.det())
    product = sp.expand(det_sym * det_antisym)
    product_matches = sp.simplify(product - full_det) == 0

    return {
        "graph": "scalar S^5 Bargmann-Segal N_max=3",
        "V": V,
        "E": E,
        "beta_1": E - V + 1,
        "dim_sym": dim_sym_expected,
        "dim_antisym": dim_antisym_expected,
        "fixed_points_count": len(fixed),
        "det_sym_degree": det_sym_poly.degree(),
        "det_antisym_degree": det_antisym_poly.degree(),
        "det_sym_factored": str(det_sym_factored),
        "det_antisym_factored": str(det_antisym_factored),
        "sym_block_contains_P12": bool(sym_contains_P12),
        "sym_block_contains_P22": bool(sym_contains_P22),
        "antisym_block_contains_P12": bool(antisym_contains_P12),
        "antisym_block_contains_P22": bool(antisym_contains_P22),
        "product_equals_full_det": bool(product_matches),
        "target_P12_degree": P12.degree(),
        "target_P22_degree": P22.degree(),
        "target_P12_expr": str(P12.as_expr()),
        "target_P22_expr": str(P22.as_expr()),
    }


def run_dirac_ruleB_experiment() -> Dict:
    """Test Hypothesis 1 analog on Dirac-S^3 Rule B n_max=3."""
    A, labels, info = dirac_ruleB_data()
    V = info["V"]
    E = info["E"]
    deg = A.sum(axis=1)
    Q = np.diag(deg - 1)

    # Build m_j-reflection
    P = dirac_mj_reflection(labels)
    comm = P @ A - A @ P
    assert np.abs(comm).sum() == 0, "A must commute with m_j-reflection"

    # No fixed points (m_j is half-integer)
    fixed = [i for i, lab in enumerate(labels) if lab.two_m_j == 0]
    assert len(fixed) == 0, "Expected no m_j=0 fixed points"
    dim_sym_expected = V // 2
    dim_antisym_expected = V // 2

    U_sym, U_antisym = build_reflection_projectors_dirac(labels)
    assert U_sym.shape == (V, dim_sym_expected)
    assert U_antisym.shape == (V, dim_antisym_expected)

    I_sym = sp.simplify(U_sym.T * U_sym)
    I_antisym = sp.simplify(U_antisym.T * U_antisym)
    cross = sp.simplify(U_sym.T * U_antisym)
    assert I_sym == sp.eye(dim_sym_expected)
    assert I_antisym == sp.eye(dim_antisym_expected)
    assert cross == sp.zeros(dim_sym_expected, dim_antisym_expected)

    print(f"  Computing det(M_sym) (dim {dim_sym_expected}) ...")
    det_sym = block_det_I_minus_sA_plus_sQ(A, U_sym, Q)
    print(f"  Computing det(M_antisym) (dim {dim_antisym_expected}) ...")
    det_antisym = block_det_I_minus_sA_plus_sQ(A, U_antisym, Q)

    s = sp.symbols("s")
    P22, P24 = target_P22_P24_dirac()

    det_sym_factored = sp.factor(det_sym)
    det_antisym_factored = sp.factor(det_antisym)
    det_sym_poly = sp.Poly(det_sym, s)
    det_antisym_poly = sp.Poly(det_antisym, s)

    sym_contains_P22 = sp.rem(det_sym, P22.as_expr(), s) == 0
    sym_contains_P24 = sp.rem(det_sym, P24.as_expr(), s) == 0
    antisym_contains_P22 = sp.rem(det_antisym, P22.as_expr(), s) == 0
    antisym_contains_P24 = sp.rem(det_antisym, P24.as_expr(), s) == 0

    full_M = sp.eye(V) - s * sp.Matrix(A.tolist()) + s * s * sp.Matrix(Q.tolist())
    full_det = sp.expand(full_M.det())
    product = sp.expand(det_sym * det_antisym)
    product_matches = sp.simplify(product - full_det) == 0

    return {
        "graph": "Dirac-S^3 Rule B n_max=3",
        "V": V,
        "E": E,
        "beta_1": E - V + 1,
        "dim_sym": dim_sym_expected,
        "dim_antisym": dim_antisym_expected,
        "fixed_points_count": len(fixed),
        "det_sym_degree": det_sym_poly.degree(),
        "det_antisym_degree": det_antisym_poly.degree(),
        "det_sym_factored": str(det_sym_factored),
        "det_antisym_factored": str(det_antisym_factored),
        "sym_block_contains_P22": bool(sym_contains_P22),
        "sym_block_contains_P24": bool(sym_contains_P24),
        "antisym_block_contains_P22": bool(antisym_contains_P22),
        "antisym_block_contains_P24": bool(antisym_contains_P24),
        "product_equals_full_det": bool(product_matches),
        "target_P22_degree": P22.degree(),
        "target_P24_degree": P24.degree(),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(output_dir: Path = None) -> Dict:
    if output_dir is None:
        output_dir = Path(__file__).parent / "data"
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results: Dict = {
        "track": "RH-E",
        "hypothesis": "Paper 29 Hypothesis 1 — Hopf-U(1) block decomposition",
        "methodology": (
            "Z_2 m-reflection decomposition at the node-level matrix "
            "(I - sA + s^2 Q) via Ihara-Bass. The Z_2 acts as m → -m "
            "(m = m_l for scalar S^5, m = m_j for Dirac-S^3). Both A and "
            "Q commute with the reflection; the block determinants are "
            "the non-trivial factor polynomials P_d."
        ),
    }

    print("=== Scalar S^5 Bargmann-Segal N_max=3 ===")
    scalar = run_scalar_s5_experiment()
    all_results["scalar_S5"] = scalar
    print(f"  V={scalar['V']}, E={scalar['E']}, beta_1={scalar['beta_1']}")
    print(f"  dim_sym={scalar['dim_sym']} (contains P12={scalar['sym_block_contains_P12']}, "
          f"P22={scalar['sym_block_contains_P22']})")
    print(f"  dim_antisym={scalar['dim_antisym']} (contains P12={scalar['antisym_block_contains_P12']}, "
          f"P22={scalar['antisym_block_contains_P22']})")
    print(f"  product = full det? {scalar['product_equals_full_det']}")

    print("\n=== Dirac-S^3 Rule B n_max=3 ===")
    dirac = run_dirac_ruleB_experiment()
    all_results["dirac_S3_ruleB"] = dirac
    print(f"  V={dirac['V']}, E={dirac['E']}, beta_1={dirac['beta_1']}")
    print(f"  dim_sym={dirac['dim_sym']} (contains P22={dirac['sym_block_contains_P22']}, "
          f"P24={dirac['sym_block_contains_P24']})")
    print(f"  dim_antisym={dirac['dim_antisym']} (contains P22={dirac['antisym_block_contains_P22']}, "
          f"P24={dirac['antisym_block_contains_P24']})")
    print(f"  product = full det? {dirac['product_equals_full_det']}")

    # Verdict
    scalar_pass = (scalar["antisym_block_contains_P12"]
                   and scalar["sym_block_contains_P22"]
                   and scalar["product_equals_full_det"])
    dirac_pass = (dirac["sym_block_contains_P22"]
                  and dirac["antisym_block_contains_P24"]
                  and dirac["product_equals_full_det"])

    all_results["verdict"] = {
        "scalar_S5_hypothesis_validated": scalar_pass,
        "dirac_ruleB_hypothesis_validated": dirac_pass,
        "overall_verdict": (
            "VALIDATED (both graphs)"
            if scalar_pass and dirac_pass
            else "PARTIAL — see details"
        ),
    }

    out_path = output_dir / "hopf_u1_block_test.json"
    with out_path.open("w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nVerdict: {all_results['verdict']['overall_verdict']}")
    print(f"Wrote: {out_path}")
    return all_results


if __name__ == "__main__":
    main()
